from functools import partial
from fastapi import FastAPI, WebSocket, WebSocketDisconnect
from fastapi.responses import FileResponse
from fastapi.staticfiles import StaticFiles
from fastapi.middleware.cors import CORSMiddleware
import asyncio
from collections import defaultdict
import json
import os
import argparse
import sys
from typing import Dict, cast

import httpx

from charge.experiments.LMOExperiment import (
    LMOExperiment as LeadMoleculeOptimization,
)
from charge.experiments.RetrosynthesisExperiment import (
    TemplateFreeRetrosynthesisExperiment as RetrosynthesisExperiment,
    TemplateFreeReactionOutputSchema as ReactionOutputSchema,
)

import charge.utils.helper_funcs as lmo_helper_funcs
from charge.servers.server_utils import try_get_public_hostname

import os
from charge.clients.Client import Client
from charge.clients.autogen import AutoGenClient
import charge.servers.AiZynthTools as aizynth_funcs
import logging
from aizynthfinder.utils.logging import setup_logger

setup_logger(console_level=logging.INFO)

from loguru import logger
from callback_logger import callback_logger

import sys
from backend_helper_funcs import (
    CallbackHandler,
    RetrosynthesisContext,
    Node,
    Edge,
    calculate_positions,
)
from retro_charge_backend_funcs import (
    generate_molecules,
    optimize_molecule_retro,
)
import copy
from lmo_charge_backend_funcs import lead_molecule
from aizynth_backend_funcs import aizynth_retro
from retro_charge_backend_funcs import (
    unconstrained_retro,
    constrained_retro,
    get_constrained_prompt,
    get_unconstrained_prompt,
)


parser = argparse.ArgumentParser()

parser.add_argument(
    "--json_file",
    type=str,
    default="known_molecules.json",
    help="Path to the JSON file containing known molecules.",
)

parser.add_argument(
    "--max_iterations",
    type=int,
    default=5,
)

parser.add_argument(
    "--config-file",
    type=str,
    default="config.yml",
    help="Path to the configuration file for AiZynthFinder.",
)
parser.add_argument("--port", type=int, default=8001, help="Port to run the server on")
parser.add_argument("--host", type=str, default=None, help="Host to run the server on")

parser.add_argument("--lmo-urls", nargs="*", type=str, default=[])
parser.add_argument("--retro-urls", nargs="*", type=str, default=[])

# Add standard CLI arguments
Client.add_std_parser_arguments(parser)

args = parser.parse_args()


app = FastAPI()

# CORS for development
app.add_middleware(
    CORSMiddleware,
    allow_origins=["*"],
    allow_credentials=True,
    allow_methods=["*"],
    allow_headers=["*"],
)

if 'FLASK_APPDIR' in os.environ:
    DIST_PATH = os.environ['FLASK_APPDIR']
else:
    DIST_PATH = os.path.join(os.path.dirname(__file__), "flask-app", "dist")
ASSETS_PATH = os.path.join(DIST_PATH, "assets")


(MODEL, BACKEND, API_KEY, MODEL_KWARGS) = AutoGenClient.configure(args.model, args.backend)

server_urls = args.server_urls
assert server_urls is not None, "Server URLs must be provided"
for url in server_urls + args.lmo_urls + args.retro_urls:
    assert url.endswith("/sse"), f"Server URL {url} must end with /sse"

LMO_URLS = args.lmo_urls + server_urls
RETRO_URLS = args.retro_urls + server_urls


if os.path.exists(ASSETS_PATH):
    # Serve the frontend
    app.mount("/assets", StaticFiles(directory=ASSETS_PATH), name="assets")
    app.mount("/rdkit", StaticFiles(directory=os.path.join(DIST_PATH, "rdkit")), name="rdkit")

    @app.get("/")
    async def root():
        return FileResponse(os.path.join(DIST_PATH, "index.html"))


def make_client(client, experiment, server_urls, websocket):
    if client is None:
        return AutoGenClient(
            experiment_type=experiment,
            model=MODEL,
            backend=BACKEND,
            api_key=API_KEY,
            model_kwargs=MODEL_KWARGS,
            server_url=server_urls,
            thoughts_callback=CallbackHandler(websocket),
        )
    else:
        client.experiment_type = experiment
        return client

@app.websocket("/ws")
async def websocket_endpoint(websocket: WebSocket):
    await websocket.accept()

    # Keep track of currently running task
    CURRENT_TASK: asyncio.Task | None = None
    # Initialize Charge experiment with a dummy lead molecule
    lmo_experiment = None
    lmo_runner = None

    clogger = callback_logger(websocket)

    retro_synth_context: RetrosynthesisContext | None = None
    try:
        while True:
            data = await websocket.receive_json()
            action = data.get("action")

            if action in [
                "compute",
                "optimize-from",
                "compute-reaction-from",
                "recompute-reaction",
                "recompute-parent-reaction",
            ]:
                # cancel any running task first
                if CURRENT_TASK and not CURRENT_TASK.done():
                    logger.info("Cancelling existing compute task...")
                    CURRENT_TASK.cancel()
                    try:
                        await CURRENT_TASK
                    except asyncio.CancelledError:
                        logger.info("Previous compute task cancelled.")

            if action == "compute":

                if data["problemType"] == "optimization":

                    # Task to optimize lead molecule using LMO
                    clogger.info("Start Optimization action received")
                    logger.info(f"Data: {data}")

                    lmo_experiment = LeadMoleculeOptimization(lead_molecule=data["smiles"])
                    if lmo_runner is None:
                        lmo_runner = AutoGenClient(
                            experiment_type=lmo_experiment,
                            model=MODEL,
                            backend=BACKEND,
                            api_key=API_KEY,
                            model_kwargs=MODEL_KWARGS,
                            server_url=LMO_URLS,
                            thoughts_callback=CallbackHandler(websocket),
                        )
                    else:
                        lmo_runner.experiment_type = lmo_experiment

                    run_func = partial(
                        lead_molecule,
                        data["smiles"],
                        lmo_experiment,
                        lmo_runner,
                        args.json_file,
                        args.max_iterations,
                        data.get("depth", 3),
                        websocket,
                    )
                elif data["problemType"] == "retrosynthesis":
                    # Set up retrosynthesis experiment to retrosynthesis
                    # to ensure the reactant is not used in the synthesis
                    clogger.info("Setting up retrosynthesis experiment...")
                    logger.info(f"Data: {data}")

                    if retro_synth_context is None:
                        retro_synth_context = RetrosynthesisContext()

                    assert retro_synth_context is not None
                    run_func = partial(
                        generate_molecules,
                        data["smiles"],
                        args.config_file,
                        retro_synth_context,
                        websocket,
                    )
                else:
                    raise ValueError(f"Unknown problem type: {data['problemType']}")

                async def run_task():
                    await run_func()

                # start a new task
                CURRENT_TASK = asyncio.create_task(run_task())

            elif action == "compute-reaction-from":
                assert retro_synth_context is not None

                # Leaf node optimization
                logger.info("Synthesize tree leaf action received")
                logger.info(f"Data: {data}")
                await optimize_molecule_retro(data["nodeId"], retro_synth_context, websocket, MODEL, BACKEND, API_KEY, MODEL_KWARGS, RETRO_URLS)
                await websocket.send_json({"type": "complete"})
            elif action == "optimize-from":
                # Leaf node optimization
                prompt = data.get("query", None)
                if prompt:
                    await websocket.send_json(
                        {
                            "type": "response",
                            "message": {
                                "source": "System",
                                "message": f"Processing optimization query: {data['query']} for node {data['nodeId']}",
                            }
                        }
                    )
                logger.info("Optimize from action received")
                logger.info(f"Data: {data}")

                pass
            elif action == "recompute-reaction":
                prompt = data.get("query", None)
                if prompt:
                    await websocket.send_json(
                        {
                            "type": "response",
                            "message": {
                                "source": "System",
                                "message": f"Processing reaction query: {data['query']} for node {data['nodeId']}",
                            }
                        }
                    )

                logger.info("Recompute reaction action received")
                logger.info(f"Data: {data}")
                await websocket.send_json({"type": "complete"})

            elif action == "custom_query":
                await websocket.send_json(
                    {
                        "type": "response",
                        "message": {
                            "source": "System",
                            "message": f"Processing query: {data['query']} for node {data['nodeId']}",
                        }
                    }
                )
                await asyncio.sleep(3)
                await websocket.send_json({"type": "complete"})

            elif action == "reset":
                if lmo_runner:
                    lmo_runner.reset()
                if retro_synth_context is not None:
                    retro_synth_context.reset()
                logger.info("Experiment state has been reset.")

            elif action == "stop":
                if CURRENT_TASK and not CURRENT_TASK.done():
                    logger.info("Stopping current task as per user request.")
                    CURRENT_TASK.cancel()
                    try:
                        await CURRENT_TASK
                    except asyncio.CancelledError:
                        logger.info("Current task cancelled successfully.")

                else:
                    logger.info(f"Is Current task done: {CURRENT_TASK}")
                    if CURRENT_TASK:
                        logger.info(f"Current Task states: {CURRENT_TASK.done()}")
            else:
                logger.warning(f"Unknown action received: {action}")

    except WebSocketDisconnect:
        pass
    except httpx.ConnectError as e:
        logger.error(f"Connection error: {e}")
    except Exception as e:
        logger.error(f"Error in WebSocket connection: {e}")


if __name__ == "__main__":
    import uvicorn

    host = args.host
    if host is None:
        _, host = try_get_public_hostname()

    uvicorn.run(app, host=host, port=args.port)
