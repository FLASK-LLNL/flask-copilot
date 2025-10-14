from functools import partial
from fastapi import FastAPI,  WebSocket, WebSocketDisconnect
from fastapi.responses import FileResponse
from fastapi.staticfiles import StaticFiles
from fastapi.middleware.cors import CORSMiddleware
import asyncio
import os
import argparse
import sys

import httpx

cur_dir = os.path.dirname(os.path.abspath(__file__))
sys.path.append(os.path.join(cur_dir, "ChARGe", "experiments", "Molecule_Generation"))
from ChARGe.experiments.Molecule_Generation.LMOExperiment import (
    LMOExperiment as LeadMoleculeOptimization,
)
from ChARGe.experiments.Retrosynthesis.RetrosynthesisExperiment import (
    TemplateFreeRetrosynthesisExperiment as RetrosynthesisExperiment,
)


from ChARGe.charge.servers.server_utils import try_get_public_hostname

import os
from charge.clients.Client import Client
from charge.clients.autogen import AutoGenClient
import charge.servers.AiZynthTools as aizynth_funcs
import logging
from aizynthfinder.utils.logging import setup_logger

setup_logger(console_level=logging.INFO)

from loguru import logger
import sys
from backend_helper_funcs import (
    CallbackHandler,
    RetroSynthesisContext,
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

# Add standard CLI arguments
Client.add_std_parser_arguments(parser)

args = parser.parse_args()
<<<<<<< HEAD

=======
>>>>>>> bfb777c (Streamlined implementation of the Charge backend server)


app = FastAPI()

# CORS for development
app.add_middleware(
    CORSMiddleware,
    allow_origins=["*"],
    allow_credentials=True,
    allow_methods=["*"],
    allow_headers=["*"],
)

BUILD_PATH = os.path.join(os.path.dirname(__file__), "flask-app", "build")
STATIC_PATH = os.path.join(BUILD_PATH, "static")


(model, backend, API_KEY, kwargs) = AutoGenClient.configure(args.model, args.backend)

server_urls = args.server_urls
assert server_urls is not None, "Server URLs must be provided"
for url in server_urls:
    assert url.endswith("/sse"), f"Server URL {url} must end with /sse"


if os.path.exists(STATIC_PATH):
    # Serve the frontend
    app.mount("/static", StaticFiles(directory=STATIC_PATH), name="static")

    @app.get("/")
    async def root():
        return FileResponse(os.path.join(BUILD_PATH, "index.html"))


def make_client(client, experiment, server_urls, websocket):
    if client is None:
        return AutoGenClient(
            experiment_type=experiment,
            model=model,
            backend=backend,
            api_key=API_KEY,
            model_kwargs=kwargs,
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

    retro_experiment = None
    retro_runner = None

    aizynthfinder_planner = aizynth_funcs.RetroPlanner(configfile=args.config_file)
    retro_synth_context = None
    try:
        while True:
            data = await websocket.receive_json()
            action = data.get("action")


            if action in [
                "compute",
                "compute-from",
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
                
                assert data["problemType"] in [
                    "optimization",
                    "retrosynthesis",
                ], f"Unknown problem type: {data['problemType']}"

                if action == "compute" and data["problemType"] == "optimization":
                    lmo_experiment = LeadMoleculeOptimization(
                        lead_molecule=data["smiles"]
                    )

                    lmo_runner = make_client(
                        lmo_runner, lmo_experiment, server_urls, websocket
                    )

                    # Task to optimize lead molecule using LMO
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

                elif action == "compute" and data["problemType"] == "retrosynthesis":
                    # Task to generate retrosynthesis tree using AiZynthFinder
                    run_func = partial(
                        aizynth_retro,
                        data["smiles"],
                        aizynthfinder_planner,
                        retro_synth_context if retro_synth_context else RetroSynthesisContext(),
                        websocket,
                    )
                else:

                    if action == "compute-from" or action == "recompute-reaction":
                        logger.info("Compute from action received")
                        logger.info(f"Data: {data}")

                        smiles = data["smiles"]
                        node_id = data["nodeId"]

                        user_prompt = get_unconstrained_prompt(smiles)
                        retro_experiment = RetrosynthesisExperiment(
                            user_prompt=user_prompt
                        )
                        retro_runner = make_client(
                            retro_runner, retro_experiment, server_urls, websocket
                        )
                        run_func = partial(unconstrained_retro,
                                            node_id,
                                            smiles,
                                            retro_runner,
                                            aizynthfinder_planner,
                                            retro_synth_context if retro_synth_context else RetroSynthesisContext(),
                                            websocket
                                            )
                    elif action == "recompute-parent-reaction":
                        logger.info("Recompute parent reaction action received")
                        logger.info(f"Data: {data}")

                        smiles = data["smiles"]
                        node_id = data["nodeId"]
                        assert retro_synth_context is not None, "Retro synthesis context is not initialized"
                        parent_node = retro_synth_context.get_parent(node_id)

                        assert (
                            parent_node is not None
                        ), f"Parent node for {node_id} not found"

                        parent_id = parent_node.id
                        parent_smiles = parent_node.smiles
                        constraint_id = node_id
                        constraint_smiles = smiles

                        user_prompt = get_constrained_prompt(
                            parent_smiles, constraint_smiles
                        )

                        retro_experiment = RetrosynthesisExperiment(
                            user_prompt=user_prompt
                        )
                        retro_runner = make_client(
                            retro_runner, retro_experiment, server_urls, websocket
                        )

                        run_func = partial(
                            constrained_retro,
                            parent_id,
                            parent_smiles,
                            constraint_id,
                            constraint_smiles,
                            retro_runner,
                            aizynthfinder_planner
                            retro_synth_context,
                            websocket,
                        )
                    else:
                        raise ValueError(f"Unknown action: {action}")

                async def run_task():
                    await run_func()

                # start a new task
                CURRENT_TASK = asyncio.create_task(run_task())

            elif action == "custom_query":
                await websocket.send_json(
                    {
                        "type": "response",
                        "message": f"Processing query: {data['query']} for node {data['nodeId']}",
                    }
                )
                await asyncio.sleep(3)
                await websocket.send_json({"type": "complete"})

            elif action == "reset":
                if lmo_runner:
                    lmo_runner.reset()
                if retro_runner:
                    retro_runner.reset()
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
