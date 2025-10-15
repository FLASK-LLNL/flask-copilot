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
import sys
from backend_helper_funcs import (
    CallbackHandler,
    RetrosynthesisContext,
    Node,
    Edge,
)
from retro_charge_backend_funcs import (
    RETROSYNTH_UNCONSTRAINED_USER_PROMPT_TEMPLATE,
    RETROSYNTH_CONSTRAINED_USER_PROMPT_TEMPLATE,
    generate_molecules,
)
import copy
from lmo_charge_backend_funcs import lead_molecule
from aizynth_backend_funcs import aizynth_retro, calculate_positions
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

BUILD_PATH = os.path.join(os.path.dirname(__file__), "flask-app", "build")
STATIC_PATH = os.path.join(BUILD_PATH, "static")


(MODEL, BACKEND, API_KEY, MODEL_KWARGS) = AutoGenClient.configure(args.model, args.backend)

server_urls = args.server_urls
assert server_urls is not None, "Server URLs must be provided"
for url in server_urls + args.lmo_urls + args.retro_urls:
    assert url.endswith("/sse"), f"Server URL {url} must end with /sse"

LMO_URLS = args.lmo_urls + server_urls
RETRO_URLS = args.retro_urls + server_urls


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


async def highlight_node(node: Node, websocket: WebSocket, highlight: bool):
    await websocket.send_json({"type": "node_update", "id": node.id, "highlight": "yellow" if highlight else "normal"})


async def optimize_molecule_retro(node_id: str, context: RetrosynthesisContext, websocket: WebSocket):
    """Optimize a molecule using retrosynthesis by node ID"""
    current_node = context.node_ids.get(node_id)
    assert current_node is not None, f"Node ID {node_id} not found"

    await websocket.send_json(
        {
            "type": "response",
            "message": f"Finding synthesis pathway to {current_node.smiles}...",
            "smiles": current_node.smiles,
        }
    )

    if node_id in context.node_id_to_charge_client:
        # Existing context
        runner = context.node_id_to_charge_client[node_id]
    else:
        # New context
        user_prompt = RETROSYNTH_UNCONSTRAINED_USER_PROMPT_TEMPLATE.format(target_molecule=current_node.smiles)
        user_prompt += "\nDouble check the reactants with the `predict_reaction_products` tool to see if the products are equivalent to the given product. If there is any inconsistency (canonicalize both sides of the equation first), log it and try some other set of reactants."
        retro_experiment = RetrosynthesisExperiment(user_prompt=user_prompt)
        runner = AutoGenClient(
            experiment_type=retro_experiment,
            model=MODEL,
            backend=BACKEND,
            api_key=API_KEY,
            model_kwargs=MODEL_KWARGS,
            server_url=RETRO_URLS,
            thoughts_callback=CallbackHandler(websocket),
        )
        context.node_id_to_charge_client[node_id] = runner

    logger.info(f"Optimizing {current_node.smiles} using retrosynthesis.")

    # Run experiment
    await highlight_node(current_node, websocket, True)
    result: ReactionOutputSchema = cast(ReactionOutputSchema, await runner.run())
    await highlight_node(current_node, websocket, False)

    level = current_node.level + 1
    num_nodes = len(context.node_ids)

    nodes: list[Node] = []
    edges: list[Edge] = []
    for i, smiles in enumerate(result.reactants_smiles_list):
        node = Node(f"node_{num_nodes+i}", smiles, smiles, "Discovered", level, current_node.id)
        nodes.append(node)
        context.node_ids[node.id] = node
        edges.append(Edge(f"edge_{node_id}_{node.id}", node_id, node.id, "complete"))

    calculate_positions(nodes, context.nodes_per_level[level])
    context.nodes_per_level[level] += len(nodes)

    for node, edge in zip(nodes, edges):
        await websocket.send_json({"type": "node", **node.json()})
        await websocket.send_json({"type": "edge", **edge.json()})


@app.websocket("/ws")
async def websocket_endpoint(websocket: WebSocket):
    await websocket.accept()

    # Keep track of currently running task
    CURRENT_TASK: asyncio.Task | None = None
    # Initialize Charge experiment with a dummy lead molecule
    lmo_experiment = None
    lmo_runner = None

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
                elif data["problemType"] == "retrosynthesis":
                    # Set up retrosynthesis experiment to retrosynthesis
                    # to ensure the reactant is not used in the synthesis
                    logger.info("Setting up retrosynthesis experiment...")

                    if retro_synth_context is None:
                        retro_synth_context = RetrosynthesisContext()
                else:
                    raise ValueError(f"Unknown action: {action}")

                async def run_task():
                    if data["problemType"] == "optimization":
                        await lead_molecule(
                            data["smiles"],
                            lmo_experiment,
                            lmo_runner,
                            args.json_file,
                            args.max_iterations,
                            data.get("depth", 3),
                            websocket,
                        )
                    elif data["problemType"] == "retrosynthesis":
                        assert retro_synth_context is not None
                        await generate_molecules(
                            data["smiles"],
                            retro_synth_context,
                            websocket,
                        )
                    else:
                        logger.error(f"Unknown problem type: {data['problemType']}")

                # start a new task
                CURRENT_TASK = asyncio.create_task(run_task())

            elif action == "compute-reaction-from":
                assert retro_synth_context is not None

                # Leaf node optimization
                logger.info("Synthesize tree leaf action received")
                logger.info(f"Data: {data}")
                await optimize_molecule_retro(data["nodeId"], retro_synth_context, websocket)
                await websocket.send_json({"type": "complete"})
            elif action == "optimize-from":
                # Leaf node optimization
                prompt = data.get("query", None)
                if prompt:
                    await websocket.send_json(
                        {
                            "type": "response",
                            "message": f"Processing optimization query: {data['query']} for node {data['nodeId']}",
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
                            "message": f"Processing reaction query: {data['query']} for node {data['nodeId']}",
                        }
                    )

                logger.info("Recompute reaction action received")
                logger.info(f"Data: {data}")
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
