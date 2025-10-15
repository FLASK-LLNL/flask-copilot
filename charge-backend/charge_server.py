from fastapi import FastAPI, Request, WebSocket, WebSocketDisconnect
from fastapi.responses import StreamingResponse, FileResponse
from fastapi.staticfiles import StaticFiles
from fastapi.middleware.cors import CORSMiddleware
import asyncio
from collections import defaultdict
import json
import os
import random
import argparse
import sys
from typing import Dict, cast

import httpx

cur_dir = os.path.dirname(os.path.abspath(__file__))
sys.path.append(os.path.join(cur_dir, "ChARGe", "experiments", "Molecule_Generation"))
from ChARGe.experiments.Molecule_Generation.LMOExperiment import (
    LMOExperiment as LeadMoleculeOptimization,
)
from ChARGe.experiments.Retrosynthesis.RetrosynthesisExperiment import (
    TemplateFreeRetrosynthesisExperiment as RetrosynthesisExperiment,
    TemplateFreeReactionOutputSchema as ReactionOutputSchema,
)

import ChARGe.experiments.Molecule_Generation.helper_funcs as lmo_helper_funcs
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
    RETROSYNTH_UNCONSTRAINED_USER_PROMPT_TEMPLATE,
    RETROSYNTH_CONSTRAINED_USER_PROMPT_TEMPLATE,
    Node,
    Edge,
)
import copy

# sys.path.append(os.path.dirname(os.path.abspath(__file__)))
# from ChARGe.experiments.LMOExperiment import LMOExperiment

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


# TODO: Convert this to a dataclass
MOLECULE_HOVER_TEMPLATE = """**SMILES:** `{}`\n
## Properties
 - Molecule Weight: {:.3f}
 - **Cost:** {:.2f}
 - **Density:** {:.3f}
 - **SA Score:** {:.3f}"""

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


class RetrosynthesisContext:
    """
    Manages a retrosynthesis experiment
    """

    def __init__(self):
        self.node_ids: dict[str, Node] = {}
        self.node_id_to_planner: dict[str, aizynth_funcs.RetroPlanner] = {}
        self.node_id_to_charge_client: dict[str, AutoGenClient] = {}
        self.azf_nodes: dict[str, aizynth_funcs.Node] = {}
        self.nodes_per_level: dict[int, int] = defaultdict(int)
        self.parents: dict[str, str] = {}

    def reset(self):
        self.node_ids.clear()
        self.node_id_to_planner.clear()
        self.node_id_to_charge_client.clear()
        self.azf_nodes.clear()
        self.nodes_per_level.clear()
        self.parents.clear()


def generate_tree_structure(
    reaction_path_dict: Dict[int, aizynth_funcs.Node],
    retro_synth_context: RetrosynthesisContext,
    start_level: int = 0,
) -> tuple[list[Node], list[Edge]]:
    """Generate nodes and edges from reaction path dict"""
    nodes = []
    edges = []

    root_id = 0
    node_queue = [(reaction_path_dict[root_id], start_level)]  # (node, level)

    while node_queue:
        current_node, level = node_queue.pop(0)
        node_id = current_node.node_id
        smiles = current_node.smiles
        purchasable = current_node.purchasable
        leaf = current_node.is_leaf
        node_id_str = f"node_{node_id}"
        hover_info = f"# Molecule \n **SMILES:** {smiles}\n"
        if leaf:
            if purchasable:
                hover_info += " - This molecule is purchasable.\n"
                # TODO: For ... chemprice
            else:
                hover_info += " - This molecule is NOT purchasable.\n"

        node = Node(
            id=node_id_str,
            smiles=smiles,
            label=smiles,
            hoverInfo=hover_info,
            level=level,
            parentId=(f"node_{current_node.parent_id}" if current_node.parent_id is not None else None),
            highlight=("red" if (leaf and not purchasable) else "normal"),
        )

        retro_synth_context.node_ids[node_id_str] = node
        retro_synth_context.azf_nodes[node_id_str] = current_node
        retro_synth_context.nodes_per_level[level] += 1
        nodes.append(node)

        if current_node.parent_id is not None:
            edge = Edge(
                id=f"edge_{current_node.parent_id}_{node_id}",
                fromNode=f"node_{current_node.parent_id}",
                toNode=node_id_str,
                status="complete",
                label=None,
            )
            edges.append(edge)
            retro_synth_context.parents[node_id_str] = f"node_{current_node.parent_id}"

        for child_id in current_node.children:
            child_node = reaction_path_dict[child_id]
            node_queue.append((child_node, level + 1))

    return nodes, edges


def calculate_positions(nodes: list[Node], y_offset: int = 0):
    """
    Calculate positions for all nodes (matching frontend logic).
    Operates in-place.
    """
    BOX_WIDTH = 270  # Must match with javascript!
    BOX_GAP = 160  # Must match with javascript!
    level_gap = BOX_WIDTH + BOX_GAP
    node_spacing = 150

    # Group by level
    levels = {}
    for node in nodes:
        level = node.level
        if level not in levels:
            levels[level] = []
        levels[level].append(node)

    # Position nodes
    for node in nodes:
        level_nodes = levels[node.level]
        index_in_level = level_nodes.index(node) + y_offset

        node.x = 100 + node.level * level_gap
        node.y = 100 + index_in_level * node_spacing


async def generate_molecules(start_smiles: str, context: RetrosynthesisContext, websocket: WebSocket):
    """Stream positioned nodes and edges"""
    logger.info(f"Planning retrosynthesis for: {start_smiles}")

    # Generate and position entire tree upfront

    root = Node(
        id="node_0",
        smiles=start_smiles,
        label=start_smiles,
        hoverInfo=f"# Root molecule \n **SMILES:** {start_smiles}",
        level=0,
        parentId=None,
        cost=None,
        bandgap=None,
        yield_=None,
        highlight="yellow",
        x=100,
        y=100,
    )
    await websocket.send_json({"type": "node", **root.json()})
    planner = context.node_id_to_planner[root.id] = aizynth_funcs.RetroPlanner(configfile=args.config_file)
    _, _, routes = planner.plan(start_smiles)
    if not routes:
        await websocket.send_json(
            {
                "type": "response",
                "message": f"No synthesis routes found for {start_smiles}.",
                "smiles": start_smiles,
            }
        )
        await websocket.send_json({"type": "complete"})
        return
    logger.info(f"Found {len(routes)} routes for {start_smiles}.")

    planner.last_route_used = 0
    reaction_path = aizynth_funcs.ReactionPath(route=routes[0])
    nodes, edges = generate_tree_structure(reaction_path.nodes, context)
    logger.info(f"Generated {len(nodes)} nodes and {len(edges)} edges.")

    calculate_positions(nodes)

    # Stream remaining nodes with edges
    for node in nodes[1:]:
        # Find edge for this node
        edge = next((e for e in edges if e.toNode == node.id), None)

        # Send node
        await websocket.send_json({"type": "node", **node.json()})
        if edge:
            await websocket.send_json({"type": "edge", **edge.json()})

    await websocket.send_json({"type": "node_update", "id": root.id, "highlight": "normal"})
    await websocket.send_json({"type": "complete"})


async def lead_molecule(
    start_smiles: str,
    experiment,
    lmo_runner,
    depth: int,
    websocket: WebSocket,
):
    """Stream positioned nodes and edges"""

    mol_file_path = args.json_file

    lead_molecule_smiles = start_smiles
    logger.info(f"Starting experiment with lead molecule: {lead_molecule_smiles}")
    parent_id = 0
    node_id = 0
    lead_molecule_data = lmo_helper_funcs.post_process_smiles(
        smiles=lead_molecule_smiles, parent_id=parent_id - 1, node_id=node_id
    )

    # Start the db with the lead molecule
    lmo_helper_funcs.save_list_to_json_file(data=[lead_molecule_data], file_path=mol_file_path)

    # Start the db with the lead molecule
    lmo_helper_funcs.save_list_to_json_file(data=[lead_molecule_data], file_path=mol_file_path)
    logger.info(f"Storing found molecules in {mol_file_path}")

    # Run the experiment in a loop
    new_molecules = lmo_helper_funcs.get_list_from_json_file(file_path=mol_file_path)  # Start with known molecules

    leader_hov = MOLECULE_HOVER_TEMPLATE.format(
        lead_molecule_smiles,
        0.0,  # TODO: Add molecule weight calculation
        0.0,  # TODO: Add cost calculation
        lead_molecule_data["density"],
        lead_molecule_data["sascore"],
    )
    node = dict(
        id=f"node_{node_id}",
        smiles=lead_molecule_smiles,
        label=f"{lead_molecule_smiles}",
        # Add property calculations here
        energy=lead_molecule_data["density"],
        level=0,
        cost=lead_molecule_data["sascore"],
        # Not sure what to put here
        hoverInfo=leader_hov,
        x=0,
        y=node_id * 150,
    )

    await websocket.send_json({"type": "node", **node})

    edge_data = {
        "type": "edge",
        "id": f"edge_{0}_{1}",
        "status": "computing",
        "label": "Optimizing",
        "fromNode": {"id": f"node_{0}", "x": 0, "y": -150},
        "toNode": {"id": f"node_{1}", "x": 200, "y": -150},
    }
    await websocket.send_json(edge_data)
    # Generate one node at a time

    mol_data = [lead_molecule_data]
    max_iterations = args.max_iterations

    for i in range(depth):
        if i > 0:
            edge_data = {
                "type": "edge",
                "id": f"edge_{0}_{1}",
                "status": "computing",
                "label": "Optimizing",
                "fromNode": {"id": f"node_{0}", "x": 0, "y": -150},
                "toNode": {"id": f"node_{1}", "x": 200, "y": -150},
            }
        await websocket.send_json(edge_data)
        # Generate new molecule

        iteration = 0
        while iteration < max_iterations:

            try:
                iteration += 1
                results = await lmo_runner.run()
                results = results.as_list()  # Convert to list of strings
                logger.info(f"New molecules generated: {results}")
                processed_mol = lmo_helper_funcs.post_process_smiles(
                    smiles=results[0], parent_id=parent_id, node_id=node_id
                )
                canonical_smiles = processed_mol["smiles"]
                if canonical_smiles not in new_molecules and canonical_smiles != "Invalid SMILES":
                    new_molecules.append(canonical_smiles)
                    mol_data.append(processed_mol)
                    lmo_helper_funcs.save_list_to_json_file(data=mol_data, file_path=mol_file_path)
                    logger.info(f"New molecule added: {canonical_smiles}")
                    node_id += 1
                    mol_hov = MOLECULE_HOVER_TEMPLATE.format(
                        canonical_smiles,
                        0.0,  # TODO: Add molecule weight calculation
                        0.0,  # TODO: Add cost calculation
                        processed_mol["density"],
                        processed_mol["sascore"],
                    )
                    node = dict(
                        id=f"node_{node_id}",
                        smiles=canonical_smiles,
                        label=f"{canonical_smiles}",
                        # Add property calculations here
                        energy=processed_mol["density"],
                        level=0,
                        cost=processed_mol["sascore"],
                        # Not sure what to put here
                        hoverInfo=mol_hov,
                        x=0,
                        y=node_id * 150,
                    )

                    await websocket.send_json({"type": "node", **node})

                    edge_data = {
                        "type": "edge",
                        "id": f"edge_{node_id}_{node_id+1}",
                        "status": "computing",
                        "label": "Optimizing",
                        "fromNode": {"id": f"node_{node_id}", "x": 0, "y": -150},
                        "toNode": {"id": f"node_{node_id+1}", "x": 200, "y": -150},
                    }

                    await websocket.send_json(edge_data)
                    experiment = LeadMoleculeOptimization(lead_molecule=canonical_smiles)
                    lmo_runner.experiment_type = experiment
                    parent_id = node_id

                    break  # Exit while loop to proceed to next node
                else:
                    logger.info(f"Duplicate molecule found: {canonical_smiles}")
                    # Continue the while loop to try generating again
            except WebSocketDisconnect:
                logger.info("WebSocket disconnected")
                raise
            except asyncio.CancelledError:
                await websocket.send_json({"type": "stopped"})
                raise  # re-raise so cancellation propagates
            except Exception as e:
                logger.error(f"Error occurred: {e}")
        if i == depth - 1:
            break
        edge_data = {
            "type": "edge",
            "id": f"edge_{0}_{1}",
            "status": "computing",
            "label": "Optimizing",
            "fromNode": {"id": f"node_{0}", "x": 0, "y": -150},
            "toNode": {"id": f"node_{1}", "x": 200, "y": -150},
        }
        await websocket.send_json(edge_data)
        # TODO: Compute here!!!
        await asyncio.sleep(0.8)

    edge_data = {
        "type": "edge",
        "id": f"edge_{0}_{1}",
        "status": "Completed",
        "label": "Optimization Complete",
        "fromNode": {"id": f"node_{0}", "x": 0, "y": -150},
        "toNode": {"id": f"node_{1}", "x": 200, "y": -150},
    }

    await websocket.send_json({"type": "complete"})


async def constrained_opt(parent_smiles, constraint_smiles, planner, websocket: WebSocket):
    """Constrained optimization using retrosynthesis"""

    await websocket.send_json(
        {
            "type": "response",
            "message": f"Finding synthesis pathway for {parent_smiles}...",
            "smiles": parent_smiles,
        }
    )
    await websocket.send_json(
        {
            "type": "response",
            "message": f"Searching for alternatives without {constraint_smiles}",
            "smiles": constraint_smiles,
        }
    )

    user_prompt = RETROSYNTH_CONSTRAINED_USER_PROMPT_TEMPLATE.format(
        target_molecule=parent_smiles, constrained_reactant=constraint_smiles
    )
    retro_experiment = RetrosynthesisExperiment(user_prompt=user_prompt)
    planner.experiment_type = retro_experiment
    logger.info(f"Optimizing {parent_smiles} without using {constraint_smiles} in the synthesis.")
    result = await planner.run()
    return result


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
                    logger.error(f"Unknown problem type: {data['problemType']}")

                async def run_task():
                    if data["problemType"] == "optimization":
                        await lead_molecule(
                            data["smiles"],
                            lmo_experiment,
                            lmo_runner,
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
                pass
            elif action == "recompute-parent-reaction":
                logger.info("Recompute parent reaction action received")
                logger.info(f"Data: {data}")
                pass

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
