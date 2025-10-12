from fastapi import FastAPI, Request, WebSocket, WebSocketDisconnect
from fastapi.responses import StreamingResponse, FileResponse
from fastapi.staticfiles import StaticFiles
from fastapi.middleware.cors import CORSMiddleware
import asyncio
import json
import os
import random
import argparse
import sys

cur_dir = os.path.dirname(os.path.abspath(__file__))
sys.path.append(os.path.join(cur_dir, "ChARGe", "experiments", "Molecule_Generation"))
from ChARGe.experiments.Molecule_Generation.LMOExperiment import (
    LMOExperiment as LeadMoleculeOptimization,
)
import ChARGe.experiments.Molecule_Generation.helper_funcs as helper_funcs
import os
from charge.clients.Client import Client
from charge.clients.autogen import AutoGenClient
import charge.servers.AiZynthTools as aizynth_funcs

from loguru import logger
import sys

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

# Add standard CLI arguments
Client.add_std_parser_arguments(parser)

args = parser.parse_args()
# Keep track of currently running task
CURRENT_TASK: asyncio.Task | None = None

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


def calculate_positions(nodes):
    """Calculate positions for all nodes (matching frontend logic)"""
    BOX_WIDTH = 270  # Must match with javascript!
    BOX_GAP = 160  # Must match with javascript!
    level_gap = BOX_WIDTH + BOX_GAP
    node_spacing = 150

    # Group by level
    levels = {}
    for node in nodes:
        level = node["level"]
        if level not in levels:
            levels[level] = []
        levels[level].append(node)

    # Position nodes
    positioned = []
    for node in nodes:
        level_nodes = levels[node["level"]]
        index_in_level = level_nodes.index(node)

        positioned_node = {
            **node,
            "x": 100 + node["level"] * level_gap,
            "y": 100 + index_in_level * node_spacing,
        }
        positioned.append(positioned_node)

    return positioned


async def generate_molecules(start_smiles: str, planner, websocket: WebSocket = None):
    """Stream positioned nodes and edges"""
    logger.info(f"Planning retrosynthesis for: {start_smiles}")

    # Generate and position entire tree upfront
    tree, stats, routes = planner.plan(start_smiles)
    reaction_path = aizynth_funcs.ReactionPath(route=routes[0])
    nodes, edges = generate_tree_structure(reaction_path.nodes)

    positioned_nodes = calculate_positions(nodes)

    # Create node map
    node_map = {node["id"]: node for node in positioned_nodes}

    # Stream root first
    root = positioned_nodes[0]
    await websocket.send_json({"type": "node", **root})
    await asyncio.sleep(0.8)

    # Stream remaining nodes with edges
    for i in range(1, len(positioned_nodes)):
        node = positioned_nodes[i]

        # Find edge for this node
        edge = next((e for e in edges if e["to"] == node["id"]), None)

        if edge:
            # Send edge with computing status
            edge_data = {
                "type": "edge",
                **edge,
                "status": "computing",
                "label": f"Computing: {edge['reactionType']}",
                "fromNode": node_map[edge["from"]],
                "toNode": node,
            }
            await websocket.send_json(edge_data)

            await asyncio.sleep(0.6)

            # Send node
            await websocket.send_json({"type": "node", **node})

            # Update edge to complete
            edge_complete = {
                "type": "edge_update",
                "id": edge["id"],
                "status": "complete",
                "label": edge["reactionType"],
                "fromNode": node_map[edge["from"]],
                "toNode": node,
            }
            await websocket.send_json(edge_complete)

            await asyncio.sleep(0.2)

    await websocket.send_json({"type": "complete"})


async def lead_molecule(
    start_smiles: str,
    experiment,
    lmo_runner,
    depth: int = 3,
    websocket: WebSocket = None,
):
    """Stream positioned nodes and edges"""

    mol_file_path = args.json_file

    lead_molecule_smiles = start_smiles
    logger.info(f"Starting experiment with lead molecule: {lead_molecule_smiles}")
    parent_id = 0
    node_id = 0
    lead_molecule_data = helper_funcs.post_process_smiles(
        smiles=lead_molecule_smiles, parent_id=parent_id - 1, node_id=node_id
    )

    # Start the db with the lead molecule
    helper_funcs.save_list_to_json_file(
        data=[lead_molecule_data], file_path=mol_file_path
    )

    # Start the db with the lead molecule
    helper_funcs.save_list_to_json_file(
        data=[lead_molecule_data], file_path=mol_file_path
    )
    logger.info(f"Storing found molecules in {mol_file_path}")

    # Run the experiment in a loop
    new_molecules = helper_funcs.get_list_from_json_file(
        file_path=mol_file_path
    )  # Start with known molecules

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
                processed_mol = helper_funcs.post_process_smiles(
                    smiles=results[0], parent_id=parent_id, node_id=node_id
                )
                canonical_smiles = processed_mol["smiles"]
                if (
                    canonical_smiles not in new_molecules
                    and canonical_smiles != "Invalid SMILES"
                ):
                    new_molecules.append(canonical_smiles)
                    mol_data.append(processed_mol)
                    helper_funcs.save_list_to_json_file(
                        data=mol_data, file_path=mol_file_path
                    )
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
                    experiment = LeadMoleculeOptimization(
                        lead_molecule=canonical_smiles
                    )
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


@app.websocket("/ws")
async def websocket_endpoint(websocket: WebSocket):
    global CURRENT_TASK
    await websocket.accept()

    # Initialize Charge experiment with a dummy lead molecule
    lmo_experiment = None
    lmo_runner = None
    try:
        while True:
            data = await websocket.receive_json()
            action = data.get("action")

            if action == "compute":
                # cancel any running task first
                if CURRENT_TASK and not CURRENT_TASK.done():
                    logger.info("Cancelling existing compute task...")
                    CURRENT_TASK.cancel()
                    try:
                        await CURRENT_TASK
                    except asyncio.CancelledError:
                        logger.info("Previous compute task cancelled.")

                lmo_experiment = LeadMoleculeOptimization(lead_molecule=data["smiles"])

                if lmo_runner is None:
                    lmo_runner = AutoGenClient(
                        experiment_type=lmo_experiment,
                        model=model,
                        backend=backend,
                        api_key=API_KEY,
                        model_kwargs=kwargs,
                        server_url=server_urls,
                    )
                else:
                    lmo_runner.experiment_type = lmo_experiment

                async def run_task():
                    if data["problemType"] == "optimization":
                        await lead_molecule(
                            data["smiles"],
                            lmo_experiment,
                            lmo_runner,
                            data.get("depth", 3),
                            websocket,
                        )
                    else:
                        planner = aizynth_funcs.RetroPlanner()
                        await generate_molecules(data["smiles"], planner, websocket)

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

    except WebSocketDisconnect:
        pass


if __name__ == "__main__":
    import uvicorn

    uvicorn.run(app, host="127.0.0.1", port=8001)
