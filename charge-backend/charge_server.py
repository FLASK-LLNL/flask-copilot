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
from typing import Dict

cur_dir = os.path.dirname(os.path.abspath(__file__))
sys.path.append(os.path.join(cur_dir, "ChARGe", "experiments", "Molecule_Generation"))
from ChARGe.experiments.Molecule_Generation.LMOExperiment import (
    LMOExperiment as LeadMoleculeOptimization,
)
from ChARGe.experiments.Retrosynthesis.RetrosynthesisExperiment import (
    TemplateFreeRetrosynthesisExperiment as RetrosynthesisExperiment,
)

import ChARGe.experiments.Molecule_Generation.helper_funcs as lmo_helper_funcs

import os
from charge.clients.Client import Client
from charge.clients.autogen import AutoGenClient
import charge.servers.AiZynthTools as aizynth_funcs

from loguru import logger
import sys
from backend_helper_funcs import (
    CallbackHandler,
    RETROSYNTH_UNCONSTRAINED_USER_PROMPT_TEMPLATE,
    RETROSYNTH_CONSTRAINED_USER_PROMPT_TEMPLATE,
    RetroSynthesisContext,
)

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
parser.add_argument(
    "--host", type=str, default="127.0.0.1", help="Host to run the server on"
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


def generate_tree_structure(reaction_path_dict: Dict[int, aizynth_funcs.Node]):
    """Generate nodes and edges from reaction path dict"""
    nodes = []
    edges = []

    root_id = 0
    node_queue = [(reaction_path_dict[root_id], 0)]  # (node, level)

    while node_queue:
        current_node, level = node_queue.pop(0)
        node_id = current_node.node_id
        smiles = current_node.smiles
        purchasable = current_node.purchasable
        intermediate = not (current_node.is_root or current_node.is_leaf)
        leaf = current_node.is_leaf
        node = {
            "id": node_id,
            "label": smiles,
            "smiles": smiles,
            "level": level,
            "purchasable": purchasable,
            "energy": random.uniform(0, 1),  # Placeholder for energy
            "cost": random.uniform(0, 1),  # Placeholder for cost
            "parent_id": current_node.parent_id,
            "hoverInfo": f"Purchasable: {purchasable}\n",  # Placeholder for hover info
        }
        RetroSynthesisContext.node_ids[node_id] = current_node
        RetroSynthesisContext.node_by_smiles[smiles] = current_node
        nodes.append(node)
        if current_node.parent_id is not None:
            edge = {
                "id": f"edge_{current_node.parent_id}_{node_id}",
                "from": current_node.parent_id,
                "to": node_id,
                "reactionType": " ",
            }
            edges.append(edge)
        for child_id in current_node.children:
            child_node = reaction_path_dict[child_id]
            node_queue.append((child_node, level + 1))

    return nodes, edges


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


async def generate_molecules(
    start_smiles: str, planner, charge_retro_planner, websocket: WebSocket = None
):
    """Stream positioned nodes and edges"""
    logger.info(f"Planning retrosynthesis for: {start_smiles}")

    # Generate and position entire tree upfront
    tree, stats, routes = planner.plan(start_smiles)
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

    reaction_path = aizynth_funcs.ReactionPath(route=routes[0])
    nodes, edges = generate_tree_structure(reaction_path.nodes)
    logger.info(f"Generated {len(nodes)} nodes and {len(edges)} edges.")

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

    # Check if all leaf nodes are purchasable
    leaf_nodes = reaction_path.leaf_nodes
    for leaf_node_id in leaf_nodes:
        leaf_node = reaction_path.nodes[leaf_node_id]
        if not leaf_node.purchasable:
            await websocket.send_json(
                {
                    "type": "response",
                    "message": f"Leaf molecule {leaf_node.smiles} is not purchasable.",
                    "smiles": leaf_node.smiles,
                }
            )
            user_prompt = RETROSYNTH_UNCONSTRAINED_USER_PROMPT_TEMPLATE.format(
                target_molecule=start_smiles
            )
            retro_experiment = RetrosynthesisExperiment(user_prompt=user_prompt)
            charge_retro_planner.experiment_type = retro_experiment
            logger.info(
                f"{leaf_node.smiles} is not purchasable. "
                "Switching to template-free retrosynthesis."
            )
            await websocket.send_json(
                {"type": "response", "message": "Searching for alternatives..."}
            )

            results = await charge_retro_planner.run()

            results = results.as_list()  # Convert to list of strings
            logger.info(f"New retrosynthetic step found: {results}")
            for res in results:
                await websocket.send_json(
                    {
                        "type": "response",
                        "message": f"Proposed step: {res}",
                        "smiles": res,
                    }
                )
            ## TODO: Call AiZynthFinder again with the new reactants

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
    lead_molecule_data = lmo_helper_funcs.post_process_smiles(
        smiles=lead_molecule_smiles, parent_id=parent_id - 1, node_id=node_id
    )

    # Start the db with the lead molecule
    lmo_helper_funcs.save_list_to_json_file(
        data=[lead_molecule_data], file_path=mol_file_path
    )

    # Start the db with the lead molecule
    lmo_helper_funcs.save_list_to_json_file(
        data=[lead_molecule_data], file_path=mol_file_path
    )
    logger.info(f"Storing found molecules in {mol_file_path}")

    # Run the experiment in a loop
    new_molecules = lmo_helper_funcs.get_list_from_json_file(
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
                processed_mol = lmo_helper_funcs.post_process_smiles(
                    smiles=results[0], parent_id=parent_id, node_id=node_id
                )
                canonical_smiles = processed_mol["smiles"]
                if (
                    canonical_smiles not in new_molecules
                    and canonical_smiles != "Invalid SMILES"
                ):
                    new_molecules.append(canonical_smiles)
                    mol_data.append(processed_mol)
                    lmo_helper_funcs.save_list_to_json_file(
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


async def unconstrained_opt(parent_smiles, planner, websocket: WebSocket):
    """Unconstrained optimization using retrosynthesis"""

    await websocket.send_json(
        {
            "type": "response",
            "message": f"Optimizing {parent_smiles}...",
            "smiles": parent_smiles,
        }
    )

    user_prompt = RETROSYNTH_UNCONSTRAINED_USER_PROMPT_TEMPLATE.format(
        target_molecule=parent_smiles
    )
    retro_experiment = RetrosynthesisExperiment(user_prompt=user_prompt)
    planner.experiment_type = retro_experiment
    logger.info(f"Optimizing {parent_smiles} using retrosynthesis.")
    result = await planner.run()
    return result


async def constrained_opt(
    parent_smiles, constraint_smiles, planner, websocket: WebSocket
):
    """Constrained optimization using retrosynthesis"""

    await websocket.send_json(
        {
            "type": "response",
            "message": f"Optimizing {parent_smiles}...",
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
    logger.info(
        f"Optimizing {parent_smiles} without using {constraint_smiles} in the synthesis."
    )
    result = await planner.run()
    return result


async def regenerate_sub_tree(result, starting_node, planner, websocket: WebSocket):
    """Regenerate sub-tree from a given node"""

    pass


async def optimize_molecule_retro(
    node_id, opt_type, experiment, planner, websocket: WebSocket
):
    """Optimize a molecule using retrosynthesis by node ID"""
    current_node = RetroSynthesisContext.node_ids.get(node_id)

    assert current_node is not None, f"Node ID {node_id} not found"

    if opt_type == "product_optimization_retro":
        result = await unconstrained_opt(current_node.smiles, planner, websocket)
        starting_node = current_node
    else:
        parent_id = current_node.parent_id
        parent_node = RetroSynthesisContext.node_by_smiles.get(parent_id)
        assert parent_node is not None, f"Parent node {parent_id} not found"
        result = await constrained_opt(
            parent_node.smiles, current_node.smiles, planner, websocket
        )
        starting_node = parent_node

    await regenerate_sub_tree(result, starting_node, planner, websocket)


@app.websocket("/ws")
async def websocket_endpoint(websocket: WebSocket):
    global CURRENT_TASK
    await websocket.accept()

    # Initialize Charge experiment with a dummy lead molecule
    lmo_experiment = None
    lmo_runner = None

    retro_experiment = None
    retro_runner = None

    aizynthfinder_planner = aizynth_funcs.RetroPlanner(configfile=args.config_file)
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

                if data["problemType"] == "optimization":
                    lmo_experiment = LeadMoleculeOptimization(
                        lead_molecule=data["smiles"]
                    )
                    if lmo_runner is None:
                        lmo_runner = AutoGenClient(
                            experiment_type=lmo_experiment,
                            model=model,
                            backend=backend,
                            api_key=API_KEY,
                            model_kwargs=kwargs,
                            server_url=server_urls,  # TODO: Change this to LMFO specific server
                            thoughts_callback=CallbackHandler(websocket),
                        )
                    else:
                        lmo_runner.experiment_type = lmo_experiment
                elif data["problemType"] in [
                    "product_optimization_retro",
                    "reactant_optimization_retro",
                    "retrosynthesis",
                ]:
                    if data["problemType"] == "product_optimization_retro":
                        # We will set this up in the optimize_molecule_retro function
                        retro_experiment = None
                    elif data["problemType"] == "reactant_optimization_retro":
                        # We will set this up in the optimize_molecule_retro function
                        retro_experiment = None
                    else:
                        # Set up retrosynthesis experiment to retrosynthesis
                        # to ensure the reactant is not used in the synthesis
                        retro_experiment = RetrosynthesisExperiment(
                            user_prompt=RETROSYNTH_UNCONSTRAINED_USER_PROMPT_TEMPLATE.format(
                                target_molecule=data["smiles"]
                            )
                        )

                    if retro_runner is None:
                        retro_runner = AutoGenClient(
                            experiment_type=retro_experiment,
                            model=model,
                            backend=backend,
                            api_key=API_KEY,
                            model_kwargs=kwargs,
                            server_url=server_urls,  # TODO: Change this to a retrosynthesis specific server
                            thoughts_callback=CallbackHandler(websocket),
                        )
                    else:
                        retro_runner.experiment_type = retro_experiment

                async def run_task():
                    if data["problemType"] == "optimization":
                        await lead_molecule(
                            data["smiles"],
                            lmo_experiment,
                            lmo_runner,
                            data.get("depth", 3),
                            websocket,
                        )
                    elif data["problemType"] in [
                        "product_optimization_retro",
                        "reactant_optimization_retro",
                    ]:
                        await optimize_molecule_retro(
                            data["id"],
                            data["problemType"],
                            retro_experiment,
                            retro_runner,
                            websocket,
                        )
                    else:
                        await generate_molecules(
                            data["smiles"],
                            aizynthfinder_planner,
                            retro_runner,
                            websocket,
                        )

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
                RetroSynthesisContext.reset()
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

    uvicorn.run(app, host=args.host, port=args.port)
