"""
WebSocket server that can send mock data to the FLASK Copilot web UI.
Will serve the web app, if exists.

Supported messages from server to frontend:
    * ``node``: New node
    * ``edge``: New edge
    * ``edge_update``: Update existing edge properties
    * ``complete``: Free up UI for user input
    * ``response``: Server response to user prompt
    * ``error``: Server error message

Supported messages from frontend to server:
    * ``compute``: Start the given computation
    * ``custom_query``: Execute custom user query
"""

from fastapi import FastAPI, Request, WebSocket, WebSocketDisconnect
from fastapi.responses import StreamingResponse, FileResponse
from fastapi.staticfiles import StaticFiles
from fastapi.middleware.cors import CORSMiddleware
import asyncio
import os
import random
from typing import Any

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

if os.path.exists(STATIC_PATH):
    # Serve the frontend
    app.mount("/static", StaticFiles(directory=STATIC_PATH), name="static")

    @app.get("/")
    async def root():
        return FileResponse(os.path.join(BUILD_PATH, "index.html"))


def generate_tree_structure(start_smiles: str, depth: int = 3):
    """
    Generate entire tree structure upfront.
    """
    nodes = []
    edges = []
    node_counter = 0

    ATOMS = ["C", "N", "O", "Br"]

    def build_subtree(parent_smiles, parent_id, level):
        nonlocal node_counter
        if level > depth:
            return

        num_children = random.choice([1, 2])

        for i in range(num_children):
            node_id = f"node_{node_counter}"
            node_counter += 1
            child_smiles = f"{parent_smiles}{ATOMS[i]}"

            node = {
                "id": node_id,
                "smiles": child_smiles,
                "label": f"Molecule-{level}{i}",
                "cost": random.uniform(10, 110),
                "energy": random.uniform(100, 600),
                "yield": random.uniform(0, 100),
                "level": level,
                "parentId": parent_id,
                "hoverInfo": f"# Molecule {level}-{i}\n**SMILES:** `{child_smiles}`\n**Level:** {level}",
            }
            nodes.append(node)

            edge = {
                "id": f"edge_{parent_id}_{node_id}",
                "from": parent_id,
                "to": node_id,
                "reactionType": random.choice(
                    ["Hydrogenation", "Oxidation", "Methylation", "Reduction", "Cyclization"]
                ),
            }
            edges.append(edge)

            build_subtree(child_smiles, node_id, level + 1)

    # Root node
    root_id = "root"
    root = {
        "id": root_id,
        "smiles": start_smiles,
        "label": "Root Molecule",
        "cost": random.uniform(10, 110),
        "energy": random.uniform(100, 600),
        "yield": 2.0,
        "level": 0,
        "parentId": None,
        "hoverInfo": f"# Root Molecule\n**SMILES:** `{start_smiles}`",
    }
    nodes.insert(0, root)

    build_subtree(start_smiles, root_id, 1)

    return nodes, edges


def calculate_positions(nodes: list[dict[str, Any]]):
    """
    Calculate positions for all nodes (matching frontend logic).
    """
    BOX_WIDTH = 220  # Must match with javascript!
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

        positioned_node = {**node, "x": 100 + node["level"] * level_gap, "y": 100 + index_in_level * node_spacing}
        positioned.append(positioned_node)

    return positioned


async def generate_molecules(start_smiles: str, depth: int = 3, websocket: WebSocket = None):
    """
    Stream positioned nodes and edges for the retrosynthesis sample.
    """

    # Generate and position entire tree upfront
    nodes, edges = generate_tree_structure(start_smiles, depth)
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


async def lead_molecule(start_smiles: str, depth: int = 3, websocket: WebSocket = None):
    """
    Stream positioned nodes and edges for the lead molecule optimization sample.
    """

    # Generate one node at a time
    for i in range(depth):
        if i > 0:
            edge_complete = {
                "type": "edge_update",
                "id": f"edge_{i-1}_{i}",
                "status": "complete",
                "label": "",
                "fromNode": {"id": f"node_{i-1}", "x": 0, "y": 0},
                "toNode": {"id": f"node_{i}", "x": 0, "y": 0},
            }
            await websocket.send_json(edge_complete)
        node = dict(
            id=f"node_{i}",
            smiles=start_smiles + "C" * i,
            label="Water",
            energy=i * random.uniform(0, 16),
            level=0,
            hoverInfo="This is some markdown\n# Hej",
            x=0,
            y=i * 150,
        )
        await websocket.send_json({"type": "node", **node})
        if i == depth - 1:
            break
        edge_data = {
            "type": "edge",
            "id": f"edge_{i}_{i+1}",
            "status": "computing",
            "label": "Optimizing",
            "fromNode": {"id": f"node_{i}", "x": 0, "y": 0},
            "toNode": {"id": f"node_{i+1}", "x": 0, "y": 0},
        }
        await websocket.send_json(edge_data)
        # TODO: Compute here
        await asyncio.sleep(0.8)

    await websocket.send_json({"type": "complete"})


@app.websocket("/ws")
async def websocket_endpoint(websocket: WebSocket):
    await websocket.accept()
    try:
        while True:
            data = await websocket.receive_json()

            if data["action"] == "compute":
                if data["problemType"] == "optimization":
                    await lead_molecule(data["smiles"], data.get("depth", 10), websocket=websocket)
                elif data["problemType"] == "retrosynthesis":
                    await generate_molecules(data["smiles"], data.get("depth", 3), websocket)
                else:
                    await websocket.send_json(
                        {"type": "error", "message": f"Unsupported problem type {data['problemType']}"}
                    )
            if data["action"] == "custom_query":
                if "water" in data['query'].lower():
                    await websocket.send_json(
                        {"type": "response", "message": f"Processing query: {data['query']} for node {data['nodeId']}", "smiles": "O"}
                    )
                else:
                    await websocket.send_json(
                        {"type": "response", "message": f"Processing query: {data['query']} for node {data['nodeId']}"}
                    )
                await asyncio.sleep(3)  # Random wait
                await websocket.send_json({"type": "complete"})
    except WebSocketDisconnect:
        pass


if __name__ == "__main__":
    import uvicorn

    uvicorn.run(app, host="127.0.0.1", port=8001)
