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
from fastapi.responses import StreamingResponse, HTMLResponse
from fastapi.staticfiles import StaticFiles
from fastapi.middleware.cors import CORSMiddleware
from dataclasses import dataclass, asdict
import asyncio
import copy
import os
import random
from typing import Any, Optional, Literal
import requests
from datetime import datetime
from loguru import logger

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

if os.path.exists(ASSETS_PATH):
    # Serve the frontend
    app.mount("/assets", StaticFiles(directory=ASSETS_PATH), name="assets")
    app.mount("/rdkit", StaticFiles(directory=os.path.join(DIST_PATH, "rdkit")), name="rdkit")

    @app.get("/")
    async def root(request: Request):
        logger.info(f"Request for Web UI received. Headers: {str(request.headers)}")
        with open(os.path.join(DIST_PATH, "index.html"), 'r') as fp:
            html = fp.read()

        html = html.replace('<!-- APP CONFIG -->', f'''
           <script>
           window.APP_CONFIG = {{
               WS_SERVER: '{os.getenv("WS_SERVER", "ws://localhost:8001/ws")}'
           }};
           </script>''')
        return HTMLResponse(html)

CACTUS = "https://cactus.nci.nih.gov/chemical/structure/{0}/{1}"

def smiles_to_iupac(smiles):
    return smiles
    try:
        rep = "iupac_name"
        url = CACTUS.format(smiles, rep)
        response = requests.get(url)
        response.raise_for_status()
        if response.text.startswith("<"):  # HTML
            return smiles
        return response.text
    except requests.exceptions.HTTPError:
        return smiles

@dataclass
class Node:
    id: str
    smiles: str
    label: str
    hoverInfo: str
    level: int
    parentId: Optional[str] = None
    x: Optional[int] = None
    y: Optional[int] = None
    # Properties
    cost: Optional[float] = None
    bandgap: Optional[float] = None
    density: Optional[float] = None
    yield_: Optional[float] = None
    highlight: Optional[str] = "normal"

    def json(self):
        ret = asdict(self)
        ret["yield"] = ret["yield_"]
        del ret["yield_"]
        return ret


@dataclass
class Edge:
    id: str
    fromNode: str
    toNode: str
    status: Literal["computing", "complete"]
    label: Optional[str] = None

    def json(self):
        return asdict(self)


@dataclass
class SidebarMessage:
    content: str
    smiles: str

    def json(self):
        return asdict(self)

@dataclass
class Tool:
    name: str
    description: Optional[str] = None

    def json(self):
        return asdict(self)


def generate_tree_structure(start_smiles: str, depth: int = 3):
    """
    Generate entire tree structure upfront.
    """

    nodes: list[Node] = []
    edges: list[Edge] = []
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

            node = Node(
                id=node_id,
                smiles=child_smiles,
                label=child_smiles,#smiles_to_iupac(child_smiles),
                cost=random.uniform(10, 110),
                bandgap=random.uniform(100, 600),
                yield_=random.uniform(0, 100),
                level=level,
                parentId=parent_id,
                hoverInfo=f"# Molecule {level}-{i}\n**SMILES:** `{child_smiles}`\n**Level:** {level}",
            )
            nodes.append(node)

            edge = Edge(
                id=f"edge_{parent_id}_{node_id}",
                fromNode=parent_id,
                toNode=node_id,
                status="computing",
                label=random.choice(["Hydrogenation", "Oxidation", "Methylation", "Reduction", "Cyclization"]),
            )
            edges.append(edge)

            build_subtree(child_smiles, node_id, level + 1)

    # Root node
    root_id = "root"
    root = Node(
        id=root_id,
        smiles=start_smiles,
        label=smiles_to_iupac(start_smiles),
        cost=random.uniform(10, 110),
        bandgap=random.uniform(100, 600),
        yield_=2.0,
        level=0,
        parentId=None,
        hoverInfo=f"# Root Molecule\n**SMILES:** `{start_smiles}`",
    )
    nodes.insert(0, root)

    build_subtree(start_smiles, root_id, 1)

    return nodes, edges


def calculate_positions(nodes: list[Node]):
    """
    Calculate positions for all nodes (matching frontend logic).
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
    positioned: list[Node] = []
    for node in nodes:
        level_nodes = levels[node.level]
        index_in_level = level_nodes.index(node)

        positioned_node = copy.deepcopy(node)
        positioned_node.x = 100 + node.level * level_gap
        positioned_node.y = 100 + index_in_level * node_spacing
        positioned.append(positioned_node)

    return positioned


async def generate_molecules(start_smiles: str, depth: int = 3, websocket: WebSocket = None):
    """
    Stream positioned nodes and edges for the retrosynthesis sample.
    """

    # Generate and position entire tree upfront
    nodes, edges = generate_tree_structure(start_smiles, depth)
    positioned_nodes = calculate_positions(nodes)

    # Stream root first
    root = positioned_nodes[0]
    await websocket.send_json({"type": "node", "node": root.json()})
    await asyncio.sleep(0.8)

    # Stream remaining nodes with edges
    for i in range(1, len(positioned_nodes)):
        node = positioned_nodes[i]

        # Find edge for this node
        edge = next((e for e in edges if e.toNode == node.id), None)

        if edge:
            # Send edge with computing status
            edge_data = {
                "type": "edge",
                "edge": edge.json(),
            }
            edge_data["edge"]["label"] = f"Computing: {edge.label}"
            edge_data["edge"]["toNode"] = node.id
            await websocket.send_json(edge_data)

            await asyncio.sleep(0.6)

            # Send node
            await websocket.send_json({"type": "node", "node": node.json()})

            # Update edge to complete
            edge_complete = {
                "type": "edge_update",
                "edge": {
                    "id": edge.id,
                    "status": "complete",
                    "label": edge.label,
                }
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
                "edge": {
                    "id": f"edge_{i-1}_{i}",
                    "status": "complete",
                    "label": "",
                    "fromNode": {"id": f"node_{i-1}", "x": 0, "y": 0},
                    "toNode": {"id": f"node_{i}", "x": 0, "y": 0},
                }
            }
            await websocket.send_json(edge_complete)
        node = dict(
            id=f"node_{i}",
            smiles=start_smiles + "C" * i,
            label="Water",
            bandgap=i * random.uniform(0, 16),
            level=0,
            hoverInfo="This is some markdown\n# Hej",
            x=0,
            y=i * 150,
        )
        await websocket.send_json({"type": "node", "node": node})
        if i == depth - 1:
            break
        edge_data = {
            "type": "edge",
            "edge": {
                "id": f"edge_{i}_{i+1}",
                "status": "computing",
                "label": "Optimizing",
                "fromNode": {"id": f"node_{i}", "x": 0, "y": 0},
                "toNode": {"id": f"node_{i+1}", "x": 0, "y": 0},
            }
        }
        await websocket.send_json(edge_data)
        # TODO: Compute here
        await asyncio.sleep(0.8)

    await websocket.send_json({"type": "complete"})


@app.websocket("/ws")
async def websocket_endpoint(websocket: WebSocket):
    logger.info(f"Request for websocket received. Headers: {str(websocket.headers)}")
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
                        {"type": "error", "error": f"Unsupported problem type {data['problemType']}"}
                    )
            elif data["action"] == "compute-reaction-from":
                await websocket.send_json(
                    {"type": "subtree_update", "node": {"id": data['nodeId'], "highlight": "yellow"}, "withNode": True}
                )
            elif data["action"] == "recompute-reaction" and "query" not in data:
                await websocket.send_json(
                        {"type": "response", "message": { "source": "Some custom source", "message": f"Hi from server", "smiles": "CCO"} }
                    )
                await websocket.send_json({"type": "complete"})
            elif data["action"] == "recompute-reaction":
                if "water" in data["query"].lower():
                    await websocket.send_json(
                        {
                            "type": "response",
                            "message": {
                                "source": "System",
                                "message": f"Processing query: {data['query']} for node {data['nodeId']}",
                                "smiles": "O",
                            }
                        }
                    )
                else:
                    await websocket.send_json(
                        {"type": "response", "message": {"message": f"Processing query: {data['query']} for node {data['nodeId']}"}}
                    )
                await asyncio.sleep(3)  # Random wait
                await websocket.send_json({"type": "complete"})
            elif data["action"] == "list-tools":
                tools = [Tool(f"tool_{i}", f"Does what tool {i} does") for i in range(1, 16)]
                tools.append(Tool("a_tool_with_no_desc"))
                await websocket.send_json({
                    "type": "available-tools-response",
                    "tools": [tool.json() for tool in tools],
                })
            elif data["action"] == "save-context":
                await websocket.send_json({
                    "type": "save-context-response",
                    "experimentContext": f"this is some sample context we are saving at {datetime.now():%Y-%m-%d %H:%M:%S}",
                })
            elif data["action"] == "load-context":
                print(f"LOADED CONTEXT: {data['experimentContext']}")
            else:
                print("WARN: Unhandled message:", data)
    except WebSocketDisconnect:
        pass


if __name__ == "__main__":
    import uvicorn

    uvicorn.run(app, host="127.0.0.1", port=8001)
