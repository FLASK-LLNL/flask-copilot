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
from fastapi.responses import JSONResponse, HTMLResponse
from fastapi.staticfiles import StaticFiles
from fastapi.middleware.cors import CORSMiddleware
from contextlib import asynccontextmanager
from dataclasses import dataclass, asdict
import asyncio
import copy
import os
import sys
import random
import uuid
from typing import Any, Optional, Literal
from dataclasses import dataclass, field
from datetime import datetime, timedelta
import requests
from loguru import logger
from charge_backend.moleculedb.molecule_naming import smiles_to_html, MolNameFormat

# Import database components
sys.path.insert(0, os.path.dirname(__file__))
from backend.database.engine import engine, Base
from backend.routers import sessions as sessionsrouter


# Session tracking for computation resume
@dataclass
class ComputationSession:
    """Track an active computation for potential resume"""
    session_id: str
    smiles: str
    problem_type: str
    depth: int
    websocket: WebSocket
    created_at: datetime = field(default_factory=datetime.now)
    is_complete: bool = False
    is_cancelled: bool = False
    sent_nodes: list = field(default_factory=list)
    sent_edges: list = field(default_factory=list)
    pending_nodes: list = field(default_factory=list)
    pending_edges: list = field(default_factory=list)
    current_index: int = 0

active_sessions: dict[str, ComputationSession] = {}
SESSION_TIMEOUT_HOURS = 24

def cleanup_old_sessions():
    """Remove sessions older than SESSION_TIMEOUT_HOURS"""
    cutoff = datetime.now() - timedelta(hours=SESSION_TIMEOUT_HOURS)
    expired = [sid for sid, sess in active_sessions.items() if sess.created_at < cutoff]
    for sid in expired:
        del active_sessions[sid]


@asynccontextmanager
async def lifespan(app: FastAPI):
    """Database initialization on startup"""
    async with engine.begin() as conn:
        await conn.run_sync(Base.metadata.create_all)
    yield


app = FastAPI(lifespan=lifespan)

# CORS for development
app.add_middleware(
    CORSMiddleware,
    allow_origins=["*"],
    allow_credentials=True,
    allow_methods=["*"],
    allow_headers=["*"],
)

# Include the sessions router for database persistence
app.include_router(sessionsrouter.router)

if "FLASK_APPDIR" in os.environ:
    DIST_PATH = os.environ["FLASK_APPDIR"]
else:
    DIST_PATH = os.path.join(os.path.dirname(__file__), "flask-app", "dist")
ASSETS_PATH = os.path.join(DIST_PATH, "assets")

if os.path.exists(ASSETS_PATH):
    # Serve the frontend
    app.mount("/assets", StaticFiles(directory=ASSETS_PATH), name="assets")
    app.mount(
        "/rdkit", StaticFiles(directory=os.path.join(DIST_PATH, "rdkit")), name="rdkit"
    )

    @app.get("/")
    async def root(request: Request):
        logger.info(f"Request for Web UI received. Headers: {str(request.headers)}")
        with open(os.path.join(DIST_PATH, "index.html"), "r") as fp:
            html = fp.read()

        html = html.replace(
            "<!-- APP CONFIG -->",
            f"""
           <script>
           window.APP_CONFIG = {{
               WS_SERVER: '{os.getenv("WS_SERVER", "ws://localhost:8001/ws")}',
               VERSION: '{os.getenv("SERVER_VERSION", "")}'
           }};
           </script>""",
        )
        return HTMLResponse(html)


@dataclass
class PathwayStep:
    smiles: list[str]
    label: list[str]


@dataclass
class ReactionAlternative:
    id: str
    name: str
    type: Literal["exact", "template", "ai"]
    status: Literal["active", "available", "computing"]
    pathway: list[PathwayStep]
    hoverInfo: str
    disabled: Optional[bool] = None
    disabledReason: Optional[str] = None


@dataclass
class Reaction:
    id: str
    hoverInfo: str
    highlight: str = "normal"
    label: Optional[str] = None
    alternatives: Optional[list[ReactionAlternative]] = None
    templatesSearched: bool = False


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
    reaction: Optional[Reaction] = None

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


def generate_tree_structure(
    start_smiles: str,
    depth: int = 3,
    molecule_name_format: MolNameFormat = "brand",
):
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
                label=smiles_to_html(child_smiles, molecule_name_format),
                cost=random.uniform(10, 110),
                bandgap=random.uniform(100, 600),
                yield_=random.uniform(0, 100),
                level=level,
                parentId=parent_id,
                hoverInfo=f"# Molecule {level}-{i}\n**SMILES:** `{child_smiles}`\n**Level:** {level}",
                reaction=(
                    None
                    if level == depth
                    else Reaction(
                        "reaction_0",
                        """# Reaction (AI-based)\nSome information here

| Column A | Column B |
| -------- | -------- |
| 1        | `def`    |
| 2        | `ghi`    |
| 3        | $a+b$    |
                        """,
                        label=random.choice(
                            [
                                "",
                                "Hydrogenation",
                                "Oxidation",
                                "Methylation",
                                "Reduction",
                                "Cyclization",
                            ]
                        ),
                        alternatives=[
                            ReactionAlternative(
                                "a0",
                                "Oxidation",
                                "template",
                                "active",
                                [PathwayStep(["O", "CCO"], ["", ""])],
                                "Info",
                            ),
                        ]
                        + [
                            ReactionAlternative(
                                f"a{2*i+1}",
                                "Something",
                                "template",
                                "available",
                                [
                                    PathwayStep(["CCO", "CCC"], ["", ""]),
                                    PathwayStep(["CCI"], [""]),
                                ]
                                * 5,
                                "Info",
                            )
                            for i in range(1, 10)
                        ],
                        templatesSearched=True,
                    )
                ),
            )
            nodes.append(node)

            edge = Edge(
                id=f"edge_{parent_id}_{node_id}",
                fromNode=parent_id,
                toNode=node_id,
                status="computing",
            )
            edges.append(edge)

            build_subtree(child_smiles, node_id, level + 1)

    # Root node
    root_id = "root"
    root = Node(
        id=root_id,
        smiles=start_smiles,
        label=smiles_to_html(start_smiles, molecule_name_format),
        cost=random.uniform(10, 110),
        bandgap=random.uniform(100, 600),
        yield_=2.0,
        level=0,
        parentId=None,
        hoverInfo=f"# Root Molecule\n**SMILES:** `{start_smiles}`",
        reaction=Reaction(
            "reaction_0",
            "# Reaction (AI-based)\nSome information here",
            label=random.choice(
                [
                    "Hydrogenation",
                    "Oxidation",
                    "Methylation",
                    "Reduction",
                    "Cyclization",
                ]
            ),
        ),
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
    TEXT_FACTOR = 8  # Must match with javascript!
    node_spacing = 150

    # Group by level
    levels: dict[int, list[Node]] = {}
    for node in nodes:
        level = node.level
        if level not in levels:
            levels[level] = []
        levels[level].append(node)

    # Compute level gaps based on reaction label length
    level_gaps = [BOX_WIDTH + BOX_GAP] * len(levels)
    for level, level_nodes in levels.items():
        for node in level_nodes:
            if node.reaction and node.reaction.label:
                level_gaps[level + 1] = max(
                    level_gaps[level + 1],
                    BOX_WIDTH + BOX_GAP + len(node.reaction.label) * TEXT_FACTOR,
                )

    # Position nodes
    positioned: list[Node] = []
    for node in nodes:
        level_nodes = levels[node.level]
        index_in_level = level_nodes.index(node)

        positioned_node = copy.deepcopy(node)
        positioned_node.x = 100 + node.level * level_gaps[node.level]
        positioned_node.y = 100 + index_in_level * node_spacing
        positioned.append(positioned_node)

    return positioned


async def generate_molecules(
    start_smiles: str,
    websocket: WebSocket,
    depth: int = 3,
    molecule_name_format: MolNameFormat = "brand",
):
    """
    Stream positioned nodes and edges for the retrosynthesis sample.
    """

    # Generate and position entire tree upfront
    nodes, edges = generate_tree_structure(start_smiles, depth, molecule_name_format)
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
                },
            }
            await websocket.send_json(edge_complete)

            await asyncio.sleep(0.2)

    await websocket.send_json({"type": "complete"})


async def lead_molecule(
    start_smiles: str,
    websocket: WebSocket,
    depth: int = 3,
    molecule_name_format: MolNameFormat = "brand",
):
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
                },
            }
            await websocket.send_json(edge_complete)
        node = dict(
            id=f"node_{i}",
            smiles=start_smiles + "C" * i,
            label=smiles_to_html(start_smiles + "C" * i, molecule_name_format),
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
            },
        }
        await websocket.send_json(edge_data)
        # TODO: Compute here
        await asyncio.sleep(0.8)

    await websocket.send_json({"type": "complete"})


@app.websocket("/ws")
async def websocket_endpoint(websocket: WebSocket):
    logger.info(f"Request for websocket received. Headers: {str(websocket.headers)}")

    molecule_format = "brand"
    username = "nobody"
    current_session_id = None
    
    if "x-forwarded-user" in websocket.headers:
        username = websocket.headers["x-forwarded-user"]
    
    # Cleanup old sessions periodically
    cleanup_old_sessions()

    await websocket.accept()
    try:
        while True:
            data = await websocket.receive_json()
            
            # Handle session resume request
            if data["action"] == "resume_session":
                session_id = data.get("sessionId")
                if session_id and session_id in active_sessions:
                    session = active_sessions[session_id]
                    if not session.is_complete and not session.is_cancelled:
                        current_session_id = session_id
                        session.websocket = websocket
                        logger.info(f"Resuming session {session_id}")
                        
                        # Send session_resumed message with current state
                        await websocket.send_json({
                            "type": "session_resumed",
                            "sessionId": session_id,
                            "sentNodes": len(session.sent_nodes),
                            "totalNodes": len(session.pending_nodes) if session.pending_nodes else 0,
                            "isComplete": session.is_complete
                        })
                        
                        # Resume the computation
                        if session.problem_type == "optimization":
                            await lead_molecule(
                                session.smiles,
                                session.depth,
                                websocket=websocket,
                                molecule_name_format=molecule_format
                            )
                        elif session.problem_type == "retrosynthesis":
                            await generate_molecules(
                                session.smiles,
                                session.depth,
                                websocket,
                                molecule_format
                            )
                    else:
                        # Session was already complete or cancelled
                        await websocket.send_json({
                            "type": "session_status",
                            "sessionId": session_id,
                            "status": "complete" if session.is_complete else "cancelled"
                        })
                else:
                    # Session not found or expired
                    await websocket.send_json({
                        "type": "session_not_found",
                        "sessionId": session_id
                    })
                continue

            if "runSettings" in data:
                molecule_format = data["runSettings"]["moleculeName"]

            if data["action"] == "compute":
                # Create new session
                session_id = data.get("sessionId") or str(uuid.uuid4())
                depth = data.get("depth", 10 if data["problemType"] == "optimization" else 3)
                
                session = ComputationSession(
                    session_id=session_id,
                    smiles=data["smiles"],
                    problem_type=data["problemType"],
                    depth=depth,
                    websocket=websocket
                )
                active_sessions[session_id] = session
                current_session_id = session_id
                
                # Send session ID to client
                await websocket.send_json({
                    "type": "session_started",
                    "sessionId": session_id
                })
                
                logger.info(f"Started new session {session_id} for {data['problemType']}")
                
                if data["problemType"] == "optimization":
                    await lead_molecule(
                        data["smiles"],
                        depth=data.get("depth", 10),
                        websocket=websocket,
                        molecule_name_format=molecule_format,
                    )
                elif data["problemType"] == "retrosynthesis":
                    await generate_molecules(
                        data["smiles"], websocket, data.get("depth", 3), molecule_format
                    )
                else:
                    await websocket.send_json(
                        {
                            "type": "error",
                            "error": f"Unsupported problem type {data['problemType']}",
                        }
                    )
            elif data["action"] == "compute-reaction-from":
                await websocket.send_json(
                    {
                        "type": "subtree_update",
                        "node": {"id": data["nodeId"], "highlight": "yellow"},
                        "withNode": True,
                    }
                )
            elif data["action"].startswith("query-"):
                message = f"Processing query: {data['query']} for {'molecule' if data['action'].endswith('molecule') else 'reaction'} {data['nodeId']}"
                if "water" in data["query"].lower():
                    await websocket.send_json(
                        {
                            "type": "response",
                            "message": {
                                "source": "System",
                                "message": message,
                                "smiles": "O",
                            },
                        }
                    )
                else:
                    await websocket.send_json(
                        {
                            "type": "response",
                            "message": {"message": message},
                        }
                    )
                await websocket.send_json({"type": "complete"})
            elif data["action"] == "list-tools":
                tools = [
                    Tool(f"tool_{i}", f"Does what tool {i} does") for i in range(1, 16)
                ]
                tools.append(Tool("a_tool_with_no_desc"))
                await websocket.send_json(
                    {
                        "type": "available-tools-response",
                        "tools": [tool.json() for tool in tools],
                    }
                )
            elif data["action"] == "save-context":
                await websocket.send_json(
                    {
                        "type": "save-context-response",
                        "experimentContext": f"this is some sample context we are saving at {datetime.now():%Y-%m-%d %H:%M:%S}",
                    }
                )
            elif data["action"] == "load-context":
                print(f"LOADED CONTEXT: {data['experimentContext']}")
            elif data["action"] == "get-username":
                await websocket.send_json(
                    {
                        "type": "get-username-response",
                        "username": username,
                    }
                )
            elif data["action"] == "ui-update-orchestrator-settings":
                molecule_format = data["moleculeName"]
            elif data["action"] == "stop":
                if current_session_id and current_session_id in active_sessions:
                    active_sessions[current_session_id].is_cancelled = True
                    logger.info(f"Session {current_session_id} stopped by user")
            elif data["action"] == "reset":
                if current_session_id and current_session_id in active_sessions:
                    active_sessions[current_session_id].is_cancelled = True
                    del active_sessions[current_session_id]
                current_session_id = None
                logger.info("Session reset by user")
            else:
                print("WARN: Unhandled message:", data)
    except WebSocketDisconnect:
        pass


if __name__ == "__main__":
    import uvicorn

    uvicorn.run(app, host="127.0.0.1", port=8001)
