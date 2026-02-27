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
from db_backend.database.engine import engine, Base
from db_backend.routers import sessions as sessionsrouter
from db_backend.routers import projects as projectsrouter


# Session tracking for computation resume
@dataclass
class ComputationSession:
    """Track an active computation for potential resume.

    The session holds an *optional* websocket reference.  When the browser
    disconnects the websocket is set to ``None`` (``ws_connected = False``)
    but the background computation task keeps running, accumulating nodes
    and edges.  A reconnecting browser can re-attach via ``resume_session``.
    """
    session_id: str
    smiles: str
    problem_type: str
    depth: int
    username: str = "nobody"
    experiment_id: str | None = None  # DB experiment row for direct updates
    websocket: WebSocket | None = None
    ws_connected: bool = True
    created_at: datetime = field(default_factory=datetime.now)
    is_complete: bool = False
    is_cancelled: bool = False
    sent_nodes: list = field(default_factory=list)
    sent_edges: list = field(default_factory=list)
    sent_messages: list = field(default_factory=list)  # reasoning/sidebar messages
    pending_nodes: list = field(default_factory=list)
    pending_edges: list = field(default_factory=list)
    current_index: int = 0
    task: Any = field(default=None, repr=False)  # asyncio.Task running the computation
    # How often (in nodes) to persist progress to DB while running headless
    _db_save_interval: int = 3
    _send_lock: asyncio.Lock = field(default_factory=asyncio.Lock, repr=False)

    async def safe_send(self, payload: dict) -> None:
        """Send *payload* over the WebSocket if still connected.

        If the send fails (client disconnected), the session is marked
        ``ws_connected = False`` so the computation can continue headless.
        Uses a lock to prevent interleaved sends from background task and
        WS read loop (e.g. during replay).

        Every outgoing message is tagged with ``experimentId`` (when known)
        so the frontend can route messages to the correct experiment and
        ignore messages arriving for background experiments.

        Reasoning messages (type=response) are ALWAYS tracked regardless
        of ws_connected, so they survive browser disconnects and get
        persisted to sidebar_state by save_session_to_db.
        """
        # Track reasoning messages unconditionally (even when headless)
        if payload.get("type") == "response" and payload.get("message"):
            import time
            msg = payload["message"]
            if isinstance(msg, dict):
                tracked_msg = {
                    "id": int(time.time() * 1000) + len(self.sent_messages),
                    "timestamp": datetime.now().isoformat(),
                    **msg,
                }
                if "smiles" not in tracked_msg:
                    tracked_msg["smiles"] = None
                self.sent_messages.append(tracked_msg)

        if not self.ws_connected or self.websocket is None:
            return
        async with self._send_lock:
            if not self.ws_connected or self.websocket is None:
                return
            try:
                # Tag every message with the experiment it belongs to
                if self.experiment_id and "experimentId" not in payload:
                    payload = {**payload, "experimentId": self.experiment_id}
                await self.websocket.send_json(payload)
            except (WebSocketDisconnect, RuntimeError, Exception) as exc:
                logger.info(
                    f"Session {self.session_id}: WebSocket send failed ({type(exc).__name__}), "
                    f"continuing headless"
                )
                self.ws_connected = False
                self.websocket = None

    def detach_ws(self) -> None:
        """Detach the websocket so the computation continues headless."""
        self.ws_connected = False
        self.websocket = None

    def attach_ws(self, ws: WebSocket) -> None:
        """Attach a new websocket (e.g. from a reconnecting browser)."""
        self.websocket = ws
        self.ws_connected = True


active_sessions: dict[str, ComputationSession] = {}
SESSION_TIMEOUT_HOURS = 24

def cleanup_old_sessions():
    """Remove sessions older than SESSION_TIMEOUT_HOURS"""
    cutoff = datetime.now() - timedelta(hours=SESSION_TIMEOUT_HOURS)
    expired = [sid for sid, sess in active_sessions.items() if sess.created_at < cutoff]
    for sid in expired:
        del active_sessions[sid]


async def save_session_to_db(session: ComputationSession) -> None:
    """Persist the tracked nodes/edges of a ComputationSession to MariaDB.

    Called on WebSocket disconnect and on computation completion so that
    a second browser instance can load the latest state from the database.

    The function locates the experiment row whose ``graph_state``
    ``serverSessionId`` matches ``session.session_id`` (the UUID assigned
    by the WebSocket handler).  If no matching row is found it falls back
    to the most-recently-modified running experiment for the user.
    """
    from db_backend.database.engine import AsyncSessionLocal
    from db_backend.database import models as db_models
    from sqlalchemy import select, and_, cast, String, func

    if AsyncSessionLocal is None:
        logger.warning("save_session_to_db: no database engine available")
        return

    if not session.sent_nodes:
        # Even with no nodes, persist problem_type and smiles so the
        # experiment record is correct if the user reloads before any
        # results arrive.
        if session.experiment_id and (session.problem_type or session.smiles):
            try:
                async with AsyncSessionLocal() as db:
                    result = await db.execute(
                        select(db_models.Experiment).where(
                            db_models.Experiment.id == session.experiment_id
                        )
                    )
                    experiment = result.scalar_one_or_none()
                    if experiment:
                        if session.problem_type:
                            experiment.problem_type = session.problem_type
                        if session.smiles and not experiment.smiles:
                            experiment.smiles = session.smiles
                        await db.commit()
                        logger.debug(
                            f"save_session_to_db: no nodes yet, saved "
                            f"problem_type/smiles for {session.session_id}"
                        )
            except Exception as exc:
                logger.error(
                    f"save_session_to_db: problem_type-only save failed "
                    f"for {session.session_id}: {exc}"
                )
        return

    try:
        async with AsyncSessionLocal() as db:
            experiment = None

            # Strategy 0: direct lookup by experiment_id (set from frontend)
            if session.experiment_id:
                result = await db.execute(
                    select(db_models.Experiment).where(
                        db_models.Experiment.id == session.experiment_id
                    )
                )
                experiment = result.scalar_one_or_none()

            # Strategy 1: find experiment by serverSessionId in graph_state JSON
            if experiment is None:
                try:
                    result = await db.execute(
                        select(db_models.Experiment).where(
                            func.json_unquote(
                                func.json_extract(
                                    db_models.Experiment.graph_state,
                                    "$.serverSessionId",
                                )
                            )
                            == session.session_id
                        ).order_by(db_models.Experiment.last_modified.desc()).limit(1)
                    )
                    experiment = result.scalar_one_or_none()
                except Exception as e:
                    logger.debug(f"save_session_to_db: JSON query failed ({e}), trying fallback")

            # Strategy 2: fall back to latest running experiment with matching smiles
            if experiment is None:
                result = await db.execute(
                    select(db_models.Experiment).where(
                        and_(
                            db_models.Experiment.is_running == True,
                            db_models.Experiment.smiles == session.smiles,
                        )
                    ).order_by(db_models.Experiment.last_modified.desc()).limit(1)
                )
                experiment = result.scalar_one_or_none()

            # Strategy 3: create a new experiment row with a default project.
            # This handles the edge case where the browser disconnected before
            # the frontend auto-saved (so no experiment row exists yet).
            if experiment is None:
                logger.info(
                    f"save_session_to_db: no existing experiment row for session "
                    f"{session.session_id}, creating one"
                )
                # Ensure a default project exists for this user
                default_project_id = f"default-{session.username}"
                result = await db.execute(
                    select(db_models.Project).where(
                        db_models.Project.id == default_project_id
                    )
                )
                project = result.scalar_one_or_none()
                if project is None:
                    project = db_models.Project(
                        id=default_project_id,
                        user=session.username,
                        name="Default Project",
                    )
                    db.add(project)
                    await db.flush()

                experiment = db_models.Experiment(
                    id=str(uuid.uuid4()),
                    project_id=default_project_id,
                    user=session.username,
                    name=f"Auto-saved ({session.problem_type})",
                    smiles=session.smiles,
                    problem_type=session.problem_type,
                    graph_state={"serverSessionId": session.session_id},
                )
                db.add(experiment)

            experiment.tree_nodes = session.sent_nodes
            experiment.edges = session.sent_edges
            experiment.is_running = not session.is_complete
            experiment.last_modified = datetime.now(tz=None)

            # Also persist SMILES and problem_type so the experiment
            # record is complete even when the frontend never saved.
            if session.smiles and not experiment.smiles:
                experiment.smiles = session.smiles
            # Always update problem_type from the session — the server
            # knows the authoritative computation type from the original
            # 'compute' request.
            if session.problem_type:
                experiment.problem_type = session.problem_type

            # Persist reasoning/sidebar messages so they survive browser
            # disconnects and are visible when the user reopens the
            # experiment later.  The server tracks every `response`
            # message, so its list is authoritative.
            if session.sent_messages:
                experiment.sidebar_state = {
                    "messages": session.sent_messages,
                    "sourceFilterOpen": False,
                    "visibleSources": {
                        s
                        for m in session.sent_messages
                        if (s := m.get("source"))
                    },
                }
                # Convert visible sources set to dict (frontend expects {source: bool})
                experiment.sidebar_state["visibleSources"] = {
                    src: True
                    for src in experiment.sidebar_state["visibleSources"]
                }

            # Ensure the serverSessionId is always stored in graph_state
            # so that a reconnecting browser can resume the computation.
            # Also mirror sidebarMessages into graph_state so that both
            # the session load path (_experiment_to_state) and the project
            # load path (experiment_to_dict) can find the data.
            gs = experiment.graph_state or {}
            gs_dirty = False
            if gs.get("serverSessionId") != session.session_id:
                gs["serverSessionId"] = session.session_id
                gs_dirty = True
            if session.sent_messages:
                gs["sidebarMessages"] = session.sent_messages
                gs_dirty = True
            if gs_dirty:
                experiment.graph_state = gs

            await db.commit()
            logger.info(
                f"save_session_to_db: saved {len(session.sent_nodes)} nodes / "
                f"{len(session.sent_edges)} edges for session {session.session_id} "
                f"(experiment={experiment.id}, complete={session.is_complete})"
            )
    except Exception as exc:
        logger.error(f"save_session_to_db failed for {session.session_id}: {exc}")


@asynccontextmanager
async def lifespan(app: FastAPI):
    """Database initialization on startup"""
    if engine is not None:
        try:
            async with engine.begin() as conn:
                await conn.run_sync(Base.metadata.create_all)
            logger.info("Database tables initialized successfully")
        except Exception as e:
            logger.warning(f"Database initialization failed (app will run without persistence): {e}")
    else:
        logger.warning("No database engine available. Running without persistence.")
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
# Include the projects router for multi-browser experiment sharing
app.include_router(projectsrouter.router)

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
    depth: int = 3,
    molecule_name_format: MolNameFormat = "brand",
    session: ComputationSession | None = None,
):
    """
    Generate and stream retrosynthesis molecules.

    All WebSocket sends go through ``session.safe_send`` so that a
    browser disconnect does **not** kill this coroutine.  The computation
    continues headless and persists results to MariaDB.
    """
    if session is None:
        logger.error("generate_molecules called without a session")
        return

    # Generate and position entire tree upfront
    nodes, edges = generate_tree_structure(start_smiles, depth, molecule_name_format)
    positioned_nodes = calculate_positions(nodes)

    # Send initial reasoning message
    await session.safe_send({
        "type": "response",
        "message": {
            "source": "Reasoning",
            "message": f"**Starting Retrosynthesis Analysis**\n\nTarget molecule: `{start_smiles}`\n\nExploring synthetic pathways with depth {depth}...",
        }
    })

    # Stream root first
    root = positioned_nodes[0]
    await session.safe_send({"type": "node", "node": root.json()})
    session.sent_nodes.append(root.json())

    await session.safe_send({
        "type": "response",
        "message": {
            "source": "Reasoning",
            "message": "Identified target molecule. Now analyzing possible disconnection strategies...",
        }
    })
    await asyncio.sleep(0.8)

    # Retrosynthesis reasoning templates
    retro_reasoning = [
        "Analyzing bond disconnection possibilities in the target structure...",
        "Identified potential synthetic precursors through retrosynthetic analysis.",
        "Evaluating reaction feasibility and yield predictions for this transformation.",
        "Cross-referencing with known chemical reactions in the database.",
        "This disconnection follows established synthetic methodology.",
    ]

    # Stream remaining nodes with edges
    for i in range(1, len(positioned_nodes)):
        # Check for cancellation
        if session.is_cancelled:
            logger.info(f"Session {session.session_id}: cancelled, stopping generation")
            break

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
            await session.safe_send(edge_data)
            # Track the computing edge; will be updated to 'complete' below.
            session.sent_edges.append({
                "id": edge.id,
                "fromNode": edge.fromNode,
                "toNode": edge.toNode,
                "status": "computing",
                "label": f"Computing: {edge.label}",
            })

            # Send reasoning about this step
            reasoning_msg = retro_reasoning[i % len(retro_reasoning)]
            await session.safe_send({
                "type": "response",
                "message": {
                    "source": "Reasoning",
                    "message": f"**Step {i}**: {reasoning_msg}",
                }
            })

            await asyncio.sleep(0.6)

            # Send node
            await session.safe_send({"type": "node", "node": node.json()})
            session.sent_nodes.append(node.json())

            # Send info about the discovered precursor
            await session.safe_send({
                "type": "response",
                "message": {
                    "source": "Logger (Info)",
                    "message": f"Found precursor: `{node.smiles}`",
                    "smiles": node.smiles,
                }
            })

            # Update edge to complete
            edge_complete = {
                "type": "edge_update",
                "edge": {
                    "id": edge.id,
                    "status": "complete",
                    "label": edge.label,
                },
            }
            await session.safe_send(edge_complete)
            # Update tracked edge status from computing -> complete
            for tracked in session.sent_edges:
                if tracked["id"] == edge.id:
                    tracked["status"] = "complete"
                    tracked["label"] = edge.label
                    break

            # Periodic DB save so headless progress is visible to other browsers
            if not session.ws_connected and len(session.sent_nodes) % session._db_save_interval == 0:
                await save_session_to_db(session)

            await asyncio.sleep(0.2)

    # Send completion reasoning (tracked in sent_messages by safe_send)
    await session.safe_send({
        "type": "response",
        "message": {
            "source": "Reasoning",
            "message": f"**Retrosynthesis Complete**\n\nSuccessfully identified {len(positioned_nodes)} molecules in the synthetic pathway.",
        }
    })

    # Persist to DB BEFORE sending the 'complete' WebSocket message.
    # This ensures that when the frontend receives 'complete' and
    # refreshes project data, the sidebar_state (including this
    # completion message) is already in the database.
    session.is_complete = True
    await save_session_to_db(session)

    await session.safe_send({"type": "complete"})
    logger.info(
        f"Session {session.session_id}: computation finished "
        f"({len(session.sent_nodes)} nodes, ws_connected={session.ws_connected})"
    )


async def lead_molecule(
    start_smiles: str,
    depth: int = 3,
    molecule_name_format: MolNameFormat = "brand",
    session: ComputationSession | None = None,
):
    """
    Generate and stream lead molecule optimization results.

    Like ``generate_molecules``, all sends go through ``session.safe_send``
    so a browser disconnect does not kill this coroutine.
    """
    if session is None:
        logger.error("lead_molecule called without a session")
        return

    # Sample reasoning messages to simulate LLM thinking
    reasoning_templates = [
        "Analyzing molecular structure and identifying potential modification sites...",
        "Evaluating functional groups for optimization. The current molecule shows promising structural features.",
        "Applying chemical intuition to extend the carbon chain. This modification should improve the target property.",
        "Cross-referencing with known structure-activity relationships in the chemical space.",
        "Validating the proposed modification against chemical feasibility constraints.",
        "Computing predicted properties for the new molecular candidate...",
        "The optimization step successfully generated a valid molecular structure with improved characteristics.",
    ]

    # Send initial reasoning message
    await session.safe_send({
        "type": "response",
        "message": {
            "source": "Reasoning",
            "message": f"**Starting Lead Molecule Optimization**\n\nInitial molecule: `{start_smiles}`\n\nPlanned optimization depth: {depth} iterations",
        }
    })
    await asyncio.sleep(0.3)

    # Generate one node at a time
    for i in range(depth):
        # Check for cancellation
        if session.is_cancelled:
            logger.info(f"Session {session.session_id}: cancelled, stopping optimization")
            break

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
            await session.safe_send(edge_complete)
            # Update tracked edge status to complete
            for tracked in session.sent_edges:
                if tracked["id"] == f"edge_{i-1}_{i}":
                    tracked["status"] = "complete"
                    break
        
        # Send reasoning message for this iteration
        reasoning_idx = i % len(reasoning_templates)
        await session.safe_send({
            "type": "response",
            "message": {
                "source": "Reasoning",
                "message": f"**Iteration {i + 1}/{depth}**\n\n{reasoning_templates[reasoning_idx]}",
            }
        })
        
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
        await session.safe_send({"type": "node", "node": node})
        session.sent_nodes.append(node)
        
        # Send a Logger (Info) message about the generated molecule
        await session.safe_send({
            "type": "response",
            "message": {
                "source": "Logger (Info)",
                "message": f"Generated molecule: `{start_smiles + 'C' * i}`",
                "smiles": start_smiles + "C" * i,
            }
        })
        
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
        await session.safe_send(edge_data)
        # Record edge immediately; its status will be updated to complete
        # on the next iteration's edge_update message.
        session.sent_edges.append({
            "id": f"edge_{i}_{i+1}",
            "fromNode": f"node_{i}",
            "toNode": f"node_{i+1}",
            "status": "computing",
            "label": "Optimizing",
        })

        # Periodic DB save so headless progress is visible to other browsers
        if not session.ws_connected and len(session.sent_nodes) % session._db_save_interval == 0:
            await save_session_to_db(session)

        await asyncio.sleep(0.8)

    # Send completion reasoning (tracked in sent_messages by safe_send)
    await session.safe_send({
        "type": "response",
        "message": {
            "source": "Reasoning",
            "message": f"**Optimization Complete**\n\nSuccessfully generated {depth} molecules in the optimization pathway.",
        }
    })

    # Persist to DB BEFORE sending the 'complete' WebSocket message.
    # This ensures that when the frontend receives 'complete' and
    # refreshes project data, the sidebar_state (including this
    # completion message) is already in the database.
    session.is_complete = True
    await save_session_to_db(session)

    await session.safe_send({"type": "complete"})
    logger.info(
        f"Session {session.session_id}: optimization finished "
        f"({len(session.sent_nodes)} nodes, ws_connected={session.ws_connected})"
    )


def _add_task_done_callback(task: asyncio.Task, session: ComputationSession) -> None:
    """Attach a done-callback that logs unhandled exceptions from the
    background computation task so they are not silently swallowed."""
    def _on_done(t: asyncio.Task) -> None:
        if t.cancelled():
            logger.info(f"Session {session.session_id}: task was cancelled")
        elif t.exception():
            logger.error(
                f"Session {session.session_id}: task raised {t.exception()!r}"
            )
    task.add_done_callback(_on_done)


@app.websocket("/ws")
async def websocket_endpoint(websocket: WebSocket):
    logger.info(f"Request for websocket received. Headers: {str(websocket.headers)}")

    molecule_format = "brand"
    username = "nobody"
    current_session_id = None
    # Track ALL sessions spawned/resumed by this WS connection so we can
    # persist every one of them when the browser disconnects.
    ws_session_ids: set[str] = set()
    
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
                    # Detach previous session if switching to a different one
                    if current_session_id and current_session_id != session_id and current_session_id in active_sessions:
                        active_sessions[current_session_id].detach_ws()
                        logger.info(f"Detached session {current_session_id} (switching to resumed session {session_id})")

                    session = active_sessions[session_id]
                    current_session_id = session_id
                    ws_session_ids.add(session_id)

                    if session.is_cancelled:
                        await websocket.send_json({
                            "type": "session_status",
                            "sessionId": session_id,
                            "experimentId": session.experiment_id,
                            "status": "cancelled",
                        })
                    elif session.is_complete:
                        # Computation already finished – replay the full
                        # result set so the new browser renders everything.
                        session.attach_ws(websocket)
                        exp_id = session.experiment_id
                        async with session._send_lock:
                            for node in list(session.sent_nodes):
                                await websocket.send_json({"type": "node", "node": node, "experimentId": exp_id})
                            for edge in list(session.sent_edges):
                                await websocket.send_json({"type": "edge", "edge": edge, "experimentId": exp_id})
                            # Replay reasoning/sidebar messages so the
                            # frontend populates the reasoning panel
                            # (including "Retrosynthesis Complete").
                            for msg in list(session.sent_messages):
                                await websocket.send_json({"type": "response", "message": msg, "experimentId": exp_id})
                            await websocket.send_json({"type": "complete", "experimentId": exp_id})
                        logger.info(
                            f"Replayed completed session {session_id} "
                            f"({len(session.sent_nodes)} nodes)"
                        )
                    elif session.task and not session.task.done():
                        # Computation still running in the background.
                        # Re-attach the websocket so new results stream to
                        # this browser, and replay everything accumulated
                        # so far (the frontend deduplicates by id).
                        session.attach_ws(websocket)
                        exp_id = session.experiment_id
                        async with session._send_lock:
                            for node in list(session.sent_nodes):
                                await websocket.send_json({"type": "node", "node": node, "experimentId": exp_id})
                            for edge in list(session.sent_edges):
                                await websocket.send_json({"type": "edge", "edge": edge, "experimentId": exp_id})
                        await websocket.send_json({
                            "type": "session_resumed",
                            "sessionId": session_id,
                            "experimentId": session.experiment_id,
                            "sentNodes": len(session.sent_nodes),
                            "isComplete": False,
                        })
                        logger.info(
                            f"Reconnected browser to running session {session_id} "
                            f"(replayed {len(session.sent_nodes)} nodes, "
                            f"task still running)"
                        )
                    else:
                        # Task reference missing or already finished
                        # without marking complete (edge case).
                        await websocket.send_json({
                            "type": "session_status",
                            "sessionId": session_id,
                            "status": "unknown",
                        })
                else:
                    # Session not found or expired
                    await websocket.send_json({
                        "type": "session_not_found",
                        "sessionId": session_id,
                    })
                continue

            if "runSettings" in data:
                molecule_format = data["runSettings"]["moleculeName"]

            if data["action"] == "compute":
                # Detach previous session so it continues headless
                # while the new computation streams to the frontend.
                if current_session_id and current_session_id in active_sessions:
                    old_session = active_sessions[current_session_id]
                    old_session.detach_ws()
                    logger.info(f"Detached session {current_session_id} (starting new computation)")

                # Create new session
                session_id = data.get("sessionId") or str(uuid.uuid4())
                depth = data.get("depth", 10 if data["problemType"] == "optimization" else 3)
                
                session = ComputationSession(
                    session_id=session_id,
                    smiles=data["smiles"],
                    problem_type=data["problemType"],
                    depth=depth,
                    username=username,
                    experiment_id=data.get("experimentId"),
                    websocket=websocket,
                )
                active_sessions[session_id] = session
                current_session_id = session_id
                ws_session_ids.add(session_id)
                
                # Send session ID to client
                await websocket.send_json({
                    "type": "session_started",
                    "sessionId": session_id,
                    "experimentId": data.get("experimentId"),
                })
                
                logger.info(f"Started new session {session_id} for {data['problemType']}")
                
                # Immediately persist problem_type to the experiment DB
                # row so it is available even if the user reloads before
                # any nodes are generated.
                if session.experiment_id and session.problem_type:
                    try:
                        from db_backend.database.engine import AsyncSessionLocal
                        from db_backend.database import models as db_models
                        from sqlalchemy import select as sa_select
                        if AsyncSessionLocal is not None:
                            async with AsyncSessionLocal() as _db:
                                _result = await _db.execute(
                                    sa_select(db_models.Experiment).where(
                                        db_models.Experiment.id == session.experiment_id
                                    )
                                )
                                _exp = _result.scalar_one_or_none()
                                if _exp:
                                    _exp.problem_type = session.problem_type
                                    if session.smiles and not _exp.smiles:
                                        _exp.smiles = session.smiles
                                    await _db.commit()
                    except Exception as exc:
                        logger.debug(f"Early problem_type persist failed: {exc}")
                
                # Launch computation as a *background task* so the WS
                # read loop keeps processing stop / reset / query messages.
                # If the browser disconnects, the task continues headless.
                if data["problemType"] == "optimization":
                    session.task = asyncio.create_task(
                        lead_molecule(
                            data["smiles"],
                            depth=depth,
                            molecule_name_format=molecule_format,
                            session=session,
                        )
                    )
                    _add_task_done_callback(session.task, session)
                elif data["problemType"] == "retrosynthesis":
                    session.task = asyncio.create_task(
                        generate_molecules(
                            data["smiles"],
                            depth,
                            molecule_format,
                            session=session,
                        )
                    )
                    _add_task_done_callback(session.task, session)
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
                # Stop supports targeting a specific session by sessionId,
                # falling back to current_session_id for backwards compat.
                target = data.get("sessionId") or current_session_id
                if target and target in active_sessions:
                    active_sessions[target].is_cancelled = True
                    logger.info(f"Session {target} stopped by user")
            elif data["action"] == "detach":
                # Client is switching to a non-computing experiment.
                # Detach the current session so it continues headless.
                if current_session_id and current_session_id in active_sessions:
                    active_sessions[current_session_id].detach_ws()
                    logger.info(f"Detached session {current_session_id} by client request")
                current_session_id = None
            elif data["action"] == "reset":
                if current_session_id and current_session_id in active_sessions:
                    active_sessions[current_session_id].is_cancelled = True
                    del active_sessions[current_session_id]
                current_session_id = None
                logger.info("Session reset by user")
            else:
                print("WARN: Unhandled message:", data)
    except WebSocketDisconnect:
        # Browser disconnected.  Detach ALL sessions owned by this WS
        # connection so background computation keeps running headless,
        # and persist every one of them to the DB immediately so a new
        # browser can load the latest state.
        saved_count = 0
        for sid in ws_session_ids:
            sess = active_sessions.get(sid)
            if sess is None:
                continue
            task_running = sess.task is not None and not sess.task.done()
            logger.info(
                f"WebSocket disconnected: saving session {sid} "
                f"({len(sess.sent_nodes)} nodes, task_running={task_running})"
            )
            sess.detach_ws()
            try:
                await save_session_to_db(sess)
                saved_count += 1
            except Exception as exc:
                logger.error(f"Failed to save session {sid} on disconnect: {exc}")
        if saved_count == 0:
            logger.info("WebSocket disconnected (no active sessions to save)")


if __name__ == "__main__":
    import uvicorn
    import logging

    # Filter out the noisy /api/sessions/save access log lines
    class SuppressSessionSave(logging.Filter):
        def filter(self, record: logging.LogRecord) -> bool:
            msg = record.getMessage()
            if "/api/sessions/save" in msg:
                return False
            if "/api/projects/" in msg and "GET" in msg:
                return False
            return True

    logging.getLogger("uvicorn.access").addFilter(SuppressSessionSave())

    uvicorn.run(app, host="127.0.0.1", port=8001)
