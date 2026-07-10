################################################################################
## Copyright 2025-2026 Lawrence Livermore National Security, LLC..
## See the top-level LICENSE file for details.
##
## SPDX-License-Identifier: Apache-2.0
################################################################################

from functools import partial
from fastapi import FastAPI, WebSocket, Request
from fastapi.responses import HTMLResponse
from fastapi.staticfiles import StaticFiles
from fastapi.middleware.cors import CORSMiddleware
import os
import argparse
import json
from lc_conductor import try_get_public_hostname
from lc_conductor.session import SessionTimedOut, UserSessionManager
from charge_backend.flask_session import FlaskUserSession
import os
from typing import cast

from loguru import logger
from charge.clients.client import Client


from lc_conductor.tool_registration import (
    register_url,
    register_post,
    reload_server_list,
    check_mcp_servers_endpoint,
    validate_mcp_server_endpoint,
    delete_mcp_server_endpoint,
    get_registered_servers,
)
from lc_conductor import (
    discover_models_endpoint,
    validate_initial_model,
    resolve_orchestrator_config,
)

parser = argparse.ArgumentParser()

parser.add_argument(
    "--json_file",
    type=str,
    default="known_molecules.json",
    help="Path to the JSON file containing known molecules.",
)

parser.add_argument(
    "--max_retries",
    type=int,
    default=5,
)

parser.add_argument(
    "--config-file",
    type=str,
    default="/data/config.yml",
    help="Path to the configuration file for AiZynthFinder.",
)
parser.add_argument("--port", type=int, default=8001, help="Port to run the server on")
parser.add_argument("--host", type=str, default=None, help="Host to run the server on")

parser.add_argument(
    "--tool-server-cache",
    type=str,
    default="flask_copilot_active_tool_servers.json",
    help="Path to the JSON file containing current list of active MCP tool servers.",
)

# Add standard CLI arguments
Client.add_std_parser_arguments(parser, defaults=dict(backend="livai", model="gpt-5.4"))

args, _ = parser.parse_known_args()

app = FastAPI()

# CORS for development
app.add_middleware(
    CORSMiddleware,
    allow_origins=["*"],
    allow_credentials=True,
    allow_methods=["*"],
    allow_headers=["*"],
)

if "FLASK_APPDIR" in os.environ:
    DIST_PATH = os.environ["FLASK_APPDIR"]
else:
    DIST_PATH = os.path.join(os.path.dirname(__file__), "flask-app", "dist")
ASSETS_PATH = os.path.join(DIST_PATH, "assets")

reload_server_list(args.tool_server_cache)

app.post("/register")(partial(register_post, args.tool_server_cache))

app.post("/validate-mcp-server")(
    partial(validate_mcp_server_endpoint, args.tool_server_cache)
)

app.post("/delete-mcp-server")(
    partial(delete_mcp_server_endpoint, args.tool_server_cache)
)


app.post("/check-mcp-servers")(check_mcp_servers_endpoint)

app.post("/api/discover-models")(discover_models_endpoint)


@app.get("/registered-mcp-servers")
async def registered_servers_endpoint(request: Request):
    """
    Get list of all registered MCP servers and their status.
    Passes request object to enable bearer token extraction.
    """
    return await get_registered_servers(args.tool_server_cache, request)


manual_mcp_servers_env = os.getenv("FLASK_MCP_SERVERS", "")
if manual_mcp_servers_env:
    manual_mcp_servers = [
        url.strip() for url in manual_mcp_servers_env.split(",") if url.strip()
    ]
    count = 0
    for url in manual_mcp_servers:
        status = register_url(args.tool_server_cache, url, f"m{count}")
        logger.info(f"{status}")
        count += 1

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

        # Use centralized config resolution - Note that the browser state is not available here
        orchestrator_config = resolve_orchestrator_config(
            requested_backend=args.backend,
            requested_model=args.model,
        )

        html = html.replace(
            "<!-- APP CONFIG -->",
            f"""
           <script>
           window.APP_CONFIG = {{
               WS_SERVER: '{os.getenv("WS_SERVER", "ws://localhost:8001/ws")}',
               VERSION: '{os.getenv("SERVER_VERSION", "")}',
               ORCHESTRATOR: {json.dumps(orchestrator_config)}
           }};
           </script>""",
        )
        return HTMLResponse(html)


@app.websocket("/ws")
async def websocket_endpoint(websocket: WebSocket):
    logger.info(f"Request for websocket received. Headers: {str(websocket.headers)}")

    username = "nobody"
    if "x-forwarded-user" in websocket.headers:
        username = websocket.headers["x-forwarded-user"]

    await websocket.accept()

    user_session: FlaskUserSession | None = cast(
        FlaskUserSession | None,
        UserSessionManager.get_latest_inactive_session(username),
    )

    # If an inactive session exists, attach to it
    if user_session is not None:
        # In a later PR, this branch will send the saved state back to the UI in
        # order to resume the session. For now, we terminate the existing
        # session.
        await user_session.terminate()
        user_session = None
        # try:
        #     await user_session.websocket.set_websocket(websocket)
        #     logger.info("User session refreshed with new websocket")
        # except SessionTimedOut:
        #     await user_session.terminate()
        #     user_session = None

    # Otherwise, create a new user session.
    if user_session is None:
        user_session = FlaskUserSession(username, args, websocket)
        logger.info("User session created")

    try:
        await user_session.event_loop()
    finally:
        logger.info("User session websocket loop exited")


if __name__ == "__main__":
    import uvicorn

    host = args.host
    if host is None:
        _, host = try_get_public_hostname()

    # Note: the agent backend is no longer registered globally here. Each user
    # session creates and owns its own agent backend (see ActionManager) so that
    # per-user model/credential configuration cannot leak across sessions.

    uvicorn.run(
        app,
        host=host,
        port=args.port,
        ws_ping_timeout=60.0,  # Increase websocket ping timeout to 60 seconds
        timeout_keep_alive=75,  # Keep connections alive for 75 seconds
    )
