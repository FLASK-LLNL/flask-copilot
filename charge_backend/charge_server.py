################################################################################
## Copyright 2025 Lawrence Livermore National Security, LLC. and Binghamton University.
## See the top-level LICENSE file for details.
##
## SPDX-License-Identifier: Apache-2.0
################################################################################

import asyncio
from concurrent.futures import ProcessPoolExecutor
from functools import partial
from fastapi import FastAPI, WebSocket, WebSocketDisconnect, Request
from fastapi.responses import HTMLResponse
from fastapi.staticfiles import StaticFiles
from fastapi.middleware.cors import CORSMiddleware
import os
import argparse
import httpx
from flask_mcp.utils.server_utils import try_get_public_hostname
import os


import logging
from aizynthfinder.utils.logging import setup_logger

from loguru import logger
from charge.clients.Client import Client
from charge.experiments.AutoGenExperiment import AutoGenExperiment
from charge.clients.autogen import AutoGenPool


from lc_conductor.tool_registration import (
    get_client_info,
    register_url,
    register_post,
    reload_server_list,
    split_url,
    validate_mcp_server_endpoint,
    delete_mcp_server_endpoint,
    get_registered_servers,
)

from lc_conductor import TaskManager
from backend_manager import FlaskActionManager
from charge_backend import prompt_debugger

# Pydantic models for new endpoints
from pydantic import BaseModel
from typing import Optional


class CheckServersRequest(BaseModel):
    urls: list[str]


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
Client.add_std_parser_arguments(
    parser, defaults=dict(backend="openai", model="gpt-5.1")
)

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


@app.post("/check-mcp-servers")
async def check_mcp_servers_endpoint(data: CheckServersRequest):
    """
    Check connectivity status of multiple MCP server URLs.
    Returns status and tools for each URL.

    Uses existing workbench utilities for validation.
    """
    from lc_conductor.tool_registration import _check_mcp_connectivity

    results = {}

    for url in data.urls:
        try:
            tools = await _check_mcp_connectivity(url, timeout=5.0)
            results[url] = {"status": "connected", "tools": tools}
        except Exception as e:
            results[url] = {"status": "disconnected", "error": str(e)}

    return {"results": results}


app.get("/registered-mcp-servers")(
    partial(get_registered_servers, args.tool_server_cache)
)

manual_mcp_servers_env = os.getenv("FLASK_MCP_SERVERS", "")
if manual_mcp_servers_env:
    manual_mcp_servers = manual_mcp_servers_env.split(",")
    count = 0
    for url in manual_mcp_servers:
        host, port, path, protocol = split_url(url)
        status = register_url(
            args.tool_server_cache, host, port, path, protocol, f"m{count}"
        )
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


async def _cancel_task_if_running(
    action, task: asyncio.Task | None, executor: ProcessPoolExecutor | None
):
    if action not in [
        "compute",
        "optimize-from",
        "compute-reaction-from",
        "recompute-reaction",
        "recompute-parent-reaction",
    ]:
        return

    if task and not task.done():
        logger.info("Cancelling existing compute task...")
        task.cancel()
        try:
            await task
        except asyncio.CancelledError:
            logger.info("Previous compute task cancelled.")

    if executor:
        executor.shutdown(wait=False, cancel_futures=True)


@app.websocket("/ws")
async def websocket_endpoint(websocket: WebSocket):
    logger.info(f"Request for websocket received. Headers: {str(websocket.headers)}")

    username = "nobody"
    if "x-forwarded-user" in websocket.headers:
        username = websocket.headers["x-forwarded-user"]

    await websocket.accept()

    API_KEY = os.getenv("FLASK_ORCHESTRATOR_API_KEY", None)
    model = os.getenv("FLASK_ORCHESTRATOR_MODEL", None)
    backend = os.getenv("FLASK_ORCHESTRATOR_BACKEND", None)
    BASE_URL = os.getenv("FLASK_ORCHESTRATOR_URL", None)

    if not model:
        model = args.model
    if not backend:
        backend = args.backend

    # set up an AutoGenAgent pool for tasks on this endpoint

    autogen_pool = AutoGenPool(
        model=model, backend=backend, api_key=API_KEY, base_url=BASE_URL
    )

    # Set up an experiment class for current endpoint
    experiment = AutoGenExperiment(task=None, agent_pool=autogen_pool)

    task_manager = TaskManager(websocket)

    action_manager = FlaskActionManager(task_manager, experiment, args, username)
    await action_manager.report_orchestrator_config()

    action_handlers = {
        # General actions
        "query-molecule": action_manager.handle_custom_query_molecule,
        "query-reaction": action_manager.handle_custom_query_reaction,
        "compute": action_manager.handle_compute,
        "reset": action_manager.handle_reset,
        "stop": action_manager.handle_stop,
        # Tools
        "list-tools": action_manager.handle_list_tools,
        "select-tools-for-task": action_manager.handle_select_tools_for_task,
        # Settings
        "ui-update-orchestrator-settings": action_manager.handle_orchestrator_settings_update,
        "get-username": action_manager.handle_get_username,
        # Context management
        "save-context": action_manager.handle_save_state,
        "load-context": action_manager.handle_load_state,
        # Lead molecule optimization
        "optimize-from": action_manager.handle_optimize_from,
        # Retrosynthesis
        "compute-reaction-from": action_manager.handle_compute_reaction_from,
        "compute-reaction-templates": action_manager.handle_template_retrosynthesis,
        "recompute-parent-reaction": action_manager.handle_recompute_reaction,
        "set-reaction-alternative": action_manager.handle_set_reaction_alternative,
    }

    try:
        while True:
            try:
                data = await websocket.receive_json()
                action = data.get("action")

                if action == "prompt-breakpoint-response":  # AI debugging
                    prompt_debugger.DEBUG_PROMPT_RESPONSES[websocket].set_result(data)
                    continue

                if action in action_handlers:

                    # Cancel any existing task before starting a new one
                    await task_manager.cancel_current_task()

                    handler_func = action_handlers[action]
                    await handler_func(data)
                else:
                    logger.warning(
                        f"Unknown action received: {action} with data {data}"
                    )
            except ValueError as e:
                logger.error(f"Error in internal loop connection: {e}")
                await task_manager.cancel_current_task()

    except WebSocketDisconnect:
        logger.info("WebSocket disconnected")
        pass
    except httpx.ConnectError as e:
        logger.error(f"Connection error: {e}")
    except Exception as e:
        logger.error(f"Error in WebSocket connection: {e}")
    finally:
        await task_manager.close()


if __name__ == "__main__":
    import uvicorn

    host = args.host
    if host is None:
        _, host = try_get_public_hostname()

    uvicorn.run(app, host=host, port=args.port)
