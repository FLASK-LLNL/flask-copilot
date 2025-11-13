################################################################################
## Copyright 2025 Lawrence Livermore National Security, LLC. and Binghamton University.
## See the top-level LICENSE file for details.
##
## SPDX-License-Identifier: Apache-2.0
################################################################################

from concurrent.futures import ProcessPoolExecutor
from functools import partial
from fastapi import FastAPI, WebSocket, WebSocketDisconnect, Request
from fastapi.responses import HTMLResponse
from fastapi.staticfiles import StaticFiles
from fastapi.middleware.cors import CORSMiddleware
import os
import argparse
import httpx
from charge.servers.server_utils import try_get_public_hostname
import os


import logging
from aizynthfinder.utils.logging import setup_logger

setup_logger(console_level=logging.INFO)

from loguru import logger
from charge.clients.Client import Client
from charge.experiments.AutoGenExperiment import AutoGenExperiment
from charge.clients.autogen import AutoGenPool

from tool_registration import (
    register_post,
    reload_server_list,
)

from backend_manager import TaskManager, ActionManager

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

parser.add_argument(
    "--tool-server-cache",
    type=str,
    default="flask_copilot_active_tool_servers.json",
    help="Path to the JSON file containing current list of active MCP tool servers.",
)

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

if "FLASK_APPDIR" in os.environ:
    DIST_PATH = os.environ["FLASK_APPDIR"]
else:
    DIST_PATH = os.path.join(os.path.dirname(__file__), "flask-app", "dist")
ASSETS_PATH = os.path.join(DIST_PATH, "assets")

reload_server_list(args.tool_server_cache)

app.post("/register")(partial(register_post, args.tool_server_cache))

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
               WS_SERVER: '{os.getenv("WS_SERVER", "ws://localhost:8001/ws")}'
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
    await websocket.accept()

    # set up an AutoGenAgent pool for tasks on this endpoint
    autogen_pool = AutoGenPool(model=args.model, backend=args.backend)
    # Set up an experiment class for current endpoint
    experiment = AutoGenExperiment(task=None, agent_pool=autogen_pool)

    task_manager = TaskManager(websocket)

    action_manager = ActionManager(task_manager, experiment, args)

    action_handlers = {
        "compute": action_manager.handle_compute,
        "compute-reaction-from": action_manager.handle_compute_reaction_from,
        "optimize-from": action_manager.handle_optimize_from,
        "recompute-reaction": action_manager.handle_recompute_reaction,
        "list-tools": action_manager.handle_list_tools,
        "select-tools-for-task": action_manager.handle_select_tools_for_task,
        "custom_query": action_manager.handle_custom_query,
        "reset": action_manager.handle_reset,
        "stop": action_manager.handle_stop,
    }

    try:
        while True:
            try:
                data = await websocket.receive_json()
                action = data.get("action")
                if action in action_handlers:
                    if action != "stop":
                        await task_manager.cancel_current_task()
                    handler_func = action_handlers[action]
                    await handler_func(data)
                else:
                    logger.warning(f"Unknown action received: {action}")
            except ValueError as e:
                logger.error(f"Error in internal loop connection: {e}")
                await task_manager.cancel_current_task()
                # await

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
