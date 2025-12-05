################################################################################
## Copyright 2025 Lawrence Livermore National Security, LLC. and Binghamton University.
## See the top-level LICENSE file for details.
##
## SPDX-License-Identifier: Apache-2.0
################################################################################

from functools import partial
from fastapi import WebSocket, WebSocketDisconnect
import os
import argparse
import httpx
import os
from loguru import logger

from charge.servers.server_utils import try_get_public_hostname
from charge.clients.Client import Client
from charge.experiments.AutoGenExperiment import AutoGenExperiment
from charge.clients.autogen import AutoGenPool


from tool_registration import (
    register_url,
    register_post,
    reload_server_list,
    split_url,
)

from backend_manager import TaskManager, ActionManager
from backend.server import app
from backend.routers import projects as projrouter

app.include_router(projrouter.router)

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
Client.add_std_parser_arguments(parser, defaults=dict(backend="openai", model="gpt-5-nano"))

args, _ = parser.parse_known_args()

reload_server_list(args.tool_server_cache)

app.post("/register")(partial(register_post, args.tool_server_cache))

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

    autogen_pool = AutoGenPool(model=model, backend=backend, api_key=API_KEY, base_url=BASE_URL)

    # Set up an experiment class for current endpoint
    experiment = AutoGenExperiment(task=None, agent_pool=autogen_pool)

    task_manager = TaskManager(websocket)

    action_manager = ActionManager(task_manager, experiment, args, username)
    await action_manager.report_orchestrator_config()

    action_handlers = {
        "compute": action_manager.handle_compute,
        "compute-reaction-from": action_manager.handle_compute_reaction_from,
        "optimize-from": action_manager.handle_optimize_from,
        "recompute-reaction": action_manager.handle_recompute_reaction,
        "list-tools": action_manager.handle_list_tools,
        "select-tools-for-task": action_manager.handle_select_tools_for_task,
        "update-profile-settings": action_manager.handle_profile_update,
        "query-retro-product": action_manager.handle_custom_query_retro_product,
        "query-retro-reactant": action_manager.handle_custom_query_retro_reactant,
        "query-retro-molecule": action_manager.handle_custom_query_retro_molecule,
        "reset": action_manager.handle_reset,
        "save-context": action_manager.handle_save_state,
        "load-context": action_manager.handle_load_state,
        "stop": action_manager.handle_stop,
        "get-username": action_manager.handle_get_username,
    }

    try:
        while True:
            try:
                data = await websocket.receive_json()
                action = data.get("action")
                if action in action_handlers:

                    # Cancel any existing task before starting a new one
                    await task_manager.cancel_current_task()

                    handler_func = action_handlers[action]
                    await handler_func(data)
                else:
                    logger.warning(f"Unknown action received: {action} with data {data}")
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
