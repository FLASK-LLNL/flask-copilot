################################################################################
## Copyright 2025 Lawrence Livermore National Security, LLC. and Binghamton University.
## See the top-level LICENSE file for details.
##
## SPDX-License-Identifier: Apache-2.0
################################################################################

from functools import partial
from fastapi import FastAPI, WebSocket, WebSocketDisconnect
from fastapi.responses import HTMLResponse
from fastapi.staticfiles import StaticFiles
from fastapi.middleware.cors import CORSMiddleware
import asyncio
import os
import argparse
import httpx
from charge.servers.server_utils import try_get_public_hostname
import os


import logging
from aizynthfinder.utils.logging import setup_logger

setup_logger(console_level=logging.INFO)

from loguru import logger
from callback_logger import CallbackLogger

from backend_helper_funcs import (
    CallbackHandler,
    RetrosynthesisContext,
    Node,
    Edge,
    calculate_positions,
)
from retro_charge_backend_funcs import (
    generate_molecules,
    optimize_molecule_retro,
)
from lmo_charge_backend_funcs import generate_lead_molecule
from charge.clients.Client import Client
from charge.experiments.AutoGenExperiment import AutoGenExperiment
from charge.clients.autogen import AutoGenPool

from tool_registration import (
    register_post,
    list_server_urls,
    list_server_tools,
    reload_server_list,
)

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
    async def root():
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


async def _cancel_task_if_running(action, task: asyncio.Task | None):
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


@app.websocket("/ws")
async def websocket_endpoint(websocket: WebSocket):
    await websocket.accept()

    # set up an AutoGenAgent pool for tasks on this endpoint

    autogen_pool = AutoGenPool(model=args.model, backend=args.backend)
    # Set up an experiment class for current endpoint
    experiment = AutoGenExperiment(task=None, agent_pool=autogen_pool)

    # Keep track of currently running task
    CURRENT_TASK: asyncio.Task | None = None

    clogger = CallbackLogger(websocket)

    retro_synth_context: RetrosynthesisContext | None = None

    while True:
        try:
            data = await websocket.receive_json()
            action = data.get("action")
            await _cancel_task_if_running(action, CURRENT_TASK)

            if action == "compute":

                if data["problemType"] == "optimization":

                    # Task to optimize lead molecule using LMO
                    clogger.info("Start Optimization action received")
                    logger.info(f"Data: {data}")

                    run_func = partial(
                        generate_lead_molecule,
                        data["smiles"],
                        experiment,
                        args.json_file,
                        args.max_iterations,
                        data.get("depth", 3),
                        list_server_urls(),  # TODO: This should be changed to available specified by the user
                        websocket,
                    )
                elif data["problemType"] == "retrosynthesis":
                    # Set up retrosynthesis task to retrosynthesis
                    # to ensure the reactant is not used in the synthesis
                    clogger.info("Setting up retrosynthesis task...")
                    logger.info(f"Data: {data}")

                    if retro_synth_context is None:
                        retro_synth_context = RetrosynthesisContext()

                    assert retro_synth_context is not None
                    run_func = partial(
                        generate_molecules,
                        data["smiles"],
                        args.config_file,
                        retro_synth_context,
                        websocket,
                    )
                else:
                    raise ValueError(f"Unknown problem type: {data['problemType']}")

                async def run_task():
                    await run_func()

                # start a new task
                CURRENT_TASK = asyncio.create_task(run_task())

            elif action == "compute-reaction-from":
                assert retro_synth_context is not None

                # Leaf node optimization
                logger.info("Synthesize tree leaf action received")
                logger.info(f"Data: {data}")
                await optimize_molecule_retro(
                    data["nodeId"],
                    retro_synth_context,
                    websocket,
                    experiment,
                    list_server_urls(),
                )
                await websocket.send_json({"type": "complete"})
            elif action == "optimize-from":
                # Leaf node optimization
                prompt = data.get("query", None)
                if prompt:
                    await websocket.send_json(
                        {
                            "type": "response",
                            "message": {
                                "source": "System",
                                "message": f"Processing optimization query: {data['query']} for node {data['nodeId']}",
                            },
                        }
                    )
                logger.info("Optimize from action received")
                logger.info(f"Data: {data}")

                pass
            elif action == "recompute-reaction":
                prompt = data.get("query", None)
                if prompt:
                    await websocket.send_json(
                        {
                            "type": "response",
                            "message": {
                                "source": "System",
                                "message": f"Processing reaction query: {data['query']} for node {data['nodeId']}",
                            },
                        }
                    )

                logger.info("Recompute reaction action received")
                logger.info(f"Data: {data}")
                await websocket.send_json({"type": "complete"})

            elif action == "list-tools":
                tools = []
                server_list = list_server_urls()
                for server in server_list:
                    tool_list = await list_server_tools([server])
                    tool_names = []
                    for name, _ in tool_list:
                        tool_names.append(name)
                    tools.append(ToolServer(server, tool_names))

                if tools == []:
                    await websocket.send_json(
                        {
                            "type": "available-tools-response",
                            "tools": [],
                        }
                    )
                else:
                    await websocket.send_json(
                        {
                            "type": "available-tools-response",
                            "tools": [tool.json() for tool in tools],
                        }
                    )
            elif action == "select-tools-for-task":
                query = data.get("query", None)
                logger.info("Select tools for task")
                logger.info(f"Data: {data}")
            elif action == "custom_query":
                await websocket.send_json(
                    {
                        "type": "response",
                        "message": {
                            "source": "System",
                            "message": f"Processing query: {data['query']} for node {data['nodeId']}",
                        },
                    }
                )
                await asyncio.sleep(3)
                await websocket.send_json({"type": "complete"})

            elif action == "reset":
                experiment.reset()

                if retro_synth_context is not None:
                    retro_synth_context.reset()
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
            else:
                logger.warning(f"Unknown action received: {action}")

        except WebSocketDisconnect:
            logger.info("WebSocket disconnected")
            break
        except httpx.ConnectError as e:
            logger.error(f"Connection error: {e}")
        except Exception as e:
            logger.error(f"Error in WebSocket connection: {e}")
        finally:
            if CURRENT_TASK and not CURRENT_TASK.done():
                logger.info("Stopping current task due to connection closure.")
                CURRENT_TASK.cancel()
                try:
                    await CURRENT_TASK
                except asyncio.CancelledError:
                    logger.info("Current task cancelled successfully.")

        clogger.unbind()


if __name__ == "__main__":
    import uvicorn

    host = args.host
    if host is None:
        _, host = try_get_public_hostname()

    uvicorn.run(app, host=host, port=args.port)
