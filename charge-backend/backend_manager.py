from typing import Optional, Tuple
from fastapi import WebSocket
import asyncio
from loguru import logger
from callback_logger import CallbackLogger
from concurrent.futures import ProcessPoolExecutor
from charge.experiments.AutoGenExperiment import AutoGenExperiment
from charge.clients.autogen_utils import chargeConnectionError
from backend_helper_funcs import (
    RetrosynthesisContext,
)
from lmo_charge_backend_funcs import generate_lead_molecule
from functools import partial
from tool_registration import (
    ToolList,
    list_server_urls,
    list_server_tools,
)
from retro_charge_backend_funcs import (
    generate_molecules,
    optimize_molecule_retro,
)


class TaskManager:
    """Manages background tasks and processes state for a websocket connection."""

    def __init__(self, websocket: WebSocket, max_workers: int = 4):
        self.websocket = websocket
        self.current_task: Optional[asyncio.Task] = None
        self.clogger = CallbackLogger(websocket)
        self.max_workers = max_workers
        self.executor = ProcessPoolExecutor(max_workers=max_workers)

    def _attach_done_callback(self, task: asyncio.Task) -> None:
        """Attach a done-callback to a background task so exceptions are observed.

        The callback forwards useful error metadata to the websocket and logs
        the exception type/module so class-identity mismatches (multiple
        installations of `charge`) can be diagnosed.
        """
        if task is None:
            return
        task.add_done_callback(lambda t: asyncio.create_task(self._handle_task_done(t)))

    async def _handle_task_done(self, task: asyncio.Task) -> None:
        try:
            exc = task.exception()
        except asyncio.CancelledError:
            logger.info("Background task was cancelled")
            return

        if exc is None:
            return

        if type(exc) == chargeConnectionError:
            # logger.error(f"Charge connection error in background task: {exc}")
            self.clogger.info(
                f"Unsupported model was selected.  \n Server encountered error: {exc}"
            )

        # Send a stopped message with error details to the websocket so the UI can react
        try:
            await self.websocket.send_json({"type": "complete"})
        except Exception:
            logger.exception("Failed to send task error to websocket")

    async def run_task(self, coro) -> None:
        await self.cancel_current_task()
        self.current_task = asyncio.create_task(coro)
        self._attach_done_callback(self.current_task)

    async def cancel_current_task(self) -> None:
        if self.current_task and not self.current_task.done():
            logger.info("Cancelling current task...")
            self.current_task.cancel()
            try:
                await self.current_task
            except asyncio.CancelledError:
                logger.info("Current task cancelled successfully.")
        await self.restart_executor()

    async def restart_executor(self) -> None:
        """Shutdown and recreate the process pool executor."""
        self.executor.shutdown(wait=False, cancel_futures=True)
        self.executor = ProcessPoolExecutor(max_workers=self.max_workers)

    async def close(self) -> None:
        await self.cancel_current_task()
        self.executor.shutdown(wait=False, cancel_futures=True)
        self.clogger.unbind()


class ActionManager:
    """Handles action state for a websocket connection."""

    def __init__(self, task_manager: TaskManager, experiment: AutoGenExperiment, args, username: str):
        self.task_manager = task_manager
        self.experiment = experiment
        self.args = args
        self.username = username
        self.websocket = task_manager.websocket

    def setup_retro_synth_context(self) -> None:
        if self.retro_synth_context is None:
            self.retro_synth_context = RetrosynthesisContext()

    def get_retro_synth_context(self) -> RetrosynthesisContext:
        self.setup_retro_synth_context()
        assert self.retro_synth_context is not None
        return self.retro_synth_context

    async def handle_compute(self, data: dict) -> None:
        problem_type = data.get("problemType")
        if problem_type == "optimization":
            await self._handle_optimization(data)
        elif problem_type == "retrosynthesis":
            await self._handle_retrosynthesis(data)
        else:
            raise ValueError(f"Unknown problem type: {problem_type}")

    async def handle_save_state(self, *args, **kwargs) -> None:
        """Handle save state action."""
        logger.info("Save state action received")

        experiment_context = await self.experiment.save_state()
        await self.websocket.send_json(
            {"type": "save-context-response", "experimentContext": experiment_context}
        )

    async def handle_load_state(self, data, *args, **kwargs) -> None:
        """Handle load state action."""
        logger.info("Load state action received")
        experiment_context = data.get("experimentContext")
        if not experiment_context:
            logger.error("No experiment context provided for loading state")
            return
        await self.experiment.load_state(experiment_context)

    async def _handle_optimization(self, data: dict) -> None:
        """Handle optimization problem type."""
        self.task_manager.clogger.info("Start Optimization action received")
        logger.info(f"Data: {data}")

        run_func = partial(
            generate_lead_molecule,
            data["smiles"],
            self.experiment,
            self.args.json_file,
            self.args.max_iterations,
            data.get("depth", 3),
            list_server_urls(),
            self.task_manager.websocket,
        )

        await self.task_manager.run_task(run_func())

    async def _handle_retrosynthesis(self, data: dict) -> None:
        """Handle retrosynthesis problem type."""
        self.task_manager.clogger.info("Setting up retrosynthesis task...")
        logger.info(f"Data: {data}")

        run_func = partial(
            generate_molecules,
            data["smiles"],
            self.args.config_file,
            self.get_retro_synth_context(),
            self.task_manager.executor,
            self.task_manager.websocket,
        )

        await self.task_manager.run_task(run_func())

    async def handle_compute_reaction_from(self, data: dict) -> None:
        """Handle compute-reaction-from action."""
        if self.retro_synth_context is None:
            raise ValueError("Retrosynthesis context not initialized")

        logger.info("Synthesize tree leaf action received")
        logger.info(f"Data: {data}")

        run_func = partial(
            optimize_molecule_retro,
            data["nodeId"],
            self.retro_synth_context,
            self.task_manager.websocket,
            self.experiment,
            list_server_urls(),
        )

        await self.task_manager.run_task(run_func())

        await self.websocket.send_json({"type": "complete"})

    async def handle_optimize_from(self, data: dict) -> None:
        """Handle optimize-from action."""
        prompt = data.get("query")
        if prompt:
            await self._send_processing_message(
                f"Processing optimization query: {prompt} for node {data['nodeId']}"
            )

        logger.info("Optimize from action received")
        logger.info(f"Data: {data}")
        # TODO: Implement optimization logic

    async def handle_recompute_reaction(self, data: dict) -> None:
        """Handle recompute-reaction action."""
        prompt = data.get("query")
        if prompt:
            await self._send_processing_message(
                f"Processing reaction query: {prompt} for node {data['nodeId']}"
            )

        logger.info("Recompute reaction action received")
        logger.info(f"Data: {data}")
        await self.websocket.send_json({"type": "complete"})

    async def _send_processing_message(self, message: str) -> None:
        """Send a processing message to the client."""
        await self.websocket.send_json(
            {
                "type": "response",
                "message": {
                    "source": "System",
                    "message": message,
                },
            }
        )

    async def handle_list_tools(self, *args, **kwargs) -> None:
        tools = []
        server_list = list_server_urls()
        for server in server_list:
            tool_list = await list_server_tools([server])
            tool_names = [name for name, _ in tool_list]
            tools.append(ToolList(server=server, names=tool_names))
        await self.websocket.send_json(
            {
                "type": "available-tools-response",
                "tools": [tool.json() for tool in tools] if tools else [],
            }
        )

    async def report_orchestrator_config(self) -> Tuple[str, str, str]:
        agent_pool = self.experiment.agent_pool
        # Access the raw config
        raw_config = agent_pool.model_client._raw_config
        # Access specific fields
        base_url = raw_config.get("base_url")
        model = raw_config.get("model")
        api_key = raw_config.get("api_key")
        await self.websocket.send_json(
            {
                "type": "update-orchestrator-profile",
                "profileSettings": {
                    "backend": agent_pool.backend,
                    "useCustomUrl": False,
                    "customUrl": base_url if base_url else "",
                    "model": model,
                    "useCustomModel": False,
                    "apiKey": "",
                },
            }
        )
        return agent_pool.backend, model, base_url

    async def handle_profile_update(self, data: dict) -> None:
        from charge.experiments.AutoGenExperiment import AutoGenExperiment
        from charge.clients.autogen import AutoGenPool

        backend = data["backend"]
        model = data["model"]
        base_url = data["customUrl"] if data["customUrl"] else None
        api_key = data["apiKey"] if data["apiKey"] else None
        await self.handle_reset()
        try:
            logger.info(f"Experiment is reset with model {model} and backend {backend}")
            autogen_pool = AutoGenPool(
                model=model, backend=backend, api_key=api_key, base_url=base_url
            )
            # Set up an experiment class for current endpoint
            self.experiment = AutoGenExperiment(task=None, agent_pool=autogen_pool)

            await self.websocket.send_json(
                {
                    "type": "response",
                    "message": {
                        "source": "System",
                        "message": f"Experiment is reset with model {model} and backend {backend}",
                    },
                }
            )
        except ValueError as e:
            logger.error(
                f"Orchestrator Profile Error: Unable to restart experiment: {e}"
            )
            backend, model, base_url = await self.report_orchestrator_config()
            await self.websocket.send_json(
                {
                    "type": "response",
                    "message": {
                        "source": "System",
                        "message": f"Orchestrator Profile Error: Unable to restart experiment: {e}. Experiment is still using backend {backend} with model {model} at {base_url}",
                    },
                }
            )

    async def handle_reset(self, *args, **kwargs) -> None:
        """Handle reset action."""
        await self.task_manager.cancel_current_task()
        self.experiment.reset()
        self.retro_synth_context = None

    async def handle_stop(self, *args, **kwargs) -> None:
        """Handle stop action."""
        if self.task_manager.current_task and not self.task_manager.current_task.done():
            logger.info("Stopping current task as per user request.")
            await self.task_manager.cancel_current_task()
        else:
            logger.info(
                f"No active task to stop. Task done: {self.task_manager.current_task.done() if self.task_manager.current_task else 'N/A'}"
            )

    async def handle_select_tools_for_task(self, data: dict) -> None:
        """Handle select-tools-for-task action."""
        query = data.get("query")
        logger.info("Select tools for task")
        logger.info(f"Data: {data}")
        # TODO: Implement tool selection logic

    async def handle_custom_query(self, data: dict) -> None:
        """Handle custom_query action."""
        await self._send_processing_message(
            f"Processing query: {data['query']} for node {data['nodeId']}"
        )
        await asyncio.sleep(3)
        await self.websocket.send_json({"type": "complete"})

    async def handle_get_username(self, _: dict) -> None:
        await self.websocket.send_json(
            {
                "type": "get-username-response",
                "username": self.username,
            }
        )
