from typing import Literal, Optional, Tuple
from fastapi import WebSocket
import asyncio
import os
from loguru import logger
from callback_logger import CallbackLogger
from concurrent.futures import ProcessPoolExecutor
from charge.experiments.AutoGenExperiment import AutoGenExperiment
from charge.clients.autogen_utils import chargeConnectionError
from charge.tasks.Task import Task
from backend_helper_funcs import (
    CallbackHandler,
    Reaction,
)
from retrosynthesis.context import RetrosynthesisContext
from lmo_charge_backend_funcs import generate_lead_molecule
from charge_backend_custom import run_custom_problem
from functools import partial
from tool_registration import (
    ToolList,
    list_server_urls,
    list_server_tools,
)
from retro_charge_backend_funcs import (
    template_based_retrosynthesis,
    ai_based_retrosynthesis,
    compute_templates_for_node,
    set_reaction_alternative,
)

# Mapping from backend name to human-readable labels. Mirrored from the frontend
BACKEND_LABELS = {
    "openai": "OpenAI",
    "livai": "LivAI",
    "llamame": "LLamaMe",
    "alcf": "ALCF Sophia",
    "gemini": "Google Gemini",
    "ollama": "Ollama",
    "vllm": "vLLM",
    "huggingface": "HuggingFace Local",
    "custom": "Custom URL",
}


class TaskManager:
    """Manages background tasks and processes state for a websocket connection."""

    def __init__(self, websocket: WebSocket, max_workers: int = 4):
        self.websocket = websocket
        self.current_task: Optional[asyncio.Task] = None
        self.clogger = CallbackLogger(websocket, source="backend_manager")
        self.max_workers = max_workers
        self.executor = ProcessPoolExecutor(max_workers=max_workers)
        self.available_tools: Optional[list[str]] = None

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

        # Log the exception details
        msg = f"Background task failed with exception: {type(exc).__name__}: {exc}"
        logger.error(msg)
        await self.websocket.send_json(
            {
                "type": "response",
                "message": {
                    "source": "system",
                    "message": msg,
                },
            }
        )

        if type(exc) == chargeConnectionError:
            # logger.error(f"Charge connection error in background task: {exc}")
            await self.clogger.info(
                f"Unsupported model was selected.  \n Server encountered error: {exc}"
            )
        else:
            # Log other exceptions for debugging
            logger.exception(
                f"Unexpected error in background task: {type(exc).__name__}: {exc}"
            )

        # Send a stopped message with error details to the websocket so the UI can react
        try:
            await self.websocket.send_json({"type": "complete"})
        except Exception as send_error:
            logger.exception(f"Failed to send task error to websocket: {send_error}")

    async def run_task(self, coro) -> None:
        await self.cancel_current_task()
        try:
            self.current_task = asyncio.create_task(coro)
            self._attach_done_callback(self.current_task)
            await self.current_task  # Await it to catch exceptions properly
        except asyncio.CancelledError:
            logger.info("Task was cancelled")
            raise
        except Exception as e:
            logger.error(f"Task failed: {e}")
            await self.websocket.send_json({"type": "complete"})
            # Optionally re-raise or handle as needed

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

    def __init__(
        self,
        task_manager: TaskManager,
        experiment: AutoGenExperiment,
        args,
        username: str,
    ):
        self.task_manager = task_manager
        self.experiment = experiment
        self.args = args
        self.username = username
        self.molecule_name_format: Literal["brand", "iupac", "formula", "smiles"] = (
            "brand"
        )
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
            asyncio.create_task(self._handle_optimization(data))
        elif problem_type == "retrosynthesis":
            asyncio.create_task(self._handle_retrosynthesis(data))
        elif problem_type == "custom":
            asyncio.create_task(self._handle_custom_problem(data))
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

    async def _handle_optimization(
        self,
        data: dict,
        initial_level: int = 0,
        initial_node_id: int = 0,
        initial_x_position: int = 50,
    ) -> None:
        """Handle optimization problem type."""
        await self.task_manager.clogger.info("Start Optimization action received")
        logger.info(f"Data: {data}")

        # Validate required field
        if "propertyType" not in data:
            error_msg = "Missing required field 'propertyType' for optimization problem"
            logger.error(error_msg)
            await self.websocket.send_json(
                {
                    "type": "response",
                    "message": {
                        "source": "System",
                        "message": f"Error: {error_msg}",
                    },
                }
            )
            await self.websocket.send_json({"type": "complete"})
            return

        property_type = data["propertyType"]

        # Property attributes for prompting
        if property_type == "custom":
            # Validate custom property fields
            if not data.get("customPropertyName") or not data.get("customPropertyDesc"):
                error_msg = "Custom property requires both 'customPropertyName' and 'customPropertyDesc'"
                logger.error(error_msg)
                await self.websocket.send_json(
                    {
                        "type": "response",
                        "message": {
                            "source": "System",
                            "message": f"Error: {error_msg}",
                        },
                    }
                )
                await self.websocket.send_json({"type": "complete"})
                return

            property_attributes = (
                data["customPropertyName"],
                data["customPropertyDesc"],
                "calculate_property_hf",
                "greater" if data["customPropertyAscending"] else "less",
            )
        else:
            property_name = property_type
            DEFAULT_PROPERTIES = {
                "density": (
                    "density",
                    "crystalline density (g/cm^3)",
                    "calculate_property_hf",
                    "greater",
                ),
                "hof": (
                    "heat of formation",
                    "Heat of formation (kcal/mol)",
                    "calculate_property_hf",
                    "greater",
                ),
                "bandgap": (
                    "band gap",
                    "HOMO-LUMO energy gap (Hartree)",
                    "calculate_property_hf",
                    "greater",
                ),
            }
            if property_name not in DEFAULT_PROPERTIES:
                error_msg = (
                    f"Property '{property_name}' not found in default properties"
                )
                logger.error(error_msg)
                await self.websocket.send_json(
                    {
                        "type": "response",
                        "message": {
                            "source": "System",
                            "message": f"Error: {error_msg}",
                        },
                    }
                )
                await self.websocket.send_json({"type": "complete"})
                return
            property_attributes = DEFAULT_PROPERTIES[property_name]

        # Extract customization parameters
        customization = data.get("customization", {})
        enable_constraints = customization.get("enableConstraints", False)
        molecular_similarity = customization.get("molecularSimilarity", 0.7)
        diversity_penalty = customization.get("diversityPenalty", 0.0)
        exploration_rate = customization.get("explorationRate", 0.5)
        additional_constraints = customization.get("additionalConstraints", [])

        run_func = partial(
            generate_lead_molecule,
            data["smiles"],
            self.experiment,
            self.args.json_file,
            self.args.max_iterations,
            data.get("depth", 3),
            self.task_manager.available_tools or list_server_urls(),
            self.task_manager.websocket,
            *property_attributes,
            data.get("query", None),
            initial_level,
            initial_node_id,
            initial_x_position,
            self.molecule_name_format,
            enable_constraints,
            molecular_similarity,
            diversity_penalty,
            exploration_rate,
            additional_constraints,
        )
        await self.task_manager.run_task(run_func())

    async def _handle_retrosynthesis(self, data: dict) -> None:
        """Handle retrosynthesis problem type."""
        available_tools = self.task_manager.available_tools or list_server_urls()
        await self.task_manager.clogger.info(
            f"Setting up retrosynthesis task... with available tools: {available_tools}."
        )
        await self.task_manager.clogger.info(f"Data: {data}")

        run_func = partial(
            template_based_retrosynthesis,
            data["smiles"],
            self.args.config_file,
            self.get_retro_synth_context(),
            self.task_manager.websocket,
            available_tools,
            self.molecule_name_format,
        )

        await self.task_manager.run_task(run_func())

    async def _handle_custom_problem(self, data: dict) -> None:
        """Handle custom problem type."""
        await self.task_manager.clogger.info("Setting up custom task...")
        logger.info(f"Data: {data}")

        run_func = partial(
            run_custom_problem,
            data["smiles"],
            data["systemPrompt"],
            data["userPrompt"],
            self.experiment,
            self.task_manager.available_tools or list_server_urls(),
            self.task_manager.websocket,
            self.molecule_name_format,
        )

        await self.task_manager.run_task(run_func())

    async def handle_compute_reaction_from(self, data: dict) -> None:
        """Handle compute-reaction-from action."""
        if self.retro_synth_context is None:
            raise ValueError("Retrosynthesis context not initialized")

        logger.info("Synthesize tree leaf action received")
        logger.info(f"Data: {data}")
        node = self.retro_synth_context.node_ids[data["nodeId"]]

        has_children = any(
            v == data["nodeId"] for v in self.retro_synth_context.parents.values()
        )

        await self.retro_synth_context.delete_subtree(data["nodeId"], self.websocket)
        await self.websocket.send_json(
            {
                "type": "node_update",
                "node": {
                    "id": data["nodeId"],
                    "reaction": Reaction(
                        "ai_reaction_0",
                        (
                            "Recomputing with AI orchestrator"
                            if has_children
                            else "Computing with AI orchestrator"
                        ),
                        highlight="red",
                        label="Recomputing" if has_children else "Computing",
                    ).json(),
                },
            }
        )

        run_func = partial(
            ai_based_retrosynthesis,
            data["nodeId"],
            self.retro_synth_context,
            data.get("query", None),
            None,  # Unconstrained
            self.task_manager.websocket,
            self.experiment,
            self.args.config_file,
            self.task_manager.available_tools or list_server_urls(),
            self.molecule_name_format,
        )

        asyncio.create_task(self.task_manager.run_task(run_func()))

    async def handle_template_retrosynthesis(self, data: dict) -> None:
        """Handle compute-reaction-templates action."""
        if self.retro_synth_context is None:
            raise ValueError("Retrosynthesis context not initialized")

        logger.info("Synthesize tree leaf action received")
        logger.info(f"Data: {data}")
        node = self.retro_synth_context.node_ids[data["nodeId"]]

        run_func = partial(
            compute_templates_for_node,
            node,
            self.args.config_file,
            self.retro_synth_context,
            self.task_manager.websocket,
            self.task_manager.available_tools or list_server_urls(),
            self.molecule_name_format,
        )

        await self.task_manager.run_task(run_func())

    async def handle_optimize_from(self, data: dict) -> None:
        """Handle optimize-from action."""
        prompt = data.get("query")
        if prompt:
            await self._send_processing_message(
                f"Processing optimization query: {prompt} for node {data['nodeId']}"
            )
        else:
            await self._send_processing_message(
                f"Processing optimization refinement for node {data['nodeId']}"
            )

        node: str = data["nodeId"]
        if "_" in node:
            node = node[node.rfind("_") + 1 :]
        node_id = int(node)
        # Level is always equal to node ID for now.
        # TODO: Keep LMO context with nodes
        level = node_id
        xpos = data["xpos"]

        asyncio.create_task(
            self._handle_optimization(
                data,
                initial_level=level,
                initial_node_id=node_id,
                initial_x_position=xpos,
            )
        )

    async def handle_recompute_reaction(self, data: dict) -> None:
        """Handle recompute-reaction action."""
        if self.retro_synth_context is None:
            raise ValueError("Retrosynthesis context not initialized")

        # Get parent
        if data["nodeId"] not in self.retro_synth_context.parents:
            await self._send_processing_message(
                f"Cannot find parent reaction for node {data['nodeId']}", source="Agent"
            )
            await self.websocket.send_json({"type": "complete"})

        parent_nodeid = self.retro_synth_context.parents[data["nodeId"]]
        parent_node = self.retro_synth_context.node_ids[parent_nodeid]
        smiles = self.retro_synth_context.node_ids[data["nodeId"]].smiles

        # Clear subtree and levels for layouting
        await self.retro_synth_context.delete_subtree(parent_nodeid, self.websocket)
        await self.websocket.send_json(
            {
                "type": "node_update",
                "node": {
                    "id": parent_nodeid,
                    "reaction": Reaction(
                        "ai_reaction_0",
                        "Recomputing with AI orchestrator",
                        highlight="red",
                        label="Recomputing",
                    ).json(),
                },
            }
        )

        run_func = partial(
            ai_based_retrosynthesis,
            parent_nodeid,
            self.retro_synth_context,
            data.get("query", None),
            smiles,
            self.task_manager.websocket,
            self.experiment,
            self.args.config_file,
            self.task_manager.available_tools or list_server_urls(),
            self.molecule_name_format,
        )

        asyncio.create_task(self.task_manager.run_task(run_func()))

    async def handle_set_reaction_alternative(self, data: dict) -> None:
        """Handle set-reaction-alternative action."""
        if self.retro_synth_context is None:
            raise ValueError("Retrosynthesis context not initialized")

        # Get node
        if data["nodeId"] not in self.retro_synth_context.node_ids:
            await self._send_processing_message(
                f"Cannot find node {data['nodeId']}", source="Agent"
            )
            await self.websocket.send_json({"type": "complete"})
            return
        node = self.retro_synth_context.node_ids[data["nodeId"]]
        alt = data["alternativeId"]
        if not node.reaction or not node.reaction.alternatives:
            await self._send_processing_message(
                f"No alternative found for {data['nodeId']}", source="Agent"
            )
            await self.websocket.send_json({"type": "complete"})
            return

        await set_reaction_alternative(
            node, alt, self.retro_synth_context, self.websocket
        )

    async def _send_processing_message(
        self, message: str, source: str | None = None, **kwargs
    ) -> None:
        """Send a processing message to the client."""
        await self.websocket.send_json(
            {
                "type": "response",
                "message": {
                    "source": source or "System",
                    "message": message,
                },
                **kwargs,
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
        if agent_pool.backend in ["livai", "livchat", "llamame", "alcf"]:
            useCustomUrl = True
        else:
            useCustomUrl = False
        await self.websocket.send_json(
            {
                "type": "server-update-orchestrator-settings",
                "orchestratorSettings": {
                    "backend": agent_pool.backend,
                    "backendLabel": BACKEND_LABELS.get(
                        agent_pool.backend, agent_pool.backend
                    ),
                    "useCustomUrl": useCustomUrl,
                    "customUrl": base_url if base_url else "",
                    "model": model,
                    "useCustomModel": False,
                    "apiKey": "",
                },
            }
        )
        return agent_pool.backend, model, base_url

    async def handle_orchestrator_settings_update(self, data: dict) -> None:
        from charge.experiments.AutoGenExperiment import AutoGenExperiment
        from charge.clients.autogen import AutoGenPool

        if "moleculeName" in data:
            self.molecule_name_format = data["moleculeName"]

        backend = data["backend"]
        model = data["model"]
        base_url = data["customUrl"] if data["customUrl"] else None
        api_key = data["apiKey"] if data["apiKey"] else None
        await self.handle_reset()

        # Default to server defaults
        if backend == os.getenv("FLASK_ORCHESTRATOR_BACKEND", None):
            if not api_key:
                api_key = os.getenv("FLASK_ORCHESTRATOR_API_KEY", None)
            if not base_url:
                base_url = os.getenv("FLASK_ORCHESTRATOR_URL", None)

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
        logger.info("Stop action received")
        if self.task_manager.current_task:
            if not self.task_manager.current_task.done():
                logger.info("Stopping current task as per user request.")
                await self.task_manager.cancel_current_task()

                # Send confirmation to frontend
                try:
                    await self.websocket.send_json({"type": "stopped"})
                    logger.info("Sent 'stopped' confirmation to frontend")
                except Exception as e:
                    logger.error(f"Failed to send stopped confirmation: {e}")
            else:
                logger.info(f"Task already done: {self.task_manager.current_task}")
                await self.websocket.send_json({"type": "stopped"})
        else:
            logger.info(
                f"No active task to stop. Task done: {self.task_manager.current_task.done() if self.task_manager.current_task else 'N/A'}"
            )
            try:
                await self.websocket.send_json({"type": "stopped"})
            except Exception as e:
                logger.error(f"Failed to send stopped confirmation: {e}")

    async def handle_select_tools_for_task(self, data: dict) -> None:
        """Handle select-tools-for-task action."""
        logger.info("Select tools for task")
        logger.info(f"Data: {data}")
        available_tools = []
        for server in data["enabledTools"]["selectedTools"]:
            available_tools.append(server["tool_server"]["server"])
        self.task_manager.available_tools = available_tools

    async def handle_custom_query_molecule(self, data: dict) -> None:
        """Handle a query on the given molecule."""
        await self._send_processing_message(
            f"Processing molecule query: {data['query']} for node {data['nodeId']}"
        )
        smiles = data["smiles"]

        task = Task(
            system_prompt=f"You are a helpful chemical assistant who answers in concise but factual responses. Answer the following query about the molecule given by the SMILES string `{smiles}`.",
            user_prompt=data["query"],
            server_urls=self.task_manager.available_tools or list_server_urls(),
        )

        # Use the full experiment state
        agent = self.experiment.create_agent_with_experiment_state(
            task=task,
            callback=CallbackHandler(self.websocket),
        )

        async def run_and_report():
            result = await agent.run()
            # Report answer
            await self._send_processing_message(result, source="Agent")
            await self.websocket.send_json({"type": "complete"})

        asyncio.create_task(self.task_manager.run_task(run_and_report()))

    async def handle_custom_query_reaction(self, data: dict) -> None:
        """Handle a query on the reaction (from nodeId to its reactants)."""
        assert self.retro_synth_context is not None

        await self._send_processing_message(
            f"Processing reaction query: {data['query']} for node {data['nodeId']}"
        )

        node = self.retro_synth_context.node_ids[data["nodeId"]]

        task = Task(
            system_prompt="",
            user_prompt=data["query"],
            server_urls=self.task_manager.available_tools or list_server_urls(),
        )

        # No charge client (likely only AZF was run), add context from tree
        # if data["nodeId"] not in self.retro_synth_context.node_id_to_charge_client:
        child_nodes = [
            nid
            for nid, p in self.retro_synth_context.parents.items()
            if p == data["nodeId"]
        ]
        reactants = [self.retro_synth_context.node_ids[nid] for nid in child_nodes]
        reactants_str = "\n".join(reactant.smiles for reactant in reactants)
        reaction_str = f"Product: {node.smiles}\nReactants:\n{reactants_str}"

        # Enrich context from prior discovery
        if data["nodeId"] in self.retro_synth_context.node_id_to_reasoning_summary:
            reaction_str += f"\nAdditionally, the following context is given: {self.retro_synth_context.node_id_to_reasoning_summary[data['nodeId']]}"

        task.system_prompt = f"You are a helpful chemical assistant who answers in concise but factual responses. Given the following reaction (as SMILES strings):\n{reaction_str}\n\nAnswer the following query."

        agent = self.experiment.create_agent_with_experiment_state(
            task=task,
            callback=CallbackHandler(self.websocket),
        )
        self.retro_synth_context.node_id_to_charge_client[data["nodeId"]] = agent

        # TODO(later): For some reason the below code does not work because memory is not maintained
        # else:
        #     task.system_prompt = "You are a helpful chemical assistant who answers in concise but factual responses."
        #     task.user_prompt = (
        #         "Given the last computed reaction, answer the following query."
        #     )
        #     agent = self.retro_synth_context.node_id_to_charge_client[data["nodeId"]]
        #     agent.task = task

        async def run_and_report():
            result = await agent.run()
            await self.experiment.add_to_context(agent, task, result)
            # Report answer
            await self._send_processing_message(result, source="Agent")
            await self.websocket.send_json({"type": "complete"})

        asyncio.create_task(self.task_manager.run_task(run_and_report()))

    async def handle_get_username(self, _: dict) -> None:
        await self.websocket.send_json(
            {
                "type": "get-username-response",
                "username": self.username,
            }
        )
