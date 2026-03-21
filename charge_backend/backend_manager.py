from typing import Any, Callable, Literal, Optional, Tuple
from fastapi import WebSocket
import asyncio
import os
from lc_conductor import ActionManager, TaskManager, CallbackLogger
from loguru import logger
from concurrent.futures import ProcessPoolExecutor
from charge.experiments.experiment import Experiment
from charge.tasks.task import Task
from backend_helper_funcs import (
    CallbackHandler,
    PathwayStep,
    Reaction,
    ReactionAlternative,
    FlaskRunSettings,
)
from retrosynthesis.context import RetrosynthesisContext
from lmo.lmo_charge_backend_funcs import generate_lead_molecule
from charge_backend_custom import run_custom_problem
from functools import partial
from lc_conductor.tool_registration import (
    ToolList,
    list_server_urls,
    list_server_tools,
)
from retrosynthesis.template import (
    template_based_retrosynthesis,
    compute_templates_for_node,
)
from retrosynthesis.ai import ai_based_retrosynthesis
from retrosynthesis.alternatives import set_reaction_alternative
from charge_backend.prompt_debugger import debug_prompt
from backend_helper_funcs import Node
from builtin_tools import BuiltinToolDefinition, resolve_builtin_tools


class FlaskActionManager(ActionManager):
    """Handles action state for a websocket connection."""

    def __init__(
        self,
        task_manager: TaskManager,
        experiment: Experiment,
        args,
        username: str,
        builtin_tool_definitions: Optional[list[BuiltinToolDefinition]] = None,
    ):
        super().__init__(task_manager, experiment, args, username)
        self.run_settings: FlaskRunSettings = FlaskRunSettings()
        self.websocket = task_manager.websocket
        self.builtin_tool_definitions = builtin_tool_definitions or []
        self.task_manager.available_builtin_tool_ids = None

    def _selected_mcp_tools(self) -> list[str]:
        if self.task_manager.available_tools is None:
            return list_server_urls()
        return self.task_manager.available_tools

    def _selected_builtin_tools(self) -> list[Callable[..., Any]]:
        return resolve_builtin_tools(
            self.task_manager.available_builtin_tool_ids,
            self.builtin_tool_definitions,
        )

    def setup_retro_synth_context(self) -> None:
        if self.retro_synth_context is None:
            self.retro_synth_context = RetrosynthesisContext()

    def get_retro_synth_context(self) -> RetrosynthesisContext:
        self.setup_retro_synth_context()
        assert self.retro_synth_context is not None
        return self.retro_synth_context

    def setup_run_settings(self, data: dict[str, Any]):
        if "runSettings" in data:
            self.run_settings = FlaskRunSettings(**data["runSettings"])

    async def log_progress(self, progress: str):
        logger.info(f"Reasoning: {progress}")
        await self._send_processing_message(progress, "Reasoning")

    async def handle_compute(self, data: dict) -> None:
        self.experiment.reset()
        self.setup_run_settings(data)
        problem_type = data.get("problemType")
        if problem_type == "optimization":
            asyncio.create_task(self._handle_optimization(data))
        elif problem_type == "retrosynthesis":
            asyncio.create_task(self._handle_retrosynthesis(data))
        elif problem_type == "custom":
            asyncio.create_task(self._handle_custom_problem(data))
        else:
            raise ValueError(f"Unknown problem type: {problem_type}")

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
        number_of_molecules = customization.get("numberOfMolecules", 10)
        num_top_candidates = customization.get("numTopCandidates", 3)
        depth = customization.get("depth", 3)

        run_func = partial(
            generate_lead_molecule,
            data["smiles"],
            self.experiment,
            self.args.json_file,
            self.args.max_retries,
            depth,
            self._selected_mcp_tools(),
            self._selected_builtin_tools(),
            self.task_manager.websocket,
            self.run_settings,
            self.log_progress,
            *property_attributes,
            data.get("query", None),
            initial_level,
            initial_node_id,
            initial_x_position,
            enable_constraints,
            molecular_similarity,
            diversity_penalty,
            exploration_rate,
            additional_constraints,
            number_of_molecules,
            num_top_candidates,
        )
        await self.task_manager.run_task(run_func())

    async def _handle_retrosynthesis(self, data: dict) -> None:
        """Handle retrosynthesis problem type."""
        run_func = partial(
            template_based_retrosynthesis,
            data["smiles"],
            self.args.config_file,
            self.get_retro_synth_context(),
            self.task_manager.websocket,
            self.run_settings,
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
            self._selected_mcp_tools(),
            self._selected_builtin_tools(),
            self.task_manager.websocket,
            self.run_settings,
            self.log_progress,
        )

        await self.task_manager.run_task(run_func())

    async def handle_compute_reaction_from(self, data: dict) -> None:
        """Handle compute-reaction-from action."""
        self.setup_run_settings(data)
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
            self.run_settings,
            self._selected_mcp_tools(),
            self._selected_builtin_tools(),
            self.log_progress,
        )

        asyncio.create_task(self.task_manager.run_task(run_func()))

    async def handle_template_retrosynthesis(self, data: dict) -> None:
        """Handle compute-reaction-templates action."""
        self.setup_run_settings(data)
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
            self.run_settings,
        )

        asyncio.create_task(self.task_manager.run_task(run_func()))

    async def handle_optimize_from(self, data: dict) -> None:
        """Handle optimize-from action."""
        self.setup_run_settings(data)
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
        self.setup_run_settings(data)
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
            self.run_settings,
            self._selected_mcp_tools(),
            self._selected_builtin_tools(),
            self.log_progress,
        )

        asyncio.create_task(self.task_manager.run_task(run_func()))

    async def handle_set_reaction_alternative(self, data: dict) -> None:
        """Handle set-reaction-alternative action."""
        self.setup_run_settings(data)
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

    async def handle_orchestrator_settings_update(self, data: dict) -> None:
        if "moleculeName" in data:
            self.run_settings.molecule_name_format = data["moleculeName"]
        await ActionManager.handle_orchestrator_settings_update(self, data)

    async def handle_custom_query_molecule(self, data: dict) -> None:
        """Handle a query on the given molecule."""
        self.setup_run_settings(data)
        await self._send_processing_message(
            f"Processing molecule query: {data['query']} for node {data['nodeId']}"
        )
        smiles = data["smiles"]

        task = Task(
            system_prompt=f"You are a helpful chemical assistant who answers in concise but factual responses. Answer the following query about the molecule given by the SMILES string `{smiles}`.",
            user_prompt=data["query"],
            server_urls=self._selected_mcp_tools(),
            builtin_tools=self._selected_builtin_tools(),
        )

        # Use the full experiment state
        agent = self.experiment.create_agent_with_experiment_state(
            task=task,
            callback=CallbackHandler(self.websocket),
        )

        async def run_and_report():
            if self.run_settings.prompt_debugging:
                await debug_prompt(agent, self.websocket)
            result = await agent.run(self.log_progress)
            self.experiment.add_to_context(agent, task, result)
            # Report answer
            await self._send_processing_message(result, source="Agent")
            await self.websocket.send_json({"type": "complete"})

        asyncio.create_task(self.task_manager.run_task(run_and_report()))

    async def handle_custom_query_reaction(self, data: dict) -> None:
        """Handle a query on the reaction (from nodeId to its reactants)."""
        self.setup_run_settings(data)
        assert self.retro_synth_context is not None

        await self._send_processing_message(
            f"Processing reaction query: {data['query']} for node {data['nodeId']}"
        )

        node = self.retro_synth_context.node_ids[data["nodeId"]]

        task = Task(
            system_prompt="",
            user_prompt=data["query"],
            server_urls=self._selected_mcp_tools(),
            builtin_tools=self._selected_builtin_tools(),
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
            if self.run_settings.prompt_debugging:
                await debug_prompt(agent, self.websocket)
            result = await agent.run(self.log_progress)
            self.experiment.add_to_context(agent, task, result)
            # Report answer
            await self._send_processing_message(result, source="Agent")
            await self.websocket.send_json({"type": "complete"})

        asyncio.create_task(self.task_manager.run_task(run_and_report()))

    async def handle_list_tools(self, *args, **kwargs) -> None:
        tools = []
        server_list = list_server_urls()
        for server in server_list:
            tool_list = await list_server_tools([server])
            tool_names = [name for name, _ in tool_list]
            tools.append(
                ToolList(
                    server=server,
                    names=tool_names,
                    kind="mcp",
                    identifier=server,
                )
            )

        builtin_tools = [
            tool_definition.to_client_tool()
            for tool_definition in self.builtin_tool_definitions
        ]

        await self.websocket.send_json(
            {
                "type": "available-tools-response",
                "tools": [tool.json() for tool in tools] + builtin_tools,
            }
        )

    async def handle_select_tools_for_task(self, data: dict) -> None:
        """Handle select-tools-for-task action."""
        logger.info("Select tools for task")
        logger.info(f"Data: {data}")

        selected_mcp_tools: list[str] = []
        selected_builtin_tool_ids: list[str] = []

        for server in data.get("enabledTools", {}).get("selectedTools", []):
            tool_server = server.get("tool_server", {})
            if tool_server.get("kind") == "builtin":
                identifier = tool_server.get("identifier")
                if identifier:
                    selected_builtin_tool_ids.append(identifier)
                continue

            server_url = tool_server.get("server")
            if server_url:
                selected_mcp_tools.append(server_url)

        self.task_manager.available_tools = selected_mcp_tools
        self.task_manager.available_builtin_tool_ids = selected_builtin_tool_ids

    async def handle_get_username(self, _: dict) -> None:
        await self.websocket.send_json(
            {
                "type": "get-username-response",
                "username": self.username,
            }
        )

    async def restore_retrosynth_context(self, ctxt: RetrosynthesisContext, data: dict):
        node_data = data.get("nodes", [])

        # Deserialize each node.
        #
        # Nodes are complex objects, and there doesn't seem to be a
        # built-in way to deserialize dict->(object of some dataclass)
        # automatically. The simplistic dict-expansion approach just
        # writes dicts to subobjects even if they're annotated as some
        # complex type. There seem to be robust library solutions that
        # properly handle subobjects ("dacite" comes up in many search
        # results), but I'm hesitant to pull in more dependencies for
        # a one-off task. So we "brute-force" it here, but it would be
        # good to keep an eye for other such situations and pursue
        # better solutions accordingly.
        for n in node_data:
            n["yield_"] = n["yield"]
            del n["yield"]

            # Deserialize the reaction (optional)
            if "reaction" in n:
                reaction_dict = n["reaction"]

                # Deserialize each reaction alternative (optional)
                if "alternatives" in reaction_dict:
                    reaction_alternatives = reaction_dict["alternatives"]

                    # Deserialize all the pathways (not optional)
                    for alternative_dict in reaction_alternatives:
                        alternative_dict["pathway"] = [
                            PathwayStep(**pws) for pws in alternative_dict["pathway"]
                        ]

                    reaction_dict["alternatives"] = [
                        ReactionAlternative(**ra) for ra in reaction_alternatives
                    ]

                n["reaction"] = Reaction(**reaction_dict)

            await ctxt.add_node(Node(**n))

    async def handle_load_state(self, data: dict, *args, **kwargs) -> None:
        await super().handle_load_state(data, *args, **kwargs)

        problem_type = data.get("problemType")
        if problem_type == "retrosynthesis":
            if not self.retro_synth_context:
                await self.restore_retrosynth_context(
                    self.get_retro_synth_context(), data
                )
