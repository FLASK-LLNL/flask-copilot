from typing import Any, Optional
from fastapi import WebSocket
import asyncio
import os
from lc_conductor import ActionManager, TaskManager, CallbackLogger
from lc_conductor import (
    BuiltinToolDefinition,
    ToolRuntime,
    resolve_builtin_tool_descriptors,
)
from loguru import logger
from charge.experiments.experiment import Experiment
from charge.tasks.task import Task
from charge_backend.backend_helper_funcs import (
    CallbackHandler,
    PathwayStep,
    Reaction,
    ReactionAlternative,
    FlaskRunSettings,
)
from charge_backend.retrosynthesis.context import RetrosynthesisContext
from charge_backend.lmo.lmo_charge_backend_funcs import generate_lead_molecule
from charge_backend.charge_backend_custom import run_custom_problem
from functools import partial
from charge_backend.retrosynthesis.template import (
    template_based_retrosynthesis,
    compute_templates_for_node,
)
from charge_backend.retrosynthesis.ai import (
    ai_based_retrosynthesis,
    db_then_ai_retrosynthesis,
)
from charge_backend.retrosynthesis.alternatives import set_reaction_alternative
from charge_backend.prompt_debugger import debug_prompt
from charge_backend.backend_helper_funcs import Node
from charge_backend.attachments import image_refs, validate_image_attachments
from charge_backend.pdf import PdfDocumentRegistry


class FlaskActionManager(ActionManager):
    """Handles action state for a websocket connection."""

    def __init__(
        self,
        task_manager: TaskManager,
        experiment: Experiment,
        args,
        username: str,
        pdf_registry: Optional[PdfDocumentRegistry] = None,
        builtin_tool_definitions: Optional[list[BuiltinToolDefinition]] = None,
    ):
        super().__init__(
            task_manager,
            experiment,
            args,
            username,
            builtin_tool_definitions=builtin_tool_definitions,
        )
        self.run_settings: FlaskRunSettings = FlaskRunSettings()
        self.websocket = task_manager.websocket
        self.retro_synth_context: Optional[RetrosynthesisContext] = None
        self.pdf_registry = pdf_registry or PdfDocumentRegistry()

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

    def reset_problem_context(self, problem_type: Optional[str]) -> None:
        """Reset backend-only state that is scoped to the active problem type."""
        self.experiment.reset()
        if problem_type == "retrosynthesis":
            self.retro_synth_context = RetrosynthesisContext()
        else:
            self.retro_synth_context = None

    def selected_tool_runtime(self) -> ToolRuntime:
        runtime = super().selected_tool_runtime()
        if not self.pdf_registry.has_active_document(self.username):
            return runtime
        if any(tool.identifier == "consult_with_document" for tool in runtime.tools):
            return runtime
        return ToolRuntime(
            bearer_token=runtime.bearer_token,
            tools=[
                *runtime.tools,
                *resolve_builtin_tool_descriptors(
                    ["consult_with_document"],
                    self.builtin_tool_definitions,
                ),
            ],
        )

    def _document_reference_context(self) -> str:
        metadata = self.pdf_registry.active_metadata(self.username)
        if metadata is None:
            return ""
        return (
            "\n\nAn uploaded PDF reference is available through the "
            "`consult_with_document` tool. Use it when the user asks about the "
            f"reference document `{metadata.title}` ({metadata.name}); cite page "
            "numbers from the tool output."
        )

    def _with_document_reference_context(self, system_prompt: str) -> str:
        return f"{system_prompt}{self._document_reference_context()}"

    async def handle_load_state(self, data, *args, **kwargs) -> None:
        self.pdf_registry.clear(self.username)
        await super().handle_load_state(data, *args, **kwargs)

    async def handle_configure_pdf_reference(self, data: dict[str, Any]) -> None:
        reference = data.get("pdfReference")
        if reference is None:
            self.pdf_registry.clear(self.username)
            if data.get("silent"):
                return
            await self.websocket.send_json(
                {
                    "type": "pdf-reference-response",
                    "reference": None,
                }
            )
            return

        try:
            metadata = await asyncio.to_thread(
                self.pdf_registry.set_from_attachment,
                self.username,
                reference,
            )
        except ValueError as exc:
            await self.websocket.send_json(
                {
                    "type": "pdf-reference-response",
                    "reference": {
                        "id": (
                            str(reference.get("id") or "pdf_reference")
                            if isinstance(reference, dict)
                            else "pdf_reference"
                        ),
                        "name": (
                            str(reference.get("name") or "Reference PDF")
                            if isinstance(reference, dict)
                            else "Reference PDF"
                        ),
                        "mimeType": "application/pdf",
                        "sizeBytes": (
                            int(reference.get("sizeBytes") or 0)
                            if isinstance(reference, dict)
                            else 0
                        ),
                        "status": "error",
                        "error": str(exc),
                    },
                }
            )
            return
        await self.websocket.send_json(
            {
                "type": "pdf-reference-response",
                "reference": metadata.json(),
            }
        )

    async def log_progress(self, progress: str, agent_key: Optional[str] = None):
        logger.info(f"Reasoning: {progress}")
        await self._send_processing_message(
            progress,
            "Reasoning",
        )

    async def handle_compute(self, data: dict) -> None:
        self.setup_run_settings(data)
        problem_type = data.get("problemType")
        self.reset_problem_context(problem_type)
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
        tool_runtime = self.selected_tool_runtime()
        attachments = validate_image_attachments(data)
        if not data.get("_agentRequestAnnounced"):
            request_text = (
                f"LMO request: {data.get('query')}"
                if data.get("query")
                else f"LMO request for {data.get('smiles', 'current molecule')}"
            )
            await self._send_processing_message(
                request_text,
                source="LMO request",
                images=image_refs(attachments) if attachments else None,
            )

        run_func = partial(
            generate_lead_molecule,
            data["smiles"],
            self.experiment,
            self.args.json_file,
            self.args.max_retries,
            depth,
            tool_runtime,
            self.task_manager.websocket,
            self.run_settings,
            partial(self.log_progress, agent_key="lmo:main"),
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
            attachments,
            partial(
                self.send_agent_update,
                "lmo:main",
                debug=bool(data.get("debug")),
            ),
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
        tool_runtime = self.selected_tool_runtime()
        attachments = validate_image_attachments(data)
        await self._send_processing_message(
            data.get("userPrompt") or "Custom prompt",
            source="User",
            images=image_refs(attachments) if attachments else None,
        )

        run_func = partial(
            run_custom_problem,
            data["smiles"],
            self._with_document_reference_context(data["systemPrompt"]),
            data["userPrompt"],
            self.experiment,
            tool_runtime,
            self.task_manager.websocket,
            self.run_settings,
            partial(self.log_progress, agent_key="custom:main"),
            attachments,
            partial(
                self.send_agent_update,
                "custom:main",
                debug=bool(data.get("debug")),
            ),
        )

        await self.task_manager.run_task(run_func())

    async def handle_compute_reaction_from(self, data: dict) -> None:
        """Handle compute-reaction-from action."""
        self.setup_run_settings(data)
        attachments = validate_image_attachments(data)
        if self.retro_synth_context is None:
            raise ValueError("Retrosynthesis context not initialized")

        logger.info("Synthesize tree leaf action received")
        logger.info(f"Data: {data}")
        node = self.retro_synth_context.node_ids[data["nodeId"]]
        request_text = (
            f"Retrosynthesis request for node {data['nodeId']}: {data.get('query', '')}".strip()
            if data.get("query")
            else f"Retrosynthesis request for node {data['nodeId']}"
        )
        await self._send_processing_message(
            request_text,
            source="Retrosynthesis Request",
            images=image_refs(attachments) if attachments else None,
        )

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
                        "searching_0",
                        (
                            "Recomputing: searching database..."
                            if has_children
                            else "Searching database..."
                        ),
                        highlight="red",
                        label="Recomputing" if has_children else "Computing",
                    ).json(),
                },
            }
        )

        run_func = partial(
            (
                ai_based_retrosynthesis
                if data.get("aiOnly", False)
                else db_then_ai_retrosynthesis
            ),
            data["nodeId"],
            self.retro_synth_context,
            data.get("query", None),
            None,  # Unconstrained
            self.task_manager.websocket,
            self.experiment,
            self.args.config_file,
            self.run_settings,
            self.selected_tool_runtime(),
            partial(self.log_progress, agent_key=f"reaction:{data['nodeId']}"),
            attachments,
            partial(
                self.send_agent_update,
                f"reaction:{data['nodeId']}",
                debug=bool(data.get("debug")),
            ),
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
        attachments = validate_image_attachments(data)
        prompt = data.get("query")
        if prompt:
            await self._send_processing_message(
                f"Processing optimization query: {prompt} for node {data['nodeId']}",
                source="LMO Request",
                images=image_refs(attachments) if attachments else None,
            )
        else:
            await self._send_processing_message(
                f"Processing optimization refinement for node {data['nodeId']}",
                source="LMO Request",
            )
        data["_agentRequestAnnounced"] = True

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
        attachments = validate_image_attachments(data)
        if self.retro_synth_context is None:
            raise ValueError("Retrosynthesis context not initialized")

        # Get parent
        if data["nodeId"] not in self.retro_synth_context.parents:
            await self._send_processing_message(
                f"Cannot find parent reaction for node {data['nodeId']}", source="Agent"
            )
            await self.websocket.send_json({"type": "complete"})
            return

        parent_nodeid = self.retro_synth_context.parents[data["nodeId"]]
        parent_node = self.retro_synth_context.node_ids[parent_nodeid]
        smiles = self.retro_synth_context.node_ids[data["nodeId"]].smiles
        request_text = (
            f"Retrosynthesis request for node {data['nodeId']}: {data.get('query', '')}".strip()
            if data.get("query")
            else f"Retrosynthesis request for node {data['nodeId']}"
        )
        await self._send_processing_message(
            request_text,
            source="Retrosynthesis Request",
            images=image_refs(attachments) if attachments else None,
        )

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
            self.selected_tool_runtime(),
            partial(self.log_progress, agent_key=f"reaction:{parent_nodeid}"),
            attachments,
            partial(
                self.send_agent_update,
                f"reaction:{parent_nodeid}",
                debug=bool(data.get("debug")),
            ),
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
        attachments = validate_image_attachments(data)
        await self._send_processing_message(
            f"Processing molecule query: {data['query']} for node {data['nodeId']}",
            source="User",
            images=image_refs(attachments) if attachments else None,
        )
        smiles = data["smiles"]
        tool_runtime = self.selected_tool_runtime()

        task = Task(
            system_prompt=self._with_document_reference_context(
                f"You are a helpful chemical assistant who answers in concise but factual responses. Answer the following query about the molecule given by the SMILES string `{smiles}`."
            ),
            user_prompt=data["query"],
            attachments=attachments,
            **tool_runtime.task_kwargs(),
        )

        # Use the full experiment state
        callback_handler = CallbackHandler(
            self.websocket,
            agent_key=f"molecule:{data['nodeId']}",
            on_agent_update=partial(
                self.send_agent_update,
                f"molecule:{data['nodeId']}",
                debug=bool(data.get("debug")),
            ),
        )
        agent = self.experiment.create_agent_with_experiment_state(
            task=task,
            agent_key=f"molecule:{data['nodeId']}",
            callback=callback_handler,
        )

        async def run_and_report():
            if self.run_settings.prompt_debugging:
                await debug_prompt(agent, self.websocket)
            result = await agent.run(
                partial(self.log_progress, agent_key=f"molecule:{data['nodeId']}")
            )
            await callback_handler.drain()
            self.experiment.add_to_context(agent, task, result)
            await self.send_agent_update(
                f"molecule:{data['nodeId']}", debug=bool(data.get("debug"))
            )
            # Report answer
            await self._send_processing_message(
                result,
                source="Agent",
            )
            await self.websocket.send_json({"type": "complete"})

        asyncio.create_task(self.task_manager.run_task(run_and_report()))

    async def handle_custom_query_reaction(self, data: dict) -> None:
        """Handle a query on the reaction (from nodeId to its reactants)."""
        self.setup_run_settings(data)
        assert self.retro_synth_context is not None
        attachments = validate_image_attachments(data)

        await self._send_processing_message(
            f"Processing reaction query: {data['query']} for node {data['nodeId']}",
            source="User",
            images=image_refs(attachments) if attachments else None,
        )

        node = self.retro_synth_context.node_ids[data["nodeId"]]
        tool_runtime = self.selected_tool_runtime()

        task = Task(
            system_prompt=self._with_document_reference_context(""),
            user_prompt=data["query"],
            attachments=attachments,
            **tool_runtime.task_kwargs(),
        )

        reaction_str = self._reaction_context_for_node(str(data["nodeId"]), data)

        task.system_prompt = f"You are a helpful chemical assistant who answers in concise but factual responses. Given the following reaction (as SMILES strings):\n{reaction_str}\n\nAnswer the following query."

        callback_handler = CallbackHandler(
            self.websocket,
            agent_key=f"reaction:{data['nodeId']}",
            on_agent_update=partial(
                self.send_agent_update,
                f"reaction:{data['nodeId']}",
                debug=bool(data.get("debug")),
            ),
        )
        agent = self.experiment.create_agent_with_experiment_state(
            task=task,
            agent_key=f"reaction:{data['nodeId']}",
            callback=callback_handler,
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
            result = await agent.run(
                partial(self.log_progress, agent_key=f"reaction:{data['nodeId']}")
            )
            await callback_handler.drain()
            self.experiment.add_to_context(agent, task, result)
            await self.send_agent_update(
                f"reaction:{data['nodeId']}", debug=bool(data.get("debug"))
            )
            # Report answer
            await self._send_processing_message(
                result,
                source="Agent",
            )
            await self.websocket.send_json({"type": "complete"})

        asyncio.create_task(self.task_manager.run_task(run_and_report()))

    def _agent_key_parts(self, agent_key: str) -> tuple[str, str]:
        if ":" not in agent_key:
            return agent_key, ""
        kind, target = agent_key.split(":", 1)
        return kind, target

    def _reaction_hover_info_for_node(
        self, node_id: str, data: Optional[dict[str, Any]] = None
    ) -> str:
        if (
            self.retro_synth_context is not None
            and node_id in self.retro_synth_context.node_ids
        ):
            reaction = self.retro_synth_context.node_ids[node_id].reaction
            hover_info = getattr(reaction, "hoverInfo", None)
            if isinstance(hover_info, str) and hover_info.strip():
                return hover_info.strip()

        metadata = data.get("metadata") if isinstance(data, dict) else None
        if isinstance(metadata, dict):
            for key in ("reactionHoverInfo", "hoverInfo"):
                hover_info = metadata.get(key)
                if isinstance(hover_info, str) and hover_info.strip():
                    return hover_info.strip()
        return ""

    def _reaction_context_for_node(
        self, node_id: str, data: Optional[dict[str, Any]] = None
    ) -> str:
        if (
            self.retro_synth_context is None
            or node_id not in self.retro_synth_context.node_ids
        ):
            hover_info = self._reaction_hover_info_for_node(node_id, data)
            return f"Reaction hover information:\n{hover_info}" if hover_info else ""
        node = self.retro_synth_context.node_ids[node_id]
        child_nodes = [
            nid
            for nid, parent in self.retro_synth_context.parents.items()
            if parent == node_id
        ]
        reactants = [self.retro_synth_context.node_ids[nid] for nid in child_nodes]
        reactants_str = "\n".join(reactant.smiles for reactant in reactants)
        reaction_str = f"Product: {node.smiles}\nReactants:\n{reactants_str}"
        if node_id in self.retro_synth_context.node_id_to_reasoning_summary:
            reaction_str += (
                "\nAdditionally, the following context is given: "
                f"{self.retro_synth_context.node_id_to_reasoning_summary[node_id]}"
            )
        hover_info = self._reaction_hover_info_for_node(node_id, data)
        if hover_info:
            reaction_str += f"\n\nReaction hover information:\n{hover_info}"
        return reaction_str

    def _with_reaction_hover_info(
        self, system_prompt: str, node_id: str, data: Optional[dict[str, Any]] = None
    ) -> str:
        hover_info = self._reaction_hover_info_for_node(node_id, data)
        if not hover_info or hover_info in system_prompt:
            return system_prompt
        return f"{system_prompt}\n\nReaction hover information:\n{hover_info}"

    def _chat_task_for_agent(
        self,
        agent_key: str,
        data: dict[str, Any],
        attachments: list[dict[str, object]],
    ) -> Task:
        kind, target = self._agent_key_parts(agent_key)
        existing = self.experiment.agent_registry.get(agent_key)
        existing_agent = existing.agent if existing else None
        existing_task = getattr(existing_agent, "task", None)
        system_prompt = (
            existing_task.get_system_prompt()
            if existing_task is not None and hasattr(existing_task, "get_system_prompt")
            else ""
        )

        if not system_prompt:
            if kind == "molecule":
                request_metadata = (
                    data.get("metadata")
                    if isinstance(data.get("metadata"), dict)
                    else {}
                )
                smiles = data.get("smiles") or request_metadata.get("smiles") or target
                system_prompt = (
                    "You are a helpful chemical assistant who answers in concise but "
                    f"factual responses. Answer questions about the molecule `{smiles}`."
                )
            elif kind == "reaction":
                node_id = data.get("nodeId") or target
                reaction_context = self._reaction_context_for_node(str(node_id), data)
                system_prompt = (
                    "You are a helpful chemical assistant who answers in concise but "
                    "factual responses. Given the following reaction as SMILES strings:\n"
                    f"{reaction_context}\n\nAnswer the user's questions."
                )
            else:
                system_prompt = (
                    "You are a helpful chemical assistant who answers in concise but "
                    "factual responses."
                )

            system_prompt = self._with_document_reference_context(system_prompt)
        elif kind == "reaction":
            node_id = str(data.get("nodeId") or target)
            system_prompt = self._with_reaction_hover_info(system_prompt, node_id, data)

        return Task(
            system_prompt=system_prompt,
            user_prompt=str(data.get("query") or ""),
            attachments=attachments,
            **self.selected_tool_runtime().task_kwargs(),
        )

    async def handle_chat_agent(self, data: dict[str, Any]) -> None:
        self.setup_run_settings(data)
        agent_key = str(data.get("agentKey") or "")
        if not agent_key:
            raise ValueError("agentKey is required")
        query = str(data.get("query") or "").strip()
        if not query:
            raise ValueError("query is required")

        attachments = validate_image_attachments(data)
        task = self._chat_task_for_agent(agent_key, data, attachments)
        await self._send_processing_message(
            query,
            source="User",
            images=image_refs(attachments) if attachments else None,
        )
        callback_handler = CallbackHandler(
            self.websocket,
            agent_key=agent_key,
            on_agent_update=partial(
                self.send_agent_update,
                agent_key,
                debug=bool(data.get("debug")),
            ),
        )
        agent = self.experiment.create_agent_with_experiment_state(
            task=task,
            agent_key=agent_key,
            callback=callback_handler,
        )
        kind, target = self._agent_key_parts(agent_key)
        if kind == "reaction" and self.retro_synth_context is not None:
            node_id = str(data.get("nodeId") or target)
            self.retro_synth_context.node_id_to_charge_client[node_id] = agent

        async def run_and_report():
            if self.run_settings.prompt_debugging:
                await debug_prompt(agent, self.websocket)
            result = await agent.run(partial(self.log_progress, agent_key=agent_key))
            await callback_handler.drain()
            self.experiment.add_to_context(agent, task, result)
            await self.send_agent_update(agent_key, debug=bool(data.get("debug")))
            await self._send_processing_message(
                result,
                source="Agent",
            )
            await self.websocket.send_json({"type": "complete"})

        asyncio.create_task(self.task_manager.run_task(run_and_report()))

    async def handle_get_username(self, _: dict) -> None:
        await self.websocket.send_json(
            {
                "type": "get-username-response",
                "username": self.username,
            }
        )

    async def handle_load_state(self, data: dict, *args, **kwargs) -> None:
        await super().handle_load_state(data, *args, **kwargs)

        problem_type = data.get("problemType")
        if problem_type == "retrosynthesis":
            self.retro_synth_context = RetrosynthesisContext()
            await self.retro_synth_context.load_state(data)
        else:
            self.retro_synth_context = None
