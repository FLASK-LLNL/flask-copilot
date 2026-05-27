from typing import Any, Optional
from fastapi import WebSocket
import asyncio
import os
import hashlib
import json
import math
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

    async def log_progress(self, progress: str):
        logger.info(f"Reasoning: {progress}")
        await self._send_processing_message(progress, "Reasoning")

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
            attachments,
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
        if attachments:
            await self._send_processing_message(
                "Processing custom prompt with attached images",
                source="User",
                images=image_refs(attachments),
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
            self.log_progress,
            attachments,
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
        if data.get("query") or attachments:
            await self._send_processing_message(
                f"Retrosynthesis request for node {data['nodeId']}: {data.get('query', '')}".strip(),
                source="Retrosynthesis Request",
                images=image_refs(attachments),
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
            self.log_progress,
            attachments,
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
        attachments = validate_image_attachments(data)
        if self.retro_synth_context is None:
            raise ValueError("Retrosynthesis context not initialized")

        if data.get("query") or attachments:
            await self._send_processing_message(
                f"Retrosynthesis request for node {data['nodeId']}: {data.get('query', '')}".strip(),
                source="Retrosynthesis Request",
                images=image_refs(attachments) if attachments else None,
            )

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
            self.log_progress,
            attachments,
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
        callback_handler = CallbackHandler(self.websocket)
        agent = self.experiment.create_agent_with_experiment_state(
            task=task,
            agent_key=f"molecule:{data['nodeId']}",
            agent_metadata={
                "kind": "molecule",
                "nodeId": data["nodeId"],
                "target": data["nodeId"],
                "smiles": smiles,
                "title": f"Molecule {data['nodeId']}",
                "subtitle": smiles,
            },
            callback=callback_handler,
        )

        async def run_and_report():
            if self.run_settings.prompt_debugging:
                await debug_prompt(agent, self.websocket)
            result = await agent.run(self.log_progress)
            await callback_handler.drain()
            self._record_latest_user_message_metadata(
                f"molecule:{data['nodeId']}", task
            )
            self.experiment.add_to_context(agent, task, result)
            # Report answer
            await self._send_processing_message(result, source="Agent")
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

        callback_handler = CallbackHandler(self.websocket)
        agent = self.experiment.create_agent_with_experiment_state(
            task=task,
            agent_key=f"reaction:{data['nodeId']}",
            agent_metadata={
                "kind": "reaction",
                "nodeId": data["nodeId"],
                "target": data["nodeId"],
                "smiles": node.smiles,
                "title": f"Reaction {data['nodeId']}",
                "subtitle": node.smiles,
            },
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
            result = await agent.run(self.log_progress)
            await callback_handler.drain()
            self._record_latest_user_message_metadata(
                f"reaction:{data['nodeId']}", task
            )
            self.experiment.add_to_context(agent, task, result)
            # Report answer
            await self._send_processing_message(result, source="Agent")
            await self.websocket.send_json({"type": "complete"})

        asyncio.create_task(self.task_manager.run_task(run_and_report()))

    def _agent_key_parts(self, agent_key: str) -> tuple[str, str]:
        if ":" not in agent_key:
            return agent_key, ""
        kind, target = agent_key.split(":", 1)
        return kind, target

    def _default_agent_metadata(
        self, agent_key: str, data: Optional[dict[str, Any]] = None
    ) -> dict[str, Any]:
        data = data or {}
        kind, target = self._agent_key_parts(agent_key)
        metadata = {
            "kind": kind,
            "target": target,
            **(data.get("metadata") if isinstance(data.get("metadata"), dict) else {}),
        }
        if data.get("nodeId") is not None:
            metadata["nodeId"] = data["nodeId"]
        if data.get("smiles") is not None:
            metadata["smiles"] = data["smiles"]
        if "title" not in metadata:
            if kind == "reaction":
                metadata["title"] = f"Reaction {target}"
            elif kind == "molecule":
                metadata["title"] = f"Molecule {target}"
            elif kind == "custom":
                metadata["title"] = "Custom prompt"
            elif kind == "lmo":
                metadata["title"] = "Lead molecule optimization"
            else:
                metadata["title"] = agent_key
        return metadata

    def _agent_session_records(self) -> dict[str, dict[str, Any]]:
        state = self.experiment.save_state()
        records = state.get("agentSessions", {})
        return records if isinstance(records, dict) else {}

    def _session_messages_from_record(
        self, record: dict[str, Any]
    ) -> tuple[list[dict[str, Any]], Optional[dict[str, Any]]]:
        memory = record.get("memory")
        if not memory:
            return [], None
        try:
            session = json.loads(memory) if isinstance(memory, str) else memory
        except json.JSONDecodeError:
            return [], None
        if not isinstance(session, dict):
            return [], None
        messages = session.get("state", {}).get("in_memory", {}).get("messages", [])
        return (messages if isinstance(messages, list) else []), session

    def _prompt_context_from_task(
        self, task: Optional[dict[str, Any]]
    ) -> list[dict[str, str]]:
        if not isinstance(task, dict):
            return []
        context: list[dict[str, str]] = []
        system_prompt = task.get("system_prompt")
        if isinstance(system_prompt, str) and system_prompt.strip():
            context.append({"title": "Instructions", "text": system_prompt})
        user_prompt = task.get("user_prompt")
        if isinstance(user_prompt, str) and user_prompt.strip():
            context.append({"title": "Task prompt", "text": user_prompt})
        return context

    def _prompt_context_from_record(
        self, record: dict[str, Any]
    ) -> list[dict[str, str]]:
        return self._prompt_context_from_task(record.get("task"))

    @staticmethod
    def _token_count(text: str, model: object = None) -> int:
        if not text:
            return 0
        return max(1, math.ceil(len(text) / 4))

    def _content_token_count(self, content: object, model: object = None) -> int:
        if content is None:
            return 0
        if not isinstance(content, dict):
            return self._token_count(str(content), model)

        text = content.get("text")
        if isinstance(text, str):
            return self._token_count(text, model)

        uri = content.get("uri") or content.get("dataUrl")
        if isinstance(uri, str) and uri.startswith("data:image/"):
            return 85

        for key in ("arguments", "result", "output"):
            value = content.get(key)
            if value is not None:
                return self._token_count(json.dumps(value, default=str), model)

        return self._token_count(json.dumps(content, default=str), model)

    def _message_token_count(
        self, raw_message: dict[str, Any], model: object = None
    ) -> int:
        token_count = 4 + self._token_count(
            self._normalized_message_role(raw_message), model
        )
        contents = raw_message.get("contents", [])
        if not isinstance(contents, list):
            contents = [contents]
        for content in contents:
            token_count += self._content_token_count(content, model)
        return token_count

    def _context_usage_from_model_info(
        self, model_info: dict[str, Any]
    ) -> Optional[dict[str, Any]]:
        raw_usage = model_info.get("lastUsage")
        if not isinstance(raw_usage, dict):
            return None

        input_tokens = raw_usage.get("inputTokens")
        output_tokens = raw_usage.get("outputTokens")
        total_tokens = raw_usage.get("totalTokens")
        used_tokens = total_tokens
        if not isinstance(used_tokens, int):
            return None

        model = model_info.get("model")
        usage: dict[str, Any] = {
            "usedTokens": used_tokens,
            "estimated": False,
        }
        if isinstance(input_tokens, int):
            usage["inputTokens"] = input_tokens
        if isinstance(output_tokens, int):
            usage["outputTokens"] = output_tokens
        if isinstance(total_tokens, int):
            usage["totalTokens"] = total_tokens
        if isinstance(model, str) and model.strip():
            usage["model"] = model
        return usage

    def _estimate_context_usage(
        self,
        record: dict[str, Any],
        raw_messages: list[dict[str, Any]],
    ) -> dict[str, Any]:
        model_info = record.get("modelInfo") if isinstance(record, dict) else {}
        if not isinstance(model_info, dict):
            model_info = {}
        model = model_info.get("model")
        task = record.get("task") if isinstance(record.get("task"), dict) else {}

        used_tokens = 0
        if isinstance(task, dict):
            system_prompt = task.get("system_prompt")
            if isinstance(system_prompt, str):
                used_tokens += self._token_count(system_prompt, model)
            if not raw_messages:
                user_prompt = task.get("user_prompt")
                if isinstance(user_prompt, str):
                    used_tokens += self._token_count(user_prompt, model)

        for raw_message in raw_messages:
            if isinstance(raw_message, dict):
                used_tokens += self._message_token_count(raw_message, model)

        usage: dict[str, Any] = {
            "usedTokens": used_tokens,
            "estimated": True,
        }
        if isinstance(model, str) and model.strip():
            usage["model"] = model
        return usage

    def _context_usage_for_history(
        self,
        record: dict[str, Any],
        raw_messages: list[dict[str, Any]],
    ) -> dict[str, Any]:
        model_info = record.get("modelInfo")
        if isinstance(model_info, dict):
            actual_usage = self._context_usage_from_model_info(model_info)
            if actual_usage is not None:
                return actual_usage
        return self._estimate_context_usage(record, raw_messages)

    @staticmethod
    def _message_label_for_task(
        task: Optional[dict[str, Any]],
        *,
        agent_key: str,
        metadata: dict[str, Any],
    ) -> Optional[str]:
        if not isinstance(task, dict):
            return None
        class_name = str(task.get("class_name") or "")
        module = str(task.get("module") or "")
        kind = str(metadata.get("kind") or agent_key.split(":", 1)[0])
        if "Retrosynthesis" in class_name or "retrosynthesis" in module:
            return "Retrosynthesis request"
        if class_name == "LMOTask" or ".lmo." in module:
            return "LMO request"
        if kind == "custom":
            return "Custom prompt"
        return None

    @staticmethod
    def _message_metadata_from_record(record: dict[str, Any]) -> dict[str, Any]:
        metadata = record.get("messageMetadata")
        return metadata if isinstance(metadata, dict) else {}

    @staticmethod
    def _normalized_message_role(raw_message: dict[str, Any]) -> str:
        role = str(raw_message.get("role") or raw_message.get("type") or "assistant")
        if role == "message":
            role = "assistant"
        return role if role in {"user", "assistant", "system", "tool"} else "assistant"

    def _record_latest_user_message_metadata(
        self,
        agent_key: str,
        task: Task,
        *,
        label: Optional[str] = None,
    ) -> None:
        registry_item = self.experiment.agent_registry.get(agent_key)
        if not registry_item:
            return
        agent = registry_item.get("agent")
        if agent is None or not hasattr(agent, "save_memory"):
            return
        memory = agent.save_memory()
        raw_messages, _ = self._session_messages_from_record({"memory": memory})
        latest_user_index = next(
            (
                index
                for index, raw_message in reversed(list(enumerate(raw_messages)))
                if isinstance(raw_message, dict)
                and self._normalized_message_role(raw_message) == "user"
            ),
            None,
        )
        if latest_user_index is None:
            return

        task_json = task.to_json() if hasattr(task, "to_json") else {}
        prompt_context = self._prompt_context_from_task(task_json)
        metadata: dict[str, Any] = {}
        if prompt_context:
            metadata["promptContext"] = prompt_context
        message_label = label or self._message_label_for_task(
            task_json,
            agent_key=agent_key,
            metadata=registry_item.get("metadata", {}),
        )
        if message_label:
            metadata["label"] = message_label
        if not metadata:
            return

        message_metadata = {
            **(
                registry_item.get("messageMetadata")
                if isinstance(registry_item.get("messageMetadata"), dict)
                else {}
            )
        }
        message_metadata[str(latest_user_index)] = metadata
        registry_item["messageMetadata"] = message_metadata

    @staticmethod
    def _image_ref_from_data_url(
        data_url: str, content: dict[str, Any]
    ) -> dict[str, Any]:
        image_id = content.get("id")
        if not image_id:
            image_id = (
                "session-image-"
                + hashlib.sha256(data_url.encode("utf-8")).hexdigest()[:16]
            )
        mime_type = (
            content.get("media_type")
            or content.get("mimeType")
            or data_url.split(";", 1)[0].removeprefix("data:")
            or "image"
        )
        return {
            "id": image_id,
            "name": content.get("name") or "Uploaded image",
            "mimeType": mime_type,
            "sizeBytes": int(content.get("sizeBytes") or 0),
            "dataUrl": data_url,
        }

    def _normalize_agent_history(
        self,
        agent_key: str,
        record: dict[str, Any],
        *,
        debug: bool = False,
    ) -> dict[str, Any]:
        metadata = {
            **self._default_agent_metadata(agent_key),
            **(
                record.get("metadata")
                if isinstance(record.get("metadata"), dict)
                else {}
            ),
        }
        raw_messages, raw_session = self._session_messages_from_record(record)
        messages: list[dict[str, Any]] = []
        prompt_context = self._prompt_context_from_record(record)
        message_metadata = self._message_metadata_from_record(record)
        user_message_indexes = [
            index
            for index, raw_message in enumerate(raw_messages)
            if isinstance(raw_message, dict)
            and self._normalized_message_role(raw_message) == "user"
        ]

        for index, raw_message in enumerate(raw_messages):
            if not isinstance(raw_message, dict):
                continue
            role = self._normalized_message_role(raw_message)
            text_parts: list[str] = []
            reasoning_parts: list[dict[str, Any]] = []
            tool_events: list[dict[str, Any]] = []
            images: list[dict[str, Any]] = []
            contents = raw_message.get("contents", [])
            if not isinstance(contents, list):
                contents = [contents]

            for content in contents:
                if not isinstance(content, dict):
                    text_parts.append(str(content))
                    continue
                content_type = str(content.get("type") or "")
                text = content.get("text")
                uri = content.get("uri") or content.get("dataUrl")
                if content_type == "text" and isinstance(text, str):
                    text_parts.append(text)
                elif "reasoning" in content_type:
                    reasoning_item = {
                        "type": content_type or "reasoning",
                        "text": text if isinstance(text, str) else "",
                    }
                    if debug:
                        reasoning_item["debug"] = content.get(
                            "additional_properties", {}
                        )
                    reasoning_parts.append(reasoning_item)
                elif isinstance(uri, str) and uri.startswith("data:image/"):
                    images.append(self._image_ref_from_data_url(uri, content))
                elif content_type in {
                    "function_call",
                    "function_result",
                    "mcp_server_tool_call",
                    "mcp_server_tool_result",
                }:
                    tool_event = {
                        "type": content_type,
                        "name": content.get("name") or content.get("tool_name"),
                        "text": json.dumps(content, indent=2, default=str),
                    }
                    if debug:
                        tool_event["raw"] = content
                    tool_events.append(tool_event)
                elif isinstance(text, str):
                    text_parts.append(text)
                elif debug:
                    tool_events.append(
                        {
                            "type": content_type or "raw",
                            "text": json.dumps(content, indent=2, default=str),
                            "raw": content,
                        }
                    )

            message: dict[str, Any] = {
                "id": f"{agent_key}:{index}",
                "role": role,
                "text": "\n\n".join(part for part in text_parts if part),
                "images": images,
                "reasoning": reasoning_parts,
                "toolEvents": tool_events,
            }
            if message["role"] == "user":
                raw_metadata = message_metadata.get(str(index))
                if not isinstance(raw_metadata, dict):
                    raw_metadata = {}
                message_context = raw_metadata.get("promptContext")
                if not isinstance(message_context, list):
                    message_context = []
                elif not all(isinstance(item, dict) for item in message_context):
                    message_context = []
                if (
                    not message_context
                    and not message_metadata
                    and len(user_message_indexes) <= 1
                ):
                    message_context = prompt_context
                if message_context:
                    message["context"] = message_context
                message_label = raw_metadata.get("label")
                if isinstance(message_label, str) and message_label.strip():
                    message["label"] = message_label
            if debug:
                message["raw"] = raw_message
            if message["text"] or images or reasoning_parts or tool_events or debug:
                messages.append(message)

        last_message = next((m for m in reversed(messages) if m.get("text")), None)
        history: dict[str, Any] = {
            "agentKey": agent_key,
            "title": metadata.get("title") or agent_key,
            "subtitle": metadata.get("subtitle") or metadata.get("smiles") or "",
            "metadata": metadata,
            "modelInfo": record.get("modelInfo", {}),
            "contextUsage": self._context_usage_for_history(record, raw_messages),
            "promptContext": prompt_context,
            "messages": messages,
            "lastMessage": (last_message or {}).get("text", ""),
        }
        if debug and raw_session is not None:
            history["rawSession"] = raw_session
        return history

    async def handle_list_agent_histories(self, data: dict[str, Any]) -> None:
        debug = bool(data.get("debug"))
        histories = [
            self._normalize_agent_history(agent_key, record, debug=debug)
            for agent_key, record in self._agent_session_records().items()
        ]
        await self.websocket.send_json(
            {
                "type": "agent-histories-response",
                "histories": histories,
            }
        )

    async def handle_get_agent_history(self, data: dict[str, Any]) -> None:
        agent_key = str(data.get("agentKey") or "")
        if not agent_key:
            raise ValueError("agentKey is required")
        debug = bool(data.get("debug"))
        records = self._agent_session_records()
        record = records.get(agent_key)
        fallback_task = self._chat_task_for_agent(agent_key, data, [])
        if record is None:
            record = {
                "agentKey": agent_key,
                "metadata": self._default_agent_metadata(agent_key, data),
                "memory": "",
                "modelInfo": {},
                "task": fallback_task.to_json(),
            }
        elif not isinstance(record.get("task"), dict):
            record = {
                **record,
                "metadata": {
                    **self._default_agent_metadata(agent_key, data),
                    **(
                        record.get("metadata")
                        if isinstance(record.get("metadata"), dict)
                        else {}
                    ),
                },
                "task": fallback_task.to_json(),
            }
        await self.websocket.send_json(
            {
                "type": "agent-history-response",
                "history": self._normalize_agent_history(
                    agent_key, record, debug=debug
                ),
            }
        )

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
        existing_agent = existing.get("agent") if existing else None
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
        metadata = self._default_agent_metadata(agent_key, data)
        task = self._chat_task_for_agent(agent_key, data, attachments)
        await self._send_processing_message(
            query,
            source="User",
            images=image_refs(attachments) if attachments else None,
        )
        callback_handler = CallbackHandler(self.websocket)
        agent = self.experiment.create_agent_with_experiment_state(
            task=task,
            agent_key=agent_key,
            agent_metadata=metadata,
            callback=callback_handler,
        )
        kind, target = self._agent_key_parts(agent_key)
        if kind == "reaction" and self.retro_synth_context is not None:
            node_id = str(data.get("nodeId") or metadata.get("nodeId") or target)
            self.retro_synth_context.node_id_to_charge_client[node_id] = agent

        async def run_and_report():
            if self.run_settings.prompt_debugging:
                await debug_prompt(agent, self.websocket)
            result = await agent.run(self.log_progress)
            await callback_handler.drain()
            self._record_latest_user_message_metadata(agent_key, task)
            self.experiment.add_to_context(agent, task, result)
            await self._send_processing_message(result, source="Agent")
            records = self._agent_session_records()
            await self.websocket.send_json(
                {
                    "type": "agent-history-response",
                    "history": self._normalize_agent_history(
                        agent_key,
                        records[agent_key],
                        debug=bool(data.get("debug")),
                    ),
                }
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
