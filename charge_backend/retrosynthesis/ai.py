import os
import asyncio
from fastapi import WebSocket
from lc_conductor.callback_logger import CallbackLogger
from typing import Awaitable, Callable, Optional

from charge_backend.backend_helper_funcs import (
    CallbackHandler,
    highlight_node,
    Node,
    Reaction,
    ReactionAlternative,
    PathwayStep,
    FlaskRunSettings,
)
from charge_backend.prompt_debugger import debug_prompt
from charge_backend.flask_experiment import FlaskExperiment, GraphContext
from charge_backend.moleculedb.molecule_naming import (
    smiles_to_html,
)
from charge_backend.retrosynthesis.template import (
    generate_nodes_for_molecular_graph,
    run_retro_planner,
)
from charge_backend.moleculedb.purchasable import is_purchasable
from charge_backend.retrosynthesis.mapping import build_mapped_reaction_dict_or_none
from charge_backend.retrosynthesis.database import find_exact_reactions

from charge_backend.retrosynthesis.retrosynthesis_task import (
    TemplateFreeRetrosynthesisTask as RetrosynthesisTask,
    TemplateFreeReactionOutputSchema as ReactionOutputSchema,
)

from lc_conductor import ToolRuntime


RETROSYNTH_UNCONSTRAINED_USER_PROMPT_TEMPLATE = (
    "Provide a retrosynthetic pathway for the target molecule `{target_molecule}`. "
    + "If there are `*`, the `*` indicate the boundaries of the polymer repeat unit."
    + "The pathway should be provided as a tuple of reactants as SMILES and the product as SMILES. "
    + "Perform only single step retrosynthesis. Make sure the SMILES strings are valid. "
    + "Use tools to verify the SMILES strings and diagnose any issues that arise."
    + "Do the evaluation step-by-step. Propose a retrosynthetic step, then evaluate it. "
    + "If the evaluation fails, propose a new retrosynthetic step and evaluate it again. "
    + "Find the best possible retrosynthetic step, and use tools to see if the "
    + "proposed reactants are synthesizable. "
)

RETROSYNTH_CONSTRAINED_USER_PROMPT_TEMPLATE = (
    "Provide a retrosynthetic pathway for the target molecule `{target_molecule}`. "
    + "If there are `*`, the `*` indicate the boundaries of the polymer repeat unit."
    + "The pathway should be provided as a tuple of reactants as SMILES and the product as SMILES. "
    + "Perform only single step retrosynthesis. Make sure the SMILES strings are valid. "
    + "Use tools to verify the SMILES strings and diagnose any issues that arise. "
    + "The following reactant cannot be used in the retrosynthetic step: {constrained_reactant}. "
    + "Do the evaluation step-by-step. Propose a retrosynthetic step, then evaluate it. "
    + "If the evaluation fails, propose a new retrosynthetic step and evaluate it again. "
)


async def ai_based_retrosynthesis(
    node_id: str,
    query: Optional[str],
    constraint: Optional[str],
    websocket: WebSocket,
    experiment: FlaskExperiment,
    config_file: str,
    run_settings: FlaskRunSettings,
    tool_runtime: ToolRuntime,
    attachments: Optional[list[dict[str, object]]] = None,
    history_callback: Optional[Callable[[], Awaitable[None]]] = None,
):
    """Performs template-free retrosynthesis using the AI orchestrator."""
    clogger = CallbackLogger(websocket, source="ai_based_retrosynthesis")
    context = experiment.graph_context
    current_node = context.node_ids.get(node_id)
    if current_node is None:
        await clogger.error(f"Node ID {node_id} not found")
        await websocket.send_json({"type": "complete"})
        return

    await clogger.info(
        f"Finding synthesis pathway to {current_node.smiles}... with available tools: {tool_runtime.tool_summary()}.",
        smiles=current_node.smiles,
    )

    agent_key = f"reaction:{node_id}"
    callback_handler: CallbackHandler | None = None
    if agent_key in experiment.agent_registry:
        # Existing context
        runner = experiment.agent_registry[agent_key]
    else:
        # New context
        callback_handler = CallbackHandler(
            websocket, agent_key=agent_key, on_agent_update=history_callback
        )
        runner = experiment.create_agent_with_experiment_state(
            task=None,
            agent_key=agent_key,
            callback=callback_handler,
        )

    if callback_handler is None and isinstance(runner.callback, CallbackHandler):
        callback_handler = runner.callback
    if callback_handler is not None:
        callback_handler.agent_key = agent_key
        callback_handler.on_agent_update = history_callback

    if constraint:
        user_prompt = RETROSYNTH_CONSTRAINED_USER_PROMPT_TEMPLATE.format(
            target_molecule=current_node.smiles,
            constrained_reactant=constraint,
        )
    else:
        user_prompt = RETROSYNTH_UNCONSTRAINED_USER_PROMPT_TEMPLATE.format(
            target_molecule=current_node.smiles
        )

    user_prompt += "\nDouble check the reactants with the `predict_reaction_products` tool to see if the products are equivalent to the given product. If there is any inconsistency (canonicalize both sides of the equation first), log it and try some other set of reactants."
    if query is not None:
        user_prompt += (
            f"\n\nAdditionally, adhere to the following requirements:\n{query}\n\n"
        )

    retro_task = RetrosynthesisTask(
        user_prompt=user_prompt,
        attachments=attachments or [],
        **tool_runtime.task_kwargs(),
    )
    runner.task = retro_task

    if os.getenv("CHARGE_DISABLE_OUTPUT_VALIDATION", "0") == "1":
        retro_task.structured_output_schema = None
        await clogger.warning(
            "Structure validation disabled for RetrosynthesisTask output schema."
        )

    await clogger.info(
        f"Finding synthesis routes for {current_node.smiles} using available tools: {tool_runtime.tool_summary()}."
    )

    # Run task
    await highlight_node(current_node, websocket, True)
    if run_settings.prompt_debugging:
        await debug_prompt(runner, websocket)
    output = await runner.run()
    if callback_handler is not None:
        await callback_handler.drain()
    experiment.add_to_context(runner, retro_task, output)

    if os.getenv("CHARGE_DISABLE_OUTPUT_VALIDATION", "0") == "1":
        await clogger.warning(
            "Structure validation disabled for RetrosynthesisTask output schema."
            "Returning text results without validation first before post-processing."
        )

    result = ReactionOutputSchema.model_validate_json(output)

    await highlight_node(current_node, websocket, False)

    level = current_node.level + 1

    reasoning_summary = result.reasoning_summary
    await clogger.info(
        f"Retrosynthesis reasoning summary for {current_node.smiles}:\n{reasoning_summary}",
    )
    await asyncio.sleep(0)

    # Add reaction alternative
    ai_alternative = ReactionAlternative(
        f"ai_reaction_{node_id}",
        "AI-Generated Reaction",
        "ai",
        "active",
        [
            PathwayStep([current_node.smiles], [current_node.label], [-1]),
            PathwayStep(
                list(result.reactants_smiles_list),
                [
                    smiles_to_html(s, run_settings.molecule_name_format)
                    for s in result.reactants_smiles_list
                ],
                [0 for _ in range(len(result.reactants_smiles_list))],
            ),
        ],
        reasoning_summary,
        disabled=False,
    )
    if current_node.reaction is not None and not isinstance(
        current_node.reaction, dict
    ):
        if current_node.reaction.alternatives is None:
            current_node.reaction.alternatives = []
        for alt in current_node.reaction.alternatives:
            if alt.status == "active":
                alt.status = "available"
        current_node.reaction.alternatives.insert(0, ai_alternative)
        alternatives = current_node.reaction.alternatives
        templates_searched = current_node.reaction.templatesSearched
    else:
        alternatives = [ai_alternative]
        templates_searched = False

    ai_reaction = Reaction(
        "ai_reaction_0",
        reasoning_summary,
        highlight="red",
        label="FLASK AI",
        alternatives=alternatives,
        templatesSearched=templates_searched,
    )
    ai_reaction.mappedReaction = build_mapped_reaction_dict_or_none(
        reactants=list(result.reactants_smiles_list),
        products=[current_node.smiles],
        log_msg="Failed to build rdkitjs mapped reaction for AI reaction node_id={node_id} product_smiles={smiles}",
        node_id=current_node.id,
        smiles=current_node.smiles,
    )
    current_node.reaction = ai_reaction

    # Update node with discovered reaction
    await context.update_node(current_node, websocket)
    await asyncio.sleep(0)

    context.recalculate_nodes_per_level()

    new_nodes: list[Node] = []
    purchasable: list[bool] = []
    for smiles in result.reactants_smiles_list:
        node_id_str = context.new_node_id()
        mol_sources = is_purchasable(smiles)
        if mol_sources:
            purchasable_str = f"Yes (via {', '.join(mol_sources)})"
        else:
            purchasable_str = "No"
        node = Node(
            node_id_str,
            smiles,
            smiles_to_html(smiles, run_settings.molecule_name_format),
            f"Discovered by {runner.model}.\n\n**Purchasable**? {purchasable_str}",
            level,
            current_node.id,
            purchasable=(len(mol_sources) > 0),
        )
        new_nodes.append(node)
        purchasable.append(len(mol_sources) > 0)
        # Add and stream node directly
        await context.add_node(node, websocket)
        await asyncio.sleep(0)

    for node, purch in zip(new_nodes, purchasable):
        if purch:  # Skip purchasable nodes unless explicitly asked for
            continue

        # Highlight node because we are looking for templates
        await highlight_node(node, websocket, True)

        # Find paths for the leaf nodes
        reaction, routes = await run_retro_planner(
            config_file, node.smiles, clogger, run_settings
        )
        if reaction is None:
            await clogger.warning(f"No routes found for {node.smiles}. Skipping...")
            continue
        node.reaction = reaction

        # Use child nodes and edges of the first route
        await generate_nodes_for_molecular_graph(
            routes[0].nodes,
            context,
            websocket,
            start_level=level,
            include_root_node=False,
            root_node_id=node.id,
        )

        # Attach mapped reaction for immediate reactants -> product.
        # This is needed so hover-highlighting works for the first template step
        # discovered after an AI-generated step.
        attach_mapped_reaction(node, context)
        await context.update_node(node, websocket)  # Also disables highlight

    await websocket.send_json({"type": "complete"})


def attach_mapped_reaction(node: Node, context: GraphContext) -> None:
    """Attach mapped reaction to the reaction object."""
    reaction = node.reaction
    if reaction is not None:
        child_smiles = [
            n.smiles
            for nid, n in context.node_ids.items()
            if context.parents.get(nid) == node.id
        ]
        reaction.mappedReaction = build_mapped_reaction_dict_or_none(
            reactants=child_smiles,
            products=[node.smiles],
            log_msg="Failed to build rdkitjs mapped reaction for template node_id={node_id} smiles={smiles}",
            node_id=node.id,
            smiles=node.smiles,
        )


async def db_then_ai_retrosynthesis(
    node_id: str,
    query: Optional[str],
    constraint: Optional[str],
    websocket: WebSocket,
    experiment: FlaskExperiment,
    config_file: str,
    run_settings: FlaskRunSettings,
    tool_runtime: ToolRuntime,
    attachments: Optional[list[dict[str, object]]] = None,
    history_callback: Optional[Callable[[], Awaitable[None]]] = None,
):
    """Searches the reaction database first; falls back to AI-based retrosynthesis if no exact match is found."""
    clogger = CallbackLogger(websocket, source="db_then_ai_retrosynthesis")
    context = experiment.graph_context
    node = context.node_ids.get(node_id)
    if node is None:
        await clogger.error(f"Node ID {node_id} not found")
        await websocket.send_json({"type": "complete"})
        return

    await clogger.info("Searching for exact matches...")
    reaction = await find_exact_reactions(
        node, context, clogger, websocket, run_settings
    )
    if reaction is not None:
        node.reaction = reaction
        attach_mapped_reaction(node, context)
        await context.update_node(node, websocket)
        await websocket.send_json({"type": "complete"})
        return

    await clogger.info("No exact matches found. Computing with AI orchestrator...")
    await ai_based_retrosynthesis(
        node_id,
        query,
        constraint,
        websocket,
        experiment,
        config_file,
        run_settings,
        tool_runtime,
        attachments,
        history_callback,
    )
