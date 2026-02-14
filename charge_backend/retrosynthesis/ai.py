import os
from fastapi import WebSocket
from lc_conductor.callback_logger import CallbackLogger
from typing import Optional, Union

from backend_helper_funcs import (
    CallbackHandler,
    highlight_node,
    Node,
    Reaction,
    ReactionAlternative,
    PathwayStep,
    RunSettings,
)
from charge_backend.prompt_debugger import debug_prompt
from retrosynthesis.context import RetrosynthesisContext
from charge_backend.moleculedb.molecule_naming import (
    smiles_to_html,
)
from retrosynthesis.template import (
    generate_nodes_for_molecular_graph,
    run_retro_planner,
)
from charge_backend.moleculedb.purchasable import is_purchasable
from retrosynthesis.mapping import build_mapped_reaction_dict_or_none

from charge.tasks.RetrosynthesisTask import (
    TemplateFreeRetrosynthesisTask as RetrosynthesisTask,
    TemplateFreeReactionOutputSchema as ReactionOutputSchema,
)

from charge.experiments.AutoGenExperiment import AutoGenExperiment

RETROSYNTH_UNCONSTRAINED_USER_PROMPT_TEMPLATE = (
    "Provide a retrosynthetic pathway for the target molecule {target_molecule}. "
    + "The pathway should be provided as a tuple of reactants as SMILES and the product as SMILES. "
    + "Perform only single step retrosynthesis. Make sure the SMILES strings are valid. "
    + "Use tools to verify the SMILES strings and diagnose any issues that arise."
    + "Do the evaluation step-by-step. Propose a retrosynthetic step, then evaluate it. "
    + "If the evaluation fails, propose a new retrosynthetic step and evaluate it again. "
    + "Find the best possible retrosynthetic step, and use tools to see if the "
    + "proposed reactants are synthesizable. "
)

RETROSYNTH_CONSTRAINED_USER_PROMPT_TEMPLATE = (
    "Provide a retrosynthetic pathway for the target molecule {target_molecule}. "
    + "The pathway should be provided as a tuple of reactants as SMILES and the product as SMILES. "
    + "Perform only single step retrosynthesis. Make sure the SMILES strings are valid. "
    + "Use tools to verify the SMILES strings and diagnose any issues that arise. "
    + "The following reactant cannot be used in the retrosynthetic step: {constrained_reactant}. "
    + "Do the evaluation step-by-step. Propose a retrosynthetic step, then evaluate it. "
    + "If the evaluation fails, propose a new retrosynthetic step and evaluate it again. "
)


async def ai_based_retrosynthesis(
    node_id: str,
    context: RetrosynthesisContext,
    query: Optional[str],
    constraint: Optional[str],
    websocket: WebSocket,
    experiment: AutoGenExperiment,
    config_file: str,
    run_settings: RunSettings,
    available_tools: Optional[Union[str, list[str]]] = None,
):
    """Performs template-free retrosynthesis using the AI orchestrator."""
    clogger = CallbackLogger(websocket, source="ai_based_retrosynthesis")
    current_node = context.node_ids.get(node_id)
    if current_node is None:
        await clogger.error(f"Node ID {node_id} not found")
        await websocket.send_json({"type": "complete"})
        return

    await clogger.info(
        f"Finding synthesis pathway to {current_node.smiles}... with available tools: {available_tools}.",
        smiles=current_node.smiles,
    )

    if node_id in context.node_id_to_charge_client:
        # Existing context
        runner = context.node_id_to_charge_client[node_id]
    else:
        # New context
        agent_name = experiment.agent_pool.create_agent_name(
            prefix=f"retrosynth_{node_id}_"
        )
        runner = experiment.create_agent_with_experiment_state(
            task=None,
            agent_name=agent_name,
            callback=CallbackHandler(websocket),
        )
        context.node_id_to_charge_client[node_id] = runner

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
            f"\n\nAdditionally, adhere to the following requirements: {query}"
        )

    retro_task = RetrosynthesisTask(
        user_prompt=user_prompt, server_urls=available_tools
    )
    runner.task = retro_task

    if os.getenv("CHARGE_DISABLE_OUTPUT_VALIDATION", "0") == "1":
        retro_task.structured_output_schema = None
        await clogger.warning(
            "Structure validation disabled for RetrosynthesisTask output schema."
        )

    await clogger.info(
        f"Finding synthesis routes for {current_node.smiles} using available tools: {available_tools}."
    )

    # Run task
    await highlight_node(current_node, websocket, True)
    if run_settings.prompt_debugging:
        await debug_prompt(runner, websocket)
    output = await runner.run()

    if os.getenv("CHARGE_DISABLE_OUTPUT_VALIDATION", "0") == "1":
        await clogger.warning(
            "Structure validation disabled for RetrosynthesisTask output schema."
            "Returning text results without validation first before post-processing."
        )
    await clogger.info(f"Results: {output}")

    result = ReactionOutputSchema.model_validate_json(output)

    await highlight_node(current_node, websocket, False)

    level = current_node.level + 1

    reasoning_summary = result.reasoning_summary
    await clogger.info(
        f"Retrosynthesis reasoning summary for {current_node.smiles}:\n{reasoning_summary}",
    )

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
    if current_node.reaction is not None:
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
    context.node_id_to_reasoning_summary[node_id] = reasoning_summary
    await context.update_node(current_node, websocket)

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
        await context.add_node(node, current_node, websocket)

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
        if node.reaction is not None:
            child_smiles = [
                n.smiles
                for nid, n in context.node_ids.items()
                if context.parents.get(nid) == node.id
            ]
            node.reaction.mappedReaction = build_mapped_reaction_dict_or_none(
                reactants=child_smiles,
                products=[node.smiles],
                log_msg="Failed to build rdkitjs mapped reaction for template node_id={node_id} smiles={smiles}",
                node_id=node.id,
                smiles=node.smiles,
            )
        await context.update_node(node, websocket)  # Also disables highlight

    await websocket.send_json({"type": "complete"})
