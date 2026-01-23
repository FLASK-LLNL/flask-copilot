import asyncio
from concurrent.futures import ProcessPoolExecutor
import os
from fastapi import WebSocket
from charge.clients.autogen import AutoGenAgent
from charge.servers import AiZynthTools as azf
from loguru import logger
from callback_logger import CallbackLogger
from typing import Any, Literal, Optional, Union, TYPE_CHECKING
from collections import deque

from backend_helper_funcs import (
    CallbackHandler,
    RetrosynthesisContext,
    calculate_positions,
    highlight_node,
    RetrosynthesisContext,
    Node,
    Edge,
    Reaction,
    ReactionAlternative,
    PathwayStep,
    recalculate_nodes_per_level,
)
from molecule_naming import smiles_to_html, MolNameFormat

from charge.tasks.RetrosynthesisTask import (
    TemplateFreeRetrosynthesisTask as RetrosynthesisTask,
    TemplateFreeReactionOutputSchema as ReactionOutputSchema,
)

from charge.experiments.AutoGenExperiment import AutoGenExperiment

if TYPE_CHECKING:
    import aizynthfinder.reactiontree

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


def generate_nodes_for_molecular_graph(
    reaction_path_dict: dict[int, azf.Node],
    retro_synth_context: RetrosynthesisContext,
    start_level: int = 0,
    molecule_name_format: Literal["brand", "iupac", "formula", "smiles"] = "brand",
    include_root_node: bool = True,
    root_node_id: str | None = None,
    id_offset: int = 0,
) -> tuple[list[Node], list[Edge]]:
    """Generate nodes and edges from reaction path dict"""
    nodes = []
    edges = []

    root_id = 0
    node_queue = [(reaction_path_dict[root_id], start_level)]  # (node, level)

    while node_queue:
        current_node, level = node_queue.pop(0)
        node_id = current_node.node_id
        smiles = current_node.smiles
        purchasable = current_node.purchasable
        leaf = current_node.is_leaf
        node_id_str = f"node_{node_id+id_offset}"
        hover_info = f"# Molecule \n **SMILES:** {smiles}\n"
        if leaf:
            if purchasable:
                hover_info += " - This molecule is purchasable.\n"
                # TODO: For ... chemprice
            else:
                hover_info += " - This molecule is NOT purchasable.\n"

        # A bit of an ugly way of dealing with this
        if current_node.parent_id is None:
            parent_node_id = None
        elif root_node_id is not None and current_node.parent_id == root_id:
            parent_node_id = root_node_id
        else:
            parent_node_id = f"node_{current_node.parent_id+id_offset}"

        if node_id != root_id or include_root_node:
            node = Node(
                id=node_id_str,
                smiles=smiles,
                label=smiles_to_html(smiles, molecule_name_format),
                hoverInfo=hover_info,
                level=level,
                parentId=parent_node_id,
                highlight=("red" if (leaf and not purchasable) else "normal"),
            )

            retro_synth_context.node_ids[node_id_str] = node
            nodes.append(node)

        if current_node.parent_id is not None:
            assert parent_node_id is not None
            edge = Edge(
                id=f"edge_{current_node.parent_id+id_offset}_{node_id+id_offset}",
                fromNode=parent_node_id,
                toNode=node_id_str,
                status="complete",
                label=None,
            )
            retro_synth_context.node_ids[parent_node_id].reaction = Reaction(
                "azf",
                "Reaction found with AiZynthFinder",
                highlight="yellow",
                label="Template",
                templatesSearched=False,
            )
            edges.append(edge)
            retro_synth_context.parents[node_id_str] = parent_node_id

        for child_id in current_node.children:
            child_node = reaction_path_dict[child_id]
            node_queue.append((child_node, level + 1))

    return nodes, edges


def make_reaction_alternative(
    rpath: azf.ReactionPath,
    id: int,
    tree: "aizynthfinder.reactiontree.ReactionTree",
    molecule_name_format: MolNameFormat,
) -> ReactionAlternative | None:
    """
    Creates a reaction alternative object from a given route
    """
    try:
        smarts = next(iter(tree.reactions())).metadata.get("template")
    except StopIteration:
        return None  # No reaction in tree
    try:
        prevalence = next(iter(tree.reactions())).metadata.get("library_occurence")
    except StopIteration:
        return None  # No reaction in tree
    pathway = []

    # BFS traversal over tree, merging precursors of each child for visualization purposes
    queue = deque([rpath.root])

    while queue:
        level_size = len(queue)
        level_smiles = []

        for _ in range(level_size):
            node = queue.popleft()
            level_smiles.append(node.smiles)

            # Add children to queue for next level
            if node.children:
                for child_id in node.children:
                    # Get the actual node from rpath.nodes using the child_id
                    if hasattr(rpath, "nodes") and child_id in rpath.nodes:
                        queue.append(rpath.nodes[child_id])

        # Add this level as a pathway step
        pathway.append(
            PathwayStep(
                level_smiles,
                [smiles_to_html(s, molecule_name_format) for s in level_smiles],
            )
        )

    return ReactionAlternative(
        f"reaction_{rpath.root.node_id}_{id}",
        f"Reaction {id+1} ({prevalence} occurences)",
        "template",
        "active" if id == 0 else "available",
        pathway,
        disabled=False,
        hoverInfo=f"Reaction SMARTS: `{smarts}`\n\nPattern appears {prevalence} times in database.",
    )


async def run_retro_planner(
    config_file: str,
    smiles: str,
    clogger: CallbackLogger,
    molecule_name_format: MolNameFormat,
    reaction_id: str = "azf",
) -> tuple[Reaction | None, list[azf.ReactionPath]]:
    """
    Runs AiZynthFinder (template-based multi-step retrosynthesis) on the given
    SMILES string and returns a Reaction object with all the alternatives found.

    :param config_file: Path to AiZynthFinder configuration yml file
    :param smiles: SMILES string to use
    :param clogger: Logger object that can return messages to the UI
    :param molecule_name_format: Desired default formatting for molecule names
    :param reaction_id: An optional string for a unique reaction ID
    :return: A 2-tuple of (Reaction object, list of routes) if routes found, or
             ``(None, [])`` if nothing was discovered.
    """
    await clogger.info(f"Running RetroPlanner for SMILES: {smiles}")
    planner = azf.RetroPlanner(configfile=config_file)

    _, _, routes = planner.plan(smiles)
    if len(routes) == 0:  # No routes found
        return None, []

    assert planner.finder is not None
    await clogger.info(f"Found {len(routes)} routes for {smiles}.")

    trees = planner.finder.routes.reaction_trees
    rpaths = [azf.ReactionPath(route) for route in routes]

    # All other routes become reaction alternatives
    alts: list[ReactionAlternative] = []
    for i, rpath in enumerate(rpaths):
        alt = make_reaction_alternative(rpath, i, trees[i], molecule_name_format)
        if alt is not None:
            alts.append(alt)

    try:
        root_smarts = next(iter(trees[0].reactions())).metadata.get("template")
    except StopIteration:
        return None, []
    try:
        root_prevalence = next(iter(trees[0].reactions())).metadata.get(
            "library_occurence"
        )
    except StopIteration:
        return None, []
    return (
        Reaction(
            reaction_id,
            f"""Reaction found with AiZynthFinder

Pattern occurrences in database: {root_prevalence}

Reaction SMARTS: `{root_smarts}`""",
            highlight="yellow",
            label="Template",
            alternatives=alts,
            templatesSearched=True,
        ),
        rpaths,
    )


async def template_based_retrosynthesis(
    start_smiles: str,
    config_file: str,
    context: RetrosynthesisContext,
    executor: ProcessPoolExecutor,
    websocket: WebSocket,
    available_tools: Optional[Union[str, list[str]]] = None,
    molecule_name_format: Literal["brand", "iupac", "formula", "smiles"] = "brand",
):
    """Stream positioned nodes and edges"""
    clogger = CallbackLogger(websocket, source="template_based_retrosynthesis")
    await clogger.info(
        f"Planning retrosynthesis for: {start_smiles} with available tools: {available_tools}."
    )

    # Generate root node
    root = Node(
        id="node_0",
        smiles=start_smiles,
        label=smiles_to_html(start_smiles, molecule_name_format),
        hoverInfo=f"# Root molecule \n **SMILES:** {start_smiles}",
        level=0,
        parentId=None,
        cost=None,
        bandgap=None,
        yield_=None,
        highlight="yellow",
        x=100,
        y=100,
    )
    await websocket.send_json({"type": "node", "node": root.json()})

    await clogger.info("Starting planning in executor...")
    reaction, routes = await run_retro_planner(
        config_file, start_smiles, clogger, molecule_name_format
    )
    if not reaction:
        await clogger.info(
            f"No synthesis routes found for {start_smiles}.",
            smiles=start_smiles,
        )
        # Send empty reaction
        await websocket.send_json(
            {
                "type": "node_update",
                "node": {
                    "id": root.id,
                    "highlight": "normal",
                    "reaction": Reaction(
                        "none",
                        'No exact or template-based reactions found.\n\nClick "Other Reactions..." to compute a path with FLASK AI.',
                        "empty",
                        label="No reaction",
                        templatesSearched=True,
                    ).json(),
                },
            }
        )
        await websocket.send_json({"type": "complete"})
        return

    # Use first route by default
    recalculate_nodes_per_level(context)
    nodes, edges = generate_nodes_for_molecular_graph(
        routes[-1].nodes, context, molecule_name_format=molecule_name_format
    )
    calculate_positions(nodes)

    # Notify frontend that computation completed
    context.node_ids[root.id].reaction = reaction
    await websocket.send_json(
        {
            "type": "node_update",
            "node": {
                "id": root.id,
                "highlight": "normal",
                "reaction": reaction.json(),
            },
        }
    )

    # Stream nodes from reaction path 0 with edges
    for node in nodes[1:]:
        # Find edge for this node
        edge = next((e for e in edges if e.toNode == node.id), None)

        # Send node
        await websocket.send_json({"type": "node", "node": node.json()})
        if edge:
            context.parents[node.id] = edge.fromNode
            await websocket.send_json({"type": "edge", "edge": edge.json()})

    await websocket.send_json({"type": "complete"})


async def ai_based_retrosynthesis(
    node_id: str,
    context: RetrosynthesisContext,
    query: Optional[str],
    constraint: Optional[str],
    websocket: WebSocket,
    experiment: AutoGenExperiment,
    config_file: str,
    available_tools: Optional[Union[str, list[str]]] = None,
    molecule_name_format: Literal["brand", "iupac", "formula", "smiles"] = "brand",
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
    num_nodes = len(context.node_ids)

    reasoning_summary = result.reasoning_summary
    await clogger.info(
        f"Retrosynthesis reasoning summary for {current_node.smiles}:\n{reasoning_summary}",
    )

    # Add reaction alternative
    ai_alternative = ReactionAlternative(
        "ai_reaction_0",
        "AI-Generated Reaction",
        "ai",
        "active",
        [
            PathwayStep([current_node.smiles], [current_node.label]),
            PathwayStep(
                list(result.reactants_smiles_list),
                [
                    smiles_to_html(s, molecule_name_format)
                    for s in result.reactants_smiles_list
                ],
            ),
        ],
        reasoning_summary,
        disabled=False,
    )
    if current_node.reaction is not None:
        if current_node.reaction.alternatives is None:
            current_node.reaction.alternatives = []
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
    current_node.reaction = ai_reaction

    # Update node with discovered reaction
    context.node_id_to_reasoning_summary[node_id] = reasoning_summary
    await websocket.send_json(
        {
            "type": "node_update",
            "node": {"id": node_id, "reaction": ai_reaction.json()},
        }
    )

    recalculate_nodes_per_level(context)

    for i, smiles in enumerate(result.reactants_smiles_list):
        node = Node(
            f"node_{num_nodes+i}",
            smiles,
            smiles_to_html(smiles, molecule_name_format),
            f"Discovered by {runner.model}",
            level,
            current_node.id,
        )
        edge = Edge(
            f"edge_{current_node.id}_{num_nodes + i}", node_id, node.id, "complete"
        )
        if node.id in context.node_ids:
            await clogger.warning(f"{node.id} already found in context: {node}")
        context.node_ids[node.id] = node
        context.parents[node.id] = node_id

        # Stream node directly
        calculate_positions([node], context.nodes_per_level[level])
        context.nodes_per_level[level] += 1
        node.highlight = "yellow"
        await websocket.send_json({"type": "node", "node": node.json()})
        await websocket.send_json({"type": "edge", "edge": edge.json()})

        # Find paths for the leaf nodes
        reaction, routes = await run_retro_planner(
            config_file, smiles, clogger, molecule_name_format
        )
        if reaction is None:
            await clogger.warning(f"No routes found for {smiles}. Skipping...")
            continue
        node.reaction = reaction

        # Use child nodes and edges of the first route
        child_nodes, child_edges = generate_nodes_for_molecular_graph(
            routes[0].nodes,
            context,
            start_level=level,
            include_root_node=False,
            root_node_id=node.id,
            id_offset=num_nodes,
        )
        node.highlight = "normal"
        await highlight_node(node, websocket, False)

        calculate_positions(child_nodes, context.nodes_per_level[level + 1])
        context.nodes_per_level[level + 1] += len(child_nodes)

        # Stream nodes back
        for cnode in child_nodes:
            await websocket.send_json({"type": "node", "node": cnode.json()})
        for cedge in child_edges:
            await websocket.send_json({"type": "edge", "edge": cedge.json()})

        num_nodes += routes[0].num_nodes

    await websocket.send_json({"type": "complete"})
