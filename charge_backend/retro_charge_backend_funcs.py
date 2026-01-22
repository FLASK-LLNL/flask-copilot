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
        node_id_str = f"node_{node_id}"
        hover_info = f"# Molecule \n **SMILES:** {smiles}\n"
        if leaf:
            if purchasable:
                hover_info += " - This molecule is purchasable.\n"
                # TODO: For ... chemprice
            else:
                hover_info += " - This molecule is NOT purchasable.\n"

        node = Node(
            id=node_id_str,
            smiles=smiles,
            label=smiles_to_html(smiles, molecule_name_format),
            hoverInfo=hover_info,
            level=level,
            parentId=(
                f"node_{current_node.parent_id}"
                if current_node.parent_id is not None
                else None
            ),
            highlight=("red" if (leaf and not purchasable) else "normal"),
        )

        retro_synth_context.node_ids[node_id_str] = node
        retro_synth_context.nodes_per_level[level] += 1
        nodes.append(node)

        if current_node.parent_id is not None:
            parent_node_id = f"node_{current_node.parent_id}"
            edge = Edge(
                id=f"edge_{current_node.parent_id}_{node_id}",
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
) -> ReactionAlternative:
    """
    Creates a reaction alternative object from a given route
    """
    smarts = next(iter(tree.reactions())).metadata["template"]
    prevalence = next(iter(tree.reactions())).metadata["library_occurence"]
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
        f"AiZynthFinder Reaction {id+1} ({prevalence} occurences in database)",
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
    alts = []
    for i, rpath in enumerate(rpaths):
        alts.append(make_reaction_alternative(rpath, i, trees[i], molecule_name_format))

    root_smarts = next(iter(trees[0].reactions())).metadata["template"]
    root_prevalence = next(iter(trees[0].reactions())).metadata["library_occurence"]
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
    nodes, edges = generate_nodes_for_molecular_graph(routes[0].nodes, context)
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


async def optimize_molecule_retro(
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
    """Optimize a molecule using retrosynthesis by node ID"""
    current_node = context.node_ids.get(node_id)
    assert current_node is not None, f"Node ID {node_id} not found"

    cur_node_id = node_id

    clogger = CallbackLogger(websocket, source="optimize_molecule_retro")

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

    # Update node with discovered reaction
    context.node_id_to_reasoning_summary[node_id] = reasoning_summary
    await websocket.send_json(
        {
            "type": "node_update",
            "node": {
                "id": node_id,
                "reaction": Reaction(
                    "ai_reaction_0",
                    reasoning_summary,
                    highlight="red",
                    label="FLASK AI",
                ).json(),
            },
        }
    )

    nodes: list[Node] = []
    edges: list[Edge] = []
    for i, smiles in enumerate(result.reactants_smiles_list):
        node = Node(
            f"node_{num_nodes+i}",
            smiles,
            smiles_to_html(smiles, molecule_name_format),
            "Discovered",
            level,
            current_node.id,
        )
        context.node_ids[node.id] = node
        context.parents[node.id] = node_id

        # Find paths for the leaf nodes
        routes, planner, trees = run_retro_planner(config_file, smiles)

        if not routes:
            await clogger.warning(f"No routes found for {smiles}. Skipping...")
            continue

        planner.last_route_used = 0

        reaction_path = azf.ReactionPath(
            route=routes[0],
            cur_num_nodes=num_nodes,
            root_parent_node_id=num_nodes + i,
        )
        child_nodes, child_edges = generate_nodes_for_molecular_graph(
            reaction_path.nodes, context, start_level=level
        )

        # Stream child nodes and edges
        nodes.extend(child_nodes)
        edges.extend(child_edges)

        num_nodes = reaction_path.num_nodes

    calculate_positions(nodes, context.nodes_per_level[level])

    for node, edge in zip(nodes, edges):
        await websocket.send_json({"type": "node", "node": node.json()})
        await websocket.send_json({"type": "edge", "edge": edge.json()})
    await websocket.send_json({"type": "complete"})
