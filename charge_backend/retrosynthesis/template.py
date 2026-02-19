from fastapi import WebSocket
from charge.servers import AiZynthTools as azf
from lc_conductor.callback_logger import CallbackLogger
from typing import Optional, Union, TYPE_CHECKING
from collections import deque

from backend_helper_funcs import (
    Node,
    Reaction,
    ReactionAlternative,
    PathwayStep,
    RunSettings,
)
from retrosynthesis.context import RetrosynthesisContext
from charge_backend.moleculedb.molecule_naming import (
    smiles_to_html,
    MolNameFormat,
)
from charge_backend.moleculedb.purchasable import is_purchasable
from charge_backend.retrosynthesis.database import find_exact_reactions
from retrosynthesis.mapping import build_mapped_reaction_dict_or_none

if TYPE_CHECKING:
    import aizynthfinder.reactiontree


async def generate_nodes_for_molecular_graph(
    reaction_path_dict: dict[int, azf.Node],
    retro_synth_context: RetrosynthesisContext,
    websocket: WebSocket,
    start_level: int = 0,
    molecule_name_format: MolNameFormat = "brand",
    include_root_node: bool = True,
    root_node_id: str | None = None,
) -> list[Node]:
    """Generate nodes and edges from reaction path dict"""
    nodes = []

    # (node, level, parent node)
    node_queue: list[tuple[azf.Node, int, Optional[Node]]] = []

    # Determine traversal start
    if not include_root_node:
        if root_node_id is None:
            raise ValueError(
                "If include_root_node is False, a root node ID must be provided"
            )
        parent = retro_synth_context.node_ids[root_node_id]

        # Start from AZF node 0's children
        for child_id in reaction_path_dict[0].children:
            child_node = reaction_path_dict[child_id]
            node_queue.append((child_node, start_level + 1, parent))
    else:
        node_queue.append((reaction_path_dict[0], start_level, None))

    while node_queue:
        current_node, level, parent = node_queue.pop(0)
        smiles = current_node.smiles
        node_id_str = retro_synth_context.new_node_id()
        hover_info = f"# Molecule\n\n**SMILES:** {smiles}\n\n"
        mol_sources = is_purchasable(smiles)
        if mol_sources:
            purchasable_str = f"Yes (via {', '.join(mol_sources)})"
        else:
            purchasable_str = "No"
        hover_info += f"**Purchasable**? {purchasable_str}\n"

        node = Node(
            id=node_id_str,
            smiles=smiles,
            label=smiles_to_html(smiles, molecule_name_format),
            hoverInfo=hover_info,
            level=level,
            parentId=parent.id if parent is not None else None,
            highlight="normal",
            purchasable=(len(mol_sources) > 0),
        )

        # If this node has children, attach a reaction + mapped reaction immediately.
        # This makes hover-highlighting work on multi-step default routes without
        # requiring the user to open/select alternatives.
        if getattr(current_node, "children", False):
            if len(current_node.children) > 0:
                node.reaction = Reaction(
                    "azf",
                    "Reaction found with AiZynthFinder",
                    highlight="yellow",
                    label="Template",
                    templatesSearched=False,
                )
                node.reaction.mappedReaction = build_mapped_reaction_dict_or_none(
                    reactants=[
                        reaction_path_dict[cid].smiles for cid in current_node.children
                    ],
                    products=[smiles],
                    log_msg="Failed to build rdkitjs mapped reaction for node_id={node_id} smiles={smiles}",
                    node_id=node_id_str,
                    smiles=smiles,
                )
        await retro_synth_context.add_node(node, parent, websocket)

        for child_id in current_node.children:
            child_node = reaction_path_dict[child_id]
            node_queue.append((child_node, level + 1, node))

    return nodes


def make_reaction_alternative(
    rpath: azf.ReactionPath,
    id: int,
    tree: "aizynthfinder.reactiontree.ReactionTree",
    molecule_name_format: MolNameFormat,
    all_inactive: bool = False,
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
    # Each queue item is a tuple: (node, parent_index_in_previous_level)
    queue = deque([(rpath.root, -1)])  # Root has no parent, use -1

    while queue:
        level_size = len(queue)
        level_smiles = []
        level_parents = []

        for i in range(level_size):
            node, parent_idx = queue.popleft()
            level_smiles.append(node.smiles)
            level_parents.append(parent_idx)

            # Add children to queue for next level
            if node.children:
                for child_id in node.children:
                    # Get the actual node from rpath.nodes using the child_id
                    if hasattr(rpath, "nodes") and child_id in rpath.nodes:
                        # Child's parent is at index i in the current level
                        queue.append((rpath.nodes[child_id], i))

        # Add this level as a pathway step
        pathway.append(
            PathwayStep(
                level_smiles,
                [smiles_to_html(s, molecule_name_format) for s in level_smiles],
                level_parents,
            )
        )

    return ReactionAlternative(
        f"reaction_{rpath.root.node_id}_{id}",
        f"Reaction {id+1} ({prevalence} occurences)",
        "template",
        "active" if id == 0 and not all_inactive else "available",
        pathway,
        disabled=False,
        hoverInfo=f"Reaction SMARTS: `{smarts}`\n\nPattern appears {prevalence} times in database.",
    )


async def run_retro_planner(
    config_file: str,
    smiles: str,
    clogger: CallbackLogger,
    run_settings: RunSettings,
    reaction_id: str = "azf",
    all_inactive: bool = False,
) -> tuple[Reaction | None, list[azf.ReactionPath]]:
    """
    Runs AiZynthFinder (template-based multi-step retrosynthesis) on the given
    SMILES string and returns a Reaction object with all the alternatives found.

    :param config_file: Path to AiZynthFinder configuration yml file
    :param smiles: SMILES string to use
    :param clogger: Logger object that can return messages to the UI
    :param run_settings: Desired run settings
    :param reaction_id: An optional string for a unique reaction ID
    :param all_inactive: If True, does not activate the first alternative
    :return: A 2-tuple of (Reaction object, list of routes) if routes found, or
             ``(None, [])`` if nothing was discovered.
    """
    await clogger.info(f"Running RetroPlanner for SMILES: {smiles}")
    report_init = False
    if azf.RetroPlanner.finder is None:
        report_init = True
        await clogger.info("Initializing AiZynthFinder")
    planner = azf.RetroPlanner(configfile=config_file)
    if report_init:
        await clogger.info("AiZynthFinder initialization complete")

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
        alt = make_reaction_alternative(
            rpath, i, trees[i], run_settings.molecule_name_format, all_inactive
        )
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
    websocket: WebSocket,
    run_settings: RunSettings,
    available_tools: Optional[Union[str, list[str]]] = None,
):
    """Stream positioned nodes and edges"""
    clogger = CallbackLogger(websocket, source="template_based_retrosynthesis")
    await clogger.info(
        f"Planning retrosynthesis for: `{start_smiles}` with available tools: {available_tools}."
    )

    # Generate root node
    mol_sources = is_purchasable(start_smiles)
    if mol_sources:
        purchasable_str = f"Yes (via {', '.join(mol_sources)})"
    else:
        purchasable_str = "No"
    root = Node(
        id="node_0",
        smiles=start_smiles,
        label=smiles_to_html(start_smiles, run_settings.molecule_name_format),
        hoverInfo=f"""# Root molecule
**SMILES:** {start_smiles}

**Purchasable**? {purchasable_str}""",
        level=0,
        parentId=None,
        cost=None,
        bandgap=None,
        yield_=None,
        purchasable=(len(mol_sources) > 0),
        highlight="yellow",
        x=100,
        y=100,
    )
    context.reset()  # Clear context
    await context.add_node(root, parent=None, websocket=websocket)

    await clogger.info("Searching for exact matches...")
    reaction = await find_exact_reactions(
        root, context, clogger, websocket, run_settings
    )
    if reaction is not None:  # Exact reactions were found in database
        root.reaction = reaction
        root.reaction.mappedReaction = build_mapped_reaction_dict_or_none(
            reactants=[
                n.smiles
                for nid, n in context.node_ids.items()
                if context.parents.get(nid) == root.id
            ],
            products=[root.smiles],
            log_msg="Failed to build rdkitjs mapped reaction for exact root_id={node_id} smiles={smiles}",
            node_id=root.id,
            smiles=root.smiles,
        )
        await context.update_node(root, websocket)
        await websocket.send_json({"type": "complete"})
        return

    await clogger.info("Running AiZynthFinder...")
    reaction, routes = await run_retro_planner(
        config_file, start_smiles, clogger, run_settings
    )
    if not reaction:
        await clogger.info(
            f"No synthesis routes found for {start_smiles}.",
            smiles=start_smiles,
        )
        # Send empty reaction
        root.reaction = Reaction(
            "none",
            'No exact or template-based reactions found.\n\nClick "Other Reactions..." to compute a path with FLASK AI.',
            "empty",
            label="No reaction",
            templatesSearched=True,
        )
        await context.update_node(root, websocket)
        await websocket.send_json({"type": "complete"})
        return

    # Use first route by default
    await generate_nodes_for_molecular_graph(
        routes[0].nodes,
        context,
        websocket=websocket,
        molecule_name_format=run_settings.molecule_name_format,
        include_root_node=False,
        root_node_id=root.id,
    )

    # Notify frontend that computation completed
    root.reaction = reaction
    root.reaction.mappedReaction = build_mapped_reaction_dict_or_none(
        reactants=[
            n.smiles
            for nid, n in context.node_ids.items()
            if context.parents.get(nid) == root.id
        ],
        products=[root.smiles],
        log_msg="Failed to build rdkitjs mapped reaction for root template reaction root_id={node_id} smiles={smiles}",
        node_id=root.id,
        smiles=root.smiles,
    )
    await context.update_node(root, websocket)
    await websocket.send_json({"type": "complete"})


async def compute_templates_for_node(
    node: Node,
    config_file: str,
    context: RetrosynthesisContext,
    websocket: WebSocket,
    run_settings: RunSettings,
    available_tools: Optional[Union[str, list[str]]] = None,
):
    """Computes all templates for node"""
    clogger = CallbackLogger(websocket, source="compute_templates_for_node")
    await clogger.info(
        f"Planning retrosynthesis for node {node.id} with available tools: {available_tools}."
    )
    await clogger.info("Running AiZynthFinder...")
    reaction, routes = await run_retro_planner(
        config_file,
        node.smiles,
        clogger,
        run_settings=run_settings,
        all_inactive=True,
    )

    if not reaction:
        await clogger.info(f"No template-based synthesis routes found for {node.id}.")
        # Send empty reaction
        if node.reaction is None:
            node.reaction = Reaction(
                "none",
                'No exact or template-based reactions found.\n\nClick "Other Reactions..." to compute a path with FLASK AI.',
                "empty",
                label="No reaction",
                templatesSearched=True,
            )
        else:
            node.reaction.templatesSearched = True
        await context.update_node(node, websocket)
        await websocket.send_json({"type": "complete"})
        return

    if node.reaction is None or node.reaction.id == "azf":
        node.reaction = reaction
    else:
        if not node.reaction.alternatives:
            node.reaction.alternatives = []
        if reaction.alternatives:
            node.reaction.alternatives.extend(reaction.alternatives)
        node.reaction.templatesSearched = True

    # Notify frontend that template reactions were computed
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
    await context.update_node(node, websocket)
    await websocket.send_json({"type": "complete"})
