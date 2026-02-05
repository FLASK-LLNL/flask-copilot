from fastapi import WebSocket
from callback_logger import CallbackLogger
from collections import defaultdict

from backend_helper_funcs import (
    Node,
    Reaction,
)
from retrosynthesis.context import RetrosynthesisContext
from charge_backend.moleculedb.purchasable import is_purchasable


async def set_reaction_alternative(
    node: Node,
    alternative_id: str,
    context: RetrosynthesisContext,
    websocket: WebSocket,
):
    """
    Sets a reaction alternative

    :param node: Parent node to reset the reaction thereof
    :param alternative_id: The unique ID of the alternative to choose
    :param context: Retrosynthesis context object
    :param websocket: Websocket to client
    """
    clogger = CallbackLogger(websocket, source="set_reaction_alternative")
    assert node.reaction is not None
    assert node.reaction.alternatives is not None
    try:
        alternative = next(
            a for a in node.reaction.alternatives if a.id == alternative_id
        )
    except StopIteration:
        await clogger.info(f"Alternative {alternative_id} not found for {node.id}")
        await websocket.send_json({"type": "complete"})
        return

    node.reaction.hoverInfo = alternative.hoverInfo

    # Clear subtree and levels for layouting
    await context.delete_subtree(node.id, websocket)

    # Update reaction information
    for alt in node.reaction.alternatives:
        if alt.status == "active":
            alt.status = "available"
    alternative.status = "active"
    if alternative.type == "ai":
        node.reaction.label = "FLASK AI"
        node.reaction.highlight = "red"
    elif alternative.type == "template":
        node.reaction.label = "Template"
        node.reaction.highlight = "yellow"
    elif alternative.type == "exact":
        node.reaction.label = "Exact"
        node.reaction.highlight = "normal"

    node.highlight = "normal"
    await context.update_node(node, websocket)

    # Loop over new nodes/edges (from the second level onwards)
    step_nodes: dict[int, list[Node]] = defaultdict(list)
    step_nodes[0].append(node)
    for i, step in enumerate(alternative.pathway):
        if i == 0:  # Skip root
            continue
        for smiles, label, parent in zip(step.smiles, step.label, step.parents):
            try:
                parent_node = step_nodes[i - 1][parent]

                mol_sources = is_purchasable(smiles)
                if mol_sources:
                    purchasable_str = f"Yes (via {', '.join(mol_sources)})"
                else:
                    purchasable_str = "No"
                child_node = Node(
                    id=context.new_node_id(),
                    smiles=smiles,
                    label=label,
                    hoverInfo=f"""# Molecule
**SMILES:** {smiles}

**Purchasable**? {purchasable_str}""",
                    level=node.level + i,
                    parentId=parent_node.id,
                    purchasable=(len(mol_sources) > 0),
                )
                step_nodes[i].append(child_node)
                await context.add_node(child_node, parent_node, websocket)

                # Create reaction for parent node if there are children
                if parent_node.reaction is None:
                    parent_node.reaction = Reaction(
                        "subreaction",
                        'Click "Other Reactions..." to compute alternative paths',
                        highlight=node.reaction.highlight,
                        label=node.reaction.label,
                        templatesSearched=False,
                    )
                    await context.update_node(parent_node, websocket)
            except IndexError as e:
                await clogger.error(
                    f"Index error when setting reaction alternative: parent index {parent} "
                    f"out of range for step {i-1}. Error: {e}"
                )

    await websocket.send_json({"type": "complete"})
