from loguru import logger

from fastapi import WebSocket
import asyncio
import copy
from typing import Dict

from backend_helper_funcs import Node, Edge, RetroSynthesisContext
from charge.servers.AiZynthTools import RetroPlanner, aizynth_funcs


def generate_tree_structure(
    reaction_path_dict: Dict[int, aizynth_funcs.Node],
    retro_synth_context: RetroSynthesisContext,
):
    """Generate nodes and edges from reaction path dict"""
    nodes = []
    edges = []

    root_id = 0
    node_queue = [(reaction_path_dict[root_id], 0)]  # (node, level)

    while node_queue:
        current_node, level = node_queue.pop(0)
        node_id = current_node.node_id
        smiles = current_node.smiles
        purchasable = current_node.purchasable
        intermediate = not (current_node.is_root or current_node.is_leaf)
        leaf = current_node.is_leaf
        node_id_str = f"node_{node_id}"
        hover_info = f"# Molecule \n **SMILES:** {smiles}\n - Level: {level}\n"
        if leaf:
            if purchasable:
                hover_info += " - This molecule is purchasable.\n"
            else:
                hover_info += " - This molecule is NOT purchasable.\n"

        if intermediate:
            hover_info += " - Reaction intermediate\n"
        node = Node(
            id=node_id_str,
            smiles=smiles,
            label=smiles,
            hoverInfo=hover_info,
            level=level,
            parentId=(
                f"node_{current_node.parent_id}"
                if current_node.parent_id is not None
                else None
            ),
            cost=None,
            bandgap=None,
            yield_=None,
            highlight=leaf and not purchasable,
        )

        retro_synth_context.node_ids[node_id_str] = node

        # Map by smiles in case ever needed
        retro_synth_context.node_by_smiles[smiles] = node

        nodes.append(node)
        if current_node.parent_id is not None:

            edge = Edge(
                id=f"edge_{current_node.parent_id}_{node_id}",
                fromNode=f"node_{current_node.parent_id}",
                toNode=node_id_str,
                status="complete",
                label=None,
            )

            edges.append(edge)
        for child_id in current_node.children:
            child_node = reaction_path_dict[child_id]
            node_queue.append((child_node, level + 1))

    return nodes, edges


def calculate_positions(nodes: list[Node]):
    """
    Calculate positions for all nodes (matching frontend logic).
    """
    BOX_WIDTH = 270  # Must match with javascript!
    BOX_GAP = 160  # Must match with javascript!
    level_gap = BOX_WIDTH + BOX_GAP
    node_spacing = 150

    # Group by level
    levels = {}
    for node in nodes:
        level = node.level
        if level not in levels:
            levels[level] = []
        levels[level].append(node)

    # Position nodes
    positioned: list[Node] = []
    for node in nodes:
        level_nodes = levels[node.level]
        index_in_level = level_nodes.index(node)

        positioned_node = copy.deepcopy(node)
        positioned_node.x = 100 + node.level * level_gap
        positioned_node.y = 100 + index_in_level * node_spacing
        positioned.append(positioned_node)

    return positioned


async def aizynth_retro(
    start_smiles: str,
    planner: RetroPlanner,
    retro_synth_context: RetroSynthesisContext,
    websocket: WebSocket,
):
    """Stream positioned nodes and edges"""
    logger.info(f"Planning retrosynthesis for: {start_smiles}")

    # Generate and position entire tree upfront

    root = Node(
        id="node_0",
        smiles=start_smiles,
        label=start_smiles,
        hoverInfo=f"# Root molecule \n **SMILES:** {start_smiles}",
        level=0,
        parentId=None,
        cost=None,
        bandgap=None,
        yield_=None,
        highlight=True,
        x=100,
        y=100,
    )
    await websocket.send_json({"type": "node", **root.json()})
    tree, stats, routes = planner.plan(start_smiles)
    if not routes:
        await websocket.send_json(
            {
                "type": "response",
                "message": f"No synthesis routes found for {start_smiles}.",
                "smiles": start_smiles,
            }
        )
        await websocket.send_json({"type": "complete"})
        return
    logger.info(f"Found {len(routes)} routes for {start_smiles}.")

    reaction_path = aizynth_funcs.ReactionPath(route=routes[0])
    nodes, edges = generate_tree_structure(reaction_path.nodes, retro_synth_context)
    logger.info(f"Generated {len(nodes)} nodes and {len(edges)} edges.")

    positioned_nodes = calculate_positions(nodes)

    # Create node map
    # Stream root first

    await asyncio.sleep(0.8)

    # Stream remaining nodes with edges
    for i in range(1, len(positioned_nodes)):
        node = positioned_nodes[i]

        # Find edge for this node
        edge = next((e for e in edges if e.toNode == node.id), None)

        if edge:
            # Send edge with computing status
            edge_data = {
                "type": "edge",
                **edge.json(),
            }
            edge_data["label"] = f"Computing: {edge.label}"
            edge_data["toNode"] = node.id
            await websocket.send_json(edge_data)

            await asyncio.sleep(0.6)

            # Send node
            await websocket.send_json({"type": "node", **node.json()})

            # Update edge to complete
            edge_complete = {
                "type": "edge_update",
                "id": edge.id,
                "status": "complete",
                "label": edge.label,
            }
            await websocket.send_json(edge_complete)

            await asyncio.sleep(0.2)
    await websocket.send_json(
        {"type": "node_update", "id": root.id, "highlight": False}
    )

    await websocket.send_json({"type": "complete"})
