from loguru import logger
from callback_logger import CallbackLogger

from fastapi import WebSocket
import asyncio
import copy
from typing import Dict

from backend_helper_funcs import Node, Edge, RetrosynthesisContext, calculate_positions
import charge.servers.AiZynthTools as AiZynthFuncs
from molecule_naming import smiles_to_html


def generate_tree_structure(
    reaction_path_dict: Dict[int, Node],
    retro_synth_context: RetrosynthesisContext,
    start_level: int = 0,
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
            label=smiles_to_html(smiles),
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
        retro_synth_context.azf_nodes[node_id_str] = current_node
        retro_synth_context.nodes_per_level[level] += 1
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
            retro_synth_context.parents[node_id_str] = f"node_{current_node.parent_id}"

        for child_id in current_node.children:
            child_node = reaction_path_dict[child_id]
            node_queue.append((child_node, level + 1))

    return nodes, edges


async def aizynth_retro(
    start_smiles: str,
    planner: AiZynthFuncs.RetroPlanner,
    retro_synth_context: RetrosynthesisContext,
    websocket: WebSocket,
):
    clogger = CallbackLogger(websocket)
    """Stream positioned nodes and edges"""
    await clogger.info(f"Planning retrosynthesis for: {start_smiles}")

    # Generate and position entire tree upfront

    root = Node(
        id="node_0",
        smiles=start_smiles,
        label=smiles_to_html(start_smiles),
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
    await websocket.send_json({"type": "node", "node": root.json()})
    tree, stats, routes = planner.plan(start_smiles)
    if not routes:
        await websocket.send_json(
            {
                "type": "response",
                "message": {
                    "source": "AiZynthFinder",
                    "message": f"No synthesis routes found for {start_smiles}.",
                    "smiles": start_smiles,
                },
            }
        )
        await websocket.send_json({"type": "complete"})
        return
    await clogger.info(f"Found {len(routes)} routes for {start_smiles}.")

    reaction_path = AiZynthFuncs.ReactionPath(route=routes[0])
    nodes, edges = generate_tree_structure(reaction_path.nodes, retro_synth_context)
    await clogger.info(f"Generated {len(nodes)} nodes and {len(edges)} edges.")

    calculate_positions(nodes)

    # Create node map
    # Stream root first

    await asyncio.sleep(0.8)

    # Stream remaining nodes with edges
    for i in range(1, len(nodes)):
        node = nodes[i]

        # Find edge for this node
        edge = next((e for e in edges if e.toNode == node.id), None)

        if edge:
            # Send edge with computing status
            edge_data = {
                "type": "edge",
                "edge": edge.json(),
            }
            edge_data["edge"]["label"] = f"Computing: {edge.label}"
            edge_data["edge"]["toNode"] = node.id
            await websocket.send_json(edge_data)

            await asyncio.sleep(0.6)

            # Send node
            await websocket.send_json({"type": "node", "node": node.json()})

            # Update edge to complete
            edge_complete = {
                "type": "edge_update",
                "edge": {
                    "id": edge.id,
                    "status": "complete",
                    "label": edge.label,
                },
            }
            await websocket.send_json(edge_complete)

            await asyncio.sleep(0.2)
    await websocket.send_json(
        {"type": "node_update", "node": {"id": root.id, "highlight": "normal"}}
    )

    await websocket.send_json({"type": "complete"})
