"""
Retrosynthesis context management
"""

from collections import defaultdict
from dataclasses import dataclass, field
from charge.clients.agent_factory import Agent
from fastapi import WebSocket
from typing import Optional, Any
from uuid import uuid4

from backend_helper_funcs import (
    Node,
    Edge,
    PathwayStep,
    Reaction,
    ReactionAlternative,
    calculate_positions,
)


@dataclass
class RetrosynthesisContext:
    """
    Manages a retrosynthesis experiment
    """

    node_ids: dict[str, Node] = field(default_factory=dict)
    node_id_to_charge_client: dict[str, Agent] = field(default_factory=dict)
    node_id_to_reasoning_summary: dict[str, str] = field(default_factory=dict)
    nodes_per_level: dict[int, int] = field(default_factory=lambda: defaultdict(int))
    parents: dict[str, str] = field(default_factory=dict)
    edges: dict[str, Edge] = field(default_factory=dict)

    def reset(self):
        self.node_ids.clear()
        self.node_id_to_charge_client.clear()
        self.node_id_to_reasoning_summary.clear()
        self.nodes_per_level.clear()
        self.parents.clear()
        self.edges.clear()

    def recalculate_nodes_per_level(self) -> None:
        """
        Recalculate the nodes_per_level counts from scratch based on current nodes.

        This ensures accurate node counts at each level after deletions and is useful
        for proper positioning of newly added nodes.
        """
        # Clear and recalculate from scratch
        self.nodes_per_level.clear()
        for node in self.node_ids.values():
            self.nodes_per_level[node.level] += 1

    async def delete_subtree(
        self, node_id: str, websocket: Optional[WebSocket] = None
    ) -> list[str]:
        """
        Delete a subtree rooted at node_id from the retrosynthesis context.

        This function removes all descendant nodes (children, grandchildren, etc.)
        of the given node_id from all data structures in the Retrosynthesisself,
        and recalculates the nodes_per_level counts to ensure proper positioning
        of future nodes.

        :param node_id: The ID of the root node whose descendants should be deleted

        :return: List of deleted node IDs (including the descendants, but not the root node itself)
        """
        # Find all descendants using BFS
        descendants = []
        queue = [node_id]

        while queue:
            current = queue.pop(0)
            # Find all children of the current node
            children = [
                child_id
                for child_id, parent_id in self.parents.items()
                if parent_id == current
            ]
            descendants.extend(children)
            queue.extend(children)

        # Remove descendants from all data structures
        for desc_id in descendants:
            # Remove from node_ids
            if desc_id in self.node_ids:
                del self.node_ids[desc_id]

            # Remove from node_id_to_charge_client
            if desc_id in self.node_id_to_charge_client:
                del self.node_id_to_charge_client[desc_id]

            # Remove from node_id_to_reasoning_summary
            if desc_id in self.node_id_to_reasoning_summary:
                del self.node_id_to_reasoning_summary[desc_id]

            # Remove from parents
            if desc_id in self.parents:
                del self.parents[desc_id]

        # Recalculate nodes_per_level to ensure accurate counts for positioning
        self.recalculate_nodes_per_level()

        if websocket is not None:
            await websocket.send_json(
                {"type": "subtree_delete", "node": {"id": node_id}}
            )

        return descendants

    async def add_node(
        self,
        node: Node,
        websocket: Optional[WebSocket] = None,
        edge: Optional[Edge] = None,
    ):
        """
        Add a new node to the context and molecule graph.

        :param node: The new node to add
        :param websocket: If not None, notifies websocket of addition.
        :param edge: If not None, use this edge instead of generating one.
        """
        if node.id in self.node_ids:
            raise FileExistsError(f"Node with ID {node.id} already exists in context")
        self.node_ids[node.id] = node
        self.nodes_per_level[node.level] += 1

        # Recalculate positions
        all_nodes = list(self.node_ids.values())
        calculate_positions(all_nodes)

        if websocket is not None:
            await websocket.send_json({"type": "node", "node": node.json()})

        if node.parentId is not None:
            assert (
                node.parentId in self.node_ids
            )  # Parents must be added before children
            self.parents[node.id] = node.parentId
            if edge is None:
                edge = Edge(f"edge_{uuid4()}", node.parentId, node.id, "complete")
            self.edges[edge.id] = edge
            if websocket is not None:
                await websocket.send_json({"type": "edge", "edge": edge.json()})

    async def update_node(self, node: Node, websocket: Optional[WebSocket] = None):
        """
        Updates a node in the molecule graph.

        :param node: The new node to add
        :param websocket: If not None, notifies websocket of addition.
        """
        if node.id not in self.node_ids:
            raise FileNotFoundError(f"Node with ID {node.id} does not exist in context")
        self.node_ids[node.id] = node

        if websocket is not None:
            await websocket.send_json({"type": "node_update", "node": node.json()})

    async def delete_node(self, node_id: str, websocket: Optional[WebSocket] = None):
        """
        Deletes a node from the context and molecule graph.
        Raises a ``KeyError`` if the node ID is not found in the context.

        :param node_id: The node ID to delete.
        :param websocket: If not None, notifies websocket of deletion.
        """
        if node_id not in self.node_ids:
            raise KeyError(f"Node {node_id} not found in context")
        del self.node_ids[node_id]
        del self.node_id_to_charge_client[node_id]
        del self.node_id_to_reasoning_summary[node_id]
        del self.parents[node_id]

        # Erase neighboring edges
        for eid, edge in list(self.edges.items()):
            if edge.fromNode == node_id or edge.toNode == node_id:
                del self.edges[eid]

        if websocket is not None:
            await websocket.send_json({"type": "node_delete", "nodeId": node_id})

        self.recalculate_nodes_per_level()

    async def load_state(self, data: dict[str, Any]) -> None:
        """
        Restore a context from a dictionary serialization.
        Existing state is clobbered.

        The dictionary must contain "nodes" for any meaningful action.

        It may also contain "edges". If "edges" is NOT included, new
        edges will be generated. These will be "private" backend edges
        that aren't reported to the frontend. This is only a problem
        if they're ever "edge-update"-ed. Edges that don't match any
        node are (silently) ignored.

        :param dict: Dictionary with serialized nodes and, optionally, edges.
        """
        self.reset()  # Start from a clean slate

        # Sort the node list so parents are guaranteed to be added
        # first. If a parent is missing, the child's assertion (if
        # __debug__) will trigger.
        node_data = data.get("nodes", [])
        node_data.sort(
            key=lambda x: (
                x["parentId"] if "parentId" in x and x["parentId"] is not None else ""
            )
        )

        # We explicitly restore edges now because of the UUID. The
        # worldview is node-centric, so any edge that doesn't match a
        # node will just be ignored.
        edges = [Edge(**e) for e in data.get("edges", [])]

        # Deserialize each node. The "yield_"/"yield" fiasco needs to
        # be handled manually, but otherwise, this is pretty
        # automatic -- thanks, pydantic.
        for n in node_data:
            n["yield_"] = n["yield"]
            del n["yield"]
            node = Node(**n)

            # Find a matching edge, if one is needed.
            edge = (
                next(
                    filter(
                        lambda e: e.fromNode == node.parentId and e.toNode == node.id,
                        edges,
                    ),
                    None,
                )
                if node.parentId is not None
                else None
            )

            await self.add_node(node, edge=edge)

    def new_node_id(self) -> str:
        """
        Returns a new unique node ID.
        """
        i = 0
        while f"node_{i}" in self.node_ids:
            i += 1
        return f"node_{i}"
