###############################################################################
## Copyright 2025-2026 Lawrence Livermore National Security, LLC.
## See the top-level LICENSE file for details.
##
## SPDX-License-Identifier: Apache-2.0
###############################################################################
"""
FLASK-specific experiment with graph-based context
"""

import asyncio

from charge.clients.agent_factory import Agent
from charge.experiments.experiment import Experiment
from collections import defaultdict
from fastapi import WebSocket
from pydantic import Field, TypeAdapter
from pydantic.dataclasses import dataclass
from typing import Optional, Any, TypeAlias
from uuid import uuid4

from charge_backend.backend_helper_funcs import (
    Node,
    Edge,
    PathwayStep,
    Reaction,
    ReactionAlternative,
    calculate_positions,
)


# Formerly RetroSynthesisContext
@dataclass
class GraphContext:
    """
    Manages nodes and edges in graph
    """

    node_ids: dict[str, Node] = Field(default_factory=dict)
    nodes_per_level: dict[int, int] = Field(default_factory=lambda: defaultdict(int))
    parents: dict[str, str] = Field(default_factory=dict)
    edges: dict[str, Edge] = Field(default_factory=dict)

    def reset(self) -> None:
        self.node_ids.clear()
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
        Delete a subtree rooted at node_id from the graph context.

        This function removes all descendant nodes (children, grandchildren, etc.)
        of the given node_id from all data structures in the graph,
        and recalculates the nodes_per_level counts to ensure proper positioning
        of future nodes. The node corresponding to "node_id" is not deleted.

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

            # Remove edge and parent
            if desc_id in self.parents:
                parent_id = self.parents[desc_id]

                # NOTE (trb): there should only be one, and technically only
                # the check for `e.to_node == desc_id` should be
                # necessary, but this could be a list comprehension if
                # necessary.
                eid = next(
                    (
                        k
                        for k, e in self.edges.items()
                        if e.fromNode == parent_id and e.toNode == desc_id
                    ),
                    None,
                )

                del self.edges[eid]
                del self.parents[desc_id]

        # Recalculate nodes_per_level to ensure accurate counts for positioning
        self.recalculate_nodes_per_level()

        if websocket is not None:
            await websocket.send_json(
                {"type": "subtree_delete", "node": {"id": node_id}}
            )

        return descendants

    # NOTE (trb): This is the synchronous part of "add_node". FWIW, I
    # find the asyncio model clunky at best, and it's very annoying
    # that it can so easily pollute the interface with everything up
    # to main() quickly getting `async`-ified. In this case, it's the
    # websocket's fault. If websocket is None, this function, and this
    # entire interface, is completely synchronous. Oh well.
    def _add_node(self, node: Node, edge: Edge | None) -> Edge | None:
        """
        Handles the synchronous component of add_node. Returns the
        edge that was added to the graph, if any.
        """

        if node.id in self.node_ids:
            raise FileExistsError(f"Node with ID {node.id} already exists in context")
        self.node_ids[node.id] = node
        self.nodes_per_level[node.level] += 1

        # Recalculate positions
        all_nodes = list(self.node_ids.values())
        calculate_positions(all_nodes)

        if node.parentId is not None:
            assert (
                node.parentId in self.node_ids
            )  # Parents must be added before children
            self.parents[node.id] = node.parentId

            if edge is None:
                edge = Edge(f"edge_{uuid4()}", node.parentId, node.id, "complete")

            # Verify that a user-input edge is sane.
            assert edge.id not in self.edges
            assert edge.toNode == node.id and edge.fromNode == node.parentId

            self.edges[edge.id] = edge

        return edge

    async def add_node(
        self,
        node: Node,
        websocket: Optional[WebSocket] = None,
        edge: Optional[Edge] = None,
    ) -> None:
        """
        Add a new node to the context and molecule graph.

        :param node: The new node to add
        :param websocket: If not None, notifies websocket of addition.
        :param edge: If not None, use this edge instead of generating one.
        """
        added_edge = self._add_node(node, edge)

        if websocket is not None:
            await websocket.send_json({"type": "node", "node": node.json()})
            if edge is None and added_edge is not None:
                await websocket.send_json({"type": "edge", "edge": added_edge.json()})

    async def update_node(
        self, node: Node, websocket: Optional[WebSocket] = None
    ) -> None:
        """Updates a node in the molecule graph.

        If 'node' is a reference to an existing node, this will check
        a few things and possibly send update messages back to the
        frontend.

        If 'node' is not a reference, it must have the same node.id as
        some other node in the graph, which it will replace.

        :param node: The node to update.
        :param websocket: If not None, notifies websocket of addition.

        """

        # Keep the old value alive (this will raise KeyError if
        # node.id is unknown)
        old_node = self.node_ids[node.id]

        # Following the order of add_node, we reset the node in
        # node_ids, then work on parent relationships, etc.
        self.node_ids[node.id] = node

        # Parental relationship management.
        #
        # NOTE (trb): I don't think this actually happens in the code
        # so far, but I'm adding this for consistency (and
        # future-proofing, or whatever).
        if node.parentId != self.parents.get(node.id):
            if node.id in self.parents:
                # Update an existing edge
                edge = next(
                    edge for _, edge in self.edges.items() if edge.toNode == node.id
                )
                edge.fromNode = node.parentId

                # Push update to client
                if websocket is not None:
                    await websocket.send_json(
                        {"type": "edge_update", "edge": self.edges[eid].json()}
                    )

            else:
                # Create a new edge
                edge = Edge(f"edge_{uuid4()}", node.parentId, node.id, "complete")
                self.edges[edge.id] = edge
                if websocket is not None:
                    await websocket.send_json({"type": "edge", "edge": edge.json()})

            # Add the parent relationship
            self.parents[node.id] = node.parentId

        # Recalculate in case anything changed
        self.recalculate_nodes_per_level()

        if websocket is not None:
            await websocket.send_json({"type": "node_update", "node": node.json()})

    async def delete_node(
        self, node_id: str, websocket: Optional[WebSocket] = None
    ) -> None:
        """
        Deletes a node from the context and molecule graph.
        Raises a ``KeyError`` if the node ID is not found in the context.

        :param node_id: The node ID to delete.
        :param websocket: If not None, notifies websocket of deletion.
        """

        # NOTE (trb): The previous implemention did not account for
        # nodes in various states (e.g., it would KeyError out if one
        # added a node and then immediately deleted it). Id didn't
        # consider that the node might not be a leaf and provided no
        # means to "stitch" the graph back together. And even if the
        # graph break is fine, it didn't clean up the 'parentId' in
        # the child nodes. Moveover, this function was not just an
        # implementation detail of `delete_subtree`, which duplicated
        # all of the work that should have been done here anyway. An
        # alternative could be to refactor that function to do a
        # depth-first delete and restore this as the implementation
        # detail handling deletion of a leaf node, but why bother?
        raise NotImplementedError(
            "The previous implementation was incomplete and underspecified."
            "See comments in code."
        )

    def save_state(self) -> dict[str, Any]:
        """Export the class data as a Python dictionary."""

        return TypeAdapter(GraphContext).dump_python(self)

    def load_state(self, data: dict[str, Any]) -> None:
        """Restore a context from a dictionary serialization.
        Existing state is clobbered, even if the dictionary
        serialization is empty or "node"-less.

        This short-circuits after "reset()" if "nodes" is not in the
        dictionary.

        Input data may also contain "edges". If "edges" is NOT
        included, new edges will be generated. These will be "private"
        backend edges that aren't reported to the frontend (no
        websocket here). This is only a problem if they're ever
        "edge-update"-ed. Edges that don't match any node are
        (silently) ignored and not added to the context.

        The implementation just calls add_node() repeatedly,
        rebuilding "parents" and "nodes_per_level" as it goes. Passing
        those explicity in the input data has no effect.

        :param data: Dictionary with serialized nodes and, optionally, edges.

        """

        # Start from a clean slate
        self.reset()

        # Short-circuit an empty object.
        if len(data) == 0 or "nodes" not in data:
            return

        # Grab the list of nodes
        node_data = data.get("nodes", [])

        # Just enough of an ordering to make sure that parents are
        # ordered before children.
        class NodeCmp:
            def __init__(self, node_as_dict: dict[str, Any]):
                self.node_dict = node_as_dict

            def __lt__(self, other):
                return self.node_dict.get("id") == other.node_dict.get("parentId")

        node_data.sort(key=NodeCmp)

        # We explicitly restore edges now because of the UUID. The
        # worldview is node-centric, so any edge that doesn't match a
        # node will just be ignored.
        edges = [Edge(**e) for e in data.get("edges", [])]

        # Deserialize each node. The "yield_"/"yield" fiasco needs to
        # be handled manually, but otherwise, this is pretty
        # automatic -- thanks, pydantic.
        for n in node_data:
            # opting for symmetry with Node.json()
            node = TypeAdapter(Node).validate_python(n)

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

            self._add_node(node, edge)

    def is_empty(self) -> bool:
        return len(self.node_ids) == 0

    def new_node_id(self) -> str:
        """
        Returns a new "unique" (enough, hopefully) node ID.
        """
        i = 0
        while f"node_{i}" in self.node_ids:
            i += 1
        return f"node_{i}"


class FlaskExperiment(Experiment):
    """Experiment subclass that manages a graph context"""

    graph_context: GraphContext

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.graph_context = GraphContext()

    def save_state(self) -> dict[str, Any]:
        state = super().save_state()
        state["graphContext"] = self.graph_context.save_state()
        return state

    def load_state(self, state: dict[str, Any]) -> None:
        # Make sure we're working with a dict
        if isinstance(state, str):
            state = json.loads(state)

        super().load_state(state)

        # If we have an "old-style" serialization, we won't get the
        # "graphContext", likely just "nodes" and "edges".
        if "graphContext" in state:
            # Preferring symmetry with GraphContext.save_state's use
            # of TypeAdapter(GraphContext).dump_python().
            graph_context = TypeAdapter(GraphContext).validate_python(
                state.get("graphContext")
            )
        else:
            self.graph_context.load_state(state)

    def reset(self):
        super().reset()
        self.graph_context.reset()
