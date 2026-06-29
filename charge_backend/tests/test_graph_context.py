###############################################################################
## Copyright 2025-2026 Lawrence Livermore National Security, LLC.
## See the top-level LICENSE file for details.
##
## SPDX-License-Identifier: Apache-2.0
###############################################################################

import asyncio

from charge_backend.experiment import GraphContext
from charge_backend.backend_helper_funcs import Node, Edge
from typing import Optional


def make_node(node_id: str, level: int, parent_id: Optional[str] = None) -> Node:
    return Node(
        node_id,
        smiles=f"{node_id}_smiles",
        label=f"{node_id}_label",
        hoverInfo=f"{node_id}_hover",
        level=level,
        parentId=parent_id,
    )


def test_adding_isolated_nodes_to_graph_context():
    g = GraphContext()

    assert len(g.node_ids) == 0

    asyncio.run(g.add_node(make_node("a", 0)))
    asyncio.run(g.add_node(make_node("b", 0)))

    assert len(g.node_ids) == 2
    assert len(g.edges) == 0
    assert len(g.parents) == 0
    assert g.nodes_per_level[0] == 2

    g.reset()

    assert len(g.node_ids) == 0
    assert len(g.edges) == 0
    assert len(g.parents) == 0
    assert len(g.nodes_per_level) == 0


def test_adding_connected_nodes_to_graph_context():
    g = GraphContext()

    asyncio.run(g.add_node(make_node("a", 0)))
    asyncio.run(g.add_node(make_node("b", 1, "a")))

    assert len(g.node_ids) == 2
    assert len(g.edges) == 1
    assert len(g.parents) == 1
    assert g.nodes_per_level[0] == 1
    assert g.nodes_per_level[1] == 1

    edge = next((edge for _, edge in g.edges.items()), None)

    assert edge is not None
    assert edge.fromNode == "a"
    assert edge.toNode == "b"

    new_node = make_node("c", 1, "a")
    new_edge = Edge("unique_edge_id", "a", "c", "complete")

    asyncio.run(g.add_node(new_node, edge=new_edge))

    assert len(g.node_ids) == 3
    assert len(g.edges) == 2
    assert len(g.parents) == 2
    assert g.nodes_per_level[0] == 1
    assert g.nodes_per_level[1] == 2
    assert "unique_edge_id" in g.edges

    g.reset()

    assert len(g.node_ids) == 0
    assert len(g.edges) == 0
    assert len(g.parents) == 0
    assert len(g.nodes_per_level) == 0


def test_deleting_subtree():

    g = GraphContext()
    asyncio.run(g.add_node(make_node("a", 0)))
    asyncio.run(g.add_node(make_node("b", 1, "a")))
    asyncio.run(g.add_node(make_node("c", 2, "b")))
    asyncio.run(g.add_node(make_node("d", 2, "b")))
    asyncio.run(g.add_node(make_node("e", 3, "c")))

    assert len(g.node_ids) == 5
    assert len(g.edges) == 4
    assert len(g.parents) == 4
    assert g.nodes_per_level[0] == 1
    assert g.nodes_per_level[1] == 1
    assert g.nodes_per_level[2] == 2
    assert g.nodes_per_level[3] == 1

    asyncio.run(g.delete_subtree("b"))

    assert "b" in g.node_ids  # root is NOT deleted
    assert all([nd_id not in g.node_ids for nd_id in ("c", "d", "e")])
    assert len(g.node_ids) == 2
    assert len(g.edges) == 1
    assert len(g.parents) == 1
    assert g.nodes_per_level[0] == 1
    assert g.nodes_per_level[1] == 1
    assert g.nodes_per_level[2] == 0
    assert g.nodes_per_level[3] == 0


def test_updating_node_with_new_data():
    g = GraphContext()

    node_a = make_node("a", 0)
    asyncio.run(g.add_node(node_a))

    assert g.node_ids["a"].bandgap is None

    # Update by reference. Since there's no change in parent
    # relationship and no websocket available here, this basically
    # just ensures things don't break.
    node_a.bandgap = 13.0
    asyncio.run(g.update_node(node_a))

    assert g.node_ids["a"].bandgap == 13.0

    # Make a new instance to show the update_node function does
    # something rather than just modifying a reference.
    tmp = make_node("a", 0)
    tmp.bandgap = 42.0
    asyncio.run(g.update_node(tmp))

    assert g.node_ids["a"].bandgap == 42.0


def test_updating_node_with_new_parent():
    g = GraphContext()
    asyncio.run(g.add_node(make_node("a", 0)))
    asyncio.run(g.add_node(make_node("b", 0)))

    assert len(g.node_ids) == 2
    assert len(g.edges) == 0
    assert len(g.parents) == 0
    assert g.nodes_per_level[0] == 2

    # Make a new instance to show the update_node function does
    # something rather than just modifying a reference.
    tmp = make_node("b", 1, "a")
    asyncio.run(g.update_node(tmp))

    assert len(g.node_ids) == 2
    assert len(g.edges) == 1
    assert len(g.parents) == 1
    assert g.nodes_per_level[0] == 1
    assert g.nodes_per_level[1] == 1

    asyncio.run(g.add_node(make_node("c", 0)))

    # Update by reference, make sure the parent change is detected.
    tmp.parentId = "c"
    asyncio.run(g.update_node(tmp))

    assert len(g.node_ids) == 3
    assert len(g.edges) == 1
    assert len(g.parents) == 1
    assert g.parents["b"] == "c"
    assert g.nodes_per_level[0] == 2
    assert g.nodes_per_level[1] == 1

    # Update with a new object
    new_b = make_node("b", 1, "a")
    asyncio.run(g.update_node(new_b))

    assert len(g.node_ids) == 3
    assert len(g.edges) == 1
    assert len(g.parents) == 1
    assert g.parents["b"] == "a"
    assert "c" not in g.parents.values()
    assert g.nodes_per_level[0] == 2
    assert g.nodes_per_level[1] == 1
