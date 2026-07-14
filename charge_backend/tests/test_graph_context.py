###############################################################################
## Copyright 2025-2026 Lawrence Livermore National Security, LLC.
## See the top-level LICENSE file for details.
##
## SPDX-License-Identifier: Apache-2.0
###############################################################################

import asyncio

from charge_backend.flask_experiment import FlaskExperiment, GraphContext
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


def test_graph_context_save_state():
    g = GraphContext()
    asyncio.run(g.add_node(make_node("a", 0)))
    asyncio.run(g.add_node(make_node("b", 1, "a")))
    asyncio.run(g.add_node(make_node("c", 2, "b")))
    asyncio.run(g.add_node(make_node("d", 2, "b")))
    asyncio.run(g.add_node(make_node("e", 3, "c")))

    g_dict = g.save_state()

    assert "node_ids" in g_dict
    assert len(g_dict["node_ids"]) == 5
    assert "edges" in g_dict
    assert len(g_dict["edges"]) == 4

    # A dumb little test, just to sanity check that the pydantic alias
    # logic for Nodes is propagated through "model_dump" correctly.
    assert "yield" in g_dict["node_ids"]["a"]
    assert "yield_" not in g_dict["node_ids"]["a"]


def test_graph_context_load_state():

    g_dict = {}
    g_empty = GraphContext()

    assert g_empty.is_empty()
    g_empty.load_state(g_dict)
    assert g_empty.is_empty()

    # "nodes" is a list of serialized Node objects
    g_dict["nodes"] = [
        {
            "id": "a",
            "smiles": "a_smiles",
            "label": "a_label",
            "hoverInfo": "a_hover",
            "level": 0,
            "parentId": None,
        },
        {
            "id": "b",
            "smiles": "b_smiles",
            "label": "b_label",
            "hoverInfo": "b_hover",
            "level": 1,
            "parentId": "a",
        },
        {
            "id": "c",
            "smiles": "c_smiles",
            "label": "c_label",
            "hoverInfo": "c_hover",
            "level": 2,
            "parentId": "b",
        },
        {
            "id": "d",
            "smiles": "d_smiles",
            "label": "d_label",
            "hoverInfo": "d_hover",
            "level": 2,
            "parentId": "b",
        },
        {
            "id": "e",
            "smiles": "e_smiles",
            "label": "e_label",
            "hoverInfo": "e_hover",
            "level": 3,
            "parentId": "c",
        },
    ]

    g_from_nodes = GraphContext()
    assert g_from_nodes.is_empty()
    g_from_nodes.load_state(g_dict)

    assert len(g_from_nodes.node_ids) == 5
    assert g_from_nodes.node_ids.keys() == {"a", "b", "c", "d", "e"}
    assert len(g_from_nodes.edges) == 4

    # "edges" is a list of serialized Edge objects
    g_dict["edges"] = [
        {
            "id": "edge_a_b",
            "fromNode": "a",
            "toNode": "b",
            "status": "complete",
            "label": None,
        },
        {
            "id": "edge_b_c",
            "fromNode": "b",
            "toNode": "c",
            "status": "complete",
            "label": None,
        },
        {
            "id": "edge_b_d",
            "fromNode": "b",
            "toNode": "d",
            "status": "complete",
            "label": None,
        },
        {
            "id": "edge_c_e",
            "fromNode": "c",
            "toNode": "e",
            "status": "complete",
            "label": None,
        },
    ]
    g_from_nodes_edges = GraphContext()
    assert g_from_nodes_edges.is_empty()
    g_from_nodes_edges.load_state(g_dict)

    assert len(g_from_nodes_edges.node_ids) == 5
    assert g_from_nodes_edges.node_ids.keys() == {"a", "b", "c", "d", "e"}
    assert len(g_from_nodes_edges.edges) == 4
    assert "edge_a_b" in g_from_nodes_edges.edges
    assert "edge_b_c" in g_from_nodes_edges.edges
    assert "edge_b_d" in g_from_nodes_edges.edges
    assert "edge_c_e" in g_from_nodes_edges.edges
    assert g_from_nodes_edges.edges.keys() == {
        "edge_a_b",
        "edge_b_c",
        "edge_b_d",
        "edge_c_e",
    }


def test_graph_context_in_flask_experiment_load_state():

    g = GraphContext()
    asyncio.run(g.add_node(make_node("a", 0)))
    asyncio.run(g.add_node(make_node("b", 1, "a")))
    asyncio.run(g.add_node(make_node("c", 2, "b")))
    asyncio.run(g.add_node(make_node("d", 2, "b")))
    asyncio.run(g.add_node(make_node("e", 3, "c")))

    experimentContext = {}
    experimentContext["graphContext"] = {}
    experimentContext["graphContext"] = g.save_state()

    e = FlaskExperiment(task=None)
    e.load_state(experimentContext)

    assert not e.graph_context.is_empty()
    assert e.graph_context == g

    # "graphContext" takes precedence over "nodes"
    experimentContext["nodes"] = [
        {
            "id": "f",
            "smiles": "f_smiles",
            "label": "f_label",
            "hoverInfo": "f_hover",
            "level": 0,
            "parentId": None,
        }
    ]
    e.load_state(experimentContext)

    assert not e.graph_context.is_empty()
    assert "f" not in e.graph_context.node_ids
    assert e.graph_context == g


def test_flask_experiment_load_state_preserves_graph_context_identity():
    # load_state must repopulate the EXISTING graph_context in place, not
    # replace the object. Callers hand out a reference (get_retro_synth_context)
    # and then mutate it; swapping the object would orphan those references and
    # cause nodes to be written to a stale graph (seen as duplicated molecules).
    e = FlaskExperiment(task=None)
    handed_out = e.graph_context

    g = GraphContext()
    asyncio.run(g.add_node(make_node("a", 0)))
    e.load_state({"graphContext": g.save_state()})

    assert e.graph_context is handed_out
    assert "a" in handed_out.node_ids

    # The reset() path (old-style / empty) must also preserve identity.
    e.load_state({"graphContext": GraphContext().save_state()})
    assert e.graph_context is handed_out
    assert handed_out.is_empty()
