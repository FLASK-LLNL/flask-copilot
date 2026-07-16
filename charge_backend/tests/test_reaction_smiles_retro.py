###############################################################################
## Copyright 2025-2026 Lawrence Livermore National Security, LLC.
## See the top-level LICENSE file for details.
##
## SPDX-License-Identifier: Apache-2.0
###############################################################################

import asyncio

from charge_backend.flask_experiment import GraphContext
from charge_backend.backend_helper_funcs import FlaskRunSettings
from charge_backend.retrosynthesis.reaction_smiles import (
    reaction_smiles_retrosynthesis,
)
from flask_tools.chemistry.smiles_utils import canonicalize_smiles


class FakeWebSocket:
    """Records send_json calls so tests can inspect what was transmitted."""

    def __init__(self):
        self.messages = []

    async def send_json(self, data):
        self.messages.append(data)


def _run(reaction_smiles: str):
    g = GraphContext()
    ws = FakeWebSocket()
    asyncio.run(
        reaction_smiles_retrosynthesis(reaction_smiles, g, ws, FlaskRunSettings())
    )
    return g, ws


# Node SMILES are the raw components from the input reaction SMILES (split on
# '>' and '.'), so tests compare against the exact strings supplied.


def test_single_product_builds_partial_graph():
    g, ws = _run("BrC1=CC=NC=C1.C=CB(O)O>>C=Cc2ccncc2")

    # One root (product) + two children (reactants). Graph nodes store
    # canonical SMILES, so expected values are canonicalized too.
    assert len(g.node_ids) == 3
    root = g.node_ids["node_0"]
    assert root.level == 0
    assert root.parentId is None
    assert root.smiles == canonicalize_smiles("C=Cc2ccncc2")

    children = [n for nid, n in g.node_ids.items() if g.parents.get(nid) == "node_0"]
    assert len(children) == 2
    assert all(c.level == 1 for c in children)
    child_smiles = {c.smiles for c in children}
    assert child_smiles == {
        canonicalize_smiles("BrC1=CC=NC=C1"),
        canonicalize_smiles("C=CB(O)O"),
    }

    # Root has a reaction with a mapped reaction for hover-highlighting
    assert root.reaction is not None
    assert root.reaction.id == "user_rxn"
    assert root.reaction.mappedReaction is not None

    assert ws.messages[-1] == {"type": "complete"}


def test_multi_product_roots_on_largest():
    # Larger product is the pyridine ring system, not the methane byproduct.
    g, _ = _run("BrC1=CC=NC=C1.C=CB(O)O>>C=Cc2ccncc2.C")
    root = g.node_ids["node_0"]
    assert root.smiles == canonicalize_smiles("C=Cc2ccncc2")
    # Byproduct "C" is dropped; only the two reactants are children.
    assert len(g.node_ids) == 3


def test_reagents_appear_as_labeled_children():
    # reactants>reagent>product form: the Pd reagent becomes a child labeled "(Reagent)".
    g, _ = _run("BrC1=CC=NC=C1.C=CB(O)O>[Pd]>C=Cc2ccncc2")
    root = g.node_ids["node_0"]
    assert root.smiles == canonicalize_smiles("C=Cc2ccncc2")
    children = [n for nid, n in g.node_ids.items() if g.parents.get(nid) == "node_0"]
    # Two reactants + one reagent
    assert len(children) == 3
    child_smiles = {c.smiles for c in children}
    assert child_smiles == {
        canonicalize_smiles("BrC1=CC=NC=C1"),
        canonicalize_smiles("C=CB(O)O"),
        canonicalize_smiles("[Pd]"),
    }
    # The reagent child carries a "(Reagent)" role in its label.
    reagent_child = next(c for c in children if c.smiles == "[Pd]")
    assert "(Reagent)" in reagent_child.label
    assert "# Reagent" in reagent_child.hoverInfo


def test_invalid_reaction_smiles_does_not_crash():
    g, ws = _run("this is not a reaction")
    assert len(g.node_ids) == 0
    assert ws.messages[-1] == {"type": "complete"}


# --- Daylight reaction SMILES spec examples (reactant '>' reagent '>' product) ---


def test_spec_example_no_reagent():
    # "C=CCBr>>C=CCI" -- allyl bromide to allyl iodide, no reagent.
    g, ws = _run("C=CCBr>>C=CCI")

    assert len(g.node_ids) == 2
    root = g.node_ids["node_0"]
    assert root.smiles == "C=CCI"

    children = [n for nid, n in g.node_ids.items() if g.parents.get(nid) == "node_0"]
    assert len(children) == 1
    assert children[0].smiles == "C=CCBr"
    assert ws.messages[-1] == {"type": "complete"}


def test_spec_example_with_ions():
    # "[I-].[Na+].C=CCBr>>[Na+].[Br-].C=CCI" -- same reaction, fully specified.
    # The organic product (most heavy atoms) is the root; the Na+/Br- salt
    # byproducts are dropped. All three reactants become children.
    g, _ = _run("[I-].[Na+].C=CCBr>>[Na+].[Br-].C=CCI")

    root = g.node_ids["node_0"]
    assert root.smiles == "C=CCI"

    children = [n for nid, n in g.node_ids.items() if g.parents.get(nid) == "node_0"]
    assert len(children) == 3
    child_smiles = {c.smiles for c in children}
    assert child_smiles == {"[I-]", "[Na+]", "C=CCBr"}


def test_spec_example_with_reagent():
    # "C=CCBr.[Na+].[I-]>CC(=O)C>C=CCI.[Na+].[Br-]" -- acetone is a reagent
    # (solvent). It is surfaced as a child node labeled "(Reagent)".
    g, _ = _run("C=CCBr.[Na+].[I-]>CC(=O)C>C=CCI.[Na+].[Br-]")

    root = g.node_ids["node_0"]
    assert root.smiles == "C=CCI"

    children = [n for nid, n in g.node_ids.items() if g.parents.get(nid) == "node_0"]
    # Three reactants + the acetone reagent.
    assert len(children) == 4
    child_smiles = {c.smiles for c in children}
    assert child_smiles == {
        canonicalize_smiles("C=CCBr"),
        canonicalize_smiles("[Na+]"),
        canonicalize_smiles("[I-]"),
        canonicalize_smiles("CC(=O)C"),
    }
    # The acetone reagent appears with a "(Reagent)" role label.
    reagent_child = next(
        c for c in children if c.smiles == canonicalize_smiles("CC(=O)C")
    )
    assert "(Reagent)" in reagent_child.label


# --- Multiple reaction SMILES (multi-step chain, one per line) ---


def _children_of(g, node_id):
    return [n for nid, n in g.node_ids.items() if g.parents.get(nid) == node_id]


def test_two_step_chain_expands_matching_leaf():
    # Line 1: build C <- A + B. Line 2: expand A (a leaf) into D + E.
    g, ws = _run("BrC1=CC=NC=C1.C=CB(O)O>>C=Cc2ccncc2\n" "BrBr.c1ccncc1>>BrC1=CC=NC=C1")

    assert len(g.node_ids) == 5
    root = g.node_ids["node_0"]
    assert root.smiles == canonicalize_smiles("C=Cc2ccncc2")

    # The pyridyl bromide leaf was matched and expanded one level deeper.
    matched = next(
        c
        for c in _children_of(g, "node_0")
        if c.smiles == canonicalize_smiles("BrC1=CC=NC=C1")
    )
    grandchildren = _children_of(g, matched.id)
    assert len(grandchildren) == 2
    assert all(gc.level == 2 for gc in grandchildren)
    assert {gc.smiles for gc in grandchildren} == {
        canonicalize_smiles("BrBr"),
        canonicalize_smiles("c1ccncc1"),
    }

    # The expanded (formerly leaf) node now carries a reaction.
    assert matched.reaction is not None
    assert ws.messages[-1] == {"type": "complete"}


def test_tab_separated_reactions():
    g, _ = _run("BrC1=CC=NC=C1.C=CB(O)O>>C=Cc2ccncc2\tBrBr.c1ccncc1>>BrC1=CC=NC=C1")
    assert len(g.node_ids) == 5


def test_unmatched_later_line_is_skipped():
    # Second reaction's product matches no leaf; it is skipped, first intact.
    g, ws = _run("BrC1=CC=NC=C1.C=CB(O)O>>C=Cc2ccncc2\n" "CCO.CCN>>CCCCCCCCCC")
    assert len(g.node_ids) == 3  # only the first reaction's graph
    assert ws.messages[-1] == {"type": "complete"}


def test_bad_first_line_produces_no_graph():
    g, ws = _run("this is not a reaction\nA>>B")
    assert len(g.node_ids) == 0
    assert ws.messages[-1] == {"type": "complete"}


def test_chain_matches_leaf_across_noncanonical_smiles():
    # Line 1 writes the pyridyl reactant aromatic (c1ccncc1); line 2 writes the
    # same molecule Kekule (C1=CC=NC=C1). Matching must canonicalize both sides
    # so the leaf is still found and expanded.
    g, ws = _run("c1ccncc1.C=CB(O)O>>C=Cc2ccncc2\nBrBr>>C1=CC=NC=C1")

    assert len(g.node_ids) == 4
    matched = next(
        c
        for c in _children_of(g, "node_0")
        if canonicalize_smiles(c.smiles) == canonicalize_smiles("c1ccncc1")
    )
    grandchildren = _children_of(g, matched.id)
    assert len(grandchildren) == 1
    assert grandchildren[0].level == 2
    assert canonicalize_smiles(grandchildren[0].smiles) == canonicalize_smiles("BrBr")
    assert ws.messages[-1] == {"type": "complete"}
