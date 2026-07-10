###############################################################################
## Copyright 2025-2026 Lawrence Livermore National Security, LLC.
## See the top-level LICENSE file for details.
##
## SPDX-License-Identifier: Apache-2.0
###############################################################################

import asyncio

from rdkit import Chem

from charge_backend.flask_experiment import GraphContext
from charge_backend.backend_helper_funcs import FlaskRunSettings
from charge_backend.retrosynthesis.reaction_smiles import (
    reaction_smiles_retrosynthesis,
)


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


def _canon(smiles: str) -> str:
    return Chem.MolToSmiles(Chem.MolFromSmiles(smiles))


def test_single_product_builds_partial_graph():
    g, ws = _run("BrC1=CC=NC=C1.C=CB(O)O>>C=Cc2ccncc2")

    # One root (product) + two children (reactants)
    assert len(g.node_ids) == 3
    root = g.node_ids["node_0"]
    assert root.level == 0
    assert root.parentId is None
    assert root.smiles == _canon("C=Cc2ccncc2")

    children = [n for nid, n in g.node_ids.items() if g.parents.get(nid) == "node_0"]
    assert len(children) == 2
    assert all(c.level == 1 for c in children)
    child_smiles = {c.smiles for c in children}
    assert child_smiles == {_canon("BrC1=CC=NC=C1"), _canon("C=CB(O)O")}

    # Root has a reaction with a mapped reaction for hover-highlighting
    assert root.reaction is not None
    assert root.reaction.id == "user_rxn"
    assert root.reaction.mappedReaction is not None

    assert ws.messages[-1] == {"type": "complete"}


def test_multi_product_roots_on_largest():
    # Larger product is the pyridine ring system, not the methane byproduct.
    g, _ = _run("BrC1=CC=NC=C1.C=CB(O)O>>C=Cc2ccncc2.C")
    root = g.node_ids["node_0"]
    assert root.smiles == _canon("C=Cc2ccncc2")
    # Byproduct "C" is dropped; only the two reactants are children.
    assert len(g.node_ids) == 3


def test_agents_are_ignored():
    # reactants>agent>product form
    g, _ = _run("BrC1=CC=NC=C1.C=CB(O)O>[Pd]>C=Cc2ccncc2")
    root = g.node_ids["node_0"]
    assert root.smiles == _canon("C=Cc2ccncc2")
    children = [n for nid, n in g.node_ids.items() if g.parents.get(nid) == "node_0"]
    assert len(children) == 2


def test_invalid_reaction_smiles_does_not_crash():
    g, ws = _run("this is not a reaction")
    assert len(g.node_ids) == 0
    assert ws.messages[-1] == {"type": "complete"}
