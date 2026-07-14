###############################################################################
## Copyright 2025-2026 Lawrence Livermore National Security, LLC.
## See the top-level LICENSE file for details.
##
## SPDX-License-Identifier: Apache-2.0
###############################################################################

import asyncio

from charge_backend.flask_experiment import FlaskExperiment
from charge_backend.backend_helper_funcs import FlaskRunSettings
from charge_backend.charge_backend_custom import run_custom_problem
from charge_backend.retrosynthesis.reaction_smiles import (
    reaction_smiles_retrosynthesis,
)


class FakeWebSocket:
    def __init__(self):
        self.messages = []

    async def send_json(self, data):
        self.messages.append(data)


class FakeAgent:
    """Returns a canned result that contains an extractable SMILES field."""

    def __init__(self, result):
        self._result = result

    async def run(self):
        return self._result


class FakeToolRuntime:
    def task_kwargs(self):
        return {}


class FakeExperiment(FlaskExperiment):
    """FlaskExperiment with the agent runtime stubbed out."""

    def __init__(self, agent_result):
        super().__init__(task=None)
        self._agent_result = agent_result

    def create_agent_with_experiment_state(
        self, task, *, agent_key=None, callback=None
    ):
        return FakeAgent(self._agent_result)

    def add_to_context(self, agent, task, result):
        pass


def _run_custom(reaction_smiles, *, build_first, agent_result='{"smiles": "CCO"}'):
    exp = FakeExperiment(agent_result)
    ws = FakeWebSocket()
    is_reaction = ">" in reaction_smiles

    async def go():
        if build_first and is_reaction:
            await reaction_smiles_retrosynthesis(
                reaction_smiles,
                exp.graph_context,
                ws,
                FlaskRunSettings(),
                send_complete=False,
            )
        await run_custom_problem(
            reaction_smiles,
            "system",
            "user",
            exp,
            FakeToolRuntime(),
            ws,
            FlaskRunSettings(),
            build_nodes_from_response=not (build_first and is_reaction),
        )

    asyncio.run(go())
    return exp, ws


def test_custom_reaction_smiles_builds_partial_graph_without_collision():
    # The custom flow builds the reaction graph first, then runs the agent.
    exp, ws = _run_custom(
        "C=CCBr.[Na+].[I-]>CC(=O)C>C=CCI.[Na+].[Br-]", build_first=True
    )
    g = exp.graph_context

    # Reaction graph: product root + 3 reactants + 1 reagent = 5 nodes.
    assert len(g.node_ids) == 5
    root = g.node_ids["node_0"]
    assert root.parentId is None
    # The agent's "CCO" response node must NOT have been added (no collision).
    assert all(n.smiles != "CCO" for n in g.node_ids.values())
    # Exactly one "complete" (from the agent run, not the graph build).
    assert ws.messages.count({"type": "complete"}) == 1


def test_custom_reaction_smiles_not_duplicated_when_agent_echoes_molecules():
    # Regression: the agent, given the full reaction SMILES, echoes the product
    # and reactant SMILES back in its JSON response. The response-node builder
    # must stay suppressed so those molecules are not added a SECOND time
    # (which previously produced product+reactants duplicated, reagent single).
    rxn = "C=CCBr.[Na+].[I-]>CC(=O)C>C=CCI.[Na+].[Br-]"
    agent_echo = (
        'Product {"smiles": "C=CCI"} from reactants '
        '{"smiles": "C=CCBr"}, {"smiles": "[Na+]"}, {"smiles": "[I-]"}'
    )
    exp, _ = _run_custom(rxn, build_first=True, agent_result=agent_echo)
    g = exp.graph_context

    # Still exactly the 5 reaction-graph nodes; no echoed duplicates.
    assert len(g.node_ids) == 5
    smiles_counts = {}
    for n in g.node_ids.values():
        smiles_counts[n.smiles] = smiles_counts.get(n.smiles, 0) + 1
    assert all(count == 1 for count in smiles_counts.values())


def test_custom_plain_smiles_still_builds_nodes_from_response():
    # Non-reaction input: agent response nodes are still created as before.
    exp, _ = _run_custom("c1ccccc1", build_first=True)
    g = exp.graph_context
    assert "node_0" in g.node_ids
    assert g.node_ids["node_0"].smiles == "CCO"
