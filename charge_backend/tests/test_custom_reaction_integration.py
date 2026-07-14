###############################################################################
## Copyright 2025-2026 Lawrence Livermore National Security, LLC.
## See the top-level LICENSE file for details.
##
## SPDX-License-Identifier: Apache-2.0
###############################################################################

"""Integration test driving the real FlaskActionManager._handle_custom_problem.

Guards against regressions in the message sequence streamed to the browser for
a custom problem whose input is a reaction SMILES -- specifically that BOTH
node and edge messages are emitted so the frontend can wire the graph.
"""

import argparse
import asyncio

from charge_backend.backend_manager import FlaskActionManager
from charge_backend import flask_experiment as fe


class CollectingWebSocket:
    """Fake websocket that records every message sent."""

    headers: dict = {}

    def __init__(self):
        self.messages = []

    async def send_json(self, payload):
        self.messages.append(payload)

    async def send_text(self, payload):
        self.messages.append(payload)


class FakeAgent:
    def __init__(self, result):
        self._result = result

    async def run(self):
        return self._result


def _make_manager():
    websocket = CollectingWebSocket()
    manager = FlaskActionManager(
        websocket=websocket,
        args=argparse.Namespace(
            backend="openai",
            model="gpt-4o-mini",
            config_file=None,
            json_file=None,
            max_retries=1,
        ),
        username="test-user",
    )
    return manager, websocket


def test_custom_reaction_smiles_streams_nodes_and_edges(monkeypatch):
    manager, websocket = _make_manager()

    # Stub the agent runtime so no real backend is needed. The agent echoes the
    # product/reactant SMILES, which must NOT be re-added as nodes.
    agent_reply = 'Analysis: product {"smiles": "C=CCI"} from {"smiles": "C=CCBr"}.'
    monkeypatch.setattr(
        fe.FlaskExperiment,
        "create_agent_with_experiment_state",
        lambda self, *a, **k: FakeAgent(agent_reply),
    )
    monkeypatch.setattr(
        fe.FlaskExperiment, "add_to_context", lambda self, *a, **k: None
    )

    # selected_tool_runtime() needs a wrapped websocket for bearer-token
    # extraction, which is irrelevant here (FakeAgent ignores the runtime).
    from lc_conductor import ToolRuntime

    monkeypatch.setattr(
        manager,
        "selected_tool_runtime",
        lambda: ToolRuntime(bearer_token=None, tools=[]),
    )

    data = {
        "smiles": "C=CCBr.[Na+].[I-]>CC(=O)C>C=CCI.[Na+].[Br-]",
        "problemType": "custom",
        "systemPrompt": "sys",
        "userPrompt": "analyze",
    }

    asyncio.run(manager._handle_custom_problem(data))

    node_msgs = [m for m in websocket.messages if m.get("type") == "node"]
    edge_msgs = [m for m in websocket.messages if m.get("type") == "edge"]

    # Product root + 3 reactants + 1 reagent = 5 nodes.
    node_ids = [m["node"]["id"] for m in node_msgs]
    assert len(node_ids) == 5
    assert len(node_ids) == len(set(node_ids)), "nodes must not be duplicated"

    # Each non-root node must have a corresponding edge from the root, so the
    # frontend can draw the reaction graph.
    edge_pairs = {(m["edge"]["fromNode"], m["edge"]["toNode"]) for m in edge_msgs}
    assert len(edge_msgs) == 4
    for nid in node_ids:
        if nid != "node_0":
            assert ("node_0", nid) in edge_pairs

    # Every edge endpoint must reference a streamed node (no dangling edges).
    streamed = set(node_ids)
    for frm, to in edge_pairs:
        assert frm in streamed and to in streamed

    # Exactly one terminal complete.
    assert sum(1 for m in websocket.messages if m.get("type") == "complete") == 1
