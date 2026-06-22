import asyncio
from types import SimpleNamespace

from charge.clients.agent_factory import AgentRuntimeConfig
from charge.experiments.experiment import AgentRegistryEntry, Experiment

from charge_backend.backend_helper_funcs import Node
from charge_backend.backend_manager import FlaskActionManager
from charge_backend.experiment import GraphContext


class FakeWebSocket:
    headers = {}

    async def send_json(self, payload):
        self.last_payload = payload


def make_manager():
    websocket = FakeWebSocket()
    return FlaskActionManager(
        websocket=websocket,
        args=SimpleNamespace(),
        username="test-user",
    )


def test_reset_problem_context_clears_retrosynthesis_state_for_lmo():
    manager = make_manager()
    manager.retro_synth_context = GraphContext()
    manager.experiment.agent_registry["reaction:node_0"] = AgentRegistryEntry(
        agent=SimpleNamespace(task=None, save_memory=lambda: "old"),
        runtime_config=AgentRuntimeConfig(),
    )

    manager.reset_problem_context("optimization")

    assert manager.retro_synth_context is None
    assert "agentSessions" not in manager.experiment.save_state()


def test_load_retrosynthesis_state_replaces_existing_graph():
    async def run() -> None:
        manager = make_manager()
        manager.retro_synth_context = GraphContext()
        await manager.retro_synth_context.add_node(
            Node(
                id="old_node",
                smiles="CCO",
                label="old",
                hoverInfo="old",
                level=0,
            )
        )

        await manager.handle_load_state(
            {
                "problemType": "retrosynthesis",
                "nodes": [
                    Node(
                        id="new_node",
                        smiles="CCN",
                        label="new",
                        hoverInfo="new",
                        level=0,
                    ).json()
                ],
                "edges": [],
            }
        )

        assert list(manager.retro_synth_context.node_ids) == ["new_node"]

    asyncio.run(run())
