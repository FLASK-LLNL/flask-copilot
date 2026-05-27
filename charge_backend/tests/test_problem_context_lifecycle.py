import asyncio
from types import SimpleNamespace

from charge.experiments.experiment import Experiment

from charge_backend.backend_helper_funcs import Node
from charge_backend.backend_manager import FlaskActionManager
from charge_backend.retrosynthesis.context import RetrosynthesisContext


class FakeWebSocket:
    headers = {}

    async def send_json(self, payload):
        self.last_payload = payload


class FakeTaskManager:
    def __init__(self, websocket):
        self.websocket = websocket
        self.configured_tool_servers = []
        self.discovered_local_mcp_tools = {}
        self.selected_tool_runtime = None


def make_manager():
    websocket = FakeWebSocket()
    task_manager = FakeTaskManager(websocket)
    return FlaskActionManager(
        task_manager=task_manager,
        experiment=Experiment(task=None),
        args=SimpleNamespace(),
        username="test-user",
    )


def test_reset_problem_context_clears_retrosynthesis_state_for_lmo():
    manager = make_manager()
    manager.retro_synth_context = RetrosynthesisContext()
    manager.experiment.saved_agent_sessions = {"reaction:node_0": {"memory": "old"}}

    manager.reset_problem_context("optimization")

    assert manager.retro_synth_context is None
    assert "agentSessions" not in manager.experiment.save_state()


def test_load_retrosynthesis_state_replaces_existing_graph():
    async def run() -> None:
        manager = make_manager()
        manager.retro_synth_context = RetrosynthesisContext()
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
