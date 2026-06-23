import asyncio
from types import SimpleNamespace
from typing import Any

from charge.clients.agent_factory import Agent, AgentRuntimeConfig
from charge.experiments.experiment import AgentRegistryEntry
from charge.tasks.task import Task

from charge_backend.backend_helper_funcs import Node, Reaction
from charge_backend.backend_manager import FlaskActionManager
from charge_backend.experiment import FlaskExperiment, GraphContext


class FakeWebSocket:
    headers = {}

    async def send_json(self, payload):
        self.last_payload = payload


class FakeAgent(Agent):
    def __init__(
        self,
        task: Task | None,
        *,
        memory: str,
        model_info: dict[str, Any],
    ):
        super().__init__(task)
        self.memory = memory
        self.model_info = model_info

    def run(self, **kwargs) -> str:
        return ""

    def load_memory(self, json_str: str) -> None:
        self.memory = json_str

    def save_memory(self) -> str:
        return self.memory

    def get_model_info(self) -> dict[str, Any]:
        return self.model_info


def make_manager():
    websocket = FakeWebSocket()
    manager = FlaskActionManager(
        websocket=websocket,
        args=SimpleNamespace(),
        username="test-user",
    )
    return manager


def test_get_agent_returns_serialized_agent_record():
    manager = make_manager()
    memory = '{"state":{"in_memory":{"messages":[{"role":"user"}]}}}'
    task = Task(system_prompt="System", user_prompt="Hello")
    manager.experiment.agent_registry["custom:main"] = AgentRegistryEntry(
        agent=FakeAgent(
            task=task,
            memory=memory,
            model_info={"backend": "dummy", "model": "dummy-model"},
        ),
        runtime_config=AgentRuntimeConfig(
            backend="dummy",
            model="dummy-model",
        ),
    )

    asyncio.run(manager.handle_get_agent({"agentKey": "custom:main"}))

    assert manager.websocket.last_payload["type"] == "agent-response"
    assert manager.websocket.last_payload["agentKey"] == "custom:main"
    agent = manager.websocket.last_payload["agent"]
    assert agent["memory"] == memory
    assert agent["task"]["system_prompt"] == "System"
    assert agent["runtimeConfig"]["backend"] == "dummy"
    assert agent["modelInfo"]["model"] == "dummy-model"


def test_list_agents_returns_serialized_agent_records():
    manager = make_manager()
    manager.experiment.agent_registry["custom:main"] = AgentRegistryEntry(
        agent=FakeAgent(
            task=None,
            memory="",
            model_info={"backend": "dummy", "model": "dummy-model"},
        ),
        runtime_config=AgentRuntimeConfig(
            backend="dummy",
            model="dummy-model",
        ),
    )

    asyncio.run(manager.handle_list_agents({}))

    assert manager.websocket.last_payload["type"] == "list-agents-response"
    assert manager.websocket.last_payload["agents"] == ["custom:main"]


def test_get_agent_returns_empty_record_for_missing_agent():
    manager = make_manager()

    asyncio.run(manager.handle_get_agent({"agentKey": "missing:main"}))

    assert manager.websocket.last_payload == {
        "type": "agent-response",
        "agentKey": "missing:main",
        "agent": {
            "memory": "",
            "modelInfo": {},
        },
    }


def test_get_agent_includes_pending_user_message_from_running_request():
    manager = make_manager()
    agent = FakeAgent(
        task=Task(system_prompt="System", user_prompt="Current question"),
        memory='{"state":{"in_memory":{"messages":[{"role":"user"}]}}}',
        model_info={"backend": "dummy", "model": "dummy-model"},
    )
    agent.pending_user_message = {
        "text": "Current question",
        "afterMessageCount": 1,
    }
    manager.experiment.agent_registry["custom:main"] = AgentRegistryEntry(
        agent=agent,
        runtime_config=AgentRuntimeConfig(
            backend="dummy",
            model="dummy-model",
        ),
    )

    asyncio.run(manager.handle_get_agent({"agentKey": "custom:main"}))

    assert manager.websocket.last_payload["type"] == "agent-response"
    assert manager.websocket.last_payload["agentKey"] == "custom:main"
    assert manager.websocket.last_payload["agent"]["pendingUserMessage"] == {
        "text": "Current question",
        "afterMessageCount": 1,
    }


def test_agent_task_lifecycle_tracks_pending_message_and_instructions():
    agent = FakeAgent(
        task=Task(system_prompt="System", user_prompt="Current question"),
        memory='{"state":{"in_memory":{"messages":[{"role":"user"}]}}}',
        model_info={},
    )

    agent.begin_task_run()
    assert agent.pending_user_message == {
        "text": "Current question",
        "afterMessageCount": 1,
    }

    agent.finish_task_run()
    assert agent.pending_user_message is None
    assert [snapshot.to_json() for snapshot in agent.instruction_history] == [
        {"messageCount": 1, "instructions": "System"}
    ]


def test_reaction_context_includes_reaction_hover_info():
    manager = make_manager()
    manager.experiment.graph_context = GraphContext(
        node_ids={
            "product": Node(
                id="product",
                smiles="CCO",
                label="Product",
                hoverInfo="Product hover",
                level=0,
                reaction=Reaction(
                    id="rxn",
                    hoverInfo="# Literature reaction\nYield: 82%",
                ),
            ),
            "reactant": Node(
                id="reactant",
                smiles="CCBr",
                label="Reactant",
                hoverInfo="Reactant hover",
                level=1,
            ),
        },
        parents={"reactant": "product"},
    )

    context = manager._reaction_context_for_node("product")

    assert "Product: CCO" in context
    assert "Reactants:\nCCBr" in context
    assert "Reaction hover information:\n# Literature reaction\nYield: 82%" in context


def test_reaction_context_uses_hover_info_metadata_fallback():
    manager = make_manager()

    context = manager._reaction_context_for_node(
        "product",
        {"metadata": {"reactionHoverInfo": "Database yield: 75%"}},
    )

    assert "Reaction hover information:\nDatabase yield: 75%" in context


def test_chat_agent_routes_molecule_to_custom_query_handler():
    manager = make_manager()
    routed = []

    async def fake_handle_custom_query_molecule(data):
        routed.append(data)

    manager.handle_custom_query_molecule = fake_handle_custom_query_molecule

    asyncio.run(
        manager.handle_chat_agent(
            {
                "agentKey": "molecule:product",
                "query": "  Why this scaffold?  ",
                "metadata": {"smiles": "CCO"},
            }
        )
    )

    assert routed == [
        {
            "agentKey": "molecule:product",
            "query": "Why this scaffold?",
            "metadata": {"smiles": "CCO"},
            "nodeId": "product",
            "smiles": "CCO",
        }
    ]


def test_chat_agent_routes_reaction_to_custom_query_handler():
    manager = make_manager()
    routed = []

    async def fake_handle_custom_query_reaction(data):
        routed.append(data)

    manager.handle_custom_query_reaction = fake_handle_custom_query_reaction

    asyncio.run(
        manager.handle_chat_agent(
            {
                "agentKey": "reaction:product",
                "query": "What is the yield?",
                "metadata": {"reactionHoverInfo": "Database yield: 75%"},
            }
        )
    )

    assert routed == [
        {
            "agentKey": "reaction:product",
            "query": "What is the yield?",
            "metadata": {"reactionHoverInfo": "Database yield: 75%"},
            "nodeId": "product",
        }
    ]
