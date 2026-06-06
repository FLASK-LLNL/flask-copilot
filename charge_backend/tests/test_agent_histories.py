import asyncio
from types import SimpleNamespace
from typing import Any

from charge.clients.agent_factory import Agent, AgentRuntimeConfig
from charge.experiments.experiment import AgentRegistryEntry, Experiment
from charge.tasks.task import Task

from charge_backend.backend_helper_funcs import Node, Reaction
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

    def run(self, reasoning_callback=None, **kwargs) -> str:
        return ""

    def load_memory(self, json_str: str) -> None:
        self.memory = json_str

    def save_memory(self) -> str:
        return self.memory

    def get_model_info(self) -> dict[str, Any]:
        return self.model_info


class FakeToolRuntime:
    def task_kwargs(self) -> dict[str, Any]:
        return {}


def make_manager():
    websocket = FakeWebSocket()
    task_manager = FakeTaskManager(websocket)
    manager = FlaskActionManager(
        task_manager=task_manager,
        experiment=Experiment(task=None),
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
    manager.retro_synth_context = RetrosynthesisContext(
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


def test_reaction_chat_task_appends_hover_info_to_existing_prompt():
    manager = make_manager()
    manager.retro_synth_context = RetrosynthesisContext(
        node_ids={
            "product": Node(
                id="product",
                smiles="CCO",
                label="Product",
                hoverInfo="Product hover",
                level=0,
                reaction=Reaction(id="rxn", hoverInfo="Template match: Suzuki"),
            )
        }
    )
    manager.experiment.agent_registry["reaction:product"] = AgentRegistryEntry(
        agent=FakeAgent(
            task=Task(system_prompt="Existing instructions", user_prompt=""),
            memory="",
            model_info={"backend": "dummy", "model": "dummy-model"},
        ),
        runtime_config=AgentRuntimeConfig(),
    )
    manager.selected_tool_runtime = FakeToolRuntime

    task = manager._chat_task_for_agent(
        "reaction:product",
        {"nodeId": "product", "query": "Why this route?"},
        [],
    )
    task_json = task.to_json()

    assert "Existing instructions" in task_json["system_prompt"]
    assert (
        "Reaction hover information:\nTemplate match: Suzuki"
        in task_json["system_prompt"]
    )


def test_reaction_chat_task_uses_hover_info_metadata_fallback():
    manager = make_manager()
    manager.selected_tool_runtime = FakeToolRuntime

    task = manager._chat_task_for_agent(
        "reaction:product",
        {
            "nodeId": "product",
            "query": "What is the yield?",
            "metadata": {"reactionHoverInfo": "Database yield: 75%"},
        },
        [],
    )
    task_json = task.to_json()

    assert (
        "Reaction hover information:\nDatabase yield: 75%" in task_json["system_prompt"]
    )
