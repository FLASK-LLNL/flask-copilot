import asyncio
from types import SimpleNamespace

from charge.clients.agent_factory import AgentRuntimeConfig
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
        agent_key="custom:main",
        agent=SimpleNamespace(
            task=task,
            save_memory=lambda: memory,
            get_model_info=lambda: {"backend": "dummy", "model": "dummy-model"},
        ),
        runtime_config=AgentRuntimeConfig(
            agent_key="custom:main",
            backend="dummy",
            model="dummy-model",
        ),
    )

    asyncio.run(manager.handle_get_agent({"agentKey": "custom:main"}))

    assert manager.websocket.last_payload["type"] == "agent-response"
    agent = manager.websocket.last_payload["agent"]
    assert agent["agentKey"] == "custom:main"
    assert agent["memory"] == memory
    assert agent["task"]["system_prompt"] == "System"
    assert agent["runtimeConfig"]["backend"] == "dummy"
    assert agent["modelInfo"]["model"] == "dummy-model"


def test_list_agents_returns_serialized_agent_records():
    manager = make_manager()
    manager.experiment.agent_registry["custom:main"] = AgentRegistryEntry(
        agent_key="custom:main",
        agent=SimpleNamespace(
            task=None,
            save_memory=lambda: "",
            get_model_info=lambda: {"backend": "dummy", "model": "dummy-model"},
        ),
        runtime_config=AgentRuntimeConfig(
            agent_key="custom:main",
            backend="dummy",
            model="dummy-model",
        ),
    )

    asyncio.run(manager.handle_list_agents({}))

    assert manager.websocket.last_payload["type"] == "list-agents-response"
    assert [
        agent["agentKey"] for agent in manager.websocket.last_payload["agents"]
    ] == ["custom:main"]


def test_get_agent_returns_empty_record_for_missing_agent():
    manager = make_manager()

    asyncio.run(manager.handle_get_agent({"agentKey": "missing:main"}))

    assert manager.websocket.last_payload == {
        "type": "agent-response",
        "agent": {
            "agentKey": "missing:main",
            "memory": "",
            "modelInfo": {},
        },
    }


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
        agent_key="reaction:product",
        agent=SimpleNamespace(
            task=SimpleNamespace(get_system_prompt=lambda: "Existing instructions")
        ),
        runtime_config=AgentRuntimeConfig(
            agent_key="reaction:product",
        ),
    )
    manager.selected_tool_runtime = lambda: SimpleNamespace(task_kwargs=lambda: {})

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
    manager.selected_tool_runtime = lambda: SimpleNamespace(task_kwargs=lambda: {})

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
