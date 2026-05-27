import json
from types import SimpleNamespace

from charge.experiments.experiment import Experiment

from charge_backend.backend_manager import FlaskActionManager


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


def test_agent_history_normalization_hides_debug_fields_by_default():
    manager = make_manager()
    record = {
        "metadata": {"title": "Reaction node_1", "subtitle": "CCO"},
        "modelInfo": {"model": "test-model"},
        "task": {
            "system_prompt": "You are given reaction context: CCO>>CC=O",
            "user_prompt": "What is this?",
        },
        "memory": json.dumps(
            {
                "state": {
                    "in_memory": {
                        "messages": [
                            {
                                "role": "user",
                                "contents": [
                                    {"type": "text", "text": "What is this?"},
                                    {
                                        "type": "data",
                                        "uri": "data:image/png;base64,ZmFrZQ==",
                                        "media_type": "image/png",
                                        "name": "sample.png",
                                    },
                                ],
                            },
                            {
                                "role": "assistant",
                                "contents": [
                                    {
                                        "type": "text_reasoning",
                                        "text": "reasoning summary",
                                        "additional_properties": {
                                            "encrypted_content": "secret"
                                        },
                                    },
                                    {"type": "text", "text": "It is ethanol."},
                                ],
                            },
                        ]
                    }
                }
            }
        ),
    }

    history = manager._normalize_agent_history("reaction:node_1", record)

    assert history["title"] == "Reaction node_1"
    assert "rawSession" not in history
    assert history["messages"][0]["text"] == "What is this?"
    assert history["messages"][0]["context"][0]["title"] == "Instructions"
    assert "CCO>>CC=O" in history["messages"][0]["context"][0]["text"]
    assert history["messages"][0]["images"][0]["dataUrl"].startswith("data:image/png")
    assert history["messages"][1]["reasoning"][0]["text"] == "reasoning summary"
    assert "debug" not in history["messages"][1]["reasoning"][0]


def test_agent_history_normalization_includes_raw_fields_in_debug_mode():
    manager = make_manager()
    record = {
        "memory": json.dumps(
            {
                "state": {
                    "in_memory": {
                        "messages": [
                            {
                                "role": "assistant",
                                "contents": [
                                    {
                                        "type": "text_reasoning",
                                        "text": "reasoning summary",
                                        "additional_properties": {
                                            "encrypted_content": "secret"
                                        },
                                    }
                                ],
                            }
                        ]
                    }
                }
            }
        )
    }

    history = manager._normalize_agent_history("custom:main", record, debug=True)

    assert history["rawSession"]["state"]["in_memory"]["messages"]
    assert history["messages"][0]["raw"]["role"] == "assistant"
    assert (
        history["messages"][0]["reasoning"][0]["debug"]["encrypted_content"] == "secret"
    )


def test_agent_history_reports_estimated_context_usage():
    manager = make_manager()
    record = {
        "modelInfo": {"model": "gpt-5.4"},
        "task": {
            "system_prompt": "You are a careful chemistry assistant.",
            "user_prompt": "Summarize this reaction.",
        },
        "memory": json.dumps(
            {
                "state": {
                    "in_memory": {
                        "messages": [
                            {
                                "role": "user",
                                "contents": [
                                    {"type": "text", "text": "Summarize this reaction."}
                                ],
                            },
                            {
                                "role": "assistant",
                                "contents": [
                                    {
                                        "type": "text",
                                        "text": "It forms the oxidized product.",
                                    }
                                ],
                            },
                        ]
                    }
                }
            }
        ),
    }

    history = manager._normalize_agent_history("reaction:node_1", record)

    assert history["contextUsage"]["usedTokens"] > 0
    assert "maxTokens" not in history["contextUsage"]
    assert "percentUsed" not in history["contextUsage"]
    assert history["contextUsage"]["estimated"] is True
    assert history["contextUsage"]["source"] == "estimate"
    assert history["contextUsage"]["model"] == "gpt-5.4"


def test_agent_history_prefers_provider_context_usage():
    manager = make_manager()
    record = {
        "modelInfo": {
            "model": "gpt-5.4",
            "lastUsage": {
                "inputTokens": 1234,
                "outputTokens": 56,
                "totalTokens": 1290,
            },
        },
        "memory": json.dumps(
            {
                "state": {
                    "in_memory": {
                        "messages": [
                            {
                                "role": "user",
                                "contents": [{"type": "text", "text": "Hello"}],
                            }
                        ]
                    }
                }
            }
        ),
    }

    history = manager._normalize_agent_history("custom:main", record)

    assert history["contextUsage"]["usedTokens"] == 1234
    assert history["contextUsage"]["outputTokens"] == 56
    assert history["contextUsage"]["totalTokens"] == 1290
    assert history["contextUsage"]["estimated"] is False
    assert history["contextUsage"]["source"] == "provider"
    assert "maxTokens" not in history["contextUsage"]


def test_agent_history_uses_per_message_metadata_for_prompt_context_and_label():
    manager = make_manager()
    record = {
        "task": {
            "system_prompt": "Latest task only",
            "user_prompt": "Latest user prompt",
        },
        "messageMetadata": {
            "0": {
                "label": "Retrosynthesis request",
                "promptContext": [
                    {"title": "Instructions", "text": "Original reaction context"}
                ],
            },
            "2": {
                "promptContext": [
                    {"title": "Instructions", "text": "Follow-up chat context"}
                ],
            },
        },
        "memory": json.dumps(
            {
                "state": {
                    "in_memory": {
                        "messages": [
                            {
                                "role": "user",
                                "contents": [{"type": "text", "text": "Compute path"}],
                            },
                            {
                                "role": "assistant",
                                "contents": [{"type": "text", "text": "Done"}],
                            },
                            {
                                "role": "user",
                                "contents": [
                                    {"type": "text", "text": "What just happened?"}
                                ],
                            },
                        ]
                    }
                }
            }
        ),
    }

    history = manager._normalize_agent_history("reaction:node_1", record)

    assert history["messages"][0]["label"] == "Retrosynthesis request"
    assert history["messages"][0]["context"][0]["text"] == "Original reaction context"
    assert history["messages"][2]["context"][0]["text"] == "Follow-up chat context"
    assert "Latest task only" not in history["messages"][0]["context"][0]["text"]
