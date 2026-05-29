import json
from typing import Any, Optional

from charge.experiments.experiment import Experiment
from charge.tasks.task import Task


def prompt_context_from_task(task: Optional[dict[str, Any]]) -> list[dict[str, str]]:
    if not isinstance(task, dict):
        return []
    context: list[dict[str, str]] = []
    system_prompt = task.get("system_prompt")
    if isinstance(system_prompt, str) and system_prompt.strip():
        context.append({"title": "Instructions", "text": system_prompt})
    user_prompt = task.get("user_prompt")
    if isinstance(user_prompt, str) and user_prompt.strip():
        context.append({"title": "Task prompt", "text": user_prompt})
    return context


def message_label_for_task(
    task: Optional[dict[str, Any]],
    *,
    agent_key: str,
    metadata: dict[str, Any],
) -> Optional[str]:
    if not isinstance(task, dict):
        return None
    class_name = str(task.get("class_name") or "")
    module = str(task.get("module") or "")
    kind = str(metadata.get("kind") or agent_key.split(":", 1)[0])
    if "Retrosynthesis" in class_name or "retrosynthesis" in module:
        return "Retrosynthesis request"
    if class_name == "LMOTask" or ".lmo." in module:
        return "LMO request"
    if kind == "custom":
        return "Custom prompt"
    return None


def normalized_message_role(raw_message: dict[str, Any]) -> str:
    role = str(raw_message.get("role") or raw_message.get("type") or "assistant")
    if role == "message":
        role = "assistant"
    return role if role in {"user", "assistant", "system", "tool"} else "assistant"


def session_messages_from_memory(memory: object) -> list[dict[str, Any]]:
    if not memory:
        return []
    try:
        session = json.loads(memory) if isinstance(memory, str) else memory
    except json.JSONDecodeError:
        return []
    if not isinstance(session, dict):
        return []
    messages = session.get("state", {}).get("in_memory", {}).get("messages", [])
    return messages if isinstance(messages, list) else []


def record_latest_user_message_metadata(
    experiment: Experiment,
    agent_key: str,
    task: Task,
    *,
    label: Optional[str] = None,
    display_text: Optional[str] = None,
) -> None:
    registry_item = experiment.agent_registry.get(agent_key)
    if not registry_item:
        return
    agent = registry_item.get("agent")
    if agent is None or not hasattr(agent, "save_memory"):
        return
    raw_messages = session_messages_from_memory(agent.save_memory())
    latest_user_index = next(
        (
            index
            for index, raw_message in reversed(list(enumerate(raw_messages)))
            if isinstance(raw_message, dict)
            and normalized_message_role(raw_message) == "user"
        ),
        None,
    )
    if latest_user_index is None:
        return

    task_json = task.to_json() if hasattr(task, "to_json") else {}
    prompt_context = prompt_context_from_task(task_json)
    metadata: dict[str, Any] = {}
    if prompt_context:
        metadata["promptContext"] = prompt_context
    message_label = label or message_label_for_task(
        task_json,
        agent_key=agent_key,
        metadata=registry_item.get("metadata", {}),
    )
    if message_label:
        metadata["label"] = message_label
    if isinstance(display_text, str) and display_text.strip():
        metadata["displayText"] = display_text
    if not metadata:
        return

    message_metadata = {
        **(
            registry_item.get("messageMetadata")
            if isinstance(registry_item.get("messageMetadata"), dict)
            else {}
        )
    }
    message_metadata[str(latest_user_index)] = metadata
    registry_item["messageMetadata"] = message_metadata
