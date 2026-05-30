"""
Prompt debugger wrapper around a ChARGe agent. Prompts the user to modify
prompt using the given websocket.
"""

import asyncio
from dataclasses import asdict, dataclass
from charge.clients.agent_factory import Agent
from charge.tasks.task import Task
from fastapi import WebSocket
from typing import Any
from charge_backend.attachments import image_refs, validate_image_attachments


@dataclass
class PromptMetadata:
    temperature: float = -1.0


@dataclass
class PromptMessage:
    prompt: str
    metadata: PromptMetadata | None
    images: dict[str, dict[str, Any]] | None = None
    attachments: list[dict[str, Any]] | None = None

    def json(self):
        result = asdict(self)
        result["type"] = "prompt-breakpoint"
        if result.get("images") is None:
            del result["images"]
        if result.get("attachments") is None:
            del result["attachments"]
        return result


DEBUG_PROMPT_RESPONSES: dict[WebSocket, asyncio.Future] = {}
"""
Dictionary mapping websockets to request/response pairs so that we can support
breakpoints for multiple clients at the same time.
"""


async def _send_prompt_breakpoint_update(
    websocket: WebSocket,
    prompt_changed: bool,
    new_prompt: str,
    attachments_changed: bool,
    new_attachments: list[dict[str, Any]],
) -> None:
    if not prompt_changed and not attachments_changed:
        return

    message_str = "Prompt breakpoint:"
    if prompt_changed:
        message_str += f"\n\nText updated to '{new_prompt}'."

    if attachments_changed:
        message_str += "\n\nAttachments updated"

    message: dict[str, Any] = {"source": "Prompt-point", "message": message_str}
    if attachments_changed:
        message["images"] = image_refs(new_attachments) if new_attachments else {}

    await websocket.send_json({"type": "response", "message": message})


async def debug_prompt_task(task: Task, websocket: WebSocket):
    """
    Prompts the Web UI to approve/edit the prompt (and potential metadata).
    Modifies the runner in place. Should be called before running.

    :param task: ChARGe Task
    :param websocket: Description
    """
    # TODO(later): Find a cleaner way than globals using objects and state
    global DEBUG_PROMPT_RESPONSES
    if websocket in DEBUG_PROMPT_RESPONSES:  # Recursive message waiting
        raise ValueError("Already waiting for a prompt breakpoint.")

    DEBUG_PROMPT_RESPONSES[websocket] = asyncio.Future()

    user_prompt = task.user_prompt or ""
    if hasattr(task, "attachments"):
        attachments = validate_image_attachments({"attachments": task.attachments})
    else:
        attachments = []
    images = image_refs(
        [attachment for attachment in attachments if isinstance(attachment, dict)]
    )
    message = PromptMessage(user_prompt, None, images or None, attachments or None)
    await websocket.send_json(message.json())

    # Wait for response
    json: dict[str, Any] = await DEBUG_PROMPT_RESPONSES[websocket]
    assert json["action"] == "prompt-breakpoint-response"
    del DEBUG_PROMPT_RESPONSES[websocket]

    # Continue with new prompt
    next_prompt = str(json["prompt"]) if "prompt" in json else user_prompt
    prompt_changed = next_prompt != user_prompt
    task.user_prompt = next_prompt

    next_attachments = attachments
    attachments_changed = False
    if "attachments" in json:
        next_attachments = validate_image_attachments(json)
        attachments_changed = next_attachments != attachments
        task.attachments = next_attachments

    await _send_prompt_breakpoint_update(
        websocket,
        prompt_changed,
        next_prompt,
        attachments_changed,
        next_attachments,
    )


async def debug_prompt(runner: Agent, websocket: WebSocket):
    """
    Prompts the Web UI to approve/edit the prompt (and potential metadata).
    Modifies the runner in place. Should be called before running.

    :param runner: A ChARGe agent.
    :param websocket: Websocket to the web UI.
    """
    assert runner.task is not None
    await debug_prompt_task(runner.task, websocket)
