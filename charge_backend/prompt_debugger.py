"""
Prompt debugger wrapper around a ChARGe agent. Prompts the user to modify
prompt using the given websocket.
"""

import asyncio
from dataclasses import dataclass, asdict
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
    attachments = getattr(task, "attachments", []) or []
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
    if json.get("prompt"):
        task.user_prompt = json["prompt"]
    if "attachments" in json:
        task.attachments = validate_image_attachments(json)
        await websocket.send_json(
            {
                "type": "response",
                "message": {
                    "source": "User",
                    "message": (
                        "Prompt breakpoint approved with attached images"
                        if task.attachments
                        else "Prompt breakpoint approved without attached images"
                    ),
                    "images": image_refs(task.attachments) if task.attachments else {},
                },
            }
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
