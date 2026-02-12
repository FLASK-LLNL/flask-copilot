"""
Prompt debugger wrapper around a ChARGe agent. Prompts the user to modify
prompt using the given websocket.
"""

import asyncio
from dataclasses import dataclass, asdict
from charge.clients.autogen import AutoGenAgent
from charge.tasks.Task import Task
from fastapi import WebSocket
from typing import Any


@dataclass
class PromptMetadata:
    temperature: float = -1.0


@dataclass
class PromptMessage:
    prompt: str
    metadata: PromptMetadata | None

    def json(self):
        result = asdict(self)
        result["type"] = "prompt-breakpoint"
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
    message = PromptMessage(user_prompt, None)
    await websocket.send_json(message.json())

    # Wait for response
    json: dict[str, Any] = await DEBUG_PROMPT_RESPONSES[websocket]
    assert json["action"] == "prompt-breakpoint-response"
    del DEBUG_PROMPT_RESPONSES[websocket]

    # Continue with new prompt
    if json.get("prompt"):
        task.user_prompt = json["prompt"]


async def debug_prompt(runner: AutoGenAgent, websocket: WebSocket):
    """
    Prompts the Web UI to approve/edit the prompt (and potential metadata).
    Modifies the runner in place. Should be called before running.

    :param runner: A ChARGe agent.
    :param websocket: Websocket to the web UI.
    """
    await debug_prompt_task(runner.task, websocket)
