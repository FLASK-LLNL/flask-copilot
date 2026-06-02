from fastapi import WebSocket
import re
from typing import Awaitable, Callable, Optional

from charge.clients.agent_factory import ReasoningCallbackType
from charge.experiments.experiment import Experiment
from charge.tasks.task import Task
from charge_backend.backend_helper_funcs import Node, CallbackHandler, FlaskRunSettings
from charge_backend.moleculedb.molecule_naming import smiles_to_html
from charge_backend.prompt_debugger import debug_prompt
from lc_conductor import ToolRuntime


async def run_custom_problem(
    start_smiles: str,
    system_prompt: str,
    user_prompt: str,
    experiment: Experiment,
    tool_runtime: ToolRuntime,
    websocket: WebSocket,
    run_settings: FlaskRunSettings,
    log_progress: ReasoningCallbackType,
    attachments: Optional[list[dict[str, object]]] = None,
    history_callback: Optional[Callable[[], Awaitable[None]]] = None,
):
    task = Task(
        system_prompt=system_prompt,
        user_prompt=user_prompt + "\n\nInitial SMILES string: " + start_smiles,
        attachments=attachments or [],
        **tool_runtime.task_kwargs(),
    )
    callback_handler = CallbackHandler(
        websocket, agent_key="custom:main", on_agent_update=history_callback
    )
    agent = experiment.create_agent_with_experiment_state(
        task=task,
        agent_key="custom:main",
        callback=callback_handler,
    )

    if run_settings.prompt_debugging:
        await debug_prompt(agent, websocket)
    result = await agent.run(log_progress)
    await callback_handler.drain()
    experiment.add_to_context(agent, task, result)
    if history_callback is not None:
        await history_callback()
    await websocket.send_json(
        {
            "type": "response",
            "message": {
                "source": "Assistant",
                "message": result,
            },
        }
    )

    # Find all SMILES values
    matches = re.findall(r'"smiles":\s*"([^"]+)"', result)

    if matches:
        for i, smiles in enumerate(matches):
            node = Node(
                id=f"node_{i}",
                smiles=smiles,
                label=smiles_to_html(smiles, run_settings.molecule_name_format),
                hoverInfo=result,
                level=0,
                x=50,
                y=100 + i * 150,
            )
            await websocket.send_json({"type": "node", "node": node.json()})

    await websocket.send_json({"type": "complete"})
