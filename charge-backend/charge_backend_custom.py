from fastapi import WebSocket
from loguru import logger
import re

from charge.experiments.AutoGenExperiment import AutoGenExperiment
from charge.tasks.Task import Task
from charge.servers.log_progress import LOG_PROGRESS_SYSTEM_PROMPT
from backend_helper_funcs import Node, CallbackHandler


async def run_custom_problem(
    start_smiles: str,
    system_prompt: str,
    user_prompt: str,
    experiment: AutoGenExperiment,
    available_tools: list[str],
    websocket: WebSocket,
):
    task = Task(
        system_prompt=system_prompt + "\n\n" + LOG_PROGRESS_SYSTEM_PROMPT,
        user_prompt=user_prompt + "\n\nInitial SMILES string: " + start_smiles,
        server_urls=available_tools,
    )
    agent = experiment.create_agent_with_experiment_state(
        task=task,
        callback=CallbackHandler(websocket),
    )

    result = await agent.run()
    await websocket.send_json(
        {
            "type": "response",
            "message": {"source": "Assistant", "message": result},
        }
    )

    # Find all SMILES values
    matches = re.findall(r'"smiles":\s*"([^"]+)"', result)

    if matches:
        for i, smiles in enumerate(matches):
            node = Node(
                id=f"node_{i}",
                smiles=smiles,
                label=smiles,
                hoverInfo=result,
                level=0,
                x=50,
                y=100 + i * 150,
            )
            await websocket.send_json({"type": "node", "node": node.json()})

    await websocket.send_json({"type": "complete"})
