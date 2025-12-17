from fastapi import WebSocket
from loguru import logger
import re
from typing import Literal

from charge.experiments.AutoGenExperiment import AutoGenExperiment
from charge.tasks.Task import Task
from charge.servers.log_progress import LOG_PROGRESS_SYSTEM_PROMPT
from backend_helper_funcs import Node, CallbackHandler
from molecule_naming import smiles_to_html


async def run_custom_problem(
    start_smiles: str,
    system_prompt: str,
    user_prompt: str,
    experiment: AutoGenExperiment,
    available_tools: list[str],
    websocket: WebSocket,
    molecule_name_format: Literal["brand", "iupac", "formula", "smiles"] = "brand",
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
                label=smiles_to_html(smiles, molecule_name_format),
                hoverInfo=result,
                level=0,
                x=50,
                y=100 + i * 150,
            )
            await websocket.send_json({"type": "node", "node": node.json()})

    await websocket.send_json({"type": "complete"})
