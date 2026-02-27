from fastapi import WebSocket
import re

from charge.experiments.experiment import Experiment
from charge.tasks.task import Task
from charge.utils.log_progress import LOG_PROGRESS_SYSTEM_PROMPT
from backend_helper_funcs import Node, CallbackHandler, RunSettings
from moleculedb.molecule_naming import smiles_to_html
from charge_backend.prompt_debugger import debug_prompt


async def run_custom_problem(
    start_smiles: str,
    system_prompt: str,
    user_prompt: str,
    experiment: Experiment,
    available_tools: list[str],
    websocket: WebSocket,
    run_settings: RunSettings,
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

    if run_settings.prompt_debugging:
        await debug_prompt(agent, websocket)
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
                label=smiles_to_html(smiles, run_settings.molecule_name_format),
                hoverInfo=result,
                level=0,
                x=50,
                y=100 + i * 150,
            )
            await websocket.send_json({"type": "node", "node": node.json()})

    await websocket.send_json({"type": "complete"})
