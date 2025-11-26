import charge.utils.helper_funcs as lmo_helper_funcs
from fastapi import WebSocket, WebSocketDisconnect
import asyncio
from loguru import logger
import sys
import os
from charge.experiments.AutoGenExperiment import AutoGenExperiment
from charge.clients.autogen_utils import chargeConnectionError
from callback_logger import CallbackLogger
from charge.tasks.LMOTask import (
    LMOTask as LeadMoleculeOptimization,
    MoleculeOutputSchema,
)
from typing import Optional
from backend_helper_funcs import (
    Node,
    Edge,
    get_bandgap,
    post_process_lmo_smiles,
    get_price,
    CallbackHandler,
)

# TODO: Convert this to a dataclass
MOLECULE_HOVER_TEMPLATE = """**SMILES:** `{smiles}`\n
## Properties
 - **Band Gap:** {bandgap:.2f}
 - **Density:** {density:.3f}
 - **Synthesizability (SA) Score:** {sascore:.3f}"""

with open("prompts/lmo_user_prompt.txt", "r") as f:
    PROPERTY_USER_PROMPT = f.read()


with open("prompts/lmo_refine_prompt.txt", "r") as f:
    FURTHER_REFINE_PROMPT = f.read()


async def generate_lead_molecule(
    start_smiles: str,
    experiment: AutoGenExperiment,
    mol_file_path: str,
    max_iterations: int,
    depth: int,
    available_tools: list[str],
    websocket: WebSocket,
    property: str = "density",
    condition: str = "greater",
    property_description: str = "molecular density (g/cc)",
    custom_prompt: Optional[str] = None,
) -> None:
    """Generate a lead molecule and stream its progress.
    Args:
        start_smiles (str): The starting SMILES string for the lead molecule.
        experiment (AutoGenExperiment): The experiment instance to run tasks.
        mol_file_path (str): Path to the file where molecules are stored.
        max_iterations (int): Maximum iterations for molecule generation.
        available_tools (list[str]): List of available tools for molecule generation.
        depth (int): Depth of the generation tree.
        websocket (WebSocket): WebSocket connection for streaming updates.

    Returns:
        None
    """

    lead_molecule_smiles = start_smiles
    clogger = CallbackLogger(websocket)

    clogger.info(
        f"Starting task with lead molecule: {lead_molecule_smiles}",
        smiles=lead_molecule_smiles,
    )

    parent_id = 0
    node_id = 0
    lead_molecule_data = post_process_lmo_smiles(
        smiles=lead_molecule_smiles, parent_id=parent_id - 1, node_id=node_id
    )

    # Start the db with the lead molecule
    lmo_helper_funcs.save_list_to_json_file(
        data=[lead_molecule_data], file_path=mol_file_path
    )

    clogger.info(f"Storing found molecules in {mol_file_path}")

    # Run the task in a loop
    new_molecules = lmo_helper_funcs.get_list_from_json_file(
        file_path=mol_file_path
    )  # Start with known molecules

    leader_hov = MOLECULE_HOVER_TEMPLATE.format(
        smiles=lead_molecule_smiles,
        bandgap=get_bandgap(lead_molecule_smiles),
        density=lead_molecule_data["density"],
        sascore=lead_molecule_data["sascore"],
    )
    node = Node(
        id=f"node_{node_id}",
        smiles=lead_molecule_smiles,
        label=f"{lead_molecule_smiles}",
        # Add property calculations here
        hoverInfo=leader_hov,
        level=0,
        x=50,
        y=100,
    )
    logger.info(f"Sending root node: {node}")

    await websocket.send_json({"type": "node", "node": node.json()})

    edge_data = Edge(
        id=f"edge_{node_id}_{node_id+1}",
        fromNode=f"node_{node_id}",
        toNode=f"node_{node_id+1}",
        status="computing",
        label="Optimizing",
    )
    await websocket.send_json({"type": "edge", "edge": edge_data.json()})
    logger.info(f"Sending initial edge: {edge_data}")

    # Generate one node at a time

    mol_data = [lead_molecule_data]

    direction = "higher" if condition == "greater" else "lower"
    ranking = "highest" if condition == "greater" else "lowest"

    formatted_user_prompt = (
        custom_prompt
        if custom_prompt is not None
        else PROPERTY_USER_PROMPT.format(
            smiles=lead_molecule_smiles,
            property=property,
            property_description=property_description,
            condition=condition,
            direction=direction,
            ranking=ranking,
        )
    )

    canonical_smiles = lead_molecule_smiles
    callback = CallbackHandler(websocket)
    lmo_task = LeadMoleculeOptimization(
        lead_molecule=lead_molecule_smiles,
        user_prompt=formatted_user_prompt + "\n",
        server_urls=available_tools,
    )

    if os.getenv("CHARGE_DISABLE_OUTPUT_VALIDATION", "0") == "1":
        lmo_task.structured_output_schema = None
        logger.warning("Structure validation disabled for LMOTask output schema.")

    experiment.add_task(lmo_task)

    generated_smiles_list = []
    generated_densities = []

    for i in range(depth):
        logger.info(f"Iteration {i}")

        # Generate new molecule

        iteration = 0
        while iteration < max_iterations:

            try:
                iteration += 1

                await experiment.run_async(callback=callback)
                finished_tasks = experiment.get_finished_tasks()
                _, results = finished_tasks[-1]

                if os.getenv("CHARGE_DISABLE_OUTPUT_VALIDATION", "0") == "1":
                    logger.warning(
                        "Structure validation disabled for LMOTask output schema."
                        "Returning text results without validation first before post-processing."
                    )
                    clogger.info(f"Results: {results}")
                results = MoleculeOutputSchema.model_validate_json(results)
                results = results.as_list()  # Convert to list of strings
                clogger.info(f"New molecules generated: {results}")

                logger.info(f"New molecules generated: {results}")
                for result in results:
                    processed_mol = post_process_lmo_smiles(
                        smiles=result, parent_id=parent_id, node_id=node_id
                    )
                    canonical_smiles = processed_mol["smiles"]
                    densities = processed_mol["density"]
                    if (
                        canonical_smiles not in new_molecules
                        and canonical_smiles != "Invalid SMILES"
                    ):
                        new_molecules.append(canonical_smiles)

                        generated_smiles_list.append(canonical_smiles)
                        generated_densities.append(densities)

                        mol_data.append(processed_mol)
                        lmo_helper_funcs.save_list_to_json_file(
                            data=mol_data, file_path=mol_file_path
                        )
                        logger.info(f"New molecule added: {canonical_smiles}")
                        node_id += 1
                        mol_hov = MOLECULE_HOVER_TEMPLATE.format(
                            smiles=canonical_smiles,
                            bandgap=0.0,  # TODO: Add cost calculation
                            density=processed_mol["density"],
                            sascore=processed_mol["sascore"],
                        )
                        node = Node(
                            id=f"node_{node_id}",
                            smiles=canonical_smiles,
                            label=f"{canonical_smiles}",
                            # Add property calculations here
                            density=processed_mol["density"],
                            bandgap=get_bandgap(canonical_smiles),
                            yield_=None,
                            level=i + 1,
                            cost=get_price(canonical_smiles),
                            # Not sure what to put here
                            hoverInfo=mol_hov,
                            x=150 + node_id * 300,
                            y=100,
                        )

                        await websocket.send_json({"type": "node", "node": node.json()})

                    break  # Exit while loop to proceed to next node
                else:
                    logger.info(f"Duplicate molecule found: {canonical_smiles}")
                    # Continue the while loop to try generating again
                if len(generated_smiles_list) > 0:

                    formatted_refine_prompt = FURTHER_REFINE_PROMPT.format(
                        previous_values=", ".join(map(str, generated_densities)),
                        previous_smiles=", ".join(generated_smiles_list),
                        property=property,
                        property_description=property_description,
                        condition=condition,
                        direction=direction,
                        ranking=ranking,
                    )

                    formatted_refine_prompt = (
                        custom_prompt
                        if custom_prompt is not None
                        else formatted_refine_prompt
                    )

                    lmo_task = LeadMoleculeOptimization(
                        lead_molecule=lead_molecule_smiles,
                        user_prompt=formatted_refine_prompt + "\n",
                        server_urls=available_tools,
                    )

                    if os.getenv("CHARGE_DISABLE_OUTPUT_VALIDATION", "0") == "1":
                        lmo_task.structured_output_schema = None
                        logger.warning(
                            "Structure validation disabled for LMOTask output schema."
                        )
                    experiment.add_task(lmo_task)
                else:
                    logger.info("No new molecules generated in this iteration.")
                    # Re-use the same task
                    experiment.add_task(lmo_task)
                    parent_id = node_id

            except WebSocketDisconnect:
                logger.info("WebSocket disconnected")
                raise
            except asyncio.CancelledError:
                await websocket.send_json({"type": "stopped"})
                raise  # re-raise so cancellation propagates
            except chargeConnectionError as e:
                logger.error(f"Charge connection error: {e}")
                await websocket.send_json({"type": "stopped"})

                raise  # re-raise so higher-level handler can deal with it
            except IndexError:
                logger.error("No finished tasks found.")
            except Exception as e:
                logger.error(f"Error occurred: {e}")
        if i == depth - 1:
            break

        # TODO: Compute here!!!
        await asyncio.sleep(0.8)
    edge_data = Edge(
        id=f"edge_{0}_{1}",
        fromNode=f"node_{0}",
        toNode=f"node_{1}",
        status="complete",
        label="Completed",
    )

    await websocket.send_json({"type": "edge_update", "edge": edge_data.json()})
    logger.info(f"Sending initial edge: {edge_data}")

    await websocket.send_json({"type": "complete"})
