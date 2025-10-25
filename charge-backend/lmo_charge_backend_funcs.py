import charge.utils.helper_funcs as lmo_helper_funcs
from fastapi import WebSocket, WebSocketDisconnect
import asyncio
from loguru import logger

from charge.clients.autogen import AutoGenClient
from callback_logger import callback_logger
from charge.experiments.Experiment import Experiment
from charge.experiments.LMOExperiment import (
    LMOExperiment as LeadMoleculeOptimization,
    MoleculeOutputSchema,
)
from charge.experiments.LMOExperiment import MoleculeOutputSchema, SCHEMA_PROMPT

from backend_helper_funcs import Node, Edge
from backend_helper_funcs import get_bandgap, post_process_lmo_smiles, get_price

# TODO: Convert this to a dataclass
MOLECULE_HOVER_TEMPLATE = """**SMILES:** `{smiles}`\n
## Properties
 - **Band Gap:** {bandgap:.2f}
 - **Density:** {density:.3f}
 - **Synthesizability (SA) Score:** {sascore:.3f}"""

with open("prompts/lmo_density_user_prompt.txt", "r") as f:
    DENSITY_USER_PROMPT = f.read()


with open("prompts/lmo_refine_prompt.txt", "r") as f:
    FURTHER_REFINE_PROMPT = f.read()


async def lead_molecule(
    start_smiles: str,
    experiment: LeadMoleculeOptimization,
    lmo_runner: AutoGenClient,
    mol_file_path: str,
    max_iterations: int,
    depth: int,
    websocket: WebSocket,
):
    """Stream positioned nodes and edges"""

    lead_molecule_smiles = start_smiles
    clogger = callback_logger(websocket)

    clogger.info(
        f"Starting experiment with lead molecule: {lead_molecule_smiles}",
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

    # Run the experiment in a loop
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

    await websocket.send_json({"type": "node", **node.json()})

    edge_data = Edge(
        id=f"edge_{node_id}_{node_id+1}",
        fromNode=f"node_{node_id}",
        toNode=f"node_{node_id+1}",
        status="computing",
        label="Optimizing",
    )
    await websocket.send_json({"type": "edge", **edge_data.json()})
    logger.info(f"Sending initial edge: {edge_data}")

    # Generate one node at a time

    mol_data = [lead_molecule_data]

    density_experiment = LeadMoleculeOptimization(
        lead_molecule=lead_molecule_smiles,
        user_prompt=DENSITY_USER_PROMPT.format(lead_molecule_smiles)
        + "\n"
        + SCHEMA_PROMPT,
    )
    lmo_runner.experiment_type = density_experiment

    for i in range(depth):
        logger.info(f"Iteration {i}")

        # Generate new molecule

        iteration = 0
        while iteration < max_iterations:

            try:
                iteration += 1
                results: MoleculeOutputSchema = await lmo_runner.run()
                results = results.as_list()  # Convert to list of strings
                clogger.info(f"New molecules generated: {results}")
                processed_mol = lmo_helper_funcs.post_process_smiles(
                    smiles=results[0], parent_id=parent_id, node_id=node_id
                )
                canonical_smiles = processed_mol["smiles"]

                generated_smiles_list = []
                generated_densities = []

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

                        await websocket.send_json({"type": "node", **node.json()})

                        parent_id = node_id
                    else:
                        logger.info(f"Duplicate molecule found: {canonical_smiles}")
                    # Continue the while loop to try generating again
                if len(generated_smiles_list) > 0:

                    density_experiment = LeadMoleculeOptimization(
                        lead_molecule=lead_molecule_smiles,
                        user_prompt=FURTHER_REFINE_PROMPT.format(
                            ", ".join(map(str, generated_densities)),
                            ", ".join(generated_smiles_list),
                        )
                        + "\n"
                        + SCHEMA_PROMPT,
                    )
                    lmo_runner.experiment_type = density_experiment
            except WebSocketDisconnect:
                logger.info("WebSocket disconnected")
                raise
            except asyncio.CancelledError:
                await websocket.send_json({"type": "stopped"})
                raise  # re-raise so cancellation propagates
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

    await websocket.send_json({"type": "edge_update", **edge_data.json()})
    logger.info(f"Sending initial edge: {edge_data}")

    await websocket.send_json({"type": "complete"})
