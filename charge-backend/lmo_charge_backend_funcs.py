import charge.utils.helper_funcs as lmo_helper_funcs
from fastapi import WebSocket, WebSocketDisconnect
import asyncio
from loguru import logger
import sys
import os
from charge.clients.autogen import AutoGenClient

from charge.experiments.LMOExperiment import (
    LMOExperiment as LeadMoleculeOptimization,
)

# TODO: Convert this to a dataclass
MOLECULE_HOVER_TEMPLATE = """**SMILES:** `{}`\n
## Properties
 - Molecule Weight: {:.3f}
 - **Cost:** {:.2f}
 - **Density:** {:.3f}
 - **SA Score:** {:.3f}"""


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
    logger.info(f"Starting experiment with lead molecule: {lead_molecule_smiles}")
    parent_id = 0
    node_id = 0
    lead_molecule_data = lmo_helper_funcs.post_process_smiles(
        smiles=lead_molecule_smiles, parent_id=parent_id - 1, node_id=node_id
    )

    # Start the db with the lead molecule
    lmo_helper_funcs.save_list_to_json_file(
        data=[lead_molecule_data], file_path=mol_file_path
    )

    # Start the db with the lead molecule
    lmo_helper_funcs.save_list_to_json_file(
        data=[lead_molecule_data], file_path=mol_file_path
    )
    logger.info(f"Storing found molecules in {mol_file_path}")

    # Run the experiment in a loop
    new_molecules = lmo_helper_funcs.get_list_from_json_file(
        file_path=mol_file_path
    )  # Start with known molecules

    leader_hov = MOLECULE_HOVER_TEMPLATE.format(
        lead_molecule_smiles,
        0.0,  # TODO: Add molecule weight calculation
        0.0,  # TODO: Add cost calculation
        lead_molecule_data["density"],
        lead_molecule_data["sascore"],
    )
    node = dict(
        id=f"node_{node_id}",
        smiles=lead_molecule_smiles,
        label=f"{lead_molecule_smiles}",
        # Add property calculations here
        energy=lead_molecule_data["density"],
        level=0,
        cost=lead_molecule_data["sascore"],
        # Not sure what to put here
        hoverInfo=leader_hov,
        x=0,
        y=node_id * 150,
    )

    await websocket.send_json({"type": "node", **node})

    edge_data = {
        "type": "edge",
        "id": f"edge_{0}_{1}",
        "status": "computing",
        "label": "Optimizing",
        "fromNode": {"id": f"node_{0}", "x": 0, "y": -150},
        "toNode": {"id": f"node_{1}", "x": 200, "y": -150},
    }
    await websocket.send_json(edge_data)
    # Generate one node at a time

    mol_data = [lead_molecule_data]

    for i in range(depth):
        if i > 0:
            edge_data = {
                "type": "edge",
                "id": f"edge_{0}_{1}",
                "status": "computing",
                "label": "Optimizing",
                "fromNode": {"id": f"node_{0}", "x": 0, "y": -150},
                "toNode": {"id": f"node_{1}", "x": 200, "y": -150},
            }
        await websocket.send_json(edge_data)
        # Generate new molecule

        iteration = 0
        while iteration < max_iterations:

            try:
                iteration += 1
                results = await lmo_runner.run()
                results = results.as_list()  # Convert to list of strings
                logger.info(f"New molecules generated: {results}")
                processed_mol = lmo_helper_funcs.post_process_smiles(
                    smiles=results[0], parent_id=parent_id, node_id=node_id
                )
                canonical_smiles = processed_mol["smiles"]
                if (
                    canonical_smiles not in new_molecules
                    and canonical_smiles != "Invalid SMILES"
                ):
                    new_molecules.append(canonical_smiles)
                    mol_data.append(processed_mol)
                    lmo_helper_funcs.save_list_to_json_file(
                        data=mol_data, file_path=mol_file_path
                    )
                    logger.info(f"New molecule added: {canonical_smiles}")
                    node_id += 1
                    mol_hov = MOLECULE_HOVER_TEMPLATE.format(
                        canonical_smiles,
                        0.0,  # TODO: Add molecule weight calculation
                        0.0,  # TODO: Add cost calculation
                        processed_mol["density"],
                        processed_mol["sascore"],
                    )
                    node = dict(
                        id=f"node_{node_id}",
                        smiles=canonical_smiles,
                        label=f"{canonical_smiles}",
                        # Add property calculations here
                        energy=processed_mol["density"],
                        level=0,
                        cost=processed_mol["sascore"],
                        # Not sure what to put here
                        hoverInfo=mol_hov,
                        x=0,
                        y=node_id * 150,
                    )

                    await websocket.send_json({"type": "node", **node})

                    edge_data = {
                        "type": "edge",
                        "id": f"edge_{node_id}_{node_id+1}",
                        "status": "computing",
                        "label": "Optimizing",
                        "fromNode": {"id": f"node_{node_id}", "x": 0, "y": -150},
                        "toNode": {"id": f"node_{node_id+1}", "x": 200, "y": -150},
                    }

                    await websocket.send_json(edge_data)
                    experiment = LeadMoleculeOptimization(
                        lead_molecule=canonical_smiles
                    )
                    lmo_runner.experiment_type = experiment
                    parent_id = node_id

                    break  # Exit while loop to proceed to next node
                else:
                    logger.info(f"Duplicate molecule found: {canonical_smiles}")
                    # Continue the while loop to try generating again
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
        edge_data = {
            "type": "edge",
            "id": f"edge_{0}_{1}",
            "status": "computing",
            "label": "Optimizing",
            "fromNode": {"id": f"node_{0}", "x": 0, "y": -150},
            "toNode": {"id": f"node_{1}", "x": 200, "y": -150},
        }
        await websocket.send_json(edge_data)
        # TODO: Compute here!!!
        await asyncio.sleep(0.8)

    edge_data = {
        "type": "edge",
        "id": f"edge_{0}_{1}",
        "status": "Completed",
        "label": "Optimization Complete",
        "fromNode": {"id": f"node_{0}", "x": 0, "y": -150},
        "toNode": {"id": f"node_{1}", "x": 200, "y": -150},
    }

    await websocket.send_json({"type": "complete"})
