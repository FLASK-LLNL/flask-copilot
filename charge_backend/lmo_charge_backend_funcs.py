import charge.utils.helper_funcs as lmo_helper_funcs
from fastapi import WebSocket, WebSocketDisconnect
import asyncio
from loguru import logger
import sys
import os
from pathlib import Path
from charge.experiments.AutoGenExperiment import AutoGenExperiment
from charge.clients.autogen_utils import chargeConnectionError
from charge.utils.mcp_workbench_utils import call_mcp_tool_directly
from charge_backend.prompt_debugger import debug_prompt_task
from callback_logger import CallbackLogger
from charge.tasks.LMOTask import (
    LMOTask as LeadMoleculeOptimization,
    MoleculeOutputSchema,
)
from typing import Literal, Optional
from backend_helper_funcs import (
    Node,
    Edge,
    get_bandgap,
    post_process_lmo_smiles,
    get_price,
    CallbackHandler,
    RunSettings,
)
from moleculedb.molecule_naming import smiles_to_html, MolNameFormat

# TODO: Convert this to a dataclass
MOLECULE_HOVER_TEMPLATE = (
    "**SMILES:** `{smiles}`\n"
    "\n"
    "## Properties\n"
    " - **Band Gap:** {bandgap:.2f}\n"
    " - **Density:** {density:.3f}\n"
    " - **Synthesizability (SA) Score:** {sascore:.3f}"
)

BASE_DIR = Path(__file__).resolve().parent
PROMPTS_DIR = BASE_DIR / "prompts"

with open(PROMPTS_DIR / "lmo_user_prompt.txt", "r") as f:
    PROPERTY_USER_PROMPT = f.read()


with open(PROMPTS_DIR / "lmo_refine_prompt.txt", "r") as f:
    FURTHER_REFINE_PROMPT = f.read()


async def generate_lead_molecule(
    start_smiles: str,
    experiment: AutoGenExperiment,
    mol_file_path: str,
    max_iterations: int,
    depth: int,
    available_tools: list[str],
    websocket: WebSocket,
    run_settings: RunSettings,
    property: str = "density",
    property_description: str = "molecular density (g/cc)",
    calculate_property_tool: str = "calculate_property_hf",
    condition: str = "greater",
    custom_prompt: Optional[str] = None,
    initial_level: int = 0,
    initial_node_id: int = 0,
    initial_x_position: int = 50,
    enable_constraints: bool = False,
    molecular_similarity: float = 0.7,
    diversity_penalty: float = 0.0,
    exploration_rate: float = 0.5,
    additional_constraints: Optional[list[str]] = None,
    number_of_molecules: int = 10,
    num_top_candidates: int = 3,
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
        property (str): Property to optimize.
        property_description (str): Description of the property.
        calculate_property_tool (str): Tool name for property calculation.
        condition (str): Optimization direction ("greater" or "less").
        custom_prompt (Optional[str]): Custom prompt for the LLM.
        initial_level (int): Initial tree level.
        initial_node_id (int): Initial node ID.
        initial_x_position (int): Initial x-position for visualization.
        run_settings: Default run settings such as molecule display format.
        enable_constraints (bool): Whether to use custom optimization strategy.
        molecular_similarity (float): Similarity threshold to parent molecule (0.0-1.0).
        diversity_penalty (float): Penalty for generating similar molecules (0.0-1.0).
        exploration_rate (float): Balance between exploitation and exploration (0.0-1.0).
        additional_constraints (list[str]): List of constraint types to apply.
        number_of_molecules (int): Number of candidate molecules to generate per iteration.
        num_top_candidates (int): Number of valid new molecules to find per depth level.
    Returns:
        None
    """

    lead_molecule_smiles = start_smiles
    clogger = CallbackLogger(websocket)

    # Only log customization if enabled
    if enable_constraints:
        constraints_str = (
            ", ".join(additional_constraints) if additional_constraints else "none"
        )
        await clogger.info(
            f"Starting optimization with custom strategy:\n"
            f"  - Lead molecule: {lead_molecule_smiles}\n"
            f"  - Molecular similarity threshold: {molecular_similarity}\n"
            f"  - Diversity penalty: {diversity_penalty}\n"
            f"  - Exploration rate: {exploration_rate}\n"
            f"  - Additional constraints: {constraints_str}\n"
            f"  - Available tools: {available_tools}",
            smiles=lead_molecule_smiles,
            source="generate_lead_molecule",
        )
    else:
        await clogger.info(
            f"Starting optimization with default strategy for lead molecule: {lead_molecule_smiles} using available tools: {available_tools}",
            smiles=lead_molecule_smiles,
            source="generate_lead_molecule",
        )

    if additional_constraints is None:
        additional_constraints = []

    parent_id = 0
    node_id = initial_node_id

    # Fix how the available tools interacts with the calculate property tool field
    property_result_msg = await call_mcp_tool_directly(
        tool_name=calculate_property_tool,
        arguments={
            "smiles": lead_molecule_smiles,
            "property": property,
        },
        urls=available_tools,
    )
    results = property_result_msg.result
    if len(results) > 1:
        property_result = float(property_result_msg.result[1].content)
    else:
        property_result = 0.0
        raise ValueError(f"{property_result_msg.result[0].content}")

    lead_molecule_data = post_process_lmo_smiles(
        smiles=lead_molecule_smiles,
        parent_id=parent_id - 1,
        node_id=node_id,
        tool_properties={property: property_result},
    )

    # Start the db with the lead molecule
    lmo_helper_funcs.save_list_to_json_file(
        data=[lead_molecule_data], file_path=mol_file_path
    )

    await clogger.info(f"Storing found molecules in {mol_file_path}")

    # Run the task in a loop
    new_molecules = lmo_helper_funcs.get_list_from_json_file(
        file_path=mol_file_path
    )  # Start with known molecules

    leader_hov = MOLECULE_HOVER_TEMPLATE.format(
        smiles=lead_molecule_smiles,
        bandgap=lead_molecule_data["bandgap"],
        density=lead_molecule_data["density"],
        sascore=lead_molecule_data["sascore"],
    )
    if initial_level == 0:
        node = Node(
            id=f"node_{node_id}",
            smiles=lead_molecule_smiles,
            label=smiles_to_html(
                lead_molecule_smiles, run_settings.molecule_name_format
            ),
            # Add property calculations here
            density=lead_molecule_data.get("density", None),
            sascore=lead_molecule_data.get("sascore", None),
            bandgap=lead_molecule_data.get("bandgap", None),
            hoverInfo=leader_hov,
            level=initial_level,
            x=50,
            y=100,
        )
        logger.info(f"Sending root node: {node}")

        await websocket.send_json({"type": "node", "node": node.json()})

    edge_data = Edge(
        id=f"edge_{node_id}_{node_id + 1}",
        fromNode=f"node_{node_id}",
        toNode=f"node_{node_id + 1}",
        status="computing",
        label="Optimizing",
    )
    await websocket.send_json({"type": "edge", "edge": edge_data.json()})
    logger.info(f"Sending initial edge: {edge_data}")

    # Generate one node at a time

    mol_data = [lead_molecule_data]

    direction = "higher" if condition == "greater" else "lower"
    ranking = "highest" if condition == "greater" else "lowest"

    # Only apply customization guidance if enabled
    customization_text = ""
    if enable_constraints:
        # Build customization guidance for the prompt
        customization_guidance = []

        if molecular_similarity < 0.5:
            customization_guidance.append(
                f"You should explore diverse chemical modifications, "
                f"as the molecular similarity threshold is low ({molecular_similarity:.2f})."
            )
        elif molecular_similarity > 0.8:
            customization_guidance.append(
                f"You should make conservative modifications, "
                f"keeping molecules very similar to the parent (similarity threshold: {molecular_similarity:.2f})."
            )

        if diversity_penalty > 0.5:
            customization_guidance.append(
                f"Prioritize generating chemically diverse molecules to explore different regions "
                f"of chemical space (diversity penalty: {diversity_penalty:.2f})."
            )

        if exploration_rate > 0.7:
            customization_guidance.append(
                f"Focus on exploration - try novel structural modifications "
                f"(exploration rate: {exploration_rate:.2f})."
            )
        elif exploration_rate < 0.3:
            customization_guidance.append(
                f"Focus on exploitation - make incremental improvements to known good structures "
                f"(exploration rate: {exploration_rate:.2f})."
            )

        # Handle additional constraints
        if additional_constraints:
            constraint_guidance = []

            if "drug-likeness" in additional_constraints:
                constraint_guidance.append(
                    "Ensure all molecules satisfy drug-likeness criteria (Lipinski's Rule of Five)"
                )

            if "synthesizability" in additional_constraints:
                constraint_guidance.append(
                    "Prioritize molecules with high synthetic accessibility scores (SA score < 3)"
                )

            if "lead-likeness" in additional_constraints:
                constraint_guidance.append(
                    "Apply lead-likeness criteria suitable for early drug discovery "
                    "(MW 200-350, LogP 1-3)"
                )

            if "pan-assay-interference" in additional_constraints:
                constraint_guidance.append(
                    "Filter out Pan-Assay Interference Compounds (PAINS) and other promiscuous binders"
                )

            if "toxicity-rules" in additional_constraints:
                constraint_guidance.append(
                    "Apply structural alerts to avoid potential toxicity issues"
                )

            if "reactive-groups" in additional_constraints:
                constraint_guidance.append(
                    "Avoid molecules with highly reactive functional groups (epoxides, acyl halides, etc.)"
                )

            if constraint_guidance:
                customization_guidance.extend(constraint_guidance)

        customization_text = (
            "\n\nOptimization Strategy:\n"
            + "\n".join(f"- {g}" for g in customization_guidance)
            + "\n"
            if customization_guidance
            else ""
        )

    # TODO: Use refinement prompt if initial_level != 0 and custom_prompt is None?
    formatted_user_prompt = (
        custom_prompt
        if custom_prompt is not None
        else PROPERTY_USER_PROMPT.format(
            smiles=lead_molecule_smiles,
            property=property,
            property_description=property_description,
            calculate_property_tool=calculate_property_tool,
            condition=condition,
            direction=direction,
            ranking=ranking,
            number_of_molecules=number_of_molecules,
        )
        + customization_text
    )

    canonical_smiles = lead_molecule_smiles
    current_best_smiles = lead_molecule_smiles  # Track the best molecule
    # Initialize best value based on the actual property being optimized
    current_best_value = lead_molecule_data.get(property, 0.0)
    callback = CallbackHandler(websocket)
    lmo_task = LeadMoleculeOptimization(
        lead_molecule=lead_molecule_smiles,
        user_prompt=formatted_user_prompt + "\n",
        property_tool_name=calculate_property_tool,
        property_name=property,
        server_urls=available_tools,
    )
    await lmo_task.get_initial_property_value()

    if os.getenv("CHARGE_DISABLE_OUTPUT_VALIDATION", "0") == "1":
        lmo_task.structured_output_schema = None
        logger.warning("Structure validation disabled for LMOTask output schema.")

    experiment.add_task(lmo_task)

    generated_smiles_list = []
    generated_properties = []

    # Track generated molecules for diversity calculation
    all_generated_smiles = {lead_molecule_smiles}

    # Determine comparison function based on optimization direction
    is_better = (
        (lambda new, old: new > old)
        if condition == "greater"
        else (lambda new, old: new < old)
    )

    for i in range(depth):
        await websocket.send_json(
            {
                "type": "response",
                "message": {
                    "source": "Logger (Info)",
                    "message": f"Iteration {i+1}/{depth}",
                },
            }
        )
        logger.info(f"Iteration {i}")

        # Generate new molecule

        iteration = 0
        iteration_found_better = False
        while iteration < max_iterations:
            try:
                iteration += 1

                if run_settings.prompt_debugging:
                    await debug_prompt_task(lmo_task, websocket)
                await experiment.run_async(callback=callback)
                finished_tasks = experiment.get_finished_tasks()
                _, results = finished_tasks[-1]

                if os.getenv("CHARGE_DISABLE_OUTPUT_VALIDATION", "0") == "1":
                    logger.warning(
                        "Structure validation disabled for LMOTask output schema."
                        "Returning text results without validation first before post-processing."
                    )
                    await clogger.info(f"Results: {results}")
                logger.debug(f"Received: {_}, {results}")
                results = MoleculeOutputSchema.model_validate_json(results)
                reasoning_summary = results.reasoning_summary
                await websocket.send_json(
                    {
                        "type": "response",
                        "message": {
                            "source": "Reasoning",
                            "message": reasoning_summary,
                        },
                    }
                )

                optimized_property = results.property_name
                results = results.as_list()
                for smiles, property_result in results:
                    # Extract properties from tool output if available
                    processed_mol = post_process_lmo_smiles(
                        smiles=smiles,
                        parent_id=parent_id,
                        node_id=node_id,
                        tool_properties={optimized_property: property_result},
                    )
                    canonical_smiles = processed_mol["smiles"]

                    # Apply diversity penalty if enabled
                    if (
                        enable_constraints
                        and diversity_penalty > 0
                        and canonical_smiles in all_generated_smiles
                    ):
                        await clogger.info(
                            f"Molecule {canonical_smiles} already generated - applying diversity penalty"
                        )
                        # Penalize the property value based on diversity_penalty
                        if condition == "greater":
                            property_result *= 1 - diversity_penalty
                        else:
                            property_result *= 1 + diversity_penalty
                        processed_mol[optimized_property] = property_result

                    if (
                        canonical_smiles not in new_molecules
                        and canonical_smiles != "Invalid SMILES"
                    ):
                        new_molecules.append(canonical_smiles)

                        generated_smiles_list.append(canonical_smiles)
                        property_value = processed_mol.get(optimized_property, 0.0)
                        generated_properties.append(property_value)
                        all_generated_smiles.add(canonical_smiles)

                        mol_data.append(processed_mol)
                        lmo_helper_funcs.save_list_to_json_file(
                            data=mol_data, file_path=mol_file_path
                        )
                        logger.info(f"New molecule added: {canonical_smiles}")
                        # Check if this is better than current best
                        if is_better(property_value, current_best_value):
                            current_best_smiles = canonical_smiles
                            current_best_value = property_value
                            iteration_found_better = True
                            await clogger.info(
                                f"New best molecule found: {canonical_smiles} with {property}={property_value}"
                            )

                        node_id += 1
                        mol_hov = MOLECULE_HOVER_TEMPLATE.format(
                            smiles=canonical_smiles,
                            bandgap=processed_mol["bandgap"],
                            density=processed_mol["density"],
                            sascore=processed_mol["sascore"],
                        )
                        node = Node(
                            id=f"node_{node_id}",
                            smiles=canonical_smiles,
                            label=smiles_to_html(
                                canonical_smiles, run_settings.molecule_name_format
                            ),
                            # Add property calculations here
                            density=processed_mol.get("density", None),
                            sascore=processed_mol.get("sascore", None),
                            bandgap=processed_mol.get("bandgap", None),
                            yield_=None,
                            level=initial_level + i + 1,
                            cost=get_price(canonical_smiles),
                            # Highlight if this is the new best
                            highlight=(
                                "yellow"
                                if canonical_smiles == current_best_smiles
                                else "normal"
                            ),
                            # Not sure what to put here
                            hoverInfo=mol_hov,
                            x=initial_x_position
                            + 100
                            + (node_id - initial_node_id) * 300,
                            y=100,
                        )

                        await websocket.send_json({"type": "node", "node": node.json()})

                    else:
                        await clogger.info(
                            f"Duplicate molecule found: {canonical_smiles}"
                        )

                # If we found molecules in this iteration, prepare for next iteration
                if len(generated_smiles_list) > 0:
                    # Break out of max_iterations loop if we found a better molecule
                    if iteration_found_better:
                        await clogger.info(
                            f"Found better molecule in iteration {iteration}, moving to next depth level"
                        )
                        break

                    formatted_refine_prompt = (
                        FURTHER_REFINE_PROMPT.format(
                            previous_values=", ".join(map(str, generated_properties)),
                            previous_smiles=", ".join(generated_smiles_list),
                            property=property,
                            property_description=property_description,
                            calculate_property_tool=calculate_property_tool,
                            condition=condition,
                            direction=direction,
                            ranking=ranking,
                            number_of_molecules=number_of_molecules,
                        )
                        + customization_text
                    )

                    formatted_refine_prompt = (
                        custom_prompt
                        if custom_prompt is not None
                        else formatted_refine_prompt
                    )

                    # Use the BEST molecule found so far as the lead molecule
                    lmo_task = LeadMoleculeOptimization(
                        lead_molecule=current_best_smiles,
                        # This needs to be fixed so that it takes a property name, value, and function for evaluating it
                        user_prompt=formatted_refine_prompt + "\n",
                        property_tool_name=calculate_property_tool,
                        property_name=property,
                        server_urls=available_tools,
                    )
                    await lmo_task.get_initial_property_value()

                    if os.getenv("CHARGE_DISABLE_OUTPUT_VALIDATION", "0") == "1":
                        lmo_task.structured_output_schema = None
                        logger.warning(
                            "Structure validation disabled for LMOTask output schema."
                        )
                    experiment.add_task(lmo_task)
                else:
                    logger.info("No new molecules generated in this iteration.")
                    # Re-use the same task but with current best molecule
                    lmo_task.lead_molecule = current_best_smiles
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

        # Send progress update about current best
        await websocket.send_json(
            {
                "type": "response",
                "message": {
                    "source": "System",
                    "message": f"Iteration {i + 1}/{depth} complete. Best molecule so far: {current_best_smiles} with {property}={current_best_value:.3f}",
                    "smiles": current_best_smiles,
                },
            }
        )

    # Final summary
    await websocket.send_json(
        {
            "type": "response",
            "message": {
                "source": "System",
                "message": f"Optimization complete! Best molecule: {current_best_smiles} with {property}={current_best_value:.3f}",
                "smiles": current_best_smiles,
            },
        }
    )

    edge_data = Edge(
        id=f"edge_{initial_node_id}_{initial_node_id + 1}",
        fromNode=f"node_{initial_node_id}",
        toNode=f"node_{initial_node_id + 1}",
        status="complete",
        label="Completed",
    )

    await websocket.send_json({"type": "edge_update", "edge": edge_data.json()})
    logger.info(f"Sending initial edge: {edge_data}")

    await websocket.send_json({"type": "complete"})
