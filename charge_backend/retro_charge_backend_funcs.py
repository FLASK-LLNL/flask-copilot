import asyncio
from concurrent.futures import ProcessPoolExecutor
import os
from fastapi import WebSocket
from charge.clients.autogen import AutoGenAgent
from charge.servers.AiZynthTools import RetroPlanner, ReactionPath
from aizynth_backend_funcs import generate_tree_structure
from loguru import logger
from callback_logger import CallbackLogger
from typing import Literal, Optional, Union

from backend_helper_funcs import (
    CallbackHandler,
    RetrosynthesisContext,
    calculate_positions,
    highlight_node,
    RetrosynthesisContext,
    Node,
    Edge,
    Reaction,
    loop_executor,
)
from molecule_naming import smiles_to_html

from charge.tasks.RetrosynthesisTask import (
    TemplateFreeRetrosynthesisTask as RetrosynthesisTask,
    TemplateFreeReactionOutputSchema as ReactionOutputSchema,
)

from charge.experiments.AutoGenExperiment import AutoGenExperiment

RETROSYNTH_UNCONSTRAINED_USER_PROMPT_TEMPLATE = (
    "Provide a retrosynthetic pathway for the target molecule {target_molecule}. "
    + "The pathway should be provided as a tuple of reactants as SMILES and the product as SMILES. "
    + "Perform only single step retrosynthesis. Make sure the SMILES strings are valid. "
    + "Use tools to verify the SMILES strings and diagnose any issues that arise."
    + "Do the evaluation step-by-step. Propose a retrosynthetic step, then evaluate it. "
    + "If the evaluation fails, propose a new retrosynthetic step and evaluate it again. "
    + "Find the best possible retrosynthetic step, and use tools to see if the "
    + "proposed reactants are synthesizable. "
)

RETROSYNTH_CONSTRAINED_USER_PROMPT_TEMPLATE = (
    "Provide a retrosynthetic pathway for the target molecule {target_molecule}. "
    + "The pathway should be provided as a tuple of reactants as SMILES and the product as SMILES. "
    + "Perform only single step retrosynthesis. Make sure the SMILES strings are valid. "
    + "Use tools to verify the SMILES strings and diagnose any issues that arise. "
    + "The following reactant cannot be used in the retrosynthetic step: {constrained_reactant}. "
    + "Do the evaluation step-by-step. Propose a retrosynthetic step, then evaluate it. "
    + "If the evaluation fails, propose a new retrosynthetic step and evaluate it again. "
)


def run_retro_planner(config_file, smiles):
    logger.info(f"Running RetroPlanner for SMILES: {smiles}")
    planner = RetroPlanner(configfile=config_file)

    _, _, routes = planner.plan(smiles)
    return routes, planner


async def generate_molecules(
    start_smiles: str,
    config_file: str,
    context: RetrosynthesisContext,
    executor: ProcessPoolExecutor,
    websocket: WebSocket,
    available_tools: Optional[Union[str, list[str]]] = None,
    molecule_name_format: Literal["brand", "iupac", "formula", "smiles"] = "brand",
):
    """Stream positioned nodes and edges"""
    clogger = CallbackLogger(websocket, source="generate_molecules")
    await clogger.info(
        f"Planning retrosynthesis for: {start_smiles} with available tools: {available_tools}."
    )

    # Generate and position entire tree upfront

    root = Node(
        id="node_0",
        smiles=start_smiles,
        label=smiles_to_html(start_smiles, molecule_name_format),
        hoverInfo=f"# Root molecule \n **SMILES:** {start_smiles}",
        level=0,
        parentId=None,
        cost=None,
        bandgap=None,
        yield_=None,
        highlight="yellow",
        x=100,
        y=100,
    )
    await websocket.send_json({"type": "node", "node": root.json()})

    await clogger.info("Starting planning in executor...")

    # Disable executor for now due to bad interaction with task_done callbacks
    # routes, planner = await loop_executor(
    #     executor, run_retro_planner, config_file, start_smiles
    # )

    routes, planner = run_retro_planner(config_file, start_smiles)
    await clogger.info(f"Running RetroPlanner for SMILES: {start_smiles}")

    context.node_id_to_planner[root.id] = planner

    if not routes:
        await clogger.info(
            f"No synthesis routes found for {start_smiles}.",
            smiles=start_smiles,
        )
        await websocket.send_json({"type": "complete"})
        return
    await clogger.info(f"Found {len(routes)} routes for {start_smiles}.")

    planner.last_route_used = 0
    reaction_path = ReactionPath(route=routes[0])
    nodes, edges = generate_tree_structure(reaction_path.nodes, context)
    logger.info(f"Generated {len(nodes)} nodes and {len(edges)} edges.")

    calculate_positions(nodes)

    await websocket.send_json(
        {
            "type": "node_update",
            "node": {
                "id": root.id,
                "highlight": "normal",
                "reaction": Reaction(
                    "azf",
                    "Reaction found with AiZynthFinder",
                    highlight="yellow",
                    label="Template",
                ).json(),
            },
        }
    )

    # Stream remaining nodes with edges
    for node in nodes[1:]:
        # Find edge for this node
        edge = next((e for e in edges if e.toNode == node.id), None)

        # Send node
        await websocket.send_json({"type": "node", "node": node.json()})
        if edge:
            context.parents[node.id] = edge.fromNode
            await websocket.send_json({"type": "edge", "edge": edge.json()})

    await websocket.send_json({"type": "complete"})


async def optimize_molecule_retro(
    node_id: str,
    context: RetrosynthesisContext,
    query: Optional[str],
    constraint: Optional[str],
    websocket: WebSocket,
    experiment: AutoGenExperiment,
    config_file: str,
    available_tools: Optional[Union[str, list[str]]] = None,
    molecule_name_format: Literal["brand", "iupac", "formula", "smiles"] = "brand",
):
    """Optimize a molecule using retrosynthesis by node ID"""
    current_node = context.node_ids.get(node_id)
    cur_node_id = context.azf_nodes.get(node_id)
    assert current_node is not None, f"Node ID {node_id} not found"
    assert cur_node_id is not None, f"AZF Node ID {node_id} not found"

    cur_node_id = cur_node_id.node_id

    clogger = CallbackLogger(websocket, source="optimize_molecule_retro")

    await clogger.info(
        f"Finding synthesis pathway to {current_node.smiles}... with available tools: {available_tools}.",
        smiles=current_node.smiles,
    )

    if node_id in context.node_id_to_charge_client:
        # Existing context
        runner = context.node_id_to_charge_client[node_id]
    else:
        # New context
        agent_name = experiment.agent_pool.create_agent_name(
            prefix=f"retrosynth_{node_id}_"
        )
        runner = experiment.create_agent_with_experiment_state(
            task=None,
            agent_name=agent_name,
            callback=CallbackHandler(websocket),
        )
        context.node_id_to_charge_client[node_id] = runner

    if constraint:
        user_prompt = RETROSYNTH_CONSTRAINED_USER_PROMPT_TEMPLATE.format(
            target_molecule=current_node.smiles,
            constrained_reactant=constraint,
        )
    else:
        user_prompt = RETROSYNTH_UNCONSTRAINED_USER_PROMPT_TEMPLATE.format(
            target_molecule=current_node.smiles
        )

    user_prompt += "\nDouble check the reactants with the `predict_reaction_products` tool to see if the products are equivalent to the given product. If there is any inconsistency (canonicalize both sides of the equation first), log it and try some other set of reactants."
    if query is not None:
        user_prompt += (
            f"\n\nAdditionally, adhere to the following requirements: {query}"
        )

    retro_task = RetrosynthesisTask(
        user_prompt=user_prompt, server_urls=available_tools
    )
    runner.task = retro_task

    if os.getenv("CHARGE_DISABLE_OUTPUT_VALIDATION", "0") == "1":
        retro_task.structured_output_schema = None
        await clogger.warning(
            "Structure validation disabled for RetrosynthesisTask output schema."
        )

    await clogger.info(
        f"Finding synthesis routes for {current_node.smiles} using available tools: {available_tools}."
    )

    # Run task
    await highlight_node(current_node, websocket, True)
    output = await runner.run()

    if os.getenv("CHARGE_DISABLE_OUTPUT_VALIDATION", "0") == "1":
        await clogger.warning(
            "Structure validation disabled for RetrosynthesisTask output schema."
            "Returning text results without validation first before post-processing."
        )
    await clogger.info(f"Results: {output}")

    result = ReactionOutputSchema.model_validate_json(output)

    await highlight_node(current_node, websocket, False)

    level = current_node.level + 1
    num_nodes = len(context.node_ids)

    reasoning_summary = result.reasoning_summary
    await clogger.info(
        f"Retrosynthesis reasoning summary for {current_node.smiles}:\n{reasoning_summary}",
    )

    # Update node with discovered reaction
    context.node_id_to_reasoning_summary[node_id] = reasoning_summary
    await websocket.send_json(
        {
            "type": "node_update",
            "node": {
                "id": node_id,
                "reaction": Reaction(
                    "ai_reaction_0",
                    reasoning_summary,
                    highlight="red",
                    label="AI",
                ).json(),
            },
        }
    )

    nodes: list[Node] = []
    edges: list[Edge] = []
    for i, smiles in enumerate(result.reactants_smiles_list):
        node = Node(
            f"node_{num_nodes+i}",
            smiles,
            smiles_to_html(smiles, molecule_name_format),
            "Discovered",
            level,
            current_node.id,
        )
        context.node_ids[node.id] = node
        context.parents[node.id] = node_id

        # Find paths for the leaf nodes
        routes, planner = run_retro_planner(config_file, smiles)

        if not routes:
            await clogger.warning(f"No routes found for {smiles}. Skipping...")
            continue

        await clogger.info(f"Found {len(routes)} routes for {smiles}.")

        context.node_id_to_planner[node.id] = planner
        planner.last_route_used = 0

        reaction_path = ReactionPath(
            route=routes[0],
            cur_num_nodes=num_nodes,
            root_parent_node_id=cur_node_id,
        )
        child_nodes, child_edges = generate_tree_structure(
            reaction_path.nodes, context, start_level=level
        )

        # Stream child nodes and edges
        nodes.extend(child_nodes)
        edges.extend(child_edges)

        num_nodes = reaction_path.num_nodes

    calculate_positions(nodes, context.nodes_per_level[level])

    for node, edge in zip(nodes, edges):
        await websocket.send_json({"type": "node", "node": node.json()})
        await websocket.send_json({"type": "edge", "edge": edge.json()})
    await websocket.send_json({"type": "complete"})
