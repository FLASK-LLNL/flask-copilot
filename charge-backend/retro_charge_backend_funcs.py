from fastapi import WebSocket
from charge.clients.autogen import AutoGenClient
from charge.servers.AiZynthTools import RetroPlanner, ReactionPath
from aizynth_backend_funcs import generate_tree_structure
from loguru import logger
from typing import Optional, Union, Dict, cast

from backend_helper_funcs import (
    CallbackHandler,
    RetrosynthesisContext,
    calculate_positions,
    highlight_node,
    RetrosynthesisContext,
    Node,
    Edge,
)

from charge.tasks.RetrosynthesisTask import (
    TemplateFreeRetrosynthesisTask as RetrosynthesisTask,
    TemplateFreeReactionOutputSchema as ReactionOutputSchema,
)

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


def get_constrained_prompt(target_molecule: str, constrained_reactant: str) -> str:
    """Generate a constrained retrosynthesis prompt."""
    return RETROSYNTH_CONSTRAINED_USER_PROMPT_TEMPLATE.format(
        target_molecule=target_molecule,
        constrained_reactant=constrained_reactant,
    )


def get_unconstrained_prompt(target_molecule: str) -> str:
    """Generate an unconstrained retrosynthesis prompt."""
    return RETROSYNTH_UNCONSTRAINED_USER_PROMPT_TEMPLATE.format(
        target_molecule=target_molecule,
    )


async def constrained_retro(
    parent_id: str,
    parent_smiles: str,
    constraint_id: str,
    constraint_smiles: str,
    charge_runner: AutoGenClient,
    aizynth_planner: RetroPlanner,
    retro_synth_context: RetrosynthesisContext,
    websocket: WebSocket,
):
    """
    Perform a constrained retrosynthetic optimization for a target molecule and
    stream progress and results over a WebSocket.

    This coroutine coordinates an automated retrosynthesis "charge" run and then
    applies a user-specified constraint to prefer or filter routes that avoid a
    given substructure/compound. It sends progress updates and final completion
    status over the provided WebSocket connection. The progress updates through
    websocket are sent by the charge_runner only and the runner must be set up
    with the CallBackHandler to send progress updates.



    Args:
            parent_id (str): Identifier for the target (parent) molecule in the system.
            parent_smiles (str): SMILES string of the target molecule to be synthesized.
            constraint_id (str): Identifier for the constrained substructure/fragment.
            constraint_smiles (str): SMILES string for the substructure/fragment to avoid
                    in proposed synthetic routes.
            charge_runner (AutoGenClient): An asynchronous runner client that, when
                    awaited (charge_runner.run()), performs the retrosynthesis / optimization
                    and returns a result object with a to_dict() method. This object also wraps
                    the task in task, which contains the prompts
            aizynth_planner (RetroPlanner): AiZynthFinder-based planner instance
            retro_synth_context (RetrosynthesisContext): Context/configuration object
                    containing Node objects representing the retrosynthetic state.
            websocket (WebSocket): Open WebSocket connection used to stream JSON messages
                    back to the caller. Messages sent by this function include at minimum
                    initial progress updates and a final completion message.

    Returns:
            None

    Side effects:
            - Sends JSON messages via websocket.send_json(...) describing progress,
                intermediate results, and completion.
            - Awaits charge_runner.run(), which may perform network or CPU-bound work.
            - Potentially invokes planner methods that mutate or read retro_synth_context.
            - Update the UI via websocket messages.
            - Update the retro_synth_context with new nodes and edges.

    Errors:
            - Exceptions raised by charge_runner.run(), aizynth_planner methods, or
                websocket.send_json() are propagated to the caller. Callers should handle
                network errors, timeouts, and planner failures.
    """

    await websocket.send_json(
        {
            "type": "response",
            "message": {
                "source": "System",
                "message": f"Finding synthesis pathway for {parent_smiles}...",
                "smiles": parent_smiles,
            }
        }
    )
    await websocket.send_json(
        {
            "type": "response",
            "message": {
                "source": "System",
                "message": f"Searching for alternatives without {constraint_smiles}",
                "smiles": constraint_smiles,
            }
        }
    )

    result = await charge_runner.run()
    retro_results = result.to_dict()

    reactants = retro_results.get("reactants_smiles_list", [])
    products = retro_results.get("products_smiles_list", [])
    # Please implement constrained retrosynthesis here

    await websocket.send_json({"type": "complete"})


async def unconstrained_retro(
    parent_id: str,
    parent_smiles: str,
    charge_runner: AutoGenClient,
    aizynth_planner: RetroPlanner,
    retro_synth_context: RetrosynthesisContext,
    websocket: WebSocket,
):
    """
    Perform a constrained retrosynthetic optimization for a target molecule and
    stream progress and results over a WebSocket.

    This coroutine coordinates an automated retrosynthesis "charge" run and then
    applies a user-specified constraint to prefer or filter routes that avoid a
    given substructure/compound. It sends progress updates and final completion
    status over the provided WebSocket connection. The progress updates through
    websocket are sent by the charge_runner only and the runner must be set up
    with the CallBackHandler to send progress updates.



    Args:
            parent_id (str): Identifier for the target (parent) molecule in the system.
            parent_smiles (str): SMILES string of the target molecule to be synthesized.
            constraint_id (str): Identifier for the constrained substructure/fragment.
            constraint_smiles (str): SMILES string for the substructure/fragment to avoid
                    in proposed synthetic routes.
            charge_runner (AutoGenClient): An asynchronous runner client that, when
                    awaited (charge_runner.run()), performs the retrosynthesis / optimization
                    and returns a result object with a to_dict() method. This object also wraps
                    the task in task, which contains the prompts
            aizynth_planner (RetroPlanner): AiZynthFinder-based planner instance
            retro_synth_context (RetrosynthesisContext): Context/configuration object
                    containing Node objects representing the retrosynthetic state.
            websocket (WebSocket): Open WebSocket connection used to stream JSON messages
                    back to the caller. Messages sent by this function include at minimum
                    initial progress updates and a final completion message.

    Returns:
            None

    Side effects:
            - Sends JSON messages via websocket.send_json(...) describing progress,
                intermediate results, and completion.
            - Awaits charge_runner.run(), which may perform network or CPU-bound work.
            - Potentially invokes planner methods that mutate or read retro_synth_context.
            - Update the UI via websocket messages.
            - Update the retro_synth_context with new nodes and edges.

    Errors:
            - Exceptions raised by charge_runner.run(), aizynth_planner methods, or
                websocket.send_json() are propagated to the caller. Callers should handle
                network errors, timeouts, and planner failures.
    """

    await websocket.send_json(
        {
            "type": "response",
            "message": {
                "source": "System",
                "message": f"Finding synthesis pathway to {parent_smiles}...",
                "smiles": parent_smiles,
            }
        }
    )

    result = await charge_runner.run()
    retro_results = result.to_dict()

    reactants = retro_results.get("reactants_smiles_list", [])
    products = retro_results.get("products_smiles_list", [])

    # Please implement unconstrained retrosynthesis here

    await websocket.send_json({"type": "complete"})


async def generate_molecules(start_smiles: str, config_file: str, context: RetrosynthesisContext, websocket: WebSocket):
    """Stream positioned nodes and edges"""
    logger.info(f"Planning retrosynthesis for: {start_smiles}")

    # Generate and position entire tree upfront

    root = Node(
        id="node_0",
        smiles=start_smiles,
        label=start_smiles,
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
    planner = context.node_id_to_planner[root.id] = RetroPlanner(configfile=config_file)
    _, _, routes = planner.plan(start_smiles)
    if not routes:
        await websocket.send_json(
            {
                "type": "response",
                "message": {
                    "source": "System",
                    "message": f"No synthesis routes found for {start_smiles}.",
                    "smiles": start_smiles,
                }
            }
        )
        await websocket.send_json({"type": "complete"})
        return
    logger.info(f"Found {len(routes)} routes for {start_smiles}.")

    planner.last_route_used = 0
    reaction_path = ReactionPath(route=routes[0])
    nodes, edges = generate_tree_structure(reaction_path.nodes, context)
    logger.info(f"Generated {len(nodes)} nodes and {len(edges)} edges.")

    calculate_positions(nodes)

    # Stream remaining nodes with edges
    for node in nodes[1:]:
        # Find edge for this node
        edge = next((e for e in edges if e.toNode == node.id), None)

        # Send node
        await websocket.send_json({"type": "node", "node": node.json()})
        if edge:
            await websocket.send_json({"type": "edge", "edge": edge.json()})

    await websocket.send_json({"type": "node_update", "node": {"id": root.id, "highlight": "normal"}})
    await websocket.send_json({"type": "complete"})




async def constrained_opt(parent_smiles, constraint_smiles, planner, websocket: WebSocket):
    """Constrained optimization using retrosynthesis"""

    await websocket.send_json(
        {
            "type": "response",
            "message": {
                "source": "System",
                "message": f"Finding synthesis pathway for {parent_smiles}...",
                "smiles": parent_smiles,
            }
        }
    )
    await websocket.send_json(
        {
            "type": "response",
            "message": {
                "source": "System",
                "message": f"Searching for alternatives without {constraint_smiles}",
                "smiles": constraint_smiles,
            }
        }
    )

    user_prompt = RETROSYNTH_CONSTRAINED_USER_PROMPT_TEMPLATE.format(
        target_molecule=parent_smiles, constrained_reactant=constraint_smiles
    )
    retro_task = RetrosynthesisTask(user_prompt=user_prompt)
    planner.task = retro_task
    logger.info(f"Optimizing {parent_smiles} without using {constraint_smiles} in the synthesis.")
    result = await planner.run()
    return result


async def optimize_molecule_retro(node_id: str,
                                  context: RetrosynthesisContext,
                                  websocket: WebSocket,
                                  model: str,
                                  backend: str,
                                  api_key: str,
                                  model_kwargs: Optional[dict] = None,
                                  server_urls: Optional[Union[str, list[str]]] = None):
    """Optimize a molecule using retrosynthesis by node ID"""
    current_node = context.node_ids.get(node_id)
    assert current_node is not None, f"Node ID {node_id} not found"

    await websocket.send_json(
        {
            "type": "response",
            "message": {
                "source": "System",
                "message": f"Finding synthesis pathway to {current_node.smiles}...",
                "smiles": current_node.smiles,
            }
        }
    )

    if node_id in context.node_id_to_charge_client:
        # Existing context
        runner = context.node_id_to_charge_client[node_id]
    else:
        # New context
        user_prompt = RETROSYNTH_UNCONSTRAINED_USER_PROMPT_TEMPLATE.format(target_molecule=current_node.smiles)
        user_prompt += "\nDouble check the reactants with the `predict_reaction_products` tool to see if the products are equivalent to the given product. If there is any inconsistency (canonicalize both sides of the equation first), log it and try some other set of reactants."
        retro_task = RetrosynthesisTask(user_prompt=user_prompt)
        runner = AutoGenClient(
            task=retro_task,
            model=model,
            backend=backend,
            api_key=api_key,
            model_kwargs=model_kwargs,
            server_url=server_urls,
            thoughts_callback=CallbackHandler(websocket),
        )
        context.node_id_to_charge_client[node_id] = runner

    logger.info(f"Optimizing {current_node.smiles} using retrosynthesis.")

    # Run task
    await highlight_node(current_node, websocket, True)
    result: ReactionOutputSchema = cast(ReactionOutputSchema, await runner.run())
    await highlight_node(current_node, websocket, False)

    level = current_node.level + 1
    num_nodes = len(context.node_ids)

    nodes: list[Node] = []
    edges: list[Edge] = []
    for i, smiles in enumerate(result.reactants_smiles_list):
        node = Node(f"node_{num_nodes+i}", smiles, smiles, "Discovered", level, current_node.id)
        nodes.append(node)
        context.node_ids[node.id] = node
        edges.append(Edge(f"edge_{node_id}_{node.id}", node_id, node.id, "complete"))

    calculate_positions(nodes, context.nodes_per_level[level])
    context.nodes_per_level[level] += len(nodes)

    for node, edge in zip(nodes, edges):
        await websocket.send_json({"type": "node", "node": node.json()})
        await websocket.send_json({"type": "edge", "edge": edge.json()})
