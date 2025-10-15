from fastapi import WebSocket
from charge.clients.autogen import AutoGenClient
from backend_helper_funcs import RetroSynthesisContext
from charge.servers.AiZynthTools import RetroPlanner


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
    retro_synth_context: RetroSynthesisContext,
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
                    the experiment in experiment_type, which contains the prompts
            aizynth_planner (RetroPlanner): AiZynthFinder-based planner instance
            retro_synth_context (RetroSynthesisContext): Context/configuration object
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
            "message": f"Finding synthesis pathway for {parent_smiles}...",
            "smiles": parent_smiles,
        }
    )
    await websocket.send_json(
        {
            "type": "response",
            "message": f"Searching for alternatives without {constraint_smiles}",
            "smiles": constraint_smiles,
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
    retro_synth_context: RetroSynthesisContext,
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
                    the experiment in experiment_type, which contains the prompts
            aizynth_planner (RetroPlanner): AiZynthFinder-based planner instance
            retro_synth_context (RetroSynthesisContext): Context/configuration object
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
            "message": f"Finding synthesis pathway to {parent_smiles}...",
            "smiles": parent_smiles,
        }
    )

    result = await charge_runner.run()
    retro_results = result.to_dict()

    reactants = retro_results.get("reactants_smiles_list", [])
    products = retro_results.get("products_smiles_list", [])

    # Please implement unconstrained retrosynthesis here

    await websocket.send_json({"type": "complete"})
