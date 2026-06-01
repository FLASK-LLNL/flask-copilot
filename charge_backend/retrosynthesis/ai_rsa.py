"""RSA-based retrosynthesis driver.

Holds the chemistry-specific orchestration that sits on top of ChARGe's
generic RSA loop: loads the retrosynthesis prompt set (standalone or RAG),
optionally queries the reaction database for RAG context, supplies the
candidate formatter and proposal validator, and dispatches the loop.

The generic N-K-T RSA algorithm itself lives in ``charge.algorithms`` as
``RSATask`` and is domain-agnostic; this module plugs the retrosynthesis
schema, prompts, and validator into it via ``RetroRSATask`` subclass.
"""

from __future__ import annotations

import json
from pathlib import Path
from fastapi import WebSocket

from lc_conductor.callback_logger import CallbackLogger
from lc_conductor import ToolRuntime

from charge.algorithms import RSATask
from charge.experiments.experiment import Experiment
from charge.clients.agent_factory import ReasoningCallbackType

from charge_backend.backend_helper_funcs import CallbackHandler, FlaskRunSettings, Node
from charge_backend.retrosynthesis.retrosynthesis_task import (
    TemplateFreeReactionOutputSchema as ReactionOutputSchema,
)


class RetroRSATask(RSATask):
    """RSATask specialized for template-free retrosynthesis.

    Uses the chemistry output schema by default and overrides the candidate
    formatter / validator to emit and check SMILES-based proposals.
    """

    def __init__(self, **kwargs):
        kwargs.setdefault("structured_output_schema", ReactionOutputSchema)
        super().__init__(**kwargs)

    def format_candidates(self, subset):
        """Format retrosynthesis proposals for the aggregation prompt."""
        text = ""
        for idx, prop in enumerate(subset, 1):
            r = prop["result"]
            text += f"\n---- Candidate {idx} ----\n"
            text += f"Reasoning: {r.reasoning_summary}\n"
            text += f"Reactants: {', '.join(r.reactants_smiles_list)}\n"
            text += f"Products: {', '.join(r.products_smiles_list)}\n"
        return text

    def validate_proposal(self, result):
        """A proposal is valid only if it has at least one reactant."""
        return (
            hasattr(result, "reactants_smiles_list")
            and len(result.reactants_smiles_list) > 0
        )


async def run_rsa_retrosynthesis(
    current_node: Node,
    runner,
    callback_handler,
    run_settings: FlaskRunSettings,
    tool_runtime: ToolRuntime,
    user_prompt: str,
    clogger: CallbackLogger,
    experiment: Experiment,
    websocket: WebSocket,
    node_id: str,
    log_progress: ReasoningCallbackType,
) -> tuple[str, ReactionOutputSchema]:
    """Run RSA-based retrosynthesis. Raises on failure so the caller can
    fall back to standard retrosynthesis (see ai.py)."""

    rsa_n = run_settings.rsa_n if hasattr(run_settings, "rsa_n") else 8
    rsa_k = run_settings.rsa_k if hasattr(run_settings, "rsa_k") else 4
    rsa_t = run_settings.rsa_t if hasattr(run_settings, "rsa_t") else 3
    rsa_mode = (
        run_settings.rsa_mode
        if hasattr(run_settings, "rsa_mode")
        else "standalone"
    )

    await clogger.info(
        f"Running RSA mode: {rsa_mode} with N={rsa_n}, K={rsa_k}, T={rsa_t}"
    )

    # Load chemistry-specific prompts for retrosynthesis (3-part structure)
    prompts_dir = Path(__file__).parent / "prompts"
    system_prompt_file = prompts_dir / f"rsa_{rsa_mode}_system.txt"
    proposal_file = prompts_dir / f"rsa_{rsa_mode}_proposal.txt"
    aggregation_file = prompts_dir / f"rsa_{rsa_mode}_aggregation.txt"

    system_prompt = (
        system_prompt_file.read_text() if system_prompt_file.exists() else None
    )
    proposal_prompt = (
        proposal_file.read_text() if proposal_file.exists() else None
    )
    aggregation_prompt = (
        aggregation_file.read_text() if aggregation_file.exists() else None
    )

    await clogger.info(f"Loading prompts for mode={rsa_mode}")
    await clogger.info(
        f"  System prompt: {system_prompt_file.exists()=}, length={len(system_prompt) if system_prompt else 0}"
    )
    await clogger.info(
        f"  Proposal prompt: {proposal_file.exists()=}, length={len(proposal_prompt) if proposal_prompt else 0}"
    )
    await clogger.info(
        f"  Aggregation prompt: {aggregation_file.exists()=}, length={len(aggregation_prompt) if aggregation_prompt else 0}"
    )

    # For RAG mode: Query database once and inject into prompts
    # For standalone mode: Remove database query tool
    user_prompt_with_rag = user_prompt
    available_tools = tool_runtime.mcp_server_urls
    builtin_tools_filtered = tool_runtime.direct_tools
    db_results_to_save = (
        None  # Store DB results to save later after log_dir is created
    )

    if rsa_mode == "rag":
        await clogger.info("RAG mode: Querying reaction database once...")
        try:
            from charge_backend.retrosynthesis.database import query_reaction_database

            db_results = query_reaction_database(current_node.smiles, top_k=10)

            if db_results and not any("error" in r for r in db_results):
                await clogger.info(
                    f"Found {len(db_results)} similar reactions in database"
                )

                summary_lines = [
                    f"**Database Query Results ({len(db_results)} reactions found):**"
                ]
                for idx, reaction in enumerate(db_results[:5], 1):
                    name = reaction.get("name", f"Reaction {idx}")
                    summary_lines.append(f"  {idx}. {name}")
                    if "components" in reaction and reaction["components"]:
                        reactants = [
                            c.get("name", c.get("smiles", "?"))
                            for c in reaction["components"]
                            if c.get("role") in ["Reactant", "Reagent"]
                        ]
                        products = [
                            c.get("name", c.get("smiles", "?"))
                            for c in reaction["components"]
                            if c.get("role") == "Product"
                        ]
                        if reactants:
                            summary_lines.append(
                                f"     Reactants: {', '.join(reactants[:3])}"
                            )
                        if products:
                            summary_lines.append(
                                f"     Products: {', '.join(products[:2])}"
                            )
                if len(db_results) > 5:
                    summary_lines.append(
                        f"  ... and {len(db_results) - 5} more reactions"
                    )
                summary_lines.append(
                    "These reactions will be injected into all proposal prompts."
                )
                await clogger.info("\n".join(summary_lines))

                rag_context = "\n\n--- REACTION DATABASE RESULTS ---\n"
                rag_context += f"These are similar reactions retrieved by comparing structural similarity to the target product ({current_node.smiles}):\n\n"

                for idx, reaction in enumerate(db_results[:10], 1):
                    rag_context += f"Reaction {idx}:\n"
                    if "reactants" in reaction:
                        rag_context += (
                            f"  Reactants: {reaction.get('reactants', 'N/A')}\n"
                        )
                    if "products" in reaction:
                        rag_context += (
                            f"  Products: {reaction.get('products', 'N/A')}\n"
                        )
                    if "text" in reaction and reaction.get("text"):
                        rag_context += f"  Description: {reaction['text']}\n"
                    rag_context += "\n"

                rag_context += "Use these reactions as supporting evidence for your retrosynthesis proposal.\n"
                rag_context += "--- END DATABASE RESULTS ---\n"

                user_prompt_with_rag = user_prompt + rag_context
                db_results_to_save = db_results
            else:
                await clogger.info("No reactions found in database")
                user_prompt_with_rag = (
                    user_prompt
                    + "\n\nNo similar reactions found in the database for this target molecule.\n"
                )
        except Exception as e:
            await clogger.warning(f"Database query failed: {str(e)}")
            user_prompt_with_rag = (
                user_prompt
                + "\n\nDatabase query failed. Proceed using chemistry knowledge only.\n"
            )

        # Filter out query_reaction_database from builtin tools (already queried once)
        builtin_tools_filtered = [
            tool
            for tool in builtin_tools_filtered
            if getattr(tool, "__name__", "") != "query_reaction_database"
        ]
        await clogger.info(
            "RAG mode: Removed query_reaction_database from tools (already queried)"
        )

    elif rsa_mode == "standalone":
        # Standalone mode: Remove database query tool entirely (no retrieval)
        builtin_tools_filtered = [
            tool
            for tool in builtin_tools_filtered
            if getattr(tool, "__name__", "") != "query_reaction_database"
        ]
        await clogger.info(
            "Standalone mode: Removed query_reaction_database from tools (no retrieval)"
        )

    # Build the retrosynthesis-specialized RSATask. All chemistry-specific
    # customization (schema, formatter, validator) is on RetroRSATask.
    rsa_task = RetroRSATask(
        n=rsa_n,
        k=rsa_k,
        t=rsa_t,
        user_prompt=user_prompt_with_rag,
        system_prompt=system_prompt,
        proposal_prompt=proposal_prompt,
        aggregation_prompt=aggregation_prompt,
        parallel=True,
        log_dir=None,  # auto-generate timestamp-based directory
        server_urls=available_tools,
        builtin_tools=builtin_tools_filtered,
        bearer_token=tool_runtime.bearer_token,
    )

    rsa_log_dir = rsa_task.log_dir
    await clogger.info(f"RSA execution logs will be saved to: {rsa_log_dir}")

    if db_results_to_save:
        try:
            with open(f"{rsa_log_dir}/database_query_results.json", "w") as f:
                json.dump(db_results_to_save, f, indent=2)
            await clogger.info(
                f"Database query results saved to {rsa_log_dir}/database_query_results.json"
            )
        except Exception as e:
            await clogger.warning(f"Failed to save database results: {str(e)}")

    # Runner factory for parallel execution: independent runner + callback per proposal.
    proposal_counter = [0]

    def create_runner():
        proposal_counter[0] += 1
        independent_callback = CallbackHandler(
            websocket=websocket, name=f"proposal_{proposal_counter[0]}"
        )
        return experiment.create_agent_with_experiment_state(
            task=None,
            agent_name=f"retrosynth_{node_id}_proposal_{proposal_counter[0]}",
            callback=independent_callback,
        )

    output, final_result = await rsa_task.run(
        runner,
        runner_factory=create_runner,
        log_progress=log_progress,
        logger_info=clogger.info,
        logger_warning=clogger.warning,
        logger_error=clogger.error,
        callback_handler=(
            callback_handler
            if isinstance(callback_handler, CallbackHandler)
            else None
        ),
    )

    await clogger.info(
        f"RSA mode completed successfully. Logs saved to: {rsa_log_dir}"
    )

    return output, final_result
