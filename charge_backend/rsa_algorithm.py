"""
Generic Recursive Self-Aggregation (RSA) Algorithm

This module provides the core N-K-T RSA loop that can be used by any task type
(retrosynthesis, LMO, etc.). It uses the existing ChARGe Agent/Task framework
without introducing new orchestration layers.

RSA Algorithm:
    Stage 1: Generate N diverse proposals
    Stages 2-T: Recursively aggregate K-subset proposals
    Output: Single best proposal from final stage
"""

import random
import json
import os
import asyncio
from typing import Any, Callable, Optional
from pathlib import Path


async def run_rsa_loop(
    n: int,
    k: int,
    t: int,
    create_proposal_task: Callable[[], Any],
    create_aggregation_task: Callable[[str, list[dict], int, int], Any],
    format_candidates: Callable[[list[dict]], str],
    runner: Any,
    log_progress: Callable,
    clogger: Any,
    log_dir: str,
    output_schema: Any,
    callback_handler: Optional[Any] = None,
    parallel: bool = True,
    runner_factory: Optional[Callable[[], Any]] = None,
) -> tuple[str, Any]:
    """
    Execute the generic N-K-T RSA algorithm.

    Args:
        n: Number of initial proposals to generate
        k: Size of subsets for aggregation (K <= N)
        t: Total number of stages (including initial proposals)
        create_proposal_task: Factory function that returns a Task for proposals
        create_aggregation_task: Factory function that returns aggregation Task
                                 Takes (candidates_text, subset, step, total_steps)
        format_candidates: Function to format proposals into text for aggregation
                          Takes list of proposals, returns formatted string
        runner: ChARGe agent runner with .task and .run() interface
        log_progress: Callback for reasoning progress
        clogger: Callback logger for UI messages
        log_dir: Directory to save execution logs
        output_schema: Pydantic schema for validating outputs
        callback_handler: Optional callback handler to drain after each task
        parallel: If True, generate initial proposals in parallel (default: True)
        runner_factory: Factory to create independent runner instances for parallel execution
                       If None and parallel=True, uses the provided runner (not truly parallel)

    Returns:
        tuple: (final_output_json, final_result_object)

    Raises:
        ValueError: If all proposals fail or invalid parameters
    """

    # Validate parameters
    if n < 1 or k < 1 or t < 1:
        raise ValueError(f"Invalid RSA parameters: N={n}, K={k}, T={t} (all must be >= 1)")
    if k > n:
        await clogger.warning(f"K ({k}) > N ({n}), adjusting K to N")
        k = n

    # Helper function to run a single proposal
    async def run_single_proposal(proposal_index: int, proposal_runner: Any):
        """Run a single proposal and return result or None if failed"""
        try:
            await clogger.info(f"Generating proposal {proposal_index+1}/{n}")

            # Create proposal task
            proposal_task = create_proposal_task()
            proposal_runner.task = proposal_task

            # Disable validation if requested
            if os.getenv("CHARGE_DISABLE_OUTPUT_VALIDATION", "0") == "1":
                proposal_task.structured_output_schema = None

            # Save proposal prompt
            proposer_log = {
                "proposal_index": proposal_index + 1,
                "system_prompt": proposal_task.get_system_prompt(),
                "user_prompt": proposal_task.get_user_prompt(),
            }
            with open(f"{log_dir}/proposer_{proposal_index+1:02d}_prompt.json", "w") as f:
                json.dump(proposer_log, f, indent=2)

            # Run proposal
            proposal_output = await proposal_runner.run(log_progress)
            if callback_handler:
                await callback_handler.drain()

            # Validate output
            proposal_result = output_schema.model_validate_json(proposal_output)

            # Check if proposal contains actual chemical information
            if not hasattr(proposal_result, 'reactants_smiles_list') or len(proposal_result.reactants_smiles_list) == 0:
                await clogger.warning(f"Proposal {proposal_index+1} has empty reactants (model refused or failed), skipping")
                return None

            # Save proposal output
            proposer_output_log = {
                "proposal_index": proposal_index + 1,
                "result": proposal_result.model_dump(),
                "full_output": json.loads(proposal_output)
            }
            with open(f"{log_dir}/proposer_{proposal_index+1:02d}_output.json", "w") as f:
                json.dump(proposer_output_log, f, indent=2)

            await clogger.info(f"Proposal {proposal_index+1} completed successfully")

            return {
                "output": proposal_output,
                "result": proposal_result,
                "index": proposal_index
            }

        except Exception as e:
            await clogger.warning(f"Proposal {proposal_index+1} failed: {str(e)}")
            return None

    # Stage 1: Generate N initial proposals
    await clogger.info(f"RSA Step 1/{t}: Generating {n} initial proposals" +
                      (" (parallel mode)" if parallel else " (sequential mode)"))

    if parallel and runner_factory:
        # Parallel mode: generate all proposals concurrently
        proposal_tasks = []
        for i in range(n):
            # Create independent runner for each proposal
            proposal_runner = runner_factory()
            task = run_single_proposal(i, proposal_runner)
            proposal_tasks.append(task)

        # Run all proposals in parallel
        proposal_results = await asyncio.gather(*proposal_tasks, return_exceptions=True)

        # Filter valid proposals and handle exceptions
        proposals = []
        for i, result in enumerate(proposal_results):
            if isinstance(result, Exception):
                await clogger.warning(f"Proposal {i+1} failed with exception: {str(result)}")
            elif result is not None:
                proposals.append(result)
    else:
        # Sequential mode: generate proposals one by one
        if parallel and not runner_factory:
            await clogger.warning("Parallel mode requested but no runner_factory provided, falling back to sequential")

        proposals = []
        for i in range(n):
            result = await run_single_proposal(i, runner)
            if result is not None:
                proposals.append(result)

    if not proposals:
        raise ValueError("All RSA proposals failed")

    await clogger.info(f"Generated {len(proposals)} valid proposals")

    # Stages 2-T: Recursive aggregation
    current_proposals = proposals

    for step in range(2, t + 1):
        await clogger.info(
            f"RSA Step {step}/{t}: Aggregating {len(current_proposals)} proposals into {k}-subsets"
        )

        # Adjust K if needed
        current_k = k
        if len(current_proposals) < k:
            await clogger.warning(
                f"Not enough proposals ({len(current_proposals)}) for K={k}, using all available"
            )
            current_k = len(current_proposals)

        # Generate aggregations
        next_proposals = []
        num_aggregations = max(n, len(current_proposals))

        for i in range(num_aggregations):
            await clogger.info(f"Aggregation {i+1}/{num_aggregations}")
            try:
                # Select K random proposals
                if len(current_proposals) <= current_k:
                    subset = current_proposals
                else:
                    subset = random.sample(current_proposals, current_k)

                # Format candidates using task-specific formatter
                candidates_text = format_candidates(subset)
                subset_indices = [prop["index"] + 1 for prop in subset]

                # Create aggregation task
                agg_task = create_aggregation_task(
                    candidates_text,
                    subset,
                    step,
                    t
                )
                runner.task = agg_task

                if os.getenv("CHARGE_DISABLE_OUTPUT_VALIDATION", "0") == "1":
                    agg_task.structured_output_schema = None

                # Save aggregation prompt
                aggregator_log = {
                    "step": step,
                    "aggregation_index": i + 1,
                    "k_subset_indices": subset_indices,
                    "system_prompt": agg_task.get_system_prompt(),
                    "user_prompt": agg_task.get_user_prompt(),
                    "candidates_text": candidates_text,
                }
                with open(f"{log_dir}/aggregator_step{step}_{i+1:02d}_prompt.json", "w") as f:
                    json.dump(aggregator_log, f, indent=2)

                # Run aggregation
                agg_output = await runner.run(log_progress)
                if callback_handler:
                    await callback_handler.drain()

                # Validate output
                agg_result = output_schema.model_validate_json(agg_output)

                # Check if aggregation contains actual chemical information
                if not hasattr(agg_result, 'reactants_smiles_list') or len(agg_result.reactants_smiles_list) == 0:
                    await clogger.warning(f"Aggregation {i+1} has empty reactants (model refused or failed), skipping")
                    continue

                # Save aggregation output
                aggregator_output_log = {
                    "step": step,
                    "aggregation_index": i + 1,
                    "k_subset_indices": subset_indices,
                    "result": agg_result.model_dump(),
                    "full_output": json.loads(agg_output)
                }
                with open(f"{log_dir}/aggregator_step{step}_{i+1:02d}_output.json", "w") as f:
                    json.dump(aggregator_output_log, f, indent=2)

                next_proposals.append({
                    "output": agg_output,
                    "result": agg_result,
                    "index": i,
                    "step": step
                })
                await clogger.info(f"Aggregation {i+1} completed successfully")

            except Exception as e:
                await clogger.warning(f"Aggregation {i+1} failed: {str(e)}")
                continue

        if not next_proposals:
            await clogger.warning(f"No successful aggregations in step {step}, using previous proposals")
            break

        current_proposals = next_proposals

    # Select final proposal (first one from final stage)
    if not current_proposals:
        raise ValueError("RSA failed to produce any valid proposals")

    final_proposal = current_proposals[0]
    final_output = final_proposal["output"]
    final_result = final_proposal["result"]

    # Save final output
    final_log = {
        "final_step": step if step <= t else t,
        "n_proposals": n,
        "k_subset_size": k,
        "t_stages": t,
        "final_result": final_result.model_dump(),
    }
    with open(f"{log_dir}/FINAL_OUTPUT.json", "w") as f:
        json.dump(final_log, f, indent=2)

    await clogger.info(f"RSA completed! Final output saved to {log_dir}")

    return final_output, final_result
