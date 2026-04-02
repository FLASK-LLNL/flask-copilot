import os
import asyncio
import random
import json
import datetime
from pathlib import Path
from fastapi import WebSocket
from lc_conductor.callback_logger import CallbackLogger
from typing import Any, Callable, Optional, Union

from backend_helper_funcs import (
    CallbackHandler,
    highlight_node,
    Node,
    Reaction,
    ReactionAlternative,
    PathwayStep,
    FlaskRunSettings,
)
from charge_backend.prompt_debugger import debug_prompt
from retrosynthesis.context import RetrosynthesisContext
from charge_backend.moleculedb.molecule_naming import (
    smiles_to_html,
)
from retrosynthesis.template import (
    generate_nodes_for_molecular_graph,
    run_retro_planner,
)
from charge_backend.moleculedb.purchasable import is_purchasable
from retrosynthesis.mapping import build_mapped_reaction_dict_or_none

from retrosynthesis.retrosynthesis_task import (
    TemplateFreeRetrosynthesisTask as RetrosynthesisTask,
    TemplateFreeReactionOutputSchema as ReactionOutputSchema,
    RSAAggregationTask,
)

from charge.experiments.experiment import Experiment
from charge.clients.agent_factory import ReasoningCallbackType


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


async def ai_based_retrosynthesis(
    node_id: str,
    context: RetrosynthesisContext,
    query: Optional[str],
    constraint: Optional[str],
    websocket: WebSocket,
    experiment: Experiment,
    config_file: str,
    run_settings: FlaskRunSettings,
    available_tools: Optional[Union[str, list[str]]],
    builtin_tools: Optional[list[Callable[..., Any]]],
    log_progress: ReasoningCallbackType,
):
    """Performs template-free retrosynthesis using the AI orchestrator."""
    clogger = CallbackLogger(websocket, source="ai_based_retrosynthesis")
    current_node = context.node_ids.get(node_id)
    if current_node is None:
        await clogger.error(f"Node ID {node_id} not found")
        await websocket.send_json({"type": "complete"})
        return

    await clogger.info(
        f"Finding synthesis pathway to {current_node.smiles}... with available tools: {available_tools}.",
        smiles=current_node.smiles,
    )

    if node_id in context.node_id_to_charge_client:
        # Existing context
        runner = context.node_id_to_charge_client[node_id]
    else:
        # New context
        callback_handler = CallbackHandler(websocket)
        runner = experiment.create_agent_with_experiment_state(
            task=None,
            agent_name=f"retrosynth_{node_id}",
            callback=callback_handler,
        )
        context.node_id_to_charge_client[node_id] = runner

    callback_handler = getattr(runner, "callback", None)

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
        user_prompt=user_prompt,
        server_urls=available_tools,
        builtin_tools=builtin_tools or [],
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

    # Run task (with optional RSA mode)
    await highlight_node(current_node, websocket, True)
    if run_settings.prompt_debugging:
        await debug_prompt(runner, websocket)

    # Check if RSA mode is enabled
    if run_settings.use_rsa:
        try:
            # RSA Mode: Recursive Self-Aggregation
            rsa_n = run_settings.rsa_n if hasattr(run_settings, 'rsa_n') else 8
            rsa_k = run_settings.rsa_k if hasattr(run_settings, 'rsa_k') else 4
            rsa_t = run_settings.rsa_t if hasattr(run_settings, 'rsa_t') else 3
            rsa_mode = run_settings.rsa_mode if hasattr(run_settings, 'rsa_mode') else "standalone"

            await clogger.info(
                f"Running RSA mode: {rsa_mode} with N={rsa_n}, K={rsa_k}, T={rsa_t}"
            )

            # Create directory for RSA execution logs
            import datetime
            timestamp = datetime.datetime.now().strftime("%Y%m%d_%H%M%S")
            rsa_log_dir = f"/tmp/rsa_execution_{timestamp}"
            os.makedirs(rsa_log_dir, exist_ok=True)
            await clogger.info(f"RSA execution logs will be saved to: {rsa_log_dir}")

            # For RAG mode: Query database once and inject into prompts
            # For standalone mode: Remove database query tool
            user_prompt_with_rag = user_prompt
            builtin_tools_filtered = builtin_tools or []

            if rsa_mode == "rag":
                await clogger.info("RAG mode: Querying reaction database once...")
                try:
                    from retrosynthesis.database import query_reaction_database
                    db_results = query_reaction_database(current_node.smiles, top_k=10)

                    if db_results and not any("error" in r for r in db_results):
                        await clogger.info(f"Found {len(db_results)} similar reactions in database")

                        # Log summary of database results to UI
                        summary_lines = [f"**Database Query Results ({len(db_results)} reactions found):**"]
                        for idx, reaction in enumerate(db_results[:5], 1):  # Show first 5 in UI
                            name = reaction.get('name', f'Reaction {idx}')
                            summary_lines.append(f"  {idx}. {name}")
                            if 'components' in reaction and reaction['components']:
                                # Extract reactants and products
                                reactants = [c.get('name', c.get('smiles', '?')) for c in reaction['components']
                                           if c.get('role') in ['Reactant', 'Reagent']]
                                products = [c.get('name', c.get('smiles', '?')) for c in reaction['components']
                                          if c.get('role') == 'Product']
                                if reactants:
                                    summary_lines.append(f"     Reactants: {', '.join(reactants[:3])}")
                                if products:
                                    summary_lines.append(f"     Products: {', '.join(products[:2])}")
                        if len(db_results) > 5:
                            summary_lines.append(f"  ... and {len(db_results) - 5} more reactions")
                        summary_lines.append("These reactions will be injected into all proposal prompts.")
                        await clogger.info("\n".join(summary_lines))

                        # Format database results with clear context
                        rag_context = "\n\n--- REACTION DATABASE RESULTS ---\n"
                        rag_context += f"These are similar reactions retrieved by comparing structural similarity to the target product ({current_node.smiles}):\n\n"

                        for idx, reaction in enumerate(db_results[:10], 1):
                            rag_context += f"Reaction {idx}:\n"
                            if "reactants" in reaction:
                                rag_context += f"  Reactants: {reaction.get('reactants', 'N/A')}\n"
                            if "products" in reaction:
                                rag_context += f"  Products: {reaction.get('products', 'N/A')}\n"
                            if "text" in reaction and reaction.get("text"):
                                rag_context += f"  Description: {reaction['text']}\n"
                            rag_context += "\n"

                        rag_context += "Use these reactions as supporting evidence for your retrosynthesis proposal.\n"
                        rag_context += "--- END DATABASE RESULTS ---\n"

                        user_prompt_with_rag = user_prompt + rag_context

                        # Save database results to log
                        with open(f"{rsa_log_dir}/database_query_results.json", "w") as f:
                            json.dump(db_results, f, indent=2)
                    else:
                        await clogger.info("No reactions found in database")
                        user_prompt_with_rag = user_prompt + "\n\nNo similar reactions found in the database for this target molecule.\n"
                except Exception as e:
                    await clogger.warning(f"Database query failed: {str(e)}")
                    user_prompt_with_rag = user_prompt + "\n\nDatabase query failed. Proceed using chemistry knowledge only.\n"

                # Filter out query_reaction_database from builtin tools (already queried once)
                builtin_tools_filtered = [
                    tool for tool in builtin_tools_filtered
                    if getattr(tool, '__name__', '') != 'query_reaction_database'
                ]
                await clogger.info("RAG mode: Removed query_reaction_database from tools (already queried)")

            elif rsa_mode == "standalone":
                # Standalone mode: Remove database query tool entirely (no retrieval)
                builtin_tools_filtered = [
                    tool for tool in builtin_tools_filtered
                    if getattr(tool, '__name__', '') != 'query_reaction_database'
                ]
                await clogger.info("Standalone mode: Removed query_reaction_database from tools (no retrieval)")

            # Step 1: Generate N initial proposals
            await clogger.info(f"RSA Step 1/{rsa_t}: Generating {rsa_n} initial proposals")
            proposals = []
            for i in range(rsa_n):
                await clogger.info(f"Generating proposal {i+1}/{rsa_n}")
                try:
                    # Create a fresh task for each proposal
                    # TODO: Add temperature=0.8 for models that support it (Claude, etc.)
                    proposal_task = RetrosynthesisTask(
                        user_prompt=user_prompt_with_rag,
                        server_urls=available_tools,
                        builtin_tools=builtin_tools_filtered,
                    )
                    runner.task = proposal_task

                    if os.getenv("CHARGE_DISABLE_OUTPUT_VALIDATION", "0") == "1":
                        proposal_task.structured_output_schema = None

                    # Save proposer prompt for examination
                    proposer_log = {
                        "proposal_index": i + 1,
                        "system_prompt": proposal_task.get_system_prompt(),
                        "user_prompt": user_prompt_with_rag,
                        "mode": rsa_mode,
                    }
                    with open(f"{rsa_log_dir}/proposer_{i+1:02d}_prompt.json", "w") as f:
                        json.dump(proposer_log, f, indent=2)

                    # Run proposal
                    proposal_output = await runner.run(log_progress)
                    if isinstance(callback_handler, CallbackHandler):
                        await callback_handler.drain()

                    # Validate and store
                    proposal_result = ReactionOutputSchema.model_validate_json(proposal_output)

                    # Save proposer output for examination
                    proposer_output_log = {
                        "proposal_index": i + 1,
                        "reasoning_summary": proposal_result.reasoning_summary,
                        "reactants_smiles": proposal_result.reactants_smiles_list,
                        "products_smiles": proposal_result.products_smiles_list,
                        "full_output": json.loads(proposal_output)
                    }
                    with open(f"{rsa_log_dir}/proposer_{i+1:02d}_output.json", "w") as f:
                        json.dump(proposer_output_log, f, indent=2)

                    proposals.append({
                        "output": proposal_output,
                        "result": proposal_result,
                        "index": i
                    })
                    await clogger.info(f"Proposal {i+1} completed successfully")
                except Exception as e:
                    await clogger.warning(f"Proposal {i+1} failed: {str(e)}")
                    continue

            if not proposals:
                raise ValueError("All RSA proposals failed, falling back to standard mode")

            await clogger.info(f"Generated {len(proposals)} valid proposals")

            # Steps 2..T: Recursive aggregation
            current_proposals = proposals
            for step in range(2, rsa_t + 1):
                await clogger.info(
                    f"RSA Step {step}/{rsa_t}: Aggregating {len(current_proposals)} proposals into {rsa_k}-subsets"
                )

                if len(current_proposals) < rsa_k:
                    await clogger.warning(
                        f"Not enough proposals ({len(current_proposals)}) for K={rsa_k}, using all available"
                    )
                    rsa_k = len(current_proposals)

                # Generate new proposals by aggregating K-subsets
                next_proposals = []
                num_aggregations = max(rsa_n, len(current_proposals))

                for i in range(num_aggregations):
                    await clogger.info(f"Aggregation {i+1}/{num_aggregations}")
                    try:
                        # Select K random proposals
                        if len(current_proposals) <= rsa_k:
                            subset = current_proposals
                        else:
                            subset = random.sample(current_proposals, rsa_k)

                        # Format candidates text
                        candidates_text = ""
                        subset_indices = []
                        for idx, prop in enumerate(subset, 1):
                            prop_result = prop["result"]
                            subset_indices.append(prop["index"] + 1)  # Convert to 1-indexed
                            candidates_text += f"\n---- Candidate {idx} ----\n"
                            candidates_text += f"Reasoning: {prop_result.reasoning_summary}\n"
                            candidates_text += f"Reactants: {', '.join(prop_result.reactants_smiles_list)}\n"
                            candidates_text += f"Products: {', '.join(prop_result.products_smiles_list)}\n"

                        # Create aggregation task
                        agg_task = RSAAggregationTask(
                            original_user_prompt=user_prompt,
                            candidates_text=candidates_text,
                            step=step,
                            total_steps=rsa_t,
                            mode=rsa_mode,
                            server_urls=available_tools,
                            builtin_tools=builtin_tools_filtered,
                        )
                        runner.task = agg_task

                        if os.getenv("CHARGE_DISABLE_OUTPUT_VALIDATION", "0") == "1":
                            agg_task.structured_output_schema = None

                        # Save aggregator prompt for examination
                        aggregator_log = {
                            "step": step,
                            "aggregation_index": i + 1,
                            "k_subset_indices": subset_indices,  # Which proposals were selected
                            "system_prompt": agg_task.get_system_prompt(),
                            "user_prompt": agg_task.get_user_prompt(),
                            "original_user_prompt": user_prompt,
                            "candidates_text": candidates_text,
                            "mode": rsa_mode,
                        }
                        with open(f"{rsa_log_dir}/aggregator_step{step}_{i+1:02d}_prompt.json", "w") as f:
                            json.dump(aggregator_log, f, indent=2)

                        # Run aggregation
                        agg_output = await runner.run(log_progress)
                        if isinstance(callback_handler, CallbackHandler):
                            await callback_handler.drain()

                        # Validate and store
                        agg_result = ReactionOutputSchema.model_validate_json(agg_output)

                        # Save aggregator output for examination
                        aggregator_output_log = {
                            "step": step,
                            "aggregation_index": i + 1,
                            "k_subset_indices": subset_indices,
                            "reasoning_summary": agg_result.reasoning_summary,
                            "reactants_smiles": agg_result.reactants_smiles_list,
                            "products_smiles": agg_result.products_smiles_list,
                            "full_output": json.loads(agg_output)
                        }
                        with open(f"{rsa_log_dir}/aggregator_step{step}_{i+1:02d}_output.json", "w") as f:
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
                    await clogger.warning(
                        f"All aggregations at step {step} failed, using previous step results"
                    )
                    break

                current_proposals = next_proposals
                await clogger.info(f"Step {step} produced {len(current_proposals)} aggregated proposals")

            # Select final output (use first/best from final step)
            if current_proposals:
                final_proposal = current_proposals[0]
                output = final_proposal["output"]

                # Save final output for examination
                final_result = json.loads(output)
                final_output_log = {
                    "final_step": rsa_t,
                    "mode": rsa_mode,
                    "n_proposals": rsa_n,
                    "k_subset_size": rsa_k,
                    "t_stages": rsa_t,
                    "final_reasoning": final_result.get("reasoning_summary", ""),
                    "final_reactants_smiles": final_result.get("reactants_smiles_list", []),
                    "final_products_smiles": final_result.get("products_smiles_list", []),
                    "full_output": final_result
                }
                with open(f"{rsa_log_dir}/FINAL_OUTPUT.json", "w") as f:
                    json.dump(final_output_log, f, indent=2)

                await clogger.info(f"RSA mode completed successfully. Logs saved to: {rsa_log_dir}")
            else:
                raise ValueError("RSA aggregation produced no valid results")

        except Exception as e:
            # Fallback to standard mode if RSA fails
            await clogger.error(f"RSA mode failed: {str(e)}, falling back to standard retrosynthesis")
            # Reset task to original
            runner.task = retro_task
            output = await runner.run(log_progress)
            if isinstance(callback_handler, CallbackHandler):
                await callback_handler.drain()
    else:
        # Standard mode (no RSA)
        output = await runner.run(log_progress)
        if isinstance(callback_handler, CallbackHandler):
            await callback_handler.drain()

    if os.getenv("CHARGE_DISABLE_OUTPUT_VALIDATION", "0") == "1":
        await clogger.warning(
            "Structure validation disabled for RetrosynthesisTask output schema."
            "Returning text results without validation first before post-processing."
        )

    result = ReactionOutputSchema.model_validate_json(output)

    await highlight_node(current_node, websocket, False)

    level = current_node.level + 1

    reasoning_summary = result.reasoning_summary
    await clogger.info(
        f"Retrosynthesis reasoning summary for {current_node.smiles}:\n{reasoning_summary}",
    )
    await asyncio.sleep(0)

    # Add reaction alternative
    ai_alternative = ReactionAlternative(
        f"ai_reaction_{node_id}",
        "AI-Generated Reaction",
        "ai",
        "active",
        [
            PathwayStep([current_node.smiles], [current_node.label], [-1]),
            PathwayStep(
                list(result.reactants_smiles_list),
                [
                    smiles_to_html(s, run_settings.molecule_name_format)
                    for s in result.reactants_smiles_list
                ],
                [0 for _ in range(len(result.reactants_smiles_list))],
            ),
        ],
        reasoning_summary,
        disabled=False,
    )
    if current_node.reaction is not None and not isinstance(
        current_node.reaction, dict
    ):
        if current_node.reaction.alternatives is None:
            current_node.reaction.alternatives = []
        for alt in current_node.reaction.alternatives:
            if alt.status == "active":
                alt.status = "available"
        current_node.reaction.alternatives.insert(0, ai_alternative)
        alternatives = current_node.reaction.alternatives
        templates_searched = current_node.reaction.templatesSearched
    else:
        alternatives = [ai_alternative]
        templates_searched = False

    ai_reaction = Reaction(
        "ai_reaction_0",
        reasoning_summary,
        highlight="red",
        label="FLASK AI",
        alternatives=alternatives,
        templatesSearched=templates_searched,
    )
    ai_reaction.mappedReaction = build_mapped_reaction_dict_or_none(
        reactants=list(result.reactants_smiles_list),
        products=[current_node.smiles],
        log_msg="Failed to build rdkitjs mapped reaction for AI reaction node_id={node_id} product_smiles={smiles}",
        node_id=current_node.id,
        smiles=current_node.smiles,
    )
    current_node.reaction = ai_reaction

    # Update node with discovered reaction
    context.node_id_to_reasoning_summary[node_id] = reasoning_summary
    await context.update_node(current_node, websocket)
    await asyncio.sleep(0)

    context.recalculate_nodes_per_level()

    new_nodes: list[Node] = []
    purchasable: list[bool] = []
    for smiles in result.reactants_smiles_list:
        node_id_str = context.new_node_id()
        mol_sources = is_purchasable(smiles)
        if mol_sources:
            purchasable_str = f"Yes (via {', '.join(mol_sources)})"
        else:
            purchasable_str = "No"
        node = Node(
            node_id_str,
            smiles,
            smiles_to_html(smiles, run_settings.molecule_name_format),
            f"Discovered by {runner.model}.\n\n**Purchasable**? {purchasable_str}",
            level,
            current_node.id,
            purchasable=(len(mol_sources) > 0),
        )
        new_nodes.append(node)
        purchasable.append(len(mol_sources) > 0)
        # Add and stream node directly
        await context.add_node(node, websocket)
        await asyncio.sleep(0)

    # Template-based expansion of reactants (optional)
    if not os.path.exists(config_file):
        await clogger.info(
            f"Template-based retrosynthesis config not found at {config_file}. "
            "Skipping template expansion. AI-based retrosynthesis completed successfully."
        )
    else:
        for node, purch in zip(new_nodes, purchasable):
            if purch:  # Skip purchasable nodes unless explicitly asked for
                continue

            # Highlight node because we are looking for templates
            await highlight_node(node, websocket, True)

            # Find paths for the leaf nodes
            try:
                reaction, routes = await run_retro_planner(
                    config_file, node.smiles, clogger, run_settings
                )
                if reaction is None:
                    await clogger.warning(f"No routes found for {node.smiles}. Skipping...")
                    continue
                node.reaction = reaction

                # Use child nodes and edges of the first route
                await generate_nodes_for_molecular_graph(
                    routes[0].nodes,
                    context,
                    websocket,
                    start_level=level,
                    include_root_node=False,
                    root_node_id=node.id,
                )
            except Exception as e:
                await clogger.warning(
                    f"Template-based expansion failed for {node.smiles}: {str(e)}. Continuing..."
                )

            # Attach mapped reaction for immediate reactants -> product.
            # This is needed so hover-highlighting works for the first template step
            # discovered after an AI-generated step.
            if node.reaction is not None:
                child_smiles = [
                    n.smiles
                    for nid, n in context.node_ids.items()
                    if context.parents.get(nid) == node.id
                ]
                node.reaction.mappedReaction = build_mapped_reaction_dict_or_none(
                    reactants=child_smiles,
                    products=[node.smiles],
                    log_msg="Failed to build rdkitjs mapped reaction for template node_id={node_id} smiles={smiles}",
                    node_id=node.id,
                    smiles=node.smiles,
                )
            await context.update_node(node, websocket)  # Also disables highlight

    await websocket.send_json({"type": "complete"})
