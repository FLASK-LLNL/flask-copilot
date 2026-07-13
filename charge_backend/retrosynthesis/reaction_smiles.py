################################################################################
## Copyright 2025-2026 Lawrence Livermore National Security, LLC..
## See the top-level LICENSE file for details.
##
## SPDX-License-Identifier: Apache-2.0
################################################################################

from fastapi import WebSocket

from rdkit import Chem

from lc_conductor.callback_logger import CallbackLogger

from charge_backend.backend_helper_funcs import (
    Node,
    Reaction,
    FlaskRunSettings,
)
from charge_backend.flask_experiment import GraphContext
from charge_backend.moleculedb.molecule_naming import smiles_to_html
from charge_backend.moleculedb.purchasable import (
    is_purchasable,
    purchasable_summary,
)
from charge_backend.retrosynthesis.mapping import build_mapped_reaction_dict_or_none


def _split_reaction_smiles(
    reaction_smiles: str,
) -> tuple[list[str], list[str], list[str]]:
    """Split a reaction SMILES into (reactants, reagents, products).

    A reaction SMILES has exactly two '>' separators
    (``reactants>reagents>products``); the reagents field may be empty. Each
    field is a '.'-separated list of molecule SMILES; empty components are
    dropped.
    """
    parts = reaction_smiles.split(">")
    if len(parts) != 3:
        raise ValueError(
            "Invalid reaction SMILES (expected reactants>reagents>products): "
            f"{reaction_smiles!r}"
        )
    reactants = [s for s in parts[0].split(".") if s]
    reagents = [s for s in parts[1].split(".") if s]
    products = [s for s in parts[2].split(".") if s]
    return reactants, reagents, products


def _select_root_product(product_smiles: list[str]) -> int:
    """Return the index of the product with the most heavy atoms.

    Additional products are treated as byproducts and ignored. Products that
    fail to parse count as zero heavy atoms.
    """

    def heavy_atom_count(smiles: str) -> int:
        mol = Chem.MolFromSmiles(smiles)
        return mol.GetNumHeavyAtoms() if mol is not None else 0

    return max(
        range(len(product_smiles)),
        key=lambda i: heavy_atom_count(product_smiles[i]),
    )


async def reaction_smiles_retrosynthesis(
    reaction_smiles: str,
    context: GraphContext,
    websocket: WebSocket,
    run_settings: FlaskRunSettings,
    send_complete: bool = True,
) -> None:
    """Build a one-step partial retrosynthesis graph from a reaction SMILES.

    The product (largest, by heavy-atom count) becomes the root target node and
    each reactant becomes a level-1 child, as if the user had performed a single
    retrosynthesis step. No route search is performed -- the user supplies the
    reaction.

    :param send_complete: If True (default), send a ``{"type": "complete"}``
        message when finished. Set False when a caller (e.g. the custom-problem
        flow) will run further work and send ``complete`` itself.
    """
    clogger = CallbackLogger(websocket, source="reaction_smiles_retrosynthesis")
    await clogger.info(f"Building partial graph from reaction: `{reaction_smiles}`.")

    context.reset()

    # Parse the reaction SMILES. Reagents (e.g. solvent, catalyst) are surfaced as
    # child nodes labeled "(Reagent)", mirroring the exact-reaction database path
    # which appends each component's role to its label.
    try:
        reactant_smiles, reagent_smiles, product_smiles = _split_reaction_smiles(
            reaction_smiles
        )
    except ValueError as e:
        await clogger.error(f"Could not parse reaction SMILES: {e}")
        if send_complete:
            await websocket.send_json({"type": "complete"})
        return

    if not reactant_smiles or not product_smiles:
        await clogger.error(
            "Reaction SMILES must contain at least one reactant and one product "
            "(format: `reactant1.reactant2>>product`)."
        )
        if send_complete:
            await websocket.send_json({"type": "complete"})
        return

    root_index = _select_root_product(product_smiles)
    product = product_smiles[root_index]

    if len(product_smiles) > 1:
        await clogger.info(
            f"Multiple products found in {product_smiles}; using `{product}` as the "
            "target and treating the rest as byproducts."
        )

    # Root node = product.
    mol_sources = is_purchasable(product)
    root = Node(
        id="node_0",
        smiles=product,
        label=smiles_to_html(product, run_settings.molecule_name_format),
        hoverInfo=f"""# Root molecule
**SMILES:** {product}

**Purchasable**? {purchasable_summary(mol_sources)}""",
        level=0,
        parentId=None,
        purchasable=(len(mol_sources) > 0),
        highlight="yellow",
        x=100,
        y=100,
    )
    await context.add_node(root, websocket=websocket)

    # Child nodes = reactants, then reagents (labeled with their role).
    async def _add_child(smiles: str, role: str) -> None:
        child_sources = is_purchasable(smiles)
        label = smiles_to_html(smiles, run_settings.molecule_name_format)
        if role != "Reactant":
            label = f"{label} ({role})"
        node = Node(
            id=context.new_node_id(),
            smiles=smiles,
            label=label,
            hoverInfo=f"""# {role}
  * SMILES: {smiles}
  * Purchasable? {purchasable_summary(child_sources)}
""",
            level=1,
            parentId=root.id,
            purchasable=(len(child_sources) > 0),
        )
        await context.add_node(node, websocket)

    for smiles in reactant_smiles:
        await _add_child(smiles, "Reactant")
    for smiles in reagent_smiles:
        await _add_child(smiles, "Reagent")

    # Attach the user-supplied reaction to the root, with a mapped reaction so
    # hover-highlighting works. templatesSearched=False so the user can still
    # expand children later with templates/AI.
    root.reaction = Reaction(
        "user_rxn",
        f"User-provided reaction\n\n**Reaction SMILES:** `{reaction_smiles}`",
        highlight="yellow",
        label="User",
        templatesSearched=False,
    )
    root.reaction.mappedReaction = build_mapped_reaction_dict_or_none(
        reactants=reactant_smiles + reagent_smiles,
        products=[product],
        log_msg="Failed to build rdkitjs mapped reaction for user reaction root_id={node_id} smiles={smiles}",
        node_id=root.id,
        smiles=root.smiles,
    )
    await context.update_node(root, websocket)
    if send_complete:
        await websocket.send_json({"type": "complete"})
