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
from charge_backend.rdkit_mol_differ import parse_reaction_smiles


def _select_root_product(products: list[Chem.Mol]) -> Chem.Mol:
    """Pick the product with the most heavy atoms as the target molecule.

    Additional products are treated as byproducts and ignored.
    """
    return max(products, key=lambda m: m.GetNumHeavyAtoms())


async def reaction_smiles_retrosynthesis(
    reaction_smiles: str,
    context: GraphContext,
    websocket: WebSocket,
    run_settings: FlaskRunSettings,
) -> None:
    """Build a one-step partial retrosynthesis graph from a reaction SMILES.

    The product (largest, by heavy-atom count) becomes the root target node and
    each reactant becomes a level-1 child, as if the user had performed a single
    retrosynthesis step. No route search is performed -- the user supplies the
    reaction.
    """
    clogger = CallbackLogger(websocket, source="reaction_smiles_retrosynthesis")
    await clogger.info(f"Building partial graph from reaction: `{reaction_smiles}`.")

    context.reset()

    # Parse the reaction SMILES (agents: e.g. solvent, catalyst, if any, are ignored).
    try:
        reactant_mols, product_mols = parse_reaction_smiles(reaction_smiles)
    except ValueError as e:
        await clogger.error(f"Could not parse reaction SMILES: {e}")
        await websocket.send_json({"type": "complete"})
        return

    if not reactant_mols or not product_mols:
        await clogger.error(
            "Reaction SMILES must contain at least one reactant and one product "
            "(format: `reactant1.reactant2>>product`)."
        )
        await websocket.send_json({"type": "complete"})
        return

    product_mol = _select_root_product(product_mols)
    product_smiles = Chem.MolToSmiles(product_mol)
    all_product_smiles = [Chem.MolToSmiles(m) for m in product_mols]
    reactant_smiles = [Chem.MolToSmiles(m) for m in reactant_mols]

    if len(product_mols) > 1:
        await clogger.info(
            f"Multiple products found in {all_product_smiles}; using `{product_smiles}` as the target "
            "and treating the rest as byproducts."
        )

    # Root node = product.
    mol_sources = is_purchasable(product_smiles)
    root = Node(
        id="node_0",
        smiles=product_smiles,
        label=smiles_to_html(product_smiles, run_settings.molecule_name_format),
        hoverInfo=f"""# Root molecule
**SMILES:** {product_smiles}

**Purchasable**? {purchasable_summary(mol_sources)}""",
        level=0,
        parentId=None,
        purchasable=(len(mol_sources) > 0),
        highlight="yellow",
        x=100,
        y=100,
    )
    await context.add_node(root, websocket=websocket)

    # Child nodes = reactants.
    for smiles in reactant_smiles:
        child_sources = is_purchasable(smiles)
        node = Node(
            id=context.new_node_id(),
            smiles=smiles,
            label=smiles_to_html(smiles, run_settings.molecule_name_format),
            hoverInfo=f"""# Reactant
  * SMILES: {smiles}
  * Purchasable? {purchasable_summary(child_sources)}
""",
            level=1,
            parentId=root.id,
            purchasable=(len(child_sources) > 0),
        )
        await context.add_node(node, websocket)

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
        reactants=reactant_smiles,
        products=[product_smiles],
        log_msg="Failed to build rdkitjs mapped reaction for user reaction root_id={node_id} smiles={smiles}",
        node_id=root.id,
        smiles=root.smiles,
    )
    await context.update_node(root, websocket)
    await websocket.send_json({"type": "complete"})
