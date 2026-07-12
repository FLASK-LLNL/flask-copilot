################################################################################
## Copyright 2025-2026 Lawrence Livermore National Security, LLC..
## See the top-level LICENSE file for details.
##
## SPDX-License-Identifier: Apache-2.0
################################################################################

import re

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


async def _expand_node(
    parent: Node,
    reactant_smiles: list[str],
    reagent_smiles: list[str],
    reaction_smiles: str,
    context: GraphContext,
    websocket: WebSocket,
    run_settings: FlaskRunSettings,
) -> list[Node]:
    """Attach a reaction's reactants/reagents as children of ``parent``.

    Adds one child node per reactant (and per reagent, labeled "(Reagent)"),
    one level below ``parent``, and attaches a user-supplied Reaction with a
    mapped reaction to ``parent`` for hover-highlighting.

    :return: The list of child nodes created.
    """
    children: list[Node] = []

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
            level=parent.level + 1,
            parentId=parent.id,
            purchasable=(len(child_sources) > 0),
        )
        await context.add_node(node, websocket)
        children.append(node)

    for smiles in reactant_smiles:
        await _add_child(smiles, "Reactant")
    for smiles in reagent_smiles:
        await _add_child(smiles, "Reagent")

    # Attach the user-supplied reaction to the parent, with a mapped reaction so
    # hover-highlighting works. templatesSearched=False so the user can still
    # expand children later with templates/AI.
    parent.reaction = Reaction(
        "user_rxn",
        f"User-provided reaction\n\n**Reaction SMILES:** `{reaction_smiles}`",
        highlight="yellow",
        label="User",
        templatesSearched=False,
    )
    parent.reaction.mappedReaction = build_mapped_reaction_dict_or_none(
        reactants=reactant_smiles + reagent_smiles,
        products=[parent.smiles],
        log_msg="Failed to build rdkitjs mapped reaction for user reaction node_id={node_id} smiles={smiles}",
        node_id=parent.id,
        smiles=parent.smiles,
    )
    await context.update_node(parent, websocket)
    return children


def _leaf_nodes_by_smiles(context: GraphContext) -> dict[str, Node]:
    """Map canonical SMILES -> leaf Node (nodes with no children)."""
    parents = set(context.parents.values())
    leaves: dict[str, Node] = {}
    for nid, node in context.node_ids.items():
        if nid in parents:
            continue
        try:
            canon = Chem.MolToSmiles(Chem.MolFromSmiles(node.smiles))
        except Exception:
            canon = node.smiles
        leaves[canon] = node
    return leaves


async def reaction_smiles_retrosynthesis(
    reaction_smiles: str,
    context: GraphContext,
    websocket: WebSocket,
    run_settings: FlaskRunSettings,
    send_complete: bool = True,
) -> None:
    """Build a partial retrosynthesis graph from one or more reaction SMILES.

    A single reaction becomes a one-step graph: the product (largest, by
    heavy-atom count) is the root target node and each reactant/reagent is a
    level-1 child, as if the user had performed a single retrosynthesis step.

    Multiple reactions may be supplied, one per line (newline- or tab-
    separated). The first line builds the root reaction. Each subsequent line
    expands the graph one step deeper: if any of its products matches an
    existing leaf node (by canonical SMILES), that leaf's reactants/reagents are
    attached as its children. Lines that fail to parse, contain no reactants, or
    whose products match no current leaf are skipped with a warning.

    No route search is performed -- the user supplies the reactions.

    :param send_complete: If True (default), send a ``{"type": "complete"}``
        message when finished. Set False when a caller (e.g. the custom-problem
        flow) will run further work and send ``complete`` itself.
    """
    clogger = CallbackLogger(websocket, source="reaction_smiles_retrosynthesis")

    # Split into individual reactions on newlines/tabs; ignore blank lines.
    lines = [ln.strip() for ln in re.split(r"[\n\r\t]+", reaction_smiles) if ln.strip()]

    async def _finish() -> None:
        if send_complete:
            await websocket.send_json({"type": "complete"})

    if not lines:
        await clogger.error("No reaction SMILES provided.")
        await _finish()
        return

    context.reset()

    # --- First line: build the root reaction. ---
    root_line = lines[0]
    await clogger.info(f"Building partial graph from reaction: `{root_line}`.")
    try:
        reactant_smiles, reagent_smiles, product_smiles = _split_reaction_smiles(
            root_line
        )
    except ValueError as e:
        await clogger.error(f"Could not parse reaction SMILES: {e}")
        await _finish()
        return

    if not reactant_smiles or not product_smiles:
        await clogger.error(
            "Reaction SMILES must contain at least one reactant and one product "
            "(format: `reactant1.reactant2>>product`)."
        )
        await _finish()
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
    await _expand_node(
        root,
        reactant_smiles,
        reagent_smiles,
        root_line,
        context,
        websocket,
        run_settings,
    )

    # --- Subsequent lines: expand a matching leaf one step deeper. ---
    for line in lines[1:]:
        try:
            r_mols, rg_mols, p_mols = parse_reaction_smiles(line, include_reagents=True)
        except ValueError as e:
            await clogger.warning(f"Skipping unparseable reaction `{line}`: {e}")
            continue
        if not r_mols or not p_mols:
            await clogger.warning(
                f"Skipping reaction `{line}`: needs at least one reactant and "
                "one product."
            )
            continue

        # Attach to the first leaf matched by ANY of this line's products.
        leaves = _leaf_nodes_by_smiles(context)
        product_canons = [Chem.MolToSmiles(m) for m in p_mols]
        target = next((leaves[pc] for pc in product_canons if pc in leaves), None)
        if target is None:
            await clogger.warning(
                f"Skipping reaction `{line}`: no product matches an existing "
                "leaf molecule in the graph."
            )
            continue

        await _expand_node(
            target,
            [Chem.MolToSmiles(m) for m in r_mols],
            [Chem.MolToSmiles(m) for m in rg_mols],
            line,
            context,
            websocket,
            run_settings,
        )

    await _finish()
