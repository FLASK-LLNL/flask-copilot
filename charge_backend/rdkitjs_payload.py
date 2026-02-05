from __future__ import annotations

from dataclasses import asdict, dataclass
from typing import Any, Dict, List, Optional, Sequence, Tuple, Union

from rdkit import Chem

import rdkit_mol_differ as mol_differ

from backend_helper_funcs import RdkitjsMolPayload


MolInput = Union[str, Chem.Mol]


@dataclass(frozen=True)
class RdkitjsReactionPayload:
    version: str
    reactants: List[RdkitjsMolPayload]
    products: List[RdkitjsMolPayload]
    main_product_index: int
    highlight_rgb: Tuple[int, int, int]
    highlight_alpha: float
    reactant_mcs_smarts: List[Optional[str]]


def _mol_payload(
    m: Chem.Mol,
    *,
    highlight_atom_idxs: Sequence[int],
) -> RdkitjsMolPayload:
    # Do not set atom-map numbers on atoms: RDKit depictions may render them as labels.
    molblock = Chem.MolToMolBlock(m)
    smiles = Chem.MolToSmiles(m)
    idxs = sorted({int(i) for i in highlight_atom_idxs})
    mapnums = [int(i) + 1 for i in idxs]
    return RdkitjsMolPayload(
        molblock=str(molblock),
        smiles=str(smiles),
        highlight_atom_idxs=idxs,
        highlight_atom_mapnums=mapnums,
    )


def build_rdkitjs_mapped_reaction(
    *,
    reactants: Sequence[MolInput],
    products: Sequence[MolInput],
    highlight_rgb: Tuple[int, int, int] = (220, 0, 0),
    highlight_alpha: float = 0.45,
    mcs_timeout_s: int = 2,
) -> RdkitjsReactionPayload:
    """Compute changed atoms and package MolBlocks + highlight indices for rdkit.js."""

    r_mols = mol_differ._ensure_mols(reactants, label="reactants")
    p_mols = mol_differ._ensure_mols(products, label="products")

    changes = mol_differ.reaction_atom_changes_from_lists(
        r_mols,
        p_mols,
        mcs_timeout_s=int(mcs_timeout_s),
    )

    r_payloads: List[RdkitjsMolPayload] = []
    for m, hl in zip(r_mols, changes.reactant_changed_atoms):
        r_payloads.append(_mol_payload(m, highlight_atom_idxs=sorted(hl)))

    p_payloads: List[RdkitjsMolPayload] = []
    for m, hl in zip(p_mols, changes.product_changed_atoms):
        p_payloads.append(_mol_payload(m, highlight_atom_idxs=sorted(hl)))

    return RdkitjsReactionPayload(
        version="rdkitjs-reaction-payload/v1",
        reactants=r_payloads,
        products=p_payloads,
        main_product_index=int(changes.main_product_index),
        highlight_rgb=(int(highlight_rgb[0]), int(highlight_rgb[1]), int(highlight_rgb[2])),
        highlight_alpha=float(highlight_alpha),
        reactant_mcs_smarts=list(changes.reactant_mcs_smarts),
    )


def reaction_payload_to_json_dict(payload: RdkitjsReactionPayload) -> Dict[str, Any]:
    """Return a JSON-serializable dict."""
    return asdict(payload)


def mapped_reaction_to_json_dict(mapped_reaction: RdkitjsReactionPayload) -> Dict[str, Any]:
    """Return a JSON-serializable dict (preferred name)."""
    return asdict(mapped_reaction)
