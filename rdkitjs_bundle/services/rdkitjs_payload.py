from __future__ import annotations

from dataclasses import asdict, dataclass
from typing import Any, Dict, List, Optional, Sequence, Tuple, Union

from rdkit import Chem

from . import mol_differ


MolInput = Union[str, Chem.Mol]


@dataclass(frozen=True)
class RdkitjsMolPayload:
    """Best-practice payload for rdkit.js rendering.

    - Use a MolBlock to preserve atom ordering across Python -> JS.
    - Provide highlight atom indices for direct rdkit.js highlighting.
    - Also include atom-map numbers (1..N per molecule) for debugging / future use.
    """

    molblock: str
    smiles: str
    highlight_atom_idxs: List[int]
    highlight_atom_mapnums: List[int]


@dataclass(frozen=True)
class RdkitjsReactionPayload:
    version: str
    reactants: List[RdkitjsMolPayload]
    products: List[RdkitjsMolPayload]
    main_product_index: int
    highlight_rgb: Tuple[int, int, int]
    highlight_alpha: float


def _ensure_mols(items: Sequence[MolInput], *, label: str) -> List[Chem.Mol]:
    mols: List[Chem.Mol] = []
    for x in items:
        if x is None:
            continue
        if isinstance(x, Chem.Mol):
            mols.append(x)
            continue
        if isinstance(x, str):
            s = x.strip()
            if not s:
                continue
            m = Chem.MolFromSmiles(s)
            if m is None:
                raise ValueError(f"Invalid SMILES in {label}: {s}")
            mols.append(m)
            continue
        raise TypeError(f"Unsupported {label} item type: {type(x)!r}")
    return mols


def _with_atom_mapnums(m: Chem.Mol) -> Chem.Mol:
    """Return a copy of m with atom-map numbers set to 1..N."""

    mc = Chem.Mol(m)
    for i, a in enumerate(mc.GetAtoms()):
        a.SetAtomMapNum(int(i) + 1)
    return mc


def _mol_payload(
    m: Chem.Mol,
    *,
    highlight_atom_idxs: Sequence[int],
) -> RdkitjsMolPayload:
    mm = _with_atom_mapnums(m)
    molblock = Chem.MolToMolBlock(mm)
    smiles = Chem.MolToSmiles(mm)
    idxs = sorted({int(i) for i in highlight_atom_idxs})
    mapnums = [int(i) + 1 for i in idxs]
    return RdkitjsMolPayload(
        molblock=str(molblock),
        smiles=str(smiles),
        highlight_atom_idxs=idxs,
        highlight_atom_mapnums=mapnums,
    )


def build_rdkitjs_reaction_payload(
    *,
    reactants: Sequence[MolInput],
    products: Sequence[MolInput],
    highlight_rgb: Tuple[int, int, int] = (220, 0, 0),
    highlight_alpha: float = 0.45,
    mcs_timeout_s: int = 2,
) -> RdkitjsReactionPayload:
    """Compute changed atoms and package MolBlocks + highlight indices for rdkit.js."""

    r_mols = _ensure_mols(reactants, label="reactants")
    p_mols = _ensure_mols(products, label="products")

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
    )


def reaction_payload_to_json_dict(payload: RdkitjsReactionPayload) -> Dict[str, Any]:
    """Return a JSON-serializable dict."""

    return asdict(payload)
