from __future__ import annotations

from dataclasses import dataclass
import io
from typing import Dict, Iterable, List, Optional, Sequence, Set, Tuple

from PIL import Image, ImageDraw
from rdkit import Chem
from rdkit.Chem import rdFMCS
from rdkit.Chem.Scaffolds import MurckoScaffold
from rdkit.Chem.Draw import rdMolDraw2D


@dataclass(frozen=True)
class ReactionAtomChanges:
    """Atom-level change highlights for a reaction.

    Indices are RDKit atom indices for each molecule in the input ordering.
    """

    reactant_changed_atoms: List[Set[int]]
    product_changed_atoms: List[Set[int]]
    main_product_index: int


def _parse_reaction_smiles(reaction_smiles: str) -> Tuple[List[Chem.Mol], List[Chem.Mol]]:
    """Parse reaction SMILES: reactants>>products.

    Reactants and products are dot-separated.
    """

    if not reaction_smiles or ">>" not in reaction_smiles:
        return [], []
    left, right = reaction_smiles.split(">>", 1)

    def parse_side(side: str) -> List[Chem.Mol]:
        mols: List[Chem.Mol] = []
        for token in [t.strip() for t in side.split(".") if t.strip()]:
            m = Chem.MolFromSmiles(token)
            if m is None:
                raise ValueError(f"Invalid SMILES in reaction: {token}")
            mols.append(m)
        return mols

    return parse_side(left.strip()), parse_side(right.strip())


def _ensure_mols(items: Sequence[object], *, label: str) -> List[Chem.Mol]:
    """Convert a list of SMILES strings / RDKit mols into RDKit mols."""

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


def _largest_mol_index(mols: Sequence[Chem.Mol]) -> int:
    if not mols:
        return -1
    best_i = 0
    best_hac = -1
    for i, m in enumerate(mols):
        hac = int(m.GetNumHeavyAtoms()) if m is not None else 0
        if hac > best_hac:
            best_hac = hac
            best_i = i
    return best_i


def _mcs_pattern(m1: Chem.Mol, m2: Chem.Mol, timeout_s: int = 2) -> Optional[Chem.Mol]:
    try:
        res = rdFMCS.FindMCS(
            [m1, m2],
            timeout=int(timeout_s),
            atomCompare=rdFMCS.AtomCompare.CompareElements,
            bondCompare=rdFMCS.BondCompare.CompareOrder,
            matchValences=False,
            ringMatchesRingOnly=True,
            completeRingsOnly=False,
            maximizeBonds=True,
        )
    except Exception:
        return None
    smarts = getattr(res, "smartsString", None)
    if not smarts:
        return None
    return Chem.MolFromSmarts(smarts)


def _pick_best_match(
    patt: Chem.Mol,
    reactant: Chem.Mol,
    product: Chem.Mol,
    product_atoms_already_mapped: Set[int],
) -> Tuple[Tuple[int, ...], Tuple[int, ...]]:
    r_matches = list(reactant.GetSubstructMatches(patt))
    p_matches = list(product.GetSubstructMatches(patt))
    if not r_matches or not p_matches:
        return (), ()

    best_p = p_matches[0]
    best_new = -1
    for pm in p_matches:
        new_cnt = sum(1 for a in pm if int(a) not in product_atoms_already_mapped)
        if new_cnt > best_new:
            best_new = new_cnt
            best_p = pm
    return r_matches[0], best_p


def _pick_greedy_product_matches(
    patt: Chem.Mol,
    product: Chem.Mol,
    product_atoms_already_mapped: Set[int],
) -> List[Tuple[int, ...]]:
    """Pick multiple product matches to cover repeated units.

    Greedily selects matches that add the most previously-unmapped product atoms.
    """

    p_matches = list(product.GetSubstructMatches(patt))
    if not p_matches:
        return []

    selected: List[Tuple[int, ...]] = []
    used: Set[int] = set(product_atoms_already_mapped)

    while True:
        best: Optional[Tuple[int, ...]] = None
        best_new = 0
        for pm in p_matches:
            new_cnt = sum(1 for a in pm if int(a) not in used)
            if new_cnt > best_new:
                best_new = new_cnt
                best = pm
        if best is None or best_new <= 0:
            break
        selected.append(best)
        for a in best:
            used.add(int(a))

    return selected


def _core_query_pattern(reactant: Chem.Mol) -> Optional[Chem.Mol]:
    """Return a smaller query to map repeated reactant cores.

    This is used to avoid marking an entire repeated unit in the product as
    "changed" just because the full reactant MCS only matches once.

    Current heuristic:
    - Use Bemis-Murcko scaffold (rings + linkers).
    - Add methyl-like substituents attached to scaffold atoms (C, degree 1).
    """

    try:
        scaf = MurckoScaffold.GetScaffoldForMol(reactant)
    except Exception:
        return None
    if scaf is None or scaf.GetNumAtoms() <= 0:
        return None

    scaf_matches = reactant.GetSubstructMatches(scaf)
    if not scaf_matches:
        return None

    core_atoms: Set[int] = set(int(x) for x in scaf_matches[0])
    # Add methyl-like substituents attached to the scaffold.
    for ai in list(core_atoms):
        a = reactant.GetAtomWithIdx(int(ai))
        for nb in a.GetNeighbors():
            ni = int(nb.GetIdx())
            if ni in core_atoms:
                continue
            if int(nb.GetAtomicNum()) == 6 and int(nb.GetDegree()) == 1:
                core_atoms.add(ni)

    try:
        smarts = Chem.MolFragmentToSmarts(reactant, atomsToUse=sorted(core_atoms))
    except Exception:
        return None
    if not smarts:
        return None
    return Chem.MolFromSmarts(smarts)


def _extend_instance_map_by_neighbors(
    reactant: Chem.Mol,
    product: Chem.Mol,
    inst_map: Dict[int, int],
    *,
    mapped_product_atoms: Set[int],
    max_new: int = 32,
) -> Dict[int, int]:
    """Greedily extend an atom map by matching immediate neighbors.

    This helps keep simple substituents (e.g., methyl) from being treated as
    "new" when the attachment atom is already mapped.
    """

    if not inst_map:
        return inst_map

    new_count = 0
    while new_count < int(max_new):
        added_any = False
        current_items = list(inst_map.items())
        used_p = set(inst_map.values()) | set(mapped_product_atoms)

        for r_ai, p_ai in current_items:
            r_atom = reactant.GetAtomWithIdx(int(r_ai))
            p_atom = product.GetAtomWithIdx(int(p_ai))

            for r_nb in r_atom.GetNeighbors():
                r_ni = int(r_nb.GetIdx())
                if r_ni in inst_map:
                    continue
                r_b = reactant.GetBondBetweenAtoms(int(r_ai), int(r_ni))
                r_bt = r_b.GetBondType() if r_b is not None else None

                candidates: List[int] = []
                for p_nb in p_atom.GetNeighbors():
                    p_ni = int(p_nb.GetIdx())
                    if p_ni in used_p:
                        continue
                    if int(p_nb.GetAtomicNum()) != int(r_nb.GetAtomicNum()):
                        continue
                    p_b = product.GetBondBetweenAtoms(int(p_ai), int(p_ni))
                    p_bt = p_b.GetBondType() if p_b is not None else None
                    if r_bt != p_bt:
                        continue
                    candidates.append(p_ni)

                if not candidates:
                    continue
                if len(candidates) == 1:
                    pick = candidates[0]
                else:
                    # Choose best by simple atom-property similarity.
                    def score(p_idx: int) -> Tuple[int, int, int, int]:
                        pa = product.GetAtomWithIdx(int(p_idx))
                        return (
                            1 if int(pa.GetIsAromatic()) == int(r_nb.GetIsAromatic()) else 0,
                            1 if int(pa.GetFormalCharge()) == int(r_nb.GetFormalCharge()) else 0,
                            1 if int(pa.GetTotalNumHs()) == int(r_nb.GetTotalNumHs()) else 0,
                            -int(pa.GetDegree()),
                        )

                    pick = max(candidates, key=score)

                inst_map[int(r_ni)] = int(pick)
                used_p.add(int(pick))
                mapped_product_atoms.add(int(pick))
                new_count += 1
                added_any = True
                if new_count >= int(max_new):
                    break
            if new_count >= int(max_new):
                break

        if not added_any:
            break

    return inst_map


def _reaction_atom_changes_mols(
    reactants: Sequence[Chem.Mol],
    products: Sequence[Chem.Mol],
    *,
    mcs_timeout_s: int = 2,
) -> ReactionAtomChanges:
    if not reactants or not products:
        return ReactionAtomChanges([set() for _ in reactants], [set() for _ in products], -1)

    main_pi = _largest_mol_index(products)
    if main_pi < 0:
        return ReactionAtomChanges([set() for _ in reactants], [set() for _ in products], -1)

    main_product = products[main_pi]

    product_atoms_mapped: Set[int] = set()
    reactant_instance_maps: List[List[Dict[int, int]]] = [[] for _ in reactants]

    # For product atoms that map back to reactant atoms, keep a single source:
    # (reactant_index, reactant_atom_index, reactant_instance_id)
    #
    # reactant_instance_id lets us explain multiple copies of the same reactant
    # appearing in the product even when the reactants list only contains one
    # instance (i.e., stoichiometry not explicitly provided).
    product_atom_source: Dict[int, Tuple[int, int, int]] = {}

    for ri, r in enumerate(reactants):
        patt = _mcs_pattern(r, main_product, timeout_s=int(mcs_timeout_s))
        if patt is None:
            continue

        r_match, _ = _pick_best_match(patt, r, main_product, product_atoms_mapped)
        if not r_match:
            continue

        # Select multiple matches in the product to explain repeated reactant units.
        p_matches = _pick_greedy_product_matches(patt, main_product, product_atoms_mapped)
        if not p_matches:
            continue

        # Record coverage + source for all selected matches.
        for inst_id, pm in enumerate(p_matches):
            if len(pm) != len(r_match):
                continue
            inst_map: Dict[int, int] = {int(r_ai): int(p_ai) for r_ai, p_ai in zip(r_match, pm)}
            inst_map = _extend_instance_map_by_neighbors(
                r,
                main_product,
                inst_map,
                mapped_product_atoms=product_atoms_mapped,
            )
            reactant_instance_maps[ri].append(inst_map)
            for r_ai, p_ai in inst_map.items():
                pai = int(p_ai)
                product_atoms_mapped.add(pai)
                if pai not in product_atom_source:
                    product_atom_source[pai] = (int(ri), int(r_ai), int(inst_id))

        # If the full-reactant mapping only explains one occurrence, try mapping a smaller core
        # to capture repeated units in the product.
        core_query = _core_query_pattern(r)
        if core_query is not None:
            r_core = r.GetSubstructMatch(core_query)
            if r_core:
                extra_pm = _pick_greedy_product_matches(core_query, main_product, product_atoms_mapped)
                if extra_pm:
                    base_inst = len(reactant_instance_maps[ri])
                    for j, pm in enumerate(extra_pm):
                        if len(pm) != len(r_core):
                            continue
                        inst_id = base_inst + j
                        inst_map = {int(r_ai): int(p_ai) for r_ai, p_ai in zip(r_core, pm)}
                        inst_map = _extend_instance_map_by_neighbors(
                            r,
                            main_product,
                            inst_map,
                            mapped_product_atoms=product_atoms_mapped,
                        )
                        reactant_instance_maps[ri].append(inst_map)
                        for r_ai, p_ai in inst_map.items():
                            pai = int(p_ai)
                            product_atoms_mapped.add(pai)
                            if pai not in product_atom_source:
                                product_atom_source[pai] = (int(ri), int(r_ai), int(inst_id))

    reactant_changed: List[Set[int]] = []
    for ri, r in enumerate(reactants):
        mapped_r_atoms: Set[int] = set()
        for inst_map in reactant_instance_maps[ri]:
            mapped_r_atoms.update(inst_map.keys())
        if not mapped_r_atoms:
            reactant_changed.append(set(range(int(r.GetNumAtoms()))))
        else:
            reactant_changed.append(set(range(int(r.GetNumAtoms()))) - mapped_r_atoms)

    product_changed: List[Set[int]] = [set() for _ in products]
    for i, pm in enumerate(products):
        if i != main_pi:
            product_changed[i] = set(range(int(pm.GetNumAtoms())))

    product_changed_main: Set[int] = set(range(int(main_product.GetNumAtoms()))) - set(product_atoms_mapped)
    mapped_product_atoms_global = set(product_atoms_mapped)

    def bond_sig_to_outside(m: Chem.Mol, atom_idx: int, inside: Set[int]) -> List[Tuple[int, int]]:
        out: List[Tuple[int, int]] = []
        a = m.GetAtomWithIdx(int(atom_idx))
        for nb in a.GetNeighbors():
            j = int(nb.GetIdx())
            if j in inside:
                continue
            b = m.GetBondBetweenAtoms(int(atom_idx), int(j))
            # Keep aromatic (1.5) distinct from single (1.0).
            bo = int(round(float(b.GetBondTypeAsDouble()) * 10)) if b is not None else 0
            out.append((int(nb.GetAtomicNum()), int(bo)))
        out.sort()
        return out

    def atom_sig(m: Chem.Mol, atom_idx: int, inside: Set[int]) -> Tuple[int, int, int, Tuple[Tuple[int, int], ...]]:
        a = m.GetAtomWithIdx(int(atom_idx))
        total_h = int(a.GetTotalNumHs())
        charge = int(a.GetFormalCharge())
        aromatic = 1 if a.GetIsAromatic() else 0
        ext = tuple(bond_sig_to_outside(m, atom_idx, inside))
        return total_h, charge, aromatic, ext

    for ri, r in enumerate(reactants):
        for inst_map in reactant_instance_maps[ri]:
            if not inst_map:
                continue
            r_inside = set(int(x) for x in inst_map.keys())
            p_inside = set(int(x) for x in inst_map.values())

            # Bond changes between mapped atoms.
            for b in r.GetBonds():
                a1 = int(b.GetBeginAtomIdx())
                a2 = int(b.GetEndAtomIdx())
                if a1 not in r_inside or a2 not in r_inside:
                    continue
                p1 = inst_map.get(a1)
                p2 = inst_map.get(a2)
                if p1 is None or p2 is None:
                    continue
                pb = main_product.GetBondBetweenAtoms(int(p1), int(p2))
                r_bt = b.GetBondType()
                p_bt = pb.GetBondType() if pb is not None else None
                if pb is None or r_bt != p_bt:
                    reactant_changed[ri].update([a1, a2])
                    product_changed_main.update([int(p1), int(p2)])

            # Atom-level signature changes.
            for a_r, a_p in inst_map.items():
                r_sig = atom_sig(r, int(a_r), r_inside)
                p_sig = atom_sig(main_product, int(a_p), p_inside)
                if r_sig != p_sig:
                    reactant_changed[ri].add(int(a_r))
                    product_changed_main.add(int(a_p))

    # New bonds formed in the product between mapped atoms.
    # This catches bond formation between atoms originating from different reactants
    # (which won't appear in any single-reactant bond list).
    for b in main_product.GetBonds():
        p1 = int(b.GetBeginAtomIdx())
        p2 = int(b.GetEndAtomIdx())
        if p1 not in mapped_product_atoms_global or p2 not in mapped_product_atoms_global:
            continue

        s1 = product_atom_source.get(p1)
        s2 = product_atom_source.get(p2)
        if not s1 or not s2:
            continue

        r1i, r1a, inst1 = s1
        r2i, r2a, inst2 = s2

        # Different reactant indices OR different instances of the same reactant
        # => treat as bond formation between distinct molecules.
        if r1i != r2i or inst1 != inst2:
            product_changed_main.update([p1, p2])
            reactant_changed[r1i].add(r1a)
            reactant_changed[r2i].add(r2a)
            continue

        rb = reactants[r1i].GetBondBetweenAtoms(int(r1a), int(r2a))
        if rb is None or rb.GetBondType() != b.GetBondType():
            product_changed_main.update([p1, p2])
            reactant_changed[r1i].update([r1a, r2a])

    product_changed[main_pi] = product_changed_main

    return ReactionAtomChanges(
        reactant_changed_atoms=reactant_changed,
        product_changed_atoms=product_changed,
        main_product_index=int(main_pi),
    )


def reaction_atom_changes(reaction_smiles: str, *, mcs_timeout_s: int = 2) -> ReactionAtomChanges:
    """Identify changed atoms for reaction SMILES: reactants>>products."""

    reactants, products = _parse_reaction_smiles(reaction_smiles)
    return _reaction_atom_changes_mols(reactants, products, mcs_timeout_s=int(mcs_timeout_s))


def reaction_atom_changes_from_lists(
    reactants: Sequence[object],
    products: Sequence[object],
    *,
    mcs_timeout_s: int = 2,
) -> ReactionAtomChanges:
    """Identify changed atoms given separate reactant/product lists.

    Inputs can be lists of SMILES strings or RDKit mols.
    """

    r_mols = _ensure_mols(reactants, label="reactants")
    p_mols = _ensure_mols(products, label="products")
    return _reaction_atom_changes_mols(r_mols, p_mols, mcs_timeout_s=int(mcs_timeout_s))


def _draw_mol(
    mol: Chem.Mol,
    *,
    size: Tuple[int, int] = (300, 300),
    highlight_atoms: Optional[Iterable[int]] = None,
    highlight_rgb: Tuple[int, int, int] = (220, 0, 0),
    highlight_alpha: float = 0.45,
) -> Image.Image:
    drawer = rdMolDraw2D.MolDraw2DCairo(int(size[0]), int(size[1]))
    hl = sorted(set(int(x) for x in (highlight_atoms or [])))
    a = float(highlight_alpha)
    if a < 0.0:
        a = 0.0
    if a > 1.0:
        a = 1.0
    # RDKit MolDraw2D highlight colours don't consistently support alpha across versions.
    # To "turn down" harsh colours, blend the highlight colour towards white.
    r = int(round(a * int(highlight_rgb[0]) + (1.0 - a) * 255.0))
    g = int(round(a * int(highlight_rgb[1]) + (1.0 - a) * 255.0))
    b = int(round(a * int(highlight_rgb[2]) + (1.0 - a) * 255.0))
    hl_colors = {
        i: (r / 255.0, g / 255.0, b / 255.0) for i in hl
    }
    rdMolDraw2D.PrepareAndDrawMolecule(
        drawer,
        mol,
        highlightAtoms=hl,
        highlightAtomColors=hl_colors,
    )
    drawer.FinishDrawing()
    png_bytes = drawer.GetDrawingText()
    return Image.open(io.BytesIO(png_bytes))


def _draw_reaction_with_atom_changes_mols(
    reactants: Sequence[Chem.Mol],
    products: Sequence[Chem.Mol],
    changes: ReactionAtomChanges,
    *,
    mol_size: Tuple[int, int] = (300, 300),
    highlight_rgb: Tuple[int, int, int] = (220, 0, 0),
    highlight_alpha: float = 0.45,
) -> Image.Image:

    def render_side(mols: Sequence[Chem.Mol], changed: Sequence[Set[int]]) -> Image.Image:
        imgs: List[Image.Image] = []
        for m, hl in zip(mols, changed):
            if m is None:
                continue
            imgs.append(
                _draw_mol(
                    m,
                    size=mol_size,
                    highlight_atoms=hl,
                    highlight_rgb=highlight_rgb,
                    highlight_alpha=highlight_alpha,
                ).convert("RGBA")
            )
        if not imgs:
            return Image.new("RGBA", (int(mol_size[0]), int(mol_size[1])), (255, 255, 255, 255))

        spacing = 18
        plus_w = 36
        insert_plus = len(imgs) > 1
        total_w = sum(im.width for im in imgs) + spacing * (len(imgs) - 1)
        if insert_plus:
            total_w += (plus_w + spacing) * (len(imgs) - 1)
        max_h = max(im.height for im in imgs)

        side = Image.new("RGBA", (int(total_w), int(max_h)), (255, 255, 255, 255))
        x = 0
        for i, im in enumerate(imgs):
            y = (max_h - im.height) // 2
            side.paste(im, (int(x), int(y)))
            x += im.width
            if i < len(imgs) - 1:
                x += spacing
                if insert_plus:
                    plus = Image.new("RGBA", (plus_w, int(max_h)), (255, 255, 255, 255))
                    d = ImageDraw.Draw(plus)
                    cx = plus_w // 2
                    cy = int(max_h) // 2
                    arm = 10
                    d.line([(cx - arm, cy), (cx + arm, cy)], fill=(0, 0, 0), width=3)
                    d.line([(cx, cy - arm), (cx, cy + arm)], fill=(0, 0, 0), width=3)
                    side.paste(plus, (int(x), 0))
                    x += plus_w + spacing
        return side

    left = render_side(reactants, changes.reactant_changed_atoms)
    right = render_side(products, changes.product_changed_atoms)

    arrow_w = 80
    total_h = max(left.height, right.height)
    arrow = Image.new("RGBA", (arrow_w, int(total_h)), (255, 255, 255, 255))
    ad = ImageDraw.Draw(arrow)
    mid = int(total_h) // 2
    ad.line([(10, mid), (arrow_w - 30, mid)], fill=(0, 0, 0), width=3)
    ad.polygon([(arrow_w - 30, mid - 8), (arrow_w - 10, mid), (arrow_w - 30, mid + 8)], fill=(0, 0, 0))

    out = Image.new("RGBA", (left.width + arrow_w + right.width, int(total_h)), (255, 255, 255, 255))
    out.paste(left, (0, (int(total_h) - left.height) // 2))
    out.paste(arrow, (left.width, 0))
    out.paste(right, (left.width + arrow_w, (int(total_h) - right.height) // 2))
    return out


def draw_reaction_with_atom_changes(
    reaction_smiles: str,
    *,
    mol_size: Tuple[int, int] = (300, 300),
    highlight_rgb: Tuple[int, int, int] = (220, 0, 0),
    highlight_alpha: float = 0.45,
) -> Image.Image:
    """Render reaction SMILES with changed atoms highlighted."""

    reactants, products = _parse_reaction_smiles(reaction_smiles)
    changes = _reaction_atom_changes_mols(reactants, products)
    return _draw_reaction_with_atom_changes_mols(
        reactants,
        products,
        changes,
        mol_size=mol_size,
        highlight_rgb=highlight_rgb,
        highlight_alpha=highlight_alpha,
    )


def draw_reaction_with_atom_changes_from_lists(
    reactants: Sequence[object],
    products: Sequence[object],
    *,
    mol_size: Tuple[int, int] = (300, 300),
    highlight_rgb: Tuple[int, int, int] = (220, 0, 0),
    highlight_alpha: float = 0.45,
    mcs_timeout_s: int = 2,
) -> Image.Image:
    """Render reactant/product lists with changed atoms highlighted."""

    r_mols = _ensure_mols(reactants, label="reactants")
    p_mols = _ensure_mols(products, label="products")
    changes = _reaction_atom_changes_mols(r_mols, p_mols, mcs_timeout_s=int(mcs_timeout_s))
    return _draw_reaction_with_atom_changes_mols(
        r_mols,
        p_mols,
        changes,
        mol_size=mol_size,
        highlight_rgb=highlight_rgb,
        highlight_alpha=highlight_alpha,
    )
