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

    def is_path_pattern(m: Chem.Mol) -> bool:
        if m is None:
            return False
        na = int(m.GetNumAtoms())
        nb = int(m.GetNumBonds())
        if na <= 1:
            return False
        if nb != na - 1:
            return False
        for a in m.GetAtoms():
            if int(a.GetDegree()) > 2:
                return False
        return True

    def outside_degree_list(m: Chem.Mol, match: Tuple[int, ...]) -> List[int]:
        inside = set(int(i) for i in match)
        out: List[int] = []
        for i in match:
            a = m.GetAtomWithIdx(int(i))
            cnt = 0
            for nb in a.GetNeighbors():
                if int(nb.GetIdx()) not in inside:
                    cnt += 1
            out.append(int(cnt))
        return out

    patt_is_path = is_path_pattern(patt)

    best_r: Tuple[int, ...] = ()
    best_p: Tuple[int, ...] = ()
    best_score = -1
    best_new = -1

    for rm in r_matches:
        r_orients = [rm]
        if patt_is_path:
            r_orients.append(tuple(reversed(rm)))

        for pm in p_matches:
            p_orients = [pm]
            if patt_is_path:
                p_orients.append(tuple(reversed(pm)))

            for r2 in r_orients:
                r_out = outside_degree_list(reactant, r2)
                for p2 in p_orients:
                    p_out = outside_degree_list(product, p2)
                    score = sum(1 for a, b in zip(r_out, p_out) if a == b)
                    new_cnt = sum(1 for a in p2 if int(a) not in product_atoms_already_mapped)
                    if score > best_score or (score == best_score and new_cnt > best_new):
                        best_score = score
                        best_new = new_cnt
                        best_r = r2
                        best_p = p2

    if not best_r or not best_p:
        return r_matches[0], p_matches[0]
    return best_r, best_p


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


def _orient_match_pair(
    patt: Chem.Mol,
    reactant: Chem.Mol,
    product: Chem.Mol,
    r_match: Tuple[int, ...],
    p_match: Tuple[int, ...],
    product_atoms_already_mapped: Set[int],
) -> Tuple[Tuple[int, ...], Tuple[int, ...]]:
    """Choose an orientation for (r_match, p_match) for path-like patterns.

    RDKit substructure matches do not encode an intrinsic direction. For
    chain/path MCS patterns, the match orientation can flip, and if we zip the
    raw tuples we may align the "reacting end" incorrectly.

    Heuristic: pick the orientation pairing that best matches each atom's
    outside-degree signature (connections to atoms outside the match).
    """

    def is_path_pattern(m: Chem.Mol) -> bool:
        if m is None:
            return False
        na = int(m.GetNumAtoms())
        nb = int(m.GetNumBonds())
        if na <= 1:
            return False
        if nb != na - 1:
            return False
        for a in m.GetAtoms():
            if int(a.GetDegree()) > 2:
                return False
        return True

    def outside_degree_list(m: Chem.Mol, match: Tuple[int, ...]) -> List[int]:
        inside = set(int(i) for i in match)
        out: List[int] = []
        for i in match:
            a = m.GetAtomWithIdx(int(i))
            cnt = 0
            for nb in a.GetNeighbors():
                if int(nb.GetIdx()) not in inside:
                    cnt += 1
            out.append(int(cnt))
        return out

    if not is_path_pattern(patt):
        return r_match, p_match

    r_orients = [r_match, tuple(reversed(r_match))]
    p_orients = [p_match, tuple(reversed(p_match))]

    best_r = r_match
    best_p = p_match
    best_score = -1
    best_new = -1

    for r2 in r_orients:
        r_out = outside_degree_list(reactant, r2)
        for p2 in p_orients:
            p_out = outside_degree_list(product, p2)
            score = sum(1 for a, b in zip(r_out, p_out) if a == b)
            new_cnt = sum(1 for a in p2 if int(a) not in product_atoms_already_mapped)
            if score > best_score or (score == best_score and new_cnt > best_new):
                best_score = score
                best_new = new_cnt
                best_r = r2
                best_p = p2

    return best_r, best_p


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

    This is intentionally conservative: it should not extend through a
    newly-formed bond to map an entire chain across the reaction.
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

                # Only map if the local outside-degree matches as well.
                # This avoids extending mappings through newly-formed bonds.
                r_out = sum(
                    1
                    for x in r_nb.GetNeighbors()
                    if int(x.GetIdx()) not in inst_map and int(x.GetIdx()) != int(r_ai)
                )
                p_atom2 = product.GetAtomWithIdx(int(pick))
                p_out = sum(
                    1
                    for x in p_atom2.GetNeighbors()
                    if int(x.GetIdx()) not in used_p and int(x.GetIdx()) != int(p_ai)
                )
                if int(r_out) != int(p_out):
                    continue

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

    # Highlight policy:
    # - Only highlight atoms that are local to the reaction center (bond changes),
    #   plus immediate leaving group atoms attached to those centers.
    # - Exclude carbon atoms unless a bond ORDER changes at that atom
    #   (single<->double, aromatic<->single, etc.). Bond formation/cleavage
    #   (bond present vs missing) does not qualify.
    reactant_bond_order_changed_atoms: List[Set[int]] = [set() for _ in reactants]
    product_bond_order_changed_atoms_main: Set[int] = set()
    reactant_bond_change_endpoints: List[Set[int]] = [set() for _ in reactants]
    product_bond_change_endpoints_main: Set[int] = set()
    reactant_sig_changed_atoms: List[Set[int]] = [set() for _ in reactants]
    product_sig_changed_atoms_main: Set[int] = set()

    # For product atoms that map back to reactant atoms, keep a single source:
    # (reactant_index, reactant_atom_index, reactant_instance_id)
    #
    # reactant_instance_id lets us explain multiple copies of the same reactant
    # appearing in the product even when the reactants list only contains one
    # instance (i.e., stoichiometry not explicitly provided).
    product_atom_source: Dict[int, Tuple[int, int, int]] = {}

    # Map larger reactants first to reduce ambiguous matches where a small
    # reactant (e.g., ethanol) can match multiple places in the product.
    reactant_order = sorted(
        range(len(reactants)),
        key=lambda i: int(reactants[i].GetNumHeavyAtoms()) if reactants[i] is not None else 0,
        reverse=True,
    )

    for ri in reactant_order:
        r = reactants[ri]
        patt = _mcs_pattern(r, main_product, timeout_s=int(mcs_timeout_s))
        if patt is None:
            continue

        r_match, p_best = _pick_best_match(patt, r, main_product, product_atoms_mapped)
        if not r_match:
            continue

        allow_multi = (
            int(r.GetNumHeavyAtoms()) >= 6
            and int(patt.GetNumAtoms()) >= 6
        )

        # Select multiple matches in the product to explain repeated reactant units.
        # For small reactants, avoid over-mapping due to symmetry/ambiguity.
        if allow_multi:
            p_matches = _pick_greedy_product_matches(patt, main_product, product_atoms_mapped)
        else:
            p_matches = [p_best] if p_best else []
        if not p_matches:
            continue

        # Record coverage + source for all selected matches.
        for inst_id, pm in enumerate(p_matches):
            if len(pm) != len(r_match):
                continue
            r_oriented, p_oriented = _orient_match_pair(
                patt,
                r,
                main_product,
                r_match,
                pm,
                product_atoms_mapped,
            )
            inst_map: Dict[int, int] = {
                int(r_ai): int(p_ai) for r_ai, p_ai in zip(r_oriented, p_oriented)
            }
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
        # to capture repeated units in the product. Only enable this for larger reactants.
        if allow_multi:
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
                    reactant_bond_change_endpoints[ri].update([a1, a2])
                    product_bond_change_endpoints_main.update([int(p1), int(p2)])

                # Track bond ORDER changes (but not bond formation/cleavage)
                if pb is not None and r_bt != p_bt:
                    reactant_bond_order_changed_atoms[ri].update([a1, a2])
                    product_bond_order_changed_atoms_main.update([int(p1), int(p2)])

            # Atom-level signature changes.
            for a_r, a_p in inst_map.items():
                r_sig = atom_sig(r, int(a_r), r_inside)
                p_sig = atom_sig(main_product, int(a_p), p_inside)
                if r_sig != p_sig:
                    reactant_changed[ri].add(int(a_r))
                    product_changed_main.add(int(a_p))
                    reactant_sig_changed_atoms[ri].add(int(a_r))
                    product_sig_changed_atoms_main.add(int(a_p))

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
            product_bond_change_endpoints_main.update([p1, p2])
            reactant_bond_change_endpoints[r1i].add(r1a)
            reactant_bond_change_endpoints[r2i].add(r2a)
            continue

        rb = reactants[r1i].GetBondBetweenAtoms(int(r1a), int(r2a))
        if rb is None or rb.GetBondType() != b.GetBondType():
            product_changed_main.update([p1, p2])
            reactant_changed[r1i].update([r1a, r2a])
            product_bond_change_endpoints_main.update([p1, p2])
            reactant_bond_change_endpoints[r1i].update([r1a, r2a])

        # Track bond ORDER changes (but not bond formation)
        if rb is not None and rb.GetBondType() != b.GetBondType():
            product_bond_order_changed_atoms_main.update([p1, p2])
            reactant_bond_order_changed_atoms[r1i].update([r1a, r2a])

    product_changed[main_pi] = product_changed_main

    # Restrict highlights to reaction-center neighborhood.
    # - Keep atoms that are endpoints of bond changes.
    # - Also keep atoms whose local signature changed (e.g., proton transfers,
    #   new bonds to previously-unmapped atoms).
    # - Also keep immediate neighbors of those endpoints, but only if those
    #   neighbors are currently marked as changed (this captures leaving groups
    #   like halides which are often unmapped).
    def neighborhood_keep(m: Chem.Mol, seeds: Set[int]) -> Set[int]:
        keep = set(int(i) for i in seeds)
        for ai in list(keep):
            a = m.GetAtomWithIdx(int(ai))
            for nb in a.GetNeighbors():
                keep.add(int(nb.GetIdx()))
        return keep

    # First, apply neighborhood restriction.
    for ri, r in enumerate(reactants):
        seeds = (
            set(reactant_bond_change_endpoints[ri])
            | set(reactant_sig_changed_atoms[ri])
            | set(reactant_bond_order_changed_atoms[ri])
        )
        keep = neighborhood_keep(r, seeds)
        reactant_changed[ri] = set(int(i) for i in reactant_changed[ri] if int(i) in keep)

    for pi, pm in enumerate(products):
        if pi == main_pi:
            seeds = (
                set(product_bond_change_endpoints_main)
                | set(product_sig_changed_atoms_main)
                | set(product_bond_order_changed_atoms_main)
            )
            keep = neighborhood_keep(pm, seeds)
            product_changed[pi] = set(int(i) for i in product_changed[pi] if int(i) in keep)
        else:
            # We don't currently compute bond-change endpoints for side products.
            # Suppress their highlights to avoid unrelated/spectator groups.
            product_changed[pi] = set()

    # Apply carbon-filter rule: exclude carbon unless bond ORDER changed.
    def filter_carbons(m: Chem.Mol, changed: Set[int], allow: Set[int]) -> Set[int]:
        out: Set[int] = set()
        for ai in changed:
            a = m.GetAtomWithIdx(int(ai))
            if int(a.GetAtomicNum()) == 6 and int(ai) not in allow:
                continue
            out.add(int(ai))
        return out

    for ri, r in enumerate(reactants):
        reactant_changed[ri] = filter_carbons(r, reactant_changed[ri], reactant_bond_order_changed_atoms[ri])

    for pi, pm in enumerate(products):
        if pi == main_pi:
            product_changed[pi] = filter_carbons(pm, product_changed[pi], product_bond_order_changed_atoms_main)
        else:
            # No bond-order-change tracking for non-main products; carbon highlights are removed.
            product_changed[pi] = filter_carbons(pm, product_changed[pi], set())

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


# Optional: Pillow drawing helpers retained for parity but not required by payload builder
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
    # Blend highlight colour towards white for softer highlight
    r = int(round(a * int(highlight_rgb[0]) + (1.0 - a) * 255.0))
    g = int(round(a * int(highlight_rgb[1]) + (1.0 - a) * 255.0))
    b = int(round(a * int(highlight_rgb[2]) + (1.0 - a) * 255.0))
    hl_colors = {i: (r / 255.0, g / 255.0, b / 255.0) for i in hl}
    rdMolDraw2D.PrepareAndDrawMolecule(
        drawer,
        mol,
        highlightAtoms=hl,
        highlightAtomColors=hl_colors,
    )
    drawer.FinishDrawing()
    png_bytes = drawer.GetDrawingText()
    return Image.open(io.BytesIO(png_bytes))
