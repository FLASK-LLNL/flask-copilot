from rdkit import Chem

import rdkit_mol_differ as md


def test_sn2_like_ether_formation_highlights_leaving_group_and_new_bond_site():
    # CCCCCCl + CCO -> CCOCCCCC
    ch = md.reaction_atom_changes_from_lists(
        reactants=["CCCCCCl", "CCO"],
        products=["CCOCCCCC"],
    )

    assert ch.main_product_index == 0

    # Ethanol: only O should be changed
    ethanol_changed = ch.reactant_changed_atoms[1]
    em = Chem.MolFromSmiles("CCO")
    assert {em.GetAtomWithIdx(i).GetSymbol() for i in ethanol_changed} == {"O"}

    # Alkyl halide: Br/Cl equivalent - leaving group and adjacent carbon
    halide_changed = ch.reactant_changed_atoms[0]
    hm = Chem.MolFromSmiles("CCCCCCl")
    assert {hm.GetAtomWithIdx(i).GetSymbol() for i in halide_changed} == {"C", "Cl"}
