
################################################################################
## Copyright 2025 Lawrence Livermore National Security, LLC. and Binghamton University.
## See the top-level LICENSE file for details.
##
## SPDX-License-Identifier: Apache-2.0
################################################################################
"""
Molecule HTML naming utilities.
"""

try:
    from rdkit import Chem
    from rdkit.Chem import rdMolDescriptors
except ImportError:
    Chem = None
    rdMolDescriptors = None

import re
import json
import os
import requests

if os.path.exists("inchi_mapping.json"):
    with open("inchi_mapping.json", "rb") as fp:
        DATABASE = json.load(fp)
else:
    DATABASE = None

def inchi_lookup(inchi: str) -> str | None:
    if DATABASE is None:
        return None
    
    if inchi not in DATABASE:
        return None
    
    if DATABASE[inchi]["name"]:  # Canonical/brand first
        return DATABASE[inchi]["name"]
    if DATABASE[inchi]["iupac_name"]:  # IUPAC second
        return DATABASE[inchi]["iupac_name"]
    return None


def smiles_to_html(smiles: str) -> str:
    if Chem is None:
        return smiles

    mol = Chem.MolFromSmiles(smiles)
    if mol is None:  # Invalid SMILES
        return smiles

    # First, try to find a canonical or IUPAC name
    inchi = Chem.MolToInchi(mol)
    name = inchi_lookup(inchi)
    if name:
       return name
    
    # Otherwise, use RDKit for a general chemical formula
    # Get the formula
    if rdMolDescriptors is None:
        return smiles
    formula = rdMolDescriptors.CalcMolFormula(mol)
    # Find charge
    charge = Chem.GetFormalCharge(mol)
    # Replace numbers with subscripts
    formula_html = re.sub(r'(\d+)', r'<sub>\1</sub>', formula)
    # Replace isotope notation (e.g., [13C]) with superscript
    # RDKit already includes isotopes as e.g. 13CH4, so we want to superscript leading numbers
    formula_html = re.sub(r'\b(\d+)([A-Z][a-z]?)', r'<sup>\1</sup>\2', formula_html)
    # Add charge as superscript at the end (if not zero)
    if charge != 0:
        sign = '+' if charge > 0 else '-'
        abs_charge = abs(charge)
        # Use only sign for single charge, number+sign for higher
        charge_str = f"{sign}" if abs_charge == 1 else f"{abs_charge}{sign}"
        formula_html += f"<sup>{charge_str}</sup>"
    return formula_html


def smiles_to_iupac_online(smiles):
    CACTUS = "https://cactus.nci.nih.gov/chemical/structure/{0}/{1}"
    try:
        rep = "iupac_name"
        url = CACTUS.format(smiles, rep)
        response = requests.get(url)
        response.raise_for_status()
        if response.text.startswith("<"):  # HTML
            return smiles
        return response.text
    except requests.exceptions.HTTPError:
        return smiles
