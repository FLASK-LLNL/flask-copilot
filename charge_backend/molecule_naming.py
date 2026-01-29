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
import pandas as pd

from typing import Literal, TypeAlias

MolNameFormat: TypeAlias = Literal["brand", "iupac", "formula", "smiles"]

_DATABASE_PATH = os.getenv("FLASK_INCHI_DB", "/data/inchi_mapping.json")
if os.path.exists(_DATABASE_PATH):
    with open(_DATABASE_PATH, "rb") as fp:
        DATABASE = json.load(fp)
else:
    DATABASE = None


def inchi_lookup(inchi: str, prefer_iupac: bool = False) -> str | None:
    if DATABASE is None:
        return None

    if inchi not in DATABASE:
        return None

    # IUPAC first, Canonical/brand second
    if prefer_iupac:
        if DATABASE[inchi]["iupac_name"]:
            return DATABASE[inchi]["iupac_name"]
        if DATABASE[inchi]["name"]:
            return DATABASE[inchi]["name"]
    else:
        if DATABASE[inchi]["name"]:  # Canonical/brand first
            return DATABASE[inchi]["name"]
        if DATABASE[inchi]["iupac_name"]:  # IUPAC second
            return DATABASE[inchi]["iupac_name"]
    return None


def smiles_to_html(
    smiles: str,
    molecule_name_format: MolNameFormat = "brand",
) -> str:
    if molecule_name_format == "smiles":
        return smiles
    if Chem is None:
        return smiles

    mol = Chem.MolFromSmiles(smiles)
    if mol is None:  # Invalid SMILES
        return smiles

    # First, try to find a canonical or IUPAC name
    if molecule_name_format in ("brand", "iupac"):
        inchi = Chem.MolToInchi(mol)
        name = inchi_lookup(inchi, molecule_name_format == "iupac")
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
    formula_html = re.sub(r"(\d+)", r"<sub>\1</sub>", formula)
    # Replace isotope notation (e.g., [13C]) with superscript
    # RDKit already includes isotopes as e.g. 13CH4, so we want to superscript leading numbers
    formula_html = re.sub(r"\b(\d+)([A-Z][a-z]?)", r"<sup>\1</sup>\2", formula_html)
    # Add charge as superscript at the end (if not zero)
    if charge != 0:
        sign = "+" if charge > 0 else "-"
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


_STOCK_DATABASE_PATH = os.getenv("FLASK_STOCK_DB", "/data/zinc_stock.hdf5")
if os.path.exists(_STOCK_DATABASE_PATH):
    STOCK_DATABASE = pd.read_hdf(_STOCK_DATABASE_PATH, "table")
    STOCK_DATABASE.set_index("inchi_key", inplace=True)  # Optimize InChI key queries
else:
    STOCK_DATABASE = None


def is_purchasable(smiles: str) -> bool:
    if Chem is None or STOCK_DATABASE is None:
        return False
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:  # Invalid SMILES
        return False
    inchi: str = str(Chem.MolToInchi(mol))
    inchi_key = Chem.InchiToInchiKey(inchi)

    return inchi_key in STOCK_DATABASE.index
