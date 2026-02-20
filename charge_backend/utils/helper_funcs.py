# from rdkit import Chem
# from rdkit.Chem import AllChem, Descriptors
from flask_mcp.chemistry import SMILES_utils
from flask_mcp.lmo.molecular_property_utils import get_density


def post_process_smiles(smiles: str, parent_id: int, node_id: int) -> dict:
    """
    Post-process a solution SMILES string, add additional properties and return
    a dictionary that can be appended to the known molecules JSON file.

    Args:
        smiles (str): The input SMILES string.
    Returns:
        dict: The post-processed dictionary.
    """
    canonical_smiles = SMILES_utils.canonicalize_smiles(smiles)
    sascore = SMILES_utils.get_synthesizability(canonical_smiles)
    density = get_density(canonical_smiles)

    return {"smiles": canonical_smiles, "sascore": sascore, "density": density}
