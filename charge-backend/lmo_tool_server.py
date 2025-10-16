################################################################################
## Copyright 2025 Lawrence Livermore National Security, LLC. and Binghamton University.
## See the top-level LICENSE file for details.
##
## SPDX-License-Identifier: Apache-2.0
################################################################################


from charge.servers.server_utils import update_mcp_network, get_hostname
from charge.servers.SMILES import SMILES_mcp
from charge.servers.molecular_property_utils import chemprop_preds_server
from charge.servers.log_progress import log_progress
import argparse
from charge.servers.server_utils import add_server_arguments
import os
from loguru import logger
from charge.servers import SMILES_utils
import json
from rdkit import Chem

# SMILES_mcp.tool()(chemprop_preds_server)
from mcp.server.fastmcp import FastMCP

parser = argparse.ArgumentParser()
add_server_arguments(parser)
args = parser.parse_args()


JSON_FILE_PATH = f"{os.getcwd()}/known_molecules.json"

mcp = FastMCP(
    "SMILES Diagnosis and retrieval MCP Server",
    port=args.port,
    website_url=f"{args.host}",
)


@mcp.tool()
def is_already_known(smiles: str) -> bool:
    """
    Check if a SMILES string provided is already known. Only provide
    valid SMILES strings. Returns True if the SMILES string is valid, and
    already in the database, False otherwise.
    Args:
        smiles (str): The input SMILES string.
    Returns:
        bool: True if the SMILES string is valid and known, False otherwise.

    Raises:
        ValueError: If the SMILES string is invalid.
    """
    if not Chem.MolFromSmiles(smiles):
        raise ValueError("Invalid SMILES string.")

    try:
        canonical_smiles = SMILES_utils.canonicalize_smiles(smiles)

        try:
            with open(JSON_FILE_PATH) as f:
                known_mols = json.load(f)
                known_smiles = [mol["smiles"] for mol in known_mols]

        except FileNotFoundError:
            logger.warning(f"{JSON_FILE_PATH} not found. Creating a new one.")
            known_mols = []

    except Exception as e:
        raise ValueError("Error in canonicalizing SMILES string.") from e

    # Check if the SMILES string is already known (in the database)
    # This is a placeholder for the actual database check
    return canonical_smiles in known_smiles


@mcp.tool()
def get_density(smiles: str) -> float:
    """
    Predict molecular properties using pre-trained Chemprop models.
    This function returns property predictions from Chemprop models. It validates the requested property name,
    constructs the appropriate model, and returns predictions for the provided SMILES input.

    Parameters
    ----------
    smiles : str
        A SMILES string representing the molecule to be evaluated.
    property : str
        The property to predict. Must be one of the valid property names listed above.

    Returns
    -------
    float
        A float representing the predicted value for the specified property.

    Raises
    ------
    ValueError
        If the smiles string is invalid
    """
    return chemprop_preds_server(smiles, property="density")


mcp.tool()(log_progress)
mcp.tool()(SMILES_utils.canonicalize_smiles)
mcp.tool()(SMILES_utils.verify_smiles)
if __name__ == "__main__":
    mcp.run(transport="sse")
