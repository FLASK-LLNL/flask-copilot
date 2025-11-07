################################################################################
## Copyright 2025 Lawrence Livermore National Security, LLC. and Binghamton University.
## See the top-level LICENSE file for details.
##
## SPDX-License-Identifier: Apache-2.0
################################################################################

import click
import sys

from charge.servers.server_utils import update_mcp_network, get_hostname
from tool_registration import register_tool_server
from loguru import logger
from rdkit import Chem


@click.command()
@click.option("--port", type=int, default=8124, help="Port to run the server on")
@click.option("--host", type=str, default=None, help="Host to run the server on")
@click.option("--name", type=str, default="lmo_tools", help="Name of the MCP server")
@click.option(
    "--copilot-port", type=int, default=8001, help="Port to the running copilot backend"
)
@click.option(
    "--copilot-host", type=str, default=None, help="Host to the running copilot backend"
)
@click.option(
    "--model",
    type=str,
    default="gpt-5-nano",
    help="Model to use for the LMO tool server diagnose functions",
)
@click.option(
    "--backend",
    type=str,
    default="openai",
    help="Backend to use for the LMO tool server diagnose functions",
)
@click.option(
    "--api-key", type=str, default=None, help="API key for the LMO tool server"
)
@click.option(
    "--base-url", type=str, default=None, help="Base URL for the LMO tool server"
)
@click.pass_context
def main(
    ctx, port, host, name, copilot_port, copilot_host, api_key, base_url, model, backend
):
    if host is None:
        _, host = get_hostname()

    register_tool_server(port, host, name, copilot_port, copilot_host)

    sys.argv = [sys.argv[0]] + ctx.args + [f"--port={port}", f"--host={host}"]

    import charge.servers.molecular_generation_server as LMO_MCP

    LMO_MCP.setup_autogen_pool(
        api_key=api_key,
        base_url=base_url,
        model=model,
        backend=backend,
    )

    mcp = LMO_MCP.mcp

    update_mcp_network(mcp, host, port)

    mcp.run(transport="sse")


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


@click.pass_context
def main(ctx, port, host, name, copilot_port, copilot_host):
    if host is None:
        _, host = get_hostname()

    register_tool_server(port, host, name, copilot_port, copilot_host)

    sys.argv = [sys.argv[0]] + ctx.args + [f"--port={port}", f"--host={host}"]

    import charge.servers.molecular_generation_server as LMO_MCP

    mcp = LMO_MCP.mcp
    mcp.tool()(is_already_known)

    update_mcp_network(mcp, host, port)

    mcp.run(transport="sse")


if __name__ == "__main__":
    main()
