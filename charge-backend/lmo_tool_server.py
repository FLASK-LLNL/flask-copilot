################################################################################
## Copyright 2025 Lawrence Livermore National Security, LLC. and Binghamton University.
## See the top-level LICENSE file for details.
##
## SPDX-License-Identifier: Apache-2.0
################################################################################


from charge.servers.server_utils import update_mcp_network, get_hostname
from charge.servers.SMILES import SMILES_mcp
from charge.servers.molecular_property_utils import chemprop_preds_server
import argparse


SMILES_mcp.tool()(chemprop_preds_server)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Run a ChARGe MCP Server")
    parser.add_argument(
        "--port", type=int, default=8124, help="Port to run the server on"
    )
    parser.add_argument(
        "--host", type=str, default=None, help="Host to run the server on"
    )
    args = parser.parse_args()
    host = args.host if args.host else get_hostname()
    port = args.port
    update_mcp_network(SMILES_mcp, host, port)

    SMILES_mcp.run(transport="sse")
