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

@click.command()
@click.option("--port", type=int, default=8124, help="Port to run the server on")
@click.option("--host", type=str, default=None, help="Host to run the server on")
@click.option("--name", type=str, default="lmo_tools", help="Name of the MCP server")
@click.option("--copilot-port", type=int, default=8001, help="Port to the running copilot backend")
@click.option("--copilot-host", type=str, default=None, help="Host to the running copilot backend")
@click.pass_context
def main(ctx, port, host, name, copilot_port, copilot_host):
    if host is None:
        _, host = get_hostname()

    register_tool_server(port, host, name, copilot_port, copilot_host)

    sys.argv = [sys.argv[0]] + ctx.args + [f"--port={port}", f"--host={host}"]

    import charge.servers.molecular_generation_server as LMO_MCP
    mcp = LMO_MCP.mcp

    update_mcp_network(mcp, host, port)

    mcp.run(transport="sse")

if __name__ == "__main__":
    main()
