################################################################################
## Copyright 2025 Lawrence Livermore National Security, LLC. and Binghamton University.
## See the top-level LICENSE file for details.
##
## SPDX-License-Identifier: Apache-2.0
################################################################################

import os
import click
import sys

from charge.servers.server_utils import update_mcp_network, get_hostname
from tool_registration import register_tool_server


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
@click.option(
    "--json-file",
    type=str,
    default="known_molecules.json",
    help="Path to known molecules JSON",
)
@click.pass_context
def main(
    ctx,
    port,
    host,
    name,
    copilot_port,
    copilot_host,
    api_key,
    base_url,
    model,
    backend,
    json_file,
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

    # This should be convierted to a shared database instance
    # to ensure the MCP and the backend use the same known molecules

    if not os.path.exists(json_file):
        abs_path = os.path.abspath(json_file)
        with open(abs_path, "w") as f:
            pass

    LMO_MCP.JSON_FILE_PATH = json_file

    mcp = LMO_MCP.mcp

    update_mcp_network(mcp, host, port)

    mcp.run(transport="sse")


if __name__ == "__main__":
    main()
