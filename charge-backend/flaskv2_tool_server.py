################################################################################
## Copyright 2025 Lawrence Livermore National Security, LLC. and Binghamton University.
## See the top-level LICENSE file for details.
##
## SPDX-License-Identifier: Apache-2.0
################################################################################

import click
import sys
from loguru import logger

import charge.servers.FLASKv2_reactions as flask
from charge.servers.server_utils import update_mcp_network, get_hostname
from tool_registration import register_tool_server

@click.command()
@click.option("--port", type=int, default=8125, help="Port to run the server on")
@click.option("--host", type=str, default=None, help="Host to run the server on")
@click.option("--name", type=str, default="retro_tools", help="Name of the MCP server")
@click.option("--copilot-port", type=int, default=8001, help="Port to the running copilot backend")
@click.option("--copilot-host", type=str, default=None, help="Host to the running copilot backend")
@click.pass_context
def main(ctx, port, host, name, copilot_port, copilot_host):
    if host is None:
        _, host = get_hostname()

    try:
        register_tool_server(port, host, name, copilot_port, copilot_host)
    except:
        logger.info(f"{name} could not connect to server for registration -- requires manual registration")

    sys.argv = [sys.argv[0]] + ctx.args + [f"--port={port}", f"--host={host}"]
    flask.main()


if __name__ == "__main__":
    main()
