################################################################################
## Copyright 2025 Lawrence Livermore National Security, LLC. and Binghamton University.
## See the top-level LICENSE file for details.
##
## SPDX-License-Identifier: Apache-2.0
################################################################################

import click
from loguru import logger
import uvicorn

from charge.servers.AiZynthTools import is_molecule_synthesizable, RetroPlanner

import charge.servers.retrosynthesis_reaction_server as RETRO_MCP
from charge.servers.server_utils import update_mcp_network, get_hostname
from tool_registration import register_tool_server, get_asgi_app


@click.command()
@click.option("--port", type=int, default=8123, help="Port to run the server on")
@click.option("--host", type=str, default=None, help="Host to run the server on")
@click.option("--name", type=str, default="retro_tools", help="Name of the MCP server")
@click.option(
    "--copilot-port", type=int, default=8001, help="Port to the running copilot backend"
)
@click.option(
    "--copilot-host", type=str, default=None, help="Host to the running copilot backend"
)
@click.option(
    "--config",
    type=str,
    default="config.yml",
    help="Path to the configuration file for AiZynthFinder",
)
def main(port, host, name, copilot_port, copilot_host, config):
    if host is None:
        _, host = get_hostname()

    try:
        register_tool_server(port, host, name, copilot_port, copilot_host)
    except:
        logger.info(
            f"{name} could not connect to server for registration -- requires manual registration"
        )

    mcp = RETRO_MCP.template_free_mcp

    RetroPlanner.initialize(configfile=config)
    mcp.tool()(is_molecule_synthesizable)

    asgi_app = get_asgi_app(mcp)
    if asgi_app:
        uvicorn.run(asgi_app, host=host or "0.0.0.0", port=port)
    else:
        logger.error("Could not access FastMCP ASGI app")


if __name__ == "__main__":
    main()
