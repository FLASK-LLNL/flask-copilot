import sys
import os


from charge.servers.AiZynthTools import is_molecule_synthesizable, RetroPlanner

import charge.experiments.Retrosynthesis.reaction_server as RETRO_MCP
from charge.servers.server_utils import update_mcp_network, get_hostname

if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description="Run a ChARGe MCP Server")
    parser.add_argument("--port", type=int, default=8123, help="Port to run the server on")
    parser.add_argument(
        "--host", type=str, default=None, help="Host to run the server on"
    )
    parser.add_argument(
        "--config",
        type=str,
        default="config.yml",
        help="Path to the configuration file for AiZynthFinder",
    )
    args = parser.parse_args()

    port = args.port
    host = args.host
    if host is None:
        _, host = get_hostname()

    mcp = RETRO_MCP.template_free_mcp

    RetroPlanner.initialize(configfile=args.config)
    mcp.tool()(is_molecule_synthesizable)

    if host is None:
        _, host = get_hostname()

    update_mcp_network(mcp, host, port)

    mcp.run(transport="sse")
