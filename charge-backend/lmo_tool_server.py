import sys
import os

import charge.servers.molecule_generation_server as LMO_MCP
from charge.servers.server_utils import update_mcp_network, get_hostname


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description="Run a ChARGe MCP Server")
    parser.add_argument("--port", type=int, default=8124, help="Port to run the server on")
    parser.add_argument(
        "--host", type=str, default=None, help="Host to run the server on"
    )
    args = parser.parse_args()

    port = args.port
    host = args.host

    mcp = LMO_MCP.mcp

    if host is None:
        _, host = get_hostname()

    update_mcp_network(mcp, host, port)

    mcp.run(transport="sse")
