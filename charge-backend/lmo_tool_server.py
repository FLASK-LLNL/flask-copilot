import sys
import os

cur_dir = os.path.dirname(os.path.abspath(__file__))
sys.path.append(os.path.join(cur_dir, "ChARGe", "experiments", "Molecule_Generation"))

import ChARGe.experiments.Molecule_Generation.mol_server as LMO_MCP

import argparse

parser = argparse.ArgumentParser(description="Run a ChARGe MCP Server")
parser.add_argument("--port", type=int, default=8000, help="Port to run the server on")
parser.add_argument(
    "--host", type=str, default="127.0.0.1", help="Host to run the server on"
)
args = parser.parse_args()


port = args.port
host = args.host

mcp = LMO_MCP.mcp

mcp.settings.port = port
mcp.settings.host = host

if __name__ == "__main__":
    mcp.run(transport="sse")
