import sys
import os

cur_dir = os.path.dirname(os.path.abspath(__file__))
sys.path.append(os.path.join(cur_dir, "ChARGe", "experiments", "Molecule_Generation"))

import ChARGe.experiments.Molecule_Generation.mol_server as LMO_MCP

mcp = LMO_MCP.mcp

if __name__ == "__main__":
    mcp.run(transport="sse")
