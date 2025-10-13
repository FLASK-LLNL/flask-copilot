import sys
import os


cur_dir = os.path.dirname(os.path.abspath(__file__))
sys.path.append(os.path.join(cur_dir, "ChARGe", "experiments", "Retrosynthesis"))

import ChARGe.experiments.Retrosynthesis.reaction_server as RETRO_MCP
import charge.servers.server_utils as server_utils
from charge.servers.AizynthTools import is_molecule_synthesizable

port = server_utils.args.port
host = server_utils.args.host

mcp = RETRO_MCP.template_free_mcp

mcp.settings.port = port
mcp.settings.host = host

mcp.tool()(is_molecule_synthesizable)

if __name__ == "__main__":
    mcp.run(transport="sse")
