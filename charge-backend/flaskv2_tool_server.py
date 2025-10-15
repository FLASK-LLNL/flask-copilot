################################################################################
## Copyright 2025 Lawrence Livermore National Security, LLC. and Binghamton University.
## See the top-level LICENSE file for details.
##
## SPDX-License-Identifier: Apache-2.0
################################################################################


import charge.servers.FLASKv2_reactions as flask
from charge.servers.server_utils import update_mcp_network, get_hostname


if __name__ == "__main__":
    flask.main()
