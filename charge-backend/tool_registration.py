################################################################################
## Copyright 2025 Lawrence Livermore National Security, LLC. and Binghamton University.
## See the top-level LICENSE file for details.
##
## SPDX-License-Identifier: Apache-2.0
################################################################################

from dataclasses import dataclass
from fastapi import Request
import socket
from loguru import logger

@dataclass(frozen=True)
class Server:
    address: str
    port: int

SERVERS: set[Server] = set()    

def get_client_info(request: Request):
    """Get client IP and hostname with fallbacks"""
    # Try to get real IP from X-Forwarded-For header
    forwarded_for = request.headers.get("X-Forwarded-For")
    if forwarded_for:
        client_ip = forwarded_for.split(",")[0].strip()
    else:
        # Fallback to direct connection IP
        client_ip = request.client.host
    
    # Try to resolve hostname
    try:
        hostname = socket.gethostbyaddr(client_ip)[0]
    except (socket.herror, socket.gaierror, OSError):
        hostname = client_ip  # Use IP if resolution fails
    
    return hostname

@dataclass
class RegistrationRequest:
    port: int

async def register_post(request: Request, data: RegistrationRequest):
    hostname = get_client_info(request)

    SERVERS.add(Server(
        address=hostname,
        port=data.port
    ))
    return {"status": f"registered MCP server at {hostname}:{data.port}"}
