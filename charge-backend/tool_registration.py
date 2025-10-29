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
import requests

@dataclass(frozen=True)
class Server:
    address: str
    port: int
    name: str

    def __str__(self):
        return f"http://{self.address}:{self.port}/sse"

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
    host: str
    port: int
    name: str

async def register_post(request: Request, data: RegistrationRequest):
    hostname = data.host
    if not hostname:
        hostname = get_client_info(request)

    SERVERS.add(Server(
        address=hostname,
        port=data.port,
        name=data.name
    ))
    return {"status": f"registered MCP server {data.name} at {hostname}:{data.port}"}

def register_tool_server(port, host, name, copilot_port, copilot_host):
    url = f"http://{copilot_host}:{copilot_port}/register"
    response = requests.post(url, json={"host": host, "port": port, "name": name})
    logger.info(response.json())

def list_server_urls() -> list[str]:
    server_urls = []
    for server in SERVERS:
        server_urls.append(f"{server}")

    assert server_urls is not None, "Server URLs must be registered"
    for url in server_urls:
        assert url.endswith("/sse"), f"Server URL {url} must end with /sse"

    return server_urls
