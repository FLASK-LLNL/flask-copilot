################################################################################
## Copyright 2025 Lawrence Livermore National Security, LLC. and Binghamton University.
## See the top-level LICENSE file for details.
##
## SPDX-License-Identifier: Apache-2.0
################################################################################

from dataclasses import dataclass, asdict
from fastapi import Request
import socket
from loguru import logger
import requests
from pydantic import BaseModel
import json
from typing import Optional

from charge.utils.system_utils import check_server_paths
from autogen_ext.tools.mcp import McpWorkbench, SseServerParams
from charge.clients.autogen_utils import (
    _list_wb_tools,
)

@dataclass
class ToolList:
    server: str
    names: Optional[str] = None

    def json(self):
        return asdict(self)

class ToolServer(BaseModel):
    address: str
    port: int
    name: str

    def __str__(self):
        return f"http://{self.address}:{self.port}/sse"
    def long_name(self):
        return f"[{self.name}] http://{self.address}:{self.port}/sse"

class ToolServerDict(BaseModel):
    servers: dict[str, ToolServer]

SERVERS: ToolServerDict = ToolServerDict(servers={})

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

def reload_server_list(filename: str):
    if filename:
        try:
            with open(filename, "r") as f:
                data: ToolServerDict = ToolServerDict.model_validate_json(f.read())
                SERVERS.servers = data.servers
        except FileNotFoundError as e:
            logger.info(e)
            return
        except json.JSONDecodeError as e:
            logger.info(e)
            return
        except Exception as e:
            logger.info(e)
            return
    else:
        return

async def register_post(filename: str, request: Request, data: RegistrationRequest):
    hostname = data.host
    if not hostname:
        hostname = get_client_info(request)

    key = f"{hostname}:{data.port}"
    new_server = ToolServer(
        address=hostname,
        port=data.port,
        name=data.name
    )

    old_server = SERVERS.servers.pop(key, None)
    if old_server:
        logger.info(f"Replacing server at {key} with new registration: {old_server.long_name()} -> {new_server.long_name()}")

    SERVERS.servers[key] = new_server
    if filename:
        try:
            with open(filename, "w") as f:
                f.write(SERVERS.model_dump_json(indent=4))
        except Exception as e:
            logger.info(e)
            pass

    return {"status": f"registered MCP server {data.name} at {hostname}:{data.port}"}

def register_tool_server(port, host, name, copilot_port, copilot_host):
    try:
        url = f"https://{copilot_host}:{copilot_port}/register"
        response = requests.post(url, json={"host": host, "port": port, "name": name})
    except:
        url = f"http://{copilot_host}:{copilot_port}/register"
        response = requests.post(url, json={"host": host, "port": port, "name": name})

    logger.info(response.json())

def list_server_urls() -> list[str]:
    server_urls = []
    invalid_keys = []
    for key, server in SERVERS.servers.items():
        validated_server = check_server_paths(f"{server}")
        if validated_server:
            server_urls.append(f"{server}")
        else:
            logger.info(f"Previously cached URL is no longer valid - removing {server.long_name()} from cache")
            invalid_keys.append(key)

    for key in invalid_keys:
        SERVERS.pop(key)

    assert server_urls is not None, "Server URLs must be registered"
    for url in server_urls:
        assert url.endswith("/sse"), f"Server URL {url} must end with /sse"

    return server_urls

async def list_server_tools(urls: list[str]):
    workbenches = [McpWorkbench(SseServerParams(url=server)) for server in urls]
    return await _list_wb_tools(workbenches)
