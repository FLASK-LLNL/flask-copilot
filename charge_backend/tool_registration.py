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
from typing import Optional, Tuple
import time
import os
import re

from charge.utils.system_utils import check_server_paths
from autogen_ext.tools.mcp import McpWorkbench, SseServerParams
from charge.clients.autogen_utils import (
    _list_wb_tools,
)


def split_url(url: str) -> Tuple[str, int, str, str]:
    # Regular expression pattern
    pattern = r"^(https?://)?([^:/]+)(?::(\d+))?(?:/(.+?))?/?$"

    match = re.match(pattern, url)

    if match:
        protocol = match.group(1) or ""
        host = match.group(2)
        port = match.group(3) or ""
        path = match.group(4) or ""

    if not port and not path:
        raise ValueError(
            f"Unusable URL provide {url} -- requires either a port or a path"
        )

    return host, port, path, protocol


@dataclass
class ToolList:
    server: str
    names: Optional[str] = None

    def json(self):
        return asdict(self)


class ToolServer(BaseModel):
    address: str
    port: int
    path: str
    name: str
    protocol: Optional[str] = None

    def __str__(self):
        path_if_valid = f"/{self.path}" if self.path else ""
        protocol_if_valid = f"{self.protocol}" if self.protocol else "http://"
        return f"{protocol_if_valid}{self.address}:{self.port}{path_if_valid}/sse"

    def long_name(self):
        path_if_valid = f"/{self.path}" if self.path else ""
        protocol_if_valid = f"{self.protocol}" if self.protocol else "http://"
        return f"[{self.name}] {protocol_if_valid}{self.address}:{self.port}{path_if_valid}/sse"


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
        # Check if file exists
        if not os.path.exists(filename):
            logger.info(f"Server list file does not exist: {filename}")
            return

        # Check if file is readable
        if not os.access(filename, os.R_OK):
            logger.error(f"Server list file is not readable: {filename}")
            return

        # Check if file has non-zero size
        try:
            file_size = os.path.getsize(filename)
            if file_size == 0:
                logger.info(f"Server list file is empty (0 bytes): {filename}")
                return
        except OSError as e:
            logger.error(f"Error getting file size for {filename}: {e}")
            return

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


def register_url(
    filename: str,
    hostname: str,
    port: int,
    path: Optional[str] = "",
    protocol: Optional[str] = "",
    name: Optional[str] = "",
):
    path_if_valid = f"/{path}" if path else ""
    protocol_if_valid = f"{protocol}" if protocol else "http://"
    key = f"{protocol_if_valid}{hostname}:{port}{path_if_valid}"
    new_server = ToolServer(
        address=hostname, port=port, path=path, name=name, protocol=protocol
    )

    old_server = SERVERS.servers.pop(key, None)
    if old_server:
        logger.info(
            f"Replacing server at {key} with new registration: {old_server.long_name()} -> {new_server.long_name()}"
        )

    SERVERS.servers[key] = new_server
    if filename:
        # Check if file exists
        file_exists = os.path.exists(filename)

        msg_base = f"registered MCP server {name} at {key}"
        # Check if file is writable (or parent directory is writable for new files)
        if file_exists:
            if not os.access(filename, os.W_OK):
                logger.error(f"Server list file is not writable: {filename}")
                return {
                    "status": f"{msg_base} (warning: could not save to disk - file not writable)"
                }
        else:
            # For new files, check if parent directory is writable
            parent_dir = os.path.dirname(filename) or "."
            if not os.access(parent_dir, os.W_OK):
                logger.error(
                    f"Cannot create server list file - parent directory not writable: {parent_dir}"
                )
                return {
                    "status": f"{msg_base} (warning: could not save to disk - directory not writable)"
                }

        try:
            with open(filename, "w") as f:
                f.write(SERVERS.model_dump_json(indent=4))
        except PermissionError as e:
            logger.error(
                f"Permission denied writing to server list file {filename}: {e}"
            )
            return {
                "status": f"{msg_base} (warning: could not save to disk - permission denied)"
            }
        except OSError as e:
            logger.error(f"OS error writing to server list file {filename}: {e}")
            return {"status": f"{msg_base} (warning: could not save to disk - {e})"}
        except Exception as e:
            logger.error(f"Error writing to server list file {filename}: {e}")
            return {"status": f"{msg_base} (warning: could not save to disk)"}

    return {"status": f"{msg_base}"}


async def register_post(filename: str, request: Request, data: RegistrationRequest):
    hostname = data.host
    if not hostname:
        hostname = get_client_info(request)
    return register_url(filename, hostname, data.port, data.name)


def register_tool_server(port, host, name, copilot_port, copilot_host):
    for i in range(5):
        try:
            try:
                url = f"https://{copilot_host}:{copilot_port}/register"
                response = requests.post(
                    url, json={"host": host, "port": port, "name": name}
                )
            except:
                url = f"http://{copilot_host}:{copilot_port}/register"
                response = requests.post(
                    url, json={"host": host, "port": port, "name": name}
                )
            logger.info(response.json())
            break
        except:
            if i == 4:
                logger.error("Could not connect to server for registration! Exiting")
                raise
            logger.info(
                "Could not connect to server for registration, retrying in 10 seconds"
            )
            time.sleep(10)
            continue


def list_server_urls() -> list[str]:
    server_urls = []
    invalid_keys = []
    for key, server in SERVERS.servers.items():
        validated_server = check_server_paths(f"{server}")
        if validated_server:
            server_urls.append(f"{server}")
        else:
            logger.info(
                f"Previously cached URL is no longer valid - removing {server.long_name()} from cache"
            )
            invalid_keys.append(key)

    for key in invalid_keys:
        SERVERS.servers.pop(key)

    assert server_urls is not None, "Server URLs must be registered"
    for url in server_urls:
        assert url.endswith("/sse"), f"Server URL {url} must end with /sse"

    return server_urls


async def list_server_tools(urls: list[str]):
    workbenches = [McpWorkbench(SseServerParams(url=server)) for server in urls]
    return await _list_wb_tools(workbenches)
