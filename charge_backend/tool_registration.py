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
from typing import Optional, Tuple, List, Dict
import time
import os
import re
import asyncio
from charge.utils.mcp_workbench_utils import (
    _setup_mcp_workbenches,
    _close_mcp_workbenches,
)

from charge.utils.system_utils import check_server_paths
from autogen_ext.tools.mcp import McpWorkbench, SseServerParams
from charge.clients.autogen_utils import (
    _list_wb_tools,
)


class ValidateMCPServerRequest(BaseModel):
    url: str
    name: Optional[str] = None


class DeleteMCPServerRequest(BaseModel):
    url: str


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
        # Ensure path ends with /sse
        sse_path_if_valid = (
            path_if_valid
            if path_if_valid.endswith("/sse")
            else f"{path_if_valid.rstrip('/')}/sse"
        )
        return f"{protocol_if_valid}{self.address}:{self.port}{sse_path_if_valid}"

    def long_name(self):
        short_name = self.__str__()
        return f"[{self.name}] {short_name}"


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


async def _check_mcp_connectivity(url: str, timeout: float) -> List[Dict]:
    """
    Connect to an MCP server and retrieve its tools list using existing workbench utilities.

    Args:
        url: MCP server URL (should end with /sse)
        timeout: Connection timeout in seconds

    Returns:
        List of tools with name and description

    Raises:
        Exception: If connection fails or server is unreachable
    """
    from charge.utils.system_utils import check_url_exists

    # Ensure URL ends with /sse
    sse_url = url if url.endswith("/sse") else f"{url.rstrip('/')}/sse"

    # First do a quick check if the URL is reachable
    if not check_url_exists(sse_url):
        raise Exception(f"Server at {sse_url} is not reachable")

    # Now use workbench utilities to connect and get tools
    try:
        workbenches = await _setup_mcp_workbenches(paths=[], urls=[sse_url])

        if not workbenches:
            raise Exception("Failed to create workbench for server")

        # Get tools from the workbench
        tools = []
        for workbench in workbenches:
            try:
                workbench_tools = await workbench.list_tools()
                for tool in workbench_tools:
                    tools.append(
                        {
                            "name": tool.get("name"),
                            "description": tool.get("description"),
                        }
                    )
            except Exception as e:
                logger.warning(f"Failed to list tools from workbench: {e}")

        # Clean up workbenches
        await _close_mcp_workbenches(workbenches)

        return tools

    except asyncio.TimeoutError:
        raise TimeoutError("Connection timeout")
    except Exception as e:
        raise Exception(f"Validation error: {str(e)}")


async def validate_and_register_mcp_server(
    filename: str, url: str, name: Optional[str] = None, timeout: float = 10.0
) -> Dict:
    """
    Validate an MCP server URL by connecting to it and listing tools.
    If valid, register it.

    Args:
        filename: Path to server cache file
        url: MCP server URL
        name: Optional display name for the server
        timeout: Connection timeout in seconds

    Returns:
        Dict with status, tools (if successful), and error message (if failed)
    """
    # Parse the URL
    try:
        host, port, path, protocol = split_url(url)
    except ValueError as e:
        return {"status": "error", "error": str(e)}

    # Ensure URL ends with /sse
    sse_url = url if url.endswith("/sse") else f"{url.rstrip('/')}/sse"

    # Validate connectivity using existing utilities
    try:
        tools = await _check_mcp_connectivity(sse_url, timeout)

        # If validation successful, register the server
        if not name:
            name = f"{host}:{port}"

        registration_result = register_url(
            filename, host, int(port) if port else 80, path, protocol, name
        )

        return {
            "status": "connected",
            "tools": tools,
            "url": sse_url,
            "registration": registration_result,
        }

    except Exception as e:
        logger.error(f"Failed to validate MCP server at {sse_url}: {e}")
        return {"status": "disconnected", "error": str(e), "url": sse_url}


async def check_registered_servers(filename: str) -> Dict[str, Dict]:
    """
    Check connectivity of all registered servers using existing utilities.

    Args:
        filename: Path to server cache file (unused but kept for API consistency)

    Returns:
        Dict mapping server URLs to their status/tools
    """
    results = {}

    for key, server in SERVERS.servers.items():
        url = str(server)
        try:
            tools = await _check_mcp_connectivity(url, timeout=5.0)
            results[url] = {"status": "connected", "tools": tools}
        except Exception as e:
            results[url] = {"status": "disconnected", "error": str(e)}

    return results


def delete_registered_server(filename: str, url: str) -> Dict:
    """Delete a registered MCP server."""
    # Find the server by URL
    key_to_delete = None
    for key, server in SERVERS.servers.items():
        if str(server) == url or key == url:
            key_to_delete = key
            break

    if not key_to_delete:
        logger.warning(f"Server not found for deletion: {url}")
        return {
            "status": "not_found",
            "message": f"Server {url} not found in registered servers",
        }

    # Delete from in-memory dict
    deleted_server = SERVERS.servers.pop(key_to_delete)
    logger.info(f"Deleted server: {deleted_server.long_name()}")

    # Save to file
    if filename:
        try:
            with open(filename, "w") as f:
                f.write(SERVERS.model_dump_json(indent=4))
            logger.info(f"Saved updated server list to {filename}")
        except Exception as e:
            logger.error(f"Failed to save server list: {e}")
            return {
                "status": "error",
                "message": f"Deleted from memory but failed to save: {str(e)}",
            }

    return {
        "status": "deleted",
        "message": f"Successfully deleted server: {deleted_server.long_name()}",
    }


async def get_registered_servers(filename: str) -> Dict:
    """
    Get list of all registered MCP servers and their status.

    This endpoint aggregates server info and checks connectivity
    using existing validation utilities.
    """
    # Get connectivity status for all servers
    statuses = await check_registered_servers(filename)

    # Build response with server info and status
    servers = []
    for key, server in SERVERS.servers.items():
        url = str(server)
        status_info = statuses.get(url, {"status": "unknown"})

        servers.append(
            {
                "id": key,
                "url": url,
                "name": server.name,
                "address": server.address,
                "port": server.port,
                "path": server.path,
                **status_info,
            }
        )

    return {"servers": servers}


async def validate_mcp_server_endpoint(
    filename: str, request: Request, data: ValidateMCPServerRequest
):
    """
    Validate and register an MCP server URL.
    Called by the frontend when a user adds a new tool server.

    This endpoint uses the existing system_utils.check_url_exists() and
    mcp_workbench_utils for validation.
    """

    # Get client info for logging
    client_info = get_client_info(request)
    logger.info(f"validate request from {client_info} for MCP server: {data.url}")

    result = await validate_and_register_mcp_server(filename, data.url, data.name)

    logger.info(f"Validate result: {result}")

    return result


async def delete_mcp_server_endpoint(
    filename: str, request: Request, data: DeleteMCPServerRequest
):
    """Delete a registered MCP server."""
    from tool_registration import delete_registered_server

    client_info = get_client_info(request)
    logger.info(f"Delete request from {client_info} for MCP server: {data.url}")

    result = delete_registered_server(filename, data.url)

    logger.info(f"Delete result: {result}")

    return result


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
