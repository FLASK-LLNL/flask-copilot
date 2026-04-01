import asyncio
import json
from dataclasses import asdict, dataclass
from typing import Any
from uuid import uuid4

from agent_framework import FunctionTool
from fastapi import WebSocket


_PENDING_LOCAL_MCP_RESPONSES: dict[WebSocket, dict[str, asyncio.Future]] = {}


@dataclass(frozen=True)
class LocalMCPToolDefinition:
    name: str
    description: str | None = None
    inputSchema: dict[str, Any] | None = None

    @classmethod
    def from_json(cls, payload: dict[str, Any]) -> "LocalMCPToolDefinition":
        return cls(
            name=str(payload["name"]),
            description=payload.get("description"),
            inputSchema=payload.get("inputSchema"),
        )

    def json(self) -> dict[str, Any]:
        return asdict(self)


def resolve_local_mcp_response(websocket: WebSocket, data: dict[str, Any]) -> bool:
    request_id = data.get("requestId")
    if not isinstance(request_id, str):
        return False

    pending_requests = _PENDING_LOCAL_MCP_RESPONSES.get(websocket)
    if not pending_requests or request_id not in pending_requests:
        return False

    future = pending_requests[request_id]
    if not future.done():
        future.set_result(data)
    return True


async def _await_local_mcp_response(
    websocket: WebSocket,
    request_kind: str,
    payload: dict[str, Any],
    timeout: float = 30.0,
) -> dict[str, Any]:
    request_id = str(uuid4())
    pending_requests = _PENDING_LOCAL_MCP_RESPONSES.setdefault(websocket, {})
    if request_id in pending_requests:
        raise ValueError(f"Duplicate local MCP request id: {request_id}")

    future: asyncio.Future = asyncio.Future()
    pending_requests[request_id] = future
    await websocket.send_json(
        {
            "type": "local-mcp-request",
            "requestId": request_id,
            "requestKind": request_kind,
            **payload,
        }
    )

    try:
        response = await asyncio.wait_for(future, timeout=timeout)
    finally:
        pending_requests.pop(request_id, None)
        if not pending_requests:
            _PENDING_LOCAL_MCP_RESPONSES.pop(websocket, None)

    if not response.get("ok", False):
        raise RuntimeError(response.get("error") or "Local MCP request failed")

    result = response.get("result")
    if not isinstance(result, dict):
        raise RuntimeError("Local MCP response was missing a result payload")
    return result


async def list_local_mcp_tools(
    websocket: WebSocket,
    server_urls: list[str],
) -> dict[str, list[LocalMCPToolDefinition]]:
    if not server_urls:
        return {}

    result = await _await_local_mcp_response(
        websocket,
        "list-tools",
        {"servers": server_urls},
    )
    tool_map: dict[str, list[LocalMCPToolDefinition]] = {}
    for server_result in result.get("servers", []):
        server_url = server_result.get("serverUrl")
        if not isinstance(server_url, str):
            continue
        tools = [
            LocalMCPToolDefinition.from_json(tool)
            for tool in server_result.get("tools", [])
            if isinstance(tool, dict) and tool.get("name")
        ]
        tool_map[server_url] = tools
    return tool_map


def _format_local_mcp_call_result(result: dict[str, Any]) -> str:
    if result.get("isError"):
        raise RuntimeError(json.dumps(result, ensure_ascii=True))

    structured = result.get("structuredContent")
    content = result.get("content")

    if structured is not None and content is None:
        return json.dumps(structured, ensure_ascii=True)

    if isinstance(content, list):
        parts: list[str] = []
        for item in content:
            if isinstance(item, dict) and item.get("type") == "text":
                parts.append(str(item.get("text", "")))
            else:
                parts.append(json.dumps(item, ensure_ascii=True))
        if structured is not None:
            parts.append(json.dumps(structured, ensure_ascii=True))
        return "\n".join(part for part in parts if part)

    if structured is not None:
        return json.dumps(structured, ensure_ascii=True)

    return json.dumps(result, ensure_ascii=True)


async def call_local_mcp_tool(
    websocket: WebSocket,
    server_url: str,
    tool_name: str,
    arguments: dict[str, Any],
) -> str:
    result = await _await_local_mcp_response(
        websocket,
        "call-tool",
        {
            "serverUrl": server_url,
            "toolName": tool_name,
            "arguments": arguments,
        },
        timeout=120.0,
    )
    return _format_local_mcp_call_result(result)


def build_local_mcp_function_tools(
    websocket: WebSocket,
    server_tool_map: dict[str, list[LocalMCPToolDefinition]] | None,
) -> list[FunctionTool]:
    if not server_tool_map:
        return []

    tools: list[FunctionTool] = []

    for server_url, tool_definitions in server_tool_map.items():
        for tool_definition in tool_definitions:

            async def _invoke_local_tool(
                _server_url: str = server_url,
                _tool_name: str = tool_definition.name,
                **kwargs,
            ) -> str:
                return await call_local_mcp_tool(
                    websocket,
                    _server_url,
                    _tool_name,
                    kwargs,
                )

            tools.append(
                FunctionTool(
                    name=tool_definition.name,
                    description=tool_definition.description
                    or f"Proxy MCP tool `{tool_definition.name}`.",
                    func=_invoke_local_tool,
                    input_model=tool_definition.inputSchema
                    or {"type": "object", "properties": {}},
                )
            )

    return tools
