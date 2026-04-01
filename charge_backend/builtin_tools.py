import inspect
from dataclasses import dataclass
from typing import Any, Callable, Iterable

from flask_tools.chemistry.smiles_utils import (
    canonicalize_smiles,
    verify_smiles,
)
from charge_backend.moleculedb.purchasable import is_purchasable
from charge_backend.retrosynthesis.database import query_reaction_database


@dataclass(frozen=True)
class BuiltinToolDefinition:
    identifier: str
    function: Callable[..., Any]
    label: str
    description: str

    def to_client_tool(self) -> dict[str, Any]:
        return {
            "kind": "builtin",
            "identifier": self.identifier,
            "server": self.label,
            "names": [self.function.__name__],
            "description": self.description,
            "executionScope": "backend",
        }


def _doc_summary(func: Callable[..., Any]) -> str:
    doc = inspect.getdoc(func)
    if not doc:
        return f"Run the backend function `{func.__name__}`."
    return doc.splitlines()[0].strip()


_BUILTIN_TOOL_DEFINITIONS = [
    BuiltinToolDefinition(
        identifier="verify_smiles",
        function=verify_smiles,
        label="Verify SMILES",
        description=_doc_summary(verify_smiles),
    ),
    BuiltinToolDefinition(
        identifier="canonicalize_smiles",
        function=canonicalize_smiles,
        label="Canonicalize SMILES",
        description=_doc_summary(canonicalize_smiles),
    ),
    BuiltinToolDefinition(
        identifier="is_purchasable",
        function=is_purchasable,
        label="Check Purchasability",
        description=_doc_summary(is_purchasable),
    ),
    BuiltinToolDefinition(
        identifier="query_reaction_database",
        function=query_reaction_database,
        label="Query Reaction Database",
        description=_doc_summary(query_reaction_database),
    ),
]


def list_builtin_tool_definitions() -> list[BuiltinToolDefinition]:
    return list(_BUILTIN_TOOL_DEFINITIONS)


def resolve_builtin_tools(
    identifiers: Iterable[str] | None,
    definitions: Iterable[BuiltinToolDefinition] | None = None,
) -> list[Callable[..., Any]]:
    tool_definitions = list(definitions or _BUILTIN_TOOL_DEFINITIONS)
    if identifiers is None:
        return [tool.function for tool in tool_definitions]

    tool_map = {tool.identifier: tool for tool in tool_definitions}
    resolved_tools: list[Callable[..., Any]] = []
    for identifier in identifiers:
        tool = tool_map.get(identifier)
        if tool is not None:
            resolved_tools.append(tool.function)
    return resolved_tools
