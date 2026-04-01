from flask_tools.chemistry.smiles_utils import (
    canonicalize_smiles,
    verify_smiles,
)
from charge_backend.moleculedb.purchasable import is_purchasable
from charge_backend.retrosynthesis.database import query_reaction_database
from lc_conductor import (
    BuiltinToolDefinition,
    doc_summary,
)


_BUILTIN_TOOL_DEFINITIONS = [
    BuiltinToolDefinition(
        identifier="verify_smiles",
        function=verify_smiles,
        label="Verify SMILES",
        description=doc_summary(verify_smiles),
    ),
    BuiltinToolDefinition(
        identifier="canonicalize_smiles",
        function=canonicalize_smiles,
        label="Canonicalize SMILES",
        description=doc_summary(canonicalize_smiles),
    ),
    BuiltinToolDefinition(
        identifier="is_purchasable",
        function=is_purchasable,
        label="Check Purchasability",
        description=doc_summary(is_purchasable),
    ),
    BuiltinToolDefinition(
        identifier="query_reaction_database",
        function=query_reaction_database,
        label="Query Reaction Database",
        description=doc_summary(query_reaction_database),
    ),
]


def list_builtin_tool_definitions() -> list[BuiltinToolDefinition]:
    return list(_BUILTIN_TOOL_DEFINITIONS)
