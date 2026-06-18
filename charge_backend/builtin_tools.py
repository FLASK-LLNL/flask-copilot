from flask_tools.chemistry.smiles_utils import (
    canonicalize_smiles,
    verify_smiles,
)
from charge_backend.moleculedb.purchasable import is_purchasable
from charge_backend.retrosynthesis.database import query_reaction_database
from charge_backend.pdf import PdfDocumentRegistry
from lc_conductor import (
    BuiltinToolDefinition,
    doc_summary,
)


def _document_consult_tool(registry: PdfDocumentRegistry):
    async def consult_with_document(
        question: str,
        document_id: str | None = None,
        max_tool_calls: int = 14,
    ) -> str:
        """Consult the active uploaded PDF reference with a private subagent."""
        return await registry.consult(question, document_id, max_tool_calls)

    return consult_with_document


def list_builtin_tool_definitions(
    pdf_registry: PdfDocumentRegistry | None = None,
) -> list[BuiltinToolDefinition]:
    definitions = [
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
    if pdf_registry is not None:
        consult_with_document = _document_consult_tool(pdf_registry)
        definitions.append(
            BuiltinToolDefinition(
                identifier="consult_with_document",
                function=consult_with_document,
                label="Consult Document",
                description=doc_summary(consult_with_document),
            )
        )
    return definitions
