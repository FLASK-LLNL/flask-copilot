from __future__ import annotations

import asyncio
import base64
import binascii
import hashlib
import json
import shutil
import tempfile
from dataclasses import asdict, dataclass, replace
from pathlib import Path
from typing import Any, Optional

from agent_framework import Content
from charge_backend.flask_experiment import FlaskExperiment
from charge.clients.agent import AgentBackend
from charge.tasks.task import Task
from loguru import logger

from charge_backend.pdf.pdf_document import (
    SubagentDocumentPdf,
    format_search_results,
    page_read_to_json,
)


MAX_PDF_ATTACHMENT_BYTES = 100 * 1024 * 1024
SUBAGENT_MAX_TOOL_CALLS = 14


@dataclass(frozen=True)
class PdfReferenceMetadata:
    id: str
    documentId: str
    name: str
    mimeType: str
    sizeBytes: int
    title: str
    pageCount: int
    createdAt: str = ""
    status: str = "available"

    def json(self) -> dict[str, Any]:
        return asdict(self)


class PdfDocument:
    def __init__(
        self,
        metadata: PdfReferenceMetadata,
        path: Path,
        document: Any,
        agent_backend: Optional[AgentBackend],
    ):
        self.metadata = metadata
        self.path = path
        self._document = document
        self._agent_backend = agent_backend

    async def consult(
        self, question: str, max_tool_calls: int = SUBAGENT_MAX_TOOL_CALLS
    ) -> str:
        if not question.strip():
            return "Please provide a non-empty question for the document."

        logger.info(
            "Document consult started for PDF '{}' ({}, {} pages): {}",
            self.metadata.title,
            self.metadata.name,
            self.metadata.pageCount,
            question.strip(),
        )
        overview = await asyncio.to_thread(self._document.overview)
        if overview.title == self.path.stem:
            overview = replace(overview, title=self.metadata.title)
        task = Task(
            system_prompt=_instructions(overview.prompt_context()),
            user_prompt=(
                "Question to answer from the uploaded PDF reference:\n"
                f"{question.strip()}\n\n"
                "Use search and read_page as needed. Cite PDF page numbers."
            ),
            builtin_tools=[self._search_tool(), self._read_page_tool()],
        )

        # Deliberately use a fresh Experiment so the subagent has no
        # orchestrator memory or task history, but reuse the session's agent
        # backend so the subagent runs on this user's configured backend.
        experiment = FlaskExperiment(task=None, backend=self._agent_backend)
        agent = experiment.create_agent_with_experiment_state(
            task=task,
            agent_key="Scholar",
            max_retries=1,
            max_tool_calls=max(1, min(int(max_tool_calls), 40)),
        )
        result = await agent.run()
        logger.info(
            "Document consult completed for PDF '{}' ({})",
            self.metadata.title,
            self.metadata.name,
        )
        return str(result)

    def _search_tool(self):
        def search(
            query: str,
            num_exact_results: int = 10,
            num_fuzzy_results: int = 3,
        ) -> str:
            """Search the PDF text layer for exact and fuzzy page-numbered matches."""
            exact_limit = max(0, min(int(num_exact_results), 25))
            fuzzy_limit = max(0, min(int(num_fuzzy_results), 15))
            logger.info(
                "PDF search '{}' query={!r} exact_limit={} fuzzy_limit={}",
                self.metadata.title,
                query,
                exact_limit,
                fuzzy_limit,
            )
            exact, fuzzy = self._document.search(
                query,
                exact_limit,
                fuzzy_limit,
            )
            logger.info(
                "PDF search '{}' query={!r} returned exact_pages={} fuzzy_pages={}",
                self.metadata.title,
                query,
                [hit.page_number for hit in exact],
                [hit.page_number for hit in fuzzy],
            )
            return format_search_results(exact, fuzzy)

        return search

    def _read_page_tool(self):
        def read_page(page_number: int, force_image: bool = False) -> list[Content]:
            """Read a PDF page; include a rendered PNG when visual/layout evidence is needed."""
            requested_page = int(page_number)
            requested_image = bool(force_image)
            logger.info(
                "PDF read_page '{}' page={} force_image={}",
                self.metadata.title,
                requested_page,
                requested_image,
            )
            page = self._document.read_page(requested_page, requested_image)
            payload = page_read_to_json(page)
            logger.info(
                "PDF read_page '{}' page={} mode={} text_chars={} image_reason={}",
                self.metadata.title,
                page.page_number,
                "image+text" if page.image_data_url else "text",
                len(page.text or ""),
                page.image_reason,
            )
            contents: list[Content] = [
                Content.from_text(json.dumps(payload, ensure_ascii=False))
            ]
            if page.image_data_url:
                contents.append(
                    Content.from_text(
                        f"Rendered image for PDF page {page.page_number}. "
                        "Use OCR/vision if the extracted text is incomplete or layout matters."
                    )
                )
                contents.append(
                    Content.from_uri(page.image_data_url, media_type="image/png")
                )
            return contents

        return read_page


class PdfDocumentRegistry:
    """
    A document registry (metadata, summary, ToC) that caches previously uploaded
    PDFs.
    """

    def __init__(
        self,
        agent_backend: Optional[AgentBackend],
        cache_dir: str | Path | None = None,
    ):
        self._agent_backend = agent_backend
        self._base_dir = Path(cache_dir) if cache_dir else Path(tempfile.mkdtemp())
        self._base_dir.mkdir(parents=True, exist_ok=True)
        self._document: PdfDocument | None = None

    def active_metadata(self) -> PdfReferenceMetadata | None:
        if self._document is None:
            return None
        return self._document.metadata

    def has_active_document(self) -> bool:
        return self._document is not None

    def set_from_attachment(self, attachment: dict[str, Any]) -> PdfReferenceMetadata:
        normalized = validate_pdf_reference(attachment)
        raw_bytes = _decode_data_url(normalized["dataUrl"])
        document_id = hashlib.sha256(raw_bytes).hexdigest()
        path = self._base_dir / f"{document_id}.pdf"
        path.write_bytes(raw_bytes)

        document = SubagentDocumentPdf(path, cache_dir=self._base_dir / "pdf_cache")
        overview = document.overview()
        title = (
            overview.title
            if overview.title != path.stem
            else Path(normalized["name"]).stem
        )
        metadata = PdfReferenceMetadata(
            id=normalized["id"],
            documentId=document_id,
            name=normalized["name"],
            mimeType=normalized["mimeType"],
            sizeBytes=normalized["sizeBytes"],
            title=title,
            pageCount=overview.page_count,
            createdAt=normalized["createdAt"],
        )
        self._document = PdfDocument(metadata, path, document, self._agent_backend)
        return metadata

    def clear(self) -> None:
        self._document = None

    def cleanup(self) -> None:
        self._document = None
        shutil.rmtree(self._base_dir, ignore_errors=True)

    async def consult(
        self,
        question: str,
        document_id: str | None = None,
        max_tool_calls: int = SUBAGENT_MAX_TOOL_CALLS,
    ) -> str:
        if not self.has_active_document():
            return (
                "No PDF reference is available in this session. Ask the user to "
                "reupload the PDF in Customize > References before consulting it."
            )
        assert self._document is not None
        active_id = self._document.metadata.documentId
        if document_id and document_id != active_id:
            return (
                f"The requested PDF `{document_id}` is not active. The only active "
                f"PDF is `{active_id}` ({self._document.metadata.name})."
            )

        return await self._document.consult(question, max_tool_calls=max_tool_calls)


def validate_pdf_reference(attachment: dict[str, Any]) -> dict[str, Any]:
    if not isinstance(attachment, dict):
        raise ValueError("PDF reference must be an object")

    data_url = attachment.get("dataUrl")
    mime_type = attachment.get("mimeType")
    if not isinstance(data_url, str) or not data_url.startswith("data:"):
        raise ValueError("PDF reference must be a data URL")
    if mime_type != "application/pdf":
        raise ValueError("PDF reference must be an application/pdf file")

    header, separator, encoded = data_url.partition(",")
    if separator != "," or ";base64" not in header:
        raise ValueError("PDF reference must be base64 encoded")

    header_mime_type = header.removeprefix("data:").split(";", 1)[0]
    if header_mime_type != mime_type:
        raise ValueError("PDF reference has inconsistent MIME metadata")

    decoded_size = _decoded_size(encoded)
    if decoded_size > MAX_PDF_ATTACHMENT_BYTES:
        raise ValueError("PDF reference exceeds 100 MB")

    return {
        "id": str(attachment.get("id") or "pdf_reference"),
        "name": str(attachment.get("name") or "Reference PDF"),
        "mimeType": mime_type,
        "sizeBytes": decoded_size,
        "dataUrl": data_url,
        "createdAt": str(attachment.get("createdAt") or ""),
    }


def _decode_data_url(data_url: str) -> bytes:
    _, _, encoded = data_url.partition(",")
    try:
        return base64.b64decode(encoded, validate=True)
    except (binascii.Error, ValueError) as exc:
        raise ValueError("PDF reference is not valid base64") from exc


def _decoded_size(encoded: str) -> int:
    try:
        return len(base64.b64decode(encoded, validate=True))
    except (binascii.Error, ValueError) as exc:
        raise ValueError("PDF reference is not valid base64") from exc


def _instructions(document_context: str) -> str:
    return f"""You are a careful, document-grounded research assistant.

You answer questions by consulting one PDF book or article. You have access to the document's table of contents,
a detected back-of-book index when present, and two tools: search and read_page.

Core rules:
- Ground the answer in this document. Do not rely on outside knowledge except for ordinary language understanding.
- Use search first for names, phrases, and concepts unless the table of contents clearly identifies the needed pages.
- Use read_page to inspect important pages before making claims, especially when a search hit is only a clue.
- When any visual aid, table layout, figure, chart, diagram, map, photo, illustration, formula layout, or page design needs inspection, call read_page with force_image=true.
- Cite PDF page numbers for substantive claims using forms like [p. 42] or [pp. 42-44].
- If the PDF page image is provided after read_page, read the image directly and treat it as the source for that page.
- If the document does not answer the question, say that and summarize the closest relevant evidence.
- Distinguish the author's claims from your interpretation.
- Prefer a detailed, structured answer with concise quotations only when the wording matters.
- Never fabricate page references. If evidence is uncertain, mark it as uncertain.

Document context:
{document_context}
"""
