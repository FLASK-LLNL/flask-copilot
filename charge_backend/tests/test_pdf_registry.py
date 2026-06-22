import base64
import asyncio
import sys
from pathlib import Path

import pytest
import pymupdf

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from pdf.registry import PdfDocumentRegistry, validate_pdf_reference


def pdf_attachment(payload: bytes = b"%PDF-1.4\n", **overrides):
    encoded = base64.b64encode(payload).decode("ascii")
    data = {
        "id": "pdf-1",
        "name": "reference.pdf",
        "mimeType": "application/pdf",
        "sizeBytes": len(payload),
        "dataUrl": f"data:application/pdf;base64,{encoded}",
        "createdAt": "2026-05-22T00:00:00Z",
    }
    data.update(overrides)
    return data


def valid_pdf_attachment(**overrides):
    pdf = pymupdf.open()
    page = pdf.new_page()
    page.insert_text((72, 72), "Subagent can search text and inspect rendered pages.")
    payload = pdf.tobytes()
    pdf.close()
    return pdf_attachment(payload, **overrides)


def test_validate_pdf_reference_accepts_pdf_data_url():
    reference = validate_pdf_reference(pdf_attachment())

    assert reference["id"] == "pdf-1"
    assert reference["mimeType"] == "application/pdf"
    assert reference["sizeBytes"] == 9


def test_validate_pdf_reference_rejects_non_pdf():
    with pytest.raises(ValueError, match="application/pdf"):
        validate_pdf_reference(
            pdf_attachment(
                mimeType="image/png",
                dataUrl="data:image/png;base64,ZmFrZQ==",
            )
        )


def test_registry_consult_without_active_document_requests_reupload():
    registry = PdfDocumentRegistry()

    assert "reupload" in asyncio.run(
        registry.consult("testusername", "What is this about?")
    )

    registry.cleanup()


def test_registry_upload_uses_filename_title_when_pdf_has_no_title():
    registry = PdfDocumentRegistry()

    metadata = registry.set_from_attachment(
        "testusername", valid_pdf_attachment(name="sample-book.pdf")
    )

    assert metadata.title == "sample-book"
    assert metadata.pageCount == 1

    registry.cleanup()


def test_read_page_tool_returns_rendered_image_without_debug_file(
    tmp_path, monkeypatch
):
    monkeypatch.chdir(tmp_path)
    registry = PdfDocumentRegistry()
    registry.set_from_attachment("testusername", valid_pdf_attachment())
    document = registry._documents.get("testusername")
    assert document is not None

    contents = document._read_page_tool()(1, True)

    assert any(
        (getattr(content, "uri", None) or "").startswith("data:image/png;base64,")
        for content in contents
    )
    assert not list(tmp_path.glob("image*.png"))

    registry.cleanup()
