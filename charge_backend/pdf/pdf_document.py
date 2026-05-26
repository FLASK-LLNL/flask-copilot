from __future__ import annotations

import base64
import hashlib
import json
import math
import os
import re
import sqlite3
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import pymupdf


TOKEN_RE = re.compile(r"[A-Za-z0-9][A-Za-z0-9'_-]*")


@dataclass(frozen=True)
class PageRead:
    page_number: int
    text: str
    image_data_url: str | None = None
    image_reason: str | None = None


@dataclass(frozen=True)
class SearchHit:
    page_number: int
    paragraph: int
    text: str
    score: float | None = None


@dataclass(frozen=True)
class TocEntry:
    level: int
    title: str
    page_number: int


@dataclass(frozen=True)
class DocumentOverview:
    title: str
    page_count: int
    toc: list[TocEntry]
    toc_excerpt: str | None
    index_start_page: int | None
    index_excerpt: str | None

    def prompt_context(
        self, max_toc_chars: int = 18_000, max_index_chars: int = 24_000
    ) -> str:
        toc_text = self._format_toc(max_toc_chars)
        if self.index_start_page and self.index_excerpt:
            index_text = _truncate(self.index_excerpt, max_index_chars)
            index_section = (
                f"Detected back-of-book index starting near page {self.index_start_page}:\n"
                f"{index_text}"
            )
        else:
            index_section = (
                "No back-of-book index was detected from the PDF text layer."
            )

        return (
            f"Document title: {self.title}\n"
            f"PDF pages: {self.page_count}\n\n"
            f"Table of contents:\n{toc_text}\n\n"
            f"{index_section}"
        )

    def _format_toc(self, max_chars: int) -> str:
        if not self.toc:
            if self.toc_excerpt:
                return (
                    "No embedded PDF table of contents was found, but this printed ToC-like "
                    f"front-matter excerpt was detected:\n{_truncate(self.toc_excerpt, max_chars)}"
                )
            return "No embedded PDF table of contents or printed ToC-like excerpt was found."

        lines = []
        for entry in self.toc:
            indent = "  " * max(entry.level - 1, 0)
            lines.append(f"{indent}- {entry.title} (p. {entry.page_number})")
        return _truncate("\n".join(lines), max_chars)


class SubagentDocumentPdf:
    """PDF text, image, and search-index access for one document."""

    def __init__(
        self, pdf_path: str | Path, cache_dir: str | Path | None = None
    ) -> None:
        self.pdf_path = Path(pdf_path).expanduser().resolve()
        if not self.pdf_path.exists():
            raise FileNotFoundError(f"PDF not found: {self.pdf_path}")
        if self.pdf_path.suffix.lower() != ".pdf":
            raise ValueError(f"Expected a PDF file, got: {self.pdf_path}")

        self.cache_dir = Path(cache_dir or default_cache_dir()).expanduser()
        self.cache_dir.mkdir(parents=True, exist_ok=True)
        self.db_path = self.cache_dir / f"{self._fingerprint()}.sqlite3"
        self._ensure_index()

    @property
    def page_count(self) -> int:
        with pymupdf.open(self.pdf_path) as doc:
            return doc.page_count

    def overview(self) -> DocumentOverview:
        with pymupdf.open(self.pdf_path) as doc:
            metadata_title = (doc.metadata or {}).get("title") or ""
            title = metadata_title.strip() or self.pdf_path.stem
            toc = [
                TocEntry(
                    level=int(level), title=str(title).strip(), page_number=int(page)
                )
                for level, title, page in doc.get_toc(simple=True)
                if str(title).strip()
            ]
            index_start, index_excerpt = self._detect_index(doc)
            toc_excerpt = None if toc else self._detect_printed_toc(doc)
            return DocumentOverview(
                title=title,
                page_count=doc.page_count,
                toc=toc,
                toc_excerpt=toc_excerpt,
                index_start_page=index_start,
                index_excerpt=index_excerpt,
            )

    def read_page(self, page_number: int, force_image: bool = False) -> PageRead:
        with pymupdf.open(self.pdf_path) as doc:
            if page_number < 1 or page_number > doc.page_count:
                raise ValueError(f"page_number must be between 1 and {doc.page_count}")

            page = doc.load_page(page_number - 1)
            text = _clean_text(page.get_text("text", sort=True))
            image_reason = None
            image_data_url = None
            if force_image:
                image_reason = (
                    "force_image was requested, so a rendered page image is included."
                )
                image_data_url = self._render_page_data_url(page)
            elif _needs_ocr_image(text):
                image_reason = "The extracted text layer is sparse, so a rendered page image is included for OCR."
                image_data_url = self._render_page_data_url(page)

            return PageRead(
                page_number=page_number,
                text=text,
                image_data_url=image_data_url,
                image_reason=image_reason,
            )

    def search(
        self,
        query: str,
        num_exact_results: int = 10,
        num_fuzzy_results: int = 3,
    ) -> tuple[list[SearchHit], list[SearchHit]]:
        normalized_query = query.strip()
        if not normalized_query:
            return [], []

        exact = self._exact_search(normalized_query, max(0, num_exact_results))
        fuzzy = self._bm25_search(
            normalized_query, max(0, num_fuzzy_results), exclude=exact
        )
        return exact, fuzzy

    def _ensure_index(self) -> None:
        if self.db_path.exists() and self._is_current_index():
            return

        conn = sqlite3.connect(self.db_path)
        try:
            conn.execute("PRAGMA journal_mode=WAL")
            conn.execute("DROP TABLE IF EXISTS meta")
            conn.execute("DROP TABLE IF EXISTS pages")
            conn.execute("DROP TABLE IF EXISTS chunks")
            conn.execute(
                "CREATE TABLE meta (key TEXT PRIMARY KEY, value TEXT NOT NULL)"
            )
            conn.execute(
                "CREATE TABLE pages (page_number INTEGER PRIMARY KEY, text TEXT NOT NULL)"
            )
            conn.execute(
                "CREATE TABLE chunks ("
                "id INTEGER PRIMARY KEY AUTOINCREMENT, "
                "page_number INTEGER NOT NULL, "
                "paragraph INTEGER NOT NULL, "
                "text TEXT NOT NULL"
                ")"
            )

            with pymupdf.open(self.pdf_path) as doc:
                for page_index in range(doc.page_count):
                    page_number = page_index + 1
                    text = _clean_text(
                        doc.load_page(page_index).get_text("text", sort=True)
                    )
                    conn.execute(
                        "INSERT INTO pages (page_number, text) VALUES (?, ?)",
                        (page_number, text),
                    )
                    for para_index, paragraph in enumerate(
                        _paragraphs_from_text(text), start=1
                    ):
                        conn.execute(
                            "INSERT INTO chunks (page_number, paragraph, text) VALUES (?, ?, ?)",
                            (page_number, para_index, paragraph),
                        )

            conn.execute(
                "INSERT INTO meta (key, value) VALUES ('fingerprint', ?)",
                (self._fingerprint(),),
            )
            conn.commit()
        finally:
            conn.close()

    def _is_current_index(self) -> bool:
        try:
            conn = sqlite3.connect(self.db_path)
            try:
                row = conn.execute(
                    "SELECT value FROM meta WHERE key = 'fingerprint'"
                ).fetchone()
            finally:
                conn.close()
            return bool(row and row[0] == self._fingerprint())
        except sqlite3.Error:
            return False

    def _exact_search(self, query: str, limit: int) -> list[SearchHit]:
        if limit == 0:
            return []

        query_lower = query.casefold()
        hits: list[SearchHit] = []
        conn = sqlite3.connect(self.db_path)
        try:
            rows = conn.execute(
                "SELECT page_number, paragraph, text FROM chunks ORDER BY page_number, paragraph"
            ).fetchall()
        finally:
            conn.close()

        for page_number, paragraph, text in rows:
            if query_lower in text.casefold():
                hits.append(
                    SearchHit(
                        page_number=int(page_number),
                        paragraph=int(paragraph),
                        text=_highlight_excerpt(text, query),
                    )
                )
                if len(hits) >= limit:
                    break
        return hits

    def _bm25_search(
        self,
        query: str,
        limit: int,
        exclude: list[SearchHit],
    ) -> list[SearchHit]:
        if limit == 0:
            return []

        query_terms = _tokenize(query)
        if not query_terms:
            return []

        conn = sqlite3.connect(self.db_path)
        try:
            rows = conn.execute(
                "SELECT page_number, paragraph, text FROM chunks WHERE length(text) > 20"
            ).fetchall()
        finally:
            conn.close()

        excluded_keys = {(hit.page_number, hit.paragraph) for hit in exclude}
        documents: list[tuple[int, int, str, list[str]]] = []
        doc_freq: dict[str, int] = {}
        total_len = 0

        for page_number, paragraph, text in rows:
            key = (int(page_number), int(paragraph))
            if key in excluded_keys:
                continue
            terms = _tokenize(text)
            if not terms:
                continue
            documents.append((key[0], key[1], text, terms))
            total_len += len(terms)
            for term in set(terms):
                doc_freq[term] = doc_freq.get(term, 0) + 1

        if not documents:
            return []

        avgdl = total_len / len(documents)
        k1 = 1.5
        b = 0.75
        scored: list[SearchHit] = []

        for page_number, paragraph, text, terms in documents:
            term_counts: dict[str, int] = {}
            for term in terms:
                term_counts[term] = term_counts.get(term, 0) + 1

            score = 0.0
            doc_len = len(terms)
            for term in query_terms:
                if term not in term_counts:
                    continue
                df = doc_freq.get(term, 0)
                idf = math.log(1 + (len(documents) - df + 0.5) / (df + 0.5))
                tf = term_counts[term]
                denom = tf + k1 * (1 - b + b * doc_len / avgdl)
                score += idf * (tf * (k1 + 1)) / denom

            if score > 0:
                scored.append(
                    SearchHit(
                        page_number=page_number,
                        paragraph=paragraph,
                        text=_truncate(text, 900),
                        score=score,
                    )
                )

        scored.sort(key=lambda hit: (hit.score or 0.0, -hit.page_number), reverse=True)
        return scored[:limit]

    def _detect_index(self, doc: pymupdf.Document) -> tuple[int | None, str | None]:
        start_scan = max(0, int(doc.page_count * 0.75) - 1)
        index_start: int | None = None

        for page_index in range(start_scan, doc.page_count):
            text = _clean_text(doc.load_page(page_index).get_text("text", sort=True))
            first_lines = "\n".join(text.splitlines()[:8])
            if re.search(
                r"(?im)^\s(index|general index|name index|subject index)\s*$",
                first_lines,
            ):
                index_start = page_index + 1
                break

        if index_start is None:
            return None, None

        excerpts = []
        for page_number in range(index_start, doc.page_count + 1):
            text = _clean_text(
                doc.load_page(page_number - 1).get_text("text", sort=True)
            )
            if text:
                excerpts.append(f"[p. {page_number}]\n{text}")

        return index_start, "\n\n".join(excerpts) if excerpts else None

    def _detect_printed_toc(self, doc: pymupdf.Document) -> str | None:
        scan_pages = min(doc.page_count, 20)
        for page_index in range(scan_pages):
            text = _clean_text(doc.load_page(page_index).get_text("text", sort=True))
            if not text:
                continue
            first_lines = "\n".join(text.splitlines()[:12])
            if not re.search(r"(?im)^\s(table of contents|contents)\s*$", first_lines):
                continue

            excerpts = []
            for toc_page_index in range(
                page_index, min(doc.page_count, page_index + 8)
            ):
                page_number = toc_page_index + 1
                page_text = _clean_text(
                    doc.load_page(toc_page_index).get_text("text", sort=True)
                )
                if page_text:
                    excerpts.append(f"[p. {page_number}]\n{page_text}")
            return "\n\n".join(excerpts) if excerpts else None
        return None

    def _render_page_data_url(self, page: pymupdf.Page) -> str:
        dpi = int(os.environ.get("SUBAGENT_RENDER_DPI", "160"))
        for candidate_dpi in (dpi, 120, 96):
            scale = candidate_dpi / 72
            pixmap = page.get_pixmap(matrix=pymupdf.Matrix(scale, scale), alpha=False)
            png_bytes = pixmap.tobytes("png")

            if len(png_bytes) < 4_500_000 or candidate_dpi == 96:
                b64 = base64.b64encode(png_bytes).decode("ascii")
                return f"data:image/png;base64,{b64}"
        raise RuntimeError("unreachable")

    def _fingerprint(self) -> str:
        stat = self.pdf_path.stat()
        raw = json.dumps(
            {
                "path": str(self.pdf_path),
                "size": stat.st_size,
                "mtime_ns": stat.st_mtime_ns,
            },
            sort_keys=True,
        ).encode("utf-8")
        return hashlib.sha256(raw).hexdigest()


def default_cache_dir() -> Path:
    return Path(
        os.environ.get(
            "SUBAGENT_CACHE_DIR", Path.home() / ".cache" / "scholar-subagent"
        )
    )


def format_search_results(exact: list[SearchHit], fuzzy: list[SearchHit]) -> str:
    parts = ["Exact paragraph matches:"]
    if exact:
        parts.extend(_format_hit(hit) for hit in exact)
    else:
        parts.append("No exact paragraph matches.")

    parts.append("\nFuzzy/RAG paragraph results:")
    if fuzzy:
        parts.extend(_format_hit(hit, include_score=True) for hit in fuzzy)
    else:
        parts.append("No fuzzy/RAG results.")

    return "\n".join(parts)


def page_read_to_json(page: PageRead) -> dict[str, Any]:
    payload: dict[str, Any] = {
        "page_number": page.page_number,
        "text": page.text or "",
    }
    if page.image_data_url:
        payload["image_included_in_next_message"] = True
        payload["image_reason"] = page.image_reason
    return payload


def _format_hit(hit: SearchHit, include_score: bool = False) -> str:
    score = f" score={hit.score:.3f}" if include_score and hit.score is not None else ""
    return f"- p. {hit.page_number}, paragraph {hit.paragraph}{score}: {hit.text}"


def _clean_text(text: str) -> str:
    text = text.replace("\x00", "")
    text = re.sub(r"[ \t]+", " ", text)
    text = re.sub(r"\n{3,}", "\n\n", text)
    return text.strip()


def _paragraphs_from_text(text: str) -> list[str]:
    text = _clean_text(text)
    if not text:
        return []

    raw_paragraphs = [
        part.strip() for part in re.split(r"\n\s*\n", text) if part.strip()
    ]
    if len(raw_paragraphs) > 1:
        return [_normalize_paragraph(part) for part in raw_paragraphs if len(part) > 10]

    lines = [line.strip() for line in text.splitlines() if line.strip()]
    paragraphs: list[str] = []
    current: list[str] = []
    for line in lines:
        current.append(line)
        line_ends_sentence = bool(re.search(r"[.!?]['\")\]]?$", line))
        if line_ends_sentence and len(" ".join(current)) >= 160:
            paragraphs.append(_normalize_paragraph(" ".join(current)))
            current = []
    if current:
        paragraphs.append(_normalize_paragraph(" ".join(current)))
    return [paragraph for paragraph in paragraphs if len(paragraph) > 10]


def _normalize_paragraph(text: str) -> str:
    return re.sub(r"\s+", " ", text).strip()


def _tokenize(text: str) -> list[str]:
    return [match.group(0).casefold() for match in TOKEN_RE.finditer(text)]


def _highlight_excerpt(text: str, query: str, radius: int = 450) -> str:
    lowered = text.casefold()
    query_lower = query.casefold()
    pos = lowered.find(query_lower)
    if pos < 0:
        return _truncate(text, radius * 2)
    start = max(0, pos - radius)
    end = min(len(text), pos + len(query) + radius)
    prefix = "..." if start > 0 else ""
    suffix = "..." if end < len(text) else ""
    return f"{prefix}{text[start:end]}{suffix}"


def _truncate(text: str, max_chars: int) -> str:
    if len(text) <= max_chars:
        return text
    return text[: max_chars - 40].rstrip() + "\n...[truncated]..."


def _needs_ocr_image(text: str) -> bool:
    if len(text) < 500:
        return True
    alnum = sum(char.isalnum() for char in text)
    return alnum < 250
