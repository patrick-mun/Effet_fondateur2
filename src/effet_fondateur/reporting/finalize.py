"""Finalisation contrôlée d'un rapport après décision humaine explicite."""

from __future__ import annotations

import argparse
import html
import re
import tempfile
from pathlib import Path
from typing import Sequence

from fpdf import FPDF

from effet_fondateur.audit import atomic_write_json, read_json, sha256_file
from effet_fondateur.contracts import validate_json_document


FORBIDDEN_CLAIMS = re.compile(
    r"\b(prouve|d[ée]montre|confirme\s+(?:un|l['’])?effet\s+fondateur|preuve\s+robuste|origine\s+unique\s+(?:confirm[ée]e?|d[ée]montr[ée]e?|certaine)|certainement)\b",
    re.IGNORECASE,
)


def _contains_forbidden_claim(value: str) -> bool:
    """Détecte une assertion forte sans rejeter les négations protectrices."""
    protected = re.sub(
        r"\b(?:ne|n['’])\s*(?:prouve|d[ée]montre|confirme)\s+pas\b",
        "formulation_prudente",
        value,
        flags=re.IGNORECASE,
    )
    return FORBIDDEN_CLAIMS.search(protected) is not None


def _approved_texts(review: dict) -> dict[str, str]:
    texts: dict[str, str] = {}
    for section in review["sections"]:
        section_id = section.get("section_id")
        approved_text = section.get("approved_text")
        if not isinstance(section_id, str) or not isinstance(approved_text, str) or not approved_text.strip():
            raise ValueError("invalid_or_empty_approved_section")
        if section_id in texts:
            raise ValueError("duplicate_approved_section")
        if _contains_forbidden_claim(approved_text):
            raise ValueError(f"forbidden_assertive_claim:{section_id}")
        texts[section_id] = approved_text.strip()
    return texts


def _locked_html(draft_html: str, review: dict, texts: dict[str, str]) -> str:
    expected_ids = set(re.findall(r'<textarea data-section="([a-z0-9_]+)">', draft_html))
    if set(texts) != expected_ids:
        raise ValueError("review_section_set_mismatch")
    locked = draft_html
    for section_id, approved_text in texts.items():
        paragraphs = "".join(f"<p>{html.escape(part)}</p>" for part in approved_text.split("\n\n") if part.strip())
        pattern = rf'<textarea data-section="{re.escape(section_id)}">.*?</textarea>'
        locked, count = re.subn(pattern, f'<div class="approved-text" data-section="{section_id}">{paragraphs}</div>', locked, count=1, flags=re.DOTALL)
        if count != 1:
            raise ValueError(f"draft_section_missing:{section_id}")
    locked = re.sub(r'<section class="card" id="human-validation">.*?</section>', "", locked, count=1, flags=re.DOTALL)
    locked = re.sub(r'<script src="report_editor\.js"></script>', "", locked, count=1)
    locked = locked.replace("script-src 'self'", "")
    locked = locked.replace("Rapport V2 révisable", "Rapport V2 approuvé")
    approval = f'<section class="card"><h2>Validation éditoriale</h2><p>Approuvé par {html.escape(review["reviewer"])} le {html.escape(review["reviewed_at"])}.</p><p>Cette approbation ne constitue pas une preuve scientifique d’effet fondateur.</p></section>'
    return locked.replace("</main>", approval + "</main>", 1)


def _pdf_text(value: str) -> str:
    return value.translate(str.maketrans({"—": "-", "–": "-", "’": "'", "′": "'", "≠": "!="})).encode("latin-1", "replace").decode("latin-1")


def _render_pdf(report_dir: Path, review: dict, texts: dict[str, str], draft_html: str, output_path: Path) -> None:
    pdf = FPDF(unit="mm", format="A4")
    pdf.set_title(_pdf_text("Rapport V2 approuve")); pdf.set_author("Effet_fondateur2 - etape 19")
    pdf.add_page(); pdf.set_font("Helvetica", "B", 19); pdf.multi_cell(0, 10, _pdf_text("Rapport V2 approuve"), new_x="LMARGIN", new_y="NEXT")
    pdf.set_font("Helvetica", "", 10); pdf.multi_cell(0, 6, _pdf_text(f"Run {review['run_id']} | Relecteur {review['reviewer']} | {review['reviewed_at']}"), new_x="LMARGIN", new_y="NEXT")
    pdf.ln(3); pdf.set_text_color(124, 45, 18); pdf.multi_cell(0, 6, _pdf_text("L'approbation editoriale ne constitue pas une preuve scientifique d'effet fondateur."), new_x="LMARGIN", new_y="NEXT")
    for section_id, approved_text in texts.items():
        pdf.add_page(); pdf.set_text_color(23, 37, 84); pdf.set_font("Helvetica", "B", 15); pdf.multi_cell(0, 8, _pdf_text(section_id.replace("_", " ").title()), new_x="LMARGIN", new_y="NEXT")
        pdf.set_font("Helvetica", "", 10); pdf.multi_cell(0, 6, _pdf_text(approved_text), new_x="LMARGIN", new_y="NEXT")
    svg_documents = []
    network = report_dir / "kinship_network.svg"
    if network.is_file():
        svg_documents.append(network.read_text(encoding="utf-8"))
    svg_documents.extend(re.findall(r'(<svg\b.*?</svg>)', draft_html, flags=re.DOTALL))
    with tempfile.TemporaryDirectory(prefix="step19_final_pdf_", dir=report_dir) as temporary:
        temporary_root = Path(temporary)
        for index, svg in enumerate(svg_documents, 1):
            svg_path = temporary_root / f"figure_{index:02d}.svg"; svg_path.write_text(_pdf_text(svg), encoding="utf-8")
            pdf.add_page(orientation="L"); pdf.image(svg_path, x=10, y=15, w=277, h=175, keep_aspect_ratio=True)
    pdf.output(output_path)


def finalize_report(report_dir: Path, review_path: Path) -> tuple[Path, Path, Path]:
    """Vérifie une relecture et publie les formats finaux sans toucher aux faits."""
    facts_path, draft_path = report_dir / "interpretation_facts.json", report_dir / "interpretation_draft.json"
    review = read_json(review_path); validate_json_document(review, "report_review.schema.json")
    facts, draft = read_json(facts_path), read_json(draft_path)
    if review["run_id"] != facts["run_id"] or review["facts_sha256"] != sha256_file(facts_path) or review["draft_sha256"] != sha256_file(draft_path):
        raise ValueError("review_provenance_mismatch")
    texts = _approved_texts(review)
    draft_html_path = report_dir / "report_draft.html"; draft_html = draft_html_path.read_text(encoding="utf-8")
    final_html_path = report_dir / "report_final.html"; final_html_path.write_text(_locked_html(draft_html, review, texts), encoding="utf-8")
    final_pdf_path = report_dir / "report_final.pdf"; _render_pdf(report_dir, review, texts, draft_html, final_pdf_path)
    validation = {"schema_version": "1.0.0", "run_id": facts["run_id"], "method_id": "reviewable_scientific_report_v1", "status": "READY_FOR_SCIENTIFIC_REVIEW", "facts_sha256": sha256_file(facts_path), "draft_sha256": sha256_file(draft_path), "human_validation_required": False, "checks": [{"check": "review_schema", "status": "PASS"}, {"check": "review_provenance", "status": "PASS"}, {"check": "assertive_language_guard", "status": "PASS"}], "blocking_reasons": []}
    validate_json_document(validation, "report_validation.schema.json")
    validation_path = report_dir / "report_final_validation.json"; atomic_write_json(validation_path, validation)
    return final_html_path, final_pdf_path, validation_path


def main(arguments: Sequence[str] | None = None) -> int:
    parser = argparse.ArgumentParser(prog="finalize-report")
    parser.add_argument("--report-dir", type=Path, required=True); parser.add_argument("--review", type=Path, required=True)
    parsed = parser.parse_args(arguments)
    try:
        finalize_report(parsed.report_dir, parsed.review)
        return 0
    except (ValueError, OSError) as error:
        parser.error(str(error))
        return 2


if __name__ == "__main__":
    raise SystemExit(main())
