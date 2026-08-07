"""Rendus HTML prioritaire et PDF secondaire depuis l'index de figures validé."""

from __future__ import annotations

import html
import importlib.metadata
import tempfile
from dataclasses import dataclass
from pathlib import Path
from typing import Any

from fpdf import FPDF

from effet_fondateur.audit import atomic_write_json, sha256_file
from effet_fondateur.contracts import validate_json_document


DOMAIN_ORDER = (
    "POPULATION_STRUCTURE",
    "FOUNDER_IBS",
    "VARIANT_AGE",
    "LOCAL_LD",
    "ROH",
    "SENSITIVITY",
)

DOMAIN_TITLES = {
    "POPULATION_STRUCTURE": "Structure populationnelle — PCA",
    "FOUNDER_IBS": "Partage IBS autour du variant cible",
    "VARIANT_AGE": "Datation du variant",
    "LOCAL_LD": "LD local secondaire",
    "ROH": "ROH et autozygotie",
    "SENSITIVITY": "Analyses de sensibilité",
}

DOMAIN_LIMITS = {
    "POPULATION_STRUCTURE": "PCA exploratoire interne à l'étude ; aucune attribution automatique d'ascendance.",
    "FOUNDER_IBS": "Le partage allélique observé est IBS et ne constitue pas une preuve IBD.",
    "VARIANT_AGE": "La datation est conditionnelle aux hypothèses Gamma, à la phase et à la carte génétique.",
    "LOCAL_LD": "Le LD est descriptif et ne démontre pas une origine fondatrice unique.",
    "ROH": "L'autozygotie individuelle reste distincte de l'IBS ou de l'IBD entre individus.",
    "SENSITIVITY": "La robustesse aux scénarios testés n'est ni une validation externe ni une preuve causale.",
}


@dataclass(frozen=True)
class RenderPublication:
    """Chemins et manifest des deux formats de consultation."""

    html_path: Path
    pdf_path: Path
    manifest_path: Path


def _ordered_figures(figure_index: dict[str, Any]) -> list[dict[str, Any]]:
    by_domain = {figure["domain"]: figure for figure in figure_index["figures"]}
    if set(by_domain) != set(DOMAIN_ORDER):
        raise ValueError("render_domain_set_mismatch")
    return [by_domain[domain] for domain in DOMAIN_ORDER]


def _validate_render_sources(
    output_dir: Path,
    figure_index: dict[str, Any],
) -> list[dict[str, Any]]:
    figures = _ordered_figures(figure_index)
    for figure in figures:
        provenance_path = output_dir / figure["provenance_path"]
        if not provenance_path.is_file():
            raise ValueError("render_provenance_missing")
        if figure["status"] == "BLOCKED":
            if figure["figure_path"] is not None:
                raise ValueError("blocked_render_has_figure")
            continue
        if figure["figure_path"] is None or not (output_dir / figure["figure_path"]).is_file():
            raise ValueError("render_figure_missing")
    return figures


def _render_html(run_id: str, figures: list[dict[str, Any]], output_path: Path) -> None:
    navigation = "".join(
        f'<a href="#{html.escape(figure["domain"].lower())}">{html.escape(DOMAIN_TITLES[figure["domain"]])}</a>'
        for figure in figures
    )
    sections: list[str] = []
    for figure in figures:
        domain = figure["domain"]
        status_class = figure["status"].lower().replace("_", "-")
        if figure["figure_path"] is None:
            visual = (
                '<div class="unavailable">Figure indisponible — '
                f'{html.escape(figure["blocked_reason"] or figure["status"])}</div>'
            )
        else:
            visual = (
                f'<img src="{html.escape(figure["figure_path"])}" '
                f'alt="{html.escape(DOMAIN_TITLES[domain])}" loading="lazy">'
            )
        sections.append(
            f'<section id="{html.escape(domain.lower())}">'
            f'<header><h2>{html.escape(DOMAIN_TITLES[domain])}</h2>'
            f'<span class="status {status_class}">{html.escape(figure["status"])}</span></header>'
            f'{visual}'
            f'<p class="limit"><strong>Limite :</strong> {html.escape(DOMAIN_LIMITS[domain])}</p>'
            f'<p class="provenance"><a href="{html.escape(figure["provenance_path"])}">Provenance versionnée</a></p>'
            '</section>'
        )
    document = f'''<!doctype html>
<html lang="fr">
<head>
  <meta charset="utf-8">
  <meta name="viewport" content="width=device-width, initial-scale=1">
  <meta name="robots" content="noindex,nofollow,noarchive">
  <meta http-equiv="Content-Security-Policy" content="default-src 'none'; img-src 'self'; style-src 'unsafe-inline'">
  <title>Visualisations consolidées — {html.escape(run_id)}</title>
  <style>
    :root {{ color-scheme: light; --ink:#172554; --muted:#475569; --paper:#f8fafc; --card:#fff; --line:#cbd5e1; }}
    * {{ box-sizing:border-box; }} body {{ margin:0; font-family:Inter,system-ui,sans-serif; color:#1f2937; background:var(--paper); }}
    .hero {{ padding:2.5rem max(2rem,6vw); color:white; background:linear-gradient(135deg,#172554,#0f766e); }}
    .hero h1 {{ margin:0 0 .5rem; font-size:clamp(2rem,5vw,3.7rem); }} .hero p {{ max-width:75ch; line-height:1.55; }}
    nav {{ display:flex; flex-wrap:wrap; gap:.65rem; padding:1rem max(2rem,6vw); position:sticky; top:0; background:#ffffffee; border-bottom:1px solid var(--line); z-index:2; }}
    nav a {{ color:#0f4c5c; text-decoration:none; padding:.45rem .7rem; border:1px solid #99c2c8; border-radius:999px; font-size:.9rem; }}
    main {{ width:min(1280px,92vw); margin:2rem auto 4rem; display:grid; gap:1.5rem; }}
    section {{ background:var(--card); border:1px solid var(--line); border-radius:18px; padding:1.35rem; box-shadow:0 10px 30px #0f172a0d; scroll-margin-top:6rem; }}
    section header {{ display:flex; gap:1rem; align-items:center; justify-content:space-between; }} h2 {{ color:var(--ink); margin:.2rem 0 1rem; }}
    img {{ display:block; width:100%; height:auto; border:1px solid #e2e8f0; border-radius:12px; background:#fbfaf7; }}
    .status {{ font:700 .75rem/1 system-ui; letter-spacing:.06em; padding:.5rem .65rem; border-radius:999px; background:#dcfce7; color:#166534; }}
    .status.not-evaluated {{ background:#ffedd5; color:#9a3412; }} .status.blocked {{ background:#fee2e2; color:#991b1b; }}
    .limit {{ color:#7c2d12; line-height:1.5; }} .provenance a {{ color:#1d4ed8; }}
    .unavailable {{ min-height:220px; display:grid; place-items:center; border:2px dashed #f59e0b; border-radius:12px; background:#fff7ed; color:#9a3412; font-weight:700; }}
    footer {{ padding:2rem max(2rem,6vw); color:#64748b; border-top:1px solid var(--line); }}
    @media print {{ nav {{ display:none; }} section {{ break-after:page; box-shadow:none; border:0; }} .hero {{ background:white; color:#172554; }} }}
  </style>
</head>
<body>
  <header class="hero"><h1>Visualisations consolidées</h1><p>Run <code>{html.escape(run_id)}</code>. Rendu pseudonymisé classé <strong>sensitive_genetic</strong>. Les figures décrivent des résultats audités et ne constituent ni une preuve causale, ni un score composite d'effet fondateur.</p></header>
  <nav aria-label="Domaines">{navigation}</nav>
  <main>{''.join(sections)}</main>
  <footer>Étape 18 — HTML prioritaire, PDF secondaire. Aucun recalcul scientifique.</footer>
</body>
</html>
'''
    output_path.write_text(document, encoding="utf-8")


def _pdf_text(value: str) -> str:
    return value.translate(str.maketrans({"—": "-", "–": "-", "’": "'", "′": "'", "≠": "!="})).encode("latin-1", "replace").decode("latin-1")


def _pdf_safe_svg(source_path: Path, target_path: Path) -> None:
    """Adapte uniquement la typographie SVG à la police interne de fpdf2."""
    source = source_path.read_text(encoding="utf-8")
    target_path.write_text(_pdf_text(source), encoding="utf-8")


def _render_pdf(
    run_id: str,
    figures: list[dict[str, Any]],
    output_dir: Path,
    output_path: Path,
) -> None:
    pdf = FPDF(orientation="L", unit="mm", format="A4")
    pdf.set_title(_pdf_text(f"Visualisations consolidees - {run_id}"))
    pdf.set_author("Effet_fondateur2 - etape 18")
    with tempfile.TemporaryDirectory(prefix="step18_pdf_svg_", dir=output_dir) as temporary_dir:
        temporary_root = Path(temporary_dir)
        for figure in figures:
            domain = figure["domain"]
            pdf.add_page()
            pdf.set_font("Helvetica", "B", 18)
            pdf.cell(0, 10, _pdf_text(DOMAIN_TITLES[domain]), new_x="LMARGIN", new_y="NEXT")
            pdf.set_font("Helvetica", "", 9)
            pdf.cell(0, 6, _pdf_text(f"Run {run_id} | statut {figure['status']} | sensitive_genetic"), new_x="LMARGIN", new_y="NEXT")
            if figure["figure_path"] is None:
                pdf.set_fill_color(255, 247, 237)
                pdf.set_text_color(154, 52, 18)
                pdf.rect(10, 35, 277, 140, style="DF")
                pdf.set_xy(20, 95)
                pdf.set_font("Helvetica", "B", 14)
                pdf.multi_cell(257, 8, _pdf_text(f"Figure indisponible : {figure['blocked_reason'] or figure['status']}"), align="C")
            else:
                source_svg = output_dir / figure["figure_path"]
                pdf_svg = temporary_root / source_svg.name
                _pdf_safe_svg(source_svg, pdf_svg)
                pdf.image(pdf_svg, x=10, y=30, w=277, h=150, keep_aspect_ratio=True)
            pdf.set_xy(10, 185)
            pdf.set_text_color(124, 45, 18)
            pdf.set_font("Helvetica", "", 9)
            pdf.multi_cell(277, 5, _pdf_text(f"Limite : {DOMAIN_LIMITS[domain]}"))
        pdf.output(output_path)


def publish_renderings(
    *,
    run_id: str,
    output_dir: Path,
    figure_index: dict[str, Any],
    figure_index_path: Path,
) -> RenderPublication:
    """Publie HTML puis PDF depuis le même index, sans lire d'artefact scientifique."""
    validate_json_document(figure_index, "figure_index.schema.json")
    figures = _validate_render_sources(output_dir, figure_index)
    html_path = output_dir / "visualization_gallery.html"
    pdf_path = output_dir / "visualization_gallery.pdf"
    _render_html(run_id, figures, html_path)
    _render_pdf(run_id, figures, output_dir, pdf_path)
    manifest = {
        "schema_version": "1.0.0",
        "run_id": run_id,
        "method_id": "validated_current_run_rendering_v1",
        "source_figure_index_sha256": sha256_file(figure_index_path),
        "domain_order": list(DOMAIN_ORDER),
        "html": {
            "path": html_path.name,
            "media_type": "text/html",
            "sha256": sha256_file(html_path),
            "renderer": "python_html_v1",
            "renderer_version": "1.0.0",
        },
        "pdf": {
            "path": pdf_path.name,
            "media_type": "application/pdf",
            "sha256": sha256_file(pdf_path),
            "renderer": "fpdf2_svg_v1",
            "renderer_version": importlib.metadata.version("fpdf2"),
        },
        "pseudonymized": True,
        "sensitivity": "sensitive_genetic",
        "scientific_recalculation_performed": False,
        "composite_founder_score_calculated": False,
    }
    validate_json_document(manifest, "visualization_render_manifest.schema.json")
    manifest_path = output_dir / "visualization_render_manifest.json"
    atomic_write_json(manifest_path, manifest)
    return RenderPublication(html_path=html_path, pdf_path=pdf_path, manifest_path=manifest_path)
