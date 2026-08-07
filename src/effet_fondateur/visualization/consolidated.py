"""Rendu SVG audité des six domaines scientifiques consolidés."""

from __future__ import annotations

import html
import re
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Callable

from effet_fondateur.audit import atomic_write_json, read_json, sha256_file
from effet_fondateur.contracts import validate_json_document, validate_tsv_table


METHOD_ID = "validated_current_run_figures_v1"
SCRIPT_VERSION = "1.0.0"


class VisualizationContractError(ValueError):
    """Signale une entrée graphique incohérente sans exposer de donnée."""


@dataclass(frozen=True)
class FigureResult:
    """Décrit une figure, sa provenance et son statut contractuel."""

    figure_id: str
    domain: str
    status: str
    figure_path: Path | None
    provenance_path: Path
    blocked_reason: str | None


def _artifact_path(run_dir: Path, artifact: dict[str, Any]) -> Path:
    relative = Path(artifact["path"])
    if relative.is_absolute() or ".." in relative.parts or relative.parts[:2] == ("data", "output"):
        raise VisualizationContractError("source_path_not_confined_to_current_run")
    path = (run_dir / relative).resolve()
    try:
        path.relative_to(run_dir.resolve())
    except ValueError as error:
        raise VisualizationContractError("source_path_not_confined_to_current_run") from error
    if not path.is_file() or sha256_file(path) != artifact["sha256"]:
        raise VisualizationContractError("source_missing_or_checksum_mismatch")
    return path


def _source_record(artifact: dict[str, Any]) -> dict[str, str]:
    return {key: artifact[key] for key in ("artifact_id", "producer_stage", "producer_signature", "schema_name", "schema_version", "sha256")}


def _validate_sources(run_dir: Path, artifacts: dict[str, dict[str, Any]], expected: dict[str, tuple[str, str]]) -> dict[str, Path]:
    paths: dict[str, Path] = {}
    signatures: dict[str, str] = {}
    descriptor_values: dict[str, set[str]] = {
        "assembly": set(), "sample_set_id": set(), "variant_set_id": set(),
    }
    for artifact_id, (producer, schema_name) in expected.items():
        artifact = artifacts.get(artifact_id)
        if artifact is None or artifact.get("producer_stage") != producer:
            raise VisualizationContractError("source_producer_mismatch")
        if artifact.get("schema_name") != schema_name or artifact.get("schema_version") != "1.0.0":
            raise VisualizationContractError("source_schema_mismatch")
        signature = artifact.get("producer_signature")
        if not isinstance(signature, str) or len(signature) != 64:
            raise VisualizationContractError("source_signature_invalid")
        previous = signatures.setdefault(producer, signature)
        if previous != signature:
            raise VisualizationContractError("source_signature_mismatch")
        for descriptor_name in descriptor_values:
            value = artifact.get(descriptor_name)
            if value is not None:
                descriptor_values[descriptor_name].add(value)
        paths[artifact_id] = _artifact_path(run_dir, artifact)
    if any(len(values) > 1 for values in descriptor_values.values()):
        raise VisualizationContractError("source_dataset_or_assembly_mismatch")
    return paths


def _svg(title: str, subtitle: str, lines: list[str], footer: str) -> str:
    height = 150 + 25 * len(lines)
    escaped = lambda value: html.escape(str(value), quote=True)
    body = [f'<svg xmlns="http://www.w3.org/2000/svg" width="1000" height="{height}" viewBox="0 0 1000 {height}">',
            '<rect width="100%" height="100%" fill="#fbfaf7"/>',
            f'<text x="40" y="45" font-family="sans-serif" font-size="25" font-weight="700" fill="#172554">{escaped(title)}</text>',
            f'<text x="40" y="75" font-family="sans-serif" font-size="15" fill="#475569">{escaped(subtitle)}</text>']
    for index, line in enumerate(lines):
        color = "#9a3412" if "NOT_" in line or "EXCLUDED" in line or "BLOCKED" in line else "#1f2937"
        body.append(f'<text x="55" y="{115 + index * 25}" font-family="monospace" font-size="14" fill="{color}">{escaped(line)}</text>')
    body.append(f'<text x="40" y="{height - 25}" font-family="sans-serif" font-size="13" fill="#7c2d12">{escaped(footer)}</text>')
    body.append("</svg>\n")
    return "".join(body)


def _population_svg(subtitle: str, lines: list[str], footer: str) -> str:
    """Rend un scree plot et une projection PC1–PC2 sans identifiant source."""
    eigen_rows = [line.split("|") for line in lines if line.startswith("EIGEN|")]
    point_rows = [line.split("|") for line in lines if line.startswith("POINT|")]
    summary_rows = [line.removeprefix("SUMMARY|") for line in lines if line.startswith("SUMMARY|")]
    body = [
        '<svg xmlns="http://www.w3.org/2000/svg" width="1100" height="650" viewBox="0 0 1100 650">',
        '<rect width="100%" height="100%" fill="#fbfaf7"/>',
        '<text x="40" y="42" font-family="sans-serif" font-size="25" font-weight="700" fill="#172554">Structure populationnelle — PCA exploratoire</text>',
        f'<text x="40" y="70" font-family="sans-serif" font-size="15" fill="#475569">{html.escape(subtitle)}</text>',
        '<text x="55" y="112" font-family="sans-serif" font-size="18" font-weight="700">Variance expliquée</text>',
        '<line x1="70" y1="310" x2="455" y2="310" stroke="#64748b"/><line x1="70" y1="135" x2="70" y2="310" stroke="#64748b"/>',
    ]
    max_ratio = max((float(row[2]) for row in eigen_rows), default=1.0) or 1.0
    for index, row in enumerate(eigen_rows):
        x = 95 + index * 65
        bar_height = 145 * float(row[2]) / max_ratio
        body.extend([
            f'<rect x="{x}" y="{310 - bar_height:.2f}" width="38" height="{bar_height:.2f}" fill="#2563eb" opacity="0.8"/>',
            f'<text x="{x}" y="330" font-family="sans-serif" font-size="12">{html.escape(row[1])}</text>',
            f'<text x="{x - 3}" y="{300 - bar_height:.2f}" font-family="sans-serif" font-size="11">{100 * float(row[2]):.1f}%</text>',
        ])
    body.extend([
        '<text x="560" y="112" font-family="sans-serif" font-size="18" font-weight="700">Projection PC1–PC2</text>',
        '<rect x="560" y="135" width="480" height="360" fill="#ffffff" stroke="#94a3b8"/>',
    ])
    point_values = [(row, float(row[2]), float(row[3])) for row in point_rows]
    pc1_values = [item[1] for item in point_values] or [0.0]
    pc2_values = [item[2] for item in point_values] or [0.0]
    pc1_min, pc1_max = min(pc1_values), max(pc1_values)
    pc2_min, pc2_max = min(pc2_values), max(pc2_values)
    pc1_span, pc2_span = pc1_max - pc1_min or 1.0, pc2_max - pc2_min or 1.0
    for row, pc1, pc2 in point_values:
        x = 590 + 420 * (pc1 - pc1_min) / pc1_span
        y = 465 - 300 * (pc2 - pc2_min) / pc2_span
        fill = "#2563eb" if row[4] == "REFERENCE" else "#f59e0b"
        stroke = "#dc2626" if row[5] == "OUTLIER" else "#ffffff"
        stroke_width = 4 if row[5] == "OUTLIER" else 1
        body.extend([
            f'<circle cx="{x:.2f}" cy="{y:.2f}" r="8" fill="{fill}" stroke="{stroke}" stroke-width="{stroke_width}"/>',
            f'<text x="{x + 10:.2f}" y="{y - 7:.2f}" font-family="monospace" font-size="10">{html.escape(row[1])}</text>',
        ])
    body.extend([
        '<text x="775" y="520" font-family="sans-serif" font-size="12">PC1</text>',
        '<text x="535" y="320" transform="rotate(-90 535 320)" font-family="sans-serif" font-size="12">PC2</text>',
        '<circle cx="590" cy="550" r="7" fill="#2563eb"/><text x="605" y="555" font-family="sans-serif" font-size="12">Référence indépendante</text>',
        '<circle cx="790" cy="550" r="7" fill="#f59e0b"/><text x="805" y="555" font-family="sans-serif" font-size="12">Individu projeté</text>',
        '<circle cx="940" cy="550" r="8" fill="none" stroke="#dc2626" stroke-width="3"/><text x="955" y="555" font-family="sans-serif" font-size="12">Outlier exploratoire</text>',
    ])
    for index, summary in enumerate(summary_rows):
        body.append(f'<text x="55" y="{385 + 24 * index}" font-family="monospace" font-size="13" fill="#334155">{html.escape(summary)}</text>')
    body.append(f'<text x="40" y="625" font-family="sans-serif" font-size="13" fill="#7c2d12">{html.escape(footer)}</text></svg>\n')
    return "".join(body)


def _chart_header(title: str, subtitle: str) -> list[str]:
    """Construit l'en-tête commun des graphiques scientifiques SVG."""
    return [
        '<svg xmlns="http://www.w3.org/2000/svg" width="1100" height="650" viewBox="0 0 1100 650">',
        '<rect width="100%" height="100%" fill="#fbfaf7"/>',
        f'<text x="40" y="42" font-family="sans-serif" font-size="25" font-weight="700" fill="#172554">{html.escape(title)}</text>',
        f'<text x="40" y="70" font-family="sans-serif" font-size="15" fill="#475569">{html.escape(subtitle)}</text>',
    ]


def _chart_footer(body: list[str], footer: str) -> str:
    body.append(f'<text x="40" y="625" font-family="sans-serif" font-size="13" fill="#7c2d12">{html.escape(footer)}</text></svg>\n')
    return "".join(body)


def _founder_svg(subtitle: str, lines: list[str], footer: str) -> str:
    """Trace les longueurs IBS gauche/droite autour du variant cible."""
    pattern = re.compile(r"^(UNIT-\d+) \| left=([0-9.]+) cM \| right=([0-9.]+) cM \| (.+)$")
    rows = [(match.group(1), float(match.group(2)), float(match.group(3)), match.group(4)) for line in lines if (match := pattern.match(line))]
    if not rows:
        return _svg("FOUNDER IBS", subtitle, lines, footer)
    body = _chart_header("Partage IBS autour du variant cible", subtitle)
    body.extend([
        '<text x="55" y="112" font-family="sans-serif" font-size="15" font-weight="700">Longueur partagée autour de la cible (cM)</text>',
        '<line x1="500" y1="135" x2="500" y2="485" stroke="#dc2626" stroke-width="2"/>',
        '<text x="478" y="515" font-family="sans-serif" font-size="12" fill="#991b1b">cible</text>',
    ])
    maximum = max((max(left, right) for _, left, right, _ in rows), default=1.0) or 1.0
    for index, (unit, left, right, status) in enumerate(rows):
        y = 165 + index * 65
        left_x = 500 - 340 * left / maximum
        right_x = 500 + 340 * right / maximum
        body.extend([
            f'<text x="55" y="{y + 5}" font-family="monospace" font-size="13">{html.escape(unit)}</text>',
            f'<line x1="{left_x:.2f}" y1="{y}" x2="{right_x:.2f}" y2="{y}" stroke="#2563eb" stroke-width="13" stroke-linecap="round"/>',
            f'<circle cx="500" cy="{y}" r="7" fill="#dc2626"/>',
            f'<text x="865" y="{y + 5}" font-family="sans-serif" font-size="12">{left:g} + {right:g} cM · {html.escape(status)}</text>',
        ])
    for index, summary in enumerate(line for line in lines if not line.startswith("UNIT-")):
        body.append(f'<text x="55" y="{545 + 20 * index}" font-family="monospace" font-size="12" fill="#475569">{html.escape(summary)}</text>')
    return _chart_footer(body, footer)


def _age_svg(subtitle: str, lines: list[str], footer: str) -> str:
    """Trace les estimations de datation et leurs intervalles en générations."""
    pattern = re.compile(r"^(PRIMARY|EXPLORATORY) (\S+) \| n=(\d+) \| estimate=([0-9.]+) generations \| CI=([0-9.]+)\.\.([0-9.]+)$")
    rows = [(match.group(1), match.group(2), int(match.group(3)), float(match.group(4)), float(match.group(5)), float(match.group(6))) for line in lines if (match := pattern.match(line))]
    if not rows:
        return _svg("VARIANT AGE", subtitle, lines, footer)
    body = _chart_header("Datation du variant", subtitle)
    body.extend([
        '<text x="55" y="112" font-family="sans-serif" font-size="15" font-weight="700">Estimation et intervalle de confiance (générations)</text>',
        '<line x1="300" y1="500" x2="990" y2="500" stroke="#64748b"/>',
    ])
    maximum = max((upper for *_, upper in rows), default=1.0) or 1.0
    for tick in range(5):
        value = maximum * tick / 4
        x = 300 + 690 * tick / 4
        body.extend([f'<line x1="{x:.2f}" y1="495" x2="{x:.2f}" y2="505" stroke="#64748b"/>', f'<text x="{x - 10:.2f}" y="525" font-family="sans-serif" font-size="11">{value:.1f}</text>'])
    for index, (role, model, count, estimate, lower, upper) in enumerate(rows):
        y = 180 + index * 105
        lower_x, estimate_x, upper_x = (300 + 690 * value / maximum for value in (lower, estimate, upper))
        color = "#0f766e" if role == "PRIMARY" else "#d97706"
        body.extend([
            f'<text x="55" y="{y - 8}" font-family="sans-serif" font-size="14" font-weight="700" fill="{color}">{role} {html.escape(model)}</text>',
            f'<text x="55" y="{y + 15}" font-family="sans-serif" font-size="12">n={count} · {estimate:g} [{lower:g}–{upper:g}] générations</text>',
            f'<line x1="{lower_x:.2f}" y1="{y}" x2="{upper_x:.2f}" y2="{y}" stroke="{color}" stroke-width="5"/>',
            f'<line x1="{lower_x:.2f}" y1="{y - 9}" x2="{lower_x:.2f}" y2="{y + 9}" stroke="{color}"/>',
            f'<line x1="{upper_x:.2f}" y1="{y - 9}" x2="{upper_x:.2f}" y2="{y + 9}" stroke="{color}"/>',
            f'<circle cx="{estimate_x:.2f}" cy="{y}" r="8" fill="{color}"/>',
        ])
    return _chart_footer(body, footer)


def _ld_svg(subtitle: str, lines: list[str], footer: str) -> str:
    """Compare r² génotypique et |D′| par classe de distance."""
    pattern = re.compile(r"^(.*?) \| n=(\d+) \| ([0-9.]+)\.\.([0-9.]+) cM \| r2=([0-9.]+) \| Dprime=([0-9.]+) \| (.+)$")
    rows = [(match.group(1), int(match.group(2)), float(match.group(3)), float(match.group(4)), float(match.group(5)), float(match.group(6)), match.group(7)) for line in lines if (match := pattern.match(line))]
    if not rows:
        return _svg("LOCAL LD", subtitle, lines, footer)
    body = _chart_header("LD local secondaire", subtitle)
    body.extend([
        '<text x="55" y="112" font-family="sans-serif" font-size="15" font-weight="700">Médianes descriptives par classe de distance</text>',
        '<line x1="220" y1="500" x2="1010" y2="500" stroke="#64748b"/>',
        '<line x1="220" y1="145" x2="220" y2="500" stroke="#64748b"/>',
    ])
    group_width = 700 / max(len(rows), 1)
    for index, (cohort, count, low, high, r2, dprime, status) in enumerate(rows):
        center = 270 + index * group_width
        body.extend([
            f'<rect x="{center - 34:.2f}" y="{500 - 320 * r2:.2f}" width="30" height="{320 * r2:.2f}" fill="#2563eb"/>',
            f'<rect x="{center + 6:.2f}" y="{500 - 320 * dprime:.2f}" width="30" height="{320 * dprime:.2f}" fill="#0f766e"/>',
            f'<text x="{center - 45:.2f}" y="525" font-family="sans-serif" font-size="11">{low:g}–{high:g} cM</text>',
            f'<text x="{center - 45:.2f}" y="545" font-family="sans-serif" font-size="10">n={count} · {html.escape(cohort)}</text>',
            f'<text x="{center - 45:.2f}" y="565" font-family="sans-serif" font-size="10">{html.escape(status)}</text>',
        ])
    body.extend([
        '<rect x="730" y="105" width="16" height="16" fill="#2563eb"/><text x="753" y="118" font-family="sans-serif" font-size="12">r² génotypique</text>',
        '<rect x="865" y="105" width="16" height="16" fill="#0f766e"/><text x="888" y="118" font-family="sans-serif" font-size="12">|D′|</text>',
    ])
    return _chart_footer(body, footer)


def _roh_svg(subtitle: str, lines: list[str], footer: str) -> str:
    """Trace la charge ROH médiane par portée et cohorte."""
    pattern = re.compile(r"^(.*?) / (.*?) \| evaluated=(\d+)/(\d+) \| median total=([0-9.]+) kb \| target in ROH=(\d+) \| (.+)$")
    rows = [(match.group(1), match.group(2), int(match.group(3)), int(match.group(4)), float(match.group(5)), int(match.group(6)), match.group(7)) for line in lines if (match := pattern.match(line))]
    if not rows:
        return _svg("ROH", subtitle, lines, footer)
    body = _chart_header("ROH et autozygotie", subtitle)
    body.append('<text x="55" y="112" font-family="sans-serif" font-size="15" font-weight="700">Charge ROH médiane (kb)</text>')
    maximum = max((total for *_, total, _, _ in rows), default=1.0) or 1.0
    for index, (scope, cohort, evaluated, count, total, target_count, status) in enumerate(rows):
        y = 165 + index * 125
        width = 650 * total / maximum
        color = "#0f766e" if "PRIMARY" in status else "#d97706"
        body.extend([
            f'<text x="55" y="{y}" font-family="sans-serif" font-size="13" font-weight="700">{html.escape(scope)} / {html.escape(cohort)}</text>',
            f'<rect x="300" y="{y - 21}" width="{width:.2f}" height="28" rx="4" fill="{color}"/>',
            f'<text x="{315 + width:.2f}" y="{y}" font-family="sans-serif" font-size="12">{total:g} kb</text>',
            f'<text x="55" y="{y + 25}" font-family="sans-serif" font-size="11">évalués {evaluated}/{count} · cible dans ROH {target_count} · {html.escape(status)}</text>',
        ])
    return _chart_footer(body, footer)


def _sensitivity_svg(subtitle: str, lines: list[str], footer: str) -> str:
    """Trace les écarts exploratoires par rapport à la référence primaire zéro."""
    pattern = re.compile(r"^EXPLORATORY (\S+) \| factor=(\S+) \| (\S+) \| (\S+) \| delta=([-0-9.]+|missing)$")
    rows = [(match.group(1), match.group(2), match.group(3), match.group(4), None if match.group(5) == "missing" else float(match.group(5))) for line in lines if (match := pattern.match(line))]
    if not rows:
        return _svg("SENSITIVITY", subtitle, lines, footer)
    body = _chart_header("Analyses de sensibilité", subtitle)
    body.extend([
        '<text x="55" y="112" font-family="sans-serif" font-size="15" font-weight="700">Variation relative des scénarios exploratoires</text>',
        '<line x1="600" y1="135" x2="600" y2="500" stroke="#0f766e" stroke-width="2"/>',
        '<text x="548" y="525" font-family="sans-serif" font-size="11">primaire = 0 %</text>',
    ])
    maximum = max((abs(delta) for *_, delta in rows if delta is not None), default=0.1) or 0.1
    for index, (domain, factor, evaluation, category, delta) in enumerate(rows):
        y = 170 + index * 75
        body.append(f'<text x="55" y="{y + 5}" font-family="sans-serif" font-size="13" font-weight="700">{html.escape(domain)}</text>')
        body.append(f'<circle cx="600" cy="{y}" r="6" fill="#0f766e"/>')
        if delta is None:
            body.append(f'<text x="655" y="{y + 5}" font-family="sans-serif" font-size="12" fill="#9a3412">NOT_EVALUATED · {html.escape(factor)}</text>')
            continue
        x = 600 + 300 * delta / maximum
        body.extend([
            f'<line x1="600" y1="{y}" x2="{x:.2f}" y2="{y}" stroke="#d97706" stroke-width="4"/>',
            f'<circle cx="{x:.2f}" cy="{y}" r="8" fill="#d97706"/>',
            f'<text x="{x + 12:.2f}" y="{y + 5}" font-family="sans-serif" font-size="11">{100 * delta:+.1f}% · {html.escape(factor)} · {html.escape(category)}</text>',
        ])
    body.append('<text x="55" y="575" font-family="sans-serif" font-size="12" fill="#7c2d12">Résultat primaire et scénarios exploratoires séparés · composite founder score: NOT CALCULATED</text>')
    return _chart_footer(body, footer)


DOMAIN_RENDERERS: dict[str, Callable[[str, list[str], str], str]] = {
    "POPULATION_STRUCTURE": _population_svg,
    "FOUNDER_IBS": _founder_svg,
    "VARIANT_AGE": _age_svg,
    "LOCAL_LD": _ld_svg,
    "ROH": _roh_svg,
    "SENSITIVITY": _sensitivity_svg,
}


def _write_result(output_dir: Path, run_id: str, figure_id: str, domain: str, status: str, lines: list[str], counts: tuple[int, int, int, int], source_artifacts: list[dict[str, Any]], warnings: list[str], limits: list[str]) -> FigureResult:
    figure_path = output_dir / f"{figure_id}.svg"
    renderer = DOMAIN_RENDERERS.get(domain)
    svg = renderer(f"run {run_id} — {status}", lines, warnings[0]) if renderer else _svg(domain.replace("_", " "), f"run {run_id} — {status}", lines, warnings[0])
    figure_path.write_text(svg, encoding="utf-8")
    provenance = {
        "schema_version": "1.0.0", "run_id": run_id, "figure_id": figure_id,
        "domain": domain, "status": status, "method_id": METHOD_ID,
        "script_version": SCRIPT_VERSION, "sources": [_source_record(item) for item in source_artifacts],
        "represented_count": counts[0], "excluded_count": counts[1],
        "missing_value_count": counts[2], "not_evaluated_count": counts[3],
        "pseudonymized": True, "sensitivity": "sensitive_genetic",
        "legend": ["PRIMARY = résultat primaire", "EXPLORATORY = scénario ou petit effectif"],
        "interpretation_warnings": warnings, "known_limits": limits,
    }
    validate_json_document(provenance, "figure_provenance.schema.json")
    provenance_path = output_dir / f"{figure_id}.figure.json"
    atomic_write_json(provenance_path, provenance)
    return FigureResult(figure_id, domain, status, figure_path, provenance_path, None)


def _population(paths: dict[str, Path]) -> tuple[str, list[str], tuple[int, int, int, int]]:
    scores = validate_tsv_table(paths["population_scores"], "population_scores.schema.json")
    eigenvalues = validate_tsv_table(paths["population_eigenvalues"], "population_eigenvalues.schema.json")
    outliers = validate_tsv_table(paths["population_outliers"], "population_outliers.schema.json")
    score_ids = {row["SAMPLE_ID"] for row in scores.rows}
    if score_ids != {row["SAMPLE_ID"] for row in outliers.rows}:
        raise VisualizationContractError("population_sample_identity_mismatch")
    outliers_by_id = {row["SAMPLE_ID"]: row for row in outliers.rows}
    reference_count = sum(row["REFERENCE_INCLUDED"] for row in scores.rows)
    if any(int(row["REFERENCE_SAMPLE_COUNT"]) != reference_count for row in eigenvalues.rows):
        raise VisualizationContractError("population_reference_count_mismatch")
    for row in scores.rows:
        outlier = outliers_by_id[row["SAMPLE_ID"]]
        if any(row[key] != outlier[key] for key in ("REFERENCE_INCLUDED", "PROJECTED", "OUTLIER_STATUS")):
            raise VisualizationContractError("population_projection_or_outlier_mismatch")
    ordered_scores = sorted(scores.rows, key=lambda row: row["SAMPLE_ID"])
    complete = [row for row in ordered_scores if row["PC1"] is not None and row["PC2"] is not None]
    missing = 2 * len(scores.rows) - sum(row[axis] is not None for row in scores.rows for axis in ("PC1", "PC2"))
    outlier_count = sum(row["OUTLIER_STATUS"] == "OUTLIER" for row in scores.rows)
    lines = [f"EIGEN|{row['COMPONENT']}|{row['EXPLAINED_VARIANCE_RATIO']}" for row in eigenvalues.rows]
    pseudonyms = {row["SAMPLE_ID"]: f"PCA-{index:03d}" for index, row in enumerate(ordered_scores, 1)}
    lines += [
        f"POINT|{pseudonyms[row['SAMPLE_ID']]}|{row['PC1']}|{row['PC2']}|{'REFERENCE' if row['REFERENCE_INCLUDED'] else 'PROJECTED'}|{row['OUTLIER_STATUS']}"
        for row in complete
    ]
    lines += [
        f"SUMMARY|reference n={reference_count}",
        f"SUMMARY|projected n={sum(row['PROJECTED'] for row in scores.rows)}",
        f"SUMMARY|outliers exploratoires={outlier_count}",
        f"SUMMARY|PC1/PC2 manquants={missing}",
    ]
    not_evaluated = int(not complete or len(eigenvalues.rows) < 2)
    return ("NOT_EVALUATED" if not_evaluated else "RENDERED", lines, (len(complete), 0, missing, not_evaluated))


def _founder(paths: dict[str, Path]) -> tuple[str, list[str], tuple[int, int, int, int]]:
    table = validate_tsv_table(paths["founder_segments"], "founder_segments.schema.json")
    summary = read_json(paths["founder_analysis_summary"]); validate_json_document(summary, "founder_analysis_summary.schema.json")
    included = [row for row in table.rows if row["SEGMENT_STATUS"] != "EXCLUDED"]
    excluded = [row for row in table.rows if row["SEGMENT_STATUS"] == "EXCLUDED"]
    if len(included) != summary["selected_carrier_count"] or len(excluded) != summary["excluded_carrier_count"]:
        raise VisualizationContractError("founder_effective_count_mismatch")
    missing = sum(row["PHASING_CONFIDENCE"] is None for row in table.rows)
    lines = [f"UNIT-{index:03d} | left={row['LEFT_LENGTH_CM']} cM | right={row['RIGHT_LENGTH_CM']} cM | {row['SEGMENT_STATUS']}" for index, row in enumerate(included, 1)]
    lines += [f"excluded units: {len(excluded)}", f"status: {summary['status']}", f"missing confidence: {missing}"]
    not_eval = int(summary["status"] != "SUPPORTED_IBS_CANDIDATE")
    return ("NOT_EVALUATED" if not_eval else "RENDERED", lines, (len(included), len(excluded), missing, not_eval))


def _age(paths: dict[str, Path]) -> tuple[str, list[str], tuple[int, int, int, int]]:
    estimates = validate_tsv_table(paths["variant_age_estimates"], "variant_age_estimates.schema.json")
    scenarios = validate_tsv_table(paths["variant_age_scenarios"], "variant_age_scenarios.schema.json")
    primary_rows = [row for row in estimates.rows if row["PRIMARY"]]
    if len(primary_rows) != 1 or primary_rows[0]["MODEL"] != "CORRELATED":
        raise VisualizationContractError("variant_age_primary_model_mismatch")
    evaluated_unit_counts = {
        row["N_UNITS"] for row in estimates.rows
        if row["ANALYSIS_STATUS"] != "NOT_ESTIMATED"
    }
    if len(evaluated_unit_counts) > 1:
        raise VisualizationContractError("variant_age_effective_count_mismatch")
    missing = sum(row["ESTIMATE_GENERATIONS"] is None for row in estimates.rows) + sum(row["ESTIMATE_GENERATIONS"] is None for row in scenarios.rows)
    not_eval = sum(row["ANALYSIS_STATUS"] == "NOT_ESTIMATED" for row in estimates.rows) + sum(row["STATUS"] == "NOT_ESTIMATED" for row in scenarios.rows)
    lines = [f"{'PRIMARY' if row['PRIMARY'] else 'EXPLORATORY'} {row['MODEL']} | n={row['N_UNITS']} | estimate={row['ESTIMATE_GENERATIONS'] or 'NOT_ESTIMATED'} generations | CI={row['CI_LOWER_GENERATIONS'] or 'missing'}..{row['CI_UPPER_GENERATIONS'] or 'missing'}" for row in estimates.rows]
    lines.append(f"exploratory scenarios: {len(scenarios.rows)} | not evaluated: {not_eval}")
    return ("NOT_EVALUATED" if not_eval == len(estimates.rows) + len(scenarios.rows) else "RENDERED", lines, (len(estimates.rows), 0, missing, not_eval))


def _ld(paths: dict[str, Path]) -> tuple[str, list[str], tuple[int, int, int, int]]:
    table = validate_tsv_table(paths["local_ld_summary"], "local_ld_summary.schema.json")
    missing = sum(row["MEDIAN_R2_GENOTYPE"] is None for row in table.rows) + sum(row["MEDIAN_D_PRIME_ABS"] is None for row in table.rows)
    not_eval = sum(row["SUMMARY_STATUS"] != "EVALUATED" for row in table.rows)
    lines = [f"{row['COHORT_ID']} | n={row['SAMPLE_COUNT']} | {row['DISTANCE_MIN_CM']}..{row['DISTANCE_MAX_CM']} cM | r2={row['MEDIAN_R2_GENOTYPE'] or 'missing'} | Dprime={row['MEDIAN_D_PRIME_ABS'] or 'missing'} | {row['SUMMARY_STATUS']}" for row in table.rows]
    return ("NOT_EVALUATED" if not_eval == len(table.rows) else "RENDERED", lines, (len(table.rows), 0, missing, not_eval))


def _roh(paths: dict[str, Path]) -> tuple[str, list[str], tuple[int, int, int, int]]:
    table = validate_tsv_table(paths["roh_cohort_summary"], "roh_cohort_summary.schema.json")
    if any(int(row["EVALUATED_SAMPLE_COUNT"]) > int(row["SAMPLE_COUNT"]) for row in table.rows):
        raise VisualizationContractError("roh_effective_count_mismatch")
    missing = sum(row[key] is None for row in table.rows for key in ("MEDIAN_N_ROH", "MEDIAN_TOTAL_ROH_KB", "MEDIAN_MAX_ROH_KB"))
    not_eval = sum(row["COHORT_STATUS"] == "NOT_EVALUATED" for row in table.rows)
    lines = [f"{row['SCOPE']} / {row['COHORT_ID']} | evaluated={row['EVALUATED_SAMPLE_COUNT']}/{row['SAMPLE_COUNT']} | median total={row['MEDIAN_TOTAL_ROH_KB'] or 'missing'} kb | target in ROH={row['TARGET_IN_ROH_COUNT']} | {row['COHORT_STATUS']}" for row in table.rows]
    return ("NOT_EVALUATED" if not_eval == len(table.rows) else "RENDERED", lines, (len(table.rows), 0, missing, not_eval))


def _sensitivity(paths: dict[str, Path]) -> tuple[str, list[str], tuple[int, int, int, int]]:
    comparisons = validate_tsv_table(paths["sensitivity_comparisons"], "sensitivity_comparisons.schema.json")
    stability = validate_tsv_table(paths["sensitivity_stability"], "sensitivity_stability.schema.json")
    primary = [row for row in comparisons.rows if row["ROLE"] == "PRIMARY"]
    if {row["DOMAIN"] for row in primary} != {row["DOMAIN"] for row in stability.rows}:
        raise VisualizationContractError("sensitivity_primary_domain_mismatch")
    for stable_row in stability.rows:
        scenario_rows = [
            row for row in comparisons.rows
            if row["DOMAIN"] == stable_row["DOMAIN"] and row["ROLE"] == "SENSITIVITY"
        ]
        if len(scenario_rows) != int(stable_row["EXPECTED_SCENARIO_COUNT"]):
            raise VisualizationContractError("sensitivity_effective_count_mismatch")
    missing = sum(row["SCENARIO_NUMERIC_VALUE"] is None for row in comparisons.rows)
    not_eval = sum(row["EVALUATION_STATUS"] == "NOT_EVALUATED" for row in comparisons.rows)
    lines = [f"PRIMARY {row['DOMAIN']} | {row['TECHNICAL_STATUS'] or 'NOT_EVALUATED'}" for row in primary]
    lines += [f"EXPLORATORY {row['DOMAIN']} | factor={row['CHANGED_FACTOR']} | {row['EVALUATION_STATUS']} | {row['CATEGORICAL_COMPARISON']} | delta={row['RELATIVE_CHANGE'] if row['RELATIVE_CHANGE'] is not None else 'missing'}" for row in comparisons.rows if row["ROLE"] == "SENSITIVITY"]
    lines.append("composite founder score: NOT CALCULATED")
    return ("NOT_EVALUATED" if not_eval == len(comparisons.rows) else "RENDERED", lines, (len(comparisons.rows), 0, missing, not_eval))


DOMAIN_SPECS: tuple[tuple[str, str, dict[str, tuple[str, str]], Callable[[dict[str, Path]], tuple[str, list[str], tuple[int, int, int, int]]], list[str], list[str]], ...] = (
    ("population_structure", "POPULATION_STRUCTURE", {"population_scores": ("analyze_population_structure", "population_scores.schema.json"), "population_eigenvalues": ("analyze_population_structure", "population_eigenvalues.schema.json"), "population_outliers": ("analyze_population_structure", "population_outliers.schema.json")}, _population, ["PCA exploratoire interne à l'étude ; aucune attribution automatique d'ascendance."], ["Les axes dépendent de la référence indépendante, des variants informatifs et de la structure présente dans l'étude.", "Les p-values d'outlier reposent sur une approximation exploratoire et aucune exclusion n'est automatique."]),
    ("founder_ibs", "FOUNDER_IBS", {"founder_segments": ("infer_founder_haplotype", "founder_segments.schema.json"), "founder_analysis_summary": ("infer_founder_haplotype", "founder_analysis_summary.schema.json")}, _founder, ["Partage IBS uniquement ; aucune preuve IBD ou causale."], ["Les limites dépendent de la phase, de la carte et de la densité des marqueurs."]),
    ("variant_age", "VARIANT_AGE", {"variant_age_estimates": ("estimate_variant_age", "variant_age_estimates.schema.json"), "variant_age_scenarios": ("estimate_variant_age", "variant_age_scenarios.schema.json")}, _age, ["Gamma conditionne sur une origine unique non démontrée par la figure."], ["Les petits effectifs et l'incertitude de phase/carte limitent l'estimation."]),
    ("local_ld", "LOCAL_LD", {"local_ld_summary": ("analyze_local_ld", "local_ld_summary.schema.json")}, _ld, ["LD descriptif secondaire ; aucune preuve d'origine fondatrice."], ["r2 et D-prime restent distincts et sensibles aux effectifs/fréquences."]),
    ("roh", "ROH", {"roh_cohort_summary": ("analyze_roh", "roh_cohort_summary.schema.json")}, _roh, ["Autozygotie individuelle distincte de l'IBS/IBD entre individus."], ["La densité ACPA et les paramètres de scan limitent les ROH observables."]),
    ("sensitivity", "SENSITIVITY", {"sensitivity_comparisons": ("run_sensitivity_analyses", "sensitivity_comparisons.schema.json"), "sensitivity_stability": ("run_sensitivity_analyses", "sensitivity_stability.schema.json")}, _sensitivity, ["Robustesse aux scénarios testés, pas validation externe ni causalité."], ["Un biais partagé par tous les runs n'est pas détecté."]),
)


def build_consolidated_figures(*, run_dir: Path, output_dir: Path, stage_inputs: dict[str, Any]) -> list[FigureResult]:
    """Construit les six figures sans jamais lire hors du manifeste d'entrée."""
    output_dir.mkdir(parents=True, exist_ok=True)
    artifacts = {item["artifact_id"]: item for item in stage_inputs["artifacts"]}
    results: list[FigureResult] = []
    for figure_id, domain, expected, builder, warnings, limits in DOMAIN_SPECS:
        source_items = [artifacts[item] for item in expected if item in artifacts]
        try:
            paths = _validate_sources(run_dir, artifacts, expected)
            status, lines, counts = builder(paths)
            results.append(_write_result(output_dir, stage_inputs["run_id"], figure_id, domain, status, lines, counts, source_items, warnings, limits))
        except (VisualizationContractError, ValueError, OSError) as error:
            reason = str(error) if isinstance(error, VisualizationContractError) else "source_contract_validation_failed"
            provenance = {"schema_version": "1.0.0", "run_id": stage_inputs["run_id"], "figure_id": figure_id, "domain": domain, "status": "BLOCKED", "method_id": METHOD_ID, "script_version": SCRIPT_VERSION, "sources": [_source_record(item) for item in source_items], "represented_count": 0, "excluded_count": 0, "missing_value_count": 0, "not_evaluated_count": 0, "pseudonymized": True, "sensitivity": "sensitive_genetic", "legend": ["BLOCKED = incohérence d'intégrité ou de contrat"], "interpretation_warnings": warnings, "known_limits": limits}
            provenance_path = output_dir / f"{figure_id}.figure.json"
            validate_json_document(provenance, "figure_provenance.schema.json"); atomic_write_json(provenance_path, provenance)
            results.append(FigureResult(figure_id, domain, "BLOCKED", None, provenance_path, reason))
    return results
