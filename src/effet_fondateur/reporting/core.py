"""Faits contrôlés, prompt prudent et rapport HTML révisable de l'étape 19."""

from __future__ import annotations

import html
import json
from dataclasses import dataclass
from pathlib import Path
from typing import Any

from effet_fondateur.audit import atomic_write_json, read_json, sha256_file
from effet_fondateur.contracts import load_pipeline_config, validate_json_document, validate_tsv_table


METHOD_ID = "reviewable_scientific_report_v1"
DOMAIN_TITLES = {
    "POPULATION_STRUCTURE": "Structure populationnelle — PCA",
    "FOUNDER_IBS": "Partage IBS autour du variant cible",
    "VARIANT_AGE": "Datation du variant",
    "LOCAL_LD": "LD local secondaire",
    "ROH": "ROH et autozygotie",
    "SENSITIVITY": "Analyses de sensibilité",
}
PARAMETER_STAGES = (
    "qc_preliminary", "build_kinship_panel", "infer_kinship",
    "analyze_population_structure", "prepare_target_region",
    "infer_founder_haplotype", "estimate_variant_age", "analyze_local_ld",
    "analyze_roh", "run_sensitivity_analyses",
)


@dataclass(frozen=True)
class ReportDraftPublication:
    facts_path: Path
    prompt_path: Path
    draft_path: Path
    html_path: Path
    editor_script_path: Path
    kinship_network_path: Path
    validation_path: Path


def _artifact_path(run_dir: Path, artifact: dict[str, Any]) -> Path:
    relative = Path(artifact["path"])
    if relative.is_absolute() or ".." in relative.parts or relative.parts[:2] == ("data", "output"):
        raise ValueError("report_source_not_confined_to_current_run")
    path = (run_dir / relative).resolve()
    try:
        path.relative_to(run_dir.resolve())
    except ValueError as error:
        raise ValueError("report_source_not_confined_to_current_run") from error
    if not path.is_file() or sha256_file(path) != artifact["sha256"]:
        raise ValueError("report_source_missing_or_checksum_mismatch")
    return path


def _artifact_map(stage_inputs: dict[str, Any]) -> dict[str, dict[str, Any]]:
    artifacts = {artifact["artifact_id"]: artifact for artifact in stage_inputs["artifacts"]}
    if len(artifacts) != len(stage_inputs["artifacts"]):
        raise ValueError("duplicate_report_artifact_id")
    return artifacts


def _study_context(run_dir: Path) -> dict[str, Any]:
    config_path = run_dir / "config.resolved.yaml"
    manifest = read_json(run_dir / "manifest.json")
    if sha256_file(config_path) != manifest["config_sha256"]:
        raise ValueError("resolved_config_checksum_mismatch")
    config = load_pipeline_config(config_path)
    parameters = {
        stage: config["stages"][stage]["parameters"]
        for stage in PARAMETER_STAGES
        if stage in config["stages"] and config["stages"][stage]["enabled"]
    }
    return {
        "project_id": config["project"]["project_id"],
        "assembly": config["project"]["assembly"],
        "target": config["target"],
        "parameters_by_stage": parameters,
    }


def _kinship_facts(run_dir: Path, artifacts: dict[str, dict[str, Any]]) -> dict[str, Any]:
    pairs_table = validate_tsv_table(_artifact_path(run_dir, artifacts["kinship_pairs"]), "kinship_pairs.schema.json")
    degree_table = validate_tsv_table(_artifact_path(run_dir, artifacts["kinship_degree_summary"]), "kinship_degree_summary.schema.json")
    report = read_json(_artifact_path(run_dir, artifacts["kinship_report"]))
    sample_ids = sorted({row[key] for row in pairs_table.rows for key in ("SAMPLE_ID_1", "SAMPLE_ID_2")})
    pseudonyms = {sample_id: f"KIN-{index:03d}" for index, sample_id in enumerate(sample_ids, 1)}
    pairs = [
        {
            "pair_id": f"PAIR-{index:04d}",
            "unit_1": pseudonyms[row["SAMPLE_ID_1"]],
            "unit_2": pseudonyms[row["SAMPLE_ID_2"]],
            "kinship": row["KINSHIP"], "ibs0": row["IBS0"],
            "inferred_degree": row["INFERRED_DEGREE"],
            "inferred_relation": row["INFERRED_RELATION"],
            "concordance_status": row["CONCORDANCE_STATUS"],
        }
        for index, row in enumerate(pairs_table.rows, 1)
    ]
    return {
        "status": "VALID" if report.get("sample_count", 0) else "NOT_EVALUATED",
        "sample_count": int(report.get("sample_count", 0)),
        "reported_pair_count": len(pairs),
        "related_pair_count": int(report.get("related_pair_count", 0)),
        "degree_counts": {row["INFERRED_DEGREE"]: int(row["PAIR_COUNT"]) for row in degree_table.rows},
        "pairs": pairs,
        "known_limits": [
            "KING estime l'apparentement genome-wide et ne prouve pas un IBD local.",
            "Les paires sous le degré demandé peuvent ne pas être matérialisées comme non apparentées.",
            "Aucune exclusion proposée n'est appliquée automatiquement.",
        ],
    }


def _section_facts(run_dir: Path, artifacts: dict[str, dict[str, Any]]) -> list[dict[str, Any]]:
    index = read_json(_artifact_path(run_dir, artifacts["figure_index"])); validate_json_document(index, "figure_index.schema.json")
    sections = []
    for figure in index["figures"]:
        provenance_id = f"figure_provenance_{figure['figure_id']}"
        provenance = read_json(_artifact_path(run_dir, artifacts[provenance_id])); validate_json_document(provenance, "figure_provenance.schema.json")
        status = {"RENDERED": "VALID", "NOT_EVALUATED": "NOT_EVALUATED", "BLOCKED": "BLOCKED"}[figure["status"]]
        sections.append({
            "section_id": figure["figure_id"], "title": DOMAIN_TITLES[figure["domain"]], "status": status,
            "facts": [
                f"Effectif représenté : {provenance['represented_count']}.",
                f"Exclusions : {provenance['excluded_count']}.",
                f"Valeurs manquantes : {provenance['missing_value_count']}.",
                f"Résultats non évalués : {provenance['not_evaluated_count']}.",
            ],
            "limits": provenance["known_limits"] + provenance["interpretation_warnings"],
            "source_artifact_ids": ["figure_index", provenance_id, f"figure_{figure['figure_id']}"] if figure["figure_path"] else ["figure_index", provenance_id],
        })
    return sections


def _build_facts(run_dir: Path, stage_inputs: dict[str, Any]) -> dict[str, Any]:
    artifacts = _artifact_map(stage_inputs)
    facts = {
        "schema_version": "1.0.0", "run_id": stage_inputs["run_id"], "method_id": METHOD_ID,
        "sensitivity": "sensitive_genetic", "scientific_recalculation_performed": False,
        "composite_founder_score_calculated": False, "study_context": _study_context(run_dir),
        "kinship": _kinship_facts(run_dir, artifacts), "sections": _section_facts(run_dir, artifacts),
    }
    validate_json_document(facts, "interpretation_facts.schema.json")
    return facts


def _prompt(facts: dict[str, Any]) -> str:
    return """Tu rédiges une proposition de commentaire scientifique uniquement à partir du JSON fourni.
N'ajoute aucun calcul, nombre, unité, individu ou hypothèse. Sépare observation, interprétation prudente et limites.
N'emploie jamais prouve, démontre, confirme l'effet fondateur, origine unique ou causalité.
NOT_EVALUATED interdit toute conclusion. Distingue primaire et exploratoire. KING n'est pas une preuve IBD locale ;
PCA/DAPC n'attribuent pas automatiquement une ascendance ; ROH, LD, IBS et datation restent des arguments distincts.
Ne calcule aucun score composite. Retourne exclusivement un JSON conforme à interpretation_draft.schema.json.

FAITS_CONTROLES:
""" + json.dumps(facts, ensure_ascii=False, sort_keys=True, indent=2) + "\n"


def _deterministic_draft(facts: dict[str, Any], facts_sha: str, prompt_sha: str) -> dict[str, Any]:
    sections = []
    kinship = facts["kinship"]
    sections.append({
        "section_id": "kinship", "observation": f"KING publie {kinship['reported_pair_count']} paire(s), dont {kinship['related_pair_count']} classée(s) apparentée(s), sur {kinship['sample_count']} échantillon(s).",
        "interpretation_prudente": "Ces résultats décrivent l'apparentement genome-wide au seuil configuré et nécessitent une revue humaine.",
        "limites": " ".join(kinship["known_limits"]), "review_required": True,
    })
    for section in facts["sections"]:
        observation = " ".join(section["facts"])
        interpretation = "Aucune conclusion scientifique n'est proposée." if section["status"] != "VALID" else "Ce résultat constitue une observation descriptive à interpréter avec les autres domaines sans synthèse causale."
        sections.append({"section_id": section["section_id"], "observation": observation, "interpretation_prudente": interpretation, "limites": " ".join(section["limits"]), "review_required": True})
    return {
        "schema_version": "1.0.0", "run_id": facts["run_id"], "method_id": METHOD_ID,
        "draft_status": "DETERMINISTIC_DRAFT", "generator": {"provider": "disabled", "model": None, "network_used": False},
        "facts_sha256": facts_sha, "prompt_sha256": prompt_sha, "sections": sections,
        "human_validation_required": True,
    }


def _network_svg(kinship: dict[str, Any]) -> str:
    pairs = [pair for pair in kinship["pairs"] if pair["inferred_degree"] != "UNRELATED"]
    nodes = sorted({pair[key] for pair in pairs for key in ("unit_1", "unit_2")})
    width, height = 1000, max(300, 90 * len(nodes) + 100)
    positions = {node: (180 + 300 * (index % 3), 90 + 90 * (index // 3)) for index, node in enumerate(nodes)}
    colors = {"FIRST": "#dc2626", "SECOND": "#d97706", "THIRD": "#2563eb", "FOURTH": "#64748b", "DUPLICATE_OR_MZ": "#7c3aed"}
    body = [f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}">', '<rect width="100%" height="100%" fill="#fbfaf7"/>']
    for pair in pairs:
        x1, y1 = positions[pair["unit_1"]]; x2, y2 = positions[pair["unit_2"]]
        body.append(f'<line x1="{x1}" y1="{y1}" x2="{x2}" y2="{y2}" stroke="{colors.get(pair["inferred_degree"], "#64748b")}" stroke-width="3"/>')
    for node, (x, y) in positions.items():
        body.extend([f'<circle cx="{x}" cy="{y}" r="25" fill="#0f766e"/>', f'<text x="{x}" y="{y + 4}" text-anchor="middle" font-family="monospace" font-size="11" fill="white">{node}</text>'])
    body.append('<text x="30" y="30" font-family="sans-serif" font-size="14" fill="#7c2d12">Réseau KING pseudonymisé — apparentement genome-wide, pas preuve IBD locale</text></svg>\n')
    return "".join(body)


def _editor_script() -> str:
    return """'use strict';
document.getElementById('approve-report').addEventListener('click', () => {
  const reviewer = document.getElementById('reviewer').value.trim();
  const checks = [...document.querySelectorAll('.approval-check')];
  if (!reviewer || checks.some(item => !item.checked)) { alert('Renseignez le relecteur et validez toute la checklist.'); return; }
  const body = document.body.dataset;
  const sections = [...document.querySelectorAll('textarea[data-section]')].map(item => ({section_id:item.dataset.section, approved_text:item.value.trim()}));
  if (sections.some(item => !item.approved_text)) { alert('Aucune section ne peut être vide.'); return; }
  const review = {schema_version:'1.0.0', run_id:body.runId, review_status:'APPROVED', reviewer,
    reviewed_at:new Date().toISOString(), facts_sha256:body.factsSha, draft_sha256:body.draftSha, sections,
    checklist:{facts_checked:true, limits_checked:true, non_causal_language_checked:true, primary_exploratory_checked:true}};
  const link=document.createElement('a'); link.href=URL.createObjectURL(new Blob([JSON.stringify(review,null,2)],{type:'application/json'}));
  link.download='report_review.json'; link.click(); URL.revokeObjectURL(link.href);
  document.getElementById('approval-result').textContent='Décision téléchargée. Faites-la valider par le finaliseur V2 pour produire le HTML verrouillé puis le PDF.';
});
"""


def _render_html(facts: dict[str, Any], draft: dict[str, Any], artifacts: dict[str, dict[str, Any]], run_dir: Path, facts_sha: str, draft_sha: str) -> str:
    esc = lambda value: html.escape(str(value), quote=True)
    context = facts["study_context"]
    sections_by_id = {item["section_id"]: item for item in draft["sections"]}
    styles = "body{margin:0;background:#f5f7fb;color:#172554;font:16px system-ui}header{padding:3rem 6%;color:white;background:linear-gradient(120deg,#172554,#0f766e)}main{max-width:1180px;margin:auto;padding:2rem}.card{background:white;border:1px solid #dbe4ee;border-radius:14px;padding:1.4rem;margin:1.4rem 0;box-shadow:0 8px 25px #0f172a10}table{width:100%;border-collapse:collapse}th,td{padding:.55rem;border:1px solid #cbd5e1;text-align:left}th{background:#e6fffb}textarea{width:100%;min-height:190px;padding:1rem;border:2px solid #94a3b8;border-radius:8px;font:15px system-ui;box-sizing:border-box}.facts{background:#f8fafc;padding:1rem;border-left:5px solid #0f766e}.limit{color:#9a3412}.badge{padding:.25rem .7rem;border-radius:2rem;background:#fef3c7}img.figure{width:100%}button{background:#0f766e;color:white;border:0;padding:1rem 1.5rem;border-radius:8px;font-weight:700;cursor:pointer}code{overflow-wrap:anywhere}"
    target = context["target"]
    content = ["<!doctype html><html lang=\"fr\"><head><meta charset=\"utf-8\"><meta name=\"robots\" content=\"noindex,nofollow\"><meta http-equiv=\"Content-Security-Policy\" content=\"default-src 'none'; img-src 'self' data:; style-src 'unsafe-inline'; script-src 'self'\"><title>Rapport V2 révisable</title>", f"<style>{styles}</style></head><body data-run-id=\"{esc(facts['run_id'])}\" data-facts-sha=\"{facts_sha}\" data-draft-sha=\"{draft_sha}\"><header><h1>Rapport V2 révisable</h1><p>Brouillon non validé · données pseudonymisées · sensitive_genetic</p></header><main>"]
    content.append(f'<section class="card"><h2>Paramètres d’entrée de l’étude</h2><table><tr><th>Run</th><td>{esc(facts["run_id"])}</td></tr><tr><th>Projet</th><td>{esc(context["project_id"])}</td></tr><tr><th>Assemblage</th><td>{esc(context["assembly"])}</td></tr><tr><th>Cible</th><td>{esc(target.get("gene"))} · chr{esc(target.get("chromosome"))}:{esc(target.get("position_bp"))} · {esc(target.get("project_variant_id"))}</td></tr></table>')
    for stage, parameters in context["parameters_by_stage"].items():
        content.append(f'<details><summary>{esc(stage)}</summary><pre>{esc(json.dumps(parameters, ensure_ascii=False, indent=2, sort_keys=True))}</pre></details>')
    content.append('</section>')
    kin = facts["kinship"]; kin_draft = sections_by_id["kinship"]
    content.append(f'<section class="card"><h2>Parenté KING <span class="badge">{esc(kin["status"])}</span></h2><p class="facts">n={kin["sample_count"]} · paires publiées={kin["reported_pair_count"]} · apparentées={kin["related_pair_count"]}</p><img class="figure" src="kinship_network.svg" alt="Réseau KING pseudonymisé"><details><summary>Tableau des paires KING pseudonymisées</summary><table><tr><th>Paire</th><th>Unités</th><th>Kinship</th><th>IBS0</th><th>Relation</th><th>Concordance</th></tr>')
    for pair in kin["pairs"]:
        content.append(f'<tr><td>{esc(pair["pair_id"])}</td><td>{esc(pair["unit_1"])} / {esc(pair["unit_2"])}</td><td>{esc(pair["kinship"])}</td><td>{esc(pair["ibs0"])}</td><td>{esc(pair["inferred_relation"])}</td><td>{esc(pair["concordance_status"])}</td></tr>')
    text = "\n\n".join((kin_draft["observation"], kin_draft["interpretation_prudente"], kin_draft["limites"]))
    content.append(f'</table></details><h3>Commentaire à relire</h3><textarea data-section="kinship">{esc(text)}</textarea></section>')
    for section in facts["sections"]:
        section_draft = sections_by_id[section["section_id"]]
        figure_id = f"figure_{section['section_id']}"
        figure_html = ""
        if figure_id in artifacts:
            svg = _artifact_path(run_dir, artifacts[figure_id]).read_text(encoding="utf-8")
            figure_html = f'<div class="figure">{svg}</div>'
        facts_html = "".join(f"<li>{esc(item)}</li>" for item in section["facts"])
        limits_html = "".join(f"<li>{esc(item)}</li>" for item in section["limits"])
        text = "\n\n".join((section_draft["observation"], section_draft["interpretation_prudente"], section_draft["limites"]))
        content.append(f'<section class="card"><h2>{esc(section["title"])} <span class="badge">{esc(section["status"])}</span></h2>{figure_html}<div class="facts"><strong>Faits contrôlés</strong><ul>{facts_html}</ul><strong>Limites</strong><ul class="limit">{limits_html}</ul></div><h3>Commentaire à relire</h3><textarea data-section="{esc(section["section_id"])}">{esc(text)}</textarea></section>')
    content.append('<section class="card" id="human-validation"><h2>Validation humaine</h2><label>Relecteur ou rôle <input id="reviewer" required></label><p><label><input class="approval-check" type="checkbox"> Faits et nombres contrôlés</label><br><label><input class="approval-check" type="checkbox"> Limites conservées</label><br><label><input class="approval-check" type="checkbox"> Aucun langage causal ou de preuve</label><br><label><input class="approval-check" type="checkbox"> Primaire et exploratoire distingués</label></p><button id="approve-report" type="button">Valider le rapport</button><p id="approval-result"></p></section><p>La validation éditoriale ne démontre pas un effet fondateur.</p></main><script src="report_editor.js"></script></body></html>')
    return "".join(content)


def build_report_draft(*, run_dir: Path, output_dir: Path, stage_inputs: dict[str, Any]) -> ReportDraftPublication:
    """Produit le paquet de faits et le rapport local révisable sans appel réseau."""
    output_dir.mkdir(parents=True, exist_ok=True)
    artifacts = _artifact_map(stage_inputs)
    facts = _build_facts(run_dir, stage_inputs)
    facts_path = output_dir / "interpretation_facts.json"; atomic_write_json(facts_path, facts)
    prompt_path = output_dir / "interpretation_prompt.txt"; prompt_path.write_text(_prompt(facts), encoding="utf-8")
    draft = _deterministic_draft(facts, sha256_file(facts_path), sha256_file(prompt_path))
    validate_json_document(draft, "interpretation_draft.schema.json")
    draft_path = output_dir / "interpretation_draft.json"; atomic_write_json(draft_path, draft)
    network_path = output_dir / "kinship_network.svg"; network_path.write_text(_network_svg(facts["kinship"]), encoding="utf-8")
    editor_path = output_dir / "report_editor.js"; editor_path.write_text(_editor_script(), encoding="utf-8")
    html_path = output_dir / "report_draft.html"; html_path.write_text(_render_html(facts, draft, artifacts, run_dir, sha256_file(facts_path), sha256_file(draft_path)), encoding="utf-8")
    blocking = [section["section_id"] for section in facts["sections"] if section["status"] == "BLOCKED"]
    validation = {"schema_version": "1.0.0", "run_id": facts["run_id"], "method_id": METHOD_ID, "status": "TECHNICAL_INCOMPLETE" if blocking else "AWAITING_HUMAN_REVIEW", "facts_sha256": sha256_file(facts_path), "draft_sha256": sha256_file(draft_path), "human_validation_required": True, "checks": [{"check": "no_external_ai_call", "status": "PASS"}, {"check": "pseudonymized_kinship", "status": "PASS"}, {"check": "no_scientific_recalculation", "status": "PASS"}], "blocking_reasons": [f"blocked_figure:{item}" for item in blocking]}
    validate_json_document(validation, "report_validation.schema.json")
    validation_path = output_dir / "report_validation.json"; atomic_write_json(validation_path, validation)
    return ReportDraftPublication(facts_path, prompt_path, draft_path, html_path, editor_path, network_path, validation_path)
