from __future__ import annotations

import json
from pathlib import Path
from typing import Any

import pytest
import yaml

from effet_fondateur.audit import atomic_write_json, sha256_file
from effet_fondateur.contracts import build_file_artifact, validate_json_document
from effet_fondateur.reporting import finalize_report
from effet_fondateur.stages.build_report import execute as execute_report
from effet_fondateur.stages.build_visualizations import execute as execute_visualizations
from tests.test_v2_visualizations import _fixture, _producer_controls, _tsv


KINSHIP_SIGNATURE = "07" * 32


def _kinship_outputs(run_dir: Path) -> list[dict[str, Any]]:
    output_dir = run_dir / "stages" / "07_infer_kinship"; output_dir.mkdir(parents=True)
    pair_path = output_dir / "kinship_pairs.tsv"
    _tsv(pair_path, ["PAIR_ID", "SAMPLE_ID_1", "SAMPLE_ID_2", "FID_1", "IID_1", "FID_2", "IID_2", "N_SNP", "HETHET", "IBS0", "KINSHIP", "INFERRED_DEGREE", "INFERRED_RELATION", "DECLARED_RELATION", "PAIR_ORIGIN", "CONCORDANCE_STATUS"], [
        ["kinship_pair_000001", "PRIVATE_KIN_A", "PRIVATE_KIN_B", "fam_a", "iid_a", "fam_b", "iid_b", 1000, 0.1, 0.001, 0.2, "FIRST", "PARENT_CHILD_COMPATIBLE", "PARENT_CHILD", "BETWEEN_FAMILY", "CONCORDANT"],
        ["kinship_pair_000002", "PRIVATE_KIN_A", "PRIVATE_KIN_C", "fam_a", "iid_a", "fam_c", "iid_c", 1000, 0.1, 0.02, 0.05, "THIRD", "THIRD_DEGREE", "NONE", "BETWEEN_FAMILY", "CRYPTIC_RELATEDNESS"],
    ])
    degree_path = output_dir / "kinship_degree_summary.tsv"
    _tsv(degree_path, ["INFERRED_DEGREE", "PAIR_COUNT", "DECLARED_PAIR_COUNT", "CRYPTIC_PAIR_COUNT", "CONCORDANT_PAIR_COUNT", "DISCORDANT_PAIR_COUNT"], [
        ["FIRST", 1, 1, 0, 1, 0], ["THIRD", 1, 0, 1, 0, 0],
    ])
    report_path = output_dir / "kinship_report.json"
    atomic_write_json(report_path, {"parameters": {"king_degree": 4}, "sample_count": 3, "reported_pair_count": 2, "related_pair_count": 2, "cryptic_pair_count": 1, "manual_validation_required": True})
    specs = [
        ("kinship_pairs", pair_path, "kinship_pairs.schema.json", "sensitive_genetic"),
        ("kinship_degree_summary", degree_path, "kinship_degree_summary.schema.json", "internal"),
        ("kinship_report", report_path, None, "internal"),
    ]
    artifacts = [build_file_artifact(
        physical_path=path, published_path=path.relative_to(run_dir).as_posix(), artifact_id=artifact_id,
        artifact_type=artifact_id, media_type="application/json" if path.suffix == ".json" else "text/tab-separated-values",
        producer_stage="infer_kinship", producer_signature=KINSHIP_SIGNATURE,
        schema_name=schema, schema_version="1.0.0" if schema else None, assembly="GRCh38",
        sample_set_id="synthetic_samples", variant_set_id=None, sensitivity=sensitivity,
    ) for artifact_id, path, schema, sensitivity in specs]
    outputs = {"schema_version": "1.0.0", "run_id": "synthetic_visual_run", "stage_id": "07", "stage_name": "infer_kinship", "signature": KINSHIP_SIGNATURE, "artifacts": artifacts}
    atomic_write_json(output_dir / "stage_outputs.json", outputs)
    audit = {"schema_version": "1.0.0", "run_id": "synthetic_visual_run", "stage_id": "07", "stage_name": "infer_kinship", "method_id": "synthetic_king_fixture_v1", "signature": KINSHIP_SIGNATURE, "started_at": "2026-08-07T00:00:00Z", "completed_at": "2026-08-07T00:00:00Z", "duration_seconds": 0.0, "inputs": [], "outputs": artifacts, "parameters": {}, "tools": [], "counts": {}, "metrics": {}, "exclusions": [], "warnings": [], "checks": [], "known_limits": [], "expected_visualizations": [], "manual_validation_required": True}
    atomic_write_json(output_dir / "audit.json", audit)
    return artifacts


def test_stage_19_builds_editable_pseudonymized_report_without_ai_call(tmp_path: Path) -> None:
    run_dir = tmp_path / "run"; visualization_dir = run_dir / "stages" / "18_build_visualizations"; visualization_dir.mkdir(parents=True)
    stage_18_inputs = _fixture(run_dir); _producer_controls(run_dir, stage_18_inputs)
    atomic_write_json(visualization_dir / "stage_inputs.json", stage_18_inputs)
    assert execute_visualizations(visualization_dir / "stage_inputs.json", visualization_dir) == 0
    kinship_artifacts = _kinship_outputs(run_dir)
    config_path = run_dir / "config.resolved.yaml"
    config = yaml.safe_load(Path("config/pipeline.example.yaml").read_text(encoding="utf-8"))
    for stage_name in ("qc_preliminary", "build_kinship_panel", "infer_kinship", "analyze_population_structure", "prepare_target_region", "infer_founder_haplotype", "estimate_variant_age", "analyze_local_ld", "analyze_roh", "run_sensitivity_analyses"):
        config["stages"][stage_name]["enabled"] = True
    config_path.write_text(yaml.safe_dump(config, sort_keys=False), encoding="utf-8")
    manifest = json.loads((run_dir / "manifest.json").read_text(encoding="utf-8"))
    stage_18_outputs = json.loads((visualization_dir / "stage_outputs.json").read_text(encoding="utf-8"))
    manifest["config_sha256"] = sha256_file(config_path)
    manifest["stages"] += [
        {"stage_id": "07", "stage_name": "infer_kinship", "state": "SUCCEEDED", "signature": KINSHIP_SIGNATURE, "audit_path": "stages/07_infer_kinship/audit.json", "audit_sha256": sha256_file(run_dir / "stages/07_infer_kinship/audit.json"), "stage_outputs_sha256": sha256_file(run_dir / "stages/07_infer_kinship/stage_outputs.json")},
        {"stage_id": "18", "stage_name": "build_visualizations", "state": "SUCCEEDED", "signature": stage_18_inputs["signature"], "audit_path": "stages/18_build_visualizations/audit.json", "audit_sha256": sha256_file(visualization_dir / "audit.json"), "stage_outputs_sha256": sha256_file(visualization_dir / "stage_outputs.json")},
    ]
    atomic_write_json(run_dir / "manifest.json", manifest)
    required_ids = {
        "figure_index", "visualization_completeness", "visualization_render_manifest",
        *[f"figure_{name}" for name in ("population_structure", "founder_ibs", "variant_age", "local_ld", "roh", "sensitivity")],
        *[f"figure_provenance_{name}" for name in ("population_structure", "founder_ibs", "variant_age", "local_ld", "roh", "sensitivity")],
    }
    inputs = kinship_artifacts + [item for item in stage_18_outputs["artifacts"] if item["artifact_id"] in required_ids]
    report_dir = run_dir / "stages" / "19_build_report"; report_dir.mkdir()
    stage_inputs = {"schema_version": "1.0.0", "run_id": "synthetic_visual_run", "stage_id": "19", "stage_name": "build_report", "signature": "19" * 32, "attempt_number": 1, "published_output_dir": "stages/19_build_report", "parameters": {"method": "reviewable_scientific_report_v1", "ai_provider": "disabled"}, "artifacts": inputs}
    atomic_write_json(report_dir / "stage_inputs.json", stage_inputs)
    assert execute_report(report_dir / "stage_inputs.json", report_dir) == 0
    facts = json.loads((report_dir / "interpretation_facts.json").read_text(encoding="utf-8"))
    draft = json.loads((report_dir / "interpretation_draft.json").read_text(encoding="utf-8"))
    validation = json.loads((report_dir / "report_validation.json").read_text(encoding="utf-8"))
    validate_json_document(facts, "interpretation_facts.schema.json")
    validate_json_document(draft, "interpretation_draft.schema.json")
    validate_json_document(validation, "report_validation.schema.json")
    combined = "".join(path.read_text(encoding="utf-8") for path in (report_dir / "interpretation_facts.json", report_dir / "interpretation_prompt.txt", report_dir / "report_draft.html"))
    assert "PRIVATE_KIN" not in combined and "KIN-001" in combined
    assert "Paramètres d’entrée de l’étude" in combined and "Tableau des paires KING" in combined
    assert "approve-report" in combined and "report_review.json" in (report_dir / "report_editor.js").read_text(encoding="utf-8")
    assert draft["draft_status"] == "DETERMINISTIC_DRAFT"
    assert draft["generator"]["network_used"] is False
    assert validation["status"] == "AWAITING_HUMAN_REVIEW"
    assert facts["scientific_recalculation_performed"] is False
    assert facts["composite_founder_score_calculated"] is False
    review = {
        "schema_version": "1.0.0", "run_id": facts["run_id"], "review_status": "APPROVED",
        "reviewer": "synthetic-reviewer", "reviewed_at": "2026-08-07T12:00:00Z",
        "facts_sha256": sha256_file(report_dir / "interpretation_facts.json"),
        "draft_sha256": sha256_file(report_dir / "interpretation_draft.json"),
        "sections": [
            {"section_id": section["section_id"], "approved_text": "\n\n".join((section["observation"], section["interpretation_prudente"], section["limites"]))}
            for section in draft["sections"]
        ],
        "checklist": {"facts_checked": True, "limits_checked": True, "non_causal_language_checked": True, "primary_exploratory_checked": True},
    }
    review_path = report_dir / "report_review.json"; atomic_write_json(review_path, review)
    final_html, final_pdf, final_validation = finalize_report(report_dir, review_path)
    assert final_html.is_file() and "<textarea" not in final_html.read_text(encoding="utf-8")
    assert final_pdf.read_bytes().startswith(b"%PDF-")
    assert json.loads(final_validation.read_text(encoding="utf-8"))["status"] == "READY_FOR_SCIENTIFIC_REVIEW"
    review["sections"][0]["approved_text"] = "Cette analyse prouve un effet fondateur."
    atomic_write_json(review_path, review)
    with pytest.raises(ValueError, match="forbidden_assertive_claim:kinship"):
        finalize_report(report_dir, review_path)
