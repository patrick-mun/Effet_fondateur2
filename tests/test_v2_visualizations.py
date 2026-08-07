from __future__ import annotations

import csv
import json
from pathlib import Path
from typing import Any

from effet_fondateur.audit import atomic_write_json, sha256_file
from effet_fondateur.contracts import build_file_artifact, validate_json_document
from effet_fondateur.stages.build_visualizations import execute
from effet_fondateur.visualization import build_consolidated_figures


SIGNATURES = {
    "analyze_population_structure": "08" * 32,
    "infer_founder_haplotype": "13" * 32,
    "estimate_variant_age": "14" * 32,
    "analyze_local_ld": "15" * 32,
    "analyze_roh": "16" * 32,
    "run_sensitivity_analyses": "17" * 32,
}


def _tsv(path: Path, columns: list[str], rows: list[list[Any]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(columns)
        writer.writerows(rows)


def _fixture(run_dir: Path, *, founder_count_mismatch: bool = False, ld_not_evaluated: bool = False, pca_reference_count_mismatch: bool = False) -> dict[str, Any]:
    definitions: list[tuple[str, str, str, list[str] | None, list[list[Any]] | dict[str, Any]]] = []
    population_score_columns = ["SAMPLE_ID", "SAMPLE_SET_ID", "FID", "IID", "FAMILY_ID", "GROUP_LABEL", "ARRAY_BATCH", "REFERENCE_INCLUDED", "PROJECTED", "PC1", "PC2", "PC3", "PC4", "PC5", "PC6", "PC7", "PC8", "PC9", "PC10", "OUTLIER_STATUS", "OUTLIER_REASON"]
    definitions.append(("population_scores", "analyze_population_structure", "population_scores.schema.json", population_score_columns, [
        ["PRIVATE_PCA_A", "synthetic_samples", "fid_a", "iid_a", "fam_a", "group_a", "batch_a", "true", "false", -1.0, 0.5, "", "", "", "", "", "", "", "", "INLIER", ""],
        ["PRIVATE_PCA_B", "synthetic_samples", "fid_b", "iid_b", "fam_b", "group_b", "batch_a", "true", "false", 0.0, -0.5, "", "", "", "", "", "", "", "", "INLIER", ""],
        ["PRIVATE_PCA_C", "synthetic_samples", "fid_c", "iid_c", "fam_c", "group_c", "batch_b", "false", "true", 1.5, 1.2, "", "", "", "", "", "", "", "", "OUTLIER", "MULTIVARIATE_PC_DISTANCE"],
    ]))
    definitions.append(("population_eigenvalues", "analyze_population_structure", "population_eigenvalues.schema.json",
        ["COMPONENT", "EIGENVALUE", "EXPLAINED_VARIANCE_RATIO", "CUMULATIVE_EXPLAINED_VARIANCE_RATIO", "SINGULAR_VALUE", "REFERENCE_SAMPLE_COUNT", "INFORMATIVE_VARIANT_COUNT"],
        [["PC1", 2.0, 0.6, 0.6, 1.414, 3 if pca_reference_count_mismatch else 2, 100], ["PC2", 1.0, 0.3, 0.9, 1.0, 3 if pca_reference_count_mismatch else 2, 100]]))
    definitions.append(("population_outliers", "analyze_population_structure", "population_outliers.schema.json",
        ["SAMPLE_ID", "SAMPLE_SET_ID", "REFERENCE_INCLUDED", "PROJECTED", "COMPONENTS_USED", "DISTANCE_SQUARED", "P_VALUE", "OUTLIER_ALPHA", "OUTLIER_STATUS", "OUTLIER_REASON", "PROPOSED_EXCLUSION", "MANUAL_VALIDATION_REQUIRED"],
        [["PRIVATE_PCA_A", "synthetic_samples", "true", "false", 2, 0.5, 0.7, 0.05, "INLIER", "", "false", "false"], ["PRIVATE_PCA_B", "synthetic_samples", "true", "false", 2, 0.8, 0.6, 0.05, "INLIER", "", "false", "false"], ["PRIVATE_PCA_C", "synthetic_samples", "false", "true", 2, 9.0, 0.01, 0.05, "OUTLIER", "MULTIVARIATE_PC_DISTANCE", "true", "true"]]))
    definitions.append(("founder_segments", "infer_founder_haplotype", "founder_segments.schema.json",
        ["INDEPENDENT_UNIT_ID", "SAMPLE_ID", "FAMILY_ID", "CARRIER_HAPLOTYPE_ID", "TARGET_VARIANT_ID", "LEFT_BOUND_BP", "TARGET_BP", "RIGHT_BOUND_BP", "LEFT_LENGTH_CM", "RIGHT_LENGTH_CM", "PHASING_CONFIDENCE", "SEGMENT_METHOD", "SEGMENT_STATUS", "EXCLUSION_CODE"],
        [["family_alpha", "PRIVATE_SAMPLE_A", "fam_a", "H1", "target_v1", 90, 100, 120, 0.1, 0.2, "", "target_centered_exact_ibs_v1", "INCLUDED", ""], ["family_beta", "PRIVATE_SAMPLE_B", "fam_b", "H2", "target_v1", 80, 100, 130, 0.2, 0.3, 0.95, "target_centered_exact_ibs_v1", "INCLUDED", ""]]))
    definitions.append(("founder_analysis_summary", "infer_founder_haplotype", "founder_analysis_summary.schema.json", None,
        {"schema_version": "1.0.0", "method_id": "target_centered_exact_ibs_v1", "interpretation": "IBS_SHARED_CANDIDATE", "status": "SUPPORTED_IBS_CANDIDATE", "selected_carrier_count": 3 if founder_count_mismatch else 2, "excluded_carrier_count": 0, "background_haplotype_count": 8, "matching_background_haplotype_count": 0, "minimum_independent_carriers": 2, "minimum_flank_markers": 1, "ibd_claimed": False}))
    definitions.append(("variant_age_estimates", "estimate_variant_age", "variant_age_estimates.schema.json",
        ["METHOD_ID", "MODEL", "ANALYSIS_STATUS", "N_UNITS", "EFFECTIVE_N", "RHO", "ESTIMATE_GENERATIONS", "CI_LOWER_GENERATIONS", "CI_UPPER_GENERATIONS", "CONFIDENCE_LEVEL", "PRIMARY", "EXCLUSION_CODE"],
        [["gamma_gandolfo_2014_v1", "CORRELATED", "EXPLORATORY", 3, 2.5, 0.1, 10, 5, 20, 0.95, "true", ""], ["gamma_gandolfo_2014_v1", "INDEPENDENT", "EXPLORATORY", 3, 3, 0, 9, 4, 18, 0.95, "false", ""]]))
    definitions.append(("variant_age_scenarios", "estimate_variant_age", "variant_age_scenarios.schema.json",
        ["SCENARIO_ID", "SCENARIO_TYPE", "MODEL", "OMITTED_INDEPENDENT_UNIT_ID", "N_UNITS", "ESTIMATE_GENERATIONS", "CI_LOWER_GENERATIONS", "CI_UPPER_GENERATIONS", "GENERATION_YEARS", "ESTIMATE_YEARS", "CI_LOWER_YEARS", "CI_UPPER_YEARS", "STATUS", "EXCLUSION_CODE"],
        [["model_independent", "MODEL", "INDEPENDENT", "", 3, 9, 4, 18, "", "", "", "", "ESTIMATED", ""]]))
    ld_status = "COHORT_TOO_SMALL" if ld_not_evaluated else "EVALUATED"
    definitions.append(("local_ld_summary", "analyze_local_ld", "local_ld_summary.schema.json",
        ["COHORT_ID", "DISTANCE_BIN_ID", "DISTANCE_MIN_CM", "DISTANCE_MAX_CM", "SAMPLE_COUNT", "COHORT_STATUS", "TOTAL_PAIR_COUNT", "EVALUATED_PAIR_COUNT", "TARGET_PAIR_COUNT", "MEDIAN_R2_GENOTYPE", "MEDIAN_D_PRIME_ABS", "SUMMARY_STATUS"],
        [["controls_unrelated", "bin_1", 0, 0.5, 4 if ld_not_evaluated else 20, "NOT_EVALUATED" if ld_not_evaluated else "DESCRIPTIVE_PRIMARY", 12, 0 if ld_not_evaluated else 12, 2, "" if ld_not_evaluated else 0.2, "" if ld_not_evaluated else 0.4, ld_status]]))
    definitions.append(("roh_cohort_summary", "analyze_roh", "roh_cohort_summary.schema.json",
        ["SCOPE", "COHORT_ID", "SAMPLE_COUNT", "EVALUATED_SAMPLE_COUNT", "MEDIAN_N_ROH", "MEDIAN_TOTAL_ROH_KB", "MEDIAN_MAX_ROH_KB", "TARGET_IN_ROH_COUNT", "COHORT_STATUS"],
        [["GENOMEWIDE_BURDEN", "controls_unrelated", 20, 20, 2, 3000, 1800, 0, "DESCRIPTIVE_PRIMARY"], ["TARGET_CHROMOSOME", "target_all", 8, 8, 1, 1600, 1600, 1, "EXPLORATORY_SMALL_N"]]))
    domains = ["FOUNDER_IBS", "VARIANT_AGE", "LOCAL_LD", "ROH"]
    comparison_columns = ["SCENARIO_ID", "SCENARIO_SIGNATURE", "ROLE", "DESIGN", "CHANGED_FACTOR", "DOMAIN", "EXPECTED", "SOURCE_RUN_ID", "SOURCE_MANIFEST_SHA256", "SOURCE_CONFIG_SHA256", "SOURCE_STAGE_SIGNATURE", "SOURCE_SUMMARY_SHA256", "EVALUATION_STATUS", "TECHNICAL_STATUS", "PRIMARY_TECHNICAL_STATUS", "CATEGORICAL_COMPARISON", "NUMERIC_METRIC", "PRIMARY_NUMERIC_VALUE", "SCENARIO_NUMERIC_VALUE", "RELATIVE_CHANGE", "QUANTITATIVE_CLASSIFICATION"]
    comparison_rows = [["primary", "aa" * 32, "PRIMARY", "BASELINE", "PRIMARY", domain, "true", "primary_run", "bb" * 32, "cc" * 32, "dd" * 32, "ee" * 32, "EVALUATED", "PRIMARY_STATUS", "PRIMARY_STATUS", "PRIMARY", "metric", 1, 1, 0, "PRIMARY"] for domain in domains]
    comparison_rows += [["window_wide", "ab" * 32, "SENSITIVITY", "SINGLE_FACTOR", "LOCAL_WINDOW", domain, "true", "scenario_run", "bc" * 32, "cd" * 32, "de" * 32, "ef" * 32, "EVALUATED", "PRIMARY_STATUS", "PRIMARY_STATUS", "STABLE", "metric", 1, 1.1, 0.1, "NOT_CLASSIFIED"] for domain in domains]
    definitions.append(("sensitivity_comparisons", "run_sensitivity_analyses", "sensitivity_comparisons.schema.json", comparison_columns, comparison_rows))
    definitions.append(("sensitivity_stability", "run_sensitivity_analyses", "sensitivity_stability.schema.json",
        ["DOMAIN", "PRIMARY_STATUS", "EXPECTED_SCENARIO_COUNT", "EVALUATED_SCENARIO_COUNT", "STABLE_SCENARIO_COUNT", "VARIABLE_SCENARIO_COUNT", "NOT_EVALUATED_SCENARIO_COUNT", "CATEGORICAL_STABILITY", "QUANTITATIVE_STABILITY", "MANUAL_REVIEW_REQUIRED"],
        [[domain, "PRIMARY_STATUS", 1, 1, 1, 0, 0, "STABLE", "NOT_CLASSIFIED", "false"] for domain in domains]))

    artifacts = []
    for artifact_id, producer, schema_name, columns, content in definitions:
        path = run_dir / "stages" / f"source_{producer}" / (f"{artifact_id}.json" if columns is None else f"{artifact_id}.tsv")
        if columns is None:
            path.parent.mkdir(parents=True, exist_ok=True); atomic_write_json(path, content)
        else:
            _tsv(path, columns, content)  # type: ignore[arg-type]
        artifacts.append(build_file_artifact(
            physical_path=path, published_path=path.relative_to(run_dir).as_posix(),
            artifact_id=artifact_id, artifact_type=artifact_id,
            media_type="application/json" if columns is None else "text/tab-separated-values",
            producer_stage=producer, producer_signature=SIGNATURES[producer],
            schema_name=schema_name, schema_version="1.0.0", assembly="GRCh38",
            sample_set_id="synthetic_samples", variant_set_id="target_v1",
            sensitivity="sensitive_genetic",
        ))
    return {"schema_version": "1.0.0", "run_id": "synthetic_visual_run", "stage_id": "18", "stage_name": "build_visualizations", "signature": "18" * 32, "attempt_number": 1, "published_output_dir": "stages/18_build_visualizations", "parameters": {"method": "validated_current_run_figures_v1"}, "artifacts": artifacts}


def _producer_controls(run_dir: Path, stage_inputs: dict[str, Any]) -> None:
    stage_ids = {"analyze_population_structure": "08", "infer_founder_haplotype": "13", "estimate_variant_age": "14", "analyze_local_ld": "15", "analyze_roh": "16", "run_sensitivity_analyses": "17"}
    records = []
    for producer, stage_id in stage_ids.items():
        artifacts = [item for item in stage_inputs["artifacts"] if item["producer_stage"] == producer]
        control_dir = run_dir / "stages" / f"source_{producer}"
        outputs = {"schema_version": "1.0.0", "run_id": stage_inputs["run_id"], "stage_id": stage_id, "stage_name": producer, "signature": SIGNATURES[producer], "artifacts": artifacts}
        outputs_path = control_dir / "stage_outputs.json"; atomic_write_json(outputs_path, outputs)
        audit = {"schema_version": "1.0.0", "run_id": stage_inputs["run_id"], "stage_id": stage_id, "stage_name": producer, "method_id": "synthetic_visual_fixture_v1", "signature": SIGNATURES[producer], "started_at": "2026-08-07T00:00:00Z", "completed_at": "2026-08-07T00:00:00Z", "duration_seconds": 0.0, "inputs": [], "outputs": artifacts, "parameters": {}, "tools": [], "counts": {}, "metrics": {}, "exclusions": [], "warnings": [], "checks": [], "known_limits": [], "expected_visualizations": [], "manual_validation_required": True}
        audit_path = control_dir / "audit.json"; atomic_write_json(audit_path, audit)
        records.append({"stage_id": stage_id, "stage_name": producer, "state": "SUCCEEDED", "signature": SIGNATURES[producer], "audit_path": audit_path.relative_to(run_dir).as_posix(), "audit_sha256": sha256_file(audit_path), "stage_outputs_sha256": sha256_file(outputs_path)})
    atomic_write_json(run_dir / "manifest.json", {"stages": records})


def test_consolidated_figures_are_separate_pseudonymized_and_non_causal(tmp_path: Path) -> None:
    stage_inputs = _fixture(tmp_path)
    results = build_consolidated_figures(run_dir=tmp_path, output_dir=tmp_path / "rendered", stage_inputs=stage_inputs)
    assert {result.domain for result in results} == {"POPULATION_STRUCTURE", "FOUNDER_IBS", "VARIANT_AGE", "LOCAL_LD", "ROH", "SENSITIVITY"}
    assert all(result.status == "RENDERED" for result in results)
    combined = "".join(result.figure_path.read_text(encoding="utf-8") for result in results if result.figure_path)
    assert "PRIVATE_SAMPLE" not in combined
    assert "PRIVATE_PCA" not in combined
    assert "UNIT-001" in combined
    assert "PCA-001" in combined and "Référence indépendante" in combined
    assert "<circle" in combined and "Variance expliquée" in combined
    assert "PRIMARY CORRELATED" in combined and "EXPLORATORY INDEPENDENT" in combined
    assert "composite founder score: NOT CALCULATED" in combined
    rendered = {result.domain: result.figure_path.read_text(encoding="utf-8") for result in results if result.figure_path}
    assert rendered["FOUNDER_IBS"].count("<line") >= 3 and "Longueur partagée" in rendered["FOUNDER_IBS"]
    assert "intervalle de confiance" in rendered["VARIANT_AGE"] and rendered["VARIANT_AGE"].count("<circle") >= 2
    assert "r² génotypique" in rendered["LOCAL_LD"] and rendered["LOCAL_LD"].count("<rect") >= 3
    assert "Charge ROH médiane" in rendered["ROH"] and rendered["ROH"].count("<rect") >= 3
    assert "Variation relative" in rendered["SENSITIVITY"] and rendered["SENSITIVITY"].count("<circle") >= 8


def test_count_incoherence_blocks_only_the_affected_figure(tmp_path: Path) -> None:
    results = build_consolidated_figures(run_dir=tmp_path, output_dir=tmp_path / "rendered", stage_inputs=_fixture(tmp_path, founder_count_mismatch=True))
    statuses = {result.domain: result.status for result in results}
    assert statuses["FOUNDER_IBS"] == "BLOCKED"
    assert all(statuses[domain] == "RENDERED" for domain in statuses if domain != "FOUNDER_IBS")
    assert not (tmp_path / "rendered" / "founder_ibs.svg").exists()


def test_checksum_or_dataset_mismatch_blocks_affected_domain(tmp_path: Path) -> None:
    stage_inputs = _fixture(tmp_path)
    age_scenarios = next(item for item in stage_inputs["artifacts"] if item["artifact_id"] == "variant_age_scenarios")
    age_scenarios["sample_set_id"] = "different_synthetic_samples"
    founder_path = tmp_path / next(item["path"] for item in stage_inputs["artifacts"] if item["artifact_id"] == "founder_segments")
    founder_path.write_text(founder_path.read_text(encoding="utf-8") + "\n", encoding="utf-8")
    results = build_consolidated_figures(run_dir=tmp_path, output_dir=tmp_path / "rendered", stage_inputs=stage_inputs)
    statuses = {result.domain: result.status for result in results}
    assert statuses["FOUNDER_IBS"] == "BLOCKED"
    assert statuses["VARIANT_AGE"] == "BLOCKED"
    assert statuses["LOCAL_LD"] == "RENDERED"


def test_pca_reference_count_mismatch_blocks_only_population_figure(tmp_path: Path) -> None:
    results = build_consolidated_figures(run_dir=tmp_path, output_dir=tmp_path / "rendered", stage_inputs=_fixture(tmp_path, pca_reference_count_mismatch=True))
    statuses = {result.domain: result.status for result in results}
    assert statuses["POPULATION_STRUCTURE"] == "BLOCKED"
    assert all(statuses[domain] == "RENDERED" for domain in statuses if domain != "POPULATION_STRUCTURE")
    assert not (tmp_path / "rendered" / "population_structure.svg").exists()


def test_not_evaluated_is_visible_and_missing_values_are_not_zeroed(tmp_path: Path) -> None:
    results = build_consolidated_figures(run_dir=tmp_path, output_dir=tmp_path / "rendered", stage_inputs=_fixture(tmp_path, ld_not_evaluated=True))
    ld = next(result for result in results if result.domain == "LOCAL_LD")
    assert ld.status == "NOT_EVALUATED"
    svg = ld.figure_path.read_text(encoding="utf-8")
    assert "missing" in svg and "COHORT_TOO_SMALL" in svg
    provenance = json.loads(ld.provenance_path.read_text(encoding="utf-8"))
    assert provenance["missing_value_count"] == 2
    assert provenance["not_evaluated_count"] == 1


def test_stage_18_publishes_versioned_index_completeness_and_audit(tmp_path: Path) -> None:
    run_dir = tmp_path / "run"; output_dir = run_dir / "stages" / ".18_build_visualizations.tmp"; output_dir.mkdir(parents=True)
    stage_inputs = _fixture(run_dir)
    _producer_controls(run_dir, stage_inputs)
    stage_inputs_path = output_dir / "stage_inputs.json"; atomic_write_json(stage_inputs_path, stage_inputs)
    assert execute(stage_inputs_path, output_dir) == 0
    index = json.loads((output_dir / "figure_index.json").read_text(encoding="utf-8"))
    completeness = json.loads((output_dir / "visualization_completeness.json").read_text(encoding="utf-8"))
    audit = json.loads((output_dir / "audit.json").read_text(encoding="utf-8"))
    outputs = json.loads((output_dir / "stage_outputs.json").read_text(encoding="utf-8"))
    render_manifest = json.loads((output_dir / "visualization_render_manifest.json").read_text(encoding="utf-8"))
    validate_json_document(index, "figure_index.schema.json")
    validate_json_document(render_manifest, "visualization_render_manifest.schema.json")
    assert completeness == {"schema_version": "1.0.0", "run_id": "synthetic_visual_run", "expected_domain_count": 6, "rendered_count": 6, "not_evaluated_count": 0, "blocked_count": 0, "complete_for_scientific_report": True}
    assert audit["metrics"]["composite_founder_score_calculated"] is False
    assert audit["metrics"]["sensitivity"] == "sensitive_genetic"
    assert audit["metrics"]["html_rendered"] is True
    assert audit["metrics"]["pdf_rendered"] is True
    artifact_ids = {artifact["artifact_id"] for artifact in outputs["artifacts"]}
    assert {"visualization_gallery_html", "visualization_gallery_pdf", "visualization_render_manifest"} <= artifact_ids
    html_document = (output_dir / "visualization_gallery.html").read_text(encoding="utf-8")
    assert html_document.startswith("<!doctype html>")
    assert html_document.count("<section id=") == 6
    assert "PRIVATE_SAMPLE" not in html_document and "PRIVATE_PCA" not in html_document
    assert html_document.index("population_structure") < html_document.index("founder_ibs")
    assert (output_dir / "visualization_gallery.pdf").read_bytes().startswith(b"%PDF-")
    assert render_manifest["scientific_recalculation_performed"] is False
    assert render_manifest["composite_founder_score_calculated"] is False
