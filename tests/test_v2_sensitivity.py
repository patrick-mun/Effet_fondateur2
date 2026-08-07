from __future__ import annotations

import csv
import json
from pathlib import Path
from typing import Any

import pytest
import yaml

from effet_fondateur.audit import sha256_file
from effet_fondateur.contracts import build_file_artifact, load_pipeline_config, validate_tsv_table
from effet_fondateur.sensitivity import SensitivityAnalysisError, publish_sensitivity_analysis
from effet_fondateur.stages.run_sensitivity_analyses import execute


TIMESTAMP = "2026-08-07T00:00:00Z"
STAGES = {
    "build_sample_registry": ("02", "samples_master", None),
    "prepare_target_variant_dataset": ("04", "target_genotype_audit", None),
    "infer_founder_haplotype": ("13", "founder_analysis_summary", "founder_analysis_summary.schema.json"),
    "estimate_variant_age": ("14", "variant_age_summary", "variant_age_summary.schema.json"),
    "analyze_local_ld": ("15", "local_ld_analysis_summary", "local_ld_analysis_summary.schema.json"),
    "analyze_roh": ("16", "roh_analysis_summary", "roh_analysis_summary.schema.json"),
}


def _summaries(founder_status: str = "SUPPORTED_IBS_CANDIDATE", age: float = 10.0) -> dict[str, dict[str, Any]]:
    return {
        "founder_analysis_summary": {
            "schema_version": "1.0.0", "method_id": "target_centered_exact_ibs_v1",
            "interpretation": "IBS_SHARED_CANDIDATE" if founder_status == "SUPPORTED_IBS_CANDIDATE" else "NO_FOUNDER_CONCLUSION",
            "status": founder_status, "selected_carrier_count": 5, "excluded_carrier_count": 0,
            "background_haplotype_count": 10, "matching_background_haplotype_count": 0,
            "minimum_independent_carriers": 3, "minimum_flank_markers": 2, "ibd_claimed": False,
        },
        "variant_age_summary": {
            "schema_version": "1.0.0", "method_id": "gamma_gandolfo_2014_v1",
            "status": "PRIMARY_ESTIMATE", "primary_model": "CORRELATED", "included_unit_count": 5,
            "minimum_units_for_estimate": 3, "minimum_units_for_primary": 5,
            "estimate_generations": age, "confidence_lower_generations": 5.0,
            "confidence_upper_generations": 20.0, "confidence_level": 0.95,
            "chance_sharing_correction_applied": False, "manual_validation_required": True,
        },
        "local_ld_analysis_summary": {
            "schema_version": "1.0.0", "method_id": "plink19_local_ld_secondary_v1",
            "analysis_role": "SECONDARY_DESCRIPTIVE", "primary_cohort": "controls_unrelated",
            "cohort_statuses": {"controls_unrelated": "DESCRIPTIVE_PRIMARY", "target_carriers_independent": "EXPLORATORY_SMALL_N"},
            "variant_count": 100, "target_variant_id": "target_GRCh38_1_100000_A_G",
            "r2_definition": "GENOTYPE_ALLELE_COUNT_CORRELATION_SQUARED",
            "dprime_definition": "ABSOLUTE_HAPLOTYPE_ML", "feeds_variant_age": False,
        },
        "roh_analysis_summary": {
            "schema_version": "1.0.0", "method_id": "plink19_array_roh_secondary_v1",
            "analysis_role": "SECONDARY_DESCRIPTIVE_AUTOZYGOSITY",
            "scope_statuses": {"genomewide_controls_unrelated": "EVALUATED", "target_chromosome_all_qc": "EVALUATED"},
            "genomewide_common_variant_count": 20000, "target_chromosome_variant_count": 500,
            "target_variant_id": "target_GRCh38_1_100000_A_G", "target_in_roh_count": 1,
            "f_roh_calculated": False, "feeds_founder_haplotype": False,
            "feeds_variant_age": False, "consumes_local_ld": False,
        },
    }


def _audit(run_id: str, stage_id: str, stage_name: str, signature: str, outputs: list[dict[str, Any]]) -> dict[str, Any]:
    return {
        "schema_version": "1.0.0", "run_id": run_id, "stage_id": stage_id,
        "stage_name": stage_name, "method_id": "synthetic_fixture_v1", "signature": signature,
        "started_at": TIMESTAMP, "completed_at": TIMESTAMP, "duration_seconds": 0.0,
        "inputs": [], "outputs": outputs, "parameters": {}, "tools": [], "counts": {},
        "metrics": {}, "exclusions": [], "warnings": [], "checks": [], "known_limits": [],
        "expected_visualizations": [], "manual_validation_required": True,
    }


def _write_source_run(
    root: Path, run_id: str, *, founder_status: str = "SUPPORTED_IBS_CANDIDATE",
    age: float = 10.0, window_delta: int = 0, disallowed_roh_change: bool = False,
    changed_anchor: bool = False,
) -> tuple[Path, str, dict[str, Any]]:
    run_dir = root / run_id
    run_dir.mkdir()
    config = yaml.safe_load(Path("config/pipeline.example.yaml").read_text(encoding="utf-8"))
    config["project"]["project_id"] = "synthetic_sensitivity"
    config["stages"]["prepare_target_region"]["parameters"]["window_left_bp"] += window_delta
    if disallowed_roh_change:
        config["stages"]["analyze_roh"]["parameters"]["minimum_roh_kb"] += 1
    config_path = run_dir / "config.resolved.yaml"
    config_path.write_text(yaml.safe_dump(config, sort_keys=False), encoding="utf-8")
    config_sha = sha256_file(config_path)
    summaries = _summaries(founder_status, age)
    stage_records = []
    for stage_name, (stage_id, artifact_id, schema_name) in STAGES.items():
        stage_dir = run_dir / "stages" / f"{stage_id}_{stage_name}"
        stage_dir.mkdir(parents=True)
        if artifact_id == "samples_master":
            artifact_path = stage_dir / "samples.master.tsv"
            artifact_path.write_text("SYNTHETIC_SAMPLE_ANCHOR\n", encoding="utf-8")
            media_type = "text/tab-separated-values"
        elif artifact_id == "target_genotype_audit":
            artifact_path = stage_dir / "target_genotype_audit.tsv"
            content = "SYNTHETIC_TARGET_ANCHOR_CHANGED\n" if changed_anchor else "SYNTHETIC_TARGET_ANCHOR\n"
            artifact_path.write_text(content, encoding="utf-8")
            media_type = "text/tab-separated-values"
        else:
            artifact_path = stage_dir / f"{artifact_id}.json"
            artifact_path.write_text(json.dumps(summaries[artifact_id], sort_keys=True) + "\n", encoding="utf-8")
            media_type = "application/json"
        signature = (stage_id * 32)[:64]
        artifact = build_file_artifact(
            physical_path=artifact_path,
            published_path=f"stages/{stage_id}_{stage_name}/{artifact_path.name}",
            artifact_id=artifact_id, artifact_type=artifact_id, media_type=media_type,
            producer_stage=stage_name, producer_signature=signature, schema_name=schema_name,
            schema_version="1.0.0" if schema_name else None, assembly="GRCh38",
            sample_set_id="synthetic_samples", variant_set_id="synthetic_variants",
            sensitivity="sensitive_genetic" if artifact_id in {"samples_master", "target_genotype_audit"} else "internal",
        )
        outputs = {
            "schema_version": "1.0.0", "run_id": run_id, "stage_id": stage_id,
            "stage_name": stage_name, "signature": signature, "artifacts": [artifact],
        }
        outputs_path = stage_dir / "stage_outputs.json"
        outputs_path.write_text(json.dumps(outputs, sort_keys=True) + "\n", encoding="utf-8")
        audit_path = stage_dir / "audit.json"
        audit_path.write_text(json.dumps(_audit(run_id, stage_id, stage_name, signature, [artifact]), sort_keys=True) + "\n", encoding="utf-8")
        stage_records.append({
            "stage_id": stage_id, "stage_name": stage_name, "state": "SUCCEEDED",
            "critical": stage_id not in {"15", "16"}, "signature": signature,
            "started_at": TIMESTAMP, "completed_at": TIMESTAMP, "duration_seconds": 0.0,
            "audit_path": f"stages/{stage_id}_{stage_name}/audit.json",
            "audit_sha256": sha256_file(audit_path), "stage_outputs_sha256": sha256_file(outputs_path),
            "attempt_count": 1, "last_error_code": None,
        })
    manifest = {
        "schema_version": "1.0.0", "run_id": run_id, "project_id": "synthetic_sensitivity",
        "created_at": TIMESTAMP, "updated_at": TIMESTAMP, "global_status": "TECHNICALLY_VALID",
        "config_path": "config.resolved.yaml", "config_sha256": config_sha,
        "environment_path": "environment.json", "stages": stage_records,
        "manual_decisions_required": [],
    }
    manifest_path = run_dir / "manifest.json"
    manifest_path.write_text(json.dumps(manifest, sort_keys=True) + "\n", encoding="utf-8")
    return run_dir, sha256_file(manifest_path), config


def _write_registry(path: Path, primary: tuple[Path, str], scenario: tuple[Path, str]) -> None:
    columns = ["SCENARIO_ID", "ROLE", "DESIGN", "CHANGED_FACTOR", "CHANGED_VALUE", "RUN_DIR", "RUN_ID", "MANIFEST_SHA256", "EXPECT_FOUNDER_IBS", "EXPECT_VARIANT_AGE", "EXPECT_LOCAL_LD", "EXPECT_ROH"]
    rows = [
        ["primary", "PRIMARY", "BASELINE", "PRIMARY", "baseline", str(primary[0]), "primary_run", primary[1], "true", "true", "true", "true"],
        ["window_wide", "SENSITIVITY", "SINGLE_FACTOR", "LOCAL_WINDOW", "left_plus_1000bp", str(scenario[0]), "scenario_run", scenario[1], "true", "true", "true", "true"],
    ]
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(columns)
        writer.writerows(rows)


def test_cross_run_sensitivity_publishes_separate_stability_domains(tmp_path: Path) -> None:
    primary_dir, primary_manifest, primary_config = _write_source_run(tmp_path, "primary_run")
    scenario_dir, scenario_manifest, _ = _write_source_run(tmp_path, "scenario_run", age=12.0, window_delta=1000)
    registry = tmp_path / "registry.tsv"
    _write_registry(registry, (primary_dir, primary_manifest), (scenario_dir, scenario_manifest))
    publication = publish_sensitivity_analysis(
        registry_path=registry, output_dir=tmp_path / "out",
        relative_change_tolerances={"FOUNDER_IBS": None, "VARIANT_AGE": 0.25, "LOCAL_LD": None, "ROH": None},
        consolidation_config=primary_config,
    )
    comparisons = validate_tsv_table(publication.comparisons_path, "sensitivity_comparisons.schema.json")
    stability = validate_tsv_table(publication.stability_path, "sensitivity_stability.schema.json")
    assert len(comparisons.rows) == 8
    assert {row["CATEGORICAL_STABILITY"] for row in stability.rows} == {"STABLE"}
    age_row = next(row for row in comparisons.rows if row["SCENARIO_ID"] == "window_wide" and row["DOMAIN"] == "VARIANT_AGE")
    assert age_row["QUANTITATIVE_CLASSIFICATION"] == "WITHIN_TOLERANCE"
    assert publication.domain_stability == {domain: "STABLE" for domain in ("FOUNDER_IBS", "VARIANT_AGE", "LOCAL_LD", "ROH")}


def test_status_change_is_variable_not_a_composite_conclusion(tmp_path: Path) -> None:
    primary_dir, primary_manifest, primary_config = _write_source_run(tmp_path, "primary_run")
    scenario_dir, scenario_manifest, _ = _write_source_run(
        tmp_path, "scenario_run", founder_status="INSUFFICIENT_CARRIERS", window_delta=1000
    )
    registry = tmp_path / "registry.tsv"
    _write_registry(registry, (primary_dir, primary_manifest), (scenario_dir, scenario_manifest))
    publication = publish_sensitivity_analysis(
        registry_path=registry, output_dir=tmp_path / "out",
        relative_change_tolerances={domain: None for domain in ("FOUNDER_IBS", "VARIANT_AGE", "LOCAL_LD", "ROH")},
        consolidation_config=primary_config,
    )
    summary = json.loads(publication.summary_path.read_text(encoding="utf-8"))
    assert summary["domain_stability"]["FOUNDER_IBS"] == "VARIABLE"
    assert summary["composite_founder_score_calculated"] is False


def test_molecular_anchor_change_blocks_sensitivity_comparison(tmp_path: Path) -> None:
    primary_dir, primary_manifest, primary_config = _write_source_run(tmp_path, "primary_run")
    scenario_dir, scenario_manifest, _ = _write_source_run(tmp_path, "scenario_run", window_delta=1000, changed_anchor=True)
    registry = tmp_path / "registry.tsv"
    _write_registry(registry, (primary_dir, primary_manifest), (scenario_dir, scenario_manifest))
    with pytest.raises(SensitivityAnalysisError, match="molecular_or_sample_anchor"):
        publish_sensitivity_analysis(
            registry_path=registry, output_dir=tmp_path / "out",
            relative_change_tolerances={}, consolidation_config=primary_config,
        )


def test_single_factor_scenario_rejects_disallowed_config_change(tmp_path: Path) -> None:
    primary_dir, primary_manifest, primary_config = _write_source_run(tmp_path, "primary_run")
    scenario_dir, scenario_manifest, _ = _write_source_run(
        tmp_path, "scenario_run", window_delta=1000, disallowed_roh_change=True
    )
    registry = tmp_path / "registry.tsv"
    _write_registry(registry, (primary_dir, primary_manifest), (scenario_dir, scenario_manifest))
    with pytest.raises(SensitivityAnalysisError, match="multiple_or_disallowed_axes"):
        publish_sensitivity_analysis(
            registry_path=registry, output_dir=tmp_path / "out",
            relative_change_tolerances={}, consolidation_config=primary_config,
        )


def test_stage_17_publishes_contracts_and_audit_from_synthetic_runs(tmp_path: Path) -> None:
    sources_root = tmp_path / "sources"
    sources_root.mkdir()
    primary_dir, primary_manifest, primary_config = _write_source_run(sources_root, "primary_run")
    scenario_dir, scenario_manifest, _ = _write_source_run(sources_root, "scenario_run", age=11.0, window_delta=1000)
    registry = tmp_path / "registry.tsv"
    _write_registry(registry, (primary_dir, primary_manifest), (scenario_dir, scenario_manifest))
    consolidation_run = tmp_path / "consolidation"
    output_dir = consolidation_run / "stages" / ".17_run_sensitivity_analyses.tmp"
    output_dir.mkdir(parents=True)
    (consolidation_run / "config.resolved.yaml").write_text(
        yaml.safe_dump(primary_config, sort_keys=False), encoding="utf-8"
    )
    signature = "17" * 32
    stage_inputs = {
        "schema_version": "1.0.0", "run_id": "consolidation_run", "stage_id": "17",
        "stage_name": "run_sensitivity_analyses", "signature": signature,
        "attempt_number": 1, "published_output_dir": "stages/17_run_sensitivity_analyses",
        "parameters": {"method": "cross_run_sensitivity_consolidation_v1", "relative_change_tolerances": {"VARIANT_AGE": 0.25}},
        "artifacts": [{
            "artifact_id": "config_input_sensitivity_scenarios", "artifact_type": "sensitivity_scenarios",
            "path": str(registry), "media_type": "text/tab-separated-values", "schema_name": None,
            "schema_version": None, "sha256": sha256_file(registry), "producer_stage": "external_source",
            "producer_signature": "ab" * 32, "assembly": "GRCh38", "sample_set_id": None,
            "variant_set_id": None, "sensitivity": "sensitive_genetic",
        }],
    }
    stage_inputs_path = output_dir / "stage_inputs.json"
    stage_inputs_path.write_text(json.dumps(stage_inputs, sort_keys=True) + "\n", encoding="utf-8")
    assert execute(stage_inputs_path, output_dir) == 0
    outputs = json.loads((output_dir / "stage_outputs.json").read_text(encoding="utf-8"))
    audit = json.loads((output_dir / "audit.json").read_text(encoding="utf-8"))
    assert {artifact["artifact_id"] for artifact in outputs["artifacts"]} == {
        "sensitivity_comparisons", "sensitivity_stability", "sensitivity_analysis_summary"
    }
    assert audit["metrics"]["composite_founder_score_calculated"] is False
    assert audit["counts"]["scenarios"] == 2
