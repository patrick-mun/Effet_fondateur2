"""Validation et comparaison de runs V2 préspécifiés pour l'étape 17."""

from __future__ import annotations

import csv
import hashlib
import json
import math
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Iterable

from effet_fondateur.audit import atomic_write_json, read_json, sha256_file
from effet_fondateur.contracts import (
    load_pipeline_config,
    validate_json_document,
    validate_tsv_table,
)


DOMAINS = ("FOUNDER_IBS", "VARIANT_AGE", "LOCAL_LD", "ROH")
DOMAIN_STAGE = {
    "FOUNDER_IBS": ("infer_founder_haplotype", "founder_analysis_summary", "founder_analysis_summary.schema.json"),
    "VARIANT_AGE": ("estimate_variant_age", "variant_age_summary", "variant_age_summary.schema.json"),
    "LOCAL_LD": ("analyze_local_ld", "local_ld_analysis_summary", "local_ld_analysis_summary.schema.json"),
    "ROH": ("analyze_roh", "roh_analysis_summary", "roh_analysis_summary.schema.json"),
}
EXPECT_COLUMNS = {
    "FOUNDER_IBS": "EXPECT_FOUNDER_IBS",
    "VARIANT_AGE": "EXPECT_VARIANT_AGE",
    "LOCAL_LD": "EXPECT_LOCAL_LD",
    "ROH": "EXPECT_ROH",
}
COMPARISON_COLUMNS = (
    "SCENARIO_ID", "SCENARIO_SIGNATURE", "ROLE", "DESIGN", "CHANGED_FACTOR",
    "DOMAIN", "EXPECTED", "SOURCE_RUN_ID", "SOURCE_MANIFEST_SHA256",
    "SOURCE_CONFIG_SHA256", "SOURCE_STAGE_SIGNATURE", "SOURCE_SUMMARY_SHA256",
    "EVALUATION_STATUS", "TECHNICAL_STATUS",
    "PRIMARY_TECHNICAL_STATUS", "CATEGORICAL_COMPARISON", "NUMERIC_METRIC",
    "PRIMARY_NUMERIC_VALUE", "SCENARIO_NUMERIC_VALUE", "RELATIVE_CHANGE",
    "QUANTITATIVE_CLASSIFICATION",
)
STABILITY_COLUMNS = (
    "DOMAIN", "PRIMARY_STATUS", "EXPECTED_SCENARIO_COUNT",
    "EVALUATED_SCENARIO_COUNT", "STABLE_SCENARIO_COUNT",
    "VARIABLE_SCENARIO_COUNT", "NOT_EVALUATED_SCENARIO_COUNT",
    "CATEGORICAL_STABILITY", "QUANTITATIVE_STABILITY",
    "MANUAL_REVIEW_REQUIRED",
)


class SensitivityAnalysisError(ValueError):
    """Signale un registre ou un run source non comparable."""


@dataclass(frozen=True)
class SensitivityPublication:
    """Chemins et compteurs publiés par le consolidateur."""

    comparisons_path: Path
    stability_path: Path
    summary_path: Path
    scenario_count: int
    evaluated_comparison_count: int
    domain_stability: dict[str, str]


@dataclass(frozen=True)
class _SourceRun:
    row: dict[str, Any]
    run_dir: Path
    manifest: dict[str, Any]
    config: dict[str, Any]
    config_sha256: str
    anchor_hashes: dict[str, str]


def _write_tsv(path: Path, columns: Iterable[str], rows: Iterable[dict[str, Any]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=tuple(columns), delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for row in rows:
            writer.writerow({key: "" if value is None else str(value).lower() if isinstance(value, bool) else value for key, value in row.items()})


def _contained_path(run_dir: Path, relative_path: str) -> Path:
    candidate = (run_dir / relative_path).resolve()
    resolved_root = run_dir.resolve()
    if candidate != resolved_root and resolved_root not in candidate.parents:
        raise SensitivityAnalysisError("source_artifact_outside_run")
    return candidate


def _stage_record(manifest: dict[str, Any], stage_name: str) -> dict[str, Any] | None:
    matches = [record for record in manifest["stages"] if record["stage_name"] == stage_name]
    if len(matches) > 1:
        raise SensitivityAnalysisError(f"ambiguous_source_stage:{stage_name}")
    return matches[0] if matches else None


def _validated_stage_outputs(source: _SourceRun, stage_name: str) -> tuple[dict[str, Any], str] | None:
    record = _stage_record(source.manifest, stage_name)
    if record is None or record["state"] not in {"SUCCEEDED", "CACHED"}:
        return None
    if not record["audit_path"] or not record["audit_sha256"] or not record["stage_outputs_sha256"] or not record["signature"]:
        raise SensitivityAnalysisError(f"incomplete_source_stage_record:{stage_name}")
    audit_path = _contained_path(source.run_dir, record["audit_path"])
    outputs_path = audit_path.parent / "stage_outputs.json"
    if not audit_path.is_file() or sha256_file(audit_path) != record["audit_sha256"]:
        raise SensitivityAnalysisError(f"source_audit_checksum_mismatch:{stage_name}")
    if not outputs_path.is_file() or sha256_file(outputs_path) != record["stage_outputs_sha256"]:
        raise SensitivityAnalysisError(f"source_outputs_checksum_mismatch:{stage_name}")
    audit = read_json(audit_path)
    outputs = read_json(outputs_path)
    validate_json_document(audit, "stage_audit.schema.json")
    validate_json_document(outputs, "stage_outputs.schema.json")
    if outputs["signature"] != record["signature"] or audit["signature"] != record["signature"]:
        raise SensitivityAnalysisError(f"source_stage_signature_mismatch:{stage_name}")
    if outputs["run_id"] != source.manifest["run_id"] or audit["run_id"] != source.manifest["run_id"]:
        raise SensitivityAnalysisError(f"source_stage_run_id_mismatch:{stage_name}")
    if outputs["stage_name"] != stage_name or audit["stage_name"] != stage_name:
        raise SensitivityAnalysisError(f"source_stage_name_mismatch:{stage_name}")
    return outputs, record["signature"]


def _artifact_from_outputs(source: _SourceRun, outputs: dict[str, Any], artifact_id: str) -> tuple[dict[str, Any], Path]:
    matches = [artifact for artifact in outputs["artifacts"] if artifact["artifact_id"] == artifact_id]
    if len(matches) != 1:
        raise SensitivityAnalysisError(f"source_artifact_missing_or_ambiguous:{artifact_id}")
    artifact = matches[0]
    path = _contained_path(source.run_dir, artifact["path"])
    if not path.is_file() or sha256_file(path) != artifact["sha256"]:
        raise SensitivityAnalysisError(f"source_artifact_checksum_mismatch:{artifact_id}")
    return artifact, path


def _anchor_hash(run_dir: Path, manifest: dict[str, Any], stage_name: str, artifact_id: str) -> str:
    placeholder = _SourceRun({}, run_dir, manifest, {}, manifest["config_sha256"], {})
    validated = _validated_stage_outputs(placeholder, stage_name)
    if validated is None:
        raise SensitivityAnalysisError(f"identity_anchor_stage_not_available:{stage_name}")
    artifact, _ = _artifact_from_outputs(placeholder, validated[0], artifact_id)
    return artifact["sha256"]


def _load_source_run(row: dict[str, Any], registry_dir: Path) -> _SourceRun:
    configured = Path(row["RUN_DIR"])
    run_dir = configured if configured.is_absolute() else registry_dir / configured
    if run_dir.is_symlink() or not run_dir.is_dir():
        raise SensitivityAnalysisError("source_run_directory_invalid")
    manifest_path = run_dir / "manifest.json"
    if not manifest_path.is_file() or sha256_file(manifest_path) != row["MANIFEST_SHA256"]:
        raise SensitivityAnalysisError("source_manifest_checksum_mismatch")
    manifest = read_json(manifest_path)
    validate_json_document(manifest, "run_manifest.schema.json")
    if manifest["run_id"] != row["RUN_ID"]:
        raise SensitivityAnalysisError("source_run_id_mismatch")
    config_path = run_dir / "config.resolved.yaml"
    if not config_path.is_file() or sha256_file(config_path) != manifest["config_sha256"]:
        raise SensitivityAnalysisError("source_config_checksum_mismatch")
    config = load_pipeline_config(config_path)
    anchors = {
        "samples_master": _anchor_hash(run_dir, manifest, "build_sample_registry", "samples_master"),
        "target_genotype_audit": _anchor_hash(run_dir, manifest, "prepare_target_variant_dataset", "target_genotype_audit"),
    }
    return _SourceRun(row, run_dir, manifest, config, manifest["config_sha256"], anchors)


def _normalized_for_diff(config: dict[str, Any]) -> dict[str, Any]:
    normalized = json.loads(json.dumps(config))
    normalized.get("inputs", {}).pop("sensitivity_scenarios", None)
    normalized.get("stages", {}).pop("run_sensitivity_analyses", None)
    return normalized


def _leaf_differences(left: Any, right: Any, prefix: str = "") -> set[str]:
    if isinstance(left, dict) and isinstance(right, dict):
        differences: set[str] = set()
        for key in set(left) | set(right):
            path = f"{prefix}.{key}" if prefix else key
            if key not in left or key not in right:
                differences.add(path)
            else:
                differences.update(_leaf_differences(left[key], right[key], path))
        return differences
    return set() if left == right else {prefix}


ALLOWED_FACTOR_PREFIXES = {
    "LEAVE_ONE_FAMILY_OUT": ("stages.freeze_cohorts.parameters", "stages.infer_founder_haplotype.parameters"),
    "LOW_PHASE_CONFIDENCE_CARRIER_EXCLUSION": ("stages.phase_target_region.parameters", "stages.freeze_cohorts.parameters", "stages.infer_founder_haplotype.parameters"),
    "HAPLOTYPE_DEFINITION": ("stages.infer_founder_haplotype.parameters",),
    "LOCAL_WINDOW": ("stages.prepare_target_region.parameters",),
    "GENETIC_MAP_INTERPOLATION": ("inputs.genetic_map", "stages.prepare_target_region.parameters"),
    "GENOMEWIDE_LD_PRUNING": ("stages.build_kinship_panel.parameters",),
    "DISTANT_RELATIVE_POLICY": ("stages.infer_kinship.parameters", "stages.freeze_cohorts.parameters"),
    "ROH_PROFILE": ("stages.analyze_roh.parameters",),
    "REFERENCE_POPULATION": ("inputs.reference_panel_catalog", "stages.phase_target_region.parameters"),
}


def _validate_comparability(primary: _SourceRun, scenario: _SourceRun) -> None:
    if scenario.config["project"]["assembly"] != primary.config["project"]["assembly"]:
        raise SensitivityAnalysisError("scenario_assembly_mismatch")
    if scenario.config["target"] != primary.config["target"]:
        raise SensitivityAnalysisError("scenario_target_definition_mismatch")
    if scenario.anchor_hashes != primary.anchor_hashes:
        raise SensitivityAnalysisError("scenario_molecular_or_sample_anchor_mismatch")
    differences = _leaf_differences(
        _normalized_for_diff(primary.config), _normalized_for_diff(scenario.config)
    )
    factor = scenario.row["CHANGED_FACTOR"]
    if scenario.row["DESIGN"] == "MULTIFACTOR":
        return
    allowed_prefixes = ALLOWED_FACTOR_PREFIXES.get(factor, ())
    disallowed = sorted(
        path for path in differences
        if not any(path == prefix or path.startswith(prefix + ".") for prefix in allowed_prefixes)
    )
    if disallowed:
        raise SensitivityAnalysisError(f"scenario_changes_multiple_or_disallowed_axes:{disallowed[0]}")


def _technical_endpoint(domain: str, summary: dict[str, Any]) -> tuple[str, str | None, float | None]:
    if domain == "FOUNDER_IBS":
        return summary["status"], None, None
    if domain == "VARIANT_AGE":
        value = summary["estimate_generations"]
        return summary["status"], "estimate_generations", None if value is None else float(value)
    if domain == "LOCAL_LD":
        statuses = summary["cohort_statuses"]
        status = ";".join(f"{key}={statuses[key]}" for key in sorted(statuses))
        return status, None, None
    statuses = summary["scope_statuses"]
    status = ";".join(f"{key}={statuses[key]}" for key in sorted(statuses))
    return status, "target_in_roh_count", float(summary["target_in_roh_count"])


def _domain_result(source: _SourceRun, domain: str, expected: bool) -> dict[str, Any]:
    if not expected:
        return {"stage_signature": None, "summary_sha256": None, "evaluation": "NOT_EVALUATED", "status": None, "metric": None, "value": None}
    stage_name, artifact_id, schema_name = DOMAIN_STAGE[domain]
    validated = _validated_stage_outputs(source, stage_name)
    if validated is None:
        return {"stage_signature": None, "summary_sha256": None, "evaluation": "NOT_EVALUATED", "status": None, "metric": None, "value": None}
    outputs, signature = validated
    artifact, summary_path = _artifact_from_outputs(source, outputs, artifact_id)
    summary = read_json(summary_path)
    validate_json_document(summary, schema_name)
    status, metric, value = _technical_endpoint(domain, summary)
    return {"stage_signature": signature, "summary_sha256": artifact["sha256"], "evaluation": "EVALUATED", "status": status, "metric": metric, "value": value}


def _scenario_signature(source: _SourceRun, domain_results: dict[str, dict[str, Any]]) -> str:
    payload = {
        "scenario_id": source.row["SCENARIO_ID"], "role": source.row["ROLE"],
        "design": source.row["DESIGN"], "changed_factor": source.row["CHANGED_FACTOR"],
        "changed_value": source.row["CHANGED_VALUE"], "run_id": source.row["RUN_ID"],
        "config_sha256": source.config_sha256,
        "stage_signatures": {domain: domain_results[domain]["stage_signature"] for domain in DOMAINS},
        "summary_sha256": {domain: domain_results[domain]["summary_sha256"] for domain in DOMAINS},
    }
    encoded = json.dumps(payload, sort_keys=True, separators=(",", ":")).encode("utf-8")
    return hashlib.sha256(encoded).hexdigest()


def _is_conclusive(domain: str, status: str | None) -> bool:
    if status is None:
        return False
    if domain == "FOUNDER_IBS":
        return status == "SUPPORTED_IBS_CANDIDATE"
    if domain == "VARIANT_AGE":
        return status == "PRIMARY_ESTIMATE"
    if domain == "LOCAL_LD":
        return "DESCRIPTIVE_PRIMARY" in status
    return "EVALUATED" in status


def _format_number(value: float | None) -> str | None:
    if value is None:
        return None
    if not math.isfinite(value):
        raise SensitivityAnalysisError("non_finite_sensitivity_metric")
    return format(value, ".12g")


def publish_sensitivity_analysis(
    *, registry_path: Path, output_dir: Path,
    relative_change_tolerances: dict[str, float | None],
    consolidation_config: dict[str, Any],
) -> SensitivityPublication:
    """Valide les runs déclarés et publie leurs comparaisons sans les recalculer."""
    table = validate_tsv_table(registry_path, "sensitivity_scenario_registry.schema.json")
    rows = list(table.rows)
    primary_rows = [row for row in rows if row["ROLE"] == "PRIMARY"]
    if len(primary_rows) != 1:
        raise SensitivityAnalysisError("exactly_one_primary_scenario_required")
    primary_row = primary_rows[0]
    if primary_row["DESIGN"] != "BASELINE" or primary_row["CHANGED_FACTOR"] != "PRIMARY":
        raise SensitivityAnalysisError("invalid_primary_scenario_definition")
    for row in rows:
        if not any(row[column] for column in EXPECT_COLUMNS.values()):
            raise SensitivityAnalysisError("scenario_without_expected_domain")
        if row["ROLE"] == "SENSITIVITY" and row["DESIGN"] == "BASELINE":
            raise SensitivityAnalysisError("sensitivity_scenario_cannot_be_baseline")
        if row["DESIGN"] == "MULTIFACTOR" and row["CHANGED_FACTOR"] != "MULTIFACTOR":
            raise SensitivityAnalysisError("multifactor_label_mismatch")
    sources = [_load_source_run(row, registry_path.parent) for row in rows]
    primary = next(source for source in sources if source.row["ROLE"] == "PRIMARY")
    if consolidation_config["project"]["assembly"] != primary.config["project"]["assembly"] or consolidation_config["target"] != primary.config["target"]:
        raise SensitivityAnalysisError("consolidation_run_not_comparable_to_primary")
    for source in sources:
        if source is not primary:
            _validate_comparability(primary, source)

    raw_results: dict[str, dict[str, dict[str, Any]]] = {}
    scenario_signatures: dict[str, str] = {}
    for source in sources:
        results = {
            domain: _domain_result(source, domain, bool(source.row[EXPECT_COLUMNS[domain]]))
            for domain in DOMAINS
        }
        raw_results[source.row["SCENARIO_ID"]] = results
        scenario_signatures[source.row["SCENARIO_ID"]] = _scenario_signature(source, results)

    primary_results = raw_results[primary.row["SCENARIO_ID"]]
    comparison_rows: list[dict[str, Any]] = []
    for source in sources:
        scenario_id = source.row["SCENARIO_ID"]
        for domain in DOMAINS:
            result = raw_results[scenario_id][domain]
            primary_result = primary_results[domain]
            expected = bool(source.row[EXPECT_COLUMNS[domain]])
            if source.row["ROLE"] == "PRIMARY":
                categorical = "PRIMARY"
            elif result["evaluation"] != "EVALUATED" or primary_result["evaluation"] != "EVALUATED":
                categorical = "NOT_COMPARABLE"
            else:
                categorical = "STABLE" if result["status"] == primary_result["status"] else "VARIABLE"
            relative_change = None
            primary_value = primary_result["value"]
            scenario_value = result["value"]
            if source.row["ROLE"] == "PRIMARY":
                quantitative = "PRIMARY"
            elif primary_value is None or scenario_value is None:
                quantitative = "NOT_AVAILABLE"
            else:
                if primary_value != 0:
                    relative_change = (scenario_value - primary_value) / abs(primary_value)
                tolerance = relative_change_tolerances.get(domain)
                if tolerance is None or relative_change is None:
                    quantitative = "NOT_CLASSIFIED"
                else:
                    quantitative = "WITHIN_TOLERANCE" if abs(relative_change) <= tolerance else "OUTSIDE_TOLERANCE"
            comparison_rows.append({
                "SCENARIO_ID": scenario_id,
                "SCENARIO_SIGNATURE": scenario_signatures[scenario_id],
                "ROLE": source.row["ROLE"], "DESIGN": source.row["DESIGN"],
                "CHANGED_FACTOR": source.row["CHANGED_FACTOR"], "DOMAIN": domain,
                "EXPECTED": expected, "SOURCE_RUN_ID": source.row["RUN_ID"],
                "SOURCE_MANIFEST_SHA256": source.row["MANIFEST_SHA256"],
                "SOURCE_CONFIG_SHA256": source.config_sha256,
                "SOURCE_STAGE_SIGNATURE": result["stage_signature"],
                "SOURCE_SUMMARY_SHA256": result["summary_sha256"],
                "EVALUATION_STATUS": result["evaluation"], "TECHNICAL_STATUS": result["status"],
                "PRIMARY_TECHNICAL_STATUS": primary_result["status"],
                "CATEGORICAL_COMPARISON": categorical, "NUMERIC_METRIC": result["metric"],
                "PRIMARY_NUMERIC_VALUE": _format_number(primary_value),
                "SCENARIO_NUMERIC_VALUE": _format_number(scenario_value),
                "RELATIVE_CHANGE": _format_number(relative_change),
                "QUANTITATIVE_CLASSIFICATION": quantitative,
            })

    stability_rows: list[dict[str, Any]] = []
    domain_stability: dict[str, str] = {}
    for domain in DOMAINS:
        primary_result = primary_results[domain]
        candidates = [
            row for row in comparison_rows
            if row["DOMAIN"] == domain and row["ROLE"] == "SENSITIVITY" and row["EXPECTED"]
        ]
        evaluated = [row for row in candidates if row["EVALUATION_STATUS"] == "EVALUATED"]
        stable_count = sum(row["CATEGORICAL_COMPARISON"] == "STABLE" for row in candidates)
        variable_count = sum(row["CATEGORICAL_COMPARISON"] == "VARIABLE" for row in candidates)
        missing_count = len(candidates) - len(evaluated)
        if not candidates or primary_result["evaluation"] != "EVALUATED":
            categorical_stability = "NOT_EVALUATED"
        elif not _is_conclusive(domain, primary_result["status"]):
            categorical_stability = "INCONCLUSIVE"
        elif variable_count:
            categorical_stability = "VARIABLE"
        elif missing_count:
            categorical_stability = "INCONCLUSIVE"
        else:
            categorical_stability = "STABLE"
        numeric_candidates = [row for row in candidates if row["PRIMARY_NUMERIC_VALUE"] is not None]
        if not numeric_candidates:
            quantitative_stability = "NOT_AVAILABLE"
        elif relative_change_tolerances.get(domain) is None:
            quantitative_stability = "NOT_CLASSIFIED"
        elif any(row["QUANTITATIVE_CLASSIFICATION"] == "OUTSIDE_TOLERANCE" for row in numeric_candidates):
            quantitative_stability = "VARIABLE"
        elif any(row["QUANTITATIVE_CLASSIFICATION"] != "WITHIN_TOLERANCE" for row in numeric_candidates):
            quantitative_stability = "INCONCLUSIVE"
        else:
            quantitative_stability = "STABLE"
        manual_review = categorical_stability != "STABLE" or quantitative_stability in {"VARIABLE", "INCONCLUSIVE"}
        domain_stability[domain] = categorical_stability
        stability_rows.append({
            "DOMAIN": domain, "PRIMARY_STATUS": primary_result["status"],
            "EXPECTED_SCENARIO_COUNT": len(candidates), "EVALUATED_SCENARIO_COUNT": len(evaluated),
            "STABLE_SCENARIO_COUNT": stable_count, "VARIABLE_SCENARIO_COUNT": variable_count,
            "NOT_EVALUATED_SCENARIO_COUNT": missing_count,
            "CATEGORICAL_STABILITY": categorical_stability,
            "QUANTITATIVE_STABILITY": quantitative_stability,
            "MANUAL_REVIEW_REQUIRED": manual_review,
        })

    comparisons_path = output_dir / "sensitivity_comparisons.tsv"
    stability_path = output_dir / "sensitivity_stability.tsv"
    summary_path = output_dir / "sensitivity_analysis_summary.json"
    _write_tsv(comparisons_path, COMPARISON_COLUMNS, comparison_rows)
    _write_tsv(stability_path, STABILITY_COLUMNS, stability_rows)
    summary = {
        "schema_version": "1.0.0", "method_id": "cross_run_sensitivity_consolidation_v1",
        "analysis_role": "ROBUSTNESS_NOT_EXTERNAL_VALIDATION",
        "primary_scenario_id": primary.row["SCENARIO_ID"], "scenario_count": len(sources),
        "domain_stability": domain_stability, "composite_founder_score_calculated": False,
        "manual_validation_required": True,
    }
    atomic_write_json(summary_path, summary)
    validate_tsv_table(comparisons_path, "sensitivity_comparisons.schema.json")
    validate_tsv_table(stability_path, "sensitivity_stability.schema.json")
    validate_json_document(summary, "sensitivity_analysis_summary.schema.json")
    return SensitivityPublication(
        comparisons_path=comparisons_path, stability_path=stability_path,
        summary_path=summary_path, scenario_count=len(sources),
        evaluated_comparison_count=sum(row["EVALUATION_STATUS"] == "EVALUATED" for row in comparison_rows),
        domain_stability=domain_stability,
    )
