"""Étape 17 : consolidation inter-runs des analyses de sensibilité."""

from __future__ import annotations

import argparse
import json
import math
import sys
from pathlib import Path
from time import monotonic
from typing import Any, Sequence

from effet_fondateur.audit import atomic_write_json, read_json, sha256_file
from effet_fondateur.contracts import (
    DocumentValidationError,
    TableValidationError,
    build_file_artifact,
    load_pipeline_config,
    validate_json_document,
)
from effet_fondateur.orchestrator.state import utc_now
from effet_fondateur.sensitivity import SensitivityAnalysisError, publish_sensitivity_analysis


class SensitivityStageInputError(ValueError):
    """Signale une entrée orchestrée invalide de l'étape 17."""


def _build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(prog="run-sensitivity-analyses")
    parser.add_argument("--stage-inputs", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    return parser


def _parameters(parameters: dict[str, Any]) -> dict[str, Any]:
    if parameters.get("method", "cross_run_sensitivity_consolidation_v1") != "cross_run_sensitivity_consolidation_v1":
        raise SensitivityStageInputError("invalid_parameter:method")
    raw_tolerances = parameters.get("relative_change_tolerances", {})
    domains = ("FOUNDER_IBS", "VARIANT_AGE", "LOCAL_LD", "ROH")
    if not isinstance(raw_tolerances, dict) or set(raw_tolerances) - set(domains):
        raise SensitivityStageInputError("invalid_parameter:relative_change_tolerances")
    tolerances: dict[str, float | None] = {}
    for domain in domains:
        value = raw_tolerances.get(domain)
        if value is None:
            tolerances[domain] = None
        elif isinstance(value, bool) or not isinstance(value, (int, float)) or not math.isfinite(value) or value < 0:
            raise SensitivityStageInputError(f"invalid_tolerance:{domain}")
        else:
            tolerances[domain] = float(value)
    return {"method": "cross_run_sensitivity_consolidation_v1", "relative_change_tolerances": tolerances}


def _registry_artifact(stage_inputs: dict[str, Any]) -> dict[str, Any]:
    matches = [artifact for artifact in stage_inputs["artifacts"] if artifact["artifact_id"] == "config_input_sensitivity_scenarios"]
    if len(matches) != 1:
        raise SensitivityStageInputError("sensitivity_registry_missing_or_ambiguous")
    return matches[0]


def execute(stage_inputs_path: Path, output_dir: Path) -> int:
    """Valide et consolide des runs signés sans relancer leurs analyses."""
    started_at = utc_now()
    started_clock = monotonic()
    stage_inputs = read_json(stage_inputs_path)
    validate_json_document(stage_inputs, "stage_inputs.schema.json")
    run_dir = output_dir.parent.parent
    config = load_pipeline_config(run_dir / "config.resolved.yaml")
    parameters = _parameters(stage_inputs["parameters"])
    registry_artifact = _registry_artifact(stage_inputs)
    registry_path = Path(registry_artifact["path"])
    if not registry_path.is_absolute():
        registry_path = Path.cwd() / registry_path
    if not registry_path.is_file() or sha256_file(registry_path) != registry_artifact["sha256"]:
        raise SensitivityStageInputError("sensitivity_registry_missing_or_modified")
    publication = publish_sensitivity_analysis(
        registry_path=registry_path,
        output_dir=output_dir / "sensitivity",
        relative_change_tolerances=parameters["relative_change_tolerances"],
        consolidation_config=config,
    )
    specifications = (
        ("sensitivity_comparisons", publication.comparisons_path, "sensitivity_comparisons.schema.json", "sensitive_genetic"),
        ("sensitivity_stability", publication.stability_path, "sensitivity_stability.schema.json", "internal"),
        ("sensitivity_analysis_summary", publication.summary_path, "sensitivity_analysis_summary.schema.json", "internal"),
    )
    output_artifacts = [
        build_file_artifact(
            physical_path=path,
            published_path=f"{stage_inputs['published_output_dir']}/{path.relative_to(output_dir).as_posix()}",
            artifact_id=artifact_id, artifact_type=artifact_id,
            media_type="application/json" if path.suffix == ".json" else "text/tab-separated-values",
            producer_stage=stage_inputs["stage_name"], producer_signature=stage_inputs["signature"],
            schema_name=schema_name, schema_version="1.0.0",
            assembly=config["project"]["assembly"], sample_set_id=None,
            variant_set_id=config["target"]["project_variant_id"], sensitivity=sensitivity,
        )
        for artifact_id, path, schema_name, sensitivity in specifications
    ]
    stage_outputs = {
        "schema_version": "1.0.0", "run_id": stage_inputs["run_id"],
        "stage_id": stage_inputs["stage_id"], "stage_name": stage_inputs["stage_name"],
        "signature": stage_inputs["signature"], "artifacts": output_artifacts,
    }
    validate_json_document(stage_outputs, "stage_outputs.schema.json")
    atomic_write_json(output_dir / "stage_outputs.json", stage_outputs)
    warnings = [
        {"code": f"{domain.lower()}_{status.lower()}", "count": 1}
        for domain, status in publication.domain_stability.items() if status != "STABLE"
    ]
    audit = {
        "schema_version": "1.0.0", "run_id": stage_inputs["run_id"],
        "stage_id": stage_inputs["stage_id"], "stage_name": stage_inputs["stage_name"],
        "method_id": "cross_run_sensitivity_consolidation_v1", "signature": stage_inputs["signature"],
        "started_at": started_at, "completed_at": utc_now(), "duration_seconds": monotonic() - started_clock,
        "inputs": [registry_artifact], "outputs": output_artifacts, "parameters": parameters,
        "tools": [{"tool": "python", "configured": Path(sys.executable).name, "version": sys.version.split()[0]}],
        "counts": {"scenarios": publication.scenario_count, "evaluated_domain_comparisons": publication.evaluated_comparison_count},
        "metrics": {"analysis_role": "ROBUSTNESS_NOT_EXTERNAL_VALIDATION", "domain_stability": publication.domain_stability, "composite_founder_score_calculated": False},
        "exclusions": [], "warnings": warnings,
        "checks": [
            {"check": "registry_contract", "status": "PASS"},
            {"check": "source_manifests_and_artifacts_integrity", "status": "PASS"},
            {"check": "molecular_and_sample_identity_anchors", "status": "PASS"},
            {"check": "single_factor_configuration_isolation", "status": "PASS"},
            {"check": "domain_separation_without_composite_score", "status": "PASS"},
        ],
        "known_limits": [
            "La robustesse ne porte que sur les scénarios préspécifiés et disponibles.",
            "Un biais commun à tous les runs n'est pas détecté par cette consolidation.",
            "La stabilité catégorielle ne prouve ni IBD, ni origine unique, ni causalité.",
            "Cette étape ne constitue pas une validation externe dans une autre population ou technologie.",
        ],
        "expected_visualizations": ["plot_sensitivity_stability", "plot_sensitivity_forest", "plot_sensitivity_status_changes"],
        "manual_validation_required": True,
    }
    validate_json_document(audit, "stage_audit.schema.json")
    atomic_write_json(output_dir / "audit.json", audit)
    (output_dir / "checksums.sha256").write_text(
        "".join(f"{artifact['sha256']}  {artifact['path'].removeprefix(stage_inputs['published_output_dir'] + '/')}\n" for artifact in output_artifacts),
        encoding="utf-8",
    )
    return 0


def _classify_error(error: Exception) -> int:
    if isinstance(error, SensitivityAnalysisError):
        return 4
    if isinstance(error, (SensitivityStageInputError, DocumentValidationError, TableValidationError, OSError, ValueError, json.JSONDecodeError)):
        return 2
    return 5


def main(arguments: Sequence[str] | None = None) -> int:
    """Point d'entrée autonome respectant les codes de retour V2."""
    parsed = _build_parser().parse_args(arguments)
    try:
        return execute(parsed.stage_inputs, parsed.output_dir)
    except Exception as error:
        sys.stderr.write(f"{error}\n")
        return _classify_error(error)


if __name__ == "__main__":
    raise SystemExit(main())
