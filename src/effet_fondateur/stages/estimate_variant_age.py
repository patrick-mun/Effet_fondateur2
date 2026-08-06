"""Étape 14 : datation Gamma du variant fondateur candidat."""

from __future__ import annotations

import argparse
import importlib.metadata
import json
import math
import sys
from pathlib import Path, PurePosixPath
from time import monotonic
from typing import Any, Sequence

from effet_fondateur.audit import atomic_write_json, read_json, sha256_file
from effet_fondateur.contracts import (
    DocumentValidationError,
    build_file_artifact,
    load_pipeline_config,
    validate_json_document,
)
from effet_fondateur.dating import GammaAgeError, publish_variant_age
from effet_fondateur.orchestrator.state import utc_now


class EstimateVariantAgeInputError(ValueError):
    """Signale une configuration ou une entrée invalide de l'étape 14."""


def _build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(prog="estimate-variant-age")
    parser.add_argument("--stage-inputs", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    return parser


def _integer(parameters: dict[str, Any], name: str, default: int) -> int:
    value = parameters.get(name, default)
    if isinstance(value, bool) or not isinstance(value, int) or value <= 0:
        raise EstimateVariantAgeInputError(f"invalid_parameter:{name}")
    return value


def _number(parameters: dict[str, Any], name: str, default: float) -> float:
    value = parameters.get(name, default)
    if isinstance(value, bool) or not isinstance(value, (int, float)) or not math.isfinite(value):
        raise EstimateVariantAgeInputError(f"invalid_parameter:{name}")
    return float(value)


def _parameters(parameters: dict[str, Any]) -> dict[str, Any]:
    method = parameters.get("method", "gamma_gandolfo_2014_v1")
    if method != "gamma_gandolfo_2014_v1":
        raise EstimateVariantAgeInputError("invalid_parameter:method")
    confidence_level = _number(parameters, "confidence_level", 0.95)
    if not 0.5 < confidence_level < 1:
        raise EstimateVariantAgeInputError("invalid_parameter:confidence_level")
    generation_years = parameters.get("generation_years", [25, 28, 30])
    if (
        not isinstance(generation_years, list) or not generation_years
        or any(isinstance(value, bool) or not isinstance(value, (int, float)) or value <= 0 for value in generation_years)
        or len(set(generation_years)) != len(generation_years)
    ):
        raise EstimateVariantAgeInputError("invalid_parameter:generation_years")
    leave_one_out = parameters.get("leave_one_family_out", True)
    chance_sharing = parameters.get("chance_sharing_correction", False)
    if not isinstance(leave_one_out, bool) or not isinstance(chance_sharing, bool):
        raise EstimateVariantAgeInputError("invalid_boolean_parameter")
    chance_parameters = None
    if chance_sharing:
        chance_parameters = {
            "median_allele_frequency": _number(parameters, "median_allele_frequency", 0.0),
            "markers_on_chromosome": _integer(parameters, "markers_on_chromosome", 0),
            "chromosome_length_cm": _number(parameters, "chromosome_length_cm", 0.0),
        }
    return {
        "method": method, "confidence_level": confidence_level,
        "minimum_units_for_estimate": _integer(parameters, "minimum_units_for_estimate", 3),
        "minimum_units_for_primary": _integer(parameters, "minimum_units_for_primary", 5),
        "generation_years": [float(value) for value in generation_years],
        "leave_one_family_out": leave_one_out,
        "chance_sharing_correction": chance_sharing,
        "chance_sharing_parameters": chance_parameters,
    }


def _artifact_by_id(stage_inputs: dict[str, Any], artifact_id: str) -> dict[str, Any]:
    matches = [artifact for artifact in stage_inputs["artifacts"] if artifact["artifact_id"] == artifact_id]
    if len(matches) != 1:
        raise EstimateVariantAgeInputError(f"{artifact_id}_missing_or_ambiguous")
    return matches[0]


def _validated_path(artifact: dict[str, Any], run_dir: Path) -> Path:
    artifact_path = Path(artifact["path"])
    if not artifact_path.is_absolute() and PurePosixPath(artifact["path"]).parts[:1] == ("stages",):
        artifact_path = run_dir / artifact_path
    elif not artifact_path.is_absolute():
        artifact_path = Path.cwd() / artifact_path
    if not artifact_path.is_file():
        raise EstimateVariantAgeInputError(f"artifact_missing:{artifact['artifact_id']}")
    if sha256_file(artifact_path) != artifact["sha256"]:
        raise EstimateVariantAgeInputError(f"artifact_checksum_mismatch:{artifact['artifact_id']}")
    return artifact_path


def execute(stage_inputs_path: Path, output_dir: Path) -> int:
    """Produit Gamma corrélé, indépendance et sensibilités sans outil externe."""
    started_at = utc_now()
    started_clock = monotonic()
    stage_inputs = read_json(stage_inputs_path)
    validate_json_document(stage_inputs, "stage_inputs.schema.json")
    run_dir = output_dir.parent.parent
    load_pipeline_config(run_dir / "config.resolved.yaml")
    parameters = _parameters(stage_inputs["parameters"])
    artifact_ids = ("founder_segments", "founder_consensus", "founder_analysis_summary")
    input_artifacts = {artifact_id: _artifact_by_id(stage_inputs, artifact_id) for artifact_id in artifact_ids}
    paths = {artifact_id: _validated_path(artifact, run_dir) for artifact_id, artifact in input_artifacts.items()}
    publication = publish_variant_age(
        founder_segments_path=paths["founder_segments"],
        founder_consensus_path=paths["founder_consensus"],
        founder_summary_path=paths["founder_analysis_summary"],
        output_dir=output_dir / "dating",
        confidence_level=parameters["confidence_level"],
        minimum_units_for_estimate=parameters["minimum_units_for_estimate"],
        minimum_units_for_primary=parameters["minimum_units_for_primary"],
        generation_years=parameters["generation_years"],
        leave_one_family_out=parameters["leave_one_family_out"],
        chance_sharing_parameters=parameters["chance_sharing_parameters"],
    )
    source = input_artifacts["founder_segments"]
    specifications = (
        ("variant_age_input", publication.input_path, "text/tab-separated-values", "variant_age_input.schema.json", "sensitive_genetic"),
        ("variant_age_estimates", publication.estimates_path, "text/tab-separated-values", "variant_age_estimates.schema.json", "internal"),
        ("variant_age_scenarios", publication.scenarios_path, "text/tab-separated-values", "variant_age_scenarios.schema.json", "sensitive_genetic"),
        ("variant_age_summary", publication.summary_path, "application/json", "variant_age_summary.schema.json", "internal"),
    )
    output_artifacts = [
        build_file_artifact(
            physical_path=path,
            published_path=f"{stage_inputs['published_output_dir']}/{path.relative_to(output_dir).as_posix()}",
            artifact_id=artifact_id, artifact_type=artifact_id, media_type=media_type,
            producer_stage=stage_inputs["stage_name"], producer_signature=stage_inputs["signature"],
            schema_name=schema_name, schema_version="1.0.0", assembly=source.get("assembly"),
            sample_set_id=source.get("sample_set_id"), variant_set_id=source.get("variant_set_id"),
            sensitivity=sensitivity,
        )
        for artifact_id, path, media_type, schema_name, sensitivity in specifications
    ]
    stage_outputs = {
        "schema_version": "1.0.0", "run_id": stage_inputs["run_id"],
        "stage_id": stage_inputs["stage_id"], "stage_name": stage_inputs["stage_name"],
        "signature": stage_inputs["signature"], "artifacts": output_artifacts,
    }
    validate_json_document(stage_outputs, "stage_outputs.schema.json")
    atomic_write_json(output_dir / "stage_outputs.json", stage_outputs)
    warning_code = None if publication.status == "PRIMARY_ESTIMATE" else publication.status.lower()
    audit = {
        "schema_version": "1.0.0", "run_id": stage_inputs["run_id"],
        "stage_id": stage_inputs["stage_id"], "stage_name": stage_inputs["stage_name"],
        "method_id": "gamma_gandolfo_2014_v1", "signature": stage_inputs["signature"],
        "started_at": started_at, "completed_at": utc_now(),
        "duration_seconds": monotonic() - started_clock,
        "inputs": list(input_artifacts.values()), "outputs": output_artifacts,
        "parameters": parameters,
        "tools": [
            {"tool": "python", "configured": Path(sys.executable).name, "version": sys.version.split()[0]},
            {"tool": "scipy", "configured": "Python package", "version": importlib.metadata.version("scipy")},
        ],
        "counts": {"included_independent_units": publication.included_unit_count, "sensitivity_scenarios": publication.scenario_count},
        "metrics": {"analysis_status": publication.status, "primary_model": "CORRELATED"},
        "exclusions": [],
        "warnings": [] if warning_code is None else [{"code": warning_code, "count": 1}],
        "checks": [
            {"check": "founder_artifact_integrity", "status": "PASS"},
            {"check": "individual_arm_lengths_cm", "status": "PASS"},
            {"check": "independent_family_units", "status": "PASS"},
            {"check": "gamma_exact_interval", "status": "PASS"},
            {"check": "non_conclusion_policy", "status": "PASS"},
        ],
        "known_limits": [
            "Gamma date le processus de recombinaison sous une origine fondatrice unique ; un segment IBS n'est pas une preuve IBD.",
            "L'intervalle Gamma conditionne sur les limites observées et ne couvre pas toute l'incertitude de phase, de carte ou de résolution des marqueurs.",
            "Le mode corrélé estime la dépendance depuis un petit nombre de longueurs et peut produire un intervalle très large.",
            "EstiAge n'est pas exécuté tant qu'un adaptateur local reproductible et respectueux des données privées n'est pas validé.",
        ],
        "expected_visualizations": ["plot_variant_age_lengths", "plot_variant_age_forest", "plot_variant_age_sensitivity"],
        "manual_validation_required": True,
    }
    validate_json_document(audit, "stage_audit.schema.json")
    atomic_write_json(output_dir / "audit.json", audit)
    (output_dir / "checksums.sha256").write_text(
        "".join(
            f"{artifact['sha256']}  {artifact['path'].removeprefix(stage_inputs['published_output_dir'] + '/')}\n"
            for artifact in output_artifacts
        ), encoding="utf-8",
    )
    return 0


def _classify_error(error: Exception) -> int:
    if isinstance(error, GammaAgeError):
        return 4
    if isinstance(error, (EstimateVariantAgeInputError, DocumentValidationError, OSError, ValueError, json.JSONDecodeError)):
        return 2
    return 5


def main(arguments: Sequence[str] | None = None) -> int:
    """Point d'entrée autonome respectant les codes de retour V2."""
    parsed_arguments = _build_parser().parse_args(arguments)
    try:
        return execute(parsed_arguments.stage_inputs, parsed_arguments.output_dir)
    except Exception as error:
        sys.stderr.write(f"{error}\n")
        return _classify_error(error)


if __name__ == "__main__":
    raise SystemExit(main())
