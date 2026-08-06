"""Étape 13 : haplotype fondateur candidat et partage IBS local."""

from __future__ import annotations

import argparse
import json
import subprocess
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
    validate_tsv_table,
)
from effet_fondateur.founder import FounderAnalysisError, infer_target_centered_ibs
from effet_fondateur.orchestrator.state import utc_now


class InferFounderInputError(ValueError):
    """Signale une configuration ou une entrée invalide de l'étape 13."""


def _build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(prog="infer-founder-haplotype")
    parser.add_argument("--stage-inputs", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    return parser


def _positive_integer(parameters: dict[str, Any], name: str, default: int) -> int:
    value = parameters.get(name, default)
    if isinstance(value, bool) or not isinstance(value, int) or value <= 0:
        raise InferFounderInputError(f"invalid_parameter:{name}")
    return value


def _parameters(parameters: dict[str, Any]) -> dict[str, Any]:
    method = parameters.get("segment_method", "target_centered_exact_ibs_v1")
    if method != "target_centered_exact_ibs_v1":
        raise InferFounderInputError("invalid_parameter:segment_method")
    return {
        "segment_method": method,
        "minimum_independent_carriers": _positive_integer(
            parameters, "minimum_independent_carriers", 3
        ),
        "minimum_flank_markers": _positive_integer(
            parameters, "minimum_flank_markers", 2
        ),
        "bcftools_timeout_seconds": _positive_integer(
            parameters, "bcftools_timeout_seconds", 300
        ),
    }


def _artifact_by_id(stage_inputs: dict[str, Any], artifact_id: str) -> dict[str, Any]:
    matches = [
        artifact for artifact in stage_inputs["artifacts"]
        if artifact["artifact_id"] == artifact_id
    ]
    if len(matches) != 1:
        raise InferFounderInputError(f"{artifact_id}_missing_or_ambiguous")
    return matches[0]


def _resolve_path(artifact: dict[str, Any], run_dir: Path) -> Path:
    artifact_path = Path(artifact["path"])
    if artifact_path.is_absolute():
        return artifact_path
    if PurePosixPath(artifact["path"]).parts[:1] == ("stages",):
        return run_dir / artifact_path
    return Path.cwd() / artifact_path


def _validated_path(artifact: dict[str, Any], run_dir: Path) -> Path:
    path = _resolve_path(artifact, run_dir)
    if not path.is_file():
        raise InferFounderInputError(f"artifact_missing:{artifact['artifact_id']}")
    if sha256_file(path) != artifact["sha256"]:
        raise InferFounderInputError(f"artifact_checksum_mismatch:{artifact['artifact_id']}")
    return path


def _bcftools_version(command: str, timeout_seconds: int) -> str:
    try:
        completed = subprocess.run(
            [command, "--version"], check=False, capture_output=True, text=True,
            timeout=timeout_seconds,
        )
    except (OSError, subprocess.TimeoutExpired) as error:
        raise FounderAnalysisError("bcftools_version_unavailable") from error
    if completed.returncode != 0 or not completed.stdout.strip():
        raise FounderAnalysisError("bcftools_version_unavailable")
    return completed.stdout.splitlines()[0].strip()


def _artifact(
    *, path: Path, output_dir: Path, stage_inputs: dict[str, Any], artifact_id: str,
    media_type: str, schema_name: str, source: dict[str, Any], sensitivity: str,
) -> dict[str, Any]:
    return build_file_artifact(
        physical_path=path,
        published_path=f"{stage_inputs['published_output_dir']}/{path.relative_to(output_dir).as_posix()}",
        artifact_id=artifact_id, artifact_type=artifact_id, media_type=media_type,
        producer_stage=stage_inputs["stage_name"], producer_signature=stage_inputs["signature"],
        schema_name=schema_name, schema_version="1.0.0", assembly=source.get("assembly"),
        sample_set_id=source.get("sample_set_id"), variant_set_id=source.get("variant_set_id"),
        sensitivity=sensitivity,
    )


def execute(stage_inputs_path: Path, output_dir: Path) -> int:
    """Exécute l'analyse IBS locale et publie un audit sans conclure automatiquement à l'IBD."""
    started_at = utc_now()
    started_clock = monotonic()
    stage_inputs = read_json(stage_inputs_path)
    validate_json_document(stage_inputs, "stage_inputs.schema.json")
    run_dir = output_dir.parent.parent
    config = load_pipeline_config(run_dir / "config.resolved.yaml")
    parameters = _parameters(stage_inputs["parameters"])
    if config["tools"]["local_ibd_adapter"] is not None:
        raise InferFounderInputError("local_ibd_adapter_not_supported_by_exact_ibs_method")
    bcftools_command = config["tools"]["bcftools"]
    if not isinstance(bcftools_command, str) or not bcftools_command:
        raise InferFounderInputError("bcftools_not_configured")
    artifact_ids = (
        "shapeit5_final_bcf", "shapeit5_final_index", "carrier_haplotypes",
        "target_genetic_map", "cohorts_frozen", "samples_master",
    )
    input_artifacts = {
        artifact_id: _artifact_by_id(stage_inputs, artifact_id)
        for artifact_id in artifact_ids
    }
    paths = {
        artifact_id: _validated_path(artifact, run_dir)
        for artifact_id, artifact in input_artifacts.items()
    }
    bcftools_version = _bcftools_version(
        bcftools_command, parameters["bcftools_timeout_seconds"]
    )
    result = infer_target_centered_ibs(
        phased_bcf_path=paths["shapeit5_final_bcf"],
        carrier_haplotypes_path=paths["carrier_haplotypes"],
        cohorts_path=paths["cohorts_frozen"], samples_master_path=paths["samples_master"],
        genetic_map_path=paths["target_genetic_map"], output_dir=output_dir / "founder",
        bcftools_command=bcftools_command,
        timeout_seconds=parameters["bcftools_timeout_seconds"],
        minimum_independent_carriers=parameters["minimum_independent_carriers"],
        minimum_flank_markers=parameters["minimum_flank_markers"],
    )
    validate_tsv_table(result.segments_path, "founder_segments.schema.json")
    validate_tsv_table(result.consensus_path, "founder_consensus.schema.json")
    validate_tsv_table(result.sharing_path, "founder_sharing_matrix.schema.json")
    validate_tsv_table(result.discordances_path, "founder_discordances.schema.json")
    validate_json_document(read_json(result.summary_path), "founder_analysis_summary.schema.json")
    source = input_artifacts["shapeit5_final_bcf"]
    specifications = (
        ("founder_segments", result.segments_path, "text/tab-separated-values", "founder_segments.schema.json", "sensitive_genetic"),
        ("founder_consensus", result.consensus_path, "text/tab-separated-values", "founder_consensus.schema.json", "sensitive_genetic"),
        ("founder_sharing_matrix", result.sharing_path, "text/tab-separated-values", "founder_sharing_matrix.schema.json", "sensitive_genetic"),
        ("founder_discordances", result.discordances_path, "text/tab-separated-values", "founder_discordances.schema.json", "sensitive_genetic"),
        ("founder_analysis_summary", result.summary_path, "application/json", "founder_analysis_summary.schema.json", "internal"),
    )
    output_artifacts = [
        _artifact(
            path=path, output_dir=output_dir, stage_inputs=stage_inputs,
            artifact_id=artifact_id, media_type=media_type, schema_name=schema_name,
            source=source, sensitivity=sensitivity,
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
    warning_code = None if result.status == "SUPPORTED_IBS_CANDIDATE" else result.status.lower()
    audit = {
        "schema_version": "1.0.0", "run_id": stage_inputs["run_id"],
        "stage_id": stage_inputs["stage_id"], "stage_name": stage_inputs["stage_name"],
        "method_id": "target_centered_exact_ibs_v1", "signature": stage_inputs["signature"],
        "started_at": started_at, "completed_at": utc_now(),
        "duration_seconds": monotonic() - started_clock,
        "inputs": list(input_artifacts.values()), "outputs": output_artifacts,
        "parameters": parameters,
        "tools": [{"tool": "bcftools", "configured": bcftools_command, "version": bcftools_version}],
        "counts": {
            "selected_carrier_units": result.selected_carrier_count,
            "background_haplotypes": result.background_haplotype_count,
            "matching_background_haplotypes": result.matching_background_haplotype_count,
        },
        "metrics": {"consensus_status": result.status, "ibd_claimed": False},
        "exclusions": [],
        "warnings": [] if warning_code is None else [{"code": warning_code, "count": 1}],
        "checks": [
            {"check": "input_artifact_integrity", "status": "PASS"},
            {"check": "phased_sample_identity", "status": "PASS"},
            {"check": "target_carrier_assignment", "status": "PASS"},
            {"check": "ibs_ibd_distinction", "status": "PASS"},
            {"check": "founder_output_contracts", "status": "PASS"},
        ],
        "known_limits": [
            "La concordance allélique démontre un partage IBS, pas une transmission IBD.",
            "La méthode exacte s'arrête au premier génotype manquant ou discordant.",
            "La résolution des limites dépend de la densité des marqueurs et de la carte génétique.",
            "Le variant cible est exclu de la signature comparée aux témoins et non-porteurs.",
        ],
        "expected_visualizations": [
            "plot_founder_haplotypes", "plot_founder_segments",
            "plot_founder_sharing_matrix", "plot_founder_background_frequency",
        ],
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
    if isinstance(error, (InferFounderInputError, FounderAnalysisError, DocumentValidationError, OSError, ValueError, json.JSONDecodeError)):
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
