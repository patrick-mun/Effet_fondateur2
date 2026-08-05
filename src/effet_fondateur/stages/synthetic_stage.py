"""Étape sans donnée génétique utilisée pour valider l'orchestration V2."""

from __future__ import annotations

import argparse
import json
from pathlib import Path, PurePosixPath
from time import monotonic
from typing import Any, Sequence

from effet_fondateur.audit import atomic_write_json, read_json, sha256_file
from effet_fondateur.contracts import DocumentValidationError, validate_json_document
from effet_fondateur.orchestrator.state import utc_now


def _build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(prog="synthetic-stage")
    parser.add_argument("--stage-inputs", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    return parser


def _artifact(
    stage_inputs: dict[str, Any], output_dir: Path, filename: str
) -> dict[str, Any]:
    published_path = str(
        PurePosixPath(stage_inputs["published_output_dir"]) / filename
    )
    return {
        "artifact_id": "synthetic_result",
        "artifact_type": "synthetic_json",
        "path": published_path,
        "media_type": "application/json",
        "schema_name": None,
        "schema_version": None,
        "sha256": sha256_file(output_dir / filename),
        "producer_stage": f"{stage_inputs['stage_id']}_{stage_inputs['stage_name']}",
        "producer_signature": stage_inputs["signature"],
        "assembly": None,
        "sample_set_id": None,
        "variant_set_id": None,
        "sensitivity": "internal",
    }


def execute(stage_inputs_path: Path, output_dir: Path) -> int:
    """Produit un artefact déterministe ou un échec contrôlé demandé par le test."""
    started_at = utc_now()
    started_clock = monotonic()
    stage_inputs = read_json(stage_inputs_path)
    validate_json_document(stage_inputs, "stage_inputs.schema.json")
    parameters = stage_inputs["parameters"]
    fail_attempts = int(parameters.get("fail_attempts", 0))
    if stage_inputs["attempt_number"] <= fail_attempts:
        return 4

    output_dir.mkdir(parents=True, exist_ok=True)
    result = {
        "schema_version": "1.0.0",
        "message": str(parameters.get("message", "synthetic-stage-ok")),
    }
    atomic_write_json(output_dir / "synthetic_result.json", result)
    artifact = _artifact(stage_inputs, output_dir, "synthetic_result.json")
    stage_outputs = {
        "schema_version": "1.0.0",
        "run_id": stage_inputs["run_id"],
        "stage_id": stage_inputs["stage_id"],
        "stage_name": stage_inputs["stage_name"],
        "signature": stage_inputs["signature"],
        "artifacts": [artifact],
    }
    validate_json_document(stage_outputs, "stage_outputs.schema.json")
    atomic_write_json(output_dir / "stage_outputs.json", stage_outputs)

    audit = {
        "schema_version": "1.0.0",
        "run_id": stage_inputs["run_id"],
        "stage_id": stage_inputs["stage_id"],
        "stage_name": stage_inputs["stage_name"],
        "method_id": "synthetic_stage_v1",
        "signature": stage_inputs["signature"],
        "started_at": started_at,
        "completed_at": utc_now(),
        "duration_seconds": monotonic() - started_clock,
        "inputs": [],
        "outputs": [artifact],
        "parameters": parameters,
        "tools": [],
        "counts": {"artifacts": 1},
        "metrics": {},
        "exclusions": [],
        "warnings": [],
        "checks": [{"check": "synthetic_output", "status": "PASS"}],
        "known_limits": ["Étape technique sans calcul scientifique."],
        "expected_visualizations": [],
        "manual_validation_required": False,
    }
    validate_json_document(audit, "stage_audit.schema.json")
    atomic_write_json(output_dir / "audit.json", audit)
    (output_dir / "checksums.sha256").write_text(
        f"{artifact['sha256']}  synthetic_result.json\n",
        encoding="utf-8",
    )
    return 0


def main(arguments: Sequence[str] | None = None) -> int:
    """Point d'entrée autonome respectant les codes de retour V2."""
    parser = _build_parser()
    parsed_arguments = parser.parse_args(arguments)
    try:
        return execute(parsed_arguments.stage_inputs, parsed_arguments.output_dir)
    except (DocumentValidationError, ValueError, OSError, json.JSONDecodeError) as error:
        parser.error(str(error))
    return 2


if __name__ == "__main__":
    raise SystemExit(main())
