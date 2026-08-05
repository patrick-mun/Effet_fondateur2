"""Exécution isolée, validation et publication des scripts d'étape."""

from __future__ import annotations

import hashlib
import importlib.util
import json
import os
import subprocess
import sys
from pathlib import Path, PurePosixPath
from time import monotonic
from typing import Any
from uuid import uuid4

from effet_fondateur.audit import atomic_write_json, read_json, sha256_file
from effet_fondateur.contracts import DocumentValidationError, validate_json_document
from effet_fondateur.orchestrator.errors import IntegrityError, StageExecutionError
from effet_fondateur.orchestrator.models import StageDefinition
from effet_fondateur.orchestrator.state import (
    find_stage_record,
    load_manifest,
    record_event,
    save_manifest,
    utc_now,
)


def _module_sha256(module_name: str) -> str:
    module_specification = importlib.util.find_spec(module_name)
    if module_specification is None or module_specification.origin is None:
        raise StageExecutionError(f"Module d'étape introuvable : {module_name}", 2)
    return sha256_file(Path(module_specification.origin))


def _stage_signature(
    definition: StageDefinition,
    parameters: dict[str, Any],
    config_sha256: str,
) -> str:
    signature_payload = {
        "stage_id": definition.stage_id,
        "stage_name": definition.stage_name,
        "module_sha256": _module_sha256(definition.module),
        "parameters": parameters,
        "config_sha256": config_sha256,
        "input_schema_version": "1.0.0",
        "output_schema_version": "1.0.0",
        "audit_schema_version": "1.0.0",
    }
    serialized_payload = json.dumps(signature_payload, sort_keys=True, separators=(",", ":"))
    return hashlib.sha256(serialized_payload.encode("utf-8")).hexdigest()


def _new_stage_record(definition: StageDefinition) -> dict[str, Any]:
    return {
        "stage_id": definition.stage_id,
        "stage_name": definition.stage_name,
        "state": "PENDING",
        "critical": definition.critical,
        "signature": None,
        "started_at": None,
        "completed_at": None,
        "duration_seconds": None,
        "audit_path": None,
        "audit_sha256": None,
        "stage_outputs_sha256": None,
        "attempt_count": 0,
        "last_error_code": None,
    }


def _validate_dependencies(
    manifest: dict[str, Any], definition: StageDefinition
) -> None:
    states_by_name = {
        stage_record["stage_name"]: stage_record["state"]
        for stage_record in manifest["stages"]
    }
    invalid_dependencies = [
        dependency
        for dependency in definition.dependencies
        if states_by_name.get(dependency) not in {"SUCCEEDED", "CACHED"}
    ]
    if invalid_dependencies:
        raise StageExecutionError(
            f"Étape {definition.stage_name} bloquée par : "
            f"{', '.join(invalid_dependencies)}",
            4,
        )


def _artifact_path_in_attempt(
    artifact_path: str,
    published_stage_dir: PurePosixPath,
    attempt_dir: Path,
) -> Path:
    declared_path = PurePosixPath(artifact_path)
    try:
        relative_path = declared_path.relative_to(published_stage_dir)
    except ValueError as error:
        raise IntegrityError(
            f"Artefact déclaré hors du dossier d'étape : {artifact_path}"
        ) from error
    if not relative_path.parts or ".." in relative_path.parts:
        raise IntegrityError(f"Chemin d'artefact non sûr : {artifact_path}")
    return attempt_dir.joinpath(*relative_path.parts)


def _validate_attempt_outputs(
    attempt_dir: Path,
    definition: StageDefinition,
    run_id: str,
    signature: str,
) -> dict[str, Any]:
    stage_outputs = read_json(attempt_dir / "stage_outputs.json")
    validate_json_document(stage_outputs, "stage_outputs.schema.json")
    expected_values = {
        "run_id": run_id,
        "stage_id": definition.stage_id,
        "stage_name": definition.stage_name,
        "signature": signature,
    }
    for field, expected_value in expected_values.items():
        if stage_outputs[field] != expected_value:
            raise IntegrityError(
                f"Sortie incohérente pour {definition.stage_name}, champ {field}."
            )

    published_stage_dir = PurePosixPath("stages") / definition.directory_name
    artifact_ids: set[str] = set()
    for artifact in stage_outputs["artifacts"]:
        if artifact["artifact_id"] in artifact_ids:
            raise IntegrityError(
                f"Identifiant d'artefact dupliqué : {artifact['artifact_id']}"
            )
        artifact_ids.add(artifact["artifact_id"])
        artifact_path = _artifact_path_in_attempt(
            artifact["path"], published_stage_dir, attempt_dir
        )
        if not artifact_path.is_file():
            raise IntegrityError(f"Artefact annoncé mais absent : {artifact['path']}")
        if sha256_file(artifact_path) != artifact["sha256"]:
            raise IntegrityError(f"Empreinte invalide : {artifact['path']}")

    audit = read_json(attempt_dir / "audit.json")
    validate_json_document(audit, "stage_audit.schema.json")
    for field, expected_value in expected_values.items():
        if audit[field] != expected_value:
            raise IntegrityError(
                f"Audit incohérent pour {definition.stage_name}, champ {field}."
            )
    return stage_outputs


def validate_published_stage(
    run_dir: Path,
    definition: StageDefinition,
    stage_record: dict[str, Any],
) -> None:
    """Recalcule les empreintes d'une étape publiée avant toute réutilisation."""
    stage_dir = run_dir / "stages" / definition.directory_name
    stage_outputs_path = stage_dir / "stage_outputs.json"
    audit_path = stage_dir / "audit.json"
    if sha256_file(stage_outputs_path) != stage_record["stage_outputs_sha256"]:
        raise IntegrityError(f"Descripteur de sorties modifié : {stage_outputs_path}")
    if sha256_file(audit_path) != stage_record["audit_sha256"]:
        raise IntegrityError(f"Audit d'étape modifié : {audit_path}")

    stage_outputs = read_json(stage_outputs_path)
    validate_json_document(stage_outputs, "stage_outputs.schema.json")
    expected_values = {
        "run_id": load_manifest(run_dir)["run_id"],
        "stage_id": definition.stage_id,
        "stage_name": definition.stage_name,
        "signature": stage_record["signature"],
    }
    for field, expected_value in expected_values.items():
        if stage_outputs[field] != expected_value:
            raise IntegrityError(
                f"Sortie publiée incohérente pour {definition.stage_name}, champ {field}."
            )
    published_stage_dir = PurePosixPath("stages") / definition.directory_name
    for artifact in stage_outputs["artifacts"]:
        artifact_path = _artifact_path_in_attempt(
            artifact["path"], published_stage_dir, stage_dir
        )
        if not artifact_path.is_file():
            raise IntegrityError(f"Artefact publié absent : {artifact['path']}")
        if sha256_file(artifact_path) != artifact["sha256"]:
            raise IntegrityError(f"Artefact publié modifié : {artifact['path']}")
    audit = read_json(audit_path)
    validate_json_document(audit, "stage_audit.schema.json")
    for field, expected_value in expected_values.items():
        if audit[field] != expected_value:
            raise IntegrityError(
                f"Audit publié incohérent pour {definition.stage_name}, champ {field}."
            )


def _record_failure(
    run_dir: Path,
    manifest: dict[str, Any],
    stage_record: dict[str, Any],
    definition: StageDefinition,
    return_code: int,
    started_clock: float,
) -> None:
    stage_record["state"] = "FAILED"
    stage_record["completed_at"] = utc_now()
    stage_record["duration_seconds"] = monotonic() - started_clock
    stage_record["last_error_code"] = return_code
    manifest["global_status"] = "BLOCKED" if definition.critical else "INCOMPLETE"
    save_manifest(run_dir, manifest)
    record_event(
        run_dir,
        definition.directory_name,
        "stage_failed",
        severity="ERROR",
        details={"return_code": return_code},
    )


def run_stage(
    run_dir: Path,
    definition: StageDefinition,
    parameters: dict[str, Any],
) -> None:
    """Lance une étape en sous-processus puis publie uniquement des sorties valides."""
    manifest = load_manifest(run_dir)
    _validate_dependencies(manifest, definition)
    stage_record = find_stage_record(manifest, definition.stage_name)
    if stage_record is None:
        stage_record = _new_stage_record(definition)
        manifest["stages"].append(stage_record)

    final_stage_dir = run_dir / "stages" / definition.directory_name
    if stage_record["state"] in {"SUCCEEDED", "CACHED"}:
        try:
            validate_published_stage(run_dir, definition, stage_record)
        except IntegrityError:
            stage_record["state"] = "FAILED"
            stage_record["last_error_code"] = 5
            manifest["global_status"] = "BLOCKED"
            save_manifest(run_dir, manifest)
            record_event(
                run_dir,
                definition.directory_name,
                "published_artifact_integrity_failed",
                severity="ERROR",
            )
            raise
        record_event(run_dir, definition.directory_name, "stage_reused")
        return
    if final_stage_dir.exists():
        raise IntegrityError(f"Dossier publié inattendu : {final_stage_dir}")

    signature = _stage_signature(definition, parameters, manifest["config_sha256"])
    attempt_identifier = uuid4().hex
    attempt_dir = run_dir / "stages" / f".{definition.directory_name}.{attempt_identifier}.tmp"
    failed_attempt_dir = (
        run_dir / "stages" / "attempts" / f"{definition.directory_name}.{attempt_identifier}.failed"
    )
    attempt_dir.mkdir(parents=True)
    stage_inputs = {
        "schema_version": "1.0.0",
        "run_id": manifest["run_id"],
        "stage_id": definition.stage_id,
        "stage_name": definition.stage_name,
        "signature": signature,
        "attempt_number": stage_record["attempt_count"] + 1,
        "published_output_dir": f"stages/{definition.directory_name}",
        "parameters": parameters,
        "artifacts": [],
    }
    validate_json_document(stage_inputs, "stage_inputs.schema.json")
    atomic_write_json(attempt_dir / "stage_inputs.json", stage_inputs)

    started_clock = monotonic()
    stage_record.update(
        {
            "state": "RUNNING",
            "signature": signature,
            "started_at": utc_now(),
            "completed_at": None,
            "duration_seconds": None,
            "audit_path": None,
            "audit_sha256": None,
            "stage_outputs_sha256": None,
            "attempt_count": stage_record["attempt_count"] + 1,
            "last_error_code": None,
        }
    )
    manifest["global_status"] = "INCOMPLETE"
    save_manifest(run_dir, manifest)
    record_event(run_dir, definition.directory_name, "stage_started")

    command = [
        sys.executable,
        "-m",
        definition.module,
        "--stage-inputs",
        str(attempt_dir / "stage_inputs.json"),
        "--output-dir",
        str(attempt_dir),
    ]
    atomic_write_json(
        attempt_dir / "command.json",
        {
            "schema_version": "1.0.0",
            "executable": Path(sys.executable).name,
            "module": definition.module,
            "arguments": [
                "--stage-inputs",
                "stage_inputs.json",
                "--output-dir",
                ".",
            ],
            "shell": False,
        },
    )
    try:
        completed_process = subprocess.run(
            command,
            capture_output=True,
            check=False,
            text=True,
        )
        (attempt_dir / "stdout.log").write_text(
            completed_process.stdout, encoding="utf-8"
        )
        (attempt_dir / "stderr.log").write_text(
            completed_process.stderr, encoding="utf-8"
        )
        if completed_process.returncode != 0:
            raise StageExecutionError(
                f"Échec de l'étape {definition.stage_name} "
                f"(code {completed_process.returncode}).",
                completed_process.returncode,
            )
        _validate_attempt_outputs(
            attempt_dir,
            definition,
            manifest["run_id"],
            signature,
        )
        os.replace(attempt_dir, final_stage_dir)
    except (DocumentValidationError, IntegrityError, OSError, ValueError) as error:
        raise StageExecutionError(
            f"Sorties invalides pour l'étape {definition.stage_name} : {error}", 5
        ) from error
    except StageExecutionError:
        raise
    finally:
        if attempt_dir.exists():
            failed_attempt_dir.parent.mkdir(parents=True, exist_ok=True)
            os.replace(attempt_dir, failed_attempt_dir)

    stage_record.update(
        {
            "state": "SUCCEEDED",
            "completed_at": utc_now(),
            "duration_seconds": monotonic() - started_clock,
            "audit_path": f"stages/{definition.directory_name}/audit.json",
            "audit_sha256": sha256_file(final_stage_dir / "audit.json"),
            "stage_outputs_sha256": sha256_file(
                final_stage_dir / "stage_outputs.json"
            ),
            "last_error_code": None,
        }
    )
    save_manifest(run_dir, manifest)
    record_event(run_dir, definition.directory_name, "stage_succeeded")


def run_stage_with_failure_audit(
    run_dir: Path,
    definition: StageDefinition,
    parameters: dict[str, Any],
) -> None:
    """Exécute une étape et garantit la mise à jour du manifest en cas d'échec."""
    started_clock = monotonic()
    try:
        run_stage(run_dir, definition, parameters)
    except StageExecutionError as error:
        manifest = load_manifest(run_dir)
        stage_record = find_stage_record(manifest, definition.stage_name)
        if stage_record is not None and stage_record["state"] == "RUNNING":
            _record_failure(
                run_dir,
                manifest,
                stage_record,
                definition,
                error.return_code,
                started_clock,
            )
        raise
