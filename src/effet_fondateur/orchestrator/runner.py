"""Exécution isolée, validation et publication des scripts d'étape."""

from __future__ import annotations

import os
import subprocess
import sys
from dataclasses import dataclass
from pathlib import Path
from time import monotonic
from typing import Any
from uuid import uuid4

from effet_fondateur.audit import atomic_write_json, sha256_file
from effet_fondateur.contracts import (
    DocumentValidationError,
    load_pipeline_config,
    validate_json_document,
)
from effet_fondateur.orchestrator.errors import IntegrityError, StageExecutionError
from effet_fondateur.orchestrator.integrity import (
    validate_attempt_outputs,
    validate_published_stage,
)
from effet_fondateur.orchestrator.inputs import resolve_stage_input_artifacts
from effet_fondateur.orchestrator.models import StageDefinition
from effet_fondateur.orchestrator.signatures import build_stage_signature
from effet_fondateur.orchestrator.state import (
    find_stage_record,
    load_manifest,
    record_event,
    save_manifest,
    utc_now,
)


@dataclass(frozen=True)
class _StageAttempt:
    """Chemins et signature immuables d'une tentative d'étape."""

    signature: str
    attempt_dir: Path
    failed_attempt_dir: Path
    final_stage_dir: Path
    started_clock: float


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


def _record_failure(
    run_dir: Path,
    manifest: dict[str, Any],
    stage_record: dict[str, Any],
    definition: StageDefinition,
    return_code: int,
    started_clock: float,
) -> None:
    """Place une étape en échec sans recopier de sortie sensible au manifest."""
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


def _reuse_published_stage(
    run_dir: Path,
    manifest: dict[str, Any],
    stage_record: dict[str, Any],
    definition: StageDefinition,
    current_signature: str,
) -> bool:
    """Réutilise une étape terminale uniquement si son intégrité est intacte."""
    if stage_record["state"] not in {"SUCCEEDED", "CACHED"}:
        return False
    if stage_record["signature"] != current_signature:
        stage_record["state"] = "FAILED"
        stage_record["last_error_code"] = 5
        manifest["global_status"] = "BLOCKED"
        save_manifest(run_dir, manifest)
        record_event(
            run_dir,
            definition.directory_name,
            "stage_input_integrity_failed",
            severity="ERROR",
        )
        raise IntegrityError(
            f"Entrées modifiées depuis la publication de {definition.stage_name}."
        )
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
    return True


def _prepare_attempt(
    run_dir: Path,
    manifest: dict[str, Any],
    stage_record: dict[str, Any],
    definition: StageDefinition,
    parameters: dict[str, Any],
    signature: str,
    input_artifacts: list[dict[str, Any]],
) -> _StageAttempt:
    """Crée les entrées temporaires et publie l'état RUNNING dans le manifest."""
    attempt_identifier = uuid4().hex
    attempt_dir = run_dir / "stages" / f".{definition.directory_name}.{attempt_identifier}.tmp"
    failed_attempt_dir = (
        run_dir
        / "stages"
        / "attempts"
        / f"{definition.directory_name}.{attempt_identifier}.failed"
    )
    final_stage_dir = run_dir / "stages" / definition.directory_name
    attempt_dir.mkdir(parents=True)
    attempt_number = stage_record["attempt_count"] + 1
    stage_inputs = {
        "schema_version": "1.0.0",
        "run_id": manifest["run_id"],
        "stage_id": definition.stage_id,
        "stage_name": definition.stage_name,
        "signature": signature,
        "attempt_number": attempt_number,
        "published_output_dir": f"stages/{definition.directory_name}",
        "parameters": parameters,
        "artifacts": input_artifacts,
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
            "attempt_count": attempt_number,
            "last_error_code": None,
        }
    )
    manifest["global_status"] = "INCOMPLETE"
    save_manifest(run_dir, manifest)
    record_event(run_dir, definition.directory_name, "stage_started")
    return _StageAttempt(
        signature=signature,
        attempt_dir=attempt_dir,
        failed_attempt_dir=failed_attempt_dir,
        final_stage_dir=final_stage_dir,
        started_clock=started_clock,
    )


def _execute_attempt(
    attempt: _StageAttempt,
    definition: StageDefinition,
    run_id: str,
) -> None:
    """Lance le sous-processus, valide ses sorties et publie son dossier."""
    command = [
        sys.executable,
        "-m",
        definition.module,
        "--stage-inputs",
        str(attempt.attempt_dir / "stage_inputs.json"),
        "--output-dir",
        str(attempt.attempt_dir),
    ]
    atomic_write_json(
        attempt.attempt_dir / "command.json",
        {
            "schema_version": "1.0.0",
            "executable": Path(sys.executable).name,
            "module": definition.module,
            "arguments": ["--stage-inputs", "stage_inputs.json", "--output-dir", "."],
            "shell": False,
        },
    )
    # Le descripteur publié masque les chemins absolus de la machine, tandis que
    # la commande réelle conserve des arguments séparés et n'utilise jamais un shell.
    try:
        completed_process = subprocess.run(
            command,
            capture_output=True,
            check=False,
            text=True,
        )
        (attempt.attempt_dir / "stdout.log").write_text(
            completed_process.stdout,
            encoding="utf-8",
        )
        (attempt.attempt_dir / "stderr.log").write_text(
            completed_process.stderr,
            encoding="utf-8",
        )
        if completed_process.returncode != 0:
            raise StageExecutionError(
                f"Échec de l'étape {definition.stage_name} "
                f"(code {completed_process.returncode}).",
                completed_process.returncode,
            )
        validate_attempt_outputs(
            attempt.attempt_dir,
            definition,
            run_id,
            attempt.signature,
        )
        os.replace(attempt.attempt_dir, attempt.final_stage_dir)
    except (DocumentValidationError, IntegrityError, OSError, ValueError) as error:
        raise StageExecutionError(
            f"Sorties invalides pour l'étape {definition.stage_name} : {error}",
            5,
        ) from error
    finally:
        # Une tentative invalide n'est jamais publiée, mais ses journaux restent
        # disponibles pour expliquer l'échec et préparer une reprise identique.
        if attempt.attempt_dir.exists():
            attempt.failed_attempt_dir.parent.mkdir(parents=True, exist_ok=True)
            os.replace(attempt.attempt_dir, attempt.failed_attempt_dir)


def _record_success(
    run_dir: Path,
    manifest: dict[str, Any],
    stage_record: dict[str, Any],
    definition: StageDefinition,
    attempt: _StageAttempt,
) -> None:
    """Enregistre les empreintes de provenance après publication atomique."""
    stage_record.update(
        {
            "state": "SUCCEEDED",
            "completed_at": utc_now(),
            "duration_seconds": monotonic() - attempt.started_clock,
            "audit_path": f"stages/{definition.directory_name}/audit.json",
            "audit_sha256": sha256_file(attempt.final_stage_dir / "audit.json"),
            "stage_outputs_sha256": sha256_file(
                attempt.final_stage_dir / "stage_outputs.json"
            ),
            "last_error_code": None,
        }
    )
    save_manifest(run_dir, manifest)
    record_event(run_dir, definition.directory_name, "stage_succeeded")


def run_stage(
    run_dir: Path,
    definition: StageDefinition,
    parameters: dict[str, Any],
) -> None:
    """Orchestre une tentative sans exécuter de logique scientifique en interne."""
    manifest = load_manifest(run_dir)
    _validate_dependencies(manifest, definition)
    stage_record = find_stage_record(manifest, definition.stage_name)
    if stage_record is None:
        stage_record = _new_stage_record(definition)
        manifest["stages"].append(stage_record)
        save_manifest(run_dir, manifest)
    config = load_pipeline_config(run_dir / "config.resolved.yaml")
    input_artifacts = resolve_stage_input_artifacts(
        config,
        definition,
        manifest["config_sha256"],
    )
    signature = build_stage_signature(
        definition,
        parameters,
        manifest["config_sha256"],
        input_artifacts,
    )
    if _reuse_published_stage(
        run_dir,
        manifest,
        stage_record,
        definition,
        signature,
    ):
        return
    final_stage_dir = run_dir / "stages" / definition.directory_name
    if final_stage_dir.exists():
        raise IntegrityError(f"Dossier publié inattendu : {final_stage_dir}")

    attempt = _prepare_attempt(
        run_dir,
        manifest,
        stage_record,
        definition,
        parameters,
        signature,
        input_artifacts,
    )
    _execute_attempt(attempt, definition, manifest["run_id"])
    _record_success(run_dir, manifest, stage_record, definition, attempt)


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
        if stage_record is not None and stage_record["state"] in {"PENDING", "RUNNING"}:
            _record_failure(
                run_dir,
                manifest,
                stage_record,
                definition,
                error.return_code,
                started_clock,
            )
        raise
