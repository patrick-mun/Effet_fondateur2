"""Contrôles d'intégrité avant publication et avant réutilisation."""

from __future__ import annotations

from pathlib import Path, PurePosixPath
from typing import Any

from effet_fondateur.audit import read_json, sha256_file
from effet_fondateur.contracts import validate_json_document
from effet_fondateur.orchestrator.errors import IntegrityError
from effet_fondateur.orchestrator.models import StageDefinition
from effet_fondateur.orchestrator.state import load_manifest


def resolve_declared_artifact_path(
    artifact_path: str,
    published_stage_dir: PurePosixPath,
    physical_stage_dir: Path,
) -> Path:
    """Traduit un chemin publié en chemin physique sans autoriser d'évasion."""
    declared_path = PurePosixPath(artifact_path)
    try:
        relative_path = declared_path.relative_to(published_stage_dir)
    except ValueError as error:
        raise IntegrityError(
            f"Artefact déclaré hors du dossier d'étape : {artifact_path}"
        ) from error
    # Une étape ne peut déclarer que ses propres sorties. Cette frontière évite
    # qu'un script déclare puis fasse publier un fichier d'un autre dossier du run.
    if not relative_path.parts or ".." in relative_path.parts:
        raise IntegrityError(f"Chemin d'artefact non sûr : {artifact_path}")
    return physical_stage_dir.joinpath(*relative_path.parts)


def _expected_identity(
    definition: StageDefinition,
    run_id: str,
    signature: str,
) -> dict[str, str]:
    return {
        "run_id": run_id,
        "stage_id": definition.stage_id,
        "stage_name": definition.stage_name,
        "signature": signature,
    }


def _validate_identity(
    document: dict[str, Any],
    expected_values: dict[str, str],
    document_label: str,
    stage_name: str,
) -> None:
    for field, expected_value in expected_values.items():
        if document[field] != expected_value:
            raise IntegrityError(
                f"{document_label} incohérent pour {stage_name}, champ {field}."
            )


def validate_attempt_outputs(
    attempt_dir: Path,
    definition: StageDefinition,
    run_id: str,
    signature: str,
) -> None:
    """Valide les documents et artefacts temporaires avant publication atomique."""
    stage_outputs = read_json(attempt_dir / "stage_outputs.json")
    validate_json_document(stage_outputs, "stage_outputs.schema.json")
    expected_values = _expected_identity(definition, run_id, signature)
    _validate_identity(stage_outputs, expected_values, "Sortie", definition.stage_name)

    published_stage_dir = PurePosixPath("stages") / definition.directory_name
    artifact_ids: set[str] = set()
    for artifact in stage_outputs["artifacts"]:
        if artifact["artifact_id"] in artifact_ids:
            raise IntegrityError(
                f"Identifiant d'artefact dupliqué : {artifact['artifact_id']}"
            )
        artifact_ids.add(artifact["artifact_id"])
        artifact_path = resolve_declared_artifact_path(
            artifact["path"], published_stage_dir, attempt_dir
        )
        if not artifact_path.is_file():
            raise IntegrityError(f"Artefact annoncé mais absent : {artifact['path']}")
        if sha256_file(artifact_path) != artifact["sha256"]:
            raise IntegrityError(f"Empreinte invalide : {artifact['path']}")

    audit = read_json(attempt_dir / "audit.json")
    validate_json_document(audit, "stage_audit.schema.json")
    _validate_identity(audit, expected_values, "Audit", definition.stage_name)


def validate_published_stage(
    run_dir: Path,
    definition: StageDefinition,
    stage_record: dict[str, Any],
) -> None:
    """Recalcule toutes les empreintes avant de réutiliser une étape publiée."""
    stage_dir = run_dir / "stages" / definition.directory_name
    stage_outputs_path = stage_dir / "stage_outputs.json"
    audit_path = stage_dir / "audit.json"
    # Le manifest protège aussi les documents de provenance. Sans ces deux
    # empreintes, un artefact intact pourrait être accompagné d'un audit altéré.
    if sha256_file(stage_outputs_path) != stage_record["stage_outputs_sha256"]:
        raise IntegrityError(f"Descripteur de sorties modifié : {stage_outputs_path}")
    if sha256_file(audit_path) != stage_record["audit_sha256"]:
        raise IntegrityError(f"Audit d'étape modifié : {audit_path}")

    run_id = load_manifest(run_dir)["run_id"]
    expected_values = _expected_identity(
        definition,
        run_id,
        stage_record["signature"],
    )
    stage_outputs = read_json(stage_outputs_path)
    validate_json_document(stage_outputs, "stage_outputs.schema.json")
    _validate_identity(
        stage_outputs,
        expected_values,
        "Sortie publiée",
        definition.stage_name,
    )

    published_stage_dir = PurePosixPath("stages") / definition.directory_name
    for artifact in stage_outputs["artifacts"]:
        artifact_path = resolve_declared_artifact_path(
            artifact["path"], published_stage_dir, stage_dir
        )
        if not artifact_path.is_file():
            raise IntegrityError(f"Artefact publié absent : {artifact['path']}")
        if sha256_file(artifact_path) != artifact["sha256"]:
            raise IntegrityError(f"Artefact publié modifié : {artifact['path']}")

    audit = read_json(audit_path)
    validate_json_document(audit, "stage_audit.schema.json")
    _validate_identity(
        audit,
        expected_values,
        "Audit publié",
        definition.stage_name,
    )
