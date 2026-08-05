"""Résolution contrôlée des sources externes déclarées pour une étape."""

from __future__ import annotations

from pathlib import Path
from typing import Any

from effet_fondateur.audit import read_json, sha256_file
from effet_fondateur.orchestrator.errors import StageExecutionError
from effet_fondateur.orchestrator.models import StageDefinition


def _resolve_configured_path(configured_path: str) -> Path:
    path = Path(configured_path)
    return path if path.is_absolute() else Path.cwd() / path


def _source_artifact(
    source_path: Path,
    configured_root: str,
    root_path: Path,
    index: int,
    config_sha256: str,
    assembly: str,
) -> dict[str, str | None]:
    relative_path = source_path.relative_to(root_path)
    configured_source_path = Path(configured_root) / relative_path
    return {
        "artifact_id": f"acpa_source_{index:06d}",
        "artifact_type": "acpa_source_file",
        "path": configured_source_path.as_posix(),
        "media_type": "text/tab-separated-values",
        "schema_name": None,
        "schema_version": None,
        "sha256": sha256_file(source_path),
        "producer_stage": "external_source",
        "producer_signature": config_sha256,
        "assembly": assembly,
        "sample_set_id": None,
        "variant_set_id": None,
        "sensitivity": "sensitive_genetic",
    }


def _config_file_artifact(
    configured_path: str,
    input_key: str,
    config_sha256: str,
    assembly: str,
) -> dict[str, str | None]:
    physical_path = _resolve_configured_path(configured_path)
    if physical_path.is_symlink() or not physical_path.is_file():
        raise StageExecutionError(
            f"Fichier source invalide pour inputs.{input_key}.",
            2,
        )
    return {
        "artifact_id": f"config_input_{input_key}",
        "artifact_type": input_key,
        "path": Path(configured_path).as_posix(),
        "media_type": "text/tab-separated-values",
        "schema_name": None,
        "schema_version": None,
        "sha256": sha256_file(physical_path),
        "producer_stage": "external_source",
        "producer_signature": config_sha256,
        "assembly": assembly,
        "sample_set_id": None,
        "variant_set_id": None,
        "sensitivity": "sensitive_genetic",
    }


def _dependency_artifacts(
    run_dir: Path,
    definition: StageDefinition,
) -> list[dict[str, Any]]:
    required_ids = set(definition.required_artifact_ids)
    if not required_ids:
        return []
    matching_artifacts: dict[str, dict[str, Any]] = {}
    for dependency in definition.dependencies:
        stage_directories = list((run_dir / "stages").glob(f"*_{dependency}"))
        if len(stage_directories) != 1:
            continue
        stage_outputs = read_json(stage_directories[0] / "stage_outputs.json")
        for artifact in stage_outputs["artifacts"]:
            artifact_id = artifact["artifact_id"]
            if artifact_id in required_ids:
                if artifact_id in matching_artifacts:
                    raise StageExecutionError(
                        f"Artefact dépendant ambigu : {artifact_id}",
                        2,
                    )
                physical_path = run_dir / artifact["path"]
                if (
                    not physical_path.is_file()
                    or sha256_file(physical_path) != artifact["sha256"]
                ):
                    raise StageExecutionError(
                        f"Artefact dépendant absent ou modifié : {artifact_id}",
                        5,
                    )
                matching_artifacts[artifact_id] = artifact
    missing_ids = required_ids - matching_artifacts.keys()
    if missing_ids:
        raise StageExecutionError(
            f"Artefacts dépendants absents : {', '.join(sorted(missing_ids))}",
            2,
        )
    return [
        matching_artifacts[artifact_id]
        for artifact_id in definition.required_artifact_ids
    ]


def resolve_stage_input_artifacts(
    config: dict[str, Any],
    definition: StageDefinition,
    config_sha256: str,
    run_dir: Path,
) -> list[dict[str, str | None]]:
    """Inventorie et empreinte les fichiers des répertoires requis par l'étape."""
    artifacts: list[dict[str, str | None]] = []
    artifact_index = 1
    for input_key in definition.config_input_directories:
        configured_root = config["inputs"].get(input_key)
        if configured_root is None:
            raise StageExecutionError(
                f"Entrée de configuration obligatoire absente : inputs.{input_key}",
                2,
            )
        root_path = _resolve_configured_path(configured_root)
        if root_path.is_symlink() or not root_path.is_dir():
            raise StageExecutionError(
                f"Répertoire source invalide pour inputs.{input_key}.",
                2,
            )
        source_entries = sorted(root_path.rglob("*"))
        if any(path.is_symlink() for path in source_entries):
            raise StageExecutionError(
                f"Lien symbolique interdit dans inputs.{input_key}.",
                2,
            )
        source_paths = [path for path in source_entries if path.is_file()]
        if not source_paths:
            raise StageExecutionError(
                f"Répertoire source vide pour inputs.{input_key}.",
                2,
            )
        for source_path in source_paths:
            artifacts.append(
                _source_artifact(
                    source_path=source_path,
                    configured_root=configured_root,
                    root_path=root_path,
                    index=artifact_index,
                    config_sha256=config_sha256,
                    assembly=config["project"]["assembly"],
                )
            )
            artifact_index += 1
    for input_key in definition.config_input_files:
        configured_path = config["inputs"].get(input_key)
        if configured_path is None:
            raise StageExecutionError(
                f"Entrée de configuration obligatoire absente : inputs.{input_key}",
                2,
            )
        artifacts.append(
            _config_file_artifact(
                configured_path=configured_path,
                input_key=input_key,
                config_sha256=config_sha256,
                assembly=config["project"]["assembly"],
            )
        )
    artifacts.extend(_dependency_artifacts(run_dir, definition))
    return artifacts
