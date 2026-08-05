"""Résolution contrôlée des sources externes déclarées pour une étape."""

from __future__ import annotations

from pathlib import Path
from typing import Any

from effet_fondateur.audit import sha256_file
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


def resolve_stage_input_artifacts(
    config: dict[str, Any],
    definition: StageDefinition,
    config_sha256: str,
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
    return artifacts
