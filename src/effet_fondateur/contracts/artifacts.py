"""Construction centralisée des descripteurs d'artefacts V2."""

from __future__ import annotations

from pathlib import Path
from typing import Literal

from effet_fondateur.audit import sha256_file


Sensitivity = Literal["public", "internal", "sensitive_genetic"]


def build_file_artifact(
    physical_path: Path,
    published_path: str,
    artifact_id: str,
    artifact_type: str,
    media_type: str,
    producer_stage: str,
    producer_signature: str,
    *,
    schema_name: str | None = None,
    schema_version: str | None = None,
    assembly: str | None = None,
    sample_set_id: str | None = None,
    variant_set_id: str | None = None,
    sensitivity: Sensitivity = "internal",
) -> dict[str, str | None]:
    """Décrit un fichier existant et calcule son empreinte au dernier moment."""
    return {
        "artifact_id": artifact_id,
        "artifact_type": artifact_type,
        "path": published_path,
        "media_type": media_type,
        "schema_name": schema_name,
        "schema_version": schema_version,
        "sha256": sha256_file(physical_path),
        "producer_stage": producer_stage,
        "producer_signature": producer_signature,
        "assembly": assembly,
        "sample_set_id": sample_set_id,
        "variant_set_id": variant_set_id,
        "sensitivity": sensitivity,
    }
