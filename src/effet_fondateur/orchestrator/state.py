"""Gestion atomique du manifest et du journal d'un run."""

from __future__ import annotations

from datetime import datetime, timezone
from pathlib import Path
from typing import Any

from effet_fondateur.audit import append_json_event, atomic_write_json, read_json
from effet_fondateur.contracts import validate_json_document


def utc_now() -> str:
    """Retourne une date UTC ISO-8601 compatible avec les schémas d'audit."""
    return datetime.now(timezone.utc).isoformat().replace("+00:00", "Z")


def load_manifest(run_dir: Path) -> dict[str, Any]:
    """Charge et valide le manifest courant d'un run."""
    manifest = read_json(run_dir / "manifest.json")
    validate_json_document(manifest, "run_manifest.schema.json")
    return manifest


def save_manifest(run_dir: Path, manifest: dict[str, Any]) -> None:
    """Valide puis remplace atomiquement le manifest courant."""
    manifest["updated_at"] = utc_now()
    validate_json_document(manifest, "run_manifest.schema.json")
    atomic_write_json(run_dir / "manifest.json", manifest)


def record_event(
    run_dir: Path,
    stage: str,
    event: str,
    severity: str = "INFO",
    details: dict[str, Any] | None = None,
) -> None:
    """Ajoute un événement agrégé sans génotype ni identifiant individuel."""
    append_json_event(
        run_dir / "events.jsonl",
        {
            "timestamp": utc_now(),
            "stage": stage,
            "event": event,
            "severity": severity,
            "details": details or {},
        },
    )


def find_stage_record(
    manifest: dict[str, Any], stage_name: str
) -> dict[str, Any] | None:
    """Retourne l'enregistrement d'une étape sans en créer implicitement."""
    for stage_record in manifest["stages"]:
        if stage_record["stage_name"] == stage_name:
            return stage_record
    return None
