"""Création atomique et auditée d'un run V2 vide."""

from __future__ import annotations

import hashlib
import json
import os
from datetime import datetime, timezone
from pathlib import Path
from time import monotonic
from typing import Any
from uuid import uuid4

import yaml

from effet_fondateur import __version__
from effet_fondateur.audit import atomic_write_json, sha256_file
from effet_fondateur.contracts import (
    build_file_artifact,
    load_pipeline_config,
    validate_json_document,
)
from effet_fondateur.orchestrator.environment import build_environment
from effet_fondateur.orchestrator.state import record_event, utc_now


def _canonical_sha256(document: dict[str, Any]) -> str:
    serialized = json.dumps(document, sort_keys=True, separators=(",", ":"))
    return hashlib.sha256(serialized.encode("utf-8")).hexdigest()


def _new_run_id(project_id: str) -> str:
    timestamp = datetime.now(timezone.utc).strftime("%Y-%m-%dT%H%M%SZ")
    return f"{timestamp}_{project_id}_{uuid4().hex[:8]}"


def _build_initialization_summary(config: dict[str, Any]) -> dict[str, Any]:
    """Résume uniquement la complétude technique, sans valeur moléculaire."""
    return {
        "schema_version": "1.0.0",
        "project_id": config["project"]["project_id"],
        "assembly": config["project"]["assembly"],
        "enabled_stage_count": sum(
            1 for stage in config["stages"].values() if stage["enabled"]
        ),
        "target_definition_complete": all(
            config["target"][field] is not None
            for field in ("chromosome", "position_bp", "ref", "alt", "project_variant_id")
        ),
    }


def _build_initialization_audit(
    run_id: str,
    signature: str,
    started_at: str,
    completed_at: str,
    duration_seconds: float,
    config_sha256: str,
    summary: dict[str, Any],
    summary_artifact: dict[str, Any],
) -> dict[str, Any]:
    """Construit l'audit non sensible de l'étape de bootstrap."""
    return {
        "schema_version": "1.0.0",
        "run_id": run_id,
        "stage_id": "00",
        "stage_name": "initialize_run",
        "method_id": "initialize_run_v1",
        "signature": signature,
        "started_at": started_at,
        "completed_at": completed_at,
        "duration_seconds": duration_seconds,
        "inputs": [{"artifact_id": "configuration", "sha256": config_sha256}],
        "outputs": [summary_artifact],
        "parameters": {},
        "tools": [],
        "counts": {"enabled_stages": summary["enabled_stage_count"]},
        "metrics": {},
        "exclusions": [],
        "warnings": [],
        "checks": [{"check": "configuration_schema", "status": "PASS"}],
        "known_limits": ["Les outils externes ne sont requis qu'au démarrage de leur étape."],
        "expected_visualizations": [],
        "manual_validation_required": not summary["target_definition_complete"],
    }


def _build_initial_manifest(
    config: dict[str, Any],
    run_id: str,
    signature: str,
    started_at: str,
    completed_at: str,
    duration_seconds: float,
    config_sha256: str,
    stage_relative_dir: str,
    stage_dir: Path,
    target_definition_complete: bool,
) -> dict[str, Any]:
    """Construit l'état initial sans créer de dépendance circulaire d'empreinte."""
    return {
        "schema_version": "1.0.0",
        "run_id": run_id,
        "project_id": config["project"]["project_id"],
        "created_at": started_at,
        "updated_at": completed_at,
        "global_status": "INCOMPLETE",
        "config_path": "config.resolved.yaml",
        "config_sha256": config_sha256,
        "environment_path": "environment.json",
        "stages": [
            {
                "stage_id": "00",
                "stage_name": "initialize_run",
                "state": "SUCCEEDED",
                "critical": True,
                "signature": signature,
                "started_at": started_at,
                "completed_at": completed_at,
                "duration_seconds": duration_seconds,
                "audit_path": f"{stage_relative_dir}/audit.json",
                "audit_sha256": sha256_file(stage_dir / "audit.json"),
                "stage_outputs_sha256": sha256_file(stage_dir / "stage_outputs.json"),
                "attempt_count": 1,
                "last_error_code": None,
            }
        ],
        "manual_decisions_required": (
            [] if target_definition_complete else ["target_variant_definition"]
        ),
    }


def _write_bootstrap_documents(
    config: dict[str, Any],
    run_id: str,
    temporary_run_dir: Path,
    stage_dir: Path,
    stage_relative_dir: str,
) -> tuple[str, str, dict[str, Any], dict[str, Any]]:
    """Écrit configuration, environnement et contrats de sortie de l'étape `00`."""
    resolved_config_path = temporary_run_dir / "config.resolved.yaml"
    resolved_config_path.write_text(
        yaml.safe_dump(config, allow_unicode=True, sort_keys=False),
        encoding="utf-8",
    )
    config_sha256 = sha256_file(resolved_config_path)
    signature = _canonical_sha256(
        {
            "method_id": "initialize_run_v1",
            "config_sha256": config_sha256,
            "pipeline_version": __version__,
        }
    )
    stage_inputs = {
        "schema_version": "1.0.0",
        "run_id": run_id,
        "stage_id": "00",
        "stage_name": "initialize_run",
        "signature": signature,
        "attempt_number": 1,
        "published_output_dir": stage_relative_dir,
        "parameters": {},
        "artifacts": [],
    }
    validate_json_document(stage_inputs, "stage_inputs.schema.json")
    atomic_write_json(stage_dir / "stage_inputs.json", stage_inputs)
    atomic_write_json(temporary_run_dir / "environment.json", build_environment(config))

    summary = _build_initialization_summary(config)
    atomic_write_json(stage_dir / "initialization_summary.json", summary)
    summary_artifact = build_file_artifact(
        physical_path=stage_dir / "initialization_summary.json",
        published_path=f"{stage_relative_dir}/initialization_summary.json",
        artifact_id="initialization_summary",
        artifact_type="initialization_summary",
        media_type="application/json",
        producer_stage="00_initialize_run",
        producer_signature=signature,
    )
    stage_outputs = {
        "schema_version": "1.0.0",
        "run_id": run_id,
        "stage_id": "00",
        "stage_name": "initialize_run",
        "signature": signature,
        "artifacts": [summary_artifact],
    }
    validate_json_document(stage_outputs, "stage_outputs.schema.json")
    atomic_write_json(stage_dir / "stage_outputs.json", stage_outputs)
    return config_sha256, signature, summary, summary_artifact


def initialize_run(config_path: Path, runs_dir: Path) -> Path:
    """Crée le run `00` ou lève une erreur avant toute publication partielle.

    La configuration est copiée et empreintée, l'environnement est inventorié,
    puis le dossier temporaire complet est renommé vers son chemin définitif.
    """
    started_at = utc_now()
    started_clock = monotonic()
    config = load_pipeline_config(config_path)
    run_id = _new_run_id(config["project"]["project_id"])
    final_run_dir = runs_dir / run_id
    temporary_run_dir = runs_dir / f".{run_id}.{uuid4().hex}.tmp"
    stage_relative_dir = "stages/00_initialize_run"
    stage_dir = temporary_run_dir / stage_relative_dir

    if final_run_dir.exists():
        raise FileExistsError(f"Le run existe déjà : {final_run_dir}")

    # Le run entier est construit hors de son chemin final : aucune autre
    # commande ne peut observer un bootstrap partiellement écrit.
    stage_dir.mkdir(parents=True)
    config_sha256, signature, summary, summary_artifact = _write_bootstrap_documents(
        config,
        run_id,
        temporary_run_dir,
        stage_dir,
        stage_relative_dir,
    )

    completed_at = utc_now()
    duration_seconds = monotonic() - started_clock
    audit = _build_initialization_audit(
        run_id=run_id,
        signature=signature,
        started_at=started_at,
        completed_at=completed_at,
        duration_seconds=duration_seconds,
        config_sha256=config_sha256,
        summary=summary,
        summary_artifact=summary_artifact,
    )
    validate_json_document(audit, "stage_audit.schema.json")
    atomic_write_json(stage_dir / "audit.json", audit)
    (stage_dir / "checksums.sha256").write_text(
        f"{summary_artifact['sha256']}  initialization_summary.json\n",
        encoding="utf-8",
    )

    # Le manifest empreinte audit et sorties, mais jamais lui-même. Cela évite
    # une dépendance circulaire tout en protégeant la provenance publiée.
    manifest = _build_initial_manifest(
        config=config,
        run_id=run_id,
        signature=signature,
        started_at=started_at,
        completed_at=completed_at,
        duration_seconds=duration_seconds,
        config_sha256=config_sha256,
        stage_relative_dir=stage_relative_dir,
        stage_dir=stage_dir,
        target_definition_complete=summary["target_definition_complete"],
    )
    validate_json_document(manifest, "run_manifest.schema.json")
    atomic_write_json(temporary_run_dir / "manifest.json", manifest)
    (temporary_run_dir / "events.jsonl").touch()
    record_event(
        temporary_run_dir,
        "00_initialize_run",
        "stage_succeeded",
        details={"manual_validation_required": audit["manual_validation_required"]},
    )

    runs_dir.mkdir(parents=True, exist_ok=True)
    os.replace(temporary_run_dir, final_run_dir)
    return final_run_dir
