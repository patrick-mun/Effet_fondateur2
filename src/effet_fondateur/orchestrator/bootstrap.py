"""Création atomique et auditée d'un run V2 vide."""

from __future__ import annotations

import hashlib
import json
import os
import platform
import shutil
import subprocess
import sys
from datetime import datetime, timezone
from pathlib import Path
from time import monotonic
from typing import Any
from uuid import uuid4

import yaml

from effet_fondateur import __version__
from effet_fondateur.audit import atomic_write_json, sha256_file
from effet_fondateur.contracts import load_pipeline_config, validate_json_document
from effet_fondateur.orchestrator.state import record_event, utc_now


TOOL_VERSION_ARGUMENTS = {
    "plink": ("--version",),
    "king": ("--version",),
    "bcftools": ("--version",),
    "rscript": ("--version",),
}


def _canonical_sha256(document: dict[str, Any]) -> str:
    serialized = json.dumps(document, sort_keys=True, separators=(",", ":"))
    return hashlib.sha256(serialized.encode("utf-8")).hexdigest()


def _capture_tool(command_name: str | None, tool_key: str) -> dict[str, Any]:
    if command_name is None:
        return {"configured": None, "available": False, "version": None}
    resolved_path = shutil.which(command_name)
    if resolved_path is None:
        return {"configured": command_name, "available": False, "version": None}

    version = None
    arguments = TOOL_VERSION_ARGUMENTS.get(tool_key, ("--version",))
    try:
        completed_process = subprocess.run(
            [resolved_path, *arguments],
            capture_output=True,
            check=False,
            text=True,
            timeout=5,
        )
        version_output = completed_process.stdout.strip() or completed_process.stderr.strip()
        if version_output:
            version = version_output.splitlines()[0][:500]
    except (OSError, subprocess.TimeoutExpired):
        version = None
    return {"configured": command_name, "available": True, "version": version}


def _build_environment(config: dict[str, Any]) -> dict[str, Any]:
    return {
        "schema_version": "1.0.0",
        "python": {
            "version": platform.python_version(),
            "implementation": platform.python_implementation(),
            "executable": Path(sys.executable).name,
        },
        "platform": {
            "system": platform.system(),
            "release": platform.release(),
            "machine": platform.machine(),
        },
        "pipeline_version": __version__,
        "tools": {
            tool_key: _capture_tool(command_name, tool_key)
            for tool_key, command_name in config["tools"].items()
        },
    }


def _new_run_id(project_id: str) -> str:
    timestamp = datetime.now(timezone.utc).strftime("%Y-%m-%dT%H%M%SZ")
    return f"{timestamp}_{project_id}_{uuid4().hex[:8]}"


def _artifact(
    run_dir: Path,
    relative_path: str,
    artifact_id: str,
    artifact_type: str,
    signature: str,
) -> dict[str, Any]:
    return {
        "artifact_id": artifact_id,
        "artifact_type": artifact_type,
        "path": relative_path,
        "media_type": "application/json",
        "schema_name": None,
        "schema_version": None,
        "sha256": sha256_file(run_dir / relative_path),
        "producer_stage": "00_initialize_run",
        "producer_signature": signature,
        "assembly": None,
        "sample_set_id": None,
        "variant_set_id": None,
        "sensitivity": "internal",
    }


def initialize_run(config_path: Path, runs_dir: Path) -> Path:
    """Valide la configuration et publie atomiquement un nouveau run audité."""
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

    stage_dir.mkdir(parents=True)
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

    environment = _build_environment(config)
    atomic_write_json(temporary_run_dir / "environment.json", environment)
    summary = {
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
    atomic_write_json(stage_dir / "initialization_summary.json", summary)
    summary_artifact = _artifact(
        temporary_run_dir,
        f"{stage_relative_dir}/initialization_summary.json",
        "initialization_summary",
        "initialization_summary",
        signature,
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

    completed_at = utc_now()
    duration_seconds = monotonic() - started_clock
    audit = {
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
    validate_json_document(audit, "stage_audit.schema.json")
    atomic_write_json(stage_dir / "audit.json", audit)
    (stage_dir / "checksums.sha256").write_text(
        f"{summary_artifact['sha256']}  initialization_summary.json\n",
        encoding="utf-8",
    )

    manifest = {
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
            ["target_variant_definition"]
            if not summary["target_definition_complete"]
            else []
        ),
    }
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
