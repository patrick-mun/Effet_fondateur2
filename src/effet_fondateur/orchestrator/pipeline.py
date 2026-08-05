"""Boucle minimale du pipeline V2 et reprise d'un run existant."""

from __future__ import annotations

from pathlib import Path
from typing import Iterable

from effet_fondateur.audit import sha256_file
from effet_fondateur.contracts import load_pipeline_config
from effet_fondateur.orchestrator.catalog import build_stage_catalog
from effet_fondateur.orchestrator.errors import IntegrityError, PipelineError
from effet_fondateur.orchestrator.models import StageDefinition
from effet_fondateur.orchestrator.runner import run_stage_with_failure_audit
from effet_fondateur.orchestrator.state import load_manifest, save_manifest
from effet_fondateur.stages.initialize_run import initialize_run


SYNTHETIC_STAGE = StageDefinition(
    stage_id="T00",
    stage_name="synthetic_stage",
    module="effet_fondateur.stages.synthetic_stage",
    critical=True,
    dependencies=("initialize_run",),
)

VALIDATE_SOURCES_STAGE = StageDefinition(
    stage_id="01",
    stage_name="validate_sources",
    module="effet_fondateur.stages.validate_sources",
    critical=True,
    dependencies=("initialize_run",),
    config_input_directories=("acpa_samples_dir",),
)

BUILD_SAMPLE_REGISTRY_STAGE = StageDefinition(
    stage_id="02",
    stage_name="build_sample_registry",
    module="effet_fondateur.stages.build_sample_registry",
    critical=True,
    dependencies=("validate_sources",),
    config_input_files=("sample_metadata",),
    required_artifact_ids=("source_inventory", "source_qc"),
    manual_decision_id="sample_registry_approval",
)

DEFAULT_STAGE_DEFINITIONS = (
    SYNTHETIC_STAGE,
    VALIDATE_SOURCES_STAGE,
    BUILD_SAMPLE_REGISTRY_STAGE,
)


def _run_enabled_stages(
    run_dir: Path,
    definitions: Iterable[StageDefinition],
) -> None:
    """Exécute les étapes activées après contrôle du catalogue et du run."""
    resolved_config_path = run_dir / "config.resolved.yaml"
    manifest = load_manifest(run_dir)
    # Une reprise doit être strictement identique. Tout changement de paramètre
    # scientifique exige un nouveau run au lieu de réécrire son histoire.
    if sha256_file(resolved_config_path) != manifest["config_sha256"]:
        manifest["global_status"] = "BLOCKED"
        save_manifest(run_dir, manifest)
        raise IntegrityError(
            f"La configuration résolue a été modifiée dans le run : {resolved_config_path}"
        )
    config = load_pipeline_config(resolved_config_path)
    definitions_by_name = build_stage_catalog(definitions)
    for stage_name, stage_config in config["stages"].items():
        if stage_name == "initialize_run" or not stage_config["enabled"]:
            continue
        definition = definitions_by_name.get(stage_name)
        if definition is None:
            raise PipelineError(f"Étape activée mais non implémentée : {stage_name}")
        run_stage_with_failure_audit(
            run_dir,
            definition,
            stage_config["parameters"],
        )

    manifest = load_manifest(run_dir)
    terminal_states = {"SUCCEEDED", "CACHED", "SKIPPED"}
    if (
        all(stage["state"] in terminal_states for stage in manifest["stages"])
        and not manifest["manual_decisions_required"]
    ):
        manifest["global_status"] = "TECHNICALLY_VALID"
        save_manifest(run_dir, manifest)


def run_pipeline(
    config_path: Path,
    runs_dir: Path,
    definitions: Iterable[StageDefinition] = DEFAULT_STAGE_DEFINITIONS,
) -> Path:
    """Crée un run puis exécute les étapes actuellement implémentées et activées."""
    run_dir = initialize_run(config_path, runs_dir)
    _run_enabled_stages(run_dir, definitions)
    return run_dir


def resume_pipeline(
    run_dir: Path,
    definitions: Iterable[StageDefinition] = DEFAULT_STAGE_DEFINITIONS,
) -> Path:
    """Reprend un run en validant les sorties déjà publiées avant réutilisation."""
    load_manifest(run_dir)
    _run_enabled_stages(run_dir, definitions)
    return run_dir
