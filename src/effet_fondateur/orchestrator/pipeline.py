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

CONVERT_ACPA_STAGE = StageDefinition(
    stage_id="03",
    stage_name="convert_acpa",
    module="effet_fondateur.stages.convert_acpa",
    critical=True,
    dependencies=("build_sample_registry",),
    config_input_directories=("acpa_samples_dir",),
    required_artifact_ids=("samples_master",),
    blocking_manual_decision_ids=("sample_registry_approval",),
)

PREPARE_TARGET_VARIANT_DATASET_STAGE = StageDefinition(
    stage_id="04",
    stage_name="prepare_target_variant_dataset",
    module="effet_fondateur.stages.prepare_target_variant_dataset",
    critical=True,
    dependencies=("build_sample_registry", "convert_acpa"),
    config_input_files=("target_variant_metadata", "target_variant_genotypes"),
    required_artifact_ids=(
        "samples_master",
        "target_chromosome_base_bed",
        "target_chromosome_base_bim",
        "target_chromosome_base_fam",
        "target_chromosome_base_dataset",
    ),
    blocking_manual_decision_ids=("sample_registry_approval",),
)

QC_PRELIMINARY_STAGE = StageDefinition(
    stage_id="05",
    stage_name="qc_preliminary",
    module="effet_fondateur.stages.qc_preliminary",
    critical=True,
    dependencies=("build_sample_registry", "convert_acpa"),
    required_artifact_ids=(
        "samples_master",
        "genomewide_base_bed",
        "genomewide_base_bim",
        "genomewide_base_fam",
        "genomewide_base_dataset",
    ),
    blocking_manual_decision_ids=("sample_registry_approval",),
)

BUILD_KINSHIP_PANEL_STAGE = StageDefinition(
    stage_id="06",
    stage_name="build_kinship_panel",
    module="effet_fondateur.stages.build_kinship_panel",
    critical=True,
    dependencies=("qc_preliminary",),
    required_artifact_ids=(
        "genomewide_pre_qc_bed",
        "genomewide_pre_qc_bim",
        "genomewide_pre_qc_fam",
        "genomewide_pre_qc_dataset",
        "qc_variant_metrics",
    ),
)

INFER_KINSHIP_STAGE = StageDefinition(
    stage_id="07",
    stage_name="infer_kinship",
    module="effet_fondateur.stages.infer_kinship",
    critical=True,
    dependencies=("build_sample_registry", "qc_preliminary", "build_kinship_panel"),
    required_artifact_ids=(
        "samples_master",
        "qc_individual_metrics",
        "kinship_panel_bed",
        "kinship_panel_bim",
        "kinship_panel_fam",
        "kinship_panel_dataset",
    ),
    manual_decision_id="kinship_exclusion_approval",
)

ANALYZE_POPULATION_STRUCTURE_STAGE = StageDefinition(
    stage_id="08",
    stage_name="analyze_population_structure",
    module="effet_fondateur.stages.analyze_population_structure",
    critical=True,
    dependencies=("build_sample_registry", "build_kinship_panel", "infer_kinship"),
    required_artifact_ids=(
        "samples_master",
        "kinship_panel_bed",
        "kinship_panel_bim",
        "kinship_panel_fam",
        "kinship_panel_dataset",
        "independent_set_proposal",
    ),
    manual_decision_id="population_structure_exclusion_approval",
)

FREEZE_COHORTS_STAGE = StageDefinition(
    stage_id="09",
    stage_name="freeze_cohorts",
    module="effet_fondateur.stages.freeze_cohorts",
    critical=True,
    dependencies=(
        "build_sample_registry",
        "prepare_target_variant_dataset",
        "qc_preliminary",
        "infer_kinship",
        "analyze_population_structure",
    ),
    config_input_files=("target_variant_metadata",),
    required_artifact_ids=(
        "samples_master",
        "target_genotype_audit",
        "qc_individual_metrics",
        "independent_set_proposal",
        "population_scores",
        "population_outliers",
    ),
    resolves_manual_decision_ids=(
        "kinship_exclusion_approval",
        "population_structure_exclusion_approval",
    ),
)

QC_FINAL_STAGE = StageDefinition(
    stage_id="10",
    stage_name="qc_final",
    module="effet_fondateur.stages.qc_final",
    critical=True,
    dependencies=(
        "build_sample_registry",
        "prepare_target_variant_dataset",
        "qc_preliminary",
        "freeze_cohorts",
    ),
    config_input_files=("target_variant_metadata",),
    required_artifact_ids=(
        "samples_master",
        "genomewide_pre_qc_bed",
        "genomewide_pre_qc_bim",
        "genomewide_pre_qc_fam",
        "genomewide_pre_qc_dataset",
        "target_variant_bed",
        "target_variant_bim",
        "target_variant_fam",
        "target_variant_dataset",
        "cohorts_frozen",
        "cohort_keep_controls_unrelated",
        "cohort_keep_target_carriers_independent",
        "cohort_keep_target_chromosome_all_qc",
    ),
)

PREPARE_TARGET_REGION_STAGE = StageDefinition(
    stage_id="11",
    stage_name="prepare_target_region",
    module="effet_fondateur.stages.prepare_target_region",
    critical=True,
    dependencies=("prepare_target_variant_dataset", "qc_final"),
    config_input_files=("target_variant_metadata", "genetic_map"),
    required_artifact_ids=(
        "target_chromosome_all_qc_bed",
        "target_chromosome_all_qc_bim",
        "target_chromosome_all_qc_fam",
        "target_chromosome_all_qc_dataset",
        "variant_qc_final",
    ),
)

PHASE_TARGET_REGION_STAGE = StageDefinition(
    stage_id="12",
    stage_name="phase_target_region",
    module="effet_fondateur.stages.phase_target_region",
    critical=True,
    dependencies=(
        "build_sample_registry",
        "prepare_target_variant_dataset",
        "prepare_target_region",
    ),
    config_input_files=("target_variant_metadata", "reference_panel_catalog"),
    required_artifact_ids=(
        "target_region_bim",
        "target_region_bed",
        "target_region_fam",
        "target_genetic_map",
        "phasing_input_manifest",
        "samples_master",
        "target_genotype_audit",
    ),
)

INFER_FOUNDER_HAPLOTYPE_STAGE = StageDefinition(
    stage_id="13",
    stage_name="infer_founder_haplotype",
    module="effet_fondateur.stages.infer_founder_haplotype",
    critical=True,
    dependencies=(
        "build_sample_registry",
        "freeze_cohorts",
        "prepare_target_region",
        "phase_target_region",
    ),
    config_input_files=(),
    required_artifact_ids=(
        "shapeit5_final_bcf",
        "shapeit5_final_index",
        "carrier_haplotypes",
        "target_genetic_map",
        "cohorts_frozen",
        "samples_master",
    ),
)

ESTIMATE_VARIANT_AGE_STAGE = StageDefinition(
    stage_id="14",
    stage_name="estimate_variant_age",
    module="effet_fondateur.stages.estimate_variant_age",
    critical=True,
    dependencies=("infer_founder_haplotype",),
    config_input_files=(),
    required_artifact_ids=(
        "founder_segments",
        "founder_consensus",
        "founder_analysis_summary",
    ),
)

ANALYZE_ROH_STAGE = StageDefinition(
    stage_id="16",
    stage_name="analyze_roh",
    module="effet_fondateur.stages.analyze_roh",
    critical=False,
    dependencies=(
        "build_sample_registry",
        "prepare_target_variant_dataset",
        "freeze_cohorts",
        "qc_final",
        "prepare_target_region",
    ),
    config_input_files=("target_variant_metadata",),
    required_artifact_ids=(
        "samples_master",
        "target_genotype_audit",
        "cohorts_frozen",
        "controls_unrelated_qc_bed",
        "controls_unrelated_qc_bim",
        "controls_unrelated_qc_fam",
        "controls_unrelated_qc_dataset",
        "target_carriers_independent_qc_bed",
        "target_carriers_independent_qc_bim",
        "target_carriers_independent_qc_fam",
        "target_carriers_independent_qc_dataset",
        "target_chromosome_all_qc_bed",
        "target_chromosome_all_qc_bim",
        "target_chromosome_all_qc_fam",
        "target_chromosome_all_qc_dataset",
        "target_genetic_map",
    ),
)

DEFAULT_STAGE_DEFINITIONS = (
    SYNTHETIC_STAGE,
    VALIDATE_SOURCES_STAGE,
    BUILD_SAMPLE_REGISTRY_STAGE,
    CONVERT_ACPA_STAGE,
    PREPARE_TARGET_VARIANT_DATASET_STAGE,
    QC_PRELIMINARY_STAGE,
    BUILD_KINSHIP_PANEL_STAGE,
    INFER_KINSHIP_STAGE,
    ANALYZE_POPULATION_STRUCTURE_STAGE,
    FREEZE_COHORTS_STAGE,
    QC_FINAL_STAGE,
    PREPARE_TARGET_REGION_STAGE,
    PHASE_TARGET_REGION_STAGE,
    INFER_FOUNDER_HAPLOTYPE_STAGE,
    ESTIMATE_VARIANT_AGE_STAGE,
    ANALYZE_ROH_STAGE,
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
