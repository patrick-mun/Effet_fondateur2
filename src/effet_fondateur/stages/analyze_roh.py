"""Étape 16 : analyse secondaire des runs of homozygosity."""

from __future__ import annotations

import argparse
import json
import math
import shutil
import subprocess
import sys
from collections import defaultdict
from pathlib import Path, PurePosixPath
from time import monotonic
from typing import Any, Sequence

import yaml

from effet_fondateur.audit import atomic_write_json, read_json, sha256_file
from effet_fondateur.contracts import (
    DocumentValidationError, TableValidationError, build_file_artifact,
    load_pipeline_config, validate_cohorts_frozen, validate_json_document,
    validate_tsv_table,
)
from effet_fondateur.orchestrator.state import utc_now
from effet_fondateur.roh import RohAnalysisError, publish_roh


GENOMEWIDE_DATASETS = {
    "controls_unrelated": "controls_unrelated_qc",
    "target_carriers_independent": "target_carriers_independent_qc",
}
TARGET_DATASET_ID = "target_chromosome_all_qc"


class AnalyzeRohInputError(ValueError):
    """Signale une configuration ou une entrée invalide de l'étape 16."""


def _build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(prog="analyze-roh")
    parser.add_argument("--stage-inputs", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    return parser


def _positive_integer(parameters: dict[str, Any], name: str, default: int) -> int:
    value = parameters.get(name, default)
    if isinstance(value, bool) or not isinstance(value, int) or value <= 0:
        raise AnalyzeRohInputError(f"invalid_parameter:{name}")
    return value


def _nonnegative_integer(parameters: dict[str, Any], name: str, default: int) -> int:
    value = parameters.get(name, default)
    if isinstance(value, bool) or not isinstance(value, int) or value < 0:
        raise AnalyzeRohInputError(f"invalid_parameter:{name}")
    return value


def _positive_number(parameters: dict[str, Any], name: str, default: float) -> float:
    value = parameters.get(name, default)
    if isinstance(value, bool) or not isinstance(value, (int, float)):
        raise AnalyzeRohInputError(f"invalid_parameter:{name}")
    result = float(value)
    if not math.isfinite(result) or result <= 0:
        raise AnalyzeRohInputError(f"invalid_parameter:{name}")
    return result


def _parameters(parameters: dict[str, Any]) -> dict[str, Any]:
    if parameters.get("method", "plink19_array_roh_secondary_v1") != "plink19_array_roh_secondary_v1":
        raise AnalyzeRohInputError("invalid_parameter:method")
    window_threshold = parameters.get("window_threshold", 0.05)
    if isinstance(window_threshold, bool) or not isinstance(window_threshold, (int, float)) or not 0 < float(window_threshold) <= 1:
        raise AnalyzeRohInputError("invalid_parameter:window_threshold")
    minimum_group = _positive_integer(parameters, "minimum_group_samples", 5)
    minimum_primary = _positive_integer(parameters, "minimum_primary_samples", 20)
    if minimum_primary < minimum_group:
        raise AnalyzeRohInputError("invalid_parameter:minimum_primary_samples")
    denominator = parameters.get("autosomal_denominator_kb")
    denominator_source = parameters.get("autosomal_denominator_source")
    if (denominator is None) != (denominator_source is None):
        raise AnalyzeRohInputError("autosomal_denominator_requires_value_and_source")
    if denominator is not None:
        denominator = _positive_number(parameters, "autosomal_denominator_kb", 1.0)
        if not isinstance(denominator_source, str) or not denominator_source.strip():
            raise AnalyzeRohInputError("invalid_parameter:autosomal_denominator_source")
    return {
        "method": "plink19_array_roh_secondary_v1",
        "minimum_roh_kb": _positive_number(parameters, "minimum_roh_kb", 1500.0),
        "minimum_roh_snps": _positive_integer(parameters, "minimum_roh_snps", 50),
        "maximum_density_kb_per_snp": _positive_number(parameters, "maximum_density_kb_per_snp", 50.0),
        "maximum_gap_kb": _positive_number(parameters, "maximum_gap_kb", 1000.0),
        "window_snps": _positive_integer(parameters, "window_snps", 50),
        "window_max_heterozygotes": _nonnegative_integer(parameters, "window_max_heterozygotes", 1),
        "window_max_missing": _nonnegative_integer(parameters, "window_max_missing", 5),
        "window_threshold": float(window_threshold),
        "minimum_genomewide_variants": _positive_integer(parameters, "minimum_genomewide_variants", 10_000),
        "minimum_genomewide_autosomes": _positive_integer(parameters, "minimum_genomewide_autosomes", 22),
        "minimum_target_chromosome_variants": _positive_integer(parameters, "minimum_target_chromosome_variants", 100),
        "minimum_group_samples": minimum_group,
        "minimum_primary_samples": minimum_primary,
        "autosomal_denominator_kb": denominator,
        "autosomal_denominator_source": denominator_source,
        "plink_timeout_seconds": _positive_integer(parameters, "plink_timeout_seconds", 300),
    }


def _artifact_by_id(stage_inputs: dict[str, Any], artifact_id: str) -> dict[str, Any]:
    matches = [artifact for artifact in stage_inputs["artifacts"] if artifact["artifact_id"] == artifact_id]
    if len(matches) != 1:
        raise AnalyzeRohInputError(f"{artifact_id}_missing_or_ambiguous")
    return matches[0]


def _validated_path(artifact: dict[str, Any], run_dir: Path) -> Path:
    path = Path(artifact["path"])
    if not path.is_absolute() and PurePosixPath(artifact["path"]).parts[:1] == ("stages",):
        path = run_dir / path
    elif not path.is_absolute():
        path = Path.cwd() / path
    if not path.is_file() or sha256_file(path) != artifact["sha256"]:
        raise AnalyzeRohInputError(f"artifact_missing_or_modified:{artifact['artifact_id']}")
    return path


def _read_fam(path: Path) -> list[tuple[str, str]]:
    rows = [line.split() for line in path.read_text(encoding="utf-8").splitlines() if line]
    if not rows or any(len(row) != 6 for row in rows):
        raise AnalyzeRohInputError(f"malformed_fam:{path.name}")
    pairs = [(row[0], row[1]) for row in rows]
    if len(pairs) != len(set(pairs)):
        raise AnalyzeRohInputError(f"duplicate_sample:{path.name}")
    return pairs


def _read_bim(path: Path) -> list[list[str]]:
    rows = [line.split() for line in path.read_text(encoding="utf-8").splitlines() if line]
    if not rows or any(len(row) != 6 for row in rows):
        raise AnalyzeRohInputError(f"malformed_bim:{path.name}")
    if len({row[1] for row in rows}) != len(rows):
        raise AnalyzeRohInputError(f"duplicate_variant:{path.name}")
    return rows


def _validate_dataset(
    dataset_id: str, paths: dict[str, Path], artifacts: dict[str, dict[str, Any]],
) -> tuple[dict[str, Any], list[tuple[str, str]], list[list[str]]]:
    descriptor = read_json(paths[f"{dataset_id}_dataset"])
    validate_json_document(descriptor, "plink_dataset.schema.json")
    if descriptor["dataset_id"] != dataset_id:
        raise AnalyzeRohInputError(f"dataset_id_mismatch:{dataset_id}")
    if paths[f"{dataset_id}_bed"].read_bytes()[:3] != b"\x6c\x1b\x01":
        raise AnalyzeRohInputError(f"invalid_bed:{dataset_id}")
    prefixes = {paths[f"{dataset_id}_{suffix}"].with_suffix("") for suffix in ("bed", "bim", "fam")}
    if len(prefixes) != 1:
        raise AnalyzeRohInputError(f"dataset_prefix_mismatch:{dataset_id}")
    fam, bim = _read_fam(paths[f"{dataset_id}_fam"]), _read_bim(paths[f"{dataset_id}_bim"])
    if descriptor["sample_count"] != len(fam) or descriptor["variant_count"] != len(bim):
        raise AnalyzeRohInputError(f"dataset_count_mismatch:{dataset_id}")
    for suffix in ("bed", "bim", "fam"):
        if descriptor["files"][suffix]["sha256"] != artifacts[f"{dataset_id}_{suffix}"]["sha256"]:
            raise AnalyzeRohInputError(f"dataset_checksum_mismatch:{dataset_id}")
    return descriptor, fam, bim


def _load_target_metadata(path: Path) -> dict[str, Any]:
    try:
        metadata = yaml.safe_load(path.read_text(encoding="utf-8"))
    except (OSError, yaml.YAMLError) as error:
        raise AnalyzeRohInputError("invalid_target_metadata") from error
    validate_json_document(metadata, "target_variant_metadata.schema.json")
    return metadata


def _sample_context(
    paths: dict[str, Path], all_pairs: set[tuple[str, str]],
) -> tuple[dict[tuple[str, str], str], dict[str, str], dict[str, set[str]]]:
    samples = validate_tsv_table(paths["samples_master"], "samples_master.schema.json")
    sample_by_pair = {(row["FID"], row["IID"]): row["SAMPLE_ID"] for row in samples.rows}
    if not all_pairs <= set(sample_by_pair):
        raise AnalyzeRohInputError("dataset_sample_absent_from_registry")
    genotypes = validate_tsv_table(paths["target_genotype_audit"], "target_genotype_audit.schema.json")
    genotype_by_sample = {row["SAMPLE_ID"]: row["GENOTYPE"] for row in genotypes.rows}
    target_sample_ids = {sample_by_pair[pair] for pair in all_pairs}
    if not target_sample_ids <= set(genotype_by_sample):
        raise AnalyzeRohInputError("target_genotype_coverage_mismatch")
    cohorts = validate_cohorts_frozen(paths["cohorts_frozen"], {row["SAMPLE_ID"] for row in samples.rows})
    roles: dict[str, set[str]] = defaultdict(set)
    for row in cohorts.rows:
        if row["INCLUDED"] and row["COHORT_ID"] in {
            "controls_unrelated", "target_carriers_independent", "family_noncarriers"
        }:
            roles[row["SAMPLE_ID"]].add(row["COHORT_ID"])
    for cohort_id in ("controls_unrelated", "target_carriers_independent"):
        units = [row["INDEPENDENT_UNIT_ID"] for row in cohorts.rows if row["INCLUDED"] and row["COHORT_ID"] == cohort_id]
        if any(unit is None for unit in units) or len(units) != len(set(units)):
            raise AnalyzeRohInputError(f"cohort_not_independent:{cohort_id}")
    return sample_by_pair, genotype_by_sample, roles


def _common_variants(first: list[list[str]], second: list[list[str]]) -> tuple[list[str], int]:
    second_by_id = {row[1]: row for row in second}
    common = []
    chromosomes = set()
    for row in first:
        other = second_by_id.get(row[1])
        if other is None:
            continue
        if (row[0], row[3], {row[4], row[5]}) != (other[0], other[3], {other[4], other[5]}):
            raise AnalyzeRohInputError("common_variant_definition_mismatch")
        common.append(row[1])
        chromosomes.add(row[0])
    return common, len(chromosomes)


def _resolve_plink(command: str | None) -> tuple[str, str]:
    executable = shutil.which(command) if command else None
    if executable is None:
        raise RohAnalysisError("plink_not_available")
    try:
        completed = subprocess.run([executable, "--version"], capture_output=True, check=False, text=True, timeout=10)
    except (OSError, subprocess.TimeoutExpired) as error:
        raise RohAnalysisError("plink_version_failed") from error
    lines = (completed.stdout or completed.stderr).strip().splitlines()
    if completed.returncode != 0 or not lines or "1.9" not in lines[0]:
        raise RohAnalysisError("plink19_required")
    return executable, lines[0][:200]


def execute(stage_inputs_path: Path, output_dir: Path) -> int:
    """Exécute les deux périmètres ROH sans consommer les étapes 13–15."""
    started_at, started_clock = utc_now(), monotonic()
    stage_inputs = read_json(stage_inputs_path)
    validate_json_document(stage_inputs, "stage_inputs.schema.json")
    run_dir = output_dir.parent.parent
    config = load_pipeline_config(run_dir / "config.resolved.yaml")
    parameters = _parameters(stage_inputs["parameters"])
    dataset_ids = (*GENOMEWIDE_DATASETS.values(), TARGET_DATASET_ID)
    artifact_ids = [
        "config_input_target_variant_metadata", "samples_master", "target_genotype_audit",
        "cohorts_frozen", "target_genetic_map",
    ]
    for dataset_id in dataset_ids:
        artifact_ids.extend(f"{dataset_id}_{suffix}" for suffix in ("bed", "bim", "fam", "dataset"))
    artifacts = {artifact_id: _artifact_by_id(stage_inputs, artifact_id) for artifact_id in artifact_ids}
    paths = {artifact_id: _validated_path(artifact, run_dir) for artifact_id, artifact in artifacts.items()}
    metadata = _load_target_metadata(paths["config_input_target_variant_metadata"])
    datasets = {dataset_id: _validate_dataset(dataset_id, paths, artifacts) for dataset_id in dataset_ids}
    descriptors = [value[0] for value in datasets.values()]
    if any(descriptor["assembly"] != metadata["assembly"] for descriptor in descriptors):
        raise AnalyzeRohInputError("dataset_assembly_mismatch")
    controls_bim = datasets["controls_unrelated_qc"][2]
    carriers_bim = datasets["target_carriers_independent_qc"][2]
    common_ids, autosome_count = _common_variants(controls_bim, carriers_bim)
    target_fam = datasets[TARGET_DATASET_ID][1]
    target_bim = datasets[TARGET_DATASET_ID][2]
    target_map = validate_tsv_table(paths["target_genetic_map"], "target_genetic_map.schema.json")
    target_rows = [row for row in target_map.rows if row["IS_TARGET_VARIANT"]]
    if len(target_rows) != 1 or target_rows[0]["VARIANT_ID"] != metadata["project_variant_id"]:
        raise AnalyzeRohInputError("target_map_mismatch")
    target_ids = {row[1] for row in target_bim}
    if metadata["project_variant_id"] not in target_ids:
        raise AnalyzeRohInputError("target_variant_absent_from_target_chromosome")
    all_pairs = set(target_fam)
    for _, fam, _ in datasets.values():
        all_pairs.update(fam)
    sample_by_pair, genotype_by_sample, roles = _sample_context(paths, all_pairs)
    executable, plink_version = _resolve_plink(config["tools"]["plink"])
    genomewide = {
        cohort_id: {
            "dataset_id": dataset_id,
            "prefix": paths[f"{dataset_id}_bed"].with_suffix(""),
            "fam_pairs": datasets[dataset_id][1],
        }
        for cohort_id, dataset_id in GENOMEWIDE_DATASETS.items()
    }
    target_dataset = {
        "dataset_id": TARGET_DATASET_ID,
        "prefix": paths[f"{TARGET_DATASET_ID}_bed"].with_suffix(""),
        "fam_pairs": target_fam,
    }
    publication = publish_roh(
        plink_executable=executable, genomewide_datasets=genomewide,
        target_dataset=target_dataset, common_variant_ids=common_ids,
        common_autosome_count=autosome_count, target_variant_count=len(target_bim),
        sample_by_pair=sample_by_pair, target_genotypes=genotype_by_sample,
        roles_by_sample=roles, target_variant_id=metadata["project_variant_id"],
        target_chromosome=str(metadata["chromosome"]), target_bp=metadata["position_bp"],
        output_dir=output_dir / "roh", parameters=parameters,
    )
    specifications = [
        ("roh_segments", publication.segments_path, "text/tab-separated-values", "roh_segments.schema.json", "sensitive_genetic"),
        ("roh_burden", publication.burden_path, "text/tab-separated-values", "roh_burden.schema.json", "sensitive_genetic"),
        ("target_in_roh", publication.target_path, "text/tab-separated-values", "target_in_roh.schema.json", "sensitive_genetic"),
        ("roh_cohort_summary", publication.cohort_summary_path, "text/tab-separated-values", "roh_cohort_summary.schema.json", "internal"),
        ("roh_analysis_summary", publication.analysis_summary_path, "application/json", "roh_analysis_summary.schema.json", "internal"),
    ]
    for native in publication.native_paths:
        specifications.append((f"roh_native_{native.name.replace('.', '_')}", native, "text/plain", None, "sensitive_genetic"))
    output_artifacts = [
        build_file_artifact(
            physical_path=path,
            published_path=f"{stage_inputs['published_output_dir']}/{path.relative_to(output_dir).as_posix()}",
            artifact_id=artifact_id, artifact_type=artifact_id, media_type=media_type,
            producer_stage=stage_inputs["stage_name"], producer_signature=stage_inputs["signature"],
            schema_name=schema_name, schema_version="1.0.0" if schema_name else None,
            assembly=metadata["assembly"], sample_set_id=None, variant_set_id=None,
            sensitivity=sensitivity,
        )
        for artifact_id, path, media_type, schema_name, sensitivity in specifications
    ]
    stage_outputs = {
        "schema_version": "1.0.0", "run_id": stage_inputs["run_id"],
        "stage_id": stage_inputs["stage_id"], "stage_name": stage_inputs["stage_name"],
        "signature": stage_inputs["signature"], "artifacts": output_artifacts,
    }
    validate_json_document(stage_outputs, "stage_outputs.schema.json")
    atomic_write_json(output_dir / "stage_outputs.json", stage_outputs)
    warnings = []
    for scope, status in publication.scope_statuses.items():
        if status != "EVALUATED":
            warnings.append({"code": f"{scope.lower()}_{status.lower()}", "count": 1})
    target_table = validate_tsv_table(publication.target_path, "target_in_roh.schema.json")
    heterozygous_warnings = sum(row["INTERPRETATION_STATUS"] == "IN_ROH_HETEROZYGOUS_WARNING" for row in target_table.rows)
    if heterozygous_warnings:
        warnings.append({"code": "heterozygous_target_inside_roh", "count": heterozygous_warnings})
    audit = {
        "schema_version": "1.0.0", "run_id": stage_inputs["run_id"],
        "stage_id": stage_inputs["stage_id"], "stage_name": stage_inputs["stage_name"],
        "method_id": "plink19_array_roh_secondary_v1", "signature": stage_inputs["signature"],
        "started_at": started_at, "completed_at": utc_now(),
        "duration_seconds": monotonic() - started_clock,
        "inputs": list(artifacts.values()), "outputs": output_artifacts,
        "parameters": parameters,
        "tools": [{"tool": "plink", "configured": config["tools"]["plink"], "version": plink_version}],
        "counts": {
            "genomewide_common_variants": len(common_ids), "genomewide_autosomes": autosome_count,
            "target_chromosome_variants": len(target_bim), "segments": publication.segment_count,
            "target_in_roh": publication.target_in_roh_count,
        },
        "metrics": {
            "analysis_role": "SECONDARY_DESCRIPTIVE_AUTOZYGOSITY",
            "scope_statuses": publication.scope_statuses,
            "f_roh_calculated": parameters["autosomal_denominator_kb"] is not None,
            "feeds_founder_haplotype": False, "feeds_variant_age": False,
            "consumes_local_ld": False,
        },
        "exclusions": [], "warnings": warnings,
        "checks": [
            {"check": "input_artifact_integrity", "status": "PASS"},
            {"check": "common_genomewide_variant_universe", "status": "PASS"},
            {"check": "independent_burden_cohorts", "status": "PASS"},
            {"check": "individual_target_intersection", "status": "PASS"},
            {"check": "stages_13_15_not_consumed", "status": "PASS"},
        ],
        "known_limits": [
            "Un ROH décrit une homozygotie individuelle et ne démontre pas un haplotype fondateur partagé entre familles.",
            "La densité ACPA, les trous entre marqueurs et les erreurs de génotypage conditionnent fortement les limites des ROH.",
            "Les paramètres PLINK sont exploratoires et doivent être soumis aux sensibilités de l'étape 17.",
            "Un génotype cible hétérozygote dans un ROH tolérant un hétérozygote exige une revue technique.",
            "F_ROH n'est calculé qu'avec un dénominateur autosomique explicite et sourcé.",
        ],
        "expected_visualizations": ["plot_roh_burden", "plot_roh_target_tracks", "plot_target_in_roh"],
        "manual_validation_required": True,
    }
    validate_json_document(audit, "stage_audit.schema.json")
    atomic_write_json(output_dir / "audit.json", audit)
    (output_dir / "checksums.sha256").write_text(
        "".join(f"{artifact['sha256']}  {artifact['path'].removeprefix(stage_inputs['published_output_dir'] + '/')}\n" for artifact in output_artifacts),
        encoding="utf-8",
    )
    return 0


def _classify_error(error: Exception) -> int:
    if isinstance(error, RohAnalysisError):
        return 3
    if isinstance(error, (AnalyzeRohInputError, DocumentValidationError, TableValidationError, OSError, ValueError, json.JSONDecodeError)):
        return 2
    return 5


def main(arguments: Sequence[str] | None = None) -> int:
    """Point d'entrée autonome respectant les codes V2."""
    parsed = _build_parser().parse_args(arguments)
    try:
        return execute(parsed.stage_inputs, parsed.output_dir)
    except Exception as error:
        sys.stderr.write(f"{error}\n")
        return _classify_error(error)


if __name__ == "__main__":
    raise SystemExit(main())
