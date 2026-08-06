"""Étape 12.1–12.5 : référence, harmonisation et entrées SHAPEIT5."""

from __future__ import annotations

import argparse
import json
import math
import sys
from collections import Counter
from pathlib import Path, PurePosixPath
from time import monotonic
from typing import Any, Sequence

from effet_fondateur.audit import atomic_write_json, read_json, sha256_file
from effet_fondateur.contracts import (
    DocumentValidationError,
    build_file_artifact,
    load_pipeline_config,
    validate_json_document,
    validate_tsv_table,
)
from effet_fondateur.orchestrator.state import utc_now
from effet_fondateur.phasing import (
    PreparedShapeit5Inputs,
    Shapeit5InputBlockError,
    Shapeit5InputError,
    Shapeit5InputExternalError,
    build_shapeit5_inputs,
)
from effet_fondateur.references import (
    ExtractedReferenceWindow,
    HarmonizedReferenceWindow,
    ReferenceCacheError,
    ReferenceCacheIntegrityError,
    ReferenceCacheOfflineMiss,
    ReferenceCatalogError,
    ReferenceHarmonizationError,
    ReferenceHarmonizationIntegrityError,
    ReferenceWindowError,
    ReferenceWindowIntegrityError,
    cache_reference_panel,
    extract_reference_window,
    harmonize_reference_window,
    resolve_reference_panel,
)


class PhaseTargetRegionInputError(ValueError):
    """Signale une entrée ou une configuration invalide de l'étape 12."""


class PhaseTargetRegionExternalError(RuntimeError):
    """Signale un échec d'acquisition ou de bcftools."""


class PhaseTargetRegionBlockError(RuntimeError):
    """Signale une référence ou une harmonisation scientifiquement inutilisable."""


def _build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(prog="phase-target-region")
    parser.add_argument("--stage-inputs", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    return parser


def _positive_integer(parameters: dict[str, Any], name: str, default: int) -> int:
    value = parameters.get(name, default)
    if isinstance(value, bool) or not isinstance(value, int) or value <= 0:
        raise PhaseTargetRegionInputError(f"invalid_parameter:{name}")
    return value


def _positive_number(parameters: dict[str, Any], name: str, default: float) -> float:
    value = parameters.get(name, default)
    if (
        isinstance(value, bool)
        or not isinstance(value, (int, float))
        or not math.isfinite(value)
        or value <= 0
    ):
        raise PhaseTargetRegionInputError(f"invalid_parameter:{name}")
    return float(value)


def _boolean(parameters: dict[str, Any], name: str, default: bool) -> bool:
    value = parameters.get(name, default)
    if not isinstance(value, bool):
        raise PhaseTargetRegionInputError(f"invalid_parameter:{name}")
    return value


def _nonempty_string(parameters: dict[str, Any], name: str) -> str:
    value = parameters.get(name)
    if not isinstance(value, str) or not value.strip():
        raise PhaseTargetRegionInputError(f"invalid_parameter:{name}")
    return value


def _parameters(parameters: dict[str, Any]) -> dict[str, Any]:
    return {
        "reference_panel_id": _nonempty_string(parameters, "reference_panel_id"),
        "reference_cache_dir": _nonempty_string(parameters, "reference_cache_dir"),
        "reference_cache_offline": _boolean(
            parameters, "reference_cache_offline", False
        ),
        "reference_lock_timeout_seconds": _positive_number(
            parameters, "reference_lock_timeout_seconds", 60
        ),
        "reference_download_timeout_seconds": _positive_integer(
            parameters, "reference_download_timeout_seconds", 7_200
        ),
        "reference_download_chunk_size": _positive_integer(
            parameters, "reference_download_chunk_size", 1_048_576
        ),
        "reference_extract_timeout_seconds": _positive_number(
            parameters, "reference_extract_timeout_seconds", 7_200
        ),
        "reference_extract_threads": _positive_integer(
            parameters, "reference_extract_threads", 1
        ),
        "minimum_common_variants": _positive_integer(
            parameters, "minimum_common_variants", 1
        ),
        "shapeit5_input_timeout_seconds": _positive_number(
            parameters, "shapeit5_input_timeout_seconds", 300
        ),
    }


def _resolve_input_path(artifact: dict[str, Any], run_dir: Path) -> Path:
    artifact_path = Path(artifact["path"])
    if artifact_path.is_absolute():
        return artifact_path
    if PurePosixPath(artifact["path"]).parts[:1] == ("stages",):
        return run_dir / artifact_path
    return Path.cwd() / artifact_path


def _artifact_by_id(stage_inputs: dict[str, Any], artifact_id: str) -> dict[str, Any]:
    matches = [
        artifact
        for artifact in stage_inputs["artifacts"]
        if artifact["artifact_id"] == artifact_id
    ]
    if len(matches) != 1:
        raise PhaseTargetRegionInputError(f"{artifact_id}_missing_or_ambiguous")
    return matches[0]


def _validated_artifact_path(artifact: dict[str, Any], run_dir: Path) -> Path:
    path = _resolve_input_path(artifact, run_dir)
    if path.is_symlink() or not path.is_file() or sha256_file(path) != artifact["sha256"]:
        raise PhaseTargetRegionInputError("declared_input_modified")
    return path


def _study_bounds(bim_path: Path, chromosome: int) -> tuple[int, int]:
    positions: list[int] = []
    for line in bim_path.read_text(encoding="utf-8").splitlines():
        if not line:
            continue
        fields = line.split()
        if len(fields) != 6:
            raise PhaseTargetRegionInputError("study_bim_malformed")
        try:
            observed_chromosome = int(fields[0].removeprefix("chr"))
            position_bp = int(fields[3])
        except ValueError as error:
            raise PhaseTargetRegionInputError("study_bim_coordinate_invalid") from error
        if observed_chromosome != chromosome or position_bp <= 0:
            raise PhaseTargetRegionInputError("study_bim_coordinate_mismatch")
        positions.append(position_bp)
    if len(positions) < 2:
        raise PhaseTargetRegionBlockError("study_region_has_fewer_than_two_variants")
    return min(positions), max(positions)


def _artifact(
    *,
    physical_path: Path,
    relative_path: Path,
    stage_inputs: dict[str, Any],
    artifact_id: str,
    artifact_type: str,
    media_type: str,
    schema_name: str | None,
    assembly: str,
    sample_set_id: str | None,
    variant_set_id: str | None,
    sensitivity: str,
) -> dict[str, Any]:
    published_dir = PurePosixPath(stage_inputs["published_output_dir"])
    return build_file_artifact(
        physical_path=physical_path,
        published_path=str(published_dir / relative_path.as_posix()),
        artifact_id=artifact_id,
        artifact_type=artifact_type,
        media_type=media_type,
        producer_stage=f"{stage_inputs['stage_id']}_{stage_inputs['stage_name']}",
        producer_signature=stage_inputs["signature"],
        schema_name=schema_name,
        schema_version="1.0.0" if schema_name else None,
        assembly=assembly,
        sample_set_id=sample_set_id,
        variant_set_id=variant_set_id,
        sensitivity=sensitivity,
    )


def _output_artifacts(
    stage_inputs: dict[str, Any],
    output_dir: Path,
    reference_window: ExtractedReferenceWindow,
    harmonized: HarmonizedReferenceWindow,
    assembly: str,
    reference_sample_set_id: str,
    reference_variant_set_id: str,
    study_sample_set_id: str,
    study_variant_set_id: str,
) -> list[dict[str, Any]]:
    specifications = (
        (
            "reference_window_vcf",
            "reference_window_vcf",
            reference_window.vcf_path.relative_to(output_dir),
            "application/gzip",
            None,
            reference_sample_set_id,
            reference_variant_set_id,
            "public",
        ),
        (
            "reference_window_index",
            "reference_window_index",
            reference_window.index_path.relative_to(output_dir),
            "application/octet-stream",
            None,
            reference_sample_set_id,
            reference_variant_set_id,
            "public",
        ),
        (
            "reference_window_manifest",
            "reference_window_manifest",
            reference_window.manifest_path.relative_to(output_dir),
            "application/json",
            "reference_window_manifest.schema.json",
            reference_sample_set_id,
            reference_variant_set_id,
            "public",
        ),
        (
            "harmonized_reference_vcf",
            "harmonized_reference_vcf",
            harmonized.vcf_path.relative_to(output_dir),
            "application/gzip",
            None,
            reference_sample_set_id,
            reference_variant_set_id,
            "public",
        ),
        (
            "harmonized_reference_index",
            "harmonized_reference_index",
            harmonized.index_path.relative_to(output_dir),
            "application/octet-stream",
            None,
            reference_sample_set_id,
            reference_variant_set_id,
            "public",
        ),
        (
            "variant_harmonization",
            "variant_harmonization",
            harmonized.harmonization_path.relative_to(output_dir),
            "text/tab-separated-values",
            "variant_harmonization.schema.json",
            study_sample_set_id,
            study_variant_set_id,
            "sensitive_genetic",
        ),
        (
            "reference_harmonization_manifest",
            "reference_harmonization_manifest",
            harmonized.manifest_path.relative_to(output_dir),
            "application/json",
            "reference_harmonization_manifest.schema.json",
            study_sample_set_id,
            study_variant_set_id,
            "internal",
        ),
    )
    return [
        _artifact(
            physical_path=output_dir / relative_path,
            relative_path=relative_path,
            stage_inputs=stage_inputs,
            artifact_id=artifact_id,
            artifact_type=artifact_type,
            media_type=media_type,
            schema_name=schema_name,
            assembly=assembly,
            sample_set_id=sample_set_id,
            variant_set_id=variant_set_id,
            sensitivity=sensitivity,
        )
        for (
            artifact_id,
            artifact_type,
            relative_path,
            media_type,
            schema_name,
            sample_set_id,
            variant_set_id,
            sensitivity,
        ) in specifications
    ]


def _shapeit5_input_artifacts(
    stage_inputs: dict[str, Any],
    output_dir: Path,
    prepared: PreparedShapeit5Inputs,
    assembly: str,
    sample_set_id: str,
    variant_set_id: str,
) -> list[dict[str, Any]]:
    specifications = (
        (
            "shapeit5_study_vcf",
            "shapeit5_study_vcf",
            prepared.study_vcf_path,
            "application/gzip",
            None,
        ),
        (
            "shapeit5_study_index",
            "shapeit5_study_index",
            prepared.study_index_path,
            "application/octet-stream",
            None,
        ),
        (
            "shapeit5_genetic_map",
            "shapeit5_genetic_map",
            prepared.genetic_map_path,
            "application/gzip",
            None,
        ),
        (
            "shapeit5_pedigree",
            "shapeit5_pedigree",
            prepared.pedigree_path,
            "text/tab-separated-values",
            None,
        ),
        (
            "shapeit5_variant_selection",
            "shapeit5_variant_selection",
            prepared.variant_selection_path,
            "text/tab-separated-values",
            "shapeit5_variant_selection.schema.json",
        ),
        (
            "shapeit5_sample_mapping",
            "shapeit5_sample_mapping",
            prepared.sample_mapping_path,
            "text/tab-separated-values",
            "shapeit5_sample_mapping.schema.json",
        ),
        (
            "shapeit5_inputs_manifest",
            "shapeit5_inputs_manifest",
            prepared.manifest_path,
            "application/json",
            "shapeit5_inputs_manifest.schema.json",
        ),
    )
    return [
        _artifact(
            physical_path=physical_path,
            relative_path=physical_path.relative_to(output_dir),
            stage_inputs=stage_inputs,
            artifact_id=artifact_id,
            artifact_type=artifact_type,
            media_type=media_type,
            schema_name=schema_name,
            assembly=assembly,
            sample_set_id=sample_set_id,
            variant_set_id=variant_set_id,
            sensitivity="sensitive_genetic",
        )
        for (
            artifact_id,
            artifact_type,
            physical_path,
            media_type,
            schema_name,
        ) in specifications
    ]


def execute(stage_inputs_path: Path, output_dir: Path) -> int:
    """Prépare la référence et les entrées validées sans lancer SHAPEIT5."""
    started_at = utc_now()
    started_clock = monotonic()
    stage_inputs = read_json(stage_inputs_path)
    validate_json_document(stage_inputs, "stage_inputs.schema.json")
    run_dir = output_dir.parent.parent
    config = load_pipeline_config(run_dir / "config.resolved.yaml")
    parameters = _parameters(stage_inputs["parameters"])
    artifact_ids = (
        "config_input_target_variant_metadata",
        "config_input_reference_panel_catalog",
        "target_region_bim",
        "target_region_bed",
        "target_region_fam",
        "target_genetic_map",
        "phasing_input_manifest",
        "samples_master",
    )
    input_artifacts = {
        artifact_id: _artifact_by_id(stage_inputs, artifact_id)
        for artifact_id in artifact_ids
    }
    paths = {
        artifact_id: _validated_artifact_path(artifact, run_dir)
        for artifact_id, artifact in input_artifacts.items()
    }
    phasing_manifest = read_json(paths["phasing_input_manifest"])
    validate_json_document(phasing_manifest, "phasing_input_manifest.schema.json")
    assembly = phasing_manifest["assembly"]
    chromosome = phasing_manifest["chromosome"]
    start_bp, end_bp = _study_bounds(paths["target_region_bim"], chromosome)
    resolved = resolve_reference_panel(
        panel_id=parameters["reference_panel_id"],
        assembly=assembly,
        chromosome=chromosome,
        catalog_path=paths["config_input_reference_panel_catalog"],
    )
    cache_root = Path(parameters["reference_cache_dir"]).expanduser()
    if not cache_root.is_absolute():
        cache_root = Path.cwd() / cache_root
    cached = cache_reference_panel(
        resolved,
        cache_root,
        offline=parameters["reference_cache_offline"],
        lock_timeout_seconds=parameters["reference_lock_timeout_seconds"],
        download_timeout_seconds=parameters["reference_download_timeout_seconds"],
        chunk_size=parameters["reference_download_chunk_size"],
    )
    bcftools_command = config["tools"]["bcftools"]
    if not isinstance(bcftools_command, str) or not bcftools_command:
        raise PhaseTargetRegionInputError("bcftools_not_configured")
    reference_window = extract_reference_window(
        resolved,
        cached,
        output_dir / "reference_window",
        start_bp=start_bp,
        end_bp=end_bp,
        bcftools_command=bcftools_command,
        threads=parameters["reference_extract_threads"],
        timeout_seconds=parameters["reference_extract_timeout_seconds"],
    )
    harmonized = harmonize_reference_window(
        reference_window,
        paths["target_region_bim"],
        paths["phasing_input_manifest"],
        paths["config_input_target_variant_metadata"],
        output_dir / "reference_harmonization",
        bcftools_command=bcftools_command,
        threads=parameters["reference_extract_threads"],
        timeout_seconds=parameters["reference_extract_timeout_seconds"],
        minimum_common_variants=parameters["minimum_common_variants"],
    )
    harmonization_manifest = read_json(harmonized.manifest_path)
    harmonization_table = validate_tsv_table(
        harmonized.harmonization_path, "variant_harmonization.schema.json"
    )
    status_counts = Counter(
        row["HARMONIZATION_STATUS"] for row in harmonization_table.rows
    )
    ineligible_counts = Counter(
        row["HARMONIZATION_STATUS"]
        for row in harmonization_table.rows
        if not row["COMMON_PHASE_ELIGIBLE"]
    )
    reference_sample_set_id = (
        f"reference_{resolved.panel_id}_{resolved.release_id}_chr{chromosome}"
    )
    reference_variant_set_id = (
        f"reference_variants_{harmonization_manifest['files']['vcf']['sha256'][:16]}"
    )
    prepared_shapeit5 = build_shapeit5_inputs(
        study_bed_path=paths["target_region_bed"],
        study_bim_path=paths["target_region_bim"],
        study_fam_path=paths["target_region_fam"],
        phasing_manifest_path=paths["phasing_input_manifest"],
        target_genetic_map_path=paths["target_genetic_map"],
        samples_master_path=paths["samples_master"],
        target_metadata_path=paths["config_input_target_variant_metadata"],
        harmonization_path=harmonized.harmonization_path,
        harmonization_manifest_path=harmonized.manifest_path,
        reference_vcf_path=harmonized.vcf_path,
        reference_index_path=harmonized.index_path,
        reference_variant_set_id=reference_variant_set_id,
        stage_root=output_dir,
        output_dir=output_dir / "shapeit5_inputs",
        plink_command=config["tools"]["plink"],
        bcftools_command=bcftools_command,
        timeout_seconds=parameters["shapeit5_input_timeout_seconds"],
    )
    output_artifacts = _output_artifacts(
        stage_inputs,
        output_dir,
        reference_window,
        harmonized,
        assembly,
        reference_sample_set_id,
        reference_variant_set_id,
        phasing_manifest["sample_set_id"],
        phasing_manifest["variant_set_id"],
    )
    shapeit5_manifest = read_json(prepared_shapeit5.manifest_path)
    output_artifacts.extend(
        _shapeit5_input_artifacts(
            stage_inputs,
            output_dir,
            prepared_shapeit5,
            assembly,
            phasing_manifest["sample_set_id"],
            shapeit5_manifest["study_variant_set_id"],
        )
    )
    stage_outputs = {
        "schema_version": "1.0.0",
        "run_id": stage_inputs["run_id"],
        "stage_id": stage_inputs["stage_id"],
        "stage_name": stage_inputs["stage_name"],
        "signature": stage_inputs["signature"],
        "artifacts": output_artifacts,
    }
    validate_json_document(stage_outputs, "stage_outputs.schema.json")
    atomic_write_json(output_dir / "stage_outputs.json", stage_outputs)
    audit = {
        "schema_version": "1.0.0",
        "run_id": stage_inputs["run_id"],
        "stage_id": stage_inputs["stage_id"],
        "stage_name": stage_inputs["stage_name"],
        "method_id": "reference_and_shapeit5_input_preparation_v1",
        "signature": stage_inputs["signature"],
        "started_at": started_at,
        "completed_at": utc_now(),
        "duration_seconds": monotonic() - started_clock,
        "inputs": list(input_artifacts.values()),
        "outputs": output_artifacts,
        "parameters": parameters,
        "tools": [
            {
                "tool": "bcftools",
                "configured": bcftools_command,
                "version": harmonization_manifest["tool"]["version"],
            },
            {
                "tool": "plink",
                "configured": config["tools"]["plink"],
                "version": shapeit5_manifest["tools"]["plink"],
            },
        ],
        "counts": {
            "reference_samples": harmonization_manifest["sample_count"],
            "reference_window_variants": reference_window.variant_count,
            "reference_harmonized_variants": harmonized.reference_variant_count,
            "study_variants": harmonized.study_variant_count,
            "common_variants": harmonized.common_variant_count,
            "harmonization_statuses": dict(sorted(status_counts.items())),
            "shapeit5_study_variants": prepared_shapeit5.study_variant_count,
            "shapeit5_pedigree_records": prepared_shapeit5.pedigree_record_count,
        },
        "metrics": {
            "panel_id": resolved.panel_id,
            "release_id": resolved.release_id,
            "cache_status": cached.status,
            "cache_entry_key": cached.cache_entry_key,
            "window_start_bp": start_bp,
            "window_end_bp": end_bp,
            "target_variant_status": harmonization_manifest["checks"][
                "target_variant"
            ],
            "shapeit5_target_role": shapeit5_manifest["target_variant_role"],
            "shapeit5_input_region": shapeit5_manifest["input_region"],
        },
        "exclusions": [
            {
                "scope": "common_phase",
                "reason": status,
                "count": count,
            }
            for status, count in sorted(ineligible_counts.items())
        ],
        "warnings": (
            [
                {
                    "code": "target_variant_absent_from_reference",
                    "count": 1,
                }
            ]
            if harmonization_manifest["checks"]["target_variant"] == "STUDY_ONLY"
            else []
        ),
        "checks": [
            {"check": "reference_catalog", "status": "PASS"},
            {"check": "reference_cache_integrity", "status": "PASS"},
            {"check": "reference_window", "status": "PASS"},
            {"check": "sample_set_and_order", "status": "PASS"},
            {"check": "canonical_variant_identity", "status": "PASS"},
            {"check": "target_variant", "status": "PASS"},
            {"check": "shapeit5_sample_identity_and_order", "status": "PASS"},
            {"check": "shapeit5_pedigree_identity", "status": "PASS"},
            {"check": "shapeit5_inputs", "status": "PASS"},
        ],
        "known_limits": [
            "Aucun gauche-alignement dépendant d'une FASTA n'est effectué ; les entrées GRCh38 doivent déjà utiliser une représentation minimale.",
            "Aucune correction automatique de brin n'est autorisée.",
            "Cette étape construit les entrées validées mais ne lance pas SHAPEIT5.",
        ],
        "expected_visualizations": [],
        "manual_validation_required": False,
    }
    validate_json_document(audit, "stage_audit.schema.json")
    atomic_write_json(output_dir / "audit.json", audit)
    (output_dir / "checksums.sha256").write_text(
        "".join(
            f"{artifact['sha256']}  "
            f"{artifact['path'].removeprefix(stage_inputs['published_output_dir'] + '/')}\n"
            for artifact in output_artifacts
        ),
        encoding="utf-8",
    )
    return 0


def _classify_error(error: Exception) -> int:
    message = str(error)
    if isinstance(
        error,
        (
            PhaseTargetRegionBlockError,
            ReferenceWindowIntegrityError,
            ReferenceHarmonizationIntegrityError,
            Shapeit5InputBlockError,
        ),
    ):
        return 4
    if isinstance(error, ReferenceCacheIntegrityError):
        return 4
    if isinstance(error, ReferenceCacheOfflineMiss):
        return 3
    if isinstance(error, (ReferenceWindowError, ReferenceHarmonizationError)) and message.startswith("bcftools_"):
        return 3
    if isinstance(error, Shapeit5InputExternalError):
        return 3
    if isinstance(error, ReferenceCacheError):
        if message.startswith("reference_download_") or message.endswith("lock_timeout"):
            return 3
        return 2
    if isinstance(
        error,
        (
            PhaseTargetRegionInputError,
            Shapeit5InputError,
            ReferenceCatalogError,
            DocumentValidationError,
            OSError,
            ValueError,
            json.JSONDecodeError,
        ),
    ):
        return 2
    return 5


def main(arguments: Sequence[str] | None = None) -> int:
    """Point d'entrée autonome respectant les codes de retour V2."""
    parser = _build_parser()
    parsed_arguments = parser.parse_args(arguments)
    try:
        return execute(parsed_arguments.stage_inputs, parsed_arguments.output_dir)
    except Exception as error:
        return_code = _classify_error(error)
        sys.stderr.write(f"{error}\n")
        return return_code


if __name__ == "__main__":
    raise SystemExit(main())
