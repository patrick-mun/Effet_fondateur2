"""Étape 11 : préparation de la région cible et de sa carte génétique."""

from __future__ import annotations

import argparse
import bisect
import csv
import hashlib
import json
import shutil
import subprocess
import sys
from decimal import Decimal, InvalidOperation, localcontext
from pathlib import Path, PurePosixPath
from time import monotonic
from typing import Any, Sequence

import yaml

from effet_fondateur.audit import atomic_write_json, read_json, sha256_file
from effet_fondateur.contracts import (
    DocumentValidationError,
    TableValidationError,
    build_file_artifact,
    load_pipeline_config,
    validate_json_document,
    validate_tsv_table,
)
from effet_fondateur.orchestrator.state import utc_now


MAP_COLUMNS = (
    "DATASET_ID",
    "VARIANT_ORDER",
    "VARIANT_ID",
    "CHROMOSOME",
    "POSITION_BP",
    "POSITION_CM",
    "MAP_STATUS",
    "LEFT_MAP_BP",
    "LEFT_MAP_CM",
    "RIGHT_MAP_BP",
    "RIGHT_MAP_CM",
    "LOCAL_RATE_CM_PER_MB",
    "IS_TARGET_VARIANT",
)


class TargetRegionInputError(ValueError):
    """Signale une entrée ou une configuration régionale invalide."""


class ExternalToolError(RuntimeError):
    """Signale un échec PLINK ou une sortie native incohérente."""


class TargetRegionBlockError(RuntimeError):
    """Signale une région scientifiquement inutilisable pour le phasage."""


def _build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(prog="prepare-target-region")
    parser.add_argument("--stage-inputs", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    return parser


def _positive_integer(parameters: dict[str, Any], name: str, default: int) -> int:
    value = parameters.get(name, default)
    if isinstance(value, bool) or not isinstance(value, int) or value <= 0:
        raise TargetRegionInputError(f"invalid_parameter:{name}")
    return value


def _parameters(parameters: dict[str, Any]) -> dict[str, int]:
    resolved = {
        "window_left_bp": _positive_integer(parameters, "window_left_bp", 5_000_000),
        "window_right_bp": _positive_integer(parameters, "window_right_bp", 5_000_000),
        "max_interpolation_gap_bp": _positive_integer(
            parameters, "max_interpolation_gap_bp", 5_000_000
        ),
        "min_region_variants": _positive_integer(parameters, "min_region_variants", 2),
        "plink_timeout_seconds": _positive_integer(
            parameters, "plink_timeout_seconds", 300
        ),
    }
    if resolved["min_region_variants"] < 2:
        raise TargetRegionInputError("invalid_parameter:min_region_variants")
    return resolved


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
        raise TargetRegionInputError(f"{artifact_id}_missing_or_ambiguous")
    return matches[0]


def _validated_artifact_path(artifact: dict[str, Any], run_dir: Path) -> Path:
    path = _resolve_input_path(artifact, run_dir)
    if not path.is_file() or sha256_file(path) != artifact["sha256"]:
        raise TargetRegionInputError("declared_input_modified")
    return path


def _load_target_metadata(path: Path) -> dict[str, Any]:
    try:
        metadata = yaml.safe_load(path.read_text(encoding="utf-8"))
    except (OSError, yaml.YAMLError) as error:
        raise TargetRegionInputError("invalid_target_variant_metadata") from error
    if not isinstance(metadata, dict):
        raise TargetRegionInputError("invalid_target_variant_metadata")
    validate_json_document(metadata, "target_variant_metadata.schema.json")
    return metadata


def _read_fam(path: Path) -> list[list[str]]:
    rows = [line.split() for line in path.read_text(encoding="utf-8").splitlines() if line]
    if not rows or any(len(row) != 6 for row in rows):
        raise TargetRegionInputError(f"malformed_fam:{path.name}")
    if len({(row[0], row[1]) for row in rows}) != len(rows):
        raise TargetRegionInputError(f"duplicate_sample_in_fam:{path.name}")
    return rows


def _read_bim(path: Path) -> list[list[str]]:
    rows = [line.split() for line in path.read_text(encoding="utf-8").splitlines() if line]
    if not rows or any(len(row) != 6 for row in rows):
        raise TargetRegionInputError(f"malformed_bim:{path.name}")
    if len({row[1] for row in rows}) != len(rows):
        raise TargetRegionInputError(f"duplicate_variant_in_bim:{path.name}")
    return rows


def _validate_source_dataset(
    paths: dict[str, Path], artifacts: dict[str, dict[str, Any]]
) -> tuple[dict[str, Any], list[list[str]], list[list[str]]]:
    descriptor = read_json(paths["target_chromosome_all_qc_dataset"])
    validate_json_document(descriptor, "plink_dataset.schema.json")
    if paths["target_chromosome_all_qc_bed"].read_bytes()[:3] != b"\x6c\x1b\x01":
        raise TargetRegionInputError("invalid_bed_header")
    stems = {
        paths[f"target_chromosome_all_qc_{suffix}"].with_suffix("")
        for suffix in ("bed", "bim", "fam")
    }
    if len(stems) != 1:
        raise TargetRegionInputError("dataset_files_do_not_share_prefix")
    for suffix in ("bed", "bim", "fam"):
        artifact_id = f"target_chromosome_all_qc_{suffix}"
        if (
            descriptor["files"][suffix]["name"] != paths[artifact_id].name
            or descriptor["files"][suffix]["sha256"] != artifacts[artifact_id]["sha256"]
            or artifacts[artifact_id]["sample_set_id"] != descriptor["sample_set_id"]
            or artifacts[artifact_id]["variant_set_id"] != descriptor["variant_set_id"]
        ):
            raise TargetRegionInputError("dataset_descriptor_mismatch")
    fam_rows = _read_fam(paths["target_chromosome_all_qc_fam"])
    bim_rows = _read_bim(paths["target_chromosome_all_qc_bim"])
    if descriptor["sample_count"] != len(fam_rows) or descriptor["variant_count"] != len(
        bim_rows
    ):
        raise TargetRegionInputError("dataset_count_mismatch")
    return descriptor, fam_rows, bim_rows


def _decimal(raw_value: str, field: str) -> Decimal:
    try:
        value = Decimal(raw_value)
    except InvalidOperation as error:
        raise TargetRegionInputError(f"invalid_decimal:{field}") from error
    if not value.is_finite() or value < 0:
        raise TargetRegionInputError(f"invalid_decimal:{field}")
    return value


def _decimal_text(value: Decimal) -> str:
    normalized = value.normalize()
    if normalized == normalized.to_integral():
        return str(normalized.quantize(Decimal(1)))
    return format(normalized, "f").rstrip("0").rstrip(".")


def _load_genetic_map(
    path: Path, assembly: str, chromosome: int
) -> tuple[str, list[tuple[int, Decimal]]]:
    table = validate_tsv_table(path, "genetic_map.schema.json")
    map_ids = {row["MAP_ID"] for row in table.rows}
    if len(map_ids) != 1:
        raise TargetRegionInputError("genetic_map_multiple_map_ids")
    if {row["ASSEMBLY"] for row in table.rows} != {assembly}:
        raise TargetRegionInputError("genetic_map_assembly_mismatch")
    chromosome_rows = [
        row for row in table.rows if int(row["CHROMOSOME"]) == chromosome
    ]
    if len(chromosome_rows) < 2:
        raise TargetRegionBlockError("genetic_map_chromosome_has_fewer_than_two_points")
    points = [
        (int(row["POSITION_BP"]), _decimal(row["POSITION_CM"], "POSITION_CM"))
        for row in chromosome_rows
    ]
    for previous, current in zip(points, points[1:], strict=False):
        if current[0] <= previous[0]:
            raise TargetRegionInputError("genetic_map_positions_not_strictly_increasing")
        if current[1] < previous[1]:
            raise TargetRegionBlockError("genetic_map_cm_not_monotonic")
    return next(iter(map_ids)), points


def _validate_target_variant_qc(path: Path, target_variant_id: str) -> None:
    table = validate_tsv_table(path, "variant_qc_final.schema.json")
    matches = [
        row
        for row in table.rows
        if row["DATASET_ID"] == "target_chromosome_all_qc"
        and row["VARIANT_ID"] == target_variant_id
    ]
    if len(matches) != 1 or not matches[0]["PRESENT_AFTER_QC"]:
        raise TargetRegionBlockError("target_variant_not_retained_by_final_qc")


def _resolve_executable(command: str | None) -> str:
    if command is None:
        raise ExternalToolError("plink_not_configured")
    executable = shutil.which(command)
    if executable is None:
        raise ExternalToolError("plink_not_available")
    return executable


def _plink_version(executable: str) -> str | None:
    try:
        completed = subprocess.run(
            [executable, "--version"],
            capture_output=True,
            check=False,
            text=True,
            timeout=10,
        )
    except (OSError, subprocess.TimeoutExpired):
        return None
    output = (completed.stdout or completed.stderr).strip().splitlines()
    return output[0][:200] if output else None


def _run_plink(executable: str, arguments: list[str], timeout_seconds: int) -> None:
    try:
        completed = subprocess.run(
            [executable, *arguments],
            capture_output=True,
            check=False,
            text=True,
            timeout=timeout_seconds,
        )
    except (OSError, subprocess.TimeoutExpired) as error:
        raise ExternalToolError("plink_target_region_failed") from error
    if completed.returncode != 0:
        raise ExternalToolError(f"plink_target_region_failed:{completed.returncode}")


def _validate_region_bim(
    bim_rows: list[list[str]], metadata: dict[str, Any], minimum_variants: int
) -> None:
    if len(bim_rows) < minimum_variants:
        raise TargetRegionBlockError("insufficient_region_variants")
    chromosome = str(metadata["chromosome"])
    if any(row[0] != chromosome for row in bim_rows):
        raise TargetRegionInputError("region_contains_unexpected_chromosome")
    positions = [int(row[3]) for row in bim_rows]
    if positions != sorted(positions) or len(positions) != len(set(positions)):
        raise TargetRegionBlockError("region_positions_not_unique_and_ordered")
    for row in bim_rows:
        if row[4] not in {"A", "C", "G", "T"} or row[5] not in {"A", "C", "G", "T"}:
            raise TargetRegionBlockError("region_contains_non_acgt_allele")
        if row[4] == row[5]:
            raise TargetRegionBlockError("region_contains_identical_alleles")
    target_rows = [row for row in bim_rows if row[1] == metadata["project_variant_id"]]
    if len(target_rows) != 1:
        raise TargetRegionBlockError("target_variant_absent_from_region")
    target = target_rows[0]
    if int(target[3]) != metadata["position_bp"]:
        raise TargetRegionInputError("target_variant_position_mismatch")
    if {target[4], target[5]} != {metadata["ref"], metadata["alt"]}:
        raise TargetRegionInputError("target_variant_allele_mismatch")


def _interpolate_variant(
    position_bp: int,
    points: list[tuple[int, Decimal]],
    maximum_gap_bp: int,
) -> tuple[Decimal, str, int, Decimal, int, Decimal, Decimal | None]:
    map_positions = [point[0] for point in points]
    insertion_index = bisect.bisect_left(map_positions, position_bp)
    if insertion_index < len(points) and points[insertion_index][0] == position_bp:
        exact_bp, exact_cm = points[insertion_index]
        return exact_cm, "EXACT", exact_bp, exact_cm, exact_bp, exact_cm, None
    if insertion_index == 0 or insertion_index == len(points):
        raise TargetRegionBlockError("variant_outside_genetic_map_bounds")
    left_bp, left_cm = points[insertion_index - 1]
    right_bp, right_cm = points[insertion_index]
    gap_bp = right_bp - left_bp
    if gap_bp > maximum_gap_bp:
        raise TargetRegionBlockError("genetic_map_interpolation_gap_too_large")
    with localcontext() as context:
        context.prec = 28
        fraction = Decimal(position_bp - left_bp) / Decimal(gap_bp)
        position_cm = left_cm + fraction * (right_cm - left_cm)
        rate = (right_cm - left_cm) * Decimal(1_000_000) / Decimal(gap_bp)
    return position_cm, "INTERPOLATED", left_bp, left_cm, right_bp, right_cm, rate


def _build_map_rows(
    bim_rows: list[list[str]],
    points: list[tuple[int, Decimal]],
    maximum_gap_bp: int,
    target_variant_id: str,
) -> list[dict[str, str]]:
    rows = []
    previous_cm: Decimal | None = None
    for variant_order, bim_row in enumerate(bim_rows, start=1):
        position_bp = int(bim_row[3])
        position_cm, status, left_bp, left_cm, right_bp, right_cm, rate = (
            _interpolate_variant(position_bp, points, maximum_gap_bp)
        )
        if previous_cm is not None and position_cm < previous_cm:
            raise TargetRegionBlockError("interpolated_variant_map_not_monotonic")
        previous_cm = position_cm
        rows.append(
            {
                "DATASET_ID": "target_region",
                "VARIANT_ORDER": str(variant_order),
                "VARIANT_ID": bim_row[1],
                "CHROMOSOME": bim_row[0],
                "POSITION_BP": bim_row[3],
                "POSITION_CM": _decimal_text(position_cm),
                "MAP_STATUS": status,
                "LEFT_MAP_BP": str(left_bp),
                "LEFT_MAP_CM": _decimal_text(left_cm),
                "RIGHT_MAP_BP": str(right_bp),
                "RIGHT_MAP_CM": _decimal_text(right_cm),
                "LOCAL_RATE_CM_PER_MB": _decimal_text(rate) if rate is not None else "",
                "IS_TARGET_VARIANT": str(bim_row[1] == target_variant_id).lower(),
            }
        )
    return rows


def _write_region_bim(path: Path, bim_rows: list[list[str]], map_rows: list[dict[str, str]]) -> None:
    cm_by_variant = {row["VARIANT_ID"]: row["POSITION_CM"] for row in map_rows}
    path.write_text(
        "".join(
            "\t".join([row[0], row[1], cm_by_variant[row[1]], row[3], row[4], row[5]])
            + "\n"
            for row in bim_rows
        ),
        encoding="utf-8",
    )


def _write_tsv(path: Path, rows: list[dict[str, str]]) -> None:
    with path.open("w", encoding="utf-8", newline="") as output_file:
        writer = csv.DictWriter(
            output_file,
            fieldnames=MAP_COLUMNS,
            delimiter="\t",
            lineterminator="\n",
        )
        writer.writeheader()
        writer.writerows({column: row[column] for column in MAP_COLUMNS} for row in rows)


def _set_id(prefix: str, values: list[str]) -> str:
    digest = hashlib.sha256("".join(f"{value}\n" for value in values).encode()).hexdigest()
    return f"{prefix}_{digest[:12]}"


def _write_descriptor(
    path: Path,
    prefix: Path,
    source_descriptor: dict[str, Any],
    fam_rows: list[list[str]],
    bim_rows: list[list[str]],
) -> dict[str, Any]:
    descriptor = {
        **source_descriptor,
        "dataset_id": "target_region",
        "sample_set_id": _set_id(
            "target_region_samples", [f"{row[0]}\t{row[1]}" for row in fam_rows]
        ),
        "variant_set_id": _set_id(
            "target_region_variants", [row[1] for row in bim_rows]
        ),
        "sample_count": len(fam_rows),
        "variant_count": len(bim_rows),
        "source_format": "PLINK_TARGET_REGION",
        "files": {
            suffix: {
                "name": prefix.with_suffix(f".{suffix}").name,
                "sha256": sha256_file(prefix.with_suffix(f".{suffix}")),
            }
            for suffix in ("bed", "bim", "fam")
        },
    }
    validate_json_document(descriptor, "plink_dataset.schema.json")
    atomic_write_json(path, descriptor)
    return descriptor


def _dataset_artifacts(
    stage_inputs: dict[str, Any],
    prefix: Path,
    descriptor_path: Path,
    descriptor: dict[str, Any],
) -> list[dict[str, Any]]:
    producer_stage = f"{stage_inputs['stage_id']}_{stage_inputs['stage_name']}"
    published_dir = PurePosixPath(stage_inputs["published_output_dir"])
    artifacts = []
    for suffix, media_type in (
        ("bed", "application/octet-stream"),
        ("bim", "text/plain"),
        ("fam", "text/plain"),
    ):
        physical_path = prefix.with_suffix(f".{suffix}")
        artifacts.append(
            build_file_artifact(
                physical_path=physical_path,
                published_path=str(published_dir / physical_path.name),
                artifact_id=f"target_region_{suffix}",
                artifact_type=f"plink_{suffix}",
                media_type=media_type,
                producer_stage=producer_stage,
                producer_signature=stage_inputs["signature"],
                assembly=descriptor["assembly"],
                sample_set_id=descriptor["sample_set_id"],
                variant_set_id=descriptor["variant_set_id"],
                sensitivity="sensitive_genetic",
            )
        )
    artifacts.append(
        build_file_artifact(
            physical_path=descriptor_path,
            published_path=str(published_dir / descriptor_path.name),
            artifact_id="target_region_dataset",
            artifact_type="plink_dataset_descriptor",
            media_type="application/json",
            producer_stage=producer_stage,
            producer_signature=stage_inputs["signature"],
            schema_name="plink_dataset.schema.json",
            schema_version="1.0.0",
            assembly=descriptor["assembly"],
            sample_set_id=descriptor["sample_set_id"],
            variant_set_id=descriptor["variant_set_id"],
            sensitivity="sensitive_genetic",
        )
    )
    return artifacts


def _clean_plink_sidecars(prefix: Path) -> None:
    for suffix in ("log", "nosex", "hh"):
        path = prefix.with_suffix(f".{suffix}")
        if path.exists():
            path.unlink()


def execute(stage_inputs_path: Path, output_dir: Path) -> int:
    """Prépare une région cible ordonnée avec positions génétiques contrôlées."""
    started_at = utc_now()
    started_clock = monotonic()
    stage_inputs = read_json(stage_inputs_path)
    validate_json_document(stage_inputs, "stage_inputs.schema.json")
    run_dir = output_dir.parent.parent
    config = load_pipeline_config(run_dir / "config.resolved.yaml")
    parameters = _parameters(stage_inputs["parameters"])
    artifact_ids = (
        "config_input_target_variant_metadata",
        "config_input_genetic_map",
        "target_chromosome_all_qc_bed",
        "target_chromosome_all_qc_bim",
        "target_chromosome_all_qc_fam",
        "target_chromosome_all_qc_dataset",
        "variant_qc_final",
    )
    input_artifacts = {
        artifact_id: _artifact_by_id(stage_inputs, artifact_id)
        for artifact_id in artifact_ids
    }
    paths = {
        artifact_id: _validated_artifact_path(artifact, run_dir)
        for artifact_id, artifact in input_artifacts.items()
    }
    metadata = _load_target_metadata(paths["config_input_target_variant_metadata"])
    descriptor, source_fam, _ = _validate_source_dataset(paths, input_artifacts)
    if descriptor["assembly"] != metadata["assembly"]:
        raise TargetRegionInputError("target_dataset_assembly_mismatch")
    if (
        descriptor["scope"] != "target_chromosome"
        or descriptor["target_chromosome"] != metadata["chromosome"]
    ):
        raise TargetRegionInputError("target_dataset_scope_mismatch")
    _validate_target_variant_qc(
        paths["variant_qc_final"], metadata["project_variant_id"]
    )
    map_id, map_points = _load_genetic_map(
        paths["config_input_genetic_map"], metadata["assembly"], metadata["chromosome"]
    )
    window_start_bp = max(1, metadata["position_bp"] - parameters["window_left_bp"])
    window_end_bp = metadata["position_bp"] + parameters["window_right_bp"]
    output_dir.mkdir(parents=True, exist_ok=True)
    region_prefix = output_dir / "target_region"
    source_prefix = paths["target_chromosome_all_qc_bed"].with_suffix("")
    executable = _resolve_executable(config["tools"]["plink"])
    plink_version = _plink_version(executable)
    _run_plink(
        executable,
        [
            "--bfile",
            str(source_prefix),
            "--chr",
            str(metadata["chromosome"]),
            "--from-bp",
            str(window_start_bp),
            "--to-bp",
            str(window_end_bp),
            "--make-bed",
            "--out",
            str(region_prefix),
        ],
        parameters["plink_timeout_seconds"],
    )
    region_fam = _read_fam(region_prefix.with_suffix(".fam"))
    region_bim = _read_bim(region_prefix.with_suffix(".bim"))
    if region_fam != source_fam:
        raise ExternalToolError("region_sample_set_or_order_mismatch")
    _validate_region_bim(region_bim, metadata, parameters["min_region_variants"])
    map_rows = _build_map_rows(
        region_bim,
        map_points,
        parameters["max_interpolation_gap_bp"],
        metadata["project_variant_id"],
    )
    _write_region_bim(region_prefix.with_suffix(".bim"), region_bim, map_rows)
    rewritten_bim = _read_bim(region_prefix.with_suffix(".bim"))
    if [row[1] for row in rewritten_bim] != [row["VARIANT_ID"] for row in map_rows]:
        raise ExternalToolError("region_variant_order_changed")
    target_map_path = output_dir / "target_genetic_map.tsv"
    _write_tsv(target_map_path, map_rows)
    validate_tsv_table(target_map_path, "target_genetic_map.schema.json")
    descriptor_path = output_dir / "target_region.dataset.json"
    output_descriptor = _write_descriptor(
        descriptor_path, region_prefix, descriptor, region_fam, rewritten_bim
    )
    target_order = next(
        int(row["VARIANT_ORDER"])
        for row in map_rows
        if row["IS_TARGET_VARIANT"] == "true"
    )
    phasing_manifest_path = output_dir / "phasing_input_manifest.json"
    phasing_manifest = {
        "schema_version": "1.0.0",
        "format": "PLINK_BED_BIM_FAM_WITH_CM",
        "assembly": metadata["assembly"],
        "chromosome": metadata["chromosome"],
        "dataset_id": output_descriptor["dataset_id"],
        "sample_set_id": output_descriptor["sample_set_id"],
        "variant_set_id": output_descriptor["variant_set_id"],
        "sample_count": output_descriptor["sample_count"],
        "variant_count": output_descriptor["variant_count"],
        "target_variant_id": metadata["project_variant_id"],
        "target_variant_order": target_order,
        "map_id": map_id,
        "files": {
            "bed": region_prefix.with_suffix(".bed").name,
            "bim": region_prefix.with_suffix(".bim").name,
            "fam": region_prefix.with_suffix(".fam").name,
            "genetic_map": target_map_path.name,
        },
    }
    validate_json_document(phasing_manifest, "phasing_input_manifest.schema.json")
    atomic_write_json(phasing_manifest_path, phasing_manifest)
    report_path = output_dir / "target_region_report.json"
    report = {
        "schema_version": "1.0.0",
        "method_id": "bounded_reference_map_interpolation_v1",
        "map_id": map_id,
        "assembly": metadata["assembly"],
        "chromosome": metadata["chromosome"],
        "target_variant_id": metadata["project_variant_id"],
        "window_start_bp": window_start_bp,
        "window_end_bp": window_end_bp,
        "observed_start_bp": int(rewritten_bim[0][3]),
        "observed_end_bp": int(rewritten_bim[-1][3]),
        "sample_count": len(region_fam),
        "variant_count": len(rewritten_bim),
        "exact_map_variant_count": sum(row["MAP_STATUS"] == "EXACT" for row in map_rows),
        "interpolated_variant_count": sum(
            row["MAP_STATUS"] == "INTERPOLATED" for row in map_rows
        ),
        "minimum_position_cm": map_rows[0]["POSITION_CM"],
        "maximum_position_cm": map_rows[-1]["POSITION_CM"],
    }
    atomic_write_json(report_path, report)
    output_artifacts = _dataset_artifacts(
        stage_inputs, region_prefix, descriptor_path, output_descriptor
    )
    producer_stage = f"{stage_inputs['stage_id']}_{stage_inputs['stage_name']}"
    published_dir = PurePosixPath(stage_inputs["published_output_dir"])
    for artifact_id, artifact_type, path, schema_name, sensitivity in (
        (
            "target_genetic_map",
            "target_genetic_map",
            target_map_path,
            "target_genetic_map.schema.json",
            "sensitive_genetic",
        ),
        (
            "phasing_input_manifest",
            "phasing_input_manifest",
            phasing_manifest_path,
            "phasing_input_manifest.schema.json",
            "sensitive_genetic",
        ),
        (
            "target_region_report",
            "target_region_report",
            report_path,
            None,
            "internal",
        ),
    ):
        output_artifacts.append(
            build_file_artifact(
                physical_path=path,
                published_path=str(published_dir / path.name),
                artifact_id=artifact_id,
                artifact_type=artifact_type,
                media_type=(
                    "application/json"
                    if path.suffix == ".json"
                    else "text/tab-separated-values"
                ),
                producer_stage=producer_stage,
                producer_signature=stage_inputs["signature"],
                schema_name=schema_name,
                schema_version="1.0.0" if schema_name else None,
                assembly=metadata["assembly"],
                sample_set_id=output_descriptor["sample_set_id"],
                variant_set_id=output_descriptor["variant_set_id"],
                sensitivity=sensitivity,
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
        "method_id": "bounded_reference_map_interpolation_v1",
        "signature": stage_inputs["signature"],
        "started_at": started_at,
        "completed_at": utc_now(),
        "duration_seconds": monotonic() - started_clock,
        "inputs": list(input_artifacts.values()),
        "outputs": output_artifacts,
        "parameters": parameters,
        "tools": [
            {
                "tool": "plink",
                "configured": config["tools"]["plink"],
                "version": plink_version,
            }
        ],
        "counts": {
            "samples": len(region_fam),
            "variants": len(rewritten_bim),
            "map_points_on_target_chromosome": len(map_points),
            "exact_map_variants": report["exact_map_variant_count"],
            "interpolated_variants": report["interpolated_variant_count"],
        },
        "metrics": {
            "map_id": map_id,
            "window_start_bp": window_start_bp,
            "window_end_bp": window_end_bp,
            "target_variant_order": target_order,
        },
        "exclusions": [],
        "warnings": [],
        "checks": [
            {"check": "target_coordinate_and_alleles", "status": "PASS"},
            {"check": "region_variant_order", "status": "PASS"},
            {"check": "reference_map_assembly", "status": "PASS"},
            {"check": "bounded_interpolation", "status": "PASS"},
            {"check": "genetic_map_monotonicity", "status": "PASS"},
            {"check": "target_variant_retained", "status": "PASS"},
            {"check": "phasing_input_manifest", "status": "PASS"},
        ],
        "known_limits": [
            "La carte fournie est utilisée telle quelle comme carte de référence ; son origine et sa population doivent être validées scientifiquement en amont.",
            "Aucune extrapolation hors des ancres de carte et aucune approximation 1 Mb = 1 cM ne sont autorisées.",
            "Le manifest reste indépendant de l'outil tant que l'adaptateur de phasage de l'étape 12 n'est pas choisi.",
        ],
        "expected_visualizations": [
            "target_physical_marker_density",
            "target_genetic_marker_density",
            "local_recombination_rate",
            "target_variant_position",
        ],
        "manual_validation_required": False,
    }
    validate_json_document(audit, "stage_audit.schema.json")
    atomic_write_json(output_dir / "audit.json", audit)
    _clean_plink_sidecars(region_prefix)
    (output_dir / "checksums.sha256").write_text(
        "".join(
            f"{artifact['sha256']}  {Path(str(artifact['path'])).name}\n"
            for artifact in output_artifacts
        ),
        encoding="utf-8",
    )
    return 0


def main(arguments: Sequence[str] | None = None) -> int:
    """Point d'entrée autonome respectant les codes de retour V2."""
    parser = _build_parser()
    parsed_arguments = parser.parse_args(arguments)
    try:
        return execute(parsed_arguments.stage_inputs, parsed_arguments.output_dir)
    except ExternalToolError as error:
        sys.stderr.write(f"{error}\n")
        return 3
    except TargetRegionBlockError as error:
        sys.stderr.write(f"{error}\n")
        return 4
    except (
        TargetRegionInputError,
        DocumentValidationError,
        TableValidationError,
        OSError,
        ValueError,
        json.JSONDecodeError,
    ) as error:
        parser.error(str(error))
    return 2


if __name__ == "__main__":
    raise SystemExit(main())
