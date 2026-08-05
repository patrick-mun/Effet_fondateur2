"""Étape 06 : construction du panel autosomique indépendant de marqueurs."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import shutil
import statistics
import subprocess
import sys
from collections import Counter, defaultdict
from pathlib import Path, PurePosixPath
from time import monotonic
from typing import Any, Sequence

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


PRUNED_VARIANT_COLUMNS = (
    "VARIANT_ID",
    "CHROMOSOME",
    "POSITION_BP",
    "MAF",
    "PRELIMINARY_ALERT_CODES",
    "IN_COMPLEX_REGION",
    "PANEL_STATUS",
    "EXCLUSION_CODE",
)
COVERAGE_COLUMNS = (
    "CHROMOSOME",
    "INPUT_VARIANT_COUNT",
    "MAF_ELIGIBLE_COUNT",
    "PRUNED_MARKER_COUNT",
    "MARKER_FRACTION",
    "SPAN_BP",
    "MEDIAN_GAP_BP",
    "MAX_GAP_BP",
    "COVERAGE_STATUS",
    "REASON_CODE",
)
LD_BIN_COLUMNS = (
    "BIN_ID",
    "R2_LOWER_BOUND",
    "R2_UPPER_BOUND",
    "UPPER_INCLUSIVE",
    "PAIR_COUNT",
)


class KinshipPanelInputError(ValueError):
    """Signale une entrée ou un paramètre invalide."""


class ExternalToolError(RuntimeError):
    """Signale un échec ou une sortie invalide de PLINK."""


class KinshipPanelBlockError(RuntimeError):
    """Signale un panel scientifiquement insuffisant."""


def _build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(prog="build-kinship-panel")
    parser.add_argument("--stage-inputs", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    return parser


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
        raise KinshipPanelInputError(f"{artifact_id}_missing_or_ambiguous")
    return matches[0]


def _validated_artifact_path(artifact: dict[str, Any], run_dir: Path) -> Path:
    physical_path = _resolve_input_path(artifact, run_dir)
    if not physical_path.is_file() or sha256_file(physical_path) != artifact["sha256"]:
        raise KinshipPanelInputError("declared_input_modified")
    return physical_path


def _read_fam(path: Path) -> list[list[str]]:
    rows = [line.split() for line in path.read_text(encoding="utf-8").splitlines() if line]
    if not rows or any(len(row) != 6 for row in rows):
        raise KinshipPanelInputError(f"malformed_fam:{path.name}")
    return rows


def _read_bim(path: Path) -> list[list[str]]:
    rows = [line.split() for line in path.read_text(encoding="utf-8").splitlines() if line]
    if not rows or any(len(row) != 6 for row in rows):
        raise KinshipPanelInputError(f"malformed_bim:{path.name}")
    if len({row[1] for row in rows}) != len(rows):
        raise KinshipPanelInputError("duplicate_variant_id_in_pre_qc_bim")
    if any(not row[0].isdigit() or not 1 <= int(row[0]) <= 22 for row in rows):
        raise KinshipPanelInputError("non_autosomal_variant_in_pre_qc_bim")
    return rows


def _read_whitespace_table(path: Path, *, allow_empty: bool = False) -> list[dict[str, str]]:
    try:
        lines = [line.split() for line in path.read_text(encoding="utf-8").splitlines() if line]
    except OSError as error:
        raise ExternalToolError(f"missing_plink_report:{path.name}") from error
    if not lines or (len(lines) < 2 and not allow_empty):
        raise ExternalToolError(f"empty_plink_report:{path.name}")
    header = lines[0]
    if any(len(row) != len(header) for row in lines[1:]):
        raise ExternalToolError(f"malformed_plink_report:{path.name}")
    return [dict(zip(header, row, strict=True)) for row in lines[1:]]


def _read_variant_list(path: Path, *, allow_empty: bool = False) -> list[str]:
    try:
        variants = [line.strip() for line in path.read_text(encoding="utf-8").splitlines() if line]
    except OSError as error:
        raise ExternalToolError(f"missing_plink_report:{path.name}") from error
    if not variants and not allow_empty:
        raise ExternalToolError(f"empty_plink_report:{path.name}")
    if len(variants) != len(set(variants)):
        raise ExternalToolError(f"duplicate_variant_in_report:{path.name}")
    return variants


def _positive_integer(parameters: dict[str, Any], name: str, default: int) -> int:
    value = parameters.get(name, default)
    if isinstance(value, bool) or not isinstance(value, int) or value <= 0:
        raise KinshipPanelInputError(f"invalid_parameter:{name}")
    return value


def _bounded_float(
    parameters: dict[str, Any],
    name: str,
    default: float,
    *,
    minimum: float = 0.0,
    maximum: float = 1.0,
) -> float:
    value = parameters.get(name, default)
    if isinstance(value, bool) or not isinstance(value, (int, float)):
        raise KinshipPanelInputError(f"invalid_parameter:{name}")
    result = float(value)
    if not math.isfinite(result) or not minimum <= result <= maximum:
        raise KinshipPanelInputError(f"invalid_parameter:{name}")
    return result


def _excluded_regions(parameters: dict[str, Any]) -> list[dict[str, int | str]]:
    raw_regions = parameters.get("excluded_regions", [])
    if not isinstance(raw_regions, list):
        raise KinshipPanelInputError("invalid_parameter:excluded_regions")
    regions: list[dict[str, int | str]] = []
    seen_ids: set[str] = set()
    for region in raw_regions:
        if not isinstance(region, dict) or set(region) != {
            "region_id",
            "chromosome",
            "start_bp",
            "end_bp",
        }:
            raise KinshipPanelInputError("invalid_parameter:excluded_regions")
        region_id = region["region_id"]
        chromosome = region["chromosome"]
        start_bp = region["start_bp"]
        end_bp = region["end_bp"]
        if (
            not isinstance(region_id, str)
            or not region_id
            or region_id in seen_ids
            or isinstance(chromosome, bool)
            or not isinstance(chromosome, int)
            or not 1 <= chromosome <= 22
            or isinstance(start_bp, bool)
            or not isinstance(start_bp, int)
            or isinstance(end_bp, bool)
            or not isinstance(end_bp, int)
            or start_bp <= 0
            or end_bp < start_bp
        ):
            raise KinshipPanelInputError("invalid_parameter:excluded_regions")
        seen_ids.add(region_id)
        regions.append(
            {
                "region_id": region_id,
                "chromosome": chromosome,
                "start_bp": start_bp,
                "end_bp": end_bp,
            }
        )
    return regions


def _parameters(parameters: dict[str, Any]) -> dict[str, Any]:
    required_autosomes = parameters.get("required_autosomes", list(range(1, 23)))
    if (
        not isinstance(required_autosomes, list)
        or not required_autosomes
        or any(
            isinstance(chromosome, bool)
            or not isinstance(chromosome, int)
            or not 1 <= chromosome <= 22
            for chromosome in required_autosomes
        )
        or len(required_autosomes) != len(set(required_autosomes))
    ):
        raise KinshipPanelInputError("invalid_parameter:required_autosomes")
    resolved = {
        "maf_min": _bounded_float(parameters, "maf_min", 0.05, maximum=0.5),
        "ld_window_variants": _positive_integer(parameters, "ld_window_variants", 50),
        "ld_step_variants": _positive_integer(parameters, "ld_step_variants", 5),
        "ld_r2_max": _bounded_float(
            parameters, "ld_r2_max", 0.2, minimum=1e-12
        ),
        "required_autosomes": sorted(required_autosomes),
        "min_markers_per_autosome": _positive_integer(
            parameters, "min_markers_per_autosome", 1
        ),
        "min_total_markers": _positive_integer(parameters, "min_total_markers", 1000),
        "max_autosome_marker_fraction": _bounded_float(
            parameters,
            "max_autosome_marker_fraction",
            0.15,
            minimum=1e-12,
        ),
        "excluded_regions": _excluded_regions(parameters),
        "residual_ld_window_variants": _positive_integer(
            parameters, "residual_ld_window_variants", 100
        ),
        "residual_ld_window_kb": _positive_integer(
            parameters, "residual_ld_window_kb", 1000
        ),
        "plink_timeout_seconds": _positive_integer(
            parameters, "plink_timeout_seconds", 300
        ),
    }
    if resolved["ld_step_variants"] > resolved["ld_window_variants"]:
        raise KinshipPanelInputError("ld_step_exceeds_window")
    return resolved


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
    output = completed.stdout.strip() or completed.stderr.strip()
    return output.splitlines()[0][:500] if output else None


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
        raise ExternalToolError("plink_execution_failed") from error
    if completed.returncode != 0:
        raise ExternalToolError("plink_command_failed")


def _validate_input_dataset(
    descriptor: dict[str, Any],
    paths: dict[str, Path],
    artifacts: dict[str, dict[str, Any]],
) -> tuple[list[list[str]], list[list[str]]]:
    validate_json_document(descriptor, "plink_dataset.schema.json")
    if descriptor["dataset_id"] != "genomewide_pre_qc":
        raise KinshipPanelInputError("unexpected_pre_qc_dataset_id")
    if descriptor["scope"] != "autosomal_genomewide":
        raise KinshipPanelInputError("kinship_panel_requires_autosomal_dataset")
    if paths["genomewide_pre_qc_bed"].read_bytes()[:3] != b"\x6c\x1b\x01":
        raise KinshipPanelInputError("invalid_pre_qc_bed_header")
    for suffix in ("bed", "bim", "fam"):
        artifact_id = f"genomewide_pre_qc_{suffix}"
        physical_path = paths[artifact_id]
        declared_file = descriptor["files"][suffix]
        if (
            declared_file["name"] != physical_path.name
            or declared_file["sha256"] != sha256_file(physical_path)
            or artifacts[artifact_id]["sample_set_id"] != descriptor["sample_set_id"]
            or artifacts[artifact_id]["variant_set_id"] != descriptor["variant_set_id"]
        ):
            raise KinshipPanelInputError("pre_qc_dataset_descriptor_mismatch")
    fam_rows = _read_fam(paths["genomewide_pre_qc_fam"])
    bim_rows = _read_bim(paths["genomewide_pre_qc_bim"])
    if descriptor["sample_count"] != len(fam_rows) or descriptor["variant_count"] != len(
        bim_rows
    ):
        raise KinshipPanelInputError("pre_qc_dataset_count_mismatch")
    return fam_rows, bim_rows


def _preliminary_alerts_by_variant(
    metrics_path: Path, bim_rows: list[list[str]]
) -> dict[str, str]:
    table = validate_tsv_table(metrics_path, "qc_variant_metrics.schema.json")
    rows_by_id = {row["VARIANT_ID"]: row for row in table.rows}
    input_ids = {row[1] for row in bim_rows}
    if not input_ids.issubset(rows_by_id):
        raise KinshipPanelInputError("preliminary_qc_variant_set_incomplete")
    if any(rows_by_id[variant_id]["QC_STATUS"] == "EXCLUDED" for variant_id in input_ids):
        raise KinshipPanelInputError("excluded_variant_present_in_pre_qc_dataset")
    return {
        variant_id: rows_by_id[variant_id]["ALERT_CODES"] or ""
        for variant_id in input_ids
    }


def _frequency_by_variant(
    rows: list[dict[str, str]], bim_rows: list[list[str]]
) -> dict[str, float | None]:
    frequency_rows = {row["SNP"]: row for row in rows}
    expected_ids = {row[1] for row in bim_rows}
    if set(frequency_rows) != expected_ids:
        raise ExternalToolError("kinship_frequency_variant_set_mismatch")
    frequencies: dict[str, float | None] = {}
    for variant_id, row in frequency_rows.items():
        try:
            maf = float(row["MAF"])
        except (KeyError, ValueError):
            maf = None
        if maf is not None and (not math.isfinite(maf) or not 0 <= maf <= 0.5):
            raise ExternalToolError("invalid_kinship_maf")
        frequencies[variant_id] = maf
    return frequencies


def _region_id_for_variant(
    chromosome: int,
    position_bp: int,
    regions: list[dict[str, int | str]],
) -> str | None:
    matching_ids = [
        str(region["region_id"])
        for region in regions
        if region["chromosome"] == chromosome
        and int(region["start_bp"]) <= position_bp <= int(region["end_bp"])
    ]
    return ",".join(matching_ids) if matching_ids else None


def _eligible_variants(
    bim_rows: list[list[str]],
    frequencies: dict[str, float | None],
    parameters: dict[str, Any],
) -> tuple[list[str], dict[str, str | None]]:
    eligible: list[str] = []
    exclusions: dict[str, str | None] = {}
    for row in bim_rows:
        chromosome = int(row[0])
        variant_id = row[1]
        position_bp = int(row[3])
        maf = frequencies[variant_id]
        region_ids = _region_id_for_variant(
            chromosome, position_bp, parameters["excluded_regions"]
        )
        if maf is None:
            exclusions[variant_id] = "maf_unavailable"
        elif maf < parameters["maf_min"]:
            exclusions[variant_id] = "low_maf"
        elif region_ids is not None:
            exclusions[variant_id] = f"complex_region:{region_ids}"
        else:
            exclusions[variant_id] = None
            eligible.append(variant_id)
    if not eligible:
        raise KinshipPanelBlockError("no_variant_eligible_before_ld_pruning")
    return eligible, exclusions


def _write_lines(path: Path, values: list[str]) -> None:
    path.write_text("".join(f"{value}\n" for value in values), encoding="utf-8")


def _validate_pruning_partition(
    eligible_ids: list[str], retained_ids: list[str], ld_excluded_ids: list[str]
) -> None:
    eligible_set = set(eligible_ids)
    retained_set = set(retained_ids)
    excluded_set = set(ld_excluded_ids)
    if retained_set & excluded_set or retained_set | excluded_set != eligible_set:
        raise ExternalToolError("invalid_ld_pruning_partition")


def _validate_panel_dataset(
    panel_prefix: Path,
    input_fam_rows: list[list[str]],
    expected_variant_ids: list[str],
) -> tuple[list[list[str]], list[list[str]]]:
    if panel_prefix.with_suffix(".bed").read_bytes()[:3] != b"\x6c\x1b\x01":
        raise ExternalToolError("invalid_kinship_panel_bed_header")
    fam_rows = _read_fam(panel_prefix.with_suffix(".fam"))
    bim_rows = _read_bim(panel_prefix.with_suffix(".bim"))
    if fam_rows != input_fam_rows:
        raise ExternalToolError("kinship_panel_sample_set_changed")
    if [row[1] for row in bim_rows] != expected_variant_ids:
        raise ExternalToolError("unexpected_kinship_panel_variant_set")
    return fam_rows, bim_rows


def _decimal(value: float | None) -> str:
    if value is None:
        return ""
    normalized = 0.0 if abs(value) < 5e-13 else value
    return f"{normalized:.8f}".rstrip("0").rstrip(".")


def _write_tsv(
    path: Path, columns: tuple[str, ...], rows: list[dict[str, str]]
) -> None:
    with path.open("w", encoding="utf-8", newline="") as output_file:
        writer = csv.DictWriter(
            output_file, fieldnames=columns, delimiter="\t", lineterminator="\n"
        )
        writer.writeheader()
        writer.writerows({column: row[column] for column in columns} for row in rows)


def _variant_audit_rows(
    bim_rows: list[list[str]],
    frequencies: dict[str, float | None],
    preliminary_alerts: dict[str, str],
    exclusions: dict[str, str | None],
    retained_ids: set[str],
) -> list[dict[str, str]]:
    rows: list[dict[str, str]] = []
    for bim_row in bim_rows:
        variant_id = bim_row[1]
        exclusion = exclusions[variant_id]
        if variant_id in retained_ids:
            status = "RETAINED"
            exclusion_code = ""
        elif exclusion == "low_maf":
            status = "EXCLUDED_LOW_MAF"
            exclusion_code = exclusion
        elif exclusion == "maf_unavailable":
            status = "EXCLUDED_MAF_UNAVAILABLE"
            exclusion_code = exclusion
        elif exclusion is not None and exclusion.startswith("complex_region:"):
            status = "EXCLUDED_COMPLEX_REGION"
            exclusion_code = exclusion
        else:
            status = "EXCLUDED_LD"
            exclusion_code = "ld_pruning"
        rows.append(
            {
                "VARIANT_ID": variant_id,
                "CHROMOSOME": bim_row[0],
                "POSITION_BP": bim_row[3],
                "MAF": _decimal(frequencies[variant_id]),
                "PRELIMINARY_ALERT_CODES": preliminary_alerts[variant_id],
                "IN_COMPLEX_REGION": str(
                    exclusion is not None and exclusion.startswith("complex_region:")
                ).lower(),
                "PANEL_STATUS": status,
                "EXCLUSION_CODE": exclusion_code,
            }
        )
    return rows


def _coverage_rows(
    input_bim_rows: list[list[str]],
    eligible_ids: set[str],
    panel_bim_rows: list[list[str]],
    parameters: dict[str, Any],
) -> tuple[list[dict[str, str]], list[str]]:
    input_counts = Counter(int(row[0]) for row in input_bim_rows)
    eligible_counts = Counter(
        int(row[0]) for row in input_bim_rows if row[1] in eligible_ids
    )
    positions_by_chromosome: dict[int, list[int]] = defaultdict(list)
    for row in panel_bim_rows:
        positions_by_chromosome[int(row[0])].append(int(row[3]))
    total_markers = len(panel_bim_rows)
    required = set(parameters["required_autosomes"])
    rows: list[dict[str, str]] = []
    blocking_reasons: list[str] = []
    for chromosome in range(1, 23):
        positions = sorted(positions_by_chromosome.get(chromosome, []))
        gaps = [right - left for left, right in zip(positions, positions[1:])]
        fraction = len(positions) / total_markers if total_markers else 0.0
        reasons: list[str] = []
        if chromosome in required and len(positions) < parameters["min_markers_per_autosome"]:
            reasons.append("insufficient_markers")
        if fraction > parameters["max_autosome_marker_fraction"]:
            reasons.append("excessive_marker_concentration")
        if reasons:
            blocking_reasons.extend(f"chr{chromosome}:{reason}" for reason in reasons)
        rows.append(
            {
                "CHROMOSOME": str(chromosome),
                "INPUT_VARIANT_COUNT": str(input_counts[chromosome]),
                "MAF_ELIGIBLE_COUNT": str(eligible_counts[chromosome]),
                "PRUNED_MARKER_COUNT": str(len(positions)),
                "MARKER_FRACTION": _decimal(fraction),
                "SPAN_BP": str(positions[-1] - positions[0]) if positions else "",
                "MEDIAN_GAP_BP": _decimal(statistics.median(gaps)) if gaps else "",
                "MAX_GAP_BP": str(max(gaps)) if gaps else "",
                "COVERAGE_STATUS": "FAIL" if reasons else "PASS",
                "REASON_CODE": ",".join(reasons),
            }
        )
    if total_markers < parameters["min_total_markers"]:
        blocking_reasons.append("insufficient_total_markers")
    return rows, blocking_reasons


def _ld_distribution(
    rows: list[dict[str, str]], panel_variant_ids: set[str]
) -> tuple[list[dict[str, str]], dict[str, Any]]:
    values: list[float] = []
    for row in rows:
        if row.get("SNP_A") not in panel_variant_ids or row.get("SNP_B") not in panel_variant_ids:
            raise ExternalToolError("residual_ld_variant_outside_panel")
        try:
            value = float(row["R2"])
        except (KeyError, ValueError) as error:
            raise ExternalToolError("invalid_residual_ld_report") from error
        if not math.isfinite(value) or not 0 <= value <= 1:
            raise ExternalToolError("invalid_residual_ld_value")
        values.append(value)
    boundaries = [index / 10 for index in range(11)]
    counts = [0] * 10
    for value in values:
        index = min(int(value * 10), 9)
        counts[index] += 1
    bins = [
        {
            "BIN_ID": f"r2_bin_{index + 1:02d}",
            "R2_LOWER_BOUND": _decimal(boundaries[index]),
            "R2_UPPER_BOUND": _decimal(boundaries[index + 1]),
            "UPPER_INCLUSIVE": str(index == 9).lower(),
            "PAIR_COUNT": str(counts[index]),
        }
        for index in range(10)
    ]
    summary = {
        "pair_count": len(values),
        "mean_r2": statistics.mean(values) if values else None,
        "median_r2": statistics.median(values) if values else None,
        "max_r2": max(values) if values else None,
    }
    return bins, summary


def _set_id(scope: str, values: list[str]) -> str:
    digest = hashlib.sha256(
        "".join(f"{value}\n" for value in values).encode("utf-8")
    ).hexdigest()
    return f"{scope}_{digest[:12]}"


def _write_dataset_descriptor(
    path: Path,
    panel_prefix: Path,
    input_descriptor: dict[str, Any],
    variant_set_id: str,
    sample_count: int,
    variant_count: int,
) -> dict[str, Any]:
    descriptor = {
        **input_descriptor,
        "dataset_id": "kinship_panel",
        "variant_set_id": variant_set_id,
        "sample_count": sample_count,
        "variant_count": variant_count,
        "source_format": "PLINK_KINSHIP_PANEL",
        "files": {
            suffix: {
                "name": panel_prefix.with_suffix(f".{suffix}").name,
                "sha256": sha256_file(panel_prefix.with_suffix(f".{suffix}")),
            }
            for suffix in ("bed", "bim", "fam")
        },
    }
    validate_json_document(descriptor, "plink_dataset.schema.json")
    atomic_write_json(path, descriptor)
    return descriptor


def _output_artifacts(
    stage_inputs: dict[str, Any],
    descriptor: dict[str, Any],
    input_variant_set_id: str,
    paths_and_contracts: list[tuple[str, str, Path, str | None, str]],
) -> list[dict[str, Any]]:
    producer_stage = f"{stage_inputs['stage_id']}_{stage_inputs['stage_name']}"
    published_dir = PurePosixPath(stage_inputs["published_output_dir"])
    artifacts: list[dict[str, Any]] = []
    for artifact_id, artifact_type, path, schema_name, sensitivity in paths_and_contracts:
        artifacts.append(
            build_file_artifact(
                physical_path=path,
                published_path=str(published_dir / path.name),
                artifact_id=artifact_id,
                artifact_type=artifact_type,
                media_type=(
                    "application/octet-stream"
                    if path.suffix == ".bed"
                    else "application/json"
                    if path.suffix == ".json"
                    else "text/tab-separated-values"
                    if path.suffix == ".tsv"
                    else "text/plain"
                ),
                producer_stage=producer_stage,
                producer_signature=stage_inputs["signature"],
                schema_name=schema_name,
                schema_version="1.0.0" if schema_name else None,
                assembly=descriptor["assembly"],
                sample_set_id=descriptor["sample_set_id"],
                variant_set_id=(
                    input_variant_set_id
                    if artifact_id == "pruned_variants"
                    else descriptor["variant_set_id"]
                ),
                sensitivity=sensitivity,
            )
        )
    return artifacts


def execute(stage_inputs_path: Path, output_dir: Path) -> int:
    """Construit et audite le panel de marqueurs destiné à KING et à la PCA."""
    started_at = utc_now()
    started_clock = monotonic()
    stage_inputs = read_json(stage_inputs_path)
    validate_json_document(stage_inputs, "stage_inputs.schema.json")
    run_dir = output_dir.parent.parent
    config = load_pipeline_config(run_dir / "config.resolved.yaml")
    parameters = _parameters(stage_inputs["parameters"])
    timeout_seconds = int(parameters["plink_timeout_seconds"])
    artifact_ids = (
        "genomewide_pre_qc_bed",
        "genomewide_pre_qc_bim",
        "genomewide_pre_qc_fam",
        "genomewide_pre_qc_dataset",
        "qc_variant_metrics",
    )
    input_artifacts = {
        artifact_id: _artifact_by_id(stage_inputs, artifact_id)
        for artifact_id in artifact_ids
    }
    paths = {
        artifact_id: _validated_artifact_path(artifact, run_dir)
        for artifact_id, artifact in input_artifacts.items()
    }
    input_descriptor = read_json(paths["genomewide_pre_qc_dataset"])
    input_fam_rows, input_bim_rows = _validate_input_dataset(
        input_descriptor, paths, input_artifacts
    )
    preliminary_alerts = _preliminary_alerts_by_variant(
        paths["qc_variant_metrics"], input_bim_rows
    )

    output_dir.mkdir(parents=True, exist_ok=True)
    plink_executable = _resolve_executable(config["tools"]["plink"])
    plink_version = _plink_version(plink_executable)
    input_prefix = paths["genomewide_pre_qc_bed"].with_suffix("")
    frequency_prefix = output_dir / "kinship_panel_frequency"
    _run_plink(
        plink_executable,
        ["--bfile", str(input_prefix), "--freq", "--out", str(frequency_prefix)],
        timeout_seconds,
    )
    frequencies = _frequency_by_variant(
        _read_whitespace_table(frequency_prefix.with_suffix(".frq")), input_bim_rows
    )
    eligible_ids, initial_exclusions = _eligible_variants(
        input_bim_rows, frequencies, parameters
    )
    eligible_path = output_dir / "kinship_panel_eligible_variants.txt"
    _write_lines(eligible_path, eligible_ids)

    pruning_prefix = output_dir / "kinship_panel_pruning"
    _run_plink(
        plink_executable,
        [
            "--bfile",
            str(input_prefix),
            "--extract",
            str(eligible_path),
            "--indep-pairwise",
            str(parameters["ld_window_variants"]),
            str(parameters["ld_step_variants"]),
            str(parameters["ld_r2_max"]),
            "--out",
            str(pruning_prefix),
        ],
        timeout_seconds,
    )
    retained_ids = _read_variant_list(pruning_prefix.with_suffix(".prune.in"))
    ld_excluded_ids = _read_variant_list(
        pruning_prefix.with_suffix(".prune.out"), allow_empty=True
    )
    _validate_pruning_partition(eligible_ids, retained_ids, ld_excluded_ids)

    retained_path = output_dir / "kinship_panel_retained_variants.txt"
    _write_lines(retained_path, retained_ids)
    panel_prefix = output_dir / "kinship_panel"
    _run_plink(
        plink_executable,
        [
            "--bfile",
            str(input_prefix),
            "--extract",
            str(retained_path),
            "--make-bed",
            "--out",
            str(panel_prefix),
        ],
        timeout_seconds,
    )
    panel_fam_rows, panel_bim_rows = _validate_panel_dataset(
        panel_prefix, input_fam_rows, retained_ids
    )

    pruned_path = output_dir / "pruned_variants.tsv"
    coverage_path = output_dir / "panel_coverage.tsv"
    pruned_rows = _variant_audit_rows(
        input_bim_rows,
        frequencies,
        preliminary_alerts,
        initial_exclusions,
        set(retained_ids),
    )
    coverage_rows, blocking_reasons = _coverage_rows(
        input_bim_rows, set(eligible_ids), panel_bim_rows, parameters
    )
    _write_tsv(pruned_path, PRUNED_VARIANT_COLUMNS, pruned_rows)
    _write_tsv(coverage_path, COVERAGE_COLUMNS, coverage_rows)
    validate_tsv_table(pruned_path, "kinship_pruned_variants.schema.json")
    validate_tsv_table(coverage_path, "kinship_panel_coverage.schema.json")
    if blocking_reasons:
        raise KinshipPanelBlockError(
            "insufficient_or_unbalanced_panel:" + ",".join(blocking_reasons)
        )

    residual_prefix = output_dir / "kinship_panel_residual_ld"
    _run_plink(
        plink_executable,
        [
            "--bfile",
            str(panel_prefix),
            "--r2",
            "--ld-window",
            str(parameters["residual_ld_window_variants"]),
            "--ld-window-kb",
            str(parameters["residual_ld_window_kb"]),
            "--ld-window-r2",
            "0",
            "--out",
            str(residual_prefix),
        ],
        timeout_seconds,
    )
    ld_rows = _read_whitespace_table(residual_prefix.with_suffix(".ld"), allow_empty=True)
    ld_bins, ld_summary = _ld_distribution(ld_rows, set(retained_ids))
    ld_bins_path = output_dir / "panel_residual_ld_bins.tsv"
    _write_tsv(ld_bins_path, LD_BIN_COLUMNS, ld_bins)
    validate_tsv_table(ld_bins_path, "kinship_panel_ld_bins.schema.json")

    variant_set_id = _set_id("kinship_variants", retained_ids)
    descriptor_path = output_dir / "kinship_panel.dataset.json"
    output_descriptor = _write_dataset_descriptor(
        descriptor_path,
        panel_prefix,
        input_descriptor,
        variant_set_id,
        len(panel_fam_rows),
        len(panel_bim_rows),
    )
    report_path = output_dir / "kinship_panel_report.json"
    status_counts = Counter(row["PANEL_STATUS"] for row in pruned_rows)
    report = {
        "schema_version": "1.0.0",
        "assembly": output_descriptor["assembly"],
        "parameters": parameters,
        "sample_count": len(panel_fam_rows),
        "input_variant_count": len(input_bim_rows),
        "maf_eligible_variant_count": sum(
            frequencies[row[1]] is not None
            and float(frequencies[row[1]]) >= parameters["maf_min"]
            for row in input_bim_rows
        ),
        "ld_eligible_variant_count": len(eligible_ids),
        "retained_marker_count": len(panel_bim_rows),
        "effective_marker_count_operational": len(panel_bim_rows),
        "status_counts": dict(sorted(status_counts.items())),
        "covered_autosome_count": sum(
            int(row["PRUNED_MARKER_COUNT"]) > 0 for row in coverage_rows
        ),
        "required_autosome_count": len(parameters["required_autosomes"]),
        "residual_ld": ld_summary,
    }
    atomic_write_json(report_path, report)

    temporary_paths = [
        eligible_path,
        retained_path,
        frequency_prefix.with_suffix(".frq"),
        pruning_prefix.with_suffix(".prune.in"),
        pruning_prefix.with_suffix(".prune.out"),
        residual_prefix.with_suffix(".ld"),
    ]
    for prefix in (frequency_prefix, pruning_prefix, panel_prefix, residual_prefix):
        temporary_paths.extend(
            prefix.with_suffix(f".{suffix}") for suffix in ("log", "nosex")
        )
    for path in temporary_paths:
        if path.exists():
            path.unlink()

    paths_and_contracts = [
        (
            f"kinship_panel_{suffix}",
            f"plink_{suffix}",
            panel_prefix.with_suffix(f".{suffix}"),
            None,
            "sensitive_genetic",
        )
        for suffix in ("bed", "bim", "fam")
    ]
    paths_and_contracts.extend(
        [
            (
                "kinship_panel_dataset",
                "plink_dataset_descriptor",
                descriptor_path,
                "plink_dataset.schema.json",
                "sensitive_genetic",
            ),
            (
                "pruned_variants",
                "kinship_pruned_variants",
                pruned_path,
                "kinship_pruned_variants.schema.json",
                "sensitive_genetic",
            ),
            (
                "panel_coverage",
                "kinship_panel_coverage",
                coverage_path,
                "kinship_panel_coverage.schema.json",
                "internal",
            ),
            (
                "panel_residual_ld_bins",
                "kinship_panel_ld_bins",
                ld_bins_path,
                "kinship_panel_ld_bins.schema.json",
                "internal",
            ),
            (
                "kinship_panel_report",
                "kinship_panel_report",
                report_path,
                None,
                "internal",
            ),
        ]
    )
    output_artifacts = _output_artifacts(
        stage_inputs,
        output_descriptor,
        input_descriptor["variant_set_id"],
        paths_and_contracts,
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
        "method_id": "plink_maf_ld_pruned_kinship_panel_v1",
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
            "samples": len(panel_fam_rows),
            "variants_before": len(input_bim_rows),
            "variants_after": len(panel_bim_rows),
            "autosomes_covered": report["covered_autosome_count"],
            "residual_ld_pairs": ld_summary["pair_count"],
        },
        "metrics": {
            "effective_marker_count_operational": len(panel_bim_rows),
            "maximum_autosome_marker_fraction": max(
                float(row["MARKER_FRACTION"]) for row in coverage_rows
            ),
            "residual_ld_mean_r2": ld_summary["mean_r2"],
            "residual_ld_max_r2": ld_summary["max_r2"],
        },
        "exclusions": [
            {"scope": "variant", "count": count, "reason": status}
            for status, count in sorted(status_counts.items())
            if status != "RETAINED"
        ],
        "warnings": (
            []
            if ld_summary["pair_count"] > 0
            else [
                {
                    "code": "residual_ld_no_pairs",
                    "message": "Aucune paire n'est évaluable dans la fenêtre de LD résiduel.",
                }
            ]
        ),
        "checks": [
            {"check": "pre_qc_input_integrity", "status": "PASS"},
            {"check": "maf_recalculated_after_sample_qc", "status": "PASS"},
            {"check": "ld_pruning_partition", "status": "PASS"},
            {"check": "required_autosome_coverage", "status": "PASS"},
            {"check": "marker_concentration", "status": "PASS"},
            {
                "check": "residual_ld_distribution",
                "status": "PASS" if ld_summary["pair_count"] > 0 else "NOT_EVALUATED",
            },
        ],
        "known_limits": [
            "Le nombre effectif publié est le nombre opérationnel de marqueurs "
            "après pruning, pas une estimation spectrale de dimension effective.",
            "Le LD résiduel publié est agrégé par classes de R2 et ne conserve "
            "pas les paires individuelles.",
        ],
        "expected_visualizations": [
            "markers_per_autosome",
            "intermarker_distances",
            "residual_ld_distribution",
        ],
        "manual_validation_required": False,
    }
    validate_json_document(audit, "stage_audit.schema.json")
    atomic_write_json(output_dir / "audit.json", audit)
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
    except KinshipPanelBlockError as error:
        sys.stderr.write(f"{error}\n")
        return 4
    except (
        KinshipPanelInputError,
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
