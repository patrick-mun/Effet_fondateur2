"""Étape 10 : QC final spécifique à chaque cohorte analytique."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import shutil
import subprocess
import sys
import tempfile
from collections import defaultdict
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
    validate_cohorts_frozen,
    validate_json_document,
    validate_tsv_table,
)
from effet_fondateur.orchestrator.state import utc_now


SAMPLE_COLUMNS = (
    "DATASET_ID",
    "COHORT_ID",
    "SAMPLE_ID",
    "FID",
    "IID",
    "ARRAY_BATCH",
    "N_MISS",
    "N_GENO",
    "F_MISS",
    "QC_STATUS",
    "EXCLUSION_CODE",
    "PRESENT_AFTER_QC",
)
VARIANT_COLUMNS = (
    "DATASET_ID",
    "COHORT_ID",
    "VARIANT_ID",
    "CHROMOSOME",
    "POSITION_BP",
    "A1",
    "A2",
    "N_MISS",
    "N_GENO",
    "F_MISS",
    "MAF",
    "HWE_P",
    "HWE_STATUS",
    "QC_STATUS",
    "FILTER_CODES",
    "IS_TARGET_VARIANT",
    "FILTER_EXCEPTION",
    "PRESENT_AFTER_QC",
)
BATCH_COLUMNS = (
    "DATASET_ID",
    "COHORT_ID",
    "ARRAY_BATCH",
    "SAMPLE_COUNT",
    "MEAN_F_MISS",
    "BATCH_DELTA_FROM_BEST",
    "QC_STATUS",
    "ALERT_CODE",
)
DATASET_POLICIES = (
    {
        "dataset_id": "controls_unrelated_qc",
        "cohort_id": "controls_unrelated",
        "source": "genomewide_pre_qc",
        "apply_hwe": True,
        "maf_parameter": "controls_maf_min",
        "variant_missingness_parameter": "controls_variant_missingness_max",
        "force_target": False,
    },
    {
        "dataset_id": "target_carriers_independent_qc",
        "cohort_id": "target_carriers_independent",
        "source": "genomewide_pre_qc",
        "apply_hwe": False,
        "maf_parameter": None,
        "variant_missingness_parameter": "carriers_variant_missingness_max",
        "force_target": False,
    },
    {
        "dataset_id": "target_chromosome_all_qc",
        "cohort_id": "target_chromosome_all_qc",
        "source": "target_variant",
        "apply_hwe": False,
        "maf_parameter": "target_region_maf_min",
        "variant_missingness_parameter": "target_region_variant_missingness_max",
        "force_target": True,
    },
)


class FinalQcInputError(ValueError):
    """Signale une entrée ou une configuration de QC final invalide."""


class ExternalToolError(RuntimeError):
    """Signale un échec PLINK ou une sortie native incohérente."""


class FinalQcBlockError(RuntimeError):
    """Signale qu'un jeu final scientifiquement requis serait vide."""


def _build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(prog="qc-final")
    parser.add_argument("--stage-inputs", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    return parser


def _bounded_float(
    parameters: dict[str, Any],
    name: str,
    default: float,
    *,
    positive: bool = False,
) -> float:
    value = parameters.get(name, default)
    if isinstance(value, bool) or not isinstance(value, (int, float)):
        raise FinalQcInputError(f"invalid_parameter:{name}")
    result = float(value)
    lower_valid = result > 0 if positive else result >= 0
    if not math.isfinite(result) or not lower_valid or result > 1:
        raise FinalQcInputError(f"invalid_parameter:{name}")
    return result


def _positive_integer(parameters: dict[str, Any], name: str, default: int) -> int:
    value = parameters.get(name, default)
    if isinstance(value, bool) or not isinstance(value, int) or value <= 0:
        raise FinalQcInputError(f"invalid_parameter:{name}")
    return value


def _parameters(parameters: dict[str, Any]) -> dict[str, Any]:
    return {
        "sample_missingness_max": _bounded_float(
            parameters, "sample_missingness_max", 0.1
        ),
        "controls_variant_missingness_max": _bounded_float(
            parameters, "controls_variant_missingness_max", 0.05
        ),
        "controls_maf_min": _bounded_float(parameters, "controls_maf_min", 0.05),
        "controls_hwe_p_min": _bounded_float(
            parameters, "controls_hwe_p_min", 1e-6, positive=True
        ),
        "carriers_variant_missingness_max": _bounded_float(
            parameters, "carriers_variant_missingness_max", 0.05
        ),
        "target_region_variant_missingness_max": _bounded_float(
            parameters, "target_region_variant_missingness_max", 0.05
        ),
        "target_region_maf_min": _bounded_float(
            parameters, "target_region_maf_min", 0.01
        ),
        "batch_missingness_delta_alert": _bounded_float(
            parameters, "batch_missingness_delta_alert", 0.02
        ),
        "plink_timeout_seconds": _positive_integer(
            parameters, "plink_timeout_seconds", 300
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
        raise FinalQcInputError(f"{artifact_id}_missing_or_ambiguous")
    return matches[0]


def _validated_artifact_path(artifact: dict[str, Any], run_dir: Path) -> Path:
    path = _resolve_input_path(artifact, run_dir)
    if not path.is_file() or sha256_file(path) != artifact["sha256"]:
        raise FinalQcInputError("declared_input_modified")
    return path


def _load_target_metadata(path: Path) -> dict[str, Any]:
    try:
        metadata = yaml.safe_load(path.read_text(encoding="utf-8"))
    except (OSError, yaml.YAMLError) as error:
        raise FinalQcInputError("invalid_target_variant_metadata") from error
    if not isinstance(metadata, dict):
        raise FinalQcInputError("invalid_target_variant_metadata")
    validate_json_document(metadata, "target_variant_metadata.schema.json")
    return metadata


def _read_fam(path: Path) -> list[list[str]]:
    rows = [line.split() for line in path.read_text(encoding="utf-8").splitlines() if line]
    if not rows or any(len(row) != 6 for row in rows):
        raise FinalQcInputError(f"malformed_fam:{path.name}")
    if len({(row[0], row[1]) for row in rows}) != len(rows):
        raise FinalQcInputError(f"duplicate_sample_in_fam:{path.name}")
    return rows


def _read_bim(path: Path) -> list[list[str]]:
    rows = [line.split() for line in path.read_text(encoding="utf-8").splitlines() if line]
    if not rows or any(len(row) != 6 for row in rows):
        raise FinalQcInputError(f"malformed_bim:{path.name}")
    if len({row[1] for row in rows}) != len(rows):
        raise FinalQcInputError(f"duplicate_variant_in_bim:{path.name}")
    return rows


def _validate_dataset(
    dataset_prefix: str,
    paths: dict[str, Path],
    artifacts: dict[str, dict[str, Any]],
) -> tuple[dict[str, Any], list[list[str]], list[list[str]]]:
    descriptor = read_json(paths[f"{dataset_prefix}_dataset"])
    validate_json_document(descriptor, "plink_dataset.schema.json")
    if paths[f"{dataset_prefix}_bed"].read_bytes()[:3] != b"\x6c\x1b\x01":
        raise FinalQcInputError(f"invalid_bed_header:{dataset_prefix}")
    stems = {
        paths[f"{dataset_prefix}_{suffix}"].with_suffix("")
        for suffix in ("bed", "bim", "fam")
    }
    if len(stems) != 1:
        raise FinalQcInputError(f"dataset_files_do_not_share_prefix:{dataset_prefix}")
    for suffix in ("bed", "bim", "fam"):
        artifact_id = f"{dataset_prefix}_{suffix}"
        path = paths[artifact_id]
        if (
            descriptor["files"][suffix]["name"] != path.name
            or descriptor["files"][suffix]["sha256"] != artifacts[artifact_id]["sha256"]
            or artifacts[artifact_id]["sample_set_id"] != descriptor["sample_set_id"]
            or artifacts[artifact_id]["variant_set_id"] != descriptor["variant_set_id"]
        ):
            raise FinalQcInputError(f"dataset_descriptor_mismatch:{dataset_prefix}")
    fam_rows = _read_fam(paths[f"{dataset_prefix}_fam"])
    bim_rows = _read_bim(paths[f"{dataset_prefix}_bim"])
    if descriptor["sample_count"] != len(fam_rows) or descriptor["variant_count"] != len(bim_rows):
        raise FinalQcInputError(f"dataset_count_mismatch:{dataset_prefix}")
    return descriptor, fam_rows, bim_rows


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
            [executable, "--version"], capture_output=True, check=False, text=True, timeout=10
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
        raise ExternalToolError("plink_final_qc_failed") from error
    if completed.returncode != 0:
        raise ExternalToolError(f"plink_final_qc_failed:{completed.returncode}")


def _read_whitespace_table(path: Path) -> list[dict[str, str]]:
    try:
        lines = [line.split() for line in path.read_text(encoding="utf-8").splitlines() if line]
    except OSError as error:
        raise ExternalToolError(f"missing_plink_report:{path.name}") from error
    if len(lines) < 2:
        raise ExternalToolError(f"empty_plink_report:{path.name}")
    header = lines[0]
    if any(len(row) != len(header) for row in lines[1:]):
        raise ExternalToolError(f"malformed_plink_report:{path.name}")
    return [dict(zip(header, row, strict=True)) for row in lines[1:]]


def _float_value(raw_value: str, field: str, *, nullable: bool = False) -> float | None:
    if nullable and raw_value in {"NA", "nan", "NaN"}:
        return None
    try:
        value = float(raw_value)
    except ValueError as error:
        raise ExternalToolError(f"invalid_plink_field:{field}") from error
    if not math.isfinite(value):
        raise ExternalToolError(f"invalid_plink_field:{field}")
    return value


def _integer_value(raw_value: str, field: str, *, minimum: int = 0) -> int:
    try:
        value = int(raw_value)
    except ValueError as error:
        raise ExternalToolError(f"invalid_plink_field:{field}") from error
    if value < minimum:
        raise ExternalToolError(f"invalid_plink_field:{field}")
    return value


def _decimal(value: float) -> str:
    if abs(value) < 5e-16:
        value = 0.0
    return format(value, ".12g")


def _keep_pairs(path: Path) -> list[tuple[str, str]]:
    rows = [line.split() for line in path.read_text(encoding="utf-8").splitlines() if line]
    if not rows or any(len(row) != 2 for row in rows):
        raise FinalQcBlockError(f"empty_or_malformed_keep:{path.name}")
    pairs = [(row[0], row[1]) for row in rows]
    if len(pairs) != len(set(pairs)):
        raise FinalQcInputError(f"duplicate_keep_sample:{path.name}")
    return pairs


def _sample_context(
    paths: dict[str, Path],
    source_fams: dict[str, list[list[str]]],
) -> tuple[dict[tuple[str, str], dict[str, Any]], dict[str, set[tuple[str, str]]]]:
    sample_table = validate_tsv_table(paths["samples_master"], "samples_master.schema.json")
    samples_by_plink = {(row["FID"], row["IID"]): row for row in sample_table.rows}
    cohorts = validate_cohorts_frozen(
        paths["cohorts_frozen"], {row["SAMPLE_ID"] for row in sample_table.rows}
    )
    included_by_cohort: dict[str, set[tuple[str, str]]] = defaultdict(set)
    samples_by_id = {row["SAMPLE_ID"]: row for row in sample_table.rows}
    for row in cohorts.rows:
        if row["INCLUDED"]:
            sample = samples_by_id[row["SAMPLE_ID"]]
            included_by_cohort[row["COHORT_ID"]].add((sample["FID"], sample["IID"]))
    source_pairs = {
        source: {(row[0], row[1]) for row in fam_rows}
        for source, fam_rows in source_fams.items()
    }
    for policy in DATASET_POLICIES:
        cohort_id = policy["cohort_id"]
        keep_artifact_id = f"cohort_keep_{cohort_id}"
        keep_pairs = set(_keep_pairs(paths[keep_artifact_id]))
        if keep_pairs != included_by_cohort[cohort_id]:
            raise FinalQcInputError(f"keep_and_cohort_mismatch:{cohort_id}")
        if not keep_pairs <= source_pairs[policy["source"]]:
            raise FinalQcBlockError(f"cohort_absent_from_source_dataset:{cohort_id}")
    return samples_by_plink, included_by_cohort


def _validate_imiss(
    rows: list[dict[str, str]],
    fam_rows: list[list[str]],
) -> dict[tuple[str, str], dict[str, str]]:
    by_sample = {(row["FID"], row["IID"]): row for row in rows}
    expected = {(row[0], row[1]) for row in fam_rows}
    if set(by_sample) != expected:
        raise ExternalToolError("imiss_sample_set_mismatch")
    return by_sample


def _sample_rows(
    dataset_id: str,
    cohort_id: str,
    fam_rows: list[list[str]],
    imiss_by_sample: dict[tuple[str, str], dict[str, str]],
    samples_by_plink: dict[tuple[str, str], dict[str, Any]],
    threshold: float,
) -> tuple[list[dict[str, str]], set[tuple[str, str]]]:
    rows = []
    excluded: set[tuple[str, str]] = set()
    for fam_row in fam_rows:
        pair = (fam_row[0], fam_row[1])
        metric = imiss_by_sample[pair]
        n_miss = _integer_value(metric["N_MISS"], "N_MISS")
        n_geno = _integer_value(metric["N_GENO"], "N_GENO", minimum=1)
        f_miss = _float_value(metric["F_MISS"], "F_MISS")
        if f_miss is None or not 0 <= f_miss <= 1:
            raise ExternalToolError("invalid_plink_field:F_MISS")
        is_excluded = f_miss > threshold
        if is_excluded:
            excluded.add(pair)
        sample = samples_by_plink.get(pair)
        if sample is None:
            raise FinalQcInputError("cohort_sample_absent_from_registry")
        rows.append(
            {
                "DATASET_ID": dataset_id,
                "COHORT_ID": cohort_id,
                "SAMPLE_ID": sample["SAMPLE_ID"],
                "FID": pair[0],
                "IID": pair[1],
                "ARRAY_BATCH": sample["ARRAY_BATCH"] or "UNASSIGNED",
                "N_MISS": str(n_miss),
                "N_GENO": str(n_geno),
                "F_MISS": _decimal(f_miss),
                "QC_STATUS": "EXCLUDED" if is_excluded else "PASS",
                "EXCLUSION_CODE": "sample_missingness" if is_excluded else "",
                "PRESENT_AFTER_QC": str(not is_excluded).lower(),
            }
        )
    return rows, excluded


def _batch_rows(
    dataset_id: str,
    cohort_id: str,
    sample_rows: list[dict[str, str]],
    delta_threshold: float,
) -> list[dict[str, str]]:
    missingness_by_batch: dict[str, list[float]] = defaultdict(list)
    for row in sample_rows:
        missingness_by_batch[row["ARRAY_BATCH"]].append(float(row["F_MISS"]))
    means = {
        batch: sum(values) / len(values)
        for batch, values in missingness_by_batch.items()
    }
    best_mean = min(means.values())
    rows = []
    for batch in sorted(means):
        delta = means[batch] - best_mean
        alert = delta > delta_threshold
        rows.append(
            {
                "DATASET_ID": dataset_id,
                "COHORT_ID": cohort_id,
                "ARRAY_BATCH": batch,
                "SAMPLE_COUNT": str(len(missingness_by_batch[batch])),
                "MEAN_F_MISS": _decimal(means[batch]),
                "BATCH_DELTA_FROM_BEST": _decimal(delta),
                "QC_STATUS": "ALERT" if alert else "PASS",
                "ALERT_CODE": "batch_missingness_delta" if alert else "",
            }
        )
    return rows


def _hwe_by_variant(rows: list[dict[str, str]]) -> dict[str, float | None]:
    all_rows = [row for row in rows if row.get("TEST") == "ALL"]
    if not all_rows:
        raise ExternalToolError("hwe_all_rows_missing")
    by_variant: dict[str, float | None] = {}
    for row in all_rows:
        variant_id = row["SNP"]
        if variant_id in by_variant:
            raise ExternalToolError("duplicate_hwe_variant")
        by_variant[variant_id] = _float_value(row["P"], "HWE_P", nullable=True)
    return by_variant


def _variant_rows(
    policy: dict[str, Any],
    bim_rows: list[list[str]],
    lmiss_rows: list[dict[str, str]],
    frequency_rows: list[dict[str, str]],
    hwe_rows: list[dict[str, str]] | None,
    parameters: dict[str, Any],
    target_variant_id: str,
) -> tuple[list[dict[str, str]], set[str]]:
    lmiss_by_variant = {row["SNP"]: row for row in lmiss_rows}
    frequency_by_variant = {row["SNP"]: row for row in frequency_rows}
    variant_ids = {row[1] for row in bim_rows}
    if set(lmiss_by_variant) != variant_ids or set(frequency_by_variant) != variant_ids:
        raise ExternalToolError("variant_metric_set_mismatch")
    hwe_by_variant = _hwe_by_variant(hwe_rows) if hwe_rows is not None else {}
    if hwe_rows is not None and set(hwe_by_variant) != variant_ids:
        raise ExternalToolError("hwe_variant_set_mismatch")
    rows = []
    excluded: set[str] = set()
    for bim_row in bim_rows:
        variant_id = bim_row[1]
        lmiss = lmiss_by_variant[variant_id]
        frequency = frequency_by_variant[variant_id]
        n_miss = _integer_value(lmiss["N_MISS"], "N_MISS")
        n_geno = _integer_value(lmiss["N_GENO"], "N_GENO", minimum=1)
        f_miss = _float_value(lmiss["F_MISS"], "F_MISS")
        maf = _float_value(frequency["MAF"], "MAF", nullable=True)
        if (
            f_miss is None
            or not 0 <= f_miss <= 1
            or (maf is not None and not 0 <= maf <= 0.5)
        ):
            raise ExternalToolError("invalid_variant_metric")
        hwe_p = hwe_by_variant.get(variant_id)
        filter_codes: list[str] = []
        if f_miss > parameters[policy["variant_missingness_parameter"]]:
            filter_codes.append("variant_missingness")
        maf_parameter = policy["maf_parameter"]
        if maf_parameter is not None and (maf is None or maf < parameters[maf_parameter]):
            filter_codes.append(
                "low_maf_controls"
                if policy["cohort_id"] == "controls_unrelated"
                else "low_maf_target_region"
            )
        if policy["apply_hwe"]:
            if hwe_p is None:
                filter_codes.append("hwe_not_evaluable")
            elif hwe_p < parameters["controls_hwe_p_min"]:
                filter_codes.append("hwe_controls")
        is_target = variant_id == target_variant_id
        forced_exception = bool(is_target and policy["force_target"] and filter_codes)
        hard_filter_codes = [code for code in filter_codes if code != "hwe_not_evaluable"]
        if forced_exception:
            status = "RETAINED_EXCEPTION"
            present = True
        elif hard_filter_codes:
            status = "EXCLUDED"
            present = False
            excluded.add(variant_id)
        elif filter_codes:
            status = "ALERT"
            present = True
        else:
            status = "PASS"
            present = True
        rows.append(
            {
                "DATASET_ID": policy["dataset_id"],
                "COHORT_ID": policy["cohort_id"],
                "VARIANT_ID": variant_id,
                "CHROMOSOME": bim_row[0],
                "POSITION_BP": bim_row[3],
                "A1": frequency["A1"],
                "A2": frequency["A2"],
                "N_MISS": str(n_miss),
                "N_GENO": str(n_geno),
                "F_MISS": _decimal(f_miss),
                "MAF": _decimal(maf) if maf is not None else "",
                "HWE_P": _decimal(hwe_p) if hwe_p is not None else "",
                "HWE_STATUS": (
                    "EVALUATED"
                    if policy["apply_hwe"] and hwe_p is not None
                    else "NOT_EVALUATED_DATA"
                    if policy["apply_hwe"]
                    else "NOT_EVALUATED_POLICY"
                ),
                "QC_STATUS": status,
                "FILTER_CODES": ",".join(filter_codes),
                "IS_TARGET_VARIANT": str(is_target).lower(),
                "FILTER_EXCEPTION": str(forced_exception).lower(),
                "PRESENT_AFTER_QC": str(present).lower(),
            }
        )
    return rows, excluded


def _write_pairs(path: Path, pairs: set[tuple[str, str]]) -> None:
    path.write_text(
        "".join(f"{fid}\t{iid}\n" for fid, iid in sorted(pairs)),
        encoding="utf-8",
    )


def _write_variants(path: Path, variant_ids: set[str]) -> None:
    path.write_text(
        "".join(f"{variant_id}\n" for variant_id in sorted(variant_ids)),
        encoding="utf-8",
    )


def _set_id(prefix: str, values: list[str]) -> str:
    digest = hashlib.sha256("".join(f"{value}\n" for value in values).encode()).hexdigest()
    return f"{prefix}_{digest[:12]}"


def _write_descriptor(
    path: Path,
    prefix: Path,
    base_descriptor: dict[str, Any],
    dataset_id: str,
    fam_rows: list[list[str]],
    bim_rows: list[list[str]],
) -> dict[str, Any]:
    descriptor = {
        **base_descriptor,
        "dataset_id": dataset_id,
        "sample_set_id": _set_id(
            f"{dataset_id}_samples", [f"{row[0]}\t{row[1]}" for row in fam_rows]
        ),
        "variant_set_id": _set_id(
            f"{dataset_id}_variants", [row[1] for row in bim_rows]
        ),
        "sample_count": len(fam_rows),
        "variant_count": len(bim_rows),
        "source_format": "PLINK_FINAL_QC",
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
        path = prefix.with_suffix(f".{suffix}")
        artifacts.append(
            build_file_artifact(
                physical_path=path,
                published_path=str(published_dir / path.name),
                artifact_id=f"{descriptor['dataset_id']}_{suffix}",
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
            artifact_id=f"{descriptor['dataset_id']}_dataset",
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


def _write_tsv(path: Path, columns: tuple[str, ...], rows: list[dict[str, str]]) -> None:
    with path.open("w", encoding="utf-8", newline="") as output_file:
        writer = csv.DictWriter(output_file, fieldnames=columns, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows({column: row[column] for column in columns} for row in rows)


def _clean_plink_sidecars(prefix: Path) -> None:
    for suffix in ("log", "nosex", "hh"):
        path = prefix.with_suffix(f".{suffix}")
        if path.exists():
            path.unlink()


def execute(stage_inputs_path: Path, output_dir: Path) -> int:
    """Produit trois jeux QC finaux sans appliquer une politique unique à tous."""
    started_at = utc_now()
    started_clock = monotonic()
    stage_inputs = read_json(stage_inputs_path)
    validate_json_document(stage_inputs, "stage_inputs.schema.json")
    run_dir = output_dir.parent.parent
    config = load_pipeline_config(run_dir / "config.resolved.yaml")
    parameters = _parameters(stage_inputs["parameters"])
    artifact_ids = (
        "config_input_target_variant_metadata",
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
    )
    input_artifacts = {artifact_id: _artifact_by_id(stage_inputs, artifact_id) for artifact_id in artifact_ids}
    paths = {artifact_id: _validated_artifact_path(artifact, run_dir) for artifact_id, artifact in input_artifacts.items()}
    metadata = _load_target_metadata(paths["config_input_target_variant_metadata"])
    source_descriptors = {}
    source_fams = {}
    source_bims = {}
    for source in ("genomewide_pre_qc", "target_variant"):
        descriptor, fam_rows, bim_rows = _validate_dataset(source, paths, input_artifacts)
        source_descriptors[source] = descriptor
        source_fams[source] = fam_rows
        source_bims[source] = bim_rows
        if descriptor["assembly"] != metadata["assembly"]:
            raise FinalQcInputError(f"assembly_mismatch:{source}")
    target_descriptor = source_descriptors["target_variant"]
    if (
        target_descriptor["scope"] != "target_chromosome"
        or target_descriptor["target_chromosome"] != metadata["chromosome"]
    ):
        raise FinalQcInputError("target_dataset_scope_mismatch")
    if metadata["project_variant_id"] not in {row[1] for row in source_bims["target_variant"]}:
        raise FinalQcBlockError("target_variant_absent_from_target_dataset")
    samples_by_plink, _ = _sample_context(paths, source_fams)
    output_dir.mkdir(parents=True, exist_ok=True)
    executable = _resolve_executable(config["tools"]["plink"])
    plink_version = _plink_version(executable)
    all_sample_rows: list[dict[str, str]] = []
    all_variant_rows: list[dict[str, str]] = []
    all_batch_rows: list[dict[str, str]] = []
    final_descriptors: dict[str, dict[str, Any]] = {}
    final_prefixes: dict[str, Path] = {}
    descriptor_paths: dict[str, Path] = {}

    with tempfile.TemporaryDirectory(prefix="final_qc_", dir=output_dir) as work_directory:
        work_dir = Path(work_directory)
        for policy in DATASET_POLICIES:
            dataset_id = policy["dataset_id"]
            cohort_id = policy["cohort_id"]
            source = policy["source"]
            source_prefix = paths[f"{source}_bed"].with_suffix("")
            keep_path = paths[f"cohort_keep_{cohort_id}"]
            subset_prefix = work_dir / f"{dataset_id}_subset"
            _run_plink(
                executable,
                ["--bfile", str(source_prefix), "--keep", str(keep_path), "--make-bed", "--out", str(subset_prefix)],
                parameters["plink_timeout_seconds"],
            )
            subset_fam = _read_fam(subset_prefix.with_suffix(".fam"))
            expected_pairs = set(_keep_pairs(keep_path))
            if {(row[0], row[1]) for row in subset_fam} != expected_pairs:
                raise ExternalToolError(f"plink_keep_mismatch:{dataset_id}")
            sample_metrics_prefix = work_dir / f"{dataset_id}_sample_metrics"
            _run_plink(
                executable,
                ["--bfile", str(subset_prefix), "--missing", "--out", str(sample_metrics_prefix)],
                parameters["plink_timeout_seconds"],
            )
            imiss = _validate_imiss(
                _read_whitespace_table(sample_metrics_prefix.with_suffix(".imiss")),
                subset_fam,
            )
            sample_rows, excluded_samples = _sample_rows(
                dataset_id,
                cohort_id,
                subset_fam,
                imiss,
                samples_by_plink,
                parameters["sample_missingness_max"],
            )
            if len(excluded_samples) == len(subset_fam):
                raise FinalQcBlockError(f"all_samples_excluded:{dataset_id}")
            all_sample_rows.extend(sample_rows)
            all_batch_rows.extend(
                _batch_rows(
                    dataset_id,
                    cohort_id,
                    sample_rows,
                    parameters["batch_missingness_delta_alert"],
                )
            )
            sample_filtered_prefix = work_dir / f"{dataset_id}_sample_filtered"
            sample_filter_arguments = ["--bfile", str(subset_prefix)]
            if excluded_samples:
                remove_path = work_dir / f"{dataset_id}.remove"
                _write_pairs(remove_path, excluded_samples)
                sample_filter_arguments.extend(["--remove", str(remove_path)])
            sample_filter_arguments.extend(["--make-bed", "--out", str(sample_filtered_prefix)])
            _run_plink(executable, sample_filter_arguments, parameters["plink_timeout_seconds"])
            sample_filtered_fam = _read_fam(sample_filtered_prefix.with_suffix(".fam"))
            sample_filtered_bim = _read_bim(sample_filtered_prefix.with_suffix(".bim"))
            retained_sample_pairs = expected_pairs - excluded_samples
            if {(row[0], row[1]) for row in sample_filtered_fam} != retained_sample_pairs:
                raise ExternalToolError(f"sample_filter_mismatch:{dataset_id}")
            variant_metrics_prefix = work_dir / f"{dataset_id}_variant_metrics"
            variant_metric_arguments = [
                "--bfile",
                str(sample_filtered_prefix),
                "--missing",
                "--freq",
            ]
            if policy["apply_hwe"]:
                variant_metric_arguments.append("--hardy")
            variant_metric_arguments.extend(["--out", str(variant_metrics_prefix)])
            _run_plink(executable, variant_metric_arguments, parameters["plink_timeout_seconds"])
            hwe_rows = (
                _read_whitespace_table(variant_metrics_prefix.with_suffix(".hwe"))
                if policy["apply_hwe"]
                else None
            )
            variant_rows, excluded_variants = _variant_rows(
                policy,
                sample_filtered_bim,
                _read_whitespace_table(variant_metrics_prefix.with_suffix(".lmiss")),
                _read_whitespace_table(variant_metrics_prefix.with_suffix(".frq")),
                hwe_rows,
                parameters,
                metadata["project_variant_id"],
            )
            if len(excluded_variants) == len(sample_filtered_bim):
                raise FinalQcBlockError(f"all_variants_excluded:{dataset_id}")
            all_variant_rows.extend(variant_rows)
            final_prefix = output_dir / dataset_id
            final_arguments = ["--bfile", str(sample_filtered_prefix)]
            if excluded_variants:
                exclude_path = work_dir / f"{dataset_id}.exclude"
                _write_variants(exclude_path, excluded_variants)
                final_arguments.extend(["--exclude", str(exclude_path)])
            final_arguments.extend(["--make-bed", "--out", str(final_prefix)])
            _run_plink(executable, final_arguments, parameters["plink_timeout_seconds"])
            final_fam = _read_fam(final_prefix.with_suffix(".fam"))
            final_bim = _read_bim(final_prefix.with_suffix(".bim"))
            if len(final_fam) != len(sample_filtered_fam):
                raise ExternalToolError(f"final_sample_count_mismatch:{dataset_id}")
            expected_variant_ids = {
                row[1] for row in sample_filtered_bim
            } - excluded_variants
            if {row[1] for row in final_bim} != expected_variant_ids:
                raise ExternalToolError(f"final_variant_set_mismatch:{dataset_id}")
            if policy["force_target"] and metadata["project_variant_id"] not in {row[1] for row in final_bim}:
                raise FinalQcBlockError("target_variant_removed_from_final_dataset")
            descriptor_path = output_dir / f"{dataset_id}.dataset.json"
            descriptor = _write_descriptor(
                descriptor_path,
                final_prefix,
                source_descriptors[source],
                dataset_id,
                final_fam,
                final_bim,
            )
            final_descriptors[dataset_id] = descriptor
            final_prefixes[dataset_id] = final_prefix
            descriptor_paths[dataset_id] = descriptor_path
            _clean_plink_sidecars(final_prefix)

    sample_path = output_dir / "sample_qc_final.tsv"
    variant_path = output_dir / "variant_qc_final.tsv"
    batch_path = output_dir / "batch_qc_final.tsv"
    _write_tsv(sample_path, SAMPLE_COLUMNS, all_sample_rows)
    _write_tsv(variant_path, VARIANT_COLUMNS, all_variant_rows)
    _write_tsv(batch_path, BATCH_COLUMNS, all_batch_rows)
    validate_tsv_table(sample_path, "sample_qc_final.schema.json")
    validate_tsv_table(variant_path, "variant_qc_final.schema.json")
    validate_tsv_table(batch_path, "batch_qc_final.schema.json")
    retained_exceptions = sum(
        row["QC_STATUS"] == "RETAINED_EXCEPTION" for row in all_variant_rows
    )
    report_path = output_dir / "qc_final_report.json"
    report = {
        "schema_version": "1.0.0",
        "method_id": "cohort_specific_plink_final_qc_v1",
        "target_variant_id": metadata["project_variant_id"],
        "dataset_counts": {
            dataset_id: {
                "samples": descriptor["sample_count"],
                "variants": descriptor["variant_count"],
            }
            for dataset_id, descriptor in final_descriptors.items()
        },
        "excluded_sample_count": sum(row["QC_STATUS"] == "EXCLUDED" for row in all_sample_rows),
        "excluded_variant_count": sum(row["QC_STATUS"] == "EXCLUDED" for row in all_variant_rows),
        "target_variant_exception_count": retained_exceptions,
        "hwe_applied_datasets": ["controls_unrelated_qc"],
        "hwe_not_applied_datasets": [
            "target_carriers_independent_qc",
            "target_chromosome_all_qc",
        ],
    }
    atomic_write_json(report_path, report)
    output_artifacts = []
    for dataset_id in [policy["dataset_id"] for policy in DATASET_POLICIES]:
        output_artifacts.extend(
            _dataset_artifacts(
                stage_inputs,
                final_prefixes[dataset_id],
                descriptor_paths[dataset_id],
                final_descriptors[dataset_id],
            )
        )
    producer_stage = f"{stage_inputs['stage_id']}_{stage_inputs['stage_name']}"
    published_dir = PurePosixPath(stage_inputs["published_output_dir"])
    for artifact_id, artifact_type, path, schema_name, sensitivity in (
        ("sample_qc_final", "sample_qc_final", sample_path, "sample_qc_final.schema.json", "sensitive_genetic"),
        ("variant_qc_final", "variant_qc_final", variant_path, "variant_qc_final.schema.json", "sensitive_genetic"),
        ("batch_qc_final", "batch_qc_final", batch_path, "batch_qc_final.schema.json", "internal"),
        ("qc_final_report", "qc_final_report", report_path, None, "internal"),
    ):
        output_artifacts.append(
            build_file_artifact(
                physical_path=path,
                published_path=str(published_dir / path.name),
                artifact_id=artifact_id,
                artifact_type=artifact_type,
                media_type="application/json" if path.suffix == ".json" else "text/tab-separated-values",
                producer_stage=producer_stage,
                producer_signature=stage_inputs["signature"],
                schema_name=schema_name,
                schema_version="1.0.0" if schema_name else None,
                assembly=metadata["assembly"],
                sample_set_id=None,
                variant_set_id=None,
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
        "method_id": "cohort_specific_plink_final_qc_v1",
        "signature": stage_inputs["signature"],
        "started_at": started_at,
        "completed_at": utc_now(),
        "duration_seconds": monotonic() - started_clock,
        "inputs": list(input_artifacts.values()),
        "outputs": output_artifacts,
        "parameters": parameters,
        "tools": [{"tool": "plink", "configured": config["tools"]["plink"], "version": plink_version}],
        "counts": {
            "datasets": len(DATASET_POLICIES),
            "sample_metric_rows": len(all_sample_rows),
            "variant_metric_rows": len(all_variant_rows),
            "batch_metric_rows": len(all_batch_rows),
            "excluded_samples": report["excluded_sample_count"],
            "excluded_variants": report["excluded_variant_count"],
            "target_variant_exceptions": retained_exceptions,
        },
        "metrics": {"dataset_counts": report["dataset_counts"]},
        "exclusions": [],
        "warnings": (
            [{"code": "target_variant_forced_retention", "count": retained_exceptions}]
            if retained_exceptions
            else []
        ),
        "checks": [
            {"check": "cohort_keep_concordance", "status": "PASS"},
            {"check": "sample_missingness", "status": "PASS"},
            {"check": "control_hwe_maf", "status": "PASS"},
            {"check": "carrier_hwe_not_applied", "status": "PASS"},
            {"check": "target_variant_retained", "status": "PASS"},
            {"check": "final_dataset_integrity", "status": "PASS"},
        ],
        "known_limits": [
            "Le HWE est évalué uniquement chez les témoins indépendants et n'est jamais utilisé chez les porteurs.",
            "Le variant cible peut être conservé par exception documentée ; cette exception ne transforme pas un appel de mauvaise qualité en génotype fiable.",
            "Les alertes de lot sont descriptives et ne déclenchent pas d'exclusion automatique.",
            "Le QC final prépare les jeux analytiques mais ne remplace pas la préparation de région et la carte génétique de l'étape 11.",
        ],
        "expected_visualizations": [
            "control_hwe_distribution",
            "cohort_maf_comparison",
            "cohort_missingness_comparison",
            "target_chromosome_marker_density",
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
    except FinalQcBlockError as error:
        sys.stderr.write(f"{error}\n")
        return 4
    except (
        FinalQcInputError,
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
