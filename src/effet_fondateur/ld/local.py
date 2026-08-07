"""Normalisation auditée des sorties PLINK de LD local."""

from __future__ import annotations

import csv
import math
import shutil
import statistics
import subprocess
import tempfile
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Iterable

from effet_fondateur.audit import atomic_write_json
from effet_fondateur.contracts import validate_json_document, validate_tsv_table


COHORT_IDS = (
    "controls_unrelated",
    "target_carriers_independent",
    "family_noncarriers",
)
VARIANT_COLUMNS = (
    "COHORT_ID", "VARIANT_ORDER", "VARIANT_ID", "CHROMOSOME", "POSITION_BP",
    "POSITION_CM", "A1", "A2", "CALLED_SAMPLE_COUNT", "MAF", "POLYMORPHIC",
    "IS_TARGET_VARIANT", "VARIANT_STATUS",
)
PAIR_COLUMNS = (
    "COHORT_ID", "VARIANT_A", "VARIANT_B", "POSITION_A_BP", "POSITION_B_BP",
    "POSITION_A_CM", "POSITION_B_CM", "DISTANCE_BP", "DISTANCE_CM",
    "CALLED_SAMPLE_COUNT", "R2_GENOTYPE", "R2_HAPLOTYPE_ML", "D_PRIME_ABS",
    "INVOLVES_TARGET", "PAIR_STATUS",
)
SUMMARY_COLUMNS = (
    "COHORT_ID", "DISTANCE_BIN_ID", "DISTANCE_MIN_CM", "DISTANCE_MAX_CM",
    "SAMPLE_COUNT", "COHORT_STATUS", "TOTAL_PAIR_COUNT", "EVALUATED_PAIR_COUNT",
    "TARGET_PAIR_COUNT", "MEDIAN_R2_GENOTYPE", "MEDIAN_D_PRIME_ABS",
    "SUMMARY_STATUS",
)


class LocalLdError(RuntimeError):
    """Signale une entrée incohérente ou une sortie PLINK invalide."""


@dataclass(frozen=True)
class LocalLdPublication:
    """Chemins et comptes publiés par l'analyse locale."""

    variants_path: Path
    pairs_path: Path
    summary_path: Path
    analysis_summary_path: Path
    native_paths: tuple[Path, ...]
    cohort_statuses: dict[str, str]
    pair_count: int
    evaluated_pair_count: int


def _decimal(value: float) -> str:
    if abs(value) < 5e-16:
        value = 0.0
    return format(value, ".12g")


def _write_tsv(path: Path, columns: tuple[str, ...], rows: Iterable[dict[str, str]]) -> None:
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows({column: row[column] for column in columns} for row in rows)


def _read_whitespace(path: Path, *, allow_header_only: bool = False) -> list[dict[str, str]]:
    try:
        lines = [line.split() for line in path.read_text(encoding="utf-8").splitlines() if line]
    except OSError as error:
        raise LocalLdError(f"missing_plink_output:{path.name}") from error
    if not lines or (len(lines) == 1 and not allow_header_only):
        raise LocalLdError(f"empty_plink_output:{path.name}")
    header = lines[0]
    if any(len(row) != len(header) for row in lines[1:]):
        raise LocalLdError(f"malformed_plink_output:{path.name}")
    return [dict(zip(header, row, strict=True)) for row in lines[1:]]


def _run_plink(executable: str, arguments: list[str], timeout_seconds: int) -> None:
    try:
        completed = subprocess.run(
            [executable, *arguments], capture_output=True, check=False, text=True,
            timeout=timeout_seconds,
        )
    except (OSError, subprocess.TimeoutExpired) as error:
        raise LocalLdError("plink_local_ld_failed") from error
    if completed.returncode != 0:
        raise LocalLdError(f"plink_local_ld_failed:{completed.returncode}")


def _pair_key(first: str, second: str, order: dict[str, int]) -> tuple[str, str]:
    if first not in order or second not in order or first == second:
        raise LocalLdError("native_ld_variant_mismatch")
    return (first, second) if order[first] < order[second] else (second, first)


def _native_ld(
    path: Path, order: dict[str, int], *, dprime: bool,
) -> dict[tuple[str, str], tuple[float, float | None]]:
    values: dict[tuple[str, str], tuple[float, float | None]] = {}
    for row in _read_whitespace(path, allow_header_only=True):
        try:
            key = _pair_key(row["SNP_A"], row["SNP_B"], order)
            r2 = float(row["R2"])
            d_value = next(
                (float(row[name]) for name in ("DP", "DPRIME", "D_PRIME") if name in row),
                None,
            )
        except (KeyError, ValueError) as error:
            raise LocalLdError(f"malformed_native_ld:{path.name}") from error
        if not math.isfinite(r2) or not 0 <= r2 <= 1:
            raise LocalLdError(f"invalid_native_r2:{path.name}")
        if dprime and (d_value is None or not math.isfinite(d_value) or not 0 <= d_value <= 1):
            raise LocalLdError(f"invalid_native_dprime:{path.name}")
        if key in values:
            raise LocalLdError(f"duplicate_native_ld_pair:{path.name}")
        values[key] = (r2, d_value)
    return values


def _dosages(path: Path, expected_variants: set[str]) -> dict[str, frozenset[int]]:
    rows = _read_whitespace(path)
    observed: dict[str, frozenset[int]] = {}
    for row in rows:
        variant_id = row.get("SNP")
        if variant_id is None or variant_id in observed:
            raise LocalLdError("malformed_plink_traw")
        fixed_columns = {"CHR", "SNP", "(C)M", "POS", "COUNTED", "ALT"}
        sample_values = [value for key, value in row.items() if key not in fixed_columns]
        observed[variant_id] = frozenset(
            index for index, value in enumerate(sample_values) if value not in {"NA", "nan"}
        )
    if set(observed) != expected_variants:
        raise LocalLdError("traw_variant_set_mismatch")
    return observed


def _frequencies(path: Path, expected_variants: set[str]) -> dict[str, float | None]:
    result: dict[str, float | None] = {}
    for row in _read_whitespace(path):
        variant_id = row.get("SNP")
        raw_maf = row.get("MAF")
        if variant_id is None or variant_id in result or raw_maf is None:
            raise LocalLdError("malformed_plink_frequency")
        if raw_maf in {"NA", "nan"}:
            result[variant_id] = None
        else:
            try:
                maf = float(raw_maf)
            except ValueError as error:
                raise LocalLdError("malformed_plink_frequency") from error
            if not math.isfinite(maf) or not 0 <= maf <= 0.5:
                raise LocalLdError("invalid_plink_frequency")
            result[variant_id] = maf
    if set(result) != expected_variants:
        raise LocalLdError("frequency_variant_set_mismatch")
    return result


def _candidate_pairs(
    variants: list[dict[str, Any]], max_bp: int, max_cm: float,
) -> Iterable[tuple[dict[str, Any], dict[str, Any]]]:
    for first_index, first in enumerate(variants):
        for second in variants[first_index + 1:]:
            distance_bp = int(second["POSITION_BP"]) - int(first["POSITION_BP"])
            if distance_bp > max_bp:
                break
            distance_cm = float(second["POSITION_CM"]) - float(first["POSITION_CM"])
            if distance_cm <= max_cm:
                yield first, second


def _bin_identifier(index: int, lower: float, upper: float) -> str:
    return f"bin_{index:02d}_{_decimal(lower).replace('.', 'p')}_{_decimal(upper).replace('.', 'p')}_cm"


def _cohort_status(sample_count: int, minimum_samples: int, primary_samples: int) -> str:
    if sample_count < minimum_samples:
        return "NOT_EVALUATED"
    if sample_count < primary_samples:
        return "EXPLORATORY_SMALL_N"
    return "DESCRIPTIVE_PRIMARY"


def _summary_rows(
    cohort_id: str,
    sample_count: int,
    cohort_status: str,
    pair_rows: list[dict[str, str]],
    distance_bins_cm: list[float],
    minimum_pairs_per_bin: int,
    non_evaluation_reason: str | None,
) -> list[dict[str, str]]:
    result = []
    lower = 0.0
    for index, upper in enumerate(distance_bins_cm, start=1):
        rows = [
            row for row in pair_rows
            if lower <= float(row["DISTANCE_CM"]) <= upper
            and (index == 1 or float(row["DISTANCE_CM"]) > lower)
        ]
        evaluated = [row for row in rows if row["PAIR_STATUS"] == "EVALUATED"]
        if non_evaluation_reason is not None:
            summary_status = non_evaluation_reason
        elif len(evaluated) < minimum_pairs_per_bin:
            summary_status = "INSUFFICIENT_PAIRS"
        else:
            summary_status = "EVALUATED"
        r2_values = [float(row["R2_GENOTYPE"]) for row in evaluated]
        dprime_values = [float(row["D_PRIME_ABS"]) for row in evaluated]
        result.append({
            "COHORT_ID": cohort_id,
            "DISTANCE_BIN_ID": _bin_identifier(index, lower, upper),
            "DISTANCE_MIN_CM": _decimal(lower),
            "DISTANCE_MAX_CM": _decimal(upper),
            "SAMPLE_COUNT": str(sample_count),
            "COHORT_STATUS": cohort_status,
            "TOTAL_PAIR_COUNT": str(len(rows)),
            "EVALUATED_PAIR_COUNT": str(len(evaluated)),
            "TARGET_PAIR_COUNT": str(sum(row["INVOLVES_TARGET"] == "true" for row in evaluated)),
            "MEDIAN_R2_GENOTYPE": _decimal(statistics.median(r2_values)) if summary_status == "EVALUATED" else "",
            "MEDIAN_D_PRIME_ABS": _decimal(statistics.median(dprime_values)) if summary_status == "EVALUATED" else "",
            "SUMMARY_STATUS": summary_status,
        })
        lower = upper
    return result


def publish_local_ld(
    *,
    plink_executable: str,
    dataset_prefix: Path,
    variants: list[dict[str, Any]],
    cohort_keep_paths: dict[str, Path],
    cohort_sample_counts: dict[str, int],
    cohort_independent: dict[str, bool],
    output_dir: Path,
    minimum_samples: int,
    primary_samples: int,
    minimum_called_samples_per_pair: int,
    minimum_variant_maf: float,
    minimum_pairs_per_bin: int,
    max_pair_distance_bp: int,
    max_pair_distance_cm: float,
    distance_bins_cm: list[float],
    max_pair_count: int,
    timeout_seconds: int,
) -> LocalLdPublication:
    """Exécute PLINK et publie des tables comparables sans inférence fondatrice."""
    if len(variants) < 2:
        raise LocalLdError("insufficient_region_variants")
    expected_ids = {row["VARIANT_ID"] for row in variants}
    if len(expected_ids) != len(variants):
        raise LocalLdError("duplicate_region_variant")
    target_ids = [row["VARIANT_ID"] for row in variants if row["IS_TARGET_VARIANT"]]
    if len(target_ids) != 1:
        raise LocalLdError("target_variant_missing_or_ambiguous")
    candidate_count = sum(1 for _ in _candidate_pairs(variants, max_pair_distance_bp, max_pair_distance_cm))
    if candidate_count > max_pair_count:
        raise LocalLdError("configured_ld_pair_limit_exceeded")

    output_dir.mkdir(parents=True, exist_ok=True)
    variant_rows: list[dict[str, str]] = []
    all_pair_rows: list[dict[str, str]] = []
    all_summary_rows: list[dict[str, str]] = []
    native_paths: list[Path] = []
    cohort_statuses: dict[str, str] = {}

    with tempfile.TemporaryDirectory(prefix="local_ld_", dir=output_dir) as temporary:
        work_dir = Path(temporary)
        for cohort_id in COHORT_IDS:
            sample_count = cohort_sample_counts[cohort_id]
            independent = cohort_independent[cohort_id]
            cohort_status = _cohort_status(sample_count, minimum_samples, primary_samples)
            non_evaluation_reason = None
            if not independent and sample_count:
                cohort_status = "NOT_EVALUATED"
                non_evaluation_reason = "COHORT_NOT_INDEPENDENT"
            elif cohort_status == "NOT_EVALUATED":
                non_evaluation_reason = "COHORT_TOO_SMALL"
            cohort_statuses[cohort_id] = cohort_status

            if non_evaluation_reason is not None:
                variant_status = non_evaluation_reason
                for variant in variants:
                    variant_rows.append({
                        "COHORT_ID": cohort_id, "VARIANT_ORDER": str(variant["VARIANT_ORDER"]),
                        "VARIANT_ID": variant["VARIANT_ID"], "CHROMOSOME": variant["CHROMOSOME"],
                        "POSITION_BP": str(variant["POSITION_BP"]), "POSITION_CM": _decimal(float(variant["POSITION_CM"])),
                        "A1": variant["A1"], "A2": variant["A2"], "CALLED_SAMPLE_COUNT": "",
                        "MAF": "", "POLYMORPHIC": "false",
                        "IS_TARGET_VARIANT": str(bool(variant["IS_TARGET_VARIANT"])).lower(),
                        "VARIANT_STATUS": variant_status,
                    })
                all_summary_rows.extend(_summary_rows(
                    cohort_id, sample_count, cohort_status, [], distance_bins_cm,
                    minimum_pairs_per_bin, non_evaluation_reason,
                ))
                continue

            cohort_prefix = work_dir / cohort_id
            report_prefix = work_dir / f"{cohort_id}_metrics"
            _run_plink(plink_executable, [
                "--bfile", str(dataset_prefix), "--keep", str(cohort_keep_paths[cohort_id]),
                "--freq", "--recode", "A-transpose", "--out", str(report_prefix),
            ], timeout_seconds)
            frequencies = _frequencies(report_prefix.with_suffix(".frq"), expected_ids)
            called = _dosages(report_prefix.with_suffix(".traw"), expected_ids)
            common_arguments = [
                "--bfile", str(dataset_prefix), "--keep", str(cohort_keep_paths[cohort_id]),
                "--ld-window", str(len(variants) + 1),
                "--ld-window-kb", str(math.ceil(max_pair_distance_bp / 1000)),
                "--ld-window-cm", _decimal(max_pair_distance_cm), "--ld-window-r2", "0",
            ]
            genotype_prefix = work_dir / f"{cohort_id}_r2_genotype"
            dprime_prefix = work_dir / f"{cohort_id}_dprime"
            _run_plink(plink_executable, [*common_arguments, "--r2", "with-freqs", "--out", str(genotype_prefix)], timeout_seconds)
            _run_plink(plink_executable, [*common_arguments, "--r2", "dprime", "with-freqs", "--out", str(dprime_prefix)], timeout_seconds)
            order = {row["VARIANT_ID"]: int(row["VARIANT_ORDER"]) for row in variants}
            genotype_ld = _native_ld(genotype_prefix.with_suffix(".ld"), order, dprime=False)
            dprime_ld = _native_ld(dprime_prefix.with_suffix(".ld"), order, dprime=True)

            for source_path, suffix in (
                (genotype_prefix.with_suffix(".ld"), "r2_genotype.ld"),
                (dprime_prefix.with_suffix(".ld"), "dprime.ld"),
            ):
                destination = output_dir / f"{cohort_id}.{suffix}"
                shutil.copyfile(source_path, destination)
                native_paths.append(destination)

            for variant in variants:
                variant_id = variant["VARIANT_ID"]
                maf = frequencies[variant_id]
                polymorphic = maf is not None and maf > 0
                variant_status = (
                    "MONOMORPHIC" if not polymorphic
                    else "LOW_MAF" if maf < minimum_variant_maf
                    else "EVALUATED"
                )
                variant_rows.append({
                    "COHORT_ID": cohort_id, "VARIANT_ORDER": str(variant["VARIANT_ORDER"]),
                    "VARIANT_ID": variant_id, "CHROMOSOME": variant["CHROMOSOME"],
                    "POSITION_BP": str(variant["POSITION_BP"]), "POSITION_CM": _decimal(float(variant["POSITION_CM"])),
                    "A1": variant["A1"], "A2": variant["A2"],
                    "CALLED_SAMPLE_COUNT": str(len(called[variant_id])),
                    "MAF": _decimal(maf) if maf is not None else "",
                    "POLYMORPHIC": str(polymorphic).lower(),
                    "IS_TARGET_VARIANT": str(bool(variant["IS_TARGET_VARIANT"])).lower(),
                    "VARIANT_STATUS": variant_status,
                })

            cohort_pair_rows: list[dict[str, str]] = []
            for first, second in _candidate_pairs(variants, max_pair_distance_bp, max_pair_distance_cm):
                first_id, second_id = first["VARIANT_ID"], second["VARIANT_ID"]
                key = (first_id, second_id)
                called_count = len(called[first_id] & called[second_id])
                first_maf = frequencies[first_id]
                second_maf = frequencies[second_id]
                polymorphic = bool(first_maf) and bool(second_maf)
                maf_eligible = (
                    first_maf is not None and second_maf is not None
                    and first_maf >= minimum_variant_maf
                    and second_maf >= minimum_variant_maf
                )
                genotype_value = genotype_ld.get(key)
                dprime_value = dprime_ld.get(key)
                if not polymorphic:
                    status = "MONOMORPHIC"
                elif not maf_eligible:
                    status = "LOW_MAF"
                elif called_count < minimum_called_samples_per_pair:
                    status = "INSUFFICIENT_CALLS"
                elif genotype_value is None or dprime_value is None:
                    status = "TOOL_NOT_REPORTED"
                else:
                    status = "EVALUATED"
                distance_bp = int(second["POSITION_BP"]) - int(first["POSITION_BP"])
                distance_cm = float(second["POSITION_CM"]) - float(first["POSITION_CM"])
                row = {
                    "COHORT_ID": cohort_id, "VARIANT_A": first_id, "VARIANT_B": second_id,
                    "POSITION_A_BP": str(first["POSITION_BP"]), "POSITION_B_BP": str(second["POSITION_BP"]),
                    "POSITION_A_CM": _decimal(float(first["POSITION_CM"])), "POSITION_B_CM": _decimal(float(second["POSITION_CM"])),
                    "DISTANCE_BP": str(distance_bp), "DISTANCE_CM": _decimal(distance_cm),
                    "CALLED_SAMPLE_COUNT": str(called_count),
                    "R2_GENOTYPE": _decimal(genotype_value[0]) if status == "EVALUATED" else "",
                    "R2_HAPLOTYPE_ML": _decimal(dprime_value[0]) if status == "EVALUATED" else "",
                    "D_PRIME_ABS": _decimal(dprime_value[1]) if status == "EVALUATED" else "",
                    "INVOLVES_TARGET": str(first_id in target_ids or second_id in target_ids).lower(),
                    "PAIR_STATUS": status,
                }
                cohort_pair_rows.append(row)
            all_pair_rows.extend(cohort_pair_rows)
            all_summary_rows.extend(_summary_rows(
                cohort_id, sample_count, cohort_status, cohort_pair_rows,
                distance_bins_cm, minimum_pairs_per_bin, None,
            ))

    variants_path = output_dir / "local_ld_variants.tsv"
    pairs_path = output_dir / "local_ld_pairs.tsv"
    summary_path = output_dir / "local_ld_summary.tsv"
    _write_tsv(variants_path, VARIANT_COLUMNS, variant_rows)
    _write_tsv(pairs_path, PAIR_COLUMNS, all_pair_rows)
    _write_tsv(summary_path, SUMMARY_COLUMNS, all_summary_rows)
    validate_tsv_table(variants_path, "local_ld_variants.schema.json")
    validate_tsv_table(pairs_path, "local_ld_pairs.schema.json")
    validate_tsv_table(summary_path, "local_ld_summary.schema.json")
    analysis_summary_path = output_dir / "local_ld_analysis_summary.json"
    analysis_summary = {
        "schema_version": "1.0.0", "method_id": "plink19_local_ld_secondary_v1",
        "analysis_role": "SECONDARY_DESCRIPTIVE", "primary_cohort": "controls_unrelated",
        "cohort_statuses": cohort_statuses, "variant_count": len(variants),
        "target_variant_id": target_ids[0],
        "r2_definition": "GENOTYPE_ALLELE_COUNT_CORRELATION_SQUARED",
        "dprime_definition": "ABSOLUTE_HAPLOTYPE_ML", "feeds_variant_age": False,
    }
    validate_json_document(analysis_summary, "local_ld_analysis_summary.schema.json")
    atomic_write_json(analysis_summary_path, analysis_summary)
    return LocalLdPublication(
        variants_path=variants_path, pairs_path=pairs_path, summary_path=summary_path,
        analysis_summary_path=analysis_summary_path, native_paths=tuple(native_paths),
        cohort_statuses=cohort_statuses, pair_count=len(all_pair_rows),
        evaluated_pair_count=sum(row["PAIR_STATUS"] == "EVALUATED" for row in all_pair_rows),
    )
