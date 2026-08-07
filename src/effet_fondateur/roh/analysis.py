"""Exécution et normalisation des ROH PLINK 1.9."""

from __future__ import annotations

import csv
import math
import shutil
import statistics
import subprocess
import tempfile
from collections import defaultdict
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Iterable

from effet_fondateur.audit import atomic_write_json
from effet_fondateur.contracts import validate_json_document, validate_tsv_table


SEGMENT_COLUMNS = (
    "SCOPE", "DATASET_ID", "COHORT_ID", "SAMPLE_ID", "FID", "IID",
    "CHROMOSOME", "SNP1", "SNP2", "POS1_BP", "POS2_BP", "LENGTH_KB",
    "NSNP", "DENSITY_KB_PER_SNP", "PHOM", "PHET", "INVOLVES_TARGET",
)
BURDEN_COLUMNS = (
    "COHORT_ID", "SAMPLE_ID", "FID", "IID", "N_ROH", "TOTAL_ROH_KB",
    "MEAN_ROH_KB", "MAX_ROH_KB", "TOTAL_1P5_TO_4MB_KB",
    "TOTAL_4_TO_8MB_KB", "TOTAL_AT_LEAST_8MB_KB",
    "AUTOZYGOSITY_DENOMINATOR_KB", "AUTOZYGOSITY_DENOMINATOR_SOURCE", "F_ROH",
    "ANALYSIS_STATUS",
)
TARGET_COLUMNS = (
    "SAMPLE_ID", "FID", "IID", "TARGET_VARIANT_ID", "TARGET_BP",
    "TARGET_GENOTYPE", "TARGET_ZYGOSITY", "ANALYTIC_ROLES", "TARGET_IN_ROH",
    "COVERING_SEGMENT_COUNT", "ROH_POS1_BP", "ROH_POS2_BP", "ROH_LENGTH_KB",
    "INTERPRETATION_STATUS",
)
SUMMARY_COLUMNS = (
    "SCOPE", "COHORT_ID", "SAMPLE_COUNT", "EVALUATED_SAMPLE_COUNT",
    "MEDIAN_N_ROH", "MEDIAN_TOTAL_ROH_KB", "MEDIAN_MAX_ROH_KB",
    "TARGET_IN_ROH_COUNT", "COHORT_STATUS",
)


class RohAnalysisError(RuntimeError):
    """Signale une entrée incohérente ou une sortie PLINK invalide."""


@dataclass(frozen=True)
class RohPublication:
    """Artefacts produits par l'analyse ROH."""

    segments_path: Path
    burden_path: Path
    target_path: Path
    cohort_summary_path: Path
    analysis_summary_path: Path
    native_paths: tuple[Path, ...]
    scope_statuses: dict[str, str]
    segment_count: int
    target_in_roh_count: int


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
        raise RohAnalysisError(f"missing_plink_roh_output:{path.name}") from error
    if not lines or (len(lines) == 1 and not allow_header_only):
        raise RohAnalysisError(f"empty_plink_roh_output:{path.name}")
    header = lines[0]
    if any(len(row) != len(header) for row in lines[1:]):
        raise RohAnalysisError(f"malformed_plink_roh_output:{path.name}")
    return [dict(zip(header, row, strict=True)) for row in lines[1:]]


def _run_plink(executable: str, arguments: list[str], timeout_seconds: int) -> None:
    try:
        completed = subprocess.run(
            [executable, *arguments], capture_output=True, check=False, text=True,
            timeout=timeout_seconds,
        )
    except (OSError, subprocess.TimeoutExpired) as error:
        raise RohAnalysisError("plink_roh_failed") from error
    if completed.returncode != 0:
        raise RohAnalysisError(f"plink_roh_failed:{completed.returncode}")


def _plink_arguments(parameters: dict[str, Any]) -> list[str]:
    return [
        "--homozyg",
        "--homozyg-snp", str(parameters["minimum_roh_snps"]),
        "--homozyg-kb", _decimal(parameters["minimum_roh_kb"]),
        "--homozyg-density", _decimal(parameters["maximum_density_kb_per_snp"]),
        "--homozyg-gap", _decimal(parameters["maximum_gap_kb"]),
        "--homozyg-window-snp", str(parameters["window_snps"]),
        "--homozyg-window-het", str(parameters["window_max_heterozygotes"]),
        "--homozyg-window-missing", str(parameters["window_max_missing"]),
        "--homozyg-window-threshold", _decimal(parameters["window_threshold"]),
    ]


def _run_scope(
    *, executable: str, dataset_prefix: Path, output_prefix: Path,
    parameters: dict[str, Any], timeout_seconds: int, extract_path: Path | None,
) -> tuple[list[dict[str, str]], list[dict[str, str]], tuple[Path, ...]]:
    arguments = ["--bfile", str(dataset_prefix)]
    if extract_path is not None:
        arguments.extend(["--extract", str(extract_path)])
    arguments.extend([*_plink_arguments(parameters), "--out", str(output_prefix)])
    _run_plink(executable, arguments, timeout_seconds)
    segments = _read_whitespace(output_prefix.with_suffix(".hom"), allow_header_only=True)
    individual = _read_whitespace(output_prefix.with_suffix(".hom.indiv"))
    native_paths = tuple(
        output_prefix.with_suffix(suffix)
        for suffix in (".hom", ".hom.indiv", ".hom.summary")
    )
    if any(not path.is_file() for path in native_paths):
        raise RohAnalysisError("missing_plink_roh_native_output")
    return segments, individual, native_paths


def _validate_number(raw: str, field: str, *, minimum: float = 0) -> float:
    try:
        value = float(raw)
    except ValueError as error:
        raise RohAnalysisError(f"invalid_plink_roh_field:{field}") from error
    if not math.isfinite(value) or value < minimum:
        raise RohAnalysisError(f"invalid_plink_roh_field:{field}")
    return value


def _normalise_segments(
    *, scope: str, dataset_id: str, cohort_id: str,
    raw_rows: list[dict[str, str]], sample_by_pair: dict[tuple[str, str], str],
    target_chromosome: str, target_bp: int,
) -> list[dict[str, str]]:
    result = []
    for raw in raw_rows:
        try:
            pair = (raw["FID"], raw["IID"])
            sample_id = sample_by_pair[pair]
            chromosome = raw["CHR"]
            pos1, pos2 = int(raw["POS1"]), int(raw["POS2"])
            nsnp = int(raw["NSNP"])
        except (KeyError, ValueError) as error:
            raise RohAnalysisError("malformed_plink_roh_segment") from error
        length = _validate_number(raw["KB"], "KB")
        density = _validate_number(raw["DENSITY"], "DENSITY")
        phom = _validate_number(raw["PHOM"], "PHOM")
        phet = _validate_number(raw["PHET"], "PHET")
        if pos1 <= 0 or pos2 < pos1 or nsnp <= 0 or phom > 1 or phet > 1:
            raise RohAnalysisError("invalid_plink_roh_segment")
        result.append({
            "SCOPE": scope, "DATASET_ID": dataset_id, "COHORT_ID": cohort_id,
            "SAMPLE_ID": sample_id, "FID": pair[0], "IID": pair[1],
            "CHROMOSOME": chromosome, "SNP1": raw["SNP1"], "SNP2": raw["SNP2"],
            "POS1_BP": str(pos1), "POS2_BP": str(pos2), "LENGTH_KB": _decimal(length),
            "NSNP": str(nsnp), "DENSITY_KB_PER_SNP": _decimal(density),
            "PHOM": _decimal(phom), "PHET": _decimal(phet),
            "INVOLVES_TARGET": str(chromosome == target_chromosome and pos1 <= target_bp <= pos2).lower(),
        })
    return result


def _burden_rows(
    *, cohort_id: str, fam_pairs: list[tuple[str, str]],
    sample_by_pair: dict[tuple[str, str], str], segments: list[dict[str, str]],
    evaluated: bool, denominator_kb: float | None, denominator_source: str | None,
) -> list[dict[str, str]]:
    by_sample: dict[str, list[float]] = defaultdict(list)
    for segment in segments:
        by_sample[segment["SAMPLE_ID"]].append(float(segment["LENGTH_KB"]))
    result = []
    for pair in fam_pairs:
        sample_id = sample_by_pair[pair]
        lengths = by_sample[sample_id]
        if not evaluated:
            values = {name: "" for name in (
                "TOTAL_ROH_KB", "MEAN_ROH_KB", "MAX_ROH_KB",
                "TOTAL_1P5_TO_4MB_KB", "TOTAL_4_TO_8MB_KB",
                "TOTAL_AT_LEAST_8MB_KB", "AUTOZYGOSITY_DENOMINATOR_KB", "F_ROH",
            )}
            values["AUTOZYGOSITY_DENOMINATOR_SOURCE"] = ""
        else:
            total = sum(lengths)
            values = {
                "TOTAL_ROH_KB": _decimal(total),
                "MEAN_ROH_KB": _decimal(statistics.mean(lengths)) if lengths else "0",
                "MAX_ROH_KB": _decimal(max(lengths)) if lengths else "0",
                "TOTAL_1P5_TO_4MB_KB": _decimal(sum(value for value in lengths if 1500 <= value < 4000)),
                "TOTAL_4_TO_8MB_KB": _decimal(sum(value for value in lengths if 4000 <= value < 8000)),
                "TOTAL_AT_LEAST_8MB_KB": _decimal(sum(value for value in lengths if value >= 8000)),
                "AUTOZYGOSITY_DENOMINATOR_KB": _decimal(denominator_kb) if denominator_kb else "",
                "AUTOZYGOSITY_DENOMINATOR_SOURCE": denominator_source or "",
                "F_ROH": _decimal(total / denominator_kb) if denominator_kb else "",
            }
        result.append({
            "COHORT_ID": cohort_id, "SAMPLE_ID": sample_id, "FID": pair[0], "IID": pair[1],
            "N_ROH": str(len(lengths)) if evaluated else "", **values,
            "ANALYSIS_STATUS": "EVALUATED" if evaluated else "NOT_EVALUATED_SCOPE",
        })
    return result


def _cohort_status(sample_count: int, minimum: int, primary: int) -> str:
    if sample_count < minimum:
        return "NOT_EVALUATED"
    return "DESCRIPTIVE_PRIMARY" if sample_count >= primary else "EXPLORATORY_SMALL_N"


def _genomewide_summary(
    cohort_id: str, burden: list[dict[str, str]], minimum: int, primary: int,
) -> dict[str, str]:
    status = _cohort_status(len(burden), minimum, primary)
    evaluated = [row for row in burden if row["ANALYSIS_STATUS"] == "EVALUATED"]
    numeric_allowed = status != "NOT_EVALUATED" and len(evaluated) == len(burden)
    return {
        "SCOPE": "GENOMEWIDE_BURDEN", "COHORT_ID": cohort_id,
        "SAMPLE_COUNT": str(len(burden)), "EVALUATED_SAMPLE_COUNT": str(len(evaluated)),
        "MEDIAN_N_ROH": _decimal(statistics.median(int(row["N_ROH"]) for row in evaluated)) if numeric_allowed else "",
        "MEDIAN_TOTAL_ROH_KB": _decimal(statistics.median(float(row["TOTAL_ROH_KB"]) for row in evaluated)) if numeric_allowed else "",
        "MEDIAN_MAX_ROH_KB": _decimal(statistics.median(float(row["MAX_ROH_KB"]) for row in evaluated)) if numeric_allowed else "",
        "TARGET_IN_ROH_COUNT": "0", "COHORT_STATUS": status,
    }


def publish_roh(
    *, plink_executable: str,
    genomewide_datasets: dict[str, dict[str, Any]], target_dataset: dict[str, Any],
    common_variant_ids: list[str], common_autosome_count: int,
    target_variant_count: int, sample_by_pair: dict[tuple[str, str], str],
    target_genotypes: dict[str, str], roles_by_sample: dict[str, set[str]],
    target_variant_id: str, target_chromosome: str, target_bp: int,
    output_dir: Path, parameters: dict[str, Any],
) -> RohPublication:
    """Publie les ROH individuels et leur intersection exacte avec la cible."""
    output_dir.mkdir(parents=True, exist_ok=True)
    genome_evaluated = (
        len(common_variant_ids) >= parameters["minimum_genomewide_variants"]
        and common_autosome_count >= parameters["minimum_genomewide_autosomes"]
    )
    target_evaluated = target_variant_count >= parameters["minimum_target_chromosome_variants"]
    scope_statuses = {
        "GENOMEWIDE_BURDEN": "EVALUATED" if genome_evaluated else "NOT_EVALUATED_DENSITY_OR_COVERAGE",
        "TARGET_CHROMOSOME": "EVALUATED" if target_evaluated else "NOT_EVALUATED_DENSITY_OR_COVERAGE",
    }
    all_segments: list[dict[str, str]] = []
    burden_rows: list[dict[str, str]] = []
    native_paths: list[Path] = []

    with tempfile.TemporaryDirectory(prefix="roh_", dir=output_dir) as temporary:
        work_dir = Path(temporary)
        extract_path = work_dir / "genomewide_common.extract"
        extract_path.write_text("".join(f"{variant_id}\n" for variant_id in common_variant_ids), encoding="utf-8")
        for cohort_id, dataset in genomewide_datasets.items():
            cohort_segments: list[dict[str, str]] = []
            if genome_evaluated:
                raw_segments, _, raw_native = _run_scope(
                    executable=plink_executable, dataset_prefix=dataset["prefix"],
                    output_prefix=work_dir / f"{cohort_id}_genomewide",
                    parameters=parameters, timeout_seconds=parameters["plink_timeout_seconds"],
                    extract_path=extract_path,
                )
                cohort_segments = _normalise_segments(
                    scope="GENOMEWIDE_BURDEN", dataset_id=dataset["dataset_id"],
                    cohort_id=cohort_id, raw_rows=raw_segments, sample_by_pair=sample_by_pair,
                    target_chromosome=target_chromosome, target_bp=target_bp,
                )
                for source in raw_native:
                    destination = output_dir / f"{cohort_id}.genomewide{''.join(source.suffixes[-2:]) if source.name.endswith('.hom.indiv') or source.name.endswith('.hom.summary') else '.hom'}"
                    shutil.copyfile(source, destination)
                    native_paths.append(destination)
            all_segments.extend(cohort_segments)
            burden_rows.extend(_burden_rows(
                cohort_id=cohort_id, fam_pairs=dataset["fam_pairs"],
                sample_by_pair=sample_by_pair, segments=cohort_segments,
                evaluated=genome_evaluated,
                denominator_kb=parameters.get("autosomal_denominator_kb"),
                denominator_source=parameters.get("autosomal_denominator_source"),
            ))

        target_segments: list[dict[str, str]] = []
        if target_evaluated:
            raw_segments, _, raw_native = _run_scope(
                executable=plink_executable, dataset_prefix=target_dataset["prefix"],
                output_prefix=work_dir / "target_chromosome", parameters=parameters,
                timeout_seconds=parameters["plink_timeout_seconds"], extract_path=None,
            )
            target_segments = _normalise_segments(
                scope="TARGET_CHROMOSOME", dataset_id=target_dataset["dataset_id"],
                cohort_id="target_chromosome_all_qc", raw_rows=raw_segments,
                sample_by_pair=sample_by_pair, target_chromosome=target_chromosome,
                target_bp=target_bp,
            )
            for source in raw_native:
                destination = output_dir / f"target_chromosome{''.join(source.suffixes[-2:]) if source.name.endswith('.hom.indiv') or source.name.endswith('.hom.summary') else '.hom'}"
                shutil.copyfile(source, destination)
                native_paths.append(destination)
        all_segments.extend(target_segments)

    covering_by_sample: dict[str, list[dict[str, str]]] = defaultdict(list)
    for segment in target_segments:
        if segment["INVOLVES_TARGET"] == "true":
            covering_by_sample[segment["SAMPLE_ID"]].append(segment)
    target_rows = []
    for pair in target_dataset["fam_pairs"]:
        sample_id = sample_by_pair[pair]
        genotype = target_genotypes[sample_id]
        alleles = genotype.split("/")
        zygosity = "HOMOZYGOUS" if alleles[0] == alleles[1] else "HETEROZYGOUS"
        covering = covering_by_sample[sample_id]
        if len(covering) > 1:
            raise RohAnalysisError("multiple_target_covering_roh_for_sample")
        segment = covering[0] if covering else None
        if not target_evaluated:
            interpretation = "NOT_EVALUATED_SCOPE"
        elif segment is None:
            interpretation = "NOT_IN_ROH"
        elif zygosity == "HETEROZYGOUS":
            interpretation = "IN_ROH_HETEROZYGOUS_WARNING"
        else:
            interpretation = "IN_ROH_HOMOZYGOUS"
        target_rows.append({
            "SAMPLE_ID": sample_id, "FID": pair[0], "IID": pair[1],
            "TARGET_VARIANT_ID": target_variant_id, "TARGET_BP": str(target_bp),
            "TARGET_GENOTYPE": genotype, "TARGET_ZYGOSITY": zygosity,
            "ANALYTIC_ROLES": ",".join(sorted(roles_by_sample.get(sample_id, set()))),
            "TARGET_IN_ROH": str(segment is not None).lower(),
            "COVERING_SEGMENT_COUNT": str(len(covering)),
            "ROH_POS1_BP": segment["POS1_BP"] if segment else "",
            "ROH_POS2_BP": segment["POS2_BP"] if segment else "",
            "ROH_LENGTH_KB": segment["LENGTH_KB"] if segment else "",
            "INTERPRETATION_STATUS": interpretation,
        })

    summary_rows = []
    for cohort_id in genomewide_datasets:
        rows = [row for row in burden_rows if row["COHORT_ID"] == cohort_id]
        summary_rows.append(_genomewide_summary(
            cohort_id, rows, parameters["minimum_group_samples"], parameters["minimum_primary_samples"]
        ))
    for cohort_id in ("controls_unrelated", "target_carriers_independent", "family_noncarriers"):
        sample_ids = {sample_id for sample_id, roles in roles_by_sample.items() if cohort_id in roles}
        rows = [row for row in target_rows if row["SAMPLE_ID"] in sample_ids]
        status = _cohort_status(len(rows), parameters["minimum_group_samples"], parameters["minimum_primary_samples"])
        summary_rows.append({
            "SCOPE": "TARGET_CHROMOSOME", "COHORT_ID": cohort_id,
            "SAMPLE_COUNT": str(len(rows)), "EVALUATED_SAMPLE_COUNT": str(len(rows) if target_evaluated else 0),
            "MEDIAN_N_ROH": "", "MEDIAN_TOTAL_ROH_KB": "", "MEDIAN_MAX_ROH_KB": "",
            "TARGET_IN_ROH_COUNT": str(sum(row["TARGET_IN_ROH"] == "true" for row in rows)),
            "COHORT_STATUS": status if target_evaluated else "NOT_EVALUATED",
        })

    segments_path = output_dir / "roh_segments.tsv"
    burden_path = output_dir / "roh_burden.tsv"
    target_path = output_dir / "target_in_roh.tsv"
    cohort_summary_path = output_dir / "roh_cohort_summary.tsv"
    _write_tsv(segments_path, SEGMENT_COLUMNS, all_segments)
    _write_tsv(burden_path, BURDEN_COLUMNS, burden_rows)
    _write_tsv(target_path, TARGET_COLUMNS, target_rows)
    _write_tsv(cohort_summary_path, SUMMARY_COLUMNS, summary_rows)
    validate_tsv_table(segments_path, "roh_segments.schema.json")
    validate_tsv_table(burden_path, "roh_burden.schema.json")
    validate_tsv_table(target_path, "target_in_roh.schema.json")
    validate_tsv_table(cohort_summary_path, "roh_cohort_summary.schema.json")

    target_count = sum(row["TARGET_IN_ROH"] == "true" for row in target_rows)
    analysis_summary_path = output_dir / "roh_analysis_summary.json"
    analysis_summary = {
        "schema_version": "1.0.0", "method_id": "plink19_array_roh_secondary_v1",
        "analysis_role": "SECONDARY_DESCRIPTIVE_AUTOZYGOSITY",
        "scope_statuses": scope_statuses,
        "genomewide_common_variant_count": len(common_variant_ids),
        "target_chromosome_variant_count": target_variant_count,
        "target_variant_id": target_variant_id, "target_in_roh_count": target_count,
        "f_roh_calculated": parameters.get("autosomal_denominator_kb") is not None,
        "feeds_founder_haplotype": False, "feeds_variant_age": False,
        "consumes_local_ld": False,
    }
    validate_json_document(analysis_summary, "roh_analysis_summary.schema.json")
    atomic_write_json(analysis_summary_path, analysis_summary)
    return RohPublication(
        segments_path=segments_path, burden_path=burden_path, target_path=target_path,
        cohort_summary_path=cohort_summary_path,
        analysis_summary_path=analysis_summary_path, native_paths=tuple(native_paths),
        scope_statuses=scope_statuses, segment_count=len(all_segments),
        target_in_roh_count=target_count,
    )
