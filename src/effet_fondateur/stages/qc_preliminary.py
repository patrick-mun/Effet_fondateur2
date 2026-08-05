"""Étape 05 : contrôle qualité technique préliminaire du jeu genome-wide."""

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
from collections import defaultdict
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


INDIVIDUAL_COLUMNS = (
    "SAMPLE_ID",
    "FID",
    "IID",
    "ARRAY_BATCH",
    "N_MISS",
    "N_GENO",
    "F_MISS",
    "OBSERVED_HETEROZYGOSITY",
    "EXPECTED_HETEROZYGOSITY",
    "INBREEDING_F",
    "HETEROZYGOSITY_Z",
    "QC_STATUS",
    "EXCLUSION_CODE",
)
VARIANT_COLUMNS = (
    "VARIANT_ID",
    "CHROMOSOME",
    "POSITION_BP",
    "A1",
    "A2",
    "N_MISS",
    "N_GENO",
    "F_MISS",
    "MAF",
    "NCHROBS",
    "BATCH_MIN_F_MISS",
    "BATCH_MAX_F_MISS",
    "BATCH_MISSINGNESS_DELTA",
    "QC_STATUS",
    "EXCLUSION_CODE",
    "ALERT_CODES",
)
BATCH_COLUMNS = (
    "ARRAY_BATCH",
    "SAMPLE_COUNT",
    "MEAN_SAMPLE_F_MISS",
    "MAX_SAMPLE_F_MISS",
    "VARIANT_COUNT",
    "MEAN_VARIANT_F_MISS",
    "MAX_VARIANT_F_MISS",
)
ALERT_COLUMNS = (
    "ALERT_ID",
    "SCOPE",
    "ENTITY_ID",
    "METRIC",
    "OBSERVED_VALUE",
    "THRESHOLD",
    "STATUS",
    "REASON_CODE",
)


class PreliminaryQcInputError(ValueError):
    """Signale une entrée ou un paramètre technique invalide."""


class ExternalToolError(RuntimeError):
    """Signale un échec ou une sortie invalide de PLINK."""


class PreliminaryQcBlockError(RuntimeError):
    """Signale qu'aucun jeu techniquement exploitable ne peut être produit."""


def _build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(prog="qc-preliminary")
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


def _validated_artifact_path(artifact: dict[str, Any], run_dir: Path) -> Path:
    physical_path = _resolve_input_path(artifact, run_dir)
    if not physical_path.is_file() or sha256_file(physical_path) != artifact["sha256"]:
        raise PreliminaryQcInputError("declared_input_modified")
    return physical_path


def _artifact_by_id(stage_inputs: dict[str, Any], artifact_id: str) -> dict[str, Any]:
    matches = [
        artifact
        for artifact in stage_inputs["artifacts"]
        if artifact["artifact_id"] == artifact_id
    ]
    if len(matches) != 1:
        raise PreliminaryQcInputError(f"{artifact_id}_missing_or_ambiguous")
    return matches[0]


def _read_whitespace_table(
    path: Path, *, allow_empty: bool = False
) -> list[dict[str, str]]:
    try:
        lines = [
            line.split()
            for line in path.read_text(encoding="utf-8").splitlines()
            if line
        ]
    except OSError as error:
        raise ExternalToolError(f"missing_plink_report:{path.name}") from error
    if not lines or (len(lines) < 2 and not allow_empty):
        raise ExternalToolError(f"empty_plink_report:{path.name}")
    header = lines[0]
    if any(len(row) != len(header) for row in lines[1:]):
        raise ExternalToolError(f"malformed_plink_report:{path.name}")
    return [dict(zip(header, row, strict=True)) for row in lines[1:]]


def _read_fam(path: Path) -> list[list[str]]:
    rows = [
        line.split()
        for line in path.read_text(encoding="utf-8").splitlines()
        if line
    ]
    if not rows or any(len(row) != 6 for row in rows):
        raise PreliminaryQcInputError("malformed_genomewide_base_fam")
    return rows


def _read_bim(path: Path) -> list[list[str]]:
    rows = [
        line.split()
        for line in path.read_text(encoding="utf-8").splitlines()
        if line
    ]
    if not rows or any(len(row) != 6 for row in rows):
        raise PreliminaryQcInputError("malformed_genomewide_base_bim")
    return rows


def _validate_base_dataset(
    descriptor: dict[str, Any],
    bed_path: Path,
    bim_path: Path,
    fam_path: Path,
    sample_rows: tuple[dict[str, Any], ...],
) -> tuple[list[list[str]], list[list[str]]]:
    validate_json_document(descriptor, "plink_dataset.schema.json")
    if descriptor["scope"] != "autosomal_genomewide":
        raise PreliminaryQcInputError("qc_requires_autosomal_genomewide_dataset")
    if bed_path.read_bytes()[:3] != b"\x6c\x1b\x01":
        raise PreliminaryQcInputError("invalid_genomewide_base_bed_header")
    for suffix, physical_path in {
        "bed": bed_path,
        "bim": bim_path,
        "fam": fam_path,
    }.items():
        declared = descriptor["files"][suffix]
        if declared["name"] != physical_path.name or declared["sha256"] != sha256_file(
            physical_path
        ):
            raise PreliminaryQcInputError("genomewide_base_descriptor_mismatch")
    fam_rows = _read_fam(fam_path)
    bim_rows = _read_bim(bim_path)
    expected_pairs = [(row["FID"], row["IID"]) for row in sample_rows]
    observed_pairs = [(row[0], row[1]) for row in fam_rows]
    if expected_pairs != observed_pairs:
        raise PreliminaryQcInputError("genomewide_base_sample_order_mismatch")
    if (
        descriptor["sample_count"] != len(fam_rows)
        or descriptor["variant_count"] != len(bim_rows)
    ):
        raise PreliminaryQcInputError("genomewide_base_descriptor_count_mismatch")
    if any(
        len({allele for allele in row[4:6] if allele != "0"}) > 2
        for row in bim_rows
    ):
        raise PreliminaryQcInputError("non_biallelic_variant_in_plink_bim")
    return fam_rows, bim_rows


def _float_parameter(
    parameters: dict[str, Any], name: str, default: float, *, maximum: float = 1.0
) -> float:
    value = float(parameters.get(name, default))
    if not 0 <= value <= maximum:
        raise PreliminaryQcInputError(f"invalid_parameter:{name}")
    return value


def _parameters(parameters: dict[str, Any]) -> dict[str, float | int]:
    resolved: dict[str, float | int] = {
        "sample_missingness_max": _float_parameter(
            parameters, "sample_missingness_max", 0.1
        ),
        "variant_missingness_max": _float_parameter(
            parameters, "variant_missingness_max", 0.05
        ),
        "kinship_maf_alert_threshold": _float_parameter(
            parameters, "kinship_maf_alert_threshold", 0.05
        ),
        "heterozygosity_z_alert": _float_parameter(
            parameters, "heterozygosity_z_alert", 3.0, maximum=100.0
        ),
        "batch_missingness_delta_alert": _float_parameter(
            parameters, "batch_missingness_delta_alert", 0.02
        ),
        "duplicate_pi_hat_alert": _float_parameter(
            parameters, "duplicate_pi_hat_alert", 0.98
        ),
    }
    minimum_variants = int(
        parameters.get("duplicate_scan_min_polymorphic_variants", 100)
    )
    timeout_seconds = int(parameters.get("plink_timeout_seconds", 300))
    if minimum_variants <= 0:
        raise PreliminaryQcInputError(
            "invalid_parameter:duplicate_scan_min_polymorphic_variants"
        )
    if timeout_seconds <= 0:
        raise PreliminaryQcInputError("invalid_parameter:plink_timeout_seconds")
    if resolved["heterozygosity_z_alert"] == 0:
        raise PreliminaryQcInputError("invalid_parameter:heterozygosity_z_alert")
    resolved["duplicate_scan_min_polymorphic_variants"] = minimum_variants
    resolved["plink_timeout_seconds"] = timeout_seconds
    return resolved


def _resolve_executable(command: str | None) -> str:
    if command is None:
        raise ExternalToolError("plink_not_configured")
    resolved = shutil.which(command)
    if resolved is None:
        raise ExternalToolError("plink_not_available")
    return resolved


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


def _write_batch_clusters(
    path: Path, sample_rows: tuple[dict[str, Any], ...]
) -> dict[tuple[str, str], str]:
    batches: dict[tuple[str, str], str] = {}
    with path.open("w", encoding="utf-8", newline="") as output_file:
        writer = csv.writer(output_file, delimiter="\t", lineterminator="\n")
        for row in sample_rows:
            batch = row["ARRAY_BATCH"] or "UNREPORTED"
            batches[(row["FID"], row["IID"])] = batch
            writer.writerow((row["FID"], row["IID"], batch))
    return batches


def _nullable_float(value: str) -> float | None:
    try:
        parsed = float(value)
    except ValueError:
        return None
    return parsed if math.isfinite(parsed) else None


def _decimal(value: float | None) -> str:
    if value is None:
        return ""
    normalized = 0.0 if abs(value) < 5e-13 else value
    return f"{normalized:.8f}".rstrip("0").rstrip(".")


def _new_alert(
    alerts: list[dict[str, str]],
    scope: str,
    entity_id: str,
    metric: str,
    observed_value: str,
    threshold: str,
    status: str,
    reason_code: str,
) -> None:
    alerts.append(
        {
            "ALERT_ID": f"qc_alert_{len(alerts) + 1:06d}",
            "SCOPE": scope,
            "ENTITY_ID": entity_id,
            "METRIC": metric,
            "OBSERVED_VALUE": observed_value,
            "THRESHOLD": threshold,
            "STATUS": status,
            "REASON_CODE": reason_code,
        }
    )


def _individual_metrics(
    imiss_rows: list[dict[str, str]],
    het_rows: list[dict[str, str]],
    sample_rows: tuple[dict[str, Any], ...],
    batches: dict[tuple[str, str], str],
    sample_missingness_max: float,
    heterozygosity_z_alert: float,
    alerts: list[dict[str, str]],
) -> list[dict[str, str]]:
    imiss_by_id = {(row["FID"], row["IID"]): row for row in imiss_rows}
    het_by_id = {(row["FID"], row["IID"]): row for row in het_rows}
    records: list[dict[str, Any]] = []
    for sample in sample_rows:
        plink_id = (sample["FID"], sample["IID"])
        if plink_id not in imiss_by_id or plink_id not in het_by_id:
            raise ExternalToolError("missing_individual_qc_metric")
        imiss = imiss_by_id[plink_id]
        het = het_by_id[plink_id]
        nonmissing_count = int(het["N(NM)"])
        observed_homozygosity = float(het["O(HOM)"])
        expected_homozygosity = float(het["E(HOM)"])
        observed_heterozygosity = (
            (nonmissing_count - observed_homozygosity) / nonmissing_count
            if nonmissing_count
            else None
        )
        expected_heterozygosity = (
            (nonmissing_count - expected_homozygosity) / nonmissing_count
            if nonmissing_count
            else None
        )
        records.append(
            {
                "sample": sample,
                "imiss": imiss,
                "batch": batches[plink_id],
                "observed_heterozygosity": observed_heterozygosity,
                "expected_heterozygosity": expected_heterozygosity,
                "inbreeding_f": _nullable_float(het["F"]),
                "heterozygosity_z": None,
            }
        )
    heterozygosities = [
        record["observed_heterozygosity"]
        for record in records
        if record["observed_heterozygosity"] is not None
    ]
    if len(heterozygosities) >= 3:
        standard_deviation = statistics.pstdev(heterozygosities)
        if standard_deviation > 0:
            mean_heterozygosity = statistics.mean(heterozygosities)
            for record in records:
                value = record["observed_heterozygosity"]
                if value is not None:
                    record["heterozygosity_z"] = (
                        value - mean_heterozygosity
                    ) / standard_deviation

    rows: list[dict[str, str]] = []
    for record in records:
        sample = record["sample"]
        imiss = record["imiss"]
        missingness = int(imiss["N_MISS"]) / int(imiss["N_GENO"])
        exclusion_code = ""
        status = "PASS"
        if missingness > sample_missingness_max:
            exclusion_code = "sample_missingness"
            status = "EXCLUDED"
            _new_alert(
                alerts,
                "SAMPLE",
                sample["SAMPLE_ID"],
                "F_MISS",
                _decimal(missingness),
                _decimal(sample_missingness_max),
                "EXCLUDED",
                exclusion_code,
            )
        heterozygosity_z = record["heterozygosity_z"]
        if (
            status != "EXCLUDED"
            and heterozygosity_z is not None
            and abs(heterozygosity_z) > heterozygosity_z_alert
        ):
            status = "ALERT"
            _new_alert(
                alerts,
                "SAMPLE",
                sample["SAMPLE_ID"],
                "HETEROZYGOSITY_Z",
                _decimal(heterozygosity_z),
                _decimal(heterozygosity_z_alert),
                "ALERT",
                "heterozygosity_outlier",
            )
        rows.append(
            {
                "SAMPLE_ID": sample["SAMPLE_ID"],
                "FID": sample["FID"],
                "IID": sample["IID"],
                "ARRAY_BATCH": record["batch"],
                "N_MISS": imiss["N_MISS"],
                "N_GENO": imiss["N_GENO"],
                "F_MISS": _decimal(missingness),
                "OBSERVED_HETEROZYGOSITY": _decimal(
                    record["observed_heterozygosity"]
                ),
                "EXPECTED_HETEROZYGOSITY": _decimal(
                    record["expected_heterozygosity"]
                ),
                "INBREEDING_F": _decimal(record["inbreeding_f"]),
                "HETEROZYGOSITY_Z": _decimal(heterozygosity_z),
                "QC_STATUS": status,
                "EXCLUSION_CODE": exclusion_code,
            }
        )
    return rows


def _batch_variant_missingness(
    rows: list[dict[str, str]],
) -> dict[str, dict[str, float]]:
    by_variant: dict[str, dict[str, float]] = defaultdict(dict)
    for row in rows:
        by_variant[row["SNP"]][row["CLST"]] = int(row["N_MISS"]) / int(
            row["N_GENO"]
        )
    return dict(by_variant)


def _variant_metrics(
    lmiss_rows: list[dict[str, str]],
    frequency_rows: list[dict[str, str]],
    batch_missingness: dict[str, dict[str, float]],
    bim_rows: list[list[str]],
    variant_missingness_max: float,
    maf_alert_threshold: float,
    batch_delta_alert: float,
    alerts: list[dict[str, str]],
) -> list[dict[str, str]]:
    lmiss_by_id = {row["SNP"]: row for row in lmiss_rows}
    frequency_by_id = {row["SNP"]: row for row in frequency_rows}
    bim_by_id = {row[1]: row for row in bim_rows}
    if set(lmiss_by_id) != set(bim_by_id) or set(frequency_by_id) != set(bim_by_id):
        raise ExternalToolError("variant_qc_metric_set_mismatch")
    output_rows: list[dict[str, str]] = []
    for variant_id in [row[1] for row in bim_rows]:
        lmiss = lmiss_by_id[variant_id]
        frequency = frequency_by_id[variant_id]
        bim = bim_by_id[variant_id]
        missingness = int(lmiss["N_MISS"]) / int(lmiss["N_GENO"])
        maf = _nullable_float(frequency["MAF"])
        batch_values = list(batch_missingness.get(variant_id, {}).values())
        batch_min = min(batch_values) if batch_values else None
        batch_max = max(batch_values) if batch_values else None
        batch_delta = (
            batch_max - batch_min
            if batch_min is not None and batch_max is not None
            else None
        )
        alert_codes: list[str] = []
        exclusion_code = ""
        status = "PASS"
        if missingness > variant_missingness_max:
            exclusion_code = "variant_missingness"
            status = "EXCLUDED"
            _new_alert(
                alerts,
                "VARIANT",
                variant_id,
                "F_MISS",
                _decimal(missingness),
                _decimal(variant_missingness_max),
                "EXCLUDED",
                exclusion_code,
            )
        if maf is not None and maf < maf_alert_threshold:
            alert_codes.append("low_maf_for_kinship")
            if status == "PASS":
                status = "ALERT"
            _new_alert(
                alerts,
                "VARIANT",
                variant_id,
                "MAF",
                _decimal(maf),
                _decimal(maf_alert_threshold),
                "ALERT",
                "low_maf_for_kinship",
            )
        if batch_delta is not None and batch_delta > batch_delta_alert:
            alert_codes.append("batch_missingness_delta")
            if status == "PASS":
                status = "ALERT"
            _new_alert(
                alerts,
                "VARIANT",
                variant_id,
                "BATCH_MISSINGNESS_DELTA",
                _decimal(batch_delta),
                _decimal(batch_delta_alert),
                "ALERT",
                "batch_missingness_delta",
            )
        output_rows.append(
            {
                "VARIANT_ID": variant_id,
                "CHROMOSOME": bim[0],
                "POSITION_BP": bim[3],
                "A1": frequency["A1"],
                "A2": frequency["A2"],
                "N_MISS": lmiss["N_MISS"],
                "N_GENO": lmiss["N_GENO"],
                "F_MISS": _decimal(missingness),
                "MAF": _decimal(maf),
                "NCHROBS": frequency["NCHROBS"],
                "BATCH_MIN_F_MISS": _decimal(batch_min),
                "BATCH_MAX_F_MISS": _decimal(batch_max),
                "BATCH_MISSINGNESS_DELTA": _decimal(batch_delta),
                "QC_STATUS": status,
                "EXCLUSION_CODE": exclusion_code,
                "ALERT_CODES": ",".join(alert_codes),
            }
        )
    return output_rows


def _batch_summary(
    individual_rows: list[dict[str, str]],
    batch_missingness: dict[str, dict[str, float]],
) -> list[dict[str, str]]:
    samples_by_batch: dict[str, list[float]] = defaultdict(list)
    for row in individual_rows:
        samples_by_batch[row["ARRAY_BATCH"]].append(float(row["F_MISS"]))
    variants_by_batch: dict[str, list[float]] = defaultdict(list)
    for values_by_batch in batch_missingness.values():
        for batch, missingness in values_by_batch.items():
            variants_by_batch[batch].append(missingness)
    rows = []
    for batch in sorted(samples_by_batch):
        sample_values = samples_by_batch[batch]
        variant_values = variants_by_batch.get(batch, [])
        if not variant_values:
            raise ExternalToolError("missing_batch_variant_metrics")
        rows.append(
            {
                "ARRAY_BATCH": batch,
                "SAMPLE_COUNT": str(len(sample_values)),
                "MEAN_SAMPLE_F_MISS": _decimal(statistics.mean(sample_values)),
                "MAX_SAMPLE_F_MISS": _decimal(max(sample_values)),
                "VARIANT_COUNT": str(len(variant_values)),
                "MEAN_VARIANT_F_MISS": _decimal(statistics.mean(variant_values)),
                "MAX_VARIANT_F_MISS": _decimal(max(variant_values)),
            }
        )
    return rows


def _duplicate_alerts(
    plink_executable: str,
    base_prefix: Path,
    output_dir: Path,
    frequency_rows: list[dict[str, str]],
    minimum_polymorphic_variants: int,
    pi_hat_threshold: float,
    timeout_seconds: int,
    alerts: list[dict[str, str]],
) -> tuple[str, int, int]:
    polymorphic_count = sum(
        (maf := _nullable_float(row["MAF"])) is not None and maf > 0
        for row in frequency_rows
    )
    if polymorphic_count < minimum_polymorphic_variants:
        _new_alert(
            alerts,
            "GLOBAL",
            "genomewide_pre_qc",
            "DUPLICATE_SCAN",
            str(polymorphic_count),
            str(minimum_polymorphic_variants),
            "NOT_EVALUATED",
            "insufficient_polymorphic_variants",
        )
        return "NOT_EVALUATED", polymorphic_count, 0
    duplicate_prefix = output_dir / "preliminary_duplicates"
    _run_plink(
        plink_executable,
        [
            "--bfile",
            str(base_prefix),
            "--genome",
            "--min",
            str(pi_hat_threshold),
            "--out",
            str(duplicate_prefix),
        ],
        timeout_seconds,
    )
    pairs = _read_whitespace_table(
        duplicate_prefix.with_suffix(".genome"), allow_empty=True
    )
    for row in pairs:
        entity_id = f"{row['FID1']}:{row['IID1']}|{row['FID2']}:{row['IID2']}"
        _new_alert(
            alerts,
            "PAIR",
            entity_id,
            "PI_HAT",
            row["PI_HAT"],
            _decimal(pi_hat_threshold),
            "ALERT",
            "potential_technical_duplicate",
        )
    return "EVALUATED", polymorphic_count, len(pairs)


def _write_tsv(
    path: Path, columns: tuple[str, ...], rows: list[dict[str, str]]
) -> None:
    with path.open("w", encoding="utf-8", newline="") as output_file:
        writer = csv.DictWriter(
            output_file, fieldnames=columns, delimiter="\t", lineterminator="\n"
        )
        writer.writeheader()
        writer.writerows({column: row[column] for column in columns} for row in rows)


def _validate_filtered_dataset(
    prefix: Path,
    expected_sample_ids: list[tuple[str, str]],
    expected_variant_ids: list[str],
) -> tuple[list[list[str]], list[list[str]]]:
    if prefix.with_suffix(".bed").read_bytes()[:3] != b"\x6c\x1b\x01":
        raise ExternalToolError("invalid_pre_qc_bed_header")
    fam_rows = _read_fam(prefix.with_suffix(".fam"))
    bim_rows = _read_bim(prefix.with_suffix(".bim"))
    if [(row[0], row[1]) for row in fam_rows] != expected_sample_ids:
        raise ExternalToolError("unexpected_pre_qc_sample_set")
    if [row[1] for row in bim_rows] != expected_variant_ids:
        raise ExternalToolError("unexpected_pre_qc_variant_set")
    return fam_rows, bim_rows


def _set_id(scope: str, values: list[str]) -> str:
    digest = hashlib.sha256(
        "".join(f"{value}\n" for value in values).encode()
    ).hexdigest()
    return f"{scope}_{digest[:12]}"


def _write_dataset_descriptor(
    path: Path,
    prefix: Path,
    base_descriptor: dict[str, Any],
    sample_set_id: str,
    variant_set_id: str,
    sample_count: int,
    variant_count: int,
) -> None:
    descriptor = {
        **base_descriptor,
        "dataset_id": "genomewide_pre_qc",
        "sample_set_id": sample_set_id,
        "variant_set_id": variant_set_id,
        "sample_count": sample_count,
        "variant_count": variant_count,
        "source_format": "PLINK_PRELIMINARY_QC",
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


def _dataset_artifacts(
    prefix: Path,
    descriptor_path: Path,
    published_dir: PurePosixPath,
    producer_stage: str,
    signature: str,
    assembly: str,
    sample_set_id: str,
    variant_set_id: str,
) -> list[dict[str, Any]]:
    artifacts = []
    media_types = {
        "bed": "application/octet-stream",
        "bim": "text/plain",
        "fam": "text/plain",
    }
    for suffix in ("bed", "bim", "fam"):
        physical_path = prefix.with_suffix(f".{suffix}")
        artifacts.append(
            build_file_artifact(
                physical_path=physical_path,
                published_path=str(published_dir / physical_path.name),
                artifact_id=f"genomewide_pre_qc_{suffix}",
                artifact_type=f"plink_{suffix}",
                media_type=media_types[suffix],
                producer_stage=producer_stage,
                producer_signature=signature,
                assembly=assembly,
                sample_set_id=sample_set_id,
                variant_set_id=variant_set_id,
                sensitivity="sensitive_genetic",
            )
        )
    artifacts.append(
        build_file_artifact(
            physical_path=descriptor_path,
            published_path=str(published_dir / descriptor_path.name),
            artifact_id="genomewide_pre_qc_dataset",
            artifact_type="plink_dataset_descriptor",
            media_type="application/json",
            producer_stage=producer_stage,
            producer_signature=signature,
            schema_name="plink_dataset.schema.json",
            schema_version="1.0.0",
            assembly=assembly,
            sample_set_id=sample_set_id,
            variant_set_id=variant_set_id,
            sensitivity="sensitive_genetic",
        )
    )
    return artifacts


def execute(stage_inputs_path: Path, output_dir: Path) -> int:
    """Produit un jeu genome-wide pré-QC et des métriques techniques auditables."""
    started_at = utc_now()
    started_clock = monotonic()
    stage_inputs = read_json(stage_inputs_path)
    validate_json_document(stage_inputs, "stage_inputs.schema.json")
    run_dir = output_dir.parent.parent
    config = load_pipeline_config(run_dir / "config.resolved.yaml")
    resolved_parameters = _parameters(stage_inputs["parameters"])
    timeout_seconds = int(resolved_parameters["plink_timeout_seconds"])

    artifact_ids = (
        "samples_master",
        "genomewide_base_bed",
        "genomewide_base_bim",
        "genomewide_base_fam",
        "genomewide_base_dataset",
    )
    input_artifacts = {
        artifact_id: _artifact_by_id(stage_inputs, artifact_id)
        for artifact_id in artifact_ids
    }
    paths = {
        artifact_id: _validated_artifact_path(artifact, run_dir)
        for artifact_id, artifact in input_artifacts.items()
    }
    samples_table = validate_tsv_table(
        paths["samples_master"], "samples_master.schema.json"
    )
    base_descriptor = read_json(paths["genomewide_base_dataset"])
    _, base_bim_rows = _validate_base_dataset(
        base_descriptor,
        paths["genomewide_base_bed"],
        paths["genomewide_base_bim"],
        paths["genomewide_base_fam"],
        samples_table.rows,
    )
    if input_artifacts["samples_master"]["sample_set_id"] != base_descriptor[
        "sample_set_id"
    ]:
        raise PreliminaryQcInputError("genomewide_base_sample_set_id_mismatch")

    output_dir.mkdir(parents=True, exist_ok=True)
    batch_path = output_dir / "sample_batches.tsv"
    batches = _write_batch_clusters(batch_path, samples_table.rows)
    plink_executable = _resolve_executable(config["tools"]["plink"])
    plink_version = _plink_version(plink_executable)
    base_prefix = paths["genomewide_base_bed"].with_suffix("")
    metrics_prefix = output_dir / "preliminary_metrics"
    _run_plink(
        plink_executable,
        [
            "--bfile",
            str(base_prefix),
            "--missing",
            "--freq",
            "--het",
            "small-sample",
            "--out",
            str(metrics_prefix),
        ],
        timeout_seconds,
    )
    batch_prefix = output_dir / "preliminary_batch_missingness"
    _run_plink(
        plink_executable,
        [
            "--bfile",
            str(base_prefix),
            "--missing",
            "--within",
            str(batch_path),
            "--out",
            str(batch_prefix),
        ],
        timeout_seconds,
    )
    imiss_rows = _read_whitespace_table(metrics_prefix.with_suffix(".imiss"))
    lmiss_rows = _read_whitespace_table(metrics_prefix.with_suffix(".lmiss"))
    frequency_rows = _read_whitespace_table(metrics_prefix.with_suffix(".frq"))
    het_rows = _read_whitespace_table(metrics_prefix.with_suffix(".het"))
    batch_lmiss_rows = _read_whitespace_table(batch_prefix.with_suffix(".lmiss"))

    alerts: list[dict[str, str]] = []
    individual_rows = _individual_metrics(
        imiss_rows,
        het_rows,
        samples_table.rows,
        batches,
        float(resolved_parameters["sample_missingness_max"]),
        float(resolved_parameters["heterozygosity_z_alert"]),
        alerts,
    )
    batch_missingness = _batch_variant_missingness(batch_lmiss_rows)
    variant_rows = _variant_metrics(
        lmiss_rows,
        frequency_rows,
        batch_missingness,
        base_bim_rows,
        float(resolved_parameters["variant_missingness_max"]),
        float(resolved_parameters["kinship_maf_alert_threshold"]),
        float(resolved_parameters["batch_missingness_delta_alert"]),
        alerts,
    )
    batch_rows = _batch_summary(individual_rows, batch_missingness)
    if len(batch_rows) < 2:
        _new_alert(
            alerts,
            "GLOBAL",
            "genomewide_pre_qc",
            "BATCH_COMPARISON",
            str(len(batch_rows)),
            "2",
            "NOT_EVALUATED",
            "single_or_unreported_batch",
        )
    duplicate_status, polymorphic_count, duplicate_pair_count = _duplicate_alerts(
        plink_executable,
        base_prefix,
        output_dir,
        frequency_rows,
        int(resolved_parameters["duplicate_scan_min_polymorphic_variants"]),
        float(resolved_parameters["duplicate_pi_hat_alert"]),
        timeout_seconds,
        alerts,
    )
    _new_alert(
        alerts,
        "GLOBAL",
        "genomewide_pre_qc",
        "SEX_DISCORDANCE",
        "",
        "",
        "NOT_EVALUATED",
        "autosomal_dataset_without_sex_chromosomes",
    )
    _new_alert(
        alerts,
        "GLOBAL",
        "genomewide_pre_qc",
        "HWE_FILTER",
        "not_applied",
        "",
        "INFO",
        "cohorts_not_frozen",
    )

    retained_sample_pairs = [
        (row["FID"], row["IID"])
        for row in individual_rows
        if row["QC_STATUS"] != "EXCLUDED"
    ]
    retained_variant_ids = [
        row["VARIANT_ID"]
        for row in variant_rows
        if row["QC_STATUS"] != "EXCLUDED"
    ]
    if not retained_sample_pairs:
        raise PreliminaryQcBlockError("all_samples_excluded_by_missingness")
    if not retained_variant_ids:
        raise PreliminaryQcBlockError("all_variants_excluded_by_missingness")
    pre_qc_prefix = output_dir / "genomewide_pre_qc"
    excluded_sample_path = output_dir / "pre_qc_excluded_samples.tsv"
    excluded_variant_path = output_dir / "pre_qc_excluded_variants.txt"
    excluded_sample_rows = [
        (row["FID"], row["IID"])
        for row in individual_rows
        if row["QC_STATUS"] == "EXCLUDED"
    ]
    excluded_variant_ids = [
        row["VARIANT_ID"]
        for row in variant_rows
        if row["QC_STATUS"] == "EXCLUDED"
    ]
    filter_arguments = ["--bfile", str(base_prefix)]
    if excluded_sample_rows:
        with excluded_sample_path.open("w", encoding="utf-8", newline="") as file:
            writer = csv.writer(file, delimiter="\t", lineterminator="\n")
            writer.writerows(excluded_sample_rows)
        filter_arguments.extend(("--remove", str(excluded_sample_path)))
    if excluded_variant_ids:
        excluded_variant_path.write_text(
            "".join(f"{variant_id}\n" for variant_id in excluded_variant_ids),
            encoding="utf-8",
        )
        filter_arguments.extend(("--exclude", str(excluded_variant_path)))
    filter_arguments.extend(("--make-bed", "--out", str(pre_qc_prefix)))
    _run_plink(
        plink_executable,
        filter_arguments,
        timeout_seconds,
    )
    filtered_fam_rows, filtered_bim_rows = _validate_filtered_dataset(
        pre_qc_prefix, retained_sample_pairs, retained_variant_ids
    )

    individual_path = output_dir / "qc_individual_metrics.tsv"
    variant_path = output_dir / "qc_variant_metrics.tsv"
    batch_summary_path = output_dir / "qc_batch_summary.tsv"
    alerts_path = output_dir / "qc_alerts.tsv"
    _write_tsv(individual_path, INDIVIDUAL_COLUMNS, individual_rows)
    _write_tsv(variant_path, VARIANT_COLUMNS, variant_rows)
    _write_tsv(batch_summary_path, BATCH_COLUMNS, batch_rows)
    _write_tsv(alerts_path, ALERT_COLUMNS, alerts)
    validate_tsv_table(individual_path, "qc_individual_metrics.schema.json")
    validate_tsv_table(variant_path, "qc_variant_metrics.schema.json")
    validate_tsv_table(batch_summary_path, "qc_batch_summary.schema.json")
    validate_tsv_table(alerts_path, "qc_alerts.schema.json")

    sample_set_id = _set_id(
        "preqc_samples", [f"{row[0]}:{row[1]}" for row in filtered_fam_rows]
    )
    variant_set_id = _set_id(
        "preqc_variants", [row[1] for row in filtered_bim_rows]
    )
    descriptor_path = output_dir / "genomewide_pre_qc.dataset.json"
    _write_dataset_descriptor(
        descriptor_path,
        pre_qc_prefix,
        base_descriptor,
        sample_set_id,
        variant_set_id,
        len(filtered_fam_rows),
        len(filtered_bim_rows),
    )
    report_path = output_dir / "qc_preliminary_report.json"
    report = {
        "schema_version": "1.0.0",
        "assembly": base_descriptor["assembly"],
        "parameters": resolved_parameters,
        "sample_count_before": len(individual_rows),
        "sample_count_after": len(filtered_fam_rows),
        "variant_count_before": len(variant_rows),
        "variant_count_after": len(filtered_bim_rows),
        "excluded_sample_count": sum(
            row["QC_STATUS"] == "EXCLUDED" for row in individual_rows
        ),
        "excluded_variant_count": sum(
            row["QC_STATUS"] == "EXCLUDED" for row in variant_rows
        ),
        "low_maf_alert_count": sum(
            "low_maf_for_kinship" in row["ALERT_CODES"] for row in variant_rows
        ),
        "batch_delta_alert_count": sum(
            "batch_missingness_delta" in row["ALERT_CODES"] for row in variant_rows
        ),
        "heterozygosity_alert_count": sum(
            row["QC_STATUS"] == "ALERT" for row in individual_rows
        ),
        "duplicate_scan_status": duplicate_status,
        "polymorphic_variant_count": polymorphic_count,
        "potential_duplicate_pair_count": duplicate_pair_count,
        "non_biallelic_variant_count": 0,
        "sex_check_status": "NOT_EVALUATED",
        "hwe_filter_applied": False,
    }
    atomic_write_json(report_path, report)

    temporary_paths = [batch_path, excluded_sample_path, excluded_variant_path]
    temporary_prefixes = (
        metrics_prefix,
        batch_prefix,
        pre_qc_prefix,
        output_dir / "preliminary_duplicates",
    )
    for prefix in temporary_prefixes:
        for suffix in (
            "imiss",
            "lmiss",
            "frq",
            "het",
            "genome",
            "log",
            "nosex",
        ):
            temporary_paths.append(prefix.with_suffix(f".{suffix}"))
    for path in temporary_paths:
        if path.exists():
            path.unlink()

    producer_stage = f"{stage_inputs['stage_id']}_{stage_inputs['stage_name']}"
    published_dir = PurePosixPath(stage_inputs["published_output_dir"])
    output_artifacts = _dataset_artifacts(
        pre_qc_prefix,
        descriptor_path,
        published_dir,
        producer_stage,
        stage_inputs["signature"],
        base_descriptor["assembly"],
        sample_set_id,
        variant_set_id,
    )
    for artifact_id, artifact_type, path, schema_name in (
        (
            "qc_individual_metrics",
            "qc_individual_metrics",
            individual_path,
            "qc_individual_metrics.schema.json",
        ),
        (
            "qc_variant_metrics",
            "qc_variant_metrics",
            variant_path,
            "qc_variant_metrics.schema.json",
        ),
        (
            "qc_batch_summary",
            "qc_batch_summary",
            batch_summary_path,
            "qc_batch_summary.schema.json",
        ),
        ("qc_alerts", "qc_alerts", alerts_path, "qc_alerts.schema.json"),
        ("qc_preliminary_report", "qc_preliminary_report", report_path, None),
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
                assembly=base_descriptor["assembly"],
                sample_set_id=sample_set_id,
                variant_set_id=variant_set_id,
                sensitivity="sensitive_genetic",
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
        "method_id": "plink_preliminary_genomewide_qc_v1",
        "signature": stage_inputs["signature"],
        "started_at": started_at,
        "completed_at": utc_now(),
        "duration_seconds": monotonic() - started_clock,
        "inputs": list(input_artifacts.values()),
        "outputs": output_artifacts,
        "parameters": resolved_parameters,
        "tools": [
            {
                "tool": "plink",
                "configured": config["tools"]["plink"],
                "version": plink_version,
            }
        ],
        "counts": {
            "samples_before": len(individual_rows),
            "samples_after": len(filtered_fam_rows),
            "variants_before": len(variant_rows),
            "variants_after": len(filtered_bim_rows),
            "alerts": len(alerts),
            "potential_duplicate_pairs": duplicate_pair_count,
        },
        "metrics": {},
        "exclusions": [
            {
                "scope": "sample",
                "count": report["excluded_sample_count"],
                "reason": "sample_missingness",
            },
            {
                "scope": "variant",
                "count": report["excluded_variant_count"],
                "reason": "variant_missingness",
            },
        ],
        "warnings": [],
        "checks": [
            {"check": "plink_input_integrity", "status": "PASS"},
            {"check": "technical_missingness_filters", "status": "PASS"},
            {"check": "hwe_filter_absent", "status": "PASS"},
            {"check": "sex_discordance", "status": "NOT_EVALUATED"},
            {"check": "duplicate_scan", "status": duplicate_status},
        ],
        "known_limits": [
            "L'hétérozygotie est une alerte non exclusive avant résolution de "
            "la structure de population.",
            "Le contrôle du sexe est impossible sur un jeu limité aux autosomes.",
            "Le scan de duplicats est préliminaire et l'étape 07 porte la "
            "décision d'apparentement.",
            "Aucun filtre HWE n'est appliqué avant le gel des cohortes "
            "indépendantes.",
        ],
        "expected_visualizations": [
            "sample_missingness",
            "heterozygosity",
            "batch_missingness",
            "maf_distribution",
            "aggregated_missingness_matrix",
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
    except PreliminaryQcBlockError as error:
        sys.stderr.write(f"{error}\n")
        return 4
    except (
        PreliminaryQcInputError,
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
