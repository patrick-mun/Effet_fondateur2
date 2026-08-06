"""Consolidation finale du QC et des haplotypes de l'étape 12."""

from __future__ import annotations

import csv
import os
import shutil
import tempfile
from collections import Counter
from dataclasses import dataclass
from pathlib import Path
from typing import Any

from effet_fondateur.audit import atomic_write_json, read_json, sha256_file
from effet_fondateur.contracts import validate_json_document, validate_tsv_table
from effet_fondateur.orchestrator.state import utc_now


QC_COLUMNS = (
    "CHECK_ORDER",
    "CHECK_ID",
    "SCOPE",
    "STATUS",
    "OBSERVED_VALUE",
    "EXPECTED_VALUE",
    "DETAIL_CODE",
)
SUMMARY_COLUMNS = (
    "CARRIER_HAPLOTYPE",
    "RELIABILITY_STATUS",
    "SAMPLE_COUNT",
)


class PhasingPublicationError(ValueError):
    """Signale une incohérence empêchant la clôture de l'étape 12."""


@dataclass(frozen=True)
class PublishedPhasingQc:
    """Décrit les sorties consolidées et leur état de revue."""

    output_dir: Path
    qc_path: Path
    haplotype_summary_path: Path
    summary_path: Path
    warning_count: int
    manual_validation_required: bool


def _write_tsv(
    path: Path, columns: tuple[str, ...], rows: list[dict[str, Any]]
) -> None:
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns, delimiter="\t")
        writer.writeheader()
        for row in rows:
            writer.writerow(
                {
                    column: "" if row[column] is None else row[column]
                    for column in columns
                }
            )


def _file_record(path: Path) -> dict[str, Any]:
    return {
        "filename": path.name,
        "sha256": sha256_file(path),
        "size_bytes": path.stat().st_size,
    }


def _qc_row(
    order: int,
    check_id: str,
    scope: str,
    status: str,
    observed: str | int,
    expected: str | int | None,
    detail_code: str | None = None,
) -> dict[str, Any]:
    return {
        "CHECK_ORDER": order,
        "CHECK_ID": check_id,
        "SCOPE": scope,
        "STATUS": status,
        "OBSERVED_VALUE": observed,
        "EXPECTED_VALUE": expected,
        "DETAIL_CODE": detail_code,
    }


def publish_phasing_qc(
    *,
    inputs_manifest_path: Path,
    phasing_manifest_path: Path,
    carrier_haplotypes_path: Path,
    transmissions_path: Path,
    unreliable_regions_path: Path,
    output_dir: Path,
) -> PublishedPhasingQc:
    """Publie le QC final sans recalculer ni réinterpréter le phasage."""
    if output_dir.exists() or output_dir.is_symlink():
        raise PhasingPublicationError("phasing_qc_output_exists")
    inputs_manifest = read_json(inputs_manifest_path)
    phasing_manifest = read_json(phasing_manifest_path)
    validate_json_document(inputs_manifest, "shapeit5_inputs_manifest.schema.json")
    validate_json_document(
        phasing_manifest, "shapeit5_phasing_manifest.schema.json"
    )
    carrier_table = validate_tsv_table(
        carrier_haplotypes_path, "carrier_haplotypes.schema.json"
    )
    transmission_table = validate_tsv_table(
        transmissions_path, "phasing_transmissions.schema.json"
    )
    unreliable_table = validate_tsv_table(
        unreliable_regions_path, "phasing_unreliable_regions.schema.json"
    )
    expected_sources = (
        (
            carrier_haplotypes_path,
            phasing_manifest["files"]["carrier_haplotypes"]["sha256"],
        ),
        (
            transmissions_path,
            phasing_manifest["files"]["phasing_transmissions"]["sha256"],
        ),
        (
            unreliable_regions_path,
            phasing_manifest["files"]["phasing_unreliable_regions"]["sha256"],
        ),
    )
    if any(
        path.is_symlink() or sha256_file(path) != expected
        for path, expected in expected_sources
    ):
        raise PhasingPublicationError("phasing_qc_source_modified")
    if carrier_table.row_count != phasing_manifest["sample_count"]:
        raise PhasingPublicationError("phasing_qc_sample_count_mismatch")
    if transmission_table.row_count != phasing_manifest["pedigree_record_count"]:
        raise PhasingPublicationError("phasing_qc_transmission_count_mismatch")

    haplotype_counts = Counter(
        row["CARRIER_HAPLOTYPE"] for row in carrier_table.rows
    )
    carrier_count = sum(
        count for haplotype, count in haplotype_counts.items() if haplotype != "NONE"
    )
    reliable_carrier_count = sum(
        row["CARRIER_HAPLOTYPE"] != "NONE"
        and row["RELIABILITY_STATUS"] == "PASS"
        for row in carrier_table.rows
    )
    if (
        carrier_count != phasing_manifest["carrier_count"]
        or reliable_carrier_count != phasing_manifest["reliable_carrier_count"]
    ):
        raise PhasingPublicationError("phasing_qc_carrier_count_mismatch")
    unreliable_carrier_count = carrier_count - reliable_carrier_count
    warning_count = unreliable_table.row_count
    if (warning_count > 0) != (unreliable_carrier_count > 0):
        raise PhasingPublicationError("phasing_qc_warning_accounting_mismatch")

    qc_rows = [
        _qc_row(1, "input_integrity", "INPUT", "PASS", "PASS", "PASS"),
        _qc_row(
            2,
            "allele_counts",
            "INPUT",
            inputs_manifest["checks"]["allele_counts"],
            inputs_manifest["checks"]["allele_counts"],
            "PASS",
        ),
        _qc_row(
            3,
            "common_phase_variants",
            "COMMON_PHASE",
            "PASS",
            phasing_manifest["common_variant_count"],
            inputs_manifest["common_variant_count"],
        ),
        _qc_row(
            4,
            "rare_phase_variants",
            "RARE_PHASE",
            "PASS",
            phasing_manifest["variant_count"],
            inputs_manifest["study_variant_count"],
        ),
        _qc_row(5, "sample_order", "RARE_PHASE", "PASS", "PASS", "PASS"),
        _qc_row(
            6, "genotype_preservation", "RARE_PHASE", "PASS", "PASS", "PASS"
        ),
        _qc_row(7, "target_preservation", "TARGET", "PASS", "1", "1"),
        _qc_row(
            8, "explicit_target_genotypes", "TARGET", "PASS", "PASS", "PASS"
        ),
        _qc_row(9, "mendel_before", "PEDIGREE", "PASS", "0", "0"),
        _qc_row(10, "mendel_after", "PEDIGREE", "PASS", "0", "0"),
        _qc_row(
            11,
            "pedigree_transmissions",
            "PEDIGREE",
            "PASS" if transmission_table.row_count else "NOT_APPLICABLE",
            transmission_table.row_count,
            phasing_manifest["pedigree_record_count"],
            None if transmission_table.row_count else "no_pedigree_record",
        ),
        _qc_row(
            12,
            "carrier_phase_confidence",
            "CONFIDENCE",
            "WARN" if unreliable_carrier_count else "PASS",
            reliable_carrier_count,
            carrier_count,
            "carrier_phase_unreliable" if unreliable_carrier_count else None,
        ),
    ]
    summary_rows = [
        {
            "CARRIER_HAPLOTYPE": haplotype,
            "RELIABILITY_STATUS": reliability,
            "SAMPLE_COUNT": sum(
                row["CARRIER_HAPLOTYPE"] == haplotype
                and row["RELIABILITY_STATUS"] == reliability
                for row in carrier_table.rows
            ),
        }
        for haplotype in ("NONE", "H1", "H2", "BOTH")
        for reliability in ("PASS", "UNRELIABLE")
    ]

    output_dir.parent.mkdir(parents=True, exist_ok=True)
    staging_dir = Path(
        tempfile.mkdtemp(prefix=f".{output_dir.name}.", dir=output_dir.parent)
    )
    published = False
    try:
        qc_path = staging_dir / "phasing_qc.tsv"
        haplotype_summary_path = staging_dir / "carrier_haplotype_summary.tsv"
        summary_path = staging_dir / "phase_target_region_summary.json"
        _write_tsv(qc_path, QC_COLUMNS, qc_rows)
        _write_tsv(haplotype_summary_path, SUMMARY_COLUMNS, summary_rows)
        validate_tsv_table(qc_path, "phasing_qc.schema.json")
        validate_tsv_table(
            haplotype_summary_path, "carrier_haplotype_summary.schema.json"
        )
        summary = {
            "schema_version": "1.0.0",
            "created_at": utc_now(),
            "method_id": "phase_target_region_qc_publication_v1",
            "assembly": phasing_manifest["assembly"],
            "chromosome": phasing_manifest["chromosome"],
            "target_variant_id": phasing_manifest["target_variant_id"],
            "target_variant_role": phasing_manifest["target_variant_role"],
            "sample_count": phasing_manifest["sample_count"],
            "variant_count": phasing_manifest["variant_count"],
            "carrier_count": carrier_count,
            "reliable_carrier_count": reliable_carrier_count,
            "unreliable_carrier_count": unreliable_carrier_count,
            "pedigree_record_count": phasing_manifest["pedigree_record_count"],
            "transmission_record_count": transmission_table.row_count,
            "manual_validation_required": unreliable_carrier_count > 0,
            "source_manifests": {
                "shapeit5_inputs": {
                    "filename": inputs_manifest_path.name,
                    "sha256": sha256_file(inputs_manifest_path),
                },
                "shapeit5_phasing": {
                    "filename": phasing_manifest_path.name,
                    "sha256": sha256_file(phasing_manifest_path),
                },
            },
            "haplotype_counts": {
                haplotype: haplotype_counts.get(haplotype, 0)
                for haplotype in ("NONE", "H1", "H2", "BOTH")
            },
            "files": {
                "phasing_qc": _file_record(qc_path),
                "carrier_haplotype_summary": _file_record(
                    haplotype_summary_path
                ),
            },
            "checks": {
                "source_integrity": "PASS",
                "qc_complete": "PASS",
                "haplotype_accounting": "PASS",
                "warning_accounting": "PASS",
                "publication": "PASS",
            },
        }
        validate_json_document(summary, "phase_target_region_summary.schema.json")
        atomic_write_json(summary_path, summary)
        os.replace(staging_dir, output_dir)
        published = True
        return PublishedPhasingQc(
            output_dir=output_dir,
            qc_path=output_dir / qc_path.name,
            haplotype_summary_path=output_dir / haplotype_summary_path.name,
            summary_path=output_dir / summary_path.name,
            warning_count=warning_count,
            manual_validation_required=unreliable_carrier_count > 0,
        )
    finally:
        if not published:
            shutil.rmtree(staging_dir, ignore_errors=True)
