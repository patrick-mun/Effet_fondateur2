"""Détection déterministe de segments IBS centrés sur un variant cible."""

from __future__ import annotations

import csv
import json
import subprocess
from dataclasses import dataclass
from pathlib import Path
from typing import Any

from effet_fondateur.audit import atomic_write_json
from effet_fondateur.contracts import validate_tsv_table


class FounderAnalysisError(ValueError):
    """Signale une entrée incompatible avec l'analyse IBS locale."""


@dataclass(frozen=True)
class Variant:
    """Variant ordonné et positionné sur la carte génétique."""

    variant_id: str
    position_bp: int
    position_cm: float
    is_target: bool
    genotypes: tuple[str, ...]


@dataclass(frozen=True)
class CarrierHaplotype:
    """Chromosome porteur fiable représentant une unité indépendante."""

    independent_unit_id: str
    sample_id: str
    family_id: str
    haplotype_index: int
    haplotype_id: str
    phase_confidence: str


@dataclass(frozen=True)
class Segment:
    """Intervalle IBS inclusif autour de la cible."""

    left_index: int
    target_index: int
    right_index: int


@dataclass(frozen=True)
class FounderAnalysisResult:
    """Chemins et métriques publiés par l'analyse."""

    segments_path: Path
    consensus_path: Path
    sharing_path: Path
    discordances_path: Path
    summary_path: Path
    status: str
    selected_carrier_count: int
    background_haplotype_count: int
    matching_background_haplotype_count: int


SEGMENT_COLUMNS = (
    "INDEPENDENT_UNIT_ID", "SAMPLE_ID", "FAMILY_ID", "CARRIER_HAPLOTYPE_ID",
    "TARGET_VARIANT_ID", "LEFT_BOUND_BP", "TARGET_BP", "RIGHT_BOUND_BP",
    "LEFT_LENGTH_CM", "RIGHT_LENGTH_CM", "PHASING_CONFIDENCE",
    "SEGMENT_METHOD", "SEGMENT_STATUS", "EXCLUSION_CODE",
)
CONSENSUS_COLUMNS = (
    "ANALYSIS_ID", "TARGET_VARIANT_ID", "INTERPRETATION", "CONSENSUS_STATUS",
    "CARRIER_UNIT_COUNT", "LEFT_BOUND_BP", "TARGET_BP", "RIGHT_BOUND_BP",
    "LEFT_LENGTH_CM", "RIGHT_LENGTH_CM", "LEFT_MARKER_COUNT",
    "RIGHT_MARKER_COUNT", "SIGNATURE_MARKER_COUNT", "BACKGROUND_COHORT",
    "BACKGROUND_HAPLOTYPE_COUNT", "MATCHING_BACKGROUND_HAPLOTYPE_COUNT",
    "BACKGROUND_MATCH_FREQUENCY", "SEGMENT_METHOD", "MAP_VERSION",
)
SHARING_COLUMNS = (
    "INDEPENDENT_UNIT_ID_1", "SAMPLE_ID_1", "CARRIER_HAPLOTYPE_ID_1",
    "INDEPENDENT_UNIT_ID_2", "SAMPLE_ID_2", "CARRIER_HAPLOTYPE_ID_2",
    "LEFT_BOUND_BP", "TARGET_BP", "RIGHT_BOUND_BP", "LEFT_LENGTH_CM",
    "RIGHT_LENGTH_CM", "LEFT_MARKER_COUNT", "RIGHT_MARKER_COUNT", "PAIR_STATUS",
)
DISCORDANCE_COLUMNS = (
    "VARIANT_ID", "POSITION_BP", "POSITION_CM", "SIDE", "CARRIER_ALLELES",
    "CALLED_CARRIER_COUNT", "DISTINCT_ALLELE_COUNT", "MARKER_STATUS",
)


def _write_tsv(path: Path, columns: tuple[str, ...], rows: list[dict[str, Any]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def _run_bcftools(command: str, arguments: list[str], timeout_seconds: float) -> str:
    try:
        completed = subprocess.run(
            [command, *arguments], check=False, capture_output=True, text=True,
            timeout=timeout_seconds,
        )
    except (OSError, subprocess.TimeoutExpired) as error:
        raise FounderAnalysisError("bcftools_query_unavailable") from error
    if completed.returncode != 0:
        raise FounderAnalysisError(f"bcftools_query_failed:{completed.returncode}")
    return completed.stdout


def _load_variants(
    bcf_path: Path, map_path: Path, bcftools_command: str, timeout_seconds: float
) -> tuple[tuple[str, ...], tuple[Variant, ...], str]:
    map_table = validate_tsv_table(map_path, "target_genetic_map.schema.json")
    map_rows = {row["VARIANT_ID"]: row for row in map_table.rows}
    target_rows = [row for row in map_table.rows if row["IS_TARGET_VARIANT"]]
    if len(target_rows) != 1:
        raise FounderAnalysisError("target_map_variant_missing_or_ambiguous")
    sample_ids = tuple(
        line for line in _run_bcftools(
            bcftools_command, ["query", "--list-samples", str(bcf_path)], timeout_seconds
        ).splitlines() if line
    )
    if not sample_ids or len(sample_ids) != len(set(sample_ids)):
        raise FounderAnalysisError("phased_sample_ids_invalid")
    query = _run_bcftools(
        bcftools_command,
        ["query", "--format", "%ID\\t%POS[\\t%GT]\\n", str(bcf_path)],
        timeout_seconds,
    )
    variants: list[Variant] = []
    observed_ids: set[str] = set()
    for line in query.splitlines():
        fields = line.split("\t")
        if len(fields) != len(sample_ids) + 2:
            raise FounderAnalysisError("phased_variant_column_count_mismatch")
        variant_id, position_text, *genotypes = fields
        if variant_id not in map_rows or variant_id in observed_ids:
            raise FounderAnalysisError("phased_variant_map_identity_mismatch")
        map_row = map_rows[variant_id]
        if int(position_text) != int(map_row["POSITION_BP"]):
            raise FounderAnalysisError("phased_variant_position_mismatch")
        variants.append(
            Variant(
                variant_id=variant_id,
                position_bp=int(position_text),
                position_cm=float(map_row["POSITION_CM"]),
                is_target=bool(map_row["IS_TARGET_VARIANT"]),
                genotypes=tuple(genotypes),
            )
        )
        observed_ids.add(variant_id)
    if observed_ids != set(map_rows):
        raise FounderAnalysisError("phased_variant_map_set_mismatch")
    variants.sort(key=lambda variant: (variant.position_bp, variant.variant_id))
    if sum(variant.is_target for variant in variants) != 1:
        raise FounderAnalysisError("target_phased_variant_missing_or_ambiguous")
    return sample_ids, tuple(variants), map_table.sha256


def _allele(genotype: str, haplotype_index: int) -> str | None:
    if "|" not in genotype:
        return None
    alleles = genotype.split("|")
    if len(alleles) != 2 or any(value not in {"0", "1", "."} for value in alleles):
        return None
    value = alleles[haplotype_index]
    return None if value == "." else value


def _shared_segment(
    variants: tuple[Variant, ...], haplotypes: tuple[tuple[int, int], ...]
) -> Segment:
    target_index = next(index for index, variant in enumerate(variants) if variant.is_target)

    def shared(index: int) -> bool:
        alleles = {
            _allele(variants[index].genotypes[sample_index], haplotype_index)
            for sample_index, haplotype_index in haplotypes
        }
        return None not in alleles and len(alleles) == 1

    if not shared(target_index):
        return Segment(target_index, target_index, target_index)
    left_index = target_index
    while left_index > 0 and shared(left_index - 1):
        left_index -= 1
    right_index = target_index
    while right_index + 1 < len(variants) and shared(right_index + 1):
        right_index += 1
    return Segment(left_index, target_index, right_index)


def _decimal(value: float) -> str:
    return f"{value:.12g}"


def _segment_values(variants: tuple[Variant, ...], segment: Segment) -> dict[str, str]:
    left = variants[segment.left_index]
    target = variants[segment.target_index]
    right = variants[segment.right_index]
    return {
        "LEFT_BOUND_BP": str(left.position_bp),
        "TARGET_BP": str(target.position_bp),
        "RIGHT_BOUND_BP": str(right.position_bp),
        "LEFT_LENGTH_CM": _decimal(target.position_cm - left.position_cm),
        "RIGHT_LENGTH_CM": _decimal(right.position_cm - target.position_cm),
        "LEFT_MARKER_COUNT": str(segment.target_index - segment.left_index),
        "RIGHT_MARKER_COUNT": str(segment.right_index - segment.target_index),
    }


def infer_target_centered_ibs(
    *, phased_bcf_path: Path, carrier_haplotypes_path: Path,
    cohorts_path: Path, samples_master_path: Path, genetic_map_path: Path,
    output_dir: Path, bcftools_command: str, timeout_seconds: float,
    minimum_independent_carriers: int, minimum_flank_markers: int,
) -> FounderAnalysisResult:
    """Publie une analyse IBS exacte, centrée sur la mutation et sans prétention IBD."""
    if minimum_independent_carriers < 2 or minimum_flank_markers < 1:
        raise FounderAnalysisError("invalid_founder_analysis_threshold")
    sample_ids, variants, map_version = _load_variants(
        phased_bcf_path, genetic_map_path, bcftools_command, timeout_seconds
    )
    sample_indexes = {sample_id: index for index, sample_id in enumerate(sample_ids)}
    carrier_rows = validate_tsv_table(
        carrier_haplotypes_path, "carrier_haplotypes.schema.json"
    ).rows
    cohort_rows = validate_tsv_table(cohorts_path, "cohorts_frozen.schema.json").rows
    sample_rows = validate_tsv_table(samples_master_path, "samples_master.schema.json").rows
    families = {row["SAMPLE_ID"]: row["FID"] for row in sample_rows}
    assignments = {row["SAMPLE_ID"]: row for row in carrier_rows}
    independent_rows = [
        row for row in cohort_rows
        if row["COHORT_ID"] == "target_carriers_independent" and row["INCLUDED"]
    ]
    selected: list[CarrierHaplotype] = []
    excluded_rows: list[dict[str, str]] = []
    target_variant = next(variant for variant in variants if variant.is_target)
    for cohort_row in independent_rows:
        sample_id = cohort_row["SAMPLE_ID"]
        assignment = assignments.get(sample_id)
        exclusion_code = None
        if sample_id not in sample_indexes or sample_id not in families or assignment is None:
            raise FounderAnalysisError("carrier_sample_identity_mismatch")
        if assignment["RELIABILITY_STATUS"] != "PASS":
            exclusion_code = "PHASING_UNRELIABLE"
        elif assignment["CARRIER_HAPLOTYPE"] == "BOTH":
            exclusion_code = "TARGET_HOMOZYGOUS_TWO_CHROMOSOMES"
        elif assignment["CARRIER_HAPLOTYPE"] not in {"H1", "H2"}:
            exclusion_code = "CARRIER_HAPLOTYPE_NOT_ASSIGNED"
        if exclusion_code:
            excluded_rows.append({
                "INDEPENDENT_UNIT_ID": cohort_row["INDEPENDENT_UNIT_ID"],
                "SAMPLE_ID": sample_id, "FAMILY_ID": families[sample_id],
                "CARRIER_HAPLOTYPE_ID": assignment["CARRIER_HAPLOTYPE"],
                "TARGET_VARIANT_ID": target_variant.variant_id,
                "LEFT_BOUND_BP": str(target_variant.position_bp),
                "TARGET_BP": str(target_variant.position_bp),
                "RIGHT_BOUND_BP": str(target_variant.position_bp),
                "LEFT_LENGTH_CM": "0", "RIGHT_LENGTH_CM": "0",
                "PHASING_CONFIDENCE": assignment["PHASE_CONFIDENCE"] or "",
                "SEGMENT_METHOD": "target_centered_exact_ibs_v1",
                "SEGMENT_STATUS": "EXCLUDED", "EXCLUSION_CODE": exclusion_code,
            })
            continue
        haplotype_index = 0 if assignment["CARRIER_HAPLOTYPE"] == "H1" else 1
        if _allele(target_variant.genotypes[sample_indexes[sample_id]], haplotype_index) != "1":
            raise FounderAnalysisError("target_carrier_haplotype_discordant")
        selected.append(CarrierHaplotype(
            independent_unit_id=cohort_row["INDEPENDENT_UNIT_ID"], sample_id=sample_id,
            family_id=families[sample_id], haplotype_index=haplotype_index,
            haplotype_id=assignment["CARRIER_HAPLOTYPE"],
            phase_confidence=assignment["PHASE_CONFIDENCE"] or "",
        ))
    selected.sort(key=lambda carrier: (carrier.independent_unit_id, carrier.sample_id))
    if not selected:
        raise FounderAnalysisError("no_reliable_independent_carrier_haplotype")
    carrier_keys = tuple(
        (sample_indexes[carrier.sample_id], carrier.haplotype_index) for carrier in selected
    )
    segment = _shared_segment(variants, carrier_keys)
    segment_values = _segment_values(variants, segment)
    enough_carriers = len(selected) >= minimum_independent_carriers
    enough_flanks = (
        int(segment_values["LEFT_MARKER_COUNT"]) >= minimum_flank_markers
        and int(segment_values["RIGHT_MARKER_COUNT"]) >= minimum_flank_markers
    )
    status = "SUPPORTED_IBS_CANDIDATE" if enough_carriers and enough_flanks else (
        "INSUFFICIENT_CARRIERS" if not enough_carriers else "INSUFFICIENT_INFORMATIVE_MARKERS"
    )
    segment_rows = excluded_rows + [
        {
            "INDEPENDENT_UNIT_ID": carrier.independent_unit_id,
            "SAMPLE_ID": carrier.sample_id, "FAMILY_ID": carrier.family_id,
            "CARRIER_HAPLOTYPE_ID": carrier.haplotype_id,
            "TARGET_VARIANT_ID": target_variant.variant_id,
            "LEFT_BOUND_BP": segment_values["LEFT_BOUND_BP"],
            "TARGET_BP": segment_values["TARGET_BP"],
            "RIGHT_BOUND_BP": segment_values["RIGHT_BOUND_BP"],
            "LEFT_LENGTH_CM": segment_values["LEFT_LENGTH_CM"],
            "RIGHT_LENGTH_CM": segment_values["RIGHT_LENGTH_CM"],
            "PHASING_CONFIDENCE": carrier.phase_confidence,
            "SEGMENT_METHOD": "target_centered_exact_ibs_v1",
            "SEGMENT_STATUS": "INCLUDED" if status == "SUPPORTED_IBS_CANDIDATE" else "INSUFFICIENT",
            "EXCLUSION_CODE": "" if status == "SUPPORTED_IBS_CANDIDATE" else status,
        }
        for carrier in selected
    ]
    sharing_rows: list[dict[str, str]] = []
    for first_index, first in enumerate(selected):
        for second in selected[first_index:]:
            pair_segment = _shared_segment(variants, (
                (sample_indexes[first.sample_id], first.haplotype_index),
                (sample_indexes[second.sample_id], second.haplotype_index),
            ))
            values = _segment_values(variants, pair_segment)
            sharing_rows.append({
                "INDEPENDENT_UNIT_ID_1": first.independent_unit_id,
                "SAMPLE_ID_1": first.sample_id, "CARRIER_HAPLOTYPE_ID_1": first.haplotype_id,
                "INDEPENDENT_UNIT_ID_2": second.independent_unit_id,
                "SAMPLE_ID_2": second.sample_id, "CARRIER_HAPLOTYPE_ID_2": second.haplotype_id,
                **values,
                "PAIR_STATUS": "SELF" if first == second else (
                    "SHARED" if int(values["LEFT_MARKER_COUNT"]) >= minimum_flank_markers and int(values["RIGHT_MARKER_COUNT"]) >= minimum_flank_markers else "INSUFFICIENT"
                ),
            })
    discordance_rows: list[dict[str, str]] = []
    for index, variant in enumerate(variants):
        alleles = [_allele(variant.genotypes[sample_index], haplotype_index) for sample_index, haplotype_index in carrier_keys]
        called = [allele for allele in alleles if allele is not None]
        side = "TARGET" if variant.is_target else ("LEFT" if index < segment.target_index else "RIGHT")
        marker_status = "MISSING" if len(called) != len(alleles) else ("CONCORDANT" if len(set(called)) <= 1 else "DISCORDANT")
        discordance_rows.append({
            "VARIANT_ID": variant.variant_id, "POSITION_BP": str(variant.position_bp),
            "POSITION_CM": _decimal(variant.position_cm), "SIDE": side,
            "CARRIER_ALLELES": ",".join(allele or "." for allele in alleles),
            "CALLED_CARRIER_COUNT": str(len(called)),
            "DISTINCT_ALLELE_COUNT": str(len(set(called))), "MARKER_STATUS": marker_status,
        })
    signature_indexes = [index for index in range(segment.left_index, segment.right_index + 1) if index != segment.target_index]
    background_ids = {
        row["SAMPLE_ID"] for row in cohort_rows
        if row["COHORT_ID"] == "controls_unrelated" and row["INCLUDED"]
    }
    background_ids.update(
        row["SAMPLE_ID"] for row in cohort_rows
        if row["COHORT_ID"] == "target_chromosome_all_qc" and row["INCLUDED"] and row["ROLE"] != "CARRIER"
    )
    background_haplotypes = 0
    matching_background = 0
    if selected and signature_indexes:
        reference_sample_index, reference_haplotype_index = carrier_keys[0]
        signature = tuple(_allele(variants[index].genotypes[reference_sample_index], reference_haplotype_index) for index in signature_indexes)
        for sample_id in sorted(background_ids):
            if sample_id not in sample_indexes:
                raise FounderAnalysisError("background_sample_identity_mismatch")
            for haplotype_index in (0, 1):
                observed = tuple(_allele(variants[index].genotypes[sample_indexes[sample_id]], haplotype_index) for index in signature_indexes)
                if None not in observed:
                    background_haplotypes += 1
                    if observed == signature:
                        matching_background += 1
    frequency = matching_background / background_haplotypes if background_haplotypes else 0.0
    consensus_rows = [{
        "ANALYSIS_ID": "primary_exact_ibs", "TARGET_VARIANT_ID": target_variant.variant_id,
        "INTERPRETATION": "IBS_SHARED_CANDIDATE" if status == "SUPPORTED_IBS_CANDIDATE" else "NO_FOUNDER_CONCLUSION",
        "CONSENSUS_STATUS": status, "CARRIER_UNIT_COUNT": str(len(selected)),
        "LEFT_BOUND_BP": segment_values["LEFT_BOUND_BP"], "TARGET_BP": segment_values["TARGET_BP"],
        "RIGHT_BOUND_BP": segment_values["RIGHT_BOUND_BP"], "LEFT_LENGTH_CM": segment_values["LEFT_LENGTH_CM"],
        "RIGHT_LENGTH_CM": segment_values["RIGHT_LENGTH_CM"], "LEFT_MARKER_COUNT": segment_values["LEFT_MARKER_COUNT"],
        "RIGHT_MARKER_COUNT": segment_values["RIGHT_MARKER_COUNT"], "SIGNATURE_MARKER_COUNT": str(len(signature_indexes)),
        "BACKGROUND_COHORT": "controls_unrelated_plus_noncarriers", "BACKGROUND_HAPLOTYPE_COUNT": str(background_haplotypes),
        "MATCHING_BACKGROUND_HAPLOTYPE_COUNT": str(matching_background), "BACKGROUND_MATCH_FREQUENCY": _decimal(frequency),
        "SEGMENT_METHOD": "target_centered_exact_ibs_v1", "MAP_VERSION": map_version,
    }]
    output_dir.mkdir(parents=True, exist_ok=True)
    paths = {
        "segments": output_dir / "founder_segments.tsv",
        "consensus": output_dir / "founder_consensus.tsv",
        "sharing": output_dir / "founder_sharing_matrix.tsv",
        "discordances": output_dir / "founder_discordances.tsv",
        "summary": output_dir / "founder_analysis_summary.json",
    }
    _write_tsv(paths["segments"], SEGMENT_COLUMNS, segment_rows)
    _write_tsv(paths["consensus"], CONSENSUS_COLUMNS, consensus_rows)
    _write_tsv(paths["sharing"], SHARING_COLUMNS, sharing_rows)
    _write_tsv(paths["discordances"], DISCORDANCE_COLUMNS, discordance_rows)
    summary = {
        "schema_version": "1.0.0", "method_id": "target_centered_exact_ibs_v1",
        "interpretation": consensus_rows[0]["INTERPRETATION"], "status": status,
        "selected_carrier_count": len(selected), "excluded_carrier_count": len(excluded_rows),
        "background_haplotype_count": background_haplotypes,
        "matching_background_haplotype_count": matching_background,
        "minimum_independent_carriers": minimum_independent_carriers,
        "minimum_flank_markers": minimum_flank_markers,
        "ibd_claimed": False,
    }
    atomic_write_json(paths["summary"], summary)
    return FounderAnalysisResult(
        segments_path=paths["segments"], consensus_path=paths["consensus"],
        sharing_path=paths["sharing"], discordances_path=paths["discordances"],
        summary_path=paths["summary"], status=status,
        selected_carrier_count=len(selected), background_haplotype_count=background_haplotypes,
        matching_background_haplotype_count=matching_background,
    )
