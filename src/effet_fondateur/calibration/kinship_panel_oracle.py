"""D21 : contrôles indépendants de MAF, partition, couverture et LD de 06.

PLINK 1.9 ne considère que les fondateurs déclarés pour fréquences et LD.
R² = corrélation de dosages au carré, sur appels conjoints ; aucune phase.
https://www.cog-genomics.org/plink/1.9/ld (consulté le 30 octobre 2026).
"""

from __future__ import annotations

from collections import Counter
from pathlib import Path
import statistics
from typing import Any

import numpy as np

from effet_fondateur.calibration.evaluate_qc import compare_metrics
from effet_fondateur.calibration.evaluate_technical import _dataset, raw_table
from effet_fondateur.calibration.materialize_ibd import _json
from effet_fondateur.calibration.qc_oracle import subset
from effet_fondateur.calibration.technical_evaluation import GenotypeMatrix, compare_genotypes


def founder_mask(matrix: GenotypeMatrix, registry: list[dict[str, str]]) -> np.ndarray:
    """Un fondateur PLINK a PID=MID=0 ; aucune ascendance réelle n'est inférée."""
    by_pair = {(row["FID"], row["IID"]): row for row in registry}
    return np.array([by_pair[pair]["PID"] == "0" and by_pair[pair]["MID"] == "0" for pair in matrix.samples])


def allele_frequencies(matrix: GenotypeMatrix, founders: np.ndarray) -> np.ndarray:
    """MAF indépendante sur A/G, manques exclus et dénominateur en chromosomes."""
    calls = matrix.alleles[:, founders]
    if np.any(~np.isin(calls, (-1, 0, 2))):
        raise ValueError("d21_non_fixture_allele")
    called = calls[:, :, 0] >= 0
    copies = 2 * called.sum(axis=1)
    frequencies = np.divide((calls == 2).sum(axis=(1, 2)), copies, out=np.full(len(matrix.variants), np.nan), where=copies > 0)
    return np.minimum(frequencies, 1 - frequencies)


def dosage_r2(first: np.ndarray, second: np.ndarray) -> np.ndarray:
    """R² par ligne sur deux tableaux (paires, fondateurs), -1 = manque."""
    valid = (first >= 0) & (second >= 0)
    first = np.where(valid, first, 0).astype(float)
    second = np.where(valid, second, 0).astype(float)
    counts = valid.sum(axis=1)
    sum_first, sum_second = first.sum(axis=1), second.sum(axis=1)
    numerator = counts * (first * second).sum(axis=1) - sum_first * sum_second
    denominator = (counts * (first * first).sum(axis=1) - sum_first ** 2) * (counts * (second * second).sum(axis=1) - sum_second ** 2)
    return np.clip(np.divide(numerator ** 2, denominator, out=np.full(len(first), np.nan), where=denominator > 0), 0, 1)


def residual_ld(matrix: GenotypeMatrix, founders: np.ndarray, parameters: dict[str, Any]) -> np.ndarray:
    """Fenêtre d'indices stricte (< --ld-window), distance bp inclusive."""
    calls = matrix.alleles[:, founders]
    dosage = np.where(calls[:, :, 0] >= 0, (calls == 2).sum(axis=2), -1)
    chromosomes = np.array([variant[0] for variant in matrix.variants])
    positions = np.array([variant[1] for variant in matrix.variants])
    values = []
    for lag in range(1, min(parameters["residual_ld_window_variants"], len(matrix.variants))):
        eligible = (chromosomes[:-lag] == chromosomes[lag:]) & (positions[lag:] - positions[:-lag] <= 1000 * parameters["residual_ld_window_kb"])
        correlations = dosage_r2(dosage[:-lag][eligible], dosage[lag:][eligible])
        values.extend(correlations[np.isfinite(correlations)].tolist())
    return np.array(values)


def evaluate_panel(preliminary: GenotypeMatrix, registry: list[dict[str, str]], qc_variants: list[dict[str, str]],
                   artifacts: dict[str, Path], parameters: dict[str, Any], config: dict[str, Any]) -> dict[str, Any]:
    """Vérifie le résultat 06 sans imposer un tie-break unique au pruning PLINK."""
    panel = _dataset(artifacts, "kinship_panel")
    retained_ids = {variant[2] for variant in panel.variants}
    if not retained_ids <= {variant[2] for variant in preliminary.variants}:
        return dict(status="FAIL", reasons=["panel_unknown_marker"], metrics={})
    expected_panel = subset(preliminary, np.ones(len(preliminary.samples), bool), np.array([variant[2] in retained_ids for variant in preliminary.variants]))
    genotypes = compare_genotypes(expected_panel, panel, stage="06")
    reasons = list(genotypes["reasons"])
    frequencies = allele_frequencies(preliminary, founder_mask(preliminary, registry))
    old_alerts = {row["VARIANT_ID"]: row["ALERT_CODES"] for row in qc_variants}
    expected_rows, eligible_ids = [], set()
    for index, (chromosome, position, variant_id) in enumerate(preliminary.variants):
        region = next((region["region_id"] for region in parameters["excluded_regions"] if region["chromosome"] == chromosome and region["start_bp"] <= position <= region["end_bp"]), None)
        maf = float(frequencies[index])
        if not np.isfinite(maf):
            status, code = "EXCLUDED_MAF_UNAVAILABLE", "maf_unavailable"
        elif maf < parameters["maf_min"]:
            status, code = "EXCLUDED_LOW_MAF", "low_maf"
        elif region is not None:
            status, code = "EXCLUDED_COMPLEX_REGION", "complex_region:" + region
        else:
            eligible_ids.add(variant_id)
            status, code = ("RETAINED", "") if variant_id in retained_ids else ("EXCLUDED_LD", "ld_pruning")
        expected_rows.append(dict(VARIANT_ID=variant_id, CHROMOSOME=chromosome, POSITION_BP=position, MAF=maf,
                                  PRELIMINARY_ALERT_CODES=old_alerts[variant_id], PANEL_STATUS=status,
                                  EXCLUSION_CODE=code, IN_COMPLEX_REGION=str(status == "EXCLUDED_COMPLEX_REGION").lower()))
    variant_comparison = compare_metrics(expected_rows, raw_table(artifacts["pruned_variants"]), "VARIANT_ID")
    reasons += ["variant_audit_" + key for key, count in variant_comparison.items() if key != "compared" and count]
    if not retained_ids <= eligible_ids:
        reasons.append("ineligible_marker_retained")
    coverage_rows = []
    for chromosome in range(1, 23):
        positions = [variant[1] for variant in panel.variants if variant[0] == chromosome]
        gaps = np.diff(positions).tolist()
        fraction = len(positions) / len(panel.variants)
        if len(positions) < parameters["min_markers_per_autosome"] or fraction > parameters["max_autosome_marker_fraction"]:
            reasons.append("panel_coverage_insufficient")
        coverage_rows.append(dict(CHROMOSOME=str(chromosome), INPUT_VARIANT_COUNT=sum(variant[0] == chromosome for variant in preliminary.variants),
                                  MAF_ELIGIBLE_COUNT=sum(variant[0] == chromosome and variant[2] in eligible_ids for variant in preliminary.variants),
                                  PRUNED_MARKER_COUNT=len(positions), MARKER_FRACTION=float(fraction), SPAN_BP=positions[-1] - positions[0] if positions else "",
                                  MEDIAN_GAP_BP=float(statistics.median(gaps)) if gaps else "", MAX_GAP_BP=max(gaps) if gaps else "", COVERAGE_STATUS="PASS", REASON_CODE=""))
    coverage = compare_metrics(coverage_rows, raw_table(artifacts["panel_coverage"]), "CHROMOSOME")
    reasons += ["coverage_" + key for key, count in coverage.items() if key != "compared" and count]
    if len(panel.variants) < parameters["min_total_markers"]:
        reasons.append("panel_total_insufficient")
    copy_pairs_checked = 0
    for chromosome in range(1, 23):
        ids = [f"kin_chr{chromosome}_{offset:04d}" for offset in config["ld_copy_offsets"]]
        if set(ids) <= eligible_ids:
            copy_pairs_checked += 1
            if set(ids) <= retained_ids:
                reasons.append("perfect_ld_copies_retained_together")
    values = residual_ld(panel, founder_mask(panel, registry), parameters)
    counts = Counter(min(int(value * 10), 9) for value in values)
    bins = [dict(BIN_ID=f"r2_bin_{index+1:02d}", R2_LOWER_BOUND=index / 10., R2_UPPER_BOUND=(index + 1) / 10.,
                 UPPER_INCLUSIVE=str(index == 9).lower(), PAIR_COUNT=counts[index]) for index in range(10)]
    ld_comparison = compare_metrics(bins, raw_table(artifacts["panel_residual_ld_bins"]), "BIN_ID")
    reasons += ["residual_ld_" + key for key, count in ld_comparison.items() if key != "compared" and count]
    report = _json(artifacts["kinship_panel_report"])
    if report["sample_count"] != len(panel.samples) or report["retained_marker_count"] != len(panel.variants) or report["effective_marker_count_operational"] != len(panel.variants) or report["ld_eligible_variant_count"] != len(eligible_ids):
        reasons.append("panel_report_count_mismatch")
    return dict(status="FAIL" if reasons else "PASS", reasons=sorted(set(reasons)), genotype_metrics=genotypes["metrics"],
                variant_comparison=variant_comparison, coverage_comparison=coverage, residual_ld_comparison=ld_comparison,
                metrics=dict(autosome_coverage=sum(bool(row["PRUNED_MARKER_COUNT"]) for row in coverage_rows), retained_markers=len(panel.variants),
                             eligible_markers=len(eligible_ids), low_maf_exclusions=sum(row["PANEL_STATUS"] == "EXCLUDED_LOW_MAF" for row in expected_rows),
                             complex_region_exclusions=sum(row["PANEL_STATUS"] == "EXCLUDED_COMPLEX_REGION" for row in expected_rows),
                             ld_exclusions=len(eligible_ids) - len(retained_ids), perfect_ld_pairs_checked=copy_pairs_checked,
                             residual_ld_pair_count=len(values), residual_ld_max_r2=float(values.max()) if len(values) else None))
