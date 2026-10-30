"""Oracle D20 indépendant des décisions et métriques de production.

Entrées : génotypes nucléotidiques (-1/-1 manquant), lots et politiques
préspecifiées. Les fractions sont calculées à partir des comptes entiers.
HWE : distribution conditionnelle exacte des hétérozygotes, sans mid-p,
énumérée en log-espace (Wigginton et al., 2005 ; PLINK 1.9 --hardy).
Référence primaire : https://pmc.ncbi.nlm.nih.gov/articles/PMC1199378/.
Les estimations PI_HAT ne sont pas réimplémentées : la paire clonée connue
est évaluée séparément comme test de détection, sans oracle d'ascendance.
"""

from __future__ import annotations

from dataclasses import dataclass
import math
from typing import Any

import numpy as np

from effet_fondateur.calibration.technical_evaluation import GenotypeMatrix


COHORTS = ("controls_unrelated", "target_carriers_independent", "target_chromosome_all_qc")


@dataclass(frozen=True)
class QcExpectation:
    """Décisions, métriques et matrice retenue ; aucune observation du pipeline."""

    samples: list[dict[str, Any]]
    variants: list[dict[str, Any]]
    batches: list[dict[str, Any]]
    retained: GenotypeMatrix


def exact_hwe(heterozygotes: int, homozygotes_a: int, homozygotes_b: int) -> float | None:
    """P-value bilatérale exacte : somme des états au plus aussi probables.

    Effectifs diploïdes non négatifs ; aucun individu appelé => None. Les
    états équiprobables sont inclus avec une tolérance numérique de 1e-12.
    La formule conditionne sur le nombre de copies de l'allèle rare, sans
    approximation chi-deux ni recours à PLINK ou à un helper de production.
    """
    counts = (heterozygotes, homozygotes_a, homozygotes_b)
    if any(type(count) is not int or count < 0 for count in counts):
        raise ValueError("d20_hwe_counts_invalid")
    total = sum(counts)
    if not total:
        return None
    rare = 2 * min(homozygotes_a, homozygotes_b) + heterozygotes
    states = range(rare % 2, rare + 1, 2)
    weights = {hets: hets * math.log(2) - math.lgamma(hets + 1)
               - math.lgamma((rare - hets) // 2 + 1)
               - math.lgamma(total - (rare + hets) // 2 + 1) for hets in states}
    maximum = max(weights.values())
    observed = weights[heterozygotes]
    denominator = sum(math.exp(weight - maximum) for weight in weights.values())
    return min(1.0, sum(math.exp(weight - maximum) for weight in weights.values()
                        if weight <= observed + 1e-12) / denominator)


def subset(matrix: GenotypeMatrix, sample_mask: np.ndarray, variant_mask: np.ndarray) -> GenotypeMatrix:
    """Sous-ensemble conservant l'ordre des axes et les allèles physiques."""
    return GenotypeMatrix(tuple(key for key, keep in zip(matrix.samples, sample_mask, strict=True) if keep),
                          tuple(key for key, keep in zip(matrix.variants, variant_mask, strict=True) if keep),
                          matrix.alleles[variant_mask][:, sample_mask].copy())


def _metrics(matrix: GenotypeMatrix) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """Manques par cellule, hétérozygotie par sujet, MAF et HWE par locus.

    Toutes les unités sont des proportions ou p-values ; une mesure sans
    appel est NaN en mémoire, puis explicitement absente dans la comparaison.
    """
    missing = matrix.alleles[:, :, 0] < 0
    called = ~missing
    heterozygous = called & (matrix.alleles[:, :, 0] != matrix.alleles[:, :, 1])
    called_counts = called.sum(axis=0)
    heterozygosity = np.divide(heterozygous.sum(axis=0), called_counts,
                              out=np.full(len(matrix.samples), np.nan), where=called_counts > 0)
    maf = np.full(len(matrix.variants), np.nan)
    hwe = np.full(len(matrix.variants), np.nan)
    for index, calls in enumerate(matrix.alleles):
        observed = calls[called[index]]
        if not len(observed):
            continue
        bases, counts = np.unique(observed, return_counts=True)
        if len(bases) > 2:
            raise ValueError("d20_non_biallelic_matrix")
        maf[index] = counts.min() / counts.sum() if len(bases) == 2 else 0.0
        het_count = int(heterozygous[index].sum())
        first_hom = int(((observed[:, 0] == bases[0]) & (observed[:, 1] == bases[0])).sum())
        hwe[index] = exact_hwe(het_count, first_hom, len(observed) - first_hom - het_count)
    return missing, heterozygosity, maf, hwe


def preliminary_expectation(matrix: GenotypeMatrix, batches: list[str], parameters: dict[str, Any]) -> QcExpectation:
    """QC 05 : exclusions simultanées sur le jeu initial ; autres défauts alertés."""
    missing, heterozygosity, maf, _ = _metrics(matrix)
    finite = heterozygosity[np.isfinite(heterozygosity)]
    std = float(finite.std()) if len(finite) >= 3 else 0.0
    z_scores = (heterozygosity - finite.mean()) / std if std > 0 else np.full(len(matrix.samples), np.nan)
    called = ~missing
    chromosome_counts = 2 * called.sum(axis=1)
    # Nei, correction small-sample explicite de PLINK --het : n copies,
    # H_e = 2p(1-p)n/(n-1). F = (O_hom-E_hom)/(N-E_hom).
    expected_het_by_variant = np.divide(2 * maf * (1 - maf) * chromosome_counts, chromosome_counts - 1,
                                      out=np.zeros(len(matrix.variants)), where=chromosome_counts > 1)
    expected_hom = called.T @ (1 - expected_het_by_variant)
    called_by_sample = called.sum(axis=0)
    expected_fraction = np.divide(called_by_sample - expected_hom, called_by_sample,
                                  out=np.full(len(matrix.samples), np.nan), where=called_by_sample > 0)
    observed_hom = (called & (matrix.alleles[:, :, 0] == matrix.alleles[:, :, 1])).sum(axis=0)
    denominator = called_by_sample - expected_hom
    inbreeding = np.divide(observed_hom - expected_hom, denominator,
                           out=np.full(len(matrix.samples), np.nan), where=denominator > 0)
    sample_rows = []
    for index, (_, sample) in enumerate(matrix.samples):
        fraction = float(missing[:, index].mean())
        excluded = fraction > parameters["sample_missingness_max"]
        het_alert = not excluded and abs(z_scores[index]) > parameters["heterozygosity_z_alert"]
        sample_rows.append(dict(IID=sample, FID=matrix.samples[index][0], SAMPLE_ID=sample, ARRAY_BATCH=batches[index],
                                N_MISS=int(missing[:, index].sum()), N_GENO=len(matrix.variants), F_MISS=fraction,
                                OBSERVED_HETEROZYGOSITY=float(heterozygosity[index]), HETEROZYGOSITY_Z=float(z_scores[index]),
                                EXPECTED_HETEROZYGOSITY=float(expected_fraction[index]), INBREEDING_F=float(inbreeding[index]),
                                QC_STATUS="EXCLUDED" if excluded else "ALERT" if het_alert else "PASS",
                                EXCLUSION_CODE="sample_missingness" if excluded else ""))
    variant_rows = []
    batch_array = np.array(batches)
    for index, (chromosome, position_bp, variant) in enumerate(matrix.variants):
        fraction = float(missing[index].mean())
        batch_values = [float(missing[index, batch_array == batch].mean()) for batch in sorted(set(batches))]
        delta = max(batch_values) - min(batch_values)
        codes = []
        if maf[index] < parameters["kinship_maf_alert_threshold"]:
            codes.append("low_maf_for_kinship")
        if delta > parameters["batch_missingness_delta_alert"]:
            codes.append("batch_missingness_delta")
        excluded = fraction > parameters["variant_missingness_max"]
        variant_rows.append(dict(VARIANT_ID=variant, CHROMOSOME=chromosome, POSITION_BP=position_bp,
                                 N_MISS=int(missing[index].sum()), N_GENO=len(matrix.samples),
                                 F_MISS=fraction, MAF=float(maf[index]), NCHROBS=2 * int((~missing[index]).sum()),
                                 BATCH_MIN_F_MISS=min(batch_values), BATCH_MAX_F_MISS=max(batch_values), BATCH_MISSINGNESS_DELTA=delta,
                                 QC_STATUS="EXCLUDED" if excluded else "ALERT" if codes else "PASS",
                                 EXCLUSION_CODE="variant_missingness" if excluded else "", ALERT_CODES=",".join(codes)))
    sample_mask = np.array([row["QC_STATUS"] != "EXCLUDED" for row in sample_rows])
    variant_mask = np.array([row["QC_STATUS"] != "EXCLUDED" for row in variant_rows])
    return QcExpectation(sample_rows, variant_rows, [], subset(matrix, sample_mask, variant_mask))


def final_expectation(matrix: GenotypeMatrix, batches: list[str], cohort: str,
                      parameters: dict[str, Any], target_id: str) -> QcExpectation:
    """QC 10 : sous-cohorte déjà formée, retrait sujets puis métriques variants."""
    if cohort not in COHORTS:
        raise ValueError("d20_unknown_cohort")
    missing = matrix.alleles[:, :, 0] < 0
    sample_rows = []
    for index, (_, sample) in enumerate(matrix.samples):
        fraction = float(missing[:, index].mean())
        excluded = fraction > parameters["sample_missingness_max"]
        sample_rows.append(dict(IID=sample, FID=matrix.samples[index][0], SAMPLE_ID=sample, ARRAY_BATCH=batches[index],
                                N_MISS=int(missing[:, index].sum()), N_GENO=len(matrix.variants), F_MISS=fraction,
                                QC_STATUS="EXCLUDED" if excluded else "PASS", EXCLUSION_CODE="sample_missingness" if excluded else "",
                                PRESENT_AFTER_QC=str(not excluded).lower()))
    means = {batch: float(np.mean([row["F_MISS"] for row, label in zip(sample_rows, batches, strict=True) if label == batch]))
             for batch in sorted(set(batches))}
    batch_rows = [dict(ARRAY_BATCH=batch, SAMPLE_COUNT=batches.count(batch), MEAN_F_MISS=mean,
                       BATCH_DELTA_FROM_BEST=mean - min(means.values()),
                       QC_STATUS="ALERT" if mean - min(means.values()) > parameters["batch_missingness_delta_alert"] else "PASS",
                       ALERT_CODE="batch_missingness_delta" if mean - min(means.values()) > parameters["batch_missingness_delta_alert"] else "")
                  for batch, mean in means.items()]
    sample_mask = np.array([row["QC_STATUS"] != "EXCLUDED" for row in sample_rows])
    retained_samples = subset(matrix, sample_mask, np.ones(len(matrix.variants), dtype=bool))
    missing, _, maf, hwe = _metrics(retained_samples)
    control = cohort == COHORTS[0]
    local = cohort == COHORTS[2]
    miss_key = "controls_variant_missingness_max" if control else "target_region_variant_missingness_max" if local else "carriers_variant_missingness_max"
    variant_rows = []
    for index, (chromosome, position_bp, variant) in enumerate(matrix.variants):
        fraction = float(missing[index].mean())
        codes = ["variant_missingness"] if fraction > parameters[miss_key] else []
        if (control or local) and (not np.isfinite(maf[index]) or maf[index] < parameters["controls_maf_min" if control else "target_region_maf_min"]):
            codes.append("low_maf_controls" if control else "low_maf_target_region")
        if control:
            if not np.isfinite(hwe[index]):
                codes.append("hwe_not_evaluable")
            elif hwe[index] < parameters["controls_hwe_p_min"]:
                codes.append("hwe_controls")
        forced = variant == target_id and local and bool(codes)
        excluded = bool(set(codes) - {"hwe_not_evaluable"}) and not forced
        variant_rows.append(dict(VARIANT_ID=variant, CHROMOSOME=chromosome, POSITION_BP=position_bp,
                                 N_MISS=int(missing[index].sum()), N_GENO=len(retained_samples.samples),
                                 F_MISS=fraction, MAF=float(maf[index]), HWE_P=float(hwe[index]) if control else None,
                                 HWE_STATUS="EVALUATED" if control and np.isfinite(hwe[index]) else "NOT_EVALUATED_DATA" if control else "NOT_EVALUATED_POLICY",
                                 QC_STATUS="RETAINED_EXCEPTION" if forced else "EXCLUDED" if excluded else "ALERT" if codes else "PASS",
                                 FILTER_CODES=",".join(codes), IS_TARGET_VARIANT=str(variant == target_id).lower(),
                                 FILTER_EXCEPTION=str(forced).lower(), PRESENT_AFTER_QC=str(not excluded).lower()))
    variant_mask = np.array([row["PRESENT_AFTER_QC"] == "true" for row in variant_rows])
    return QcExpectation(sample_rows, variant_rows, batch_rows,
                         subset(retained_samples, np.ones(len(retained_samples.samples), dtype=bool), variant_mask))


def confusion(expected: set[str], actual: set[str], universe: set[str]) -> dict[str, int | float | None]:
    """Matrice par unité homogène ; aucun dénominateur nul transformé en succès."""
    if not expected <= universe or not actual <= universe:
        raise ValueError("d20_confusion_unknown_unit")
    true_positive = len(expected & actual)
    false_positive = len(actual - expected)
    false_negative = len(expected - actual)
    true_negative = len(universe - expected - actual)
    return dict(tp=true_positive, fp=false_positive, fn=false_negative, tn=true_negative,
                sensitivity=true_positive / len(expected) if expected else None,
                specificity=true_negative / (true_negative + false_positive) if true_negative + false_positive else None)
