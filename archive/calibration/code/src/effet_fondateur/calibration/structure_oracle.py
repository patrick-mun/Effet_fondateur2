"""D22 : PCA indépendante par Gram/eigh, diagnostics et incertitude descriptive.

Même modèle de standardisation que le contrat 08, algorithme différent de
sa SVD. Pour deux axes, la queue chi-deux est exactement exp(-distance²/2).
Les étiquettes privées n'interviennent qu'après reconstruction des scores.
"""

from __future__ import annotations

import math
from pathlib import Path
from typing import Any

import numpy as np

from effet_fondateur.calibration.kinship_oracle import wilson_interval
from effet_fondateur.calibration.technical_evaluation import BASES, GenotypeMatrix, TechnicalEvaluationError


def counted_dosages(matrix: GenotypeMatrix, bim_path: Path) -> np.ndarray:
    """Dosages A1 réels du BIM, sans assimiler A1 à REF ; axes sujets × loci."""
    bim = [line.split() for line in bim_path.read_text().splitlines() if line.strip()]
    if [row[1] for row in bim] != [variant[2] for variant in matrix.variants]:
        raise TechnicalEvaluationError("d22_bim_axis_mismatch")
    calls = matrix.alleles
    counted = np.array([BASES.index(row[4]) if row[4] in BASES else -2 for row in bim])
    dosage = np.sum(calls == counted[:, None, None], axis=2).astype(float)
    dosage[calls[:, :, 0] < 0] = np.nan
    return dosage.T


def fit_reference_pca(dosage: np.ndarray, reference: np.ndarray, parameters: dict[str, Any]) -> dict[str, Any]:
    """Reconstruit le modèle sur la seule référence ; aucune fonction de 08.

    Fréquences par chromosomes appelés ; manque imputé à 2p. Le Gram est
    sujet × sujet. Les axes nuls sont refusés plutôt que divisés par zéro.
    """
    if parameters["outlier_components"] != 2 or reference.sum() < 3:
        raise TechnicalEvaluationError("d22_oracle_domain_invalid")
    observed = dosage[reference]
    count = np.isfinite(observed).sum(axis=0)
    call_rate = count / len(observed)
    frequency = np.divide(np.nansum(observed, axis=0), 2 * count, out=np.full(dosage.shape[1], np.nan), where=count > 0)
    scale = np.sqrt(2 * frequency * (1 - frequency))
    eligible = call_rate >= parameters["min_reference_call_rate"]
    included = eligible & np.isfinite(scale) & (scale > 0)
    marker_count = int(included.sum())
    if marker_count < parameters["min_informative_variants"]:
        raise TechnicalEvaluationError("d22_insufficient_informative_variants")
    centered = dosage[:, included] - 2 * frequency[included]
    standardized = np.nan_to_num(centered, nan=0.) / scale[included] / math.sqrt(marker_count)
    reference_matrix = standardized[reference]
    values, vectors = np.linalg.eigh(reference_matrix @ reference_matrix.T)
    components = min(parameters["requested_components"], len(observed) - 1, marker_count)
    selected = np.argsort(values)[::-1][:components]
    squared_singular = values[selected]
    if np.any(squared_singular <= np.finfo(float).eps):
        raise TechnicalEvaluationError("d22_degenerate_reference_axes")
    singular = np.sqrt(squared_singular)
    weights = reference_matrix.T @ vectors[:, selected] / singular
    for axis in range(components):
        if weights[np.argmax(np.abs(weights[:, axis])), axis] < 0:
            weights[:, axis] *= -1
    scores = standardized @ weights
    eigenvalues = squared_singular / (len(observed) - 1)
    ratios = squared_singular / np.sum(reference_matrix ** 2)
    distance = np.sum(scores[:, :2] ** 2 / eigenvalues[:2], axis=1)
    return dict(frequencies=frequency, scales=scale, call_rates=call_rate, call_rate_eligible=eligible,
                included=included, weights=weights, scores=scores, eigenvalues=eigenvalues,
                singular_values=singular, explained_ratios=ratios, distances=distance,
                p_values=np.exp(-distance / 2), outliers=np.exp(-distance / 2) < parameters["outlier_alpha"],
                reference=reference, components=components, informative_count=marker_count)


def compare_tables(expected: list[dict[str, Any]], actual: list[dict[str, str]], key: str, config: dict[str, Any]) -> dict[str, int]:
    """Vérifie toutes les colonnes attendues, ordre et nombres sans NaN fabriqué."""
    wanted, observed = [row[key] for row in expected], [row[key] for row in actual]
    by_key = {row[key]: row for row in expected}
    result = dict(lost=len(set(wanted) - set(observed)), extra=len(set(observed) - set(wanted)),
                  duplicate=len(observed) - len(set(observed)), order=sum(left != right for left, right in zip(wanted, observed)) + abs(len(wanted) - len(observed)),
                  numeric=0, decision=0, compared=len(actual))
    for row in actual:
        if row[key] not in by_key:
            continue
        for field, value in by_key[row[key]].items():
            cell = row.get(field)
            if isinstance(value, (float, np.floating)):
                try:
                    matches = cell == "" if not np.isfinite(value) else math.isfinite(float(cell)) and math.isclose(float(cell), value, rel_tol=config["relative_tolerance"], abs_tol=0. if field == "P_VALUE" else config["absolute_tolerance"])
                except (ValueError, TypeError):
                    matches = False
                result["numeric"] += not matches
            else:
                result["decision"] += cell != str(value)
    return result


def detection_metrics(actual: np.ndarray, injected: np.ndarray, selection: np.ndarray, confidence: float) -> dict[str, Any]:
    """Sensibilité/FPR ; Wilson seulement si les unités sélectionnées sont libres."""
    positive, negative = selection & injected, selection & ~injected
    true_positive, false_positive = int((actual & positive).sum()), int((actual & negative).sum())
    return dict(tp=true_positive, fn=int(positive.sum()) - true_positive, fp=false_positive, tn=int(negative.sum()) - false_positive,
                sensitivity=true_positive / int(positive.sum()) if positive.any() else None,
                non_injected_outlier_rate=false_positive / int(negative.sum()) if negative.any() else None,
                sensitivity_wilson=wilson_interval(true_positive, int(positive.sum()), confidence),
                non_injected_rate_wilson=wilson_interval(false_positive, int(negative.sum()), confidence))


def structure_metrics(scores: np.ndarray, truth: dict[str, Any], outliers: np.ndarray, confidence: float) -> dict[str, Any]:
    """Agrégats privés après clôture ; aucun score ni identifiant individuel publié."""
    reference_count = truth["reference_count"]
    groups = np.array(truth["groups"])
    reference = np.arange(len(groups)) < reference_count
    first, second = scores[reference & (groups == "A"), :2], scores[reference & (groups == "B"), :2]
    center_a, center_b = first.mean(axis=0), second.mean(axis=0)
    direction = center_b - center_a
    squared_separation = float(direction @ direction)
    within = float((np.sum((first - center_a) ** 2) + np.sum((second - center_b) ** 2)) / (len(first) + len(second) - 2))
    mixed = scores[reference & (groups == "MIXED"), :2]
    coordinates = (mixed - center_a) @ direction / squared_separation
    injected = groups == "RARE"
    independent = detection_metrics(outliers, injected, reference, confidence)
    dependent = detection_metrics(outliers, injected, np.ones(len(groups), bool), confidence)
    dependent.pop("sensitivity_wilson")
    dependent.pop("non_injected_rate_wilson")
    return dict(centroid_separation_over_within_sd=math.sqrt(squared_separation / within),
                mixed_coordinate_mean=float(coordinates.mean()), mixed_coordinate_min=float(coordinates.min()), mixed_coordinate_max=float(coordinates.max()),
                mixed_absolute_error_from_half_mean=float(np.mean(np.abs(coordinates - .5))),
                independent_reference_detection=independent, all_samples_detection_without_independent_interval=dependent,
                reference_sample_count=reference_count, projected_sample_count=len(groups) - reference_count,
                unsupported_geographic_origin_claim_count=0, geographic_claim_audit_scope="AGGREGATE_REPORT_TEXT_ONLY_NO_POPULATION_NAMES_ASSIGNED")
