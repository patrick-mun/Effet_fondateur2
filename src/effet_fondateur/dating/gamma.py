"""Estimateur Gamma de Gandolfo, Bahlo et Speed (2014)."""

from __future__ import annotations

import math
from dataclasses import dataclass
from typing import Sequence

import numpy as np
from scipy.stats import gamma as gamma_distribution


class GammaAgeError(ValueError):
    """Signale des longueurs ou paramètres incompatibles avec Gamma."""


@dataclass(frozen=True)
class GammaEstimate:
    """Estimation en générations et intervalle associé."""

    model: str
    sample_size: int
    effective_sample_size: float
    correlation: float | None
    generations: float
    confidence_lower: float
    confidence_upper: float
    confidence_level: float


def chance_sharing_correction_cm(
    median_allele_frequency: float,
    markers_on_chromosome: int,
    chromosome_length_cm: float,
) -> float:
    """Calcule la correction de partage fortuit décrite par l'implémentation R."""
    if not 0 < median_allele_frequency < 1:
        raise GammaAgeError("gamma_median_allele_frequency_invalid")
    if markers_on_chromosome <= 0 or chromosome_length_cm <= 0:
        raise GammaAgeError("gamma_chance_sharing_parameters_invalid")
    homozygosity = (
        median_allele_frequency**2 + (1.0 - median_allele_frequency) ** 2
    )
    matching_loci = math.log(0.01) / math.log(homozygosity)
    return matching_loci * chromosome_length_cm / markers_on_chromosome


def _lengths_morgans(
    left_lengths_cm: Sequence[float], right_lengths_cm: Sequence[float]
) -> tuple[np.ndarray, np.ndarray]:
    if len(left_lengths_cm) != len(right_lengths_cm):
        raise GammaAgeError("gamma_arm_count_mismatch")
    if len(left_lengths_cm) < 2:
        raise GammaAgeError("gamma_requires_two_segments")
    left = np.asarray(left_lengths_cm, dtype=float)
    right = np.asarray(right_lengths_cm, dtype=float)
    if (
        not np.all(np.isfinite(left))
        or not np.all(np.isfinite(right))
        or np.any(left <= 0)
        or np.any(right <= 0)
    ):
        raise GammaAgeError("gamma_lengths_must_be_positive_finite_cm")
    return left / 100.0, right / 100.0


def _chance_correction_morgans(chance_sharing_correction_cm: float) -> float:
    if (
        not math.isfinite(chance_sharing_correction_cm)
        or chance_sharing_correction_cm < 0
    ):
        raise GammaAgeError("gamma_chance_sharing_correction_invalid")
    return chance_sharing_correction_cm / 100.0


def _confidence_interval(
    estimate: float, shape: float, bias_corrected_rate: float,
    confidence_level: float,
) -> tuple[float, float]:
    if not 0.5 < confidence_level < 1.0:
        raise GammaAgeError("gamma_confidence_level_invalid")
    lower_probability = (1.0 - confidence_level) / 2.0
    lower_multiplier = gamma_distribution.ppf(
        lower_probability, a=shape, scale=1.0 / bias_corrected_rate
    )
    upper_multiplier = gamma_distribution.ppf(
        1.0 - lower_probability, a=shape, scale=1.0 / bias_corrected_rate
    )
    return estimate * float(lower_multiplier), estimate * float(upper_multiplier)


def estimate_independent_gamma(
    left_lengths_cm: Sequence[float], right_lengths_cm: Sequence[float],
    *, confidence_level: float = 0.95,
    chance_sharing_correction_cm: float = 0.0,
) -> GammaEstimate:
    """Calcule l'estimateur sans biais et son intervalle exact sous indépendance."""
    left, right = _lengths_morgans(left_lengths_cm, right_lengths_cm)
    sample_size = len(left)
    chance_correction = _chance_correction_morgans(
        chance_sharing_correction_cm if sample_size >= 10 else 0.0
    )
    arm_extension = (
        float(np.sum(left) + np.sum(right))
        - 2.0 * (sample_size - 1) * chance_correction
    ) / (2.0 * sample_size)
    corrected_sum = (
        float(np.sum(left) + np.sum(right)) + 2.0 * arm_extension
        - 2.0 * (sample_size - 1) * chance_correction
    )
    if corrected_sum <= 0:
        raise GammaAgeError("gamma_corrected_length_nonpositive")
    generations = (2.0 * sample_size - 1.0) / corrected_sum
    confidence_lower, confidence_upper = _confidence_interval(
        generations,
        shape=2.0 * sample_size,
        bias_corrected_rate=2.0 * sample_size - 1.0,
        confidence_level=confidence_level,
    )
    return GammaEstimate(
        model="INDEPENDENT", sample_size=sample_size,
        effective_sample_size=float(sample_size), correlation=None,
        generations=generations, confidence_lower=confidence_lower,
        confidence_upper=confidence_upper, confidence_level=confidence_level,
    )


def _effective_sample_size(total_lengths: np.ndarray) -> tuple[float, float]:
    sample_size = len(total_lengths)
    mean_length = float(np.mean(total_lengths))
    variance = float(np.var(total_lengths, ddof=1))
    denominator = sample_size * mean_length**2 + variance * (sample_size - 1)
    if denominator <= 0:
        raise GammaAgeError("gamma_correlation_denominator_invalid")
    correlation = (
        sample_size * mean_length**2 - variance * (1 + 2 * sample_size)
    ) / denominator
    if correlation < -2.0 / (sample_size - 1):
        effective_sample_size = sample_size / (
            1.0 + (sample_size - 1) * abs(correlation)
        )
    elif correlation < -1.0 / (sample_size - 1):
        effective_sample_size = float(sample_size)
    else:
        effective_sample_size = sample_size / (
            1.0 + (sample_size - 1) * correlation
        )
        effective_sample_size = min(float(sample_size), effective_sample_size)
    if not math.isfinite(effective_sample_size) or effective_sample_size <= 0.5:
        raise GammaAgeError("gamma_effective_sample_size_invalid")
    return correlation, effective_sample_size


def estimate_correlated_gamma(
    left_lengths_cm: Sequence[float], right_lengths_cm: Sequence[float],
    *, confidence_level: float = 0.95,
    chance_sharing_correction_cm: float = 0.0,
) -> GammaEstimate:
    """Calcule Gamma en remplaçant n par l'effectif effectif estimé des longueurs."""
    left, right = _lengths_morgans(left_lengths_cm, right_lengths_cm)
    sample_size = len(left)
    if sample_size < 3:
        raise GammaAgeError("gamma_correlated_requires_three_segments")
    chance_correction = _chance_correction_morgans(chance_sharing_correction_cm)
    arm_extension = (
        float(np.sum(left) + np.sum(right))
        - 2.0 * (sample_size - 1) * chance_correction
    ) / (2.0 * sample_size)
    corrected_left = left.copy()
    corrected_right = right.copy()
    corrected_left[int(np.argmax(corrected_left))] += arm_extension + chance_correction
    corrected_right[int(np.argmax(corrected_right))] += arm_extension + chance_correction
    total_lengths = corrected_left + corrected_right - 2.0 * chance_correction
    if np.any(total_lengths <= 0):
        raise GammaAgeError("gamma_corrected_length_nonpositive")
    correlation, effective_sample_size = _effective_sample_size(total_lengths)
    bias_correction = (2.0 * effective_sample_size - 1.0) / (
        2.0 * effective_sample_size
    )
    generations = (
        bias_correction * 2.0 * sample_size / float(np.sum(total_lengths))
    )
    confidence_lower, confidence_upper = _confidence_interval(
        generations,
        shape=2.0 * effective_sample_size,
        bias_corrected_rate=2.0 * effective_sample_size - 1.0,
        confidence_level=confidence_level,
    )
    return GammaEstimate(
        model="CORRELATED", sample_size=sample_size,
        effective_sample_size=effective_sample_size, correlation=correlation,
        generations=generations, confidence_lower=confidence_lower,
        confidence_upper=confidence_upper, confidence_level=confidence_level,
    )
