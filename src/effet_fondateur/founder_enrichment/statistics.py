"""Probabilités empiriques, quantiles et intervalle Monte-Carlo."""

from __future__ import annotations

import math
from collections.abc import Sequence

from .model import NullDraw, NullSummary


def empirical_probability(exceedance_count: int, evaluable_count: int) -> float:
    """Applique la correction conservatrice préspécifiée ``(1+k)/(1+N)``."""
    if evaluable_count < 1 or exceedance_count < 0 or exceedance_count > evaluable_count:
        raise ValueError("invalid_empirical_probability_counts")
    return (1.0 + exceedance_count) / (1.0 + evaluable_count)


def _wilson_interval(successes: int, trials: int, z: float = 1.959963984540054) -> tuple[float, float]:
    if trials < 1:
        raise ValueError("wilson_interval_requires_trials")
    proportion = successes / trials
    denominator = 1.0 + z * z / trials
    centre = (proportion + z * z / (2.0 * trials)) / denominator
    half_width = z * math.sqrt(proportion * (1.0 - proportion) / trials + z * z / (4.0 * trials * trials)) / denominator
    return max(0.0, centre - half_width), min(1.0, centre + half_width)


def summarize_null(
    draws: Sequence[NullDraw], observed_total_cm: float, *, exhaustive: bool
) -> NullSummary:
    """Résume les seules valeurs évaluables tout en comptant les exclusions."""
    evaluable = [draw for draw in draws if draw.statistic.evaluation_status == "EVALUATED"]
    if not evaluable:
        return NullSummary(len(draws), 0, len(draws), 0, None, None, None, None)
    exceedances = sum(
        float(draw.statistic.total_shared_cm) >= observed_total_cm for draw in evaluable
    )
    exact = exceedances / len(evaluable) if exhaustive else None
    interval_low, interval_high = _wilson_interval(exceedances, len(evaluable))
    return NullSummary(
        attempted_draws=len(draws), evaluable_draws=len(evaluable),
        non_evaluable_draws=len(draws) - len(evaluable),
        exceedance_count=exceedances,
        empirical_probability=empirical_probability(exceedances, len(evaluable)),
        exact_probability=exact, interval_low=interval_low, interval_high=interval_high,
    )
