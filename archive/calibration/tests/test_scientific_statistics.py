from __future__ import annotations

import math

import numpy as np
import pytest

from effet_fondateur.calibration.scientific_statistics import (
    ScientificStatisticsError,
    cluster_bootstrap_interval,
    derive_evaluation_seed,
    exact_binomial_one_sided,
    holm_adjust,
    monte_carlo_plus_one,
    wilson_one_sided,
)


def test_wilson_uses_one_sided_95_percent_bounds() -> None:
    interval = wilson_one_sided(450, 500)

    assert interval.estimate == 0.9
    assert interval.method == "WILSON_ONE_SIDED"
    assert interval.lcb95 < interval.estimate < interval.ucb95
    assert interval.lcb95 == pytest.approx(0.8757336142)
    assert interval.ucb95 == pytest.approx(0.9199608141)


def test_exact_zero_event_bound_matches_closed_form() -> None:
    interval = exact_binomial_one_sided(0, 3000)

    expected_upper = 1.0 - math.pow(0.05, 1.0 / 3000)
    assert interval.estimate == 0.0
    assert interval.lcb95 == 0.0
    assert interval.ucb95 == pytest.approx(expected_upper)
    assert interval.ucb95 == pytest.approx(0.000998079012)


def test_cluster_bootstrap_is_deterministic_and_resamples_clusters() -> None:
    values = np.asarray([1.0, 2.0, 10.0, 11.0])
    clusters = ["FAMILY-A", "FAMILY-A", "FAMILY-B", "FAMILY-B"]

    first = cluster_bootstrap_interval(
        values,
        clusters,
        statistic=lambda sample: float(np.mean(sample)),
        seed=12345,
        replicates=1000,
    )
    second = cluster_bootstrap_interval(
        values,
        clusters,
        statistic=lambda sample: float(np.mean(sample)),
        seed=12345,
        replicates=1000,
    )

    assert first == second
    assert first.estimate == 6.0
    assert first.independent_cluster_count == 2
    assert first.observation_count == 4
    assert first.lcb95 <= first.estimate <= first.ucb95


def test_evaluation_seed_is_private_stable_and_context_specific() -> None:
    secret = b"scientific-evaluation-secret-32-bytes-minimum"
    first = derive_evaluation_seed(
        secret,
        release_id="release",
        stage_id="12",
        metric_id="common_ser",
        context_id="dense",
    )
    repeated = derive_evaluation_seed(
        secret,
        release_id="release",
        stage_id="12",
        metric_id="common_ser",
        context_id="dense",
    )
    other = derive_evaluation_seed(
        secret,
        release_id="release",
        stage_id="12",
        metric_id="common_ser",
        context_id="low",
    )

    assert first == repeated
    assert first != other
    assert 0 <= first < 2**128

    with pytest.raises(ScientificStatisticsError, match="too_short"):
        derive_evaluation_seed(
            b"short",
            release_id="release",
            stage_id="12",
            metric_id="common_ser",
        )


def test_holm_adjustment_is_monotone_and_step_down() -> None:
    adjusted = holm_adjust({"A": 0.01, "B": 0.03, "C": 0.04})

    assert adjusted["A"]["adjusted_p_value"] == pytest.approx(0.03)
    assert adjusted["B"]["adjusted_p_value"] == pytest.approx(0.06)
    assert adjusted["C"]["adjusted_p_value"] == pytest.approx(0.06)
    assert adjusted["A"]["rejected"] is True
    assert adjusted["B"]["rejected"] is False
    assert adjusted["C"]["rejected"] is False


def test_monte_carlo_plus_one_records_resolution_and_uncertainty() -> None:
    result = monte_carlo_plus_one(0, 9999)

    assert result["p_mc"] == pytest.approx(0.0001)
    assert result["draws"] == 9999
    assert result["binomial_uncertainty"]["ucb95"] > 0.0


@pytest.mark.parametrize(
    ("successes", "trials"),
    [(-1, 10), (11, 10), (0, 0)],
)
def test_binomial_intervals_reject_invalid_counts(successes: int, trials: int) -> None:
    with pytest.raises(ScientificStatisticsError, match="counts_invalid"):
        wilson_one_sided(successes, trials)
