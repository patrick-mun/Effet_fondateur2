from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest

from effet_fondateur.calibration.r13_scientific_evaluator import (
    evaluate_metric_observation,
    load_metric_rule_bundle,
)
from effet_fondateur.calibration.scientific_metric_producers import (
    ScientificMetricProducerError,
    audit_metric_producer_methods,
    cluster_statistic_observation,
    exact_observation,
    proportion_observation,
)


REPOSITORY_ROOT = Path(__file__).resolve().parents[2]


def test_exact_producer_does_not_invent_an_interval() -> None:
    observation = exact_observation(
        metric_id="chromosome_completeness",
        estimate=22,
        independent_n=1,
        independent_unit="release",
    )

    assert observation["estimate"] == 22.0
    assert "lcb95" not in observation
    assert "ucb95" not in observation
    assert observation["producer"]["interval_generated"] is False


def test_wilson_producer_is_directly_consumable_by_r08_engine() -> None:
    observation = proportion_observation(
        metric_id="sensitivity",
        successes=480,
        trials=500,
        independent_unit="founder_event",
        method="WILSON_ONE_SIDED",
    )
    decision = evaluate_metric_observation(
        {
            "metric_id": "sensitivity",
            "operator": "min",
            "value": 0.90,
            "decision_value": "lcb95",
            "minimum_independent_n": 500,
        },
        observation,
    )

    assert observation["producer"]["successes"] == 480
    assert decision["status"] == "PASS"


def test_exact_binomial_producer_preserves_rare_zero_event_bound() -> None:
    observation = proportion_observation(
        metric_id="negative_pair_fpr",
        successes=0,
        trials=3000,
        independent_unit="negative_pair",
        method="CLOPPER_PEARSON_ONE_SIDED",
    )
    decision = evaluate_metric_observation(
        {
            "metric_id": "negative_pair_fpr",
            "operator": "max",
            "value": 0.001,
            "decision_value": "ucb95",
            "minimum_independent_n": 3000,
        },
        observation,
    )

    assert observation["ucb95"] == pytest.approx(0.000998079012)
    assert decision["status"] == "PASS"


def test_cluster_producer_bootstraps_independent_clusters_and_is_reproducible() -> None:
    arguments = {
        "metric_id": "boundary_error",
        "values": np.asarray([0.1, 0.2, 0.3, 0.4]),
        "cluster_ids": ["WORLD-1", "WORLD-1", "WORLD-2", "WORLD-2"],
        "independent_unit": "scenario_world",
        "statistic_id": "MEAN_ABSOLUTE",
        "evaluation_secret": b"private-evaluation-secret-at-least-32-bytes",
        "release_id": "release",
        "stage_id": "16C",
        "replicates": 1000,
    }

    first = cluster_statistic_observation(**arguments)
    second = cluster_statistic_observation(**arguments)

    assert first == second
    assert first["estimate"] == pytest.approx(0.25)
    assert first["independent_n"] == 2
    assert first["producer"]["bootstrap_replicates"] == 1000
    assert len(first["producer"]["seed_fingerprint"]) == 64


def test_context_requires_both_identifier_and_public_condition() -> None:
    with pytest.raises(ScientificMetricProducerError, match="context_invalid"):
        proportion_observation(
            metric_id="common_ser",
            successes=1,
            trials=100,
            independent_unit="informative_transition",
            method="WILSON_ONE_SIDED",
            context_id="dense",
            context=None,
        )


def test_all_frozen_rules_have_a_prespecified_inference_method() -> None:
    audit = audit_metric_producer_methods(load_metric_rule_bundle(REPOSITORY_ROOT))

    assert audit["state"] == "METHODS_AND_STAGE_BINDINGS_FROZEN"
    assert audit["stage_count"] == 25
    assert audit["rule_binding_count"] == 135
    assert audit["method_counts"]["CLOPPER_PEARSON_ONE_SIDED"] == 4
    assert audit["method_counts"]["WILSON_ONE_SIDED"] > 20
    assert audit["method_counts"]["DETERMINISTIC_CLUSTER_BOOTSTRAP_10000"] > 20
    assert audit["stage_producers_complete"] is True
    assert audit["scientific_evaluation_authorized"] is False
