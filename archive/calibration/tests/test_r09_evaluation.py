"""Ordre de décision et intersection–union R08."""

import pytest

from effet_fondateur.calibration.evaluation import (
    EvaluationStatus,
    decide_metric,
    decide_stage,
)


def test_minimum_uses_point_estimate_then_lower_bound() -> None:
    rule = {"metric_id": "sensitivity", "operator": "min", "value": 0.9}
    assert decide_metric(rule, {"estimate": 0.89, "lcb95": 0.85}).status is EvaluationStatus.FAIL
    assert decide_metric(rule, {"estimate": 0.93, "lcb95": 0.89}).status is EvaluationStatus.NOT_EVALUABLE
    assert decide_metric(rule, {"estimate": 0.94, "lcb95": 0.90}).status is EvaluationStatus.PASS


def test_maximum_uses_point_estimate_then_upper_bound() -> None:
    rule = {"metric_id": "fpr", "operator": "max", "value": 0.01}
    assert decide_metric(rule, {"estimate": 0.02, "ucb95": 0.03}).status is EvaluationStatus.FAIL
    assert decide_metric(rule, {"estimate": 0.005, "ucb95": 0.02}).status is EvaluationStatus.NOT_EVALUABLE
    assert decide_metric(rule, {"estimate": 0.004, "ucb95": 0.01}).status is EvaluationStatus.PASS


def test_equivalence_requires_the_complete_interval_inside_margin() -> None:
    rule = {"metric_id": "delta", "operator": "equivalence", "value": 0.2}
    assert decide_metric(
        rule, {"estimate": 0.05, "ci95_low": -0.2, "ci95_high": 0.2}
    ).status is EvaluationStatus.PASS
    assert decide_metric(
        rule, {"estimate": 0.05, "ci95_low": -0.25, "ci95_high": 0.2}
    ).status is EvaluationStatus.NOT_EVALUABLE


def test_structural_reason_is_not_evaluable_but_unknown_reason_needs_review() -> None:
    rule = {"metric_id": "phase", "operator": "min", "value": 0.9}
    assert decide_metric(
        rule, {"structural_reason": "INSUFFICIENT_INFORMATIVE_TRANSITIONS"}
    ).status is EvaluationStatus.NOT_EVALUABLE
    assert decide_metric(
        rule, {"structural_reason": "UNDECLARED_REASON"}
    ).status is EvaluationStatus.REVIEW_REQUIRED


def test_hard_integrity_failure_precedes_all_metrics() -> None:
    decision = decide_stage(
        stage_id="05",
        rules=[{"metric_id": "sensitivity", "operator": "min", "value": 0.9}],
        observations={"sensitivity": {"estimate": 1.0, "lcb95": 1.0}},
        hard_fail_reasons=["truth_leak"],
    )
    assert decision.status is EvaluationStatus.FAIL
    assert decision.metric_decisions == ()


def test_all_primary_metrics_must_pass() -> None:
    decision = decide_stage(
        stage_id="07",
        rules=[
            {"metric_id": "sensitivity", "operator": "min", "value": 0.9},
            {"metric_id": "mae", "operator": "max", "value": 0.03},
        ],
        observations={
            "sensitivity": {"estimate": 0.95, "lcb95": 0.90},
            "mae": {"estimate": 0.04, "ucb95": 0.05},
        },
    )
    assert decision.status is EvaluationStatus.FAIL


def test_missing_required_metric_is_fail_not_not_evaluable() -> None:
    decision = decide_stage(
        stage_id="10",
        rules=[{"metric_id": "specificity", "operator": "min", "value": 0.99}],
        observations={},
    )
    assert decision.status is EvaluationStatus.FAIL
    assert decision.reasons == ("required_metric_missing:specificity",)


def test_manual_review_is_requested_only_after_metrics_pass() -> None:
    decision = decide_stage(
        stage_id="09",
        rules=[{"metric_id": "cohort_exact", "operator": "exact", "value": 1}],
        observations={"cohort_exact": {"estimate": 1}},
        manual_review_pending=True,
    )
    assert decision.status is EvaluationStatus.REVIEW_REQUIRED


@pytest.mark.parametrize("value,status", [(1, EvaluationStatus.PASS), (0, EvaluationStatus.FAIL)])
def test_exact_contracts_are_not_diluted_by_intervals(value: int, status: EvaluationStatus) -> None:
    rule = {"metric_id": "integrity", "operator": "exact", "value": 1}
    assert decide_metric(rule, {"estimate": value}).status is status
