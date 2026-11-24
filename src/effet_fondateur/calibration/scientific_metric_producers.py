"""Producteurs génériques d'observations métriques scientifiques R13.

Ce module transforme des unités indépendantes déjà construites par un
évaluateur d'étape en observations non arrondies. Il ne choisit jamais la
population, la strate, le seuil ou la méthode après lecture du résultat.
"""

from __future__ import annotations

from collections.abc import Mapping, Sequence
from dataclasses import asdict
import math
from typing import Any

import numpy as np
from numpy.typing import NDArray

from effet_fondateur.calibration.scientific_statistics import (
    DEFAULT_BOOTSTRAP_REPLICATES,
    ScientificStatisticsError,
    cluster_bootstrap_interval,
    derive_evaluation_seed,
    exact_binomial_one_sided,
    wilson_one_sided,
)


PROPORTION_METHODS = frozenset({"WILSON_ONE_SIDED", "CLOPPER_PEARSON_ONE_SIDED"})
CLUSTER_STATISTICS = frozenset(
    {
        "MEAN",
        "MEDIAN",
        "MEAN_ABSOLUTE",
        "ABSOLUTE_MEAN",
        "MEDIAN_ABSOLUTE",
        "ROOT_MEAN_SQUARE",
        "QUANTILE_95",
        "STANDARD_DEVIATION",
    }
)
EXACT_BINOMIAL_METRICS = frozenset(
    {
        ("13", "null_fpr"),
        ("16", "roh_fpr"),
        ("16B", "null_fpr"),
        ("16C", "negative_pair_fpr"),
    }
)
PRODUCER_METHOD_STATE = "METHODS_AND_STAGE_BINDINGS_FROZEN"


class ScientificMetricProducerError(ValueError):
    """Observation impossible à produire sans violer le contrat scientifique."""

    def __init__(self, code: str, detail: str | None = None) -> None:
        self.code = code
        self.detail = detail
        super().__init__(code if detail is None else f"{code}: {detail}")


def _identity(metric_id: str, independent_unit: str) -> None:
    if not isinstance(metric_id, str) or not metric_id:
        raise ScientificMetricProducerError("metric_id_invalid")
    if not isinstance(independent_unit, str) or not independent_unit:
        raise ScientificMetricProducerError("independent_unit_invalid")


def _context_fields(
    context_id: str | None,
    context: Mapping[str, Any] | None,
) -> dict[str, Any]:
    if context_id is None and context is None:
        return {}
    if (
        not isinstance(context_id, str)
        or not context_id
        or not isinstance(context, Mapping)
        or not context
    ):
        raise ScientificMetricProducerError("metric_context_invalid")
    return {"context_id": context_id, "context": dict(context)}


def inference_method_for_rule(
    stage_id: str,
    rule: Mapping[str, Any],
) -> str:
    """Résout la méthode R08 sans dépendre d'une valeur observée."""
    metric_id = rule.get("metric_id")
    operator = rule.get("operator")
    unit = rule.get("unit")
    if not isinstance(stage_id, str) or not isinstance(metric_id, str):
        raise ScientificMetricProducerError("producer_rule_identity_invalid")
    deterministic = operator == "exact" or rule.get("deterministic") is True
    if operator == "contextual":
        contexts = rule.get("contexts")
        if not isinstance(contexts, list) or not contexts:
            raise ScientificMetricProducerError(
                "producer_contextual_rule_invalid", f"{stage_id}:{metric_id}"
            )
        deterministic = all(
            isinstance(context, Mapping) and context.get("deterministic") is True
            for context in contexts
        )
    if deterministic:
        return "DETERMINISTIC_EXACT"
    if (stage_id, metric_id) in EXACT_BINOMIAL_METRICS:
        if unit != "proportion":
            raise ScientificMetricProducerError(
                "exact_binomial_unit_mismatch", f"{stage_id}:{metric_id}"
            )
        return "CLOPPER_PEARSON_ONE_SIDED"
    if unit == "proportion":
        return "WILSON_ONE_SIDED"
    if operator in {"min", "max", "equivalence", "contextual"}:
        return "DETERMINISTIC_CLUSTER_BOOTSTRAP_10000"
    raise ScientificMetricProducerError(
        "producer_method_unresolved", f"{stage_id}:{metric_id}"
    )


def audit_metric_producer_methods(bundle: Mapping[str, Any]) -> dict[str, Any]:
    """Prouve qu'une méthode préspecifiée couvre chaque règle primaire."""
    stages = bundle.get("stages")
    if not isinstance(stages, Mapping):
        raise ScientificMetricProducerError("producer_bundle_stages_invalid")
    method_counts: dict[str, int] = {}
    binding_count = 0
    for stage_id, stage in stages.items():
        if not isinstance(stage_id, str) or not isinstance(stage, Mapping):
            raise ScientificMetricProducerError("producer_stage_invalid")
        rules = stage.get("rules")
        if not isinstance(rules, list) or not rules:
            raise ScientificMetricProducerError("producer_rules_missing", stage_id)
        for rule in rules:
            if not isinstance(rule, Mapping):
                raise ScientificMetricProducerError("producer_rule_invalid", stage_id)
            method = inference_method_for_rule(stage_id, rule)
            method_counts[method] = method_counts.get(method, 0) + 1
            binding_count += 1
    return {
        "state": PRODUCER_METHOD_STATE,
        "stage_count": len(stages),
        "rule_binding_count": binding_count,
        "method_counts": dict(sorted(method_counts.items())),
        "stage_producers_complete": True,
        "scientific_evaluation_authorized": False,
    }


def exact_observation(
    *,
    metric_id: str,
    estimate: float,
    independent_n: int,
    independent_unit: str,
    context_id: str | None = None,
    context: Mapping[str, Any] | None = None,
) -> dict[str, Any]:
    """Produit un contrôle déterministe sans inventer d'intervalle."""
    _identity(metric_id, independent_unit)
    if (
        isinstance(estimate, bool)
        or not isinstance(estimate, (int, float))
        or not math.isfinite(float(estimate))
    ):
        raise ScientificMetricProducerError("exact_estimate_invalid", metric_id)
    if (
        isinstance(independent_n, bool)
        or not isinstance(independent_n, int)
        or independent_n < 0
    ):
        raise ScientificMetricProducerError("independent_n_invalid", metric_id)
    return {
        "metric_id": metric_id,
        "independent_n": independent_n,
        "independent_unit": independent_unit,
        "integrity_passed": True,
        "estimate": float(estimate),
        "producer": {
            "method": "DETERMINISTIC_EXACT",
            "interval_generated": False,
            "decision_precision": "FULL_BEFORE_DISPLAY_ROUNDING",
        },
        **_context_fields(context_id, context),
    }


def proportion_observation(
    *,
    metric_id: str,
    successes: int,
    trials: int,
    independent_unit: str,
    method: str,
    context_id: str | None = None,
    context: Mapping[str, Any] | None = None,
) -> dict[str, Any]:
    """Produit une proportion Wilson ou exacte selon la méthode préspecifiée."""
    _identity(metric_id, independent_unit)
    if method not in PROPORTION_METHODS:
        raise ScientificMetricProducerError("proportion_method_invalid", method)
    try:
        interval = (
            wilson_one_sided(successes, trials)
            if method == "WILSON_ONE_SIDED"
            else exact_binomial_one_sided(successes, trials)
        )
    except ScientificStatisticsError as error:
        raise ScientificMetricProducerError(
            "proportion_inputs_invalid", metric_id
        ) from error
    return {
        "metric_id": metric_id,
        "independent_n": trials,
        "independent_unit": independent_unit,
        "integrity_passed": True,
        "estimate": interval.estimate,
        "lcb95": interval.lcb95,
        "ucb95": interval.ucb95,
        "producer": {
            "method": interval.method,
            "confidence": interval.confidence,
            "successes": successes,
            "trials": trials,
            "decision_precision": "FULL_BEFORE_DISPLAY_ROUNDING",
        },
        **_context_fields(context_id, context),
    }


def _statistic(statistic_id: str):
    if statistic_id not in CLUSTER_STATISTICS:
        raise ScientificMetricProducerError(
            "cluster_statistic_invalid", statistic_id
        )
    return {
        "MEAN": lambda values: float(np.mean(values)),
        "MEDIAN": lambda values: float(np.median(values)),
        "MEAN_ABSOLUTE": lambda values: float(np.mean(np.abs(values))),
        "ABSOLUTE_MEAN": lambda values: float(abs(np.mean(values))),
        "MEDIAN_ABSOLUTE": lambda values: float(np.median(np.abs(values))),
        "ROOT_MEAN_SQUARE": lambda values: float(
            np.sqrt(np.mean(np.square(values)))
        ),
        "QUANTILE_95": lambda values: float(np.quantile(values, 0.95)),
        "STANDARD_DEVIATION": lambda values: float(np.std(values, ddof=0)),
    }[statistic_id]


def cluster_statistic_observation(
    *,
    metric_id: str,
    values: Sequence[float] | NDArray[np.floating[Any]],
    cluster_ids: Sequence[str],
    independent_unit: str,
    statistic_id: str,
    evaluation_secret: bytes,
    release_id: str,
    stage_id: str,
    replicates: int = DEFAULT_BOOTSTRAP_REPLICATES,
    context_id: str | None = None,
    context: Mapping[str, Any] | None = None,
) -> dict[str, Any]:
    """Produit une statistique continue avec bootstrap des grappes R07."""
    _identity(metric_id, independent_unit)
    statistic = _statistic(statistic_id)
    try:
        seed = derive_evaluation_seed(
            evaluation_secret,
            release_id=release_id,
            stage_id=stage_id,
            metric_id=metric_id,
            context_id=context_id or "GLOBAL",
        )
        interval = cluster_bootstrap_interval(
            values,
            cluster_ids,
            statistic=statistic,
            seed=seed,
            replicates=replicates,
        )
    except ScientificStatisticsError as error:
        raise ScientificMetricProducerError(
            "cluster_statistic_inputs_invalid", metric_id
        ) from error
    interval_document = asdict(interval)
    return {
        "metric_id": metric_id,
        "independent_n": interval.independent_cluster_count,
        "independent_unit": independent_unit,
        "integrity_passed": True,
        "estimate": interval.estimate,
        "lcb95": interval.lcb95,
        "ucb95": interval.ucb95,
        "ci95": list(interval.ci95),
        "producer": {
            "method": interval.method,
            "statistic": statistic_id,
            "confidence": interval.confidence,
            "bootstrap_replicates": interval.bootstrap_replicates,
            "observation_count": interval.observation_count,
            "independent_cluster_count": interval.independent_cluster_count,
            "seed_fingerprint": interval.seed_fingerprint,
            "decision_precision": "FULL_BEFORE_DISPLAY_ROUNDING",
            "interval": interval_document,
        },
        **_context_fields(context_id, context),
    }
