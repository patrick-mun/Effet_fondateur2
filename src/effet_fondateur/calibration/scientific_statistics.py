"""Statistiques préspecifiées communes aux producteurs scientifiques R13.

Les fonctions de ce module travaillent exclusivement sur les unités
indépendantes ou grappes déclarées par R07. Elles ne choisissent ni seuil ni
stratification après observation et retournent des valeurs non arrondies.
"""

from __future__ import annotations

from collections.abc import Callable, Mapping, Sequence
from dataclasses import asdict, dataclass
from hashlib import sha256
import hmac
import math
from typing import Any

import numpy as np
from numpy.typing import NDArray
from scipy.stats import beta, norm


DEFAULT_CONFIDENCE = 0.95
DEFAULT_BOOTSTRAP_REPLICATES = 10_000


class ScientificStatisticsError(ValueError):
    """Entrées statistiques invalides ou incompatibles avec R08."""

    def __init__(self, code: str, detail: str | None = None) -> None:
        self.code = code
        self.detail = detail
        super().__init__(code if detail is None else f"{code}: {detail}")


@dataclass(frozen=True)
class ProportionInterval:
    """Estimation et bornes unilatérales, sans arrondi d'affichage."""

    estimate: float
    lcb95: float
    ucb95: float
    successes: int
    trials: int
    method: str
    confidence: float


@dataclass(frozen=True)
class BootstrapInterval:
    """Bornes bootstrap par grappes et intervalle bilatéral à 95 %."""

    estimate: float
    lcb95: float
    ucb95: float
    ci95: tuple[float, float]
    independent_cluster_count: int
    observation_count: int
    bootstrap_replicates: int
    seed_fingerprint: str
    method: str
    confidence: float


def _validate_confidence(confidence: float) -> float:
    if (
        isinstance(confidence, bool)
        or not isinstance(confidence, (int, float))
        or not math.isfinite(float(confidence))
        or not 0.5 < float(confidence) < 1.0
    ):
        raise ScientificStatisticsError("confidence_invalid")
    return float(confidence)


def _validate_binomial(successes: int, trials: int) -> None:
    if (
        isinstance(successes, bool)
        or isinstance(trials, bool)
        or not isinstance(successes, int)
        or not isinstance(trials, int)
        or trials < 1
        or successes < 0
        or successes > trials
    ):
        raise ScientificStatisticsError(
            "binomial_counts_invalid", f"successes={successes}:trials={trials}"
        )


def wilson_one_sided(
    successes: int,
    trials: int,
    *,
    confidence: float = DEFAULT_CONFIDENCE,
) -> ProportionInterval:
    """Calcule les bornes de Wilson unilatérales sur des unités indépendantes."""
    _validate_binomial(successes, trials)
    confidence = _validate_confidence(confidence)
    estimate = successes / trials
    z = float(norm.ppf(confidence))
    z_squared = z * z
    denominator = 1.0 + z_squared / trials
    center = (estimate + z_squared / (2.0 * trials)) / denominator
    half_width = (
        z
        * math.sqrt(
            estimate * (1.0 - estimate) / trials
            + z_squared / (4.0 * trials * trials)
        )
        / denominator
    )
    return ProportionInterval(
        estimate=estimate,
        lcb95=max(0.0, center - half_width),
        ucb95=min(1.0, center + half_width),
        successes=successes,
        trials=trials,
        method="WILSON_ONE_SIDED",
        confidence=confidence,
    )


def exact_binomial_one_sided(
    successes: int,
    trials: int,
    *,
    confidence: float = DEFAULT_CONFIDENCE,
) -> ProportionInterval:
    """Calcule les bornes exactes unilatérales de Clopper–Pearson."""
    _validate_binomial(successes, trials)
    confidence = _validate_confidence(confidence)
    alpha = 1.0 - confidence
    lower = (
        0.0
        if successes == 0
        else float(beta.ppf(alpha, successes, trials - successes + 1))
    )
    upper = (
        1.0
        if successes == trials
        else float(beta.ppf(confidence, successes + 1, trials - successes))
    )
    return ProportionInterval(
        estimate=successes / trials,
        lcb95=lower,
        ucb95=upper,
        successes=successes,
        trials=trials,
        method="CLOPPER_PEARSON_ONE_SIDED",
        confidence=confidence,
    )


def derive_evaluation_seed(
    evaluation_secret: bytes,
    *,
    release_id: str,
    stage_id: str,
    metric_id: str,
    context_id: str = "GLOBAL",
) -> int:
    """Dérive une graine PCG64DXSM privée et stable par métrique."""
    if len(evaluation_secret) < 32:
        raise ScientificStatisticsError("evaluation_secret_too_short")
    fields = (release_id, stage_id, metric_id, context_id)
    if any(not isinstance(field, str) or not field for field in fields):
        raise ScientificStatisticsError("evaluation_seed_identity_invalid")
    payload = "\x1f".join(fields).encode("utf-8")
    digest = hmac.new(evaluation_secret, payload, sha256).digest()
    return int.from_bytes(digest[:16], "big")


def cluster_bootstrap_interval(
    values: Sequence[float] | NDArray[np.floating[Any]],
    cluster_ids: Sequence[str],
    *,
    statistic: Callable[[NDArray[np.float64]], float],
    seed: int,
    replicates: int = DEFAULT_BOOTSTRAP_REPLICATES,
    confidence: float = DEFAULT_CONFIDENCE,
) -> BootstrapInterval:
    """Rééchantillonne les grappes avec remise en conservant leurs observations.

    ``values`` doit contenir une valeur finie par ligne et ``cluster_ids``
    l'identifiant ``isolation_unit_id`` ou ``scenario_world_id`` correspondant.
    Les grappes, et non les lignes, constituent l'unité de rééchantillonnage.
    """
    confidence = _validate_confidence(confidence)
    numeric_values = np.asarray(values, dtype=np.float64)
    if numeric_values.ndim != 1 or numeric_values.size == 0:
        raise ScientificStatisticsError("bootstrap_values_invalid")
    if not np.isfinite(numeric_values).all():
        raise ScientificStatisticsError("bootstrap_values_non_finite")
    if len(cluster_ids) != numeric_values.size:
        raise ScientificStatisticsError("bootstrap_cluster_length_mismatch")
    if (
        isinstance(replicates, bool)
        or not isinstance(replicates, int)
        or replicates < 1
    ):
        raise ScientificStatisticsError("bootstrap_replicates_invalid")
    if (
        isinstance(seed, bool)
        or not isinstance(seed, int)
        or seed < 0
        or seed >= 2**128
    ):
        raise ScientificStatisticsError("bootstrap_seed_invalid")

    cluster_positions: dict[str, list[int]] = {}
    for index, cluster_id in enumerate(cluster_ids):
        if not isinstance(cluster_id, str) or not cluster_id:
            raise ScientificStatisticsError("bootstrap_cluster_id_invalid")
        cluster_positions.setdefault(cluster_id, []).append(index)
    clusters = tuple(cluster_positions)
    if not clusters:
        raise ScientificStatisticsError("bootstrap_clusters_missing")
    positions = tuple(
        np.asarray(cluster_positions[cluster_id], dtype=np.int64)
        for cluster_id in clusters
    )
    estimate = float(statistic(numeric_values))
    if not math.isfinite(estimate):
        raise ScientificStatisticsError("bootstrap_statistic_non_finite")

    generator = np.random.Generator(np.random.PCG64DXSM(seed))
    bootstrap_values = np.empty(replicates, dtype=np.float64)
    cluster_count = len(clusters)
    for replicate in range(replicates):
        selected = generator.integers(0, cluster_count, size=cluster_count)
        sampled_positions = np.concatenate([positions[index] for index in selected])
        bootstrap_values[replicate] = float(statistic(numeric_values[sampled_positions]))
    if not np.isfinite(bootstrap_values).all():
        raise ScientificStatisticsError("bootstrap_replicate_non_finite")

    alpha = 1.0 - confidence
    lcb, ucb = np.quantile(bootstrap_values, [alpha, confidence])
    ci_low, ci_high = np.quantile(
        bootstrap_values,
        [alpha / 2.0, 1.0 - alpha / 2.0],
    )
    seed_digest = sha256(seed.to_bytes(16, "big", signed=False)).hexdigest()
    return BootstrapInterval(
        estimate=estimate,
        lcb95=float(lcb),
        ucb95=float(ucb),
        ci95=(float(ci_low), float(ci_high)),
        independent_cluster_count=cluster_count,
        observation_count=int(numeric_values.size),
        bootstrap_replicates=replicates,
        seed_fingerprint=seed_digest,
        method="DETERMINISTIC_CLUSTER_PERCENTILE_PCG64DXSM",
        confidence=confidence,
    )


def holm_adjust(
    p_values: Mapping[str, float],
    *,
    alpha: float = 0.05,
) -> dict[str, dict[str, float | bool]]:
    """Applique Holm pas à pas et retourne les p-values ajustées monotones."""
    if not p_values:
        raise ScientificStatisticsError("holm_family_empty")
    if (
        isinstance(alpha, bool)
        or not isinstance(alpha, (int, float))
        or not math.isfinite(float(alpha))
        or not 0.0 < float(alpha) < 1.0
    ):
        raise ScientificStatisticsError("holm_alpha_invalid")
    checked = []
    for identifier, value in p_values.items():
        if not isinstance(identifier, str) or not identifier:
            raise ScientificStatisticsError("holm_identifier_invalid")
        if (
            isinstance(value, bool)
            or not isinstance(value, (int, float))
            or not math.isfinite(float(value))
            or not 0.0 <= float(value) <= 1.0
        ):
            raise ScientificStatisticsError("holm_p_value_invalid", identifier)
        checked.append((identifier, float(value)))
    ordered = sorted(checked, key=lambda item: (item[1], item[0]))
    family_size = len(ordered)
    previous_adjusted = 0.0
    rejection_open = True
    results: dict[str, dict[str, float | bool]] = {}
    for rank, (identifier, value) in enumerate(ordered, start=1):
        multiplier = family_size - rank + 1
        adjusted = min(1.0, max(previous_adjusted, multiplier * value))
        previous_adjusted = adjusted
        critical = float(alpha) / multiplier
        rejected = rejection_open and value <= critical
        if not rejected:
            rejection_open = False
        results[identifier] = {
            "p_value": value,
            "adjusted_p_value": adjusted,
            "holm_critical_value": critical,
            "rejected": rejected,
        }
    return {identifier: results[identifier] for identifier in p_values}


def monte_carlo_plus_one(
    exceedances: int,
    draws: int,
    *,
    confidence: float = DEFAULT_CONFIDENCE,
) -> dict[str, Any]:
    """Calcule ``(b+1)/(N+1)`` et l'incertitude binomiale des tirages nuls."""
    _validate_binomial(exceedances, draws)
    uncertainty = exact_binomial_one_sided(
        exceedances,
        draws,
        confidence=confidence,
    )
    return {
        "exceedances": exceedances,
        "draws": draws,
        "p_mc": (exceedances + 1) / (draws + 1),
        "binomial_uncertainty": asdict(uncertainty),
        "method": "MONTE_CARLO_PLUS_ONE",
    }
