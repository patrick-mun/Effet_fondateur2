"""Types immuables du test d'enrichissement haplotypique."""

from __future__ import annotations

from dataclasses import dataclass


@dataclass(frozen=True)
class Marker:
    """Marqueur ordonné, exprimé en bp et cM."""

    variant_id: str
    position_bp: int
    position_cm: float
    is_target: bool = False


@dataclass(frozen=True)
class HaplotypeProfile:
    """Allèles 0/1 d'une copie chromosomique ; ``None`` représente un manque."""

    individual_id: str
    haplotype_id: str
    alleles: tuple[str | None, ...]
    family_id: str | None = None


@dataclass(frozen=True)
class SharingStatistic:
    """Longueurs du partage exact, mesurées séparément autour de la cible."""

    evaluation_status: str
    left_shared_cm: float | None
    right_shared_cm: float | None
    total_shared_cm: float | None
    left_shared_bp: int | None
    right_shared_bp: int | None
    left_marker_count: int | None
    right_marker_count: int | None
    left_bound_index: int | None
    right_bound_index: int | None
    non_evaluable_reason: str | None = None


@dataclass(frozen=True)
class NullDraw:
    """Un tirage nul, y compris lorsqu'il n'est pas évaluable."""

    source: str
    stratum: str
    draw_index: int
    attempt_index: int
    individual_ids: tuple[str, ...]
    haplotype_ids: tuple[str, ...]
    statistic: SharingStatistic


@dataclass(frozen=True)
class NullSummary:
    """Résumé numérique d'une distribution nulle."""

    attempted_draws: int
    evaluable_draws: int
    non_evaluable_draws: int
    exceedance_count: int
    empirical_probability: float | None
    exact_probability: float | None
    interval_low: float | None
    interval_high: float | None
