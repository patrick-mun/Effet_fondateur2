"""Comparaison exploratoire du LD ciblé à un fond génomique apparié.

Les calculs portent sur des haplotypes phasés et utilisent un seul représentant
par unité familiale. Le module ne prend aucune décision d'effet fondateur : il
produit une direction de preuve secondaire, explicitement limitée par le petit
nombre de familles porteuses.
"""

from __future__ import annotations

import math
import statistics
from bisect import bisect_left, bisect_right
from dataclasses import dataclass
from typing import Iterable, Sequence


class MatchedLdError(ValueError):
    """Signale une entrée phasée ou un paramètre de comparaison incohérent."""


@dataclass(frozen=True)
class PhasedVariant:
    """Variant biallélique phasé ; la position génétique est exprimée en cM."""

    chromosome: str
    variant_id: str
    position_bp: int
    position_cm: float
    haplotypes: tuple[int | None, ...]

    @property
    def called_haplotype_count(self) -> int:
        return sum(value is not None for value in self.haplotypes)

    @property
    def alt_allele_frequency(self) -> float | None:
        called = [value for value in self.haplotypes if value is not None]
        if not called:
            return None
        return sum(called) / len(called)

    @property
    def minor_allele_frequency(self) -> float | None:
        frequency = self.alt_allele_frequency
        return None if frequency is None else min(frequency, 1.0 - frequency)


@dataclass(frozen=True)
class PairMetric:
    """Mesure r² entre un variant ancre et un voisin."""

    anchor_id: str
    neighbor_id: str
    chromosome: str
    anchor_position_bp: int
    neighbor_position_bp: int
    anchor_position_cm: float
    neighbor_position_cm: float
    signed_distance_cm: float
    distance_bin: str
    anchor_maf: float
    neighbor_maf: float
    called_haplotype_count: int
    r2: float


def parse_phased_genotype(value: str) -> tuple[int | None, int | None]:
    """Convertit un GT diploïde phasé en deux allèles 0/1."""

    genotype = value.split(":", 1)[0]
    if genotype in {".|.", "."}:
        return None, None
    fields = genotype.split("|")
    if len(fields) != 2:
        raise MatchedLdError(f"unphased_or_malformed_genotype:{genotype}")
    result: list[int | None] = []
    for field in fields:
        if field == ".":
            result.append(None)
        elif field in {"0", "1"}:
            result.append(int(field))
        else:
            raise MatchedLdError(f"non_biallelic_genotype:{genotype}")
    return result[0], result[1]


def r_squared(
    first: Sequence[int | None], second: Sequence[int | None], *, minimum_called: int = 20,
) -> tuple[float | None, int]:
    """Calcule le carré de la corrélation allélique sur les copies appelées."""

    if len(first) != len(second):
        raise MatchedLdError("haplotype_length_mismatch")
    pairs = [(left, right) for left, right in zip(first, second, strict=True)
             if left is not None and right is not None]
    if len(pairs) < minimum_called:
        return None, len(pairs)
    left_mean = sum(left for left, _ in pairs) / len(pairs)
    right_mean = sum(right for _, right in pairs) / len(pairs)
    left_variance = sum((left - left_mean) ** 2 for left, _ in pairs)
    right_variance = sum((right - right_mean) ** 2 for _, right in pairs)
    if left_variance <= 0 or right_variance <= 0:
        return None, len(pairs)
    covariance = sum(
        (left - left_mean) * (right - right_mean) for left, right in pairs
    )
    value = covariance * covariance / (left_variance * right_variance)
    return min(1.0, max(0.0, value)), len(pairs)


def bin_labels(distance_bins_cm: Sequence[float]) -> tuple[str, ...]:
    """Retourne des identifiants stables pour les classes de distance."""

    previous = 0.0
    labels: list[str] = []
    for upper in distance_bins_cm:
        if not math.isfinite(upper) or upper <= previous:
            raise MatchedLdError("invalid_distance_bins")
        labels.append(f"{previous:g}-{upper:g}_cM")
        previous = upper
    return tuple(labels)


def distance_bin(distance_cm: float, distance_bins_cm: Sequence[float]) -> str | None:
    """Classe une distance absolue, avec borne supérieure incluse."""

    labels = bin_labels(distance_bins_cm)
    lower = 0.0
    for label, upper in zip(labels, distance_bins_cm, strict=True):
        if lower < distance_cm <= upper or (lower == 0 and distance_cm == 0):
            return label
        lower = upper
    return None


def centered_pair_metrics(
    anchor: PhasedVariant,
    variants: Sequence[PhasedVariant],
    *,
    distance_bins_cm: Sequence[float],
    minimum_neighbor_maf: float,
    minimum_called_haplotypes: int,
) -> list[PairMetric]:
    """Calcule le LD entre une ancre rare et ses voisins communs."""

    anchor_maf = anchor.minor_allele_frequency
    if anchor_maf is None or anchor_maf <= 0:
        return []
    result: list[PairMetric] = []
    positions = [variant.position_cm for variant in variants]
    maximum_distance = distance_bins_cm[-1]
    lower_index = bisect_left(positions, anchor.position_cm - maximum_distance)
    upper_index = bisect_right(positions, anchor.position_cm + maximum_distance)
    for neighbor in variants[lower_index:upper_index]:
        if neighbor.variant_id == anchor.variant_id:
            continue
        signed_distance = neighbor.position_cm - anchor.position_cm
        label = distance_bin(abs(signed_distance), distance_bins_cm)
        neighbor_maf = neighbor.minor_allele_frequency
        if label is None or neighbor_maf is None or neighbor_maf < minimum_neighbor_maf:
            continue
        r2, called = r_squared(
            anchor.haplotypes, neighbor.haplotypes,
            minimum_called=minimum_called_haplotypes,
        )
        if r2 is None:
            continue
        result.append(PairMetric(
            anchor_id=anchor.variant_id,
            neighbor_id=neighbor.variant_id,
            chromosome=anchor.chromosome,
            anchor_position_bp=anchor.position_bp,
            neighbor_position_bp=neighbor.position_bp,
            anchor_position_cm=anchor.position_cm,
            neighbor_position_cm=neighbor.position_cm,
            signed_distance_cm=signed_distance,
            distance_bin=label,
            anchor_maf=anchor_maf,
            neighbor_maf=neighbor_maf,
            called_haplotype_count=called,
            r2=r2,
        ))
    return result


def local_background_metrics(
    variants: Sequence[PhasedVariant],
    *,
    excluded_variant_id: str,
    distance_bins_cm: Sequence[float],
    minimum_maf: float,
    minimum_called_haplotypes: int,
) -> dict[str, list[float]]:
    """Calcule le fond LD commun–commun dans la région cible."""

    result = {label: [] for label in bin_labels(distance_bins_cm)}
    eligible = [
        variant for variant in variants
        if variant.variant_id != excluded_variant_id
        and variant.minor_allele_frequency is not None
        and variant.minor_allele_frequency >= minimum_maf
    ]
    max_distance = distance_bins_cm[-1]
    for first_index, first in enumerate(eligible):
        for second in eligible[first_index + 1:]:
            distance = second.position_cm - first.position_cm
            if distance > max_distance:
                break
            label = distance_bin(abs(distance), distance_bins_cm)
            if label is None:
                continue
            r2, _ = r_squared(
                first.haplotypes, second.haplotypes,
                minimum_called=minimum_called_haplotypes,
            )
            if r2 is not None:
                result[label].append(r2)
    return result


def median_by_bin(
    metrics: Iterable[PairMetric], distance_bins_cm: Sequence[float],
) -> dict[str, float | None]:
    """Résume des paires par médiane dans chaque classe de distance."""

    grouped: dict[str, list[float]] = {label: [] for label in bin_labels(distance_bins_cm)}
    for metric in metrics:
        grouped[metric.distance_bin].append(metric.r2)
    return {
        label: statistics.median(values) if values else None
        for label, values in grouped.items()
    }


def count_neighbors_by_bin(
    anchor: PhasedVariant,
    variants: Sequence[PhasedVariant],
    *,
    distance_bins_cm: Sequence[float],
    minimum_neighbor_maf: float,
    minimum_called_haplotypes: int | None = None,
) -> dict[str, int]:
    """Compte les voisins éligibles et, si demandé, leur r² évaluable."""

    result = {label: 0 for label in bin_labels(distance_bins_cm)}
    if minimum_called_haplotypes is not None:
        metrics = centered_pair_metrics(
            anchor,
            variants,
            distance_bins_cm=distance_bins_cm,
            minimum_neighbor_maf=minimum_neighbor_maf,
            minimum_called_haplotypes=minimum_called_haplotypes,
        )
        for metric in metrics:
            result[metric.distance_bin] += 1
        return result

    positions = [variant.position_cm for variant in variants]
    maximum_distance = distance_bins_cm[-1]
    lower_index = bisect_left(positions, anchor.position_cm - maximum_distance)
    upper_index = bisect_right(positions, anchor.position_cm + maximum_distance)
    for neighbor in variants[lower_index:upper_index]:
        if neighbor.variant_id == anchor.variant_id:
            continue
        label = distance_bin(abs(neighbor.position_cm - anchor.position_cm), distance_bins_cm)
        maf = neighbor.minor_allele_frequency
        if label is not None and maf is not None and maf >= minimum_neighbor_maf:
            result[label] += 1
    return result


def density_is_matched(
    candidate_counts: dict[str, int], target_counts: dict[str, int], *, ratio: float,
) -> bool:
    """Vérifie un appariement préspécifié de densité par classe."""

    if ratio < 1:
        raise MatchedLdError("density_ratio_below_one")
    for label, target_count in target_counts.items():
        candidate_count = candidate_counts.get(label, 0)
        if target_count == 0:
            # Une classe sans paire cible reste NOT_EVALUATED. Elle ne doit pas
            # imposer artificiellement une absence de marqueur aux fenêtres de
            # référence, car aucune comparaison de r² n'y est possible.
            continue
        lower = max(1, math.floor(target_count / ratio))
        upper = math.ceil(target_count * ratio)
        if not lower <= candidate_count <= upper:
            return False
    return True
