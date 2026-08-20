"""Énumération exhaustive du fond interne non porteur."""

from __future__ import annotations

from collections.abc import Sequence
from itertools import combinations

from .model import HaplotypeProfile, Marker, NullDraw
from .observed import EnrichmentAnalysisError, evaluate_exact_sharing


def enumerate_internal_null(
    markers: Sequence[Marker],
    haplotypes: Sequence[HaplotypeProfile],
    *,
    unit_count: int,
    minimum_flank_markers: int = 0,
) -> tuple[NullDraw, ...]:
    """Énumère chaque combinaison de copies provenant d'individus distincts."""
    if unit_count < 2:
        raise EnrichmentAnalysisError("null_unit_count_invalid")
    if len({profile.individual_id for profile in haplotypes}) < unit_count:
        return ()
    draws: list[NullDraw] = []
    attempt_index = 0
    for selected in combinations(haplotypes, unit_count):
        if len({profile.individual_id for profile in selected}) != unit_count:
            continue
        attempt_index += 1
        statistic = evaluate_exact_sharing(
            markers, selected, minimum_flank_markers=minimum_flank_markers
        )
        draws.append(NullDraw(
            source="INTERNAL", stratum="ALL", draw_index=len(draws) + 1,
            attempt_index=attempt_index,
            individual_ids=tuple(profile.individual_id for profile in selected),
            haplotype_ids=tuple(profile.haplotype_id for profile in selected),
            statistic=statistic,
        ))
    return tuple(draws)
