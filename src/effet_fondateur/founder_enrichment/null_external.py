"""Monte-Carlo reproductible sur les haplotypes externes."""

from __future__ import annotations

import random
from collections import defaultdict
from collections.abc import Sequence

from .model import HaplotypeProfile, Marker, NullDraw
from .observed import EnrichmentAnalysisError, evaluate_exact_sharing


def sample_external_null(
    markers: Sequence[Marker],
    haplotypes: Sequence[HaplotypeProfile],
    *,
    unit_count: int,
    evaluable_draws: int,
    max_attempts: int,
    random_seed: int,
    minimum_flank_markers: int = 0,
    stratum: str = "ALL",
) -> tuple[NullDraw, ...]:
    """Tire des individus sans remise puis une seule copie par individu."""
    if unit_count < 2 or evaluable_draws < 1 or max_attempts < evaluable_draws:
        raise EnrichmentAnalysisError("external_null_parameters_invalid")
    by_individual: dict[str, list[HaplotypeProfile]] = defaultdict(list)
    for profile in haplotypes:
        by_individual[profile.individual_id].append(profile)
    individual_ids = sorted(by_individual)
    if len(individual_ids) < unit_count:
        return ()
    generator = random.Random(random_seed)
    draws: list[NullDraw] = []
    evaluated = 0
    for attempt_index in range(1, max_attempts + 1):
        selected_individuals = tuple(generator.sample(individual_ids, unit_count))
        selected = tuple(generator.choice(by_individual[individual_id]) for individual_id in selected_individuals)
        statistic = evaluate_exact_sharing(
            markers, selected, minimum_flank_markers=minimum_flank_markers
        )
        draws.append(NullDraw(
            source="EXTERNAL", stratum=stratum, draw_index=len(draws) + 1,
            attempt_index=attempt_index, individual_ids=selected_individuals,
            haplotype_ids=tuple(profile.haplotype_id for profile in selected),
            statistic=statistic,
        ))
        if statistic.evaluation_status == "EVALUATED":
            evaluated += 1
            if evaluated == evaluable_draws:
                break
    return tuple(draws)
