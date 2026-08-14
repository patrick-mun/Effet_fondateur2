"""Statistique observée et règle exacte commune aux distributions nulles."""

from __future__ import annotations

from collections.abc import Sequence

from .model import HaplotypeProfile, Marker, SharingStatistic


class EnrichmentAnalysisError(ValueError):
    """Signale une entrée incompatible avec la méthode préspécifiée."""


def validate_independent_family_profiles(
    haplotypes: Sequence[HaplotypeProfile], minimum_units: int,
) -> None:
    """Vérifie les unités présélectionnées sans choisir de représentant."""
    if len(haplotypes) < minimum_units:
        raise EnrichmentAnalysisError("insufficient_preselected_family_units")
    family_ids = [profile.family_id for profile in haplotypes]
    if None in family_ids or len(set(family_ids)) != len(family_ids):
        raise EnrichmentAnalysisError("preselected_family_units_not_independent")


def distinct_background_count(
    haplotypes: Sequence[HaplotypeProfile], signature_indexes: Sequence[int],
) -> int:
    """Compte les signatures mutantes distinctes sans créer de réplications."""
    return len({tuple(profile.alleles[index] for index in signature_indexes) for profile in haplotypes})


def _validate_inputs(
    markers: Sequence[Marker], haplotypes: Sequence[HaplotypeProfile]
) -> int:
    if not markers or not haplotypes:
        raise EnrichmentAnalysisError("markers_or_haplotypes_empty")
    target_indexes = [index for index, marker in enumerate(markers) if marker.is_target]
    if len(target_indexes) != 1:
        raise EnrichmentAnalysisError("target_marker_missing_or_ambiguous")
    if any(len(profile.alleles) != len(markers) for profile in haplotypes):
        raise EnrichmentAnalysisError("haplotype_marker_count_mismatch")
    if any(markers[index].position_bp >= markers[index + 1].position_bp for index in range(len(markers) - 1)):
        raise EnrichmentAnalysisError("marker_bp_order_invalid")
    if any(markers[index].position_cm > markers[index + 1].position_cm for index in range(len(markers) - 1)):
        raise EnrichmentAnalysisError("marker_cm_order_invalid")
    return target_indexes[0]


def evaluate_exact_sharing(
    markers: Sequence[Marker],
    haplotypes: Sequence[HaplotypeProfile],
    *,
    minimum_flank_markers: int = 0,
) -> SharingStatistic:
    """Étend un IBS exact depuis la cible, exclue elle-même de la signature.

    Chaque bras s'arrête au premier manque ou désaccord. La fonction est la
    primitive unique utilisée pour l'observation et tous les tirages nuls.
    """
    if minimum_flank_markers < 0:
        raise EnrichmentAnalysisError("minimum_flank_markers_invalid")
    target_index = _validate_inputs(markers, haplotypes)

    def is_shared(marker_index: int) -> bool:
        alleles = {profile.alleles[marker_index] for profile in haplotypes}
        return None not in alleles and len(alleles) == 1

    left_index = target_index
    while left_index > 0 and is_shared(left_index - 1):
        left_index -= 1
    right_index = target_index
    while right_index + 1 < len(markers) and is_shared(right_index + 1):
        right_index += 1
    left_count = target_index - left_index
    right_count = right_index - target_index
    if left_count < minimum_flank_markers or right_count < minimum_flank_markers:
        return SharingStatistic(
            evaluation_status="NOT_EVALUATED", left_shared_cm=None,
            right_shared_cm=None, total_shared_cm=None, left_shared_bp=None,
            right_shared_bp=None, left_marker_count=left_count,
            right_marker_count=right_count, left_bound_index=left_index,
            right_bound_index=right_index,
            non_evaluable_reason="INSUFFICIENT_FLANK_MARKERS",
        )
    left_cm = markers[target_index].position_cm - markers[left_index].position_cm
    right_cm = markers[right_index].position_cm - markers[target_index].position_cm
    return SharingStatistic(
        evaluation_status="EVALUATED", left_shared_cm=left_cm,
        right_shared_cm=right_cm, total_shared_cm=left_cm + right_cm,
        left_shared_bp=markers[target_index].position_bp - markers[left_index].position_bp,
        right_shared_bp=markers[right_index].position_bp - markers[target_index].position_bp,
        left_marker_count=left_count, right_marker_count=right_count,
        left_bound_index=left_index, right_bound_index=right_index,
    )
