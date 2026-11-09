"""Noyau scientifique pur de 16C.

Les fonctions de ce module ne lancent aucun outil et ne connaissent aucun gène.
Elles normalisent les appels, imposent une attribution haplotypique globale et
séparent strictement scénario primaire, sensibilités et non-évaluation.
"""

from __future__ import annotations

import csv
import gzip
from dataclasses import dataclass
from itertools import combinations, product
from pathlib import Path
from typing import Iterable, Mapping, Sequence


TOOLS = ("HAP_IBD", "REFINED_IBD")


class ExplicitIbdError(ValueError):
    """Signale une entrée IBD ou une règle scientifique incohérente."""


@dataclass(frozen=True)
class Scenario:
    """Seuils préspécifiés d'un appel IBD."""

    scenario_id: str
    role: str
    minimum_cm: float
    minimum_markers: int
    calibrated: bool
    calibration_specificity_acceptable: bool
    calibration_control_frequency: float | None

    def validate(self) -> None:
        if self.role not in {"PRIMARY", "SENSITIVITY"}:
            raise ExplicitIbdError("scenario_role_invalid")
        if self.minimum_cm <= 0 or self.minimum_markers < 1:
            raise ExplicitIbdError("scenario_threshold_invalid")
        if self.role == "PRIMARY" and (
            self.minimum_cm < 2.0 or self.minimum_markers < 100
        ):
            raise ExplicitIbdError("primary_threshold_post_hoc_or_too_permissive")
        if self.role == "SENSITIVITY" and self.minimum_cm < 1.0:
            raise ExplicitIbdError("sensitivity_threshold_too_permissive")
        if self.calibration_control_frequency is not None and not (
            0.0 <= self.calibration_control_frequency <= 1.0
        ):
            raise ExplicitIbdError("calibration_control_frequency_invalid")
        if self.calibrated and self.calibration_control_frequency is None:
            raise ExplicitIbdError("calibration_control_frequency_missing")
        if not self.calibrated and (
            self.calibration_specificity_acceptable
            or self.calibration_control_frequency is not None
        ):
            raise ExplicitIbdError("uncalibrated_scenario_has_calibration_results")

    def calibration_is_acceptable(self, maximum_control_frequency: float) -> bool:
        """Vérifie la calibration externe préspécifiée du scénario."""
        return (
            self.calibrated
            and self.calibration_specificity_acceptable
            and self.calibration_control_frequency is not None
            and self.calibration_control_frequency <= maximum_control_frequency
        )


@dataclass(frozen=True)
class IbdSegment:
    """Segment normalisé, bornes inclusives en bp et cM."""

    tool: str
    scenario_id: str
    sample_1: str
    haplotype_1: str
    family_1: str
    sample_2: str
    haplotype_2: str
    family_2: str
    chromosome: str
    start_bp: int
    end_bp: int
    start_cm: float
    end_cm: float
    marker_count: int

    @property
    def length_cm(self) -> float:
        return self.end_cm - self.start_cm

    @property
    def family_pair(self) -> tuple[str, str]:
        return tuple(sorted((self.family_1, self.family_2)))

    def contains(self, position_bp: int) -> bool:
        return self.start_bp <= position_bp <= self.end_bp

    def contains_target(self, chromosome: str | int, position_bp: int) -> bool:
        """Indique si le segment contient une cible sur le même chromosome."""
        target_chromosome = str(chromosome).removeprefix("chr")
        return self.chromosome == target_chromosome and self.contains(position_bp)


@dataclass(frozen=True)
class ScenarioEvaluation:
    scenario_id: str
    evaluable: bool
    all_pairs_by_both_tools: bool
    globally_coherent: bool
    boundaries_concordant: bool
    control_frequency_acceptable: bool
    common_start_bp: int | None
    common_end_bp: int | None
    selected_segments: tuple[IbdSegment, ...]
    reason: str | None


def _open_text(path: Path):
    return gzip.open(path, "rt", encoding="utf-8", newline="") if path.suffix == ".gz" else path.open("r", encoding="utf-8", newline="")


def _hap(value: str) -> str:
    normalized = value.strip().upper()
    if normalized in {"1", "H1"}:
        return "H1"
    if normalized in {"2", "H2"}:
        return "H2"
    raise ExplicitIbdError("ibd_haplotype_invalid")


def _segment(
    tool: str,
    scenario_id: str,
    fields: Sequence[str],
    family_by_sample: Mapping[str, str],
    cm_at_bp: Mapping[int, float],
    marker_positions: Sequence[int],
) -> IbdSegment:
    if len(fields) < 7:
        raise ExplicitIbdError(f"{tool.lower()}_row_invalid")
    sample_1, hap_1, sample_2, hap_2, chromosome = fields[:5]
    try:
        start_bp, end_bp = int(fields[5]), int(fields[6])
        family_1, family_2 = family_by_sample[sample_1], family_by_sample[sample_2]
        start_cm, end_cm = cm_at_bp[start_bp], cm_at_bp[end_bp]
    except (KeyError, ValueError) as error:
        raise ExplicitIbdError(f"{tool.lower()}_coordinate_or_sample_invalid") from error
    if start_bp >= end_bp or start_cm > end_cm:
        raise ExplicitIbdError(f"{tool.lower()}_segment_invalid")
    marker_count = sum(start_bp <= position <= end_bp for position in marker_positions)
    return IbdSegment(tool, scenario_id, sample_1, _hap(hap_1), family_1, sample_2, _hap(hap_2), family_2, chromosome.removeprefix("chr"), start_bp, end_bp, start_cm, end_cm, marker_count)


def _parse(
    path: Path,
    tool: str,
    scenario_id: str,
    family_by_sample: Mapping[str, str],
    cm_at_bp: Mapping[int, float],
    marker_positions: Sequence[int],
) -> tuple[IbdSegment, ...]:
    """Lit les sept premières colonnes officielles; commentaires/en-tête permis."""
    result: list[IbdSegment] = []
    with _open_text(path) as handle:
        for raw in handle:
            line = raw.strip()
            if not line or line.startswith("#"):
                continue
            fields = next(csv.reader([line], delimiter="\t"))
            if fields[0].lower() in {"id1", "sample_1", "sample1"}:
                continue
            result.append(_segment(tool, scenario_id, fields, family_by_sample, cm_at_bp, marker_positions))
    return tuple(result)


def parse_hap_ibd(path: Path, scenario_id: str, family_by_sample: Mapping[str, str], cm_at_bp: Mapping[int, float], marker_positions: Sequence[int]) -> tuple[IbdSegment, ...]:
    """Normalise une sortie Hap-IBD `.ibd[.gz]`."""
    return _parse(path, "HAP_IBD", scenario_id, family_by_sample, cm_at_bp, marker_positions)


def parse_refined_ibd(path: Path, scenario_id: str, family_by_sample: Mapping[str, str], cm_at_bp: Mapping[int, float], marker_positions: Sequence[int]) -> tuple[IbdSegment, ...]:
    """Normalise une sortie Refined IBD `.ibd[.gz]`."""
    return _parse(path, "REFINED_IBD", scenario_id, family_by_sample, cm_at_bp, marker_positions)


def _mutant_compatible(segment: IbdSegment, mutant_haplotypes: Mapping[str, frozenset[str]], assignment: Mapping[str, str]) -> bool:
    observed = {segment.family_1: segment.haplotype_1, segment.family_2: segment.haplotype_2}
    return all(observed[family] == assignment[family] and observed[family] in mutant_haplotypes.get(family, frozenset()) for family in observed)


def _assignments(families: Sequence[str], mutant_haplotypes: Mapping[str, frozenset[str]]) -> Iterable[dict[str, str]]:
    choices = [tuple(sorted(mutant_haplotypes.get(family, frozenset()))) for family in families]
    if any(not choice for choice in choices):
        return ()
    return (dict(zip(families, values, strict=True)) for values in product(*choices))


def boundaries_are_concordant(
    segments: Sequence[IbdSegment],
    tolerance_bp: int,
) -> bool:
    """Compare les bornes d'un ensemble de segments avec une tolérance en pb."""
    if not segments or tolerance_bp < 0:
        return False
    starts = [segment.start_bp for segment in segments]
    ends = [segment.end_bp for segment in segments]
    return (
        max(starts) - min(starts) <= tolerance_bp
        and max(ends) - min(ends) <= tolerance_bp
    )


def _selection_key(segments: Sequence[IbdSegment]) -> tuple[object, ...]:
    """Ordonne les sélections de façon stable en privilégiant leur longueur."""
    return (
        -sum(segment.length_cm for segment in segments),
        tuple(
            (
                segment.tool,
                segment.family_pair,
                segment.start_bp,
                segment.end_bp,
                segment.sample_1,
                segment.haplotype_1,
                segment.sample_2,
                segment.haplotype_2,
            )
            for segment in segments
        ),
    )


def _select_for_assignment(
    match_groups: Sequence[Sequence[IbdSegment]],
    boundary_tolerance_bp: int,
) -> tuple[tuple[IbdSegment, ...], bool]:
    """Retient la meilleure combinaison concordante pour une assignation.

    Chaque borne minimale d'une combinaison réalisable appartient à un segment
    candidat. Les couples de bornes candidates permettent donc de chercher une
    combinaison sans énumérer le produit cartésien de tous les segments.
    """
    fallback = tuple(
        sorted(
            matches,
            key=lambda item: (-item.length_cm, item.start_bp, item.end_bp),
        )[0]
        for matches in match_groups
    )
    all_matches = tuple(segment for matches in match_groups for segment in matches)
    best: tuple[IbdSegment, ...] | None = None
    for minimum_start in sorted({segment.start_bp for segment in all_matches}):
        maximum_start = minimum_start + boundary_tolerance_bp
        start_filtered = tuple(
            tuple(
                segment
                for segment in matches
                if minimum_start <= segment.start_bp <= maximum_start
            )
            for matches in match_groups
        )
        if any(not matches for matches in start_filtered):
            continue
        possible_minimum_ends = sorted(
            {segment.end_bp for matches in start_filtered for segment in matches}
        )
        for minimum_end in possible_minimum_ends:
            maximum_end = minimum_end + boundary_tolerance_bp
            eligible = tuple(
                tuple(
                    segment
                    for segment in matches
                    if minimum_end <= segment.end_bp <= maximum_end
                )
                for matches in start_filtered
            )
            if any(not matches for matches in eligible):
                continue
            selection = tuple(
                sorted(
                    matches,
                    key=lambda item: (
                        -item.length_cm,
                        item.start_bp,
                        item.end_bp,
                    ),
                )[0]
                for matches in eligible
            )
            if best is None or _selection_key(selection) < _selection_key(best):
                best = selection
    return (best, True) if best is not None else (fallback, False)


def evaluate_scenario(
    scenario: Scenario,
    segments: Sequence[IbdSegment],
    families: Sequence[str],
    mutant_haplotypes: Mapping[str, frozenset[str]],
    target_chromosome: str | int,
    target_position_bp: int,
    boundary_tolerance_bp: int,
    empirical_control_frequency: float | None,
    maximum_control_frequency: float,
    control_pair_count_sufficient: bool,
) -> ScenarioEvaluation:
    """Évalue les deux outils et toutes les paires avec une assignation globale."""
    scenario.validate()
    families = tuple(sorted(set(families)))
    required_pairs = tuple(combinations(families, 2))
    if len(families) < 3:
        return ScenarioEvaluation(scenario.scenario_id, False, False, False, False, False, None, None, (), "INSUFFICIENT_INDEPENDENT_FAMILIES")
    if not scenario.calibration_is_acceptable(maximum_control_frequency):
        return ScenarioEvaluation(scenario.scenario_id, False, False, False, False, False, None, None, (), "CALIBRATION_UNACCEPTABLE_OR_MISSING")
    if not control_pair_count_sufficient:
        return ScenarioEvaluation(scenario.scenario_id, False, False, False, False, False, None, None, (), "CONTROL_BACKGROUND_INSUFFICIENT")
    if empirical_control_frequency is None:
        return ScenarioEvaluation(scenario.scenario_id, False, False, False, False, False, None, None, (), "CONTROL_BACKGROUND_NOT_EVALUATED")
    candidates = tuple(
        segment
        for segment in segments
        if segment.scenario_id == scenario.scenario_id
        and segment.contains_target(target_chromosome, target_position_bp)
        and segment.length_cm >= scenario.minimum_cm
        and segment.marker_count >= scenario.minimum_markers
    )
    best_concordant: tuple[IbdSegment, ...] | None = None
    best_complete: tuple[IbdSegment, ...] | None = None
    for assignment in _assignments(families, mutant_haplotypes):
        match_groups: list[tuple[IbdSegment, ...]] = []
        complete_assignment = True
        for tool in TOOLS:
            for pair in required_pairs:
                matches = tuple(
                    segment
                    for segment in candidates
                    if segment.tool == tool
                    and segment.family_pair == pair
                    and _mutant_compatible(
                        segment,
                        mutant_haplotypes,
                        assignment,
                    )
                )
                if not matches:
                    complete_assignment = False
                    break
                match_groups.append(matches)
            if not complete_assignment:
                break
        if not complete_assignment:
            continue
        current, concordant = _select_for_assignment(
            match_groups,
            boundary_tolerance_bp,
        )
        if best_complete is None or _selection_key(current) < _selection_key(
            best_complete
        ):
            best_complete = current
        if concordant and (
            best_concordant is None
            or _selection_key(current) < _selection_key(best_concordant)
        ):
            best_concordant = current
    selected = best_concordant or best_complete
    if selected is None:
        return ScenarioEvaluation(scenario.scenario_id, True, False, False, False, empirical_control_frequency <= maximum_control_frequency, None, None, (), "PAIR_TOOL_OR_MUTANT_ASSIGNMENT_MISSING")
    starts, ends = [item.start_bp for item in selected], [item.end_bp for item in selected]
    boundaries = boundaries_are_concordant(selected, boundary_tolerance_bp)
    common_start, common_end = max(starts), min(ends)
    if common_start > target_position_bp or common_end < target_position_bp:
        raise ExplicitIbdError("target_outside_common_interval_invariant")
    return ScenarioEvaluation(scenario.scenario_id, True, True, True, boundaries, empirical_control_frequency <= maximum_control_frequency, common_start, common_end, selected, None if boundaries else "BOUNDARIES_DISPERSED")


def classify_explicit_ibd(primary: ScenarioEvaluation, sensitivities: Sequence[ScenarioEvaluation]) -> str:
    """Classe 16C sans jamais déclarer l'effet fondateur comme prouvé."""
    if not primary.evaluable:
        return "NOT_EVALUABLE"
    if primary.all_pairs_by_both_tools and primary.globally_coherent and primary.boundaries_concordant and primary.control_frequency_acceptable:
        return "PRIMARY_CONCORDANT_IBD_SUPPORT"
    if any(item.evaluable and item.all_pairs_by_both_tools and item.globally_coherent and item.boundaries_concordant and item.control_frequency_acceptable for item in sensitivities):
        return "SENSITIVITY_ONLY_IBD_SUPPORT"
    if primary.all_pairs_by_both_tools and not primary.boundaries_concordant:
        return "METHOD_DISCORDANT"
    if primary.reason == "PAIR_TOOL_OR_MUTANT_ASSIGNMENT_MISSING":
        return "NO_PRIMARY_IBD_CALL"
    return "NO_PRIMARY_IBD_CALL"
