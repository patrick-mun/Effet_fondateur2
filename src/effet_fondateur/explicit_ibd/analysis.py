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


def evaluate_scenario(
    scenario: Scenario,
    segments: Sequence[IbdSegment],
    families: Sequence[str],
    mutant_haplotypes: Mapping[str, frozenset[str]],
    target_position_bp: int,
    boundary_tolerance_bp: int,
    control_frequency: float | None,
    maximum_control_frequency: float,
    calibration_specificity_acceptable: bool,
) -> ScenarioEvaluation:
    """Évalue les deux outils et toutes les paires avec une assignation globale."""
    scenario.validate()
    families = tuple(sorted(set(families)))
    required_pairs = set(combinations(families, 2))
    if len(families) < 3:
        return ScenarioEvaluation(scenario.scenario_id, False, False, False, False, False, None, None, (), "INSUFFICIENT_INDEPENDENT_FAMILIES")
    if not scenario.calibrated or not calibration_specificity_acceptable:
        return ScenarioEvaluation(scenario.scenario_id, False, False, False, False, False, None, None, (), "CALIBRATION_UNACCEPTABLE_OR_MISSING")
    if control_frequency is None:
        return ScenarioEvaluation(scenario.scenario_id, False, False, False, False, False, None, None, (), "CONTROL_BACKGROUND_NOT_EVALUATED")
    candidates = tuple(segment for segment in segments if segment.scenario_id == scenario.scenario_id and segment.contains(target_position_bp) and segment.length_cm >= scenario.minimum_cm and segment.marker_count >= scenario.minimum_markers)
    selected: tuple[IbdSegment, ...] = ()
    coherent = False
    for assignment in _assignments(families, mutant_haplotypes):
        current: list[IbdSegment] = []
        valid = True
        for tool in TOOLS:
            for pair in required_pairs:
                matches = [segment for segment in candidates if segment.tool == tool and segment.family_pair == pair and _mutant_compatible(segment, mutant_haplotypes, assignment)]
                if not matches:
                    valid = False
                    break
                current.append(sorted(matches, key=lambda item: (-item.length_cm, item.start_bp, item.end_bp))[0])
            if not valid:
                break
        if valid:
            selected, coherent = tuple(current), True
            break
    all_pairs = coherent and len(selected) == len(required_pairs) * len(TOOLS)
    if not all_pairs:
        return ScenarioEvaluation(scenario.scenario_id, True, False, False, False, control_frequency <= maximum_control_frequency, None, None, (), "PAIR_TOOL_OR_MUTANT_ASSIGNMENT_MISSING")
    starts, ends = [item.start_bp for item in selected], [item.end_bp for item in selected]
    boundaries = max(starts) - min(starts) <= boundary_tolerance_bp and max(ends) - min(ends) <= boundary_tolerance_bp
    common_start, common_end = max(starts), min(ends)
    if common_start > target_position_bp or common_end < target_position_bp:
        boundaries = False
    return ScenarioEvaluation(scenario.scenario_id, True, True, True, boundaries, control_frequency <= maximum_control_frequency, common_start, common_end, selected, None if boundaries else "BOUNDARIES_DISCORDANT")


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
