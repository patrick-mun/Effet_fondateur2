"""Dimensionnement D14b des fenêtres de simulation sur une carte génétique.

La demi-largeur primaire est le quantile qui borne la probabilité de censure
d'un bras sous le modèle exponentiel D14. Les conversions entre ``bp`` et
``cM`` utilisent uniquement une interpolation linéaire de la carte fournie ;
aucune extrapolation ni vitesse moyenne n'est autorisée.
"""

from __future__ import annotations

from collections.abc import Iterable, Mapping, Sequence
import csv
from dataclasses import dataclass
import math
import re


SIMULATION_WINDOW_METHOD_ID = "map_quantile_window_design_v1"
SIMULATION_WINDOW_METHOD_VERSION = "1.0.0"
_IDENTIFIER = re.compile(r"^[A-Za-z0-9][A-Za-z0-9._-]*$")
_BASES = frozenset("ACGT")


class SimulationWindowError(RuntimeError):
    """Signale une violation D14b avec un code stable et non sensible."""

    def __init__(self, code: str, detail: str | None = None) -> None:
        self.code = code
        self.detail = detail
        super().__init__(code if detail is None else f"{code}: {detail}")


@dataclass(frozen=True, order=True)
class WindowMapPoint:
    """Point de carte publique, exprimé en bp et cM."""

    chromosome: int
    position_bp: int
    position_cm: float


@dataclass(frozen=True, order=True)
class IbdMarkerRequirement:
    """Seuil observable exigé pour une longueur IBD centrée sur la cible."""

    scenario_id: str
    segment_length_cm: float
    minimum_observed_markers: int


@dataclass(frozen=True)
class WindowTarget:
    """SNV public présélectionné avant toute simulation."""

    variant_id: str
    position_bp: int
    reference_allele: str
    alternate_allele: str


@dataclass(frozen=True)
class WindowSpecification:
    """Demande primaire ou défi de censure pour un âge donné."""

    window_id: str
    purpose: str
    age_generations: int
    chromosome: int
    target: WindowTarget
    challenge_half_width_cm: float | None


@dataclass(frozen=True)
class SimulationWindowConfiguration:
    """Configuration finalisée et typée du dimensionnement D14b."""

    configuration_id: str
    panel_role: str
    maximum_arm_censor_probability: float
    minimum_flank_markers_each_side: int
    ibd_requirements: tuple[IbdMarkerRequirement, ...]
    window_specs: tuple[WindowSpecification, ...]


@dataclass(frozen=True)
class DesignedWindow:
    """Domaine physique réalisé depuis une largeur génétique préspecifiée."""

    specification: WindowSpecification
    start_bp: int
    end_bp: int
    start_cm: float
    target_cm: float
    end_cm: float
    requested_half_width_cm: float
    theoretical_arm_censor_probability: float

    @property
    def realized_left_width_cm(self) -> float:
        return self.target_cm - self.start_cm

    @property
    def realized_right_width_cm(self) -> float:
        return self.end_cm - self.target_cm


def _fail(code: str, detail: str | None = None) -> None:
    raise SimulationWindowError(code, detail)


def _identifier(value: object, namespace: str) -> str:
    if not isinstance(value, str) or _IDENTIFIER.fullmatch(value) is None:
        _fail("simulation_window_identifier_invalid", namespace)
    return value


def minimum_half_width_cm(
    age_generations: int,
    maximum_arm_censor_probability: float,
) -> float:
    """Calcule ``-100 ln(p) / g`` en cM pour le modèle D14.

    ``g`` est l'âge en générations et ``p`` la probabilité maximale de censure
    d'un bras. Les booléens sont refusés malgré leur héritage Python de ``int``.
    """
    if type(age_generations) is not int or age_generations < 1:
        _fail("simulation_window_age_invalid")
    if (
        isinstance(maximum_arm_censor_probability, bool)
        or not isinstance(maximum_arm_censor_probability, (int, float))
        or not math.isfinite(float(maximum_arm_censor_probability))
        or not 0 < float(maximum_arm_censor_probability) < 1
    ):
        _fail("simulation_window_censor_probability_invalid")
    return -100.0 * math.log(float(maximum_arm_censor_probability)) / age_generations


def parse_window_configuration(
    document: Mapping[str, object],
) -> SimulationWindowConfiguration:
    """Valide les invariants croisés non exprimables simplement par le schéma."""
    if set(document) != {
        "schema_version",
        "simulation_window_config_id",
        "config_state",
        "panel_role",
        "method",
        "marker_requirements",
        "window_specs",
    }:
        _fail("simulation_window_config_invalid")
    if (
        document.get("schema_version") != "1.0.0"
        or document.get("config_state") != "FINALIZED"
        or document.get("panel_role") not in {"CALIBRATION", "VALIDATION"}
    ):
        _fail("simulation_window_config_invalid")
    configuration_id = _identifier(
        document.get("simulation_window_config_id"), "configuration"
    )
    method = document.get("method")
    marker_requirements = document.get("marker_requirements")
    raw_specs = document.get("window_specs")
    if (
        not isinstance(method, Mapping)
        or not isinstance(marker_requirements, Mapping)
        or not isinstance(raw_specs, list)
    ):
        _fail("simulation_window_config_invalid")
    if set(method) != {
        "method_id",
        "method_version",
        "maximum_arm_censor_probability",
        "genetic_map_interpolation",
        "candidate_selection",
        "outcome_filtering",
    } or (
        method.get("method_id") != SIMULATION_WINDOW_METHOD_ID
        or method.get("method_version") != SIMULATION_WINDOW_METHOD_VERSION
        or method.get("genetic_map_interpolation")
        != "PIECEWISE_LINEAR_NO_EXTRAPOLATION"
        or method.get("candidate_selection")
        != "PRESELECTED_PUBLIC_REFERENCE_TARGET"
        or method.get("outcome_filtering") != "FORBIDDEN"
    ):
        _fail("simulation_window_method_invalid")
    probability = method.get("maximum_arm_censor_probability")
    if (
        isinstance(probability, bool)
        or not isinstance(probability, (int, float))
        or not math.isfinite(float(probability))
        or not 0 < float(probability) <= 0.05
    ):
        _fail("simulation_window_censor_probability_invalid")
    if set(marker_requirements) != {
        "minimum_flank_markers_each_side",
        "ibd_scenarios",
    }:
        _fail("simulation_window_marker_requirements_invalid")
    minimum_flanks = marker_requirements.get("minimum_flank_markers_each_side")
    raw_scenarios = marker_requirements.get("ibd_scenarios")
    if type(minimum_flanks) is not int or minimum_flanks < 1:
        _fail("simulation_window_flank_requirement_invalid")
    if not isinstance(raw_scenarios, list) or not raw_scenarios:
        _fail("simulation_window_ibd_requirements_invalid")
    scenarios: list[IbdMarkerRequirement] = []
    scenario_ids: set[str] = set()
    segment_lengths: set[float] = set()
    for raw_scenario in raw_scenarios:
        if not isinstance(raw_scenario, Mapping) or set(raw_scenario) != {
            "scenario_id",
            "segment_length_cm",
            "minimum_observed_markers",
        }:
            _fail("simulation_window_ibd_requirement_invalid")
        scenario_id = _identifier(raw_scenario.get("scenario_id"), "scenario")
        length = raw_scenario.get("segment_length_cm")
        minimum_markers = raw_scenario.get("minimum_observed_markers")
        if (
            isinstance(length, bool)
            or not isinstance(length, (int, float))
            or float(length) not in {1.0, 1.5, 2.0, 3.0, 5.0}
            or type(minimum_markers) is not int
            or minimum_markers < 1
            or scenario_id in scenario_ids
            or float(length) in segment_lengths
        ):
            _fail("simulation_window_ibd_requirement_invalid")
        scenario_ids.add(scenario_id)
        segment_lengths.add(float(length))
        scenarios.append(
            IbdMarkerRequirement(scenario_id, float(length), minimum_markers)
        )
    specs: list[WindowSpecification] = []
    spec_ids: set[str] = set()
    purposes: set[str] = set()
    chromosomes: set[int] = set()
    for raw_spec in raw_specs:
        if not isinstance(raw_spec, Mapping):
            _fail("simulation_window_spec_invalid")
        allowed = {
            "window_id",
            "purpose",
            "age_generations",
            "chromosome",
            "target",
            "challenge_half_width_cm",
        }
        if not set(raw_spec) <= allowed or not {
            "window_id",
            "purpose",
            "age_generations",
            "chromosome",
            "target",
        } <= set(raw_spec):
            _fail("simulation_window_spec_invalid")
        window_id = _identifier(raw_spec.get("window_id"), "window")
        purpose = raw_spec.get("purpose")
        age = raw_spec.get("age_generations")
        chromosome = raw_spec.get("chromosome")
        target = raw_spec.get("target")
        challenge_width = raw_spec.get("challenge_half_width_cm")
        if (
            window_id in spec_ids
            or purpose not in {"PRIMARY", "CENSORING_CHALLENGE"}
            or type(age) is not int
            or age < 1
            or type(chromosome) is not int
            or not 1 <= chromosome <= 22
            or not isinstance(target, Mapping)
        ):
            _fail("simulation_window_spec_invalid")
        if set(target) != {
            "variant_id",
            "position_bp",
            "reference_allele",
            "alternate_allele",
        }:
            _fail("simulation_window_target_invalid")
        variant_id = _identifier(target.get("variant_id"), "target")
        position_bp = target.get("position_bp")
        reference = target.get("reference_allele")
        alternate = target.get("alternate_allele")
        if (
            type(position_bp) is not int
            or position_bp < 1
            or reference not in _BASES
            or alternate not in _BASES
            or reference == alternate
        ):
            _fail("simulation_window_target_invalid")
        if purpose == "PRIMARY" and "challenge_half_width_cm" in raw_spec:
            _fail("simulation_window_primary_width_override_forbidden")
        if purpose == "CENSORING_CHALLENGE":
            if (
                isinstance(challenge_width, bool)
                or not isinstance(challenge_width, (int, float))
                or not math.isfinite(float(challenge_width))
                or float(challenge_width) <= 0
            ):
                _fail("simulation_window_challenge_width_invalid")
            challenge_width = float(challenge_width)
        else:
            challenge_width = None
        spec_ids.add(window_id)
        purposes.add(str(purpose))
        chromosomes.add(chromosome)
        specs.append(
            WindowSpecification(
                window_id,
                str(purpose),
                age,
                chromosome,
                WindowTarget(variant_id, position_bp, str(reference), str(alternate)),
                challenge_width,
            )
        )
    if purposes != {"PRIMARY", "CENSORING_CHALLENGE"} or len(chromosomes) != 1:
        _fail("simulation_window_spec_set_invalid")
    return SimulationWindowConfiguration(
        configuration_id=configuration_id,
        panel_role=str(document["panel_role"]),
        maximum_arm_censor_probability=float(probability),
        minimum_flank_markers_each_side=minimum_flanks,
        ibd_requirements=tuple(sorted(scenarios, key=lambda item: item.segment_length_cm)),
        window_specs=tuple(specs),
    )


def parse_public_genetic_map(
    lines: Iterable[str],
    *,
    expected_chromosome: int,
) -> tuple[WindowMapPoint, ...]:
    """Lit la carte D09 normalisée et impose la monotonie bp/cM."""
    reader = csv.DictReader(lines, delimiter="\t")
    required = {"ASSEMBLY", "CHROMOSOME", "POSITION_BP", "POSITION_CM"}
    if reader.fieldnames is None or not required <= set(reader.fieldnames):
        _fail("simulation_window_map_columns_invalid")
    points: list[WindowMapPoint] = []
    for row in reader:
        try:
            chromosome = int(row["CHROMOSOME"])
            position_bp = int(row["POSITION_BP"])
            position_cm = float(row["POSITION_CM"])
        except (TypeError, ValueError):
            _fail("simulation_window_map_row_invalid")
        if (
            row["ASSEMBLY"] != "GRCh38"
            or chromosome != expected_chromosome
            or position_bp < 1
            or not math.isfinite(position_cm)
            or position_cm < 0
            or (
                points
                and (
                    position_bp <= points[-1].position_bp
                    or position_cm < points[-1].position_cm
                )
            )
        ):
            _fail("simulation_window_map_not_monotonic")
        points.append(WindowMapPoint(chromosome, position_bp, position_cm))
    if len(points) < 2:
        _fail("simulation_window_map_insufficient")
    return tuple(points)


def interpolate_window_cm(
    points: Sequence[WindowMapPoint], position_bp: int
) -> float:
    """Interpole une coordonnée physique couverte, sans extrapolation."""
    if not points or position_bp < points[0].position_bp or position_bp > points[-1].position_bp:
        _fail("simulation_window_map_extrapolation_forbidden")
    for left, right in zip(points, points[1:], strict=False):
        if left.position_bp <= position_bp <= right.position_bp:
            fraction = (position_bp - left.position_bp) / (
                right.position_bp - left.position_bp
            )
            return left.position_cm + fraction * (
                right.position_cm - left.position_cm
            )
    return points[-1].position_cm


def interpolate_window_bp(
    points: Sequence[WindowMapPoint],
    position_cm: float,
    *,
    boundary: str,
) -> int:
    """Inverse la carte et arrondit vers l'extérieur du domaine demandé."""
    if boundary not in {"LOWER", "UPPER"}:
        _fail("simulation_window_boundary_invalid")
    if (
        not points
        or position_cm < points[0].position_cm
        or position_cm > points[-1].position_cm
    ):
        _fail("simulation_window_map_extrapolation_forbidden")
    candidates: list[int] = []
    for left, right in zip(points, points[1:], strict=False):
        if left.position_cm <= position_cm <= right.position_cm:
            if right.position_cm == left.position_cm:
                candidates.extend((left.position_bp, right.position_bp))
                continue
            fraction = (position_cm - left.position_cm) / (
                right.position_cm - left.position_cm
            )
            raw_position = left.position_bp + fraction * (
                right.position_bp - left.position_bp
            )
            candidates.append(
                math.floor(raw_position) if boundary == "LOWER" else math.ceil(raw_position)
            )
    if not candidates:
        _fail("simulation_window_map_inverse_failed")
    return min(candidates) if boundary == "LOWER" else max(candidates)


def design_window(
    specification: WindowSpecification,
    map_points: Sequence[WindowMapPoint],
    *,
    maximum_arm_censor_probability: float,
) -> DesignedWindow:
    """Matérialise les bornes d'une spécification et contrôle sa faisabilité."""
    if any(point.chromosome != specification.chromosome for point in map_points):
        _fail("simulation_window_map_chromosome_mismatch")
    target_cm = interpolate_window_cm(map_points, specification.target.position_bp)
    primary_width = minimum_half_width_cm(
        specification.age_generations, maximum_arm_censor_probability
    )
    if specification.purpose == "PRIMARY":
        requested_width = primary_width
    elif specification.purpose == "CENSORING_CHALLENGE":
        if specification.challenge_half_width_cm is None:
            _fail("simulation_window_challenge_width_invalid")
        requested_width = specification.challenge_half_width_cm
        if requested_width >= primary_width:
            _fail("simulation_window_challenge_not_truncated")
    else:
        _fail("simulation_window_purpose_invalid")
    requested_start_cm = target_cm - requested_width
    requested_end_cm = target_cm + requested_width
    if (
        requested_start_cm < map_points[0].position_cm
        or requested_end_cm > map_points[-1].position_cm
    ):
        _fail("simulation_window_infeasible_window", specification.window_id)
    start_bp = interpolate_window_bp(
        map_points, requested_start_cm, boundary="LOWER"
    )
    end_bp = interpolate_window_bp(map_points, requested_end_cm, boundary="UPPER")
    start_cm = interpolate_window_cm(map_points, start_bp)
    end_cm = interpolate_window_cm(map_points, end_bp)
    tolerance = 1e-9
    if (
        start_bp >= specification.target.position_bp
        or end_bp <= specification.target.position_bp
        or target_cm - start_cm + tolerance < requested_width
        or end_cm - target_cm + tolerance < requested_width
    ):
        _fail("simulation_window_realized_width_insufficient")
    censor_probability = math.exp(
        -specification.age_generations * requested_width / 100.0
    )
    if (
        specification.purpose == "PRIMARY"
        and censor_probability > maximum_arm_censor_probability + tolerance
    ):
        _fail("simulation_window_primary_censor_probability_exceeded")
    return DesignedWindow(
        specification=specification,
        start_bp=start_bp,
        end_bp=end_bp,
        start_cm=start_cm,
        target_cm=target_cm,
        end_cm=end_cm,
        requested_half_width_cm=requested_width,
        theoretical_arm_censor_probability=censor_probability,
    )
