"""Simulation D14 d'un haplotype fondateur autour d'un variant synthétique.

Les limites gauche et droite sont les premiers événements d'un processus de
recombinaison exponentiel sur la carte génétique. Une limite hors domaine est
censurée au bord et n'est jamais présentée comme une recombinaison observée.
Le segment mosaïqué est appliqué à une origine haplotypique D12 complète afin
de préserver les transmissions mendéliennes déjà établies par D13.
"""

from __future__ import annotations

from collections import defaultdict
from collections.abc import Iterable, Mapping, Sequence
import csv
from dataclasses import dataclass
import io
import math
import random
import re
from typing import TextIO


FOUNDER_RECOMBINATION_METHOD_ID = "map_first_breakpoint_founder_mosaic_v1"
FOUNDER_RECOMBINATION_METHOD_VERSION = "1.0.0"
_IDENTIFIER = re.compile(r"^[A-Za-z0-9][A-Za-z0-9._-]*$")


class FounderRecombinationError(RuntimeError):
    """Signale une violation D14 avec un code stable et non sensible."""

    def __init__(self, code: str, detail: str | None = None) -> None:
        self.code = code
        self.detail = detail
        super().__init__(code if detail is None else f"{code}: {detail}")


@dataclass(frozen=True, order=True)
class GeneticMapPoint:
    """Coordonnée physique et génétique d'un marqueur, en bp et cM."""

    position_bp: int
    position_cm: float


@dataclass(frozen=True, order=True)
class HaplotypeOrigin:
    """Origine D12 stable d'une copie haplotypique D13."""

    d12_sample_id: str
    d12_haplotype_index: int


@dataclass(frozen=True)
class PedigreeMember:
    """Membre D13 nécessaire à la sélection et à la propagation D14."""

    sample_id: str
    family_unit_id: str
    generation: int
    is_founder: bool


@dataclass(frozen=True)
class FounderSimulationConfiguration:
    """Paramètres scientifiques finalisés d'un run D14."""

    configuration_id: str
    random_seed: int
    replication_id: str
    age_generations: int
    target_variant_id: str
    chromosome: int
    target_position_bp: int
    reference_allele: str
    alternate_allele: str


@dataclass(frozen=True)
class FounderFamilySegment:
    """Segment fondateur réalisé pour une famille indépendante."""

    family_unit_id: str
    carrier_origin: HaplotypeOrigin
    left_boundary_bp: int
    right_boundary_bp: int
    left_boundary_cm: float
    target_position_cm: float
    right_boundary_cm: float
    left_censored: bool
    right_censored: bool
    recombination_event_count: int

    @property
    def left_length_cm(self) -> float:
        return self.target_position_cm - self.left_boundary_cm

    @property
    def right_length_cm(self) -> float:
        return self.right_boundary_cm - self.target_position_cm

    @property
    def total_length_cm(self) -> float:
        return self.right_boundary_cm - self.left_boundary_cm


@dataclass(frozen=True, order=True)
class CarrierHaplotype:
    """Copie présente portant la mutation et le segment fondateur."""

    family_unit_id: str
    sample_id: str
    haplotype_index: int
    carrier_origin: HaplotypeOrigin


@dataclass(frozen=True)
class FounderSimulationPlan:
    """Vérité D14 indépendante des génotypes observés."""

    configuration: FounderSimulationConfiguration
    founder_haplotype_id: str
    founder_template_origin: HaplotypeOrigin
    founder_source_stratum_id: str
    segments: tuple[FounderFamilySegment, ...]
    carriers: tuple[CarrierHaplotype, ...]
    sample_origins: Mapping[str, tuple[HaplotypeOrigin, HaplotypeOrigin]]
    background_strata: Mapping[HaplotypeOrigin, str]


@dataclass(frozen=True)
class FounderVcfStats:
    """Contrôles agrégés de la transformation VCF D14."""

    dense_variant_count: int
    observed_variant_count: int
    carrier_sample_count: int
    carrier_haplotype_count: int
    changed_haplotype_allele_count: int
    mendelian_error_count: int
    target_variant_count: int


def _fail(code: str, detail: str | None = None) -> None:
    raise FounderRecombinationError(code, detail)


def _identifier(value: object, namespace: str) -> str:
    if not isinstance(value, str) or _IDENTIFIER.fullmatch(value) is None:
        _fail("founder_identifier_invalid", namespace)
    return value


def parse_founder_configuration(
    document: Mapping[str, object],
) -> FounderSimulationConfiguration:
    """Valide le contrat finalisé D14 et retourne ses paramètres typés."""
    if set(document) != {
        "schema_version",
        "founder_simulation_config_id",
        "config_state",
        "replication_id",
        "method",
        "target_variant",
    }:
        _fail("founder_config_invalid")
    if document.get("schema_version") != "1.0.0" or document.get(
        "config_state"
    ) != "FINALIZED":
        _fail("founder_config_invalid")
    configuration_id = _identifier(
        document.get("founder_simulation_config_id"), "configuration"
    )
    replication_id = _identifier(document.get("replication_id"), "replication")
    method = document.get("method")
    target = document.get("target_variant")
    if not isinstance(method, Mapping) or not isinstance(target, Mapping):
        _fail("founder_config_invalid")
    if set(method) != {
        "method_id",
        "method_version",
        "random_seed",
        "age_generations",
        "family_selection",
        "carrier_lineage_selection",
        "founder_haplotype_selection",
        "breakpoint_process",
        "edge_policy",
    }:
        _fail("founder_method_invalid")
    if (
        method.get("method_id") != FOUNDER_RECOMBINATION_METHOD_ID
        or method.get("method_version") != FOUNDER_RECOMBINATION_METHOD_VERSION
        or method.get("family_selection") != "ALL_D13_PEDIGREE_FAMILIES"
        or method.get("carrier_lineage_selection")
        != "RANDOM_TERMINAL_DESCENDANT_HAPLOTYPE"
        or method.get("founder_haplotype_selection")
        != "RANDOM_UNRELATED_D13_HAPLOTYPE"
        or method.get("breakpoint_process")
        != "INDEPENDENT_EXPONENTIAL_FIRST_EVENT_EACH_SIDE"
        or method.get("edge_policy") != "RIGHT_CENSOR_AT_MAP_DOMAIN"
    ):
        _fail("founder_method_invalid")
    seed = method.get("random_seed")
    age = method.get("age_generations")
    if type(seed) is not int or not 0 <= seed <= 4_294_967_295:
        _fail("founder_seed_invalid")
    if type(age) is not int or age < 1:
        _fail("founder_age_invalid")
    if set(target) != {
        "variant_id",
        "chromosome",
        "position_bp",
        "reference_allele",
        "alternate_allele",
        "position_policy",
    }:
        _fail("founder_target_invalid")
    variant_id = _identifier(target.get("variant_id"), "target_variant")
    chromosome = target.get("chromosome")
    position_bp = target.get("position_bp")
    reference = target.get("reference_allele")
    alternate = target.get("alternate_allele")
    if (
        type(chromosome) is not int
        or not 1 <= chromosome <= 22
        or type(position_bp) is not int
        or position_bp < 1
        or target.get("position_policy")
        != "EXISTING_BIALLELIC_SNV_RESET_TO_SYNTHETIC_MUTATION"
        or reference not in {"A", "C", "G", "T"}
        or alternate not in {"A", "C", "G", "T"}
        or reference == alternate
    ):
        _fail("founder_target_invalid")
    return FounderSimulationConfiguration(
        configuration_id,
        seed,
        replication_id,
        age,
        variant_id,
        chromosome,
        position_bp,
        reference,
        alternate,
    )


def parse_genetic_map(
    lines: Iterable[str], chromosome: int
) -> tuple[GeneticMapPoint, ...]:
    """Lit une carte D11 strictement croissante et sans coordonnées dupliquées."""
    reader = csv.DictReader(lines, delimiter="\t")
    if tuple(reader.fieldnames or ()) != (
        "CHROMOSOME",
        "POSITION_BP",
        "POSITION_CM",
    ):
        _fail("founder_genetic_map_columns_invalid")
    points: list[GeneticMapPoint] = []
    for row in reader:
        try:
            row_chromosome = int(row["CHROMOSOME"])
            position_bp = int(row["POSITION_BP"])
            position_cm = float(row["POSITION_CM"])
        except ValueError:
            _fail("founder_genetic_map_value_invalid")
        if (
            row_chromosome != chromosome
            or position_bp < 1
            or not math.isfinite(position_cm)
            or (points and position_bp < points[-1].position_bp)
            or (points and position_cm < points[-1].position_cm)
        ):
            _fail("founder_genetic_map_order_invalid")
        if points and position_bp == points[-1].position_bp:
            if position_cm != points[-1].position_cm:
                _fail("founder_genetic_map_duplicate_inconsistent")
            continue
        points.append(GeneticMapPoint(position_bp, position_cm))
    if len(points) < 2:
        _fail("founder_genetic_map_insufficient")
    return tuple(points)


def interpolate_cm(points: Sequence[GeneticMapPoint], position_bp: int) -> float:
    """Interpole linéairement une position physique sans extrapolation."""
    if position_bp < points[0].position_bp or position_bp > points[-1].position_bp:
        _fail("founder_position_outside_map")
    for left, right in zip(points, points[1:], strict=False):
        if left.position_bp <= position_bp <= right.position_bp:
            fraction = (position_bp - left.position_bp) / (
                right.position_bp - left.position_bp
            )
            return left.position_cm + fraction * (
                right.position_cm - left.position_cm
            )
    return points[-1].position_cm


def interpolate_bp(points: Sequence[GeneticMapPoint], position_cm: float) -> int:
    """Inverse la carte par interpolation et arrondit au bp le plus proche."""
    if position_cm < points[0].position_cm or position_cm > points[-1].position_cm:
        _fail("founder_cm_outside_map")
    for left, right in zip(points, points[1:], strict=False):
        if left.position_cm <= position_cm <= right.position_cm:
            fraction = (position_cm - left.position_cm) / (
                right.position_cm - left.position_cm
            )
            value = left.position_bp + fraction * (
                right.position_bp - left.position_bp
            )
            return int(math.floor(value + 0.5))
    return points[-1].position_bp


def parse_pedigree_truth(lines: Iterable[str]) -> tuple[PedigreeMember, ...]:
    """Lit la structure D13 privée requise par la sélection des lignées."""
    reader = csv.DictReader(lines, delimiter="\t")
    if tuple(reader.fieldnames or ()) != (
        "sample_id",
        "family_unit_id",
        "father_id",
        "mother_id",
        "sex",
        "generation",
        "is_founder",
    ):
        _fail("founder_pedigree_columns_invalid")
    records: list[PedigreeMember] = []
    seen: set[str] = set()
    for row in reader:
        sample_id = _identifier(row["sample_id"], "sample")
        family_id = _identifier(row["family_unit_id"], "family")
        try:
            generation = int(row["generation"])
        except ValueError:
            _fail("founder_pedigree_generation_invalid")
        if sample_id in seen or generation < 0 or row["is_founder"] not in {
            "true",
            "false",
        }:
            _fail("founder_pedigree_record_invalid")
        seen.add(sample_id)
        records.append(
            PedigreeMember(
                sample_id,
                family_id,
                generation,
                row["is_founder"] == "true",
            )
        )
    if not records:
        _fail("founder_pedigree_empty")
    return tuple(records)


def parse_d13_origins(
    lines: Iterable[str],
) -> dict[str, tuple[HaplotypeOrigin, HaplotypeOrigin]]:
    """Lit et valide les deux origines D12 de chaque individu D13."""
    reader = csv.DictReader(lines, delimiter="\t")
    if tuple(reader.fieldnames or ()) != (
        "sample_id",
        "d13_haplotype_index",
        "d12_sample_id",
        "d12_haplotype_index",
    ):
        _fail("founder_provenance_columns_invalid")
    values: defaultdict[str, dict[int, HaplotypeOrigin]] = defaultdict(dict)
    for row in reader:
        sample_id = _identifier(row["sample_id"], "sample")
        d12_sample = _identifier(row["d12_sample_id"], "d12_sample")
        try:
            output_index = int(row["d13_haplotype_index"])
            source_index = int(row["d12_haplotype_index"])
        except ValueError:
            _fail("founder_provenance_haplotype_invalid")
        if (
            output_index not in (1, 2)
            or source_index not in (1, 2)
            or output_index in values[sample_id]
        ):
            _fail("founder_provenance_haplotype_invalid")
        values[sample_id][output_index] = HaplotypeOrigin(d12_sample, source_index)
    if not values or any(set(copies) != {1, 2} for copies in values.values()):
        _fail("founder_provenance_incomplete")
    return {sample: (copies[1], copies[2]) for sample, copies in values.items()}


def parse_d12_haplotype_strata(
    lines: Iterable[str],
) -> dict[HaplotypeOrigin, str]:
    """Lie chaque origine D12 à sa strate haplotypique publique."""
    reader = csv.DictReader(lines, delimiter="\t")
    if tuple(reader.fieldnames or ()) != (
        "sample_id",
        "output_haplotype_index",
        "d11_sample_id",
        "d11_haplotype_index",
        "source_stratum_id",
    ):
        _fail("founder_d12_provenance_columns_invalid")
    strata: dict[HaplotypeOrigin, str] = {}
    for row in reader:
        try:
            index = int(row["output_haplotype_index"])
        except ValueError:
            _fail("founder_d12_provenance_haplotype_invalid")
        origin = HaplotypeOrigin(
            _identifier(row["sample_id"], "d12_sample"), index
        )
        stratum = _identifier(row["source_stratum_id"], "source_stratum")
        if index not in (1, 2) or origin in strata:
            _fail("founder_d12_provenance_haplotype_invalid")
        strata[origin] = stratum
    if not strata:
        _fail("founder_d12_provenance_empty")
    return strata


def build_founder_plan(
    *,
    configuration: FounderSimulationConfiguration,
    genetic_map: Sequence[GeneticMapPoint],
    pedigree: Sequence[PedigreeMember],
    sample_origins: Mapping[str, tuple[HaplotypeOrigin, HaplotypeOrigin]],
    background_strata: Mapping[HaplotypeOrigin, str],
) -> FounderSimulationPlan:
    """Sélectionne les lignées et simule les premiers breakpoints D14."""
    if {member.sample_id for member in pedigree} != set(sample_origins):
        _fail("founder_sample_universe_mismatch")
    all_origins = {origin for copies in sample_origins.values() for origin in copies}
    if not all_origins <= set(background_strata):
        _fail("founder_stratum_origin_missing")
    target_cm = interpolate_cm(genetic_map, configuration.target_position_bp)
    family_members: defaultdict[str, list[PedigreeMember]] = defaultdict(list)
    unrelated_members: list[PedigreeMember] = []
    for member in pedigree:
        if member.family_unit_id.startswith("UNRELATED_"):
            unrelated_members.append(member)
        else:
            family_members[member.family_unit_id].append(member)
    if not family_members or not unrelated_members:
        _fail("founder_family_universe_invalid")
    origin_families: dict[HaplotypeOrigin, str] = {}
    for family_id, members in family_members.items():
        for member in members:
            for origin in sample_origins[member.sample_id]:
                previous_family = origin_families.setdefault(origin, family_id)
                if previous_family != family_id:
                    _fail("founder_family_origin_overlap")
    generator = random.Random(configuration.random_seed)
    template_candidates = sorted(
        origin
        for member in unrelated_members
        for origin in sample_origins[member.sample_id]
    )
    founder_template = generator.choice(template_candidates)
    if founder_template in origin_families:
        _fail("founder_template_family_overlap")
    segments: list[FounderFamilySegment] = []
    carriers: list[CarrierHaplotype] = []
    used_origins: set[HaplotypeOrigin] = set()
    for family_id in sorted(family_members):
        members = family_members[family_id]
        terminal_generation = max(member.generation for member in members)
        terminal_members = sorted(
            (
                member
                for member in members
                if member.generation == terminal_generation
            ),
            key=lambda member: member.sample_id,
        )
        selected_member = generator.choice(terminal_members)
        selected_haplotype = generator.choice((1, 2))
        carrier_origin = sample_origins[selected_member.sample_id][
            selected_haplotype - 1
        ]
        if carrier_origin in used_origins or carrier_origin == founder_template:
            _fail("founder_carrier_origin_collision")
        used_origins.add(carrier_origin)
        left_draw_cm = generator.expovariate(
            configuration.age_generations / 100.0
        )
        right_draw_cm = generator.expovariate(
            configuration.age_generations / 100.0
        )
        raw_left_cm = target_cm - left_draw_cm
        raw_right_cm = target_cm + right_draw_cm
        left_censored = raw_left_cm < genetic_map[0].position_cm
        right_censored = raw_right_cm > genetic_map[-1].position_cm
        left_bp = (
            genetic_map[0].position_bp
            if left_censored
            else interpolate_bp(genetic_map, raw_left_cm)
        )
        right_bp = (
            genetic_map[-1].position_bp
            if right_censored
            else interpolate_bp(genetic_map, raw_right_cm)
        )
        left_cm = interpolate_cm(genetic_map, left_bp)
        right_cm = interpolate_cm(genetic_map, right_bp)
        segment = FounderFamilySegment(
            family_id,
            carrier_origin,
            left_bp,
            right_bp,
            left_cm,
            target_cm,
            right_cm,
            left_censored,
            right_censored,
            int(not left_censored) + int(not right_censored),
        )
        segments.append(segment)
        for member in members:
            for index, origin in enumerate(sample_origins[member.sample_id], start=1):
                if origin == carrier_origin:
                    carriers.append(
                        CarrierHaplotype(
                            family_id,
                            member.sample_id,
                            index,
                            carrier_origin,
                        )
                    )
    if any(not any(row.family_unit_id == segment.family_unit_id for row in carriers)
           for segment in segments):
        _fail("founder_carrier_missing")
    return FounderSimulationPlan(
        configuration,
        f"{configuration.replication_id}.founder_haplotype",
        founder_template,
        background_strata[founder_template],
        tuple(segments),
        tuple(sorted(carriers)),
        dict(sample_origins),
        dict(background_strata),
    )


def canonical_founder_segments_bytes(plan: FounderSimulationPlan) -> bytes:
    """Encode les limites vraies, longueurs et censures par famille."""
    output = io.StringIO(newline="")
    writer = csv.writer(output, delimiter="\t", lineterminator="\n")
    writer.writerow(
        (
            "replication_id",
            "family_id",
            "age_generations",
            "left_boundary_bp",
            "right_boundary_bp",
            "left_boundary_cm",
            "target_position_cm",
            "right_boundary_cm",
            "left_length_cm",
            "right_length_cm",
            "total_length_cm",
            "left_censored",
            "right_censored",
            "recombination_event_count",
            "founder_haplotype_id",
        )
    )
    for segment in plan.segments:
        writer.writerow(
            (
                plan.configuration.replication_id,
                segment.family_unit_id,
                plan.configuration.age_generations,
                segment.left_boundary_bp,
                segment.right_boundary_bp,
                f"{segment.left_boundary_cm:.12g}",
                f"{segment.target_position_cm:.12g}",
                f"{segment.right_boundary_cm:.12g}",
                f"{segment.left_length_cm:.12g}",
                f"{segment.right_length_cm:.12g}",
                f"{segment.total_length_cm:.12g}",
                str(segment.left_censored).lower(),
                str(segment.right_censored).lower(),
                segment.recombination_event_count,
                plan.founder_haplotype_id,
            )
        )
    return output.getvalue().encode("utf-8")


def canonical_carrier_bytes(plan: FounderSimulationPlan) -> bytes:
    """Encode les copies présentes qui portent le segment et la mutation."""
    output = io.StringIO(newline="")
    writer = csv.writer(output, delimiter="\t", lineterminator="\n")
    writer.writerow(
        (
            "family_id",
            "sample_id",
            "carrier_haplotype_index",
            "founder_haplotype_id",
        )
    )
    for carrier in plan.carriers:
        writer.writerow(
            (
                carrier.family_unit_id,
                carrier.sample_id,
                carrier.haplotype_index,
                plan.founder_haplotype_id,
            )
        )
    return output.getvalue().encode("utf-8")


def canonical_target_truth_bytes(plan: FounderSimulationPlan) -> bytes:
    """Encode le génotype cible vrai de tous les individus D14."""
    carrier_indexes: defaultdict[str, list[int]] = defaultdict(list)
    for carrier in plan.carriers:
        carrier_indexes[carrier.sample_id].append(carrier.haplotype_index)
    output = io.StringIO(newline="")
    writer = csv.writer(output, delimiter="\t", lineterminator="\n")
    writer.writerow(
        ("sample_id", "variant_id", "genotype", "mutant_haplotype")
    )
    for sample_id in sorted(plan.sample_origins):
        indexes = carrier_indexes[sample_id]
        genotype = "0/1" if indexes else "0/0"
        mutant = f"H{indexes[0]}" if indexes else "NONE"
        writer.writerow(
            (sample_id, plan.configuration.target_variant_id, genotype, mutant)
        )
    return output.getvalue().encode("utf-8")


def canonical_local_ancestry_bytes(plan: FounderSimulationPlan) -> bytes:
    """Encode le remplacement local d'ascendance de chaque origine porteuse."""
    output = io.StringIO(newline="")
    writer = csv.writer(output, delimiter="\t", lineterminator="\n")
    writer.writerow(
        (
            "family_id",
            "background_stratum_id",
            "founder_stratum_id",
            "left_boundary_bp",
            "right_boundary_bp",
        )
    )
    for segment in plan.segments:
        writer.writerow(
            (
                segment.family_unit_id,
                plan.background_strata[segment.carrier_origin],
                plan.founder_source_stratum_id,
                segment.left_boundary_bp,
                segment.right_boundary_bp,
            )
        )
    return output.getvalue().encode("utf-8")


def _parse_gt(value: str, gt_index: int) -> tuple[int, int]:
    fields = value.split(":")
    if gt_index >= len(fields):
        _fail("founder_vcf_sample_incomplete")
    alleles = fields[gt_index].split("|")
    if len(alleles) != 2 or any(allele not in {"0", "1"} for allele in alleles):
        _fail("founder_vcf_genotype_invalid")
    return int(alleles[0]), int(alleles[1])


def transform_founder_vcf(
    source_lines: Iterable[str],
    dense_output: TextIO,
    observed_output: TextIO,
    *,
    plan: FounderSimulationPlan,
    observed_variant_keys: Mapping[tuple[str, int, str, str], str],
    transmissions: Sequence[tuple[str, str, int, int]],
) -> FounderVcfStats:
    """Mosaïque le fondateur, injecte la mutation et contrôle Mendel."""
    sample_ids = sorted(plan.sample_origins)
    carrier_lookup = {
        (carrier.sample_id, carrier.haplotype_index): carrier.family_unit_id
        for carrier in plan.carriers
    }
    segments = {segment.family_unit_id: segment for segment in plan.segments}
    representative: dict[HaplotypeOrigin, tuple[str, int]] = {}
    for sample_id, origins in plan.sample_origins.items():
        for index, origin in enumerate(origins, start=1):
            representative.setdefault(origin, (sample_id, index))
    template_sample, template_index = representative[plan.founder_template_origin]
    header = (
        "##fileformat=VCFv4.2\n"
        f"##source={FOUNDER_RECOMBINATION_METHOD_ID}_"
        f"{FOUNDER_RECOMBINATION_METHOD_VERSION}\n"
        f"##contig=<ID=chr{plan.configuration.chromosome}>\n"
        '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n'
    )
    columns = "\t".join(
        (
            "#CHROM",
            "POS",
            "ID",
            "REF",
            "ALT",
            "QUAL",
            "FILTER",
            "INFO",
            "FORMAT",
            *sample_ids,
        )
    )
    dense_output.write(header + columns + "\n")
    observed_output.write(header + columns + "\n")
    source_index: dict[str, int] | None = None
    dense_count = 0
    observed_count = 0
    target_count = 0
    changed = 0
    mendelian_errors = 0
    observed_seen: set[tuple[str, int, str, str]] = set()
    for raw_line in source_lines:
        if not raw_line.strip() or raw_line.startswith("##"):
            continue
        if raw_line.startswith("#CHROM\t"):
            source_samples = raw_line.rstrip("\r\n").split("\t")[9:]
            if set(source_samples) != set(sample_ids) or len(source_samples) != len(
                set(source_samples)
            ):
                _fail("founder_vcf_sample_universe_mismatch")
            source_index = {
                sample_id: index for index, sample_id in enumerate(source_samples)
            }
            continue
        if source_index is None:
            _fail("founder_vcf_header_missing")
        fields = raw_line.rstrip("\r\n").split("\t")
        if len(fields) != 9 + len(source_index):
            _fail("founder_vcf_column_count_invalid")
        try:
            position_bp = int(fields[1])
        except ValueError:
            _fail("founder_vcf_position_invalid")
        if fields[0] not in {
            str(plan.configuration.chromosome),
            f"chr{plan.configuration.chromosome}",
        }:
            _fail("founder_vcf_chromosome_mismatch")
        format_fields = fields[8].split(":")
        if "GT" not in format_fields:
            _fail("founder_vcf_gt_missing")
        gt_index = format_fields.index("GT")
        original = {
            sample_id: _parse_gt(fields[9 + index], gt_index)
            for sample_id, index in source_index.items()
        }
        if position_bp == plan.configuration.target_position_bp:
            if (
                fields[3] != plan.configuration.reference_allele
                or fields[4] != plan.configuration.alternate_allele
            ):
                _fail("founder_target_alleles_mismatch")
            target_count += 1
        template_allele = original[template_sample][template_index - 1]
        output_genotypes: dict[str, tuple[int, int]] = {}
        phased: list[str] = []
        unphased: list[str] = []
        for sample_id in sample_ids:
            alleles = list(original[sample_id])
            for haplotype_index in (1, 2):
                family_id = carrier_lookup.get((sample_id, haplotype_index))
                if family_id is not None:
                    segment = segments[family_id]
                    if (
                        segment.left_boundary_bp
                        <= position_bp
                        <= segment.right_boundary_bp
                    ):
                        if alleles[haplotype_index - 1] != template_allele:
                            changed += 1
                        alleles[haplotype_index - 1] = template_allele
            if position_bp == plan.configuration.target_position_bp:
                alleles = [
                    int((sample_id, index) in carrier_lookup) for index in (1, 2)
                ]
            genotype = (alleles[0], alleles[1])
            output_genotypes[sample_id] = genotype
            phased.append(f"{genotype[0]}|{genotype[1]}")
            ordered = sorted(genotype)
            unphased.append(f"{ordered[0]}/{ordered[1]}")
        for child_id, parent_id, parental_index, child_index in transmissions:
            if output_genotypes[child_id][child_index - 1] != output_genotypes[
                parent_id
            ][parental_index - 1]:
                mendelian_errors += 1
        variant_id = (
            plan.configuration.target_variant_id
            if position_bp == plan.configuration.target_position_bp
            else fields[2]
        )
        fixed = (
            f"chr{plan.configuration.chromosome}",
            fields[1],
            variant_id,
            fields[3],
            fields[4],
            ".",
            "PASS",
            ".",
            "GT",
        )
        dense_output.write("\t".join((*fixed, *phased)) + "\n")
        dense_count += 1
        key = (fixed[0], position_bp, fields[3], fields[4])
        if key in observed_variant_keys:
            observed_fixed = (
                fixed[0],
                fixed[1],
                observed_variant_keys[key],
                fixed[3],
                fixed[4],
                ".",
                "PASS",
                ".",
                "GT",
            )
            observed_output.write("\t".join((*observed_fixed, *unphased)) + "\n")
            observed_seen.add(key)
            observed_count += 1
    if source_index is None:
        _fail("founder_vcf_header_missing")
    if target_count != 1:
        _fail("founder_target_variant_count_invalid")
    if observed_seen != set(observed_variant_keys):
        _fail("founder_observed_variant_mismatch")
    return FounderVcfStats(
        dense_count,
        observed_count,
        len({carrier.sample_id for carrier in plan.carriers}),
        len(plan.carriers),
        changed,
        mendelian_errors,
        target_count,
    )
