"""Génération D13 de pedigrees et transmissions mendéliennes exactes.

Chaque famille de cousins germains est construite à partir de quatre fondateurs
D12 distincts d'une même strate. Les quatre descendants reçoivent une copie
haplotypique entière de chaque parent. D13 ne simule aucune recombinaison : ce
choix rend les transmissions vérifiables exactement et réserve les mosaïques
chromosomiques au modèle D14.
"""

from __future__ import annotations

from collections import Counter, defaultdict
from collections.abc import Iterable, Mapping, Sequence
import csv
from dataclasses import dataclass, replace
from fractions import Fraction
import io
import json
import random
import re
from typing import TextIO


PEDIGREE_METHOD_ID = "first_cousin_whole_haplotype_transmission_v1"
PEDIGREE_METHOD_VERSION = "1.0.0"
_STABLE_IDENTIFIER = re.compile(r"^[A-Za-z0-9][A-Za-z0-9._-]*$")
_RELATIONSHIP_BY_KINSHIP = {
    Fraction(0, 1): "UNRELATED",
    Fraction(1, 16): "THIRD_DEGREE",
    Fraction(1, 8): "SECOND_DEGREE",
    Fraction(1, 4): "FIRST_DEGREE",
}


class PedigreeSimulationError(RuntimeError):
    """Signale une violation D13 avec un code stable et non sensible."""

    def __init__(self, code: str, detail: str | None = None) -> None:
        self.code = code
        self.detail = detail
        super().__init__(code if detail is None else f"{code}: {detail}")


@dataclass(frozen=True, order=True)
class SourceIndividual:
    """Individu D12 et vérité d'ascendance nécessaire à la stratification."""

    sample_id: str
    stratum_id: str
    mixture_proportions: tuple[tuple[str, float], ...]
    haplotype_1_stratum_id: str
    haplotype_2_stratum_id: str


@dataclass(frozen=True, order=True)
class PedigreeFamilySpecification:
    """Nombre de familles indépendantes à construire dans une strate D12."""

    family_group_id: str
    stratum_id: str
    family_count: int


@dataclass(frozen=True, order=True)
class PedigreeHaplotypeOrigin:
    """Copie D12 ultime transmise à une copie haplotypique D13."""

    d12_sample_id: str
    d12_haplotype_index: int


@dataclass(frozen=True)
class PedigreeIndividual:
    """Individu D13 ordonné topologiquement dans son pedigree."""

    sample_id: str
    family_unit_id: str
    father_id: str | None
    mother_id: str | None
    sex: str
    generation: int
    is_founder: bool
    stratum_id: str
    mixture_proportions: tuple[tuple[str, float], ...]
    haplotype_1: PedigreeHaplotypeOrigin
    haplotype_2: PedigreeHaplotypeOrigin


@dataclass(frozen=True, order=True)
class TransmissionRecord:
    """Copie parentale transmise à une copie enfant, sans recombinaison."""

    family_unit_id: str
    child_id: str
    parent_id: str
    parental_haplotype_index: int
    child_haplotype_index: int


@dataclass(frozen=True, order=True)
class KinshipRecord:
    """Vérité de parenté théorique dérivée du pedigree."""

    sample_1_id: str
    sample_2_id: str
    relationship_class: str
    kinship_coefficient: float
    independence_unit_id: str


@dataclass(frozen=True)
class PedigreePanel:
    """Panel D13 complet et ses vérités privées."""

    individuals: tuple[PedigreeIndividual, ...]
    transmissions: tuple[TransmissionRecord, ...]
    kinship_records: tuple[KinshipRecord, ...]
    source_individual_count: int
    family_count: int
    trio_count: int
    unrelated_pair_count: int


@dataclass(frozen=True)
class PedigreeVcfStats:
    """Contrôles agrégés de la transformation VCF D13."""

    dense_variant_count: int
    observed_variant_count: int
    source_individual_count: int
    output_individual_count: int
    descendant_count: int
    founder_haplotype_copy_error_count: int
    founder_allele_count_difference: int
    mendelian_error_count: int


def _fail(code: str, detail: str | None = None) -> None:
    raise PedigreeSimulationError(code, detail)


def _identifier(value: object, namespace: str) -> str:
    if not isinstance(value, str) or _STABLE_IDENTIFIER.fullmatch(value) is None:
        _fail("pedigree_identifier_invalid", namespace)
    return value


def parse_d12_source_truth(
    ancestry_lines: Iterable[str],
    provenance_lines: Iterable[str],
) -> tuple[SourceIndividual, ...]:
    """Lie l'ascendance diploïde D12 à la strate de chaque copie haplotypique."""
    reader = csv.DictReader(ancestry_lines, delimiter="\t")
    if tuple(reader.fieldnames or ()) != (
        "sample_id",
        "stratum_id",
        "mixture_proportions",
    ):
        _fail("pedigree_ancestry_columns_invalid")
    ancestry: dict[str, tuple[str, tuple[tuple[str, float], ...]]] = {}
    identifiers: set[str] = set()
    for row in reader:
        sample_id = _identifier(row["sample_id"], "d12_sample")
        stratum_id = _identifier(row["stratum_id"], "d12_stratum")
        if sample_id.casefold() in identifiers:
            _fail("pedigree_ancestry_sample_duplicate")
        identifiers.add(sample_id.casefold())
        try:
            decoded = json.loads(row["mixture_proportions"])
        except json.JSONDecodeError:
            _fail("pedigree_ancestry_proportions_invalid")
        if not isinstance(decoded, dict) or not decoded:
            _fail("pedigree_ancestry_proportions_invalid")
        proportions: list[tuple[str, float]] = []
        for component, raw_value in decoded.items():
            component_id = _identifier(component, "ancestry_component")
            if type(raw_value) not in (int, float):
                _fail("pedigree_ancestry_proportions_invalid")
            value = float(raw_value)
            if value <= 0 or value > 1:
                _fail("pedigree_ancestry_proportions_invalid")
            proportions.append((component_id, value))
        proportions.sort()
        if abs(sum(value for _, value in proportions) - 1.0) > 1e-12:
            _fail("pedigree_ancestry_proportion_sum_invalid")
        ancestry[sample_id] = (stratum_id, tuple(proportions))
    if not ancestry:
        _fail("pedigree_ancestry_empty")

    provenance_reader = csv.DictReader(provenance_lines, delimiter="\t")
    if tuple(provenance_reader.fieldnames or ()) != (
        "sample_id",
        "output_haplotype_index",
        "d11_sample_id",
        "d11_haplotype_index",
        "source_stratum_id",
    ):
        _fail("pedigree_d12_provenance_columns_invalid")
    haplotype_strata: defaultdict[str, dict[int, str]] = defaultdict(dict)
    for row in provenance_reader:
        sample_id = _identifier(row["sample_id"], "d12_provenance_sample")
        _identifier(row["d11_sample_id"], "d11_provenance_sample")
        stratum_id = _identifier(
            row["source_stratum_id"], "d12_haplotype_stratum"
        )
        try:
            haplotype_index = int(row["output_haplotype_index"])
            d11_haplotype_index = int(row["d11_haplotype_index"])
        except ValueError:
            _fail("pedigree_d12_provenance_haplotype_invalid")
        if (
            haplotype_index not in (1, 2)
            or d11_haplotype_index not in (1, 2)
            or haplotype_index in haplotype_strata[sample_id]
        ):
            _fail("pedigree_d12_provenance_haplotype_invalid")
        haplotype_strata[sample_id][haplotype_index] = stratum_id
    if set(haplotype_strata) != set(ancestry):
        _fail("pedigree_d12_provenance_sample_mismatch")

    records: list[SourceIndividual] = []
    for sample_id in sorted(ancestry):
        strata = haplotype_strata[sample_id]
        if set(strata) != {1, 2}:
            _fail("pedigree_d12_provenance_haplotype_missing")
        realized = Counter((strata[1], strata[2]))
        realized_proportions = tuple(
            sorted((key, value / 2) for key, value in realized.items())
        )
        stratum_id, declared_proportions = ancestry[sample_id]
        if realized_proportions != declared_proportions:
            _fail("pedigree_d12_ancestry_provenance_mismatch")
        records.append(
            SourceIndividual(
                sample_id,
                stratum_id,
                declared_proportions,
                strata[1],
                strata[2],
            )
        )
    return tuple(records)


def parse_pedigree_configuration(
    config: Mapping[str, object],
) -> tuple[int, tuple[PedigreeFamilySpecification, ...], int]:
    """Valide la configuration D13 finalisée et retourne ses paramètres."""
    if set(config) != {
        "schema_version",
        "pedigree_config_id",
        "config_state",
        "method",
        "family_groups",
        "unrelated_pair_count",
    }:
        _fail("pedigree_config_invalid")
    if config.get("schema_version") != "1.0.0" or config.get(
        "config_state"
    ) != "FINALIZED":
        _fail("pedigree_config_invalid")
    _identifier(config.get("pedigree_config_id"), "configuration")
    method = config.get("method")
    groups = config.get("family_groups")
    unrelated_pair_count = config.get("unrelated_pair_count")
    if not isinstance(method, Mapping) or not isinstance(groups, list):
        _fail("pedigree_config_invalid")
    if set(method) != {
        "method_id",
        "method_version",
        "random_seed",
        "family_template",
        "source_selection",
        "recombination_policy",
    }:
        _fail("pedigree_method_invalid")
    if (
        method.get("method_id") != PEDIGREE_METHOD_ID
        or method.get("method_version") != PEDIGREE_METHOD_VERSION
        or method.get("family_template") != "FIRST_COUSIN_PEDIGREE"
        or method.get("source_selection")
        != "SAME_STRATUM_WITHOUT_REPLACEMENT"
        or method.get("recombination_policy")
        != "FORBIDDEN_DEFERRED_TO_D14"
    ):
        _fail("pedigree_method_invalid")
    random_seed = method.get("random_seed")
    if (
        type(random_seed) is not int
        or random_seed < 0
        or random_seed > 4_294_967_295
    ):
        _fail("pedigree_seed_invalid")
    if type(unrelated_pair_count) is not int or unrelated_pair_count < 1:
        _fail("pedigree_unrelated_pair_count_invalid")
    specifications: list[PedigreeFamilySpecification] = []
    group_ids: set[str] = set()
    strata: set[str] = set()
    for group in groups:
        if not isinstance(group, Mapping) or set(group) != {
            "family_group_id",
            "stratum_id",
            "family_count",
        }:
            _fail("pedigree_family_group_invalid")
        group_id = _identifier(group.get("family_group_id"), "family_group")
        stratum_id = _identifier(group.get("stratum_id"), "family_stratum")
        family_count = group.get("family_count")
        if type(family_count) is not int or family_count < 1:
            _fail("pedigree_family_count_invalid")
        if group_id.casefold() in group_ids or stratum_id in strata:
            _fail("pedigree_family_group_duplicate")
        group_ids.add(group_id.casefold())
        strata.add(stratum_id)
        specifications.append(
            PedigreeFamilySpecification(group_id, stratum_id, family_count)
        )
    if not specifications:
        _fail("pedigree_family_group_empty")
    return random_seed, tuple(specifications), unrelated_pair_count


def _realized_ancestry(
    origins: Sequence[PedigreeHaplotypeOrigin],
    haplotype_strata: Mapping[tuple[str, int], str],
) -> tuple[str, tuple[tuple[str, float], ...]]:
    strata = [
        haplotype_strata[(origin.d12_sample_id, origin.d12_haplotype_index)]
        for origin in origins
    ]
    counts = Counter(strata)
    proportions = tuple(sorted((key, value / 2) for key, value in counts.items()))
    stratum_id = strata[0] if len(counts) == 1 else "ADMIX_" + "_".join(sorted(counts))
    return stratum_id, proportions


def _transmitted_child(
    *,
    sample_id: str,
    family_unit_id: str,
    father: PedigreeIndividual,
    mother: PedigreeIndividual,
    sex: str,
    generation: int,
    generator: random.Random,
    haplotype_strata: Mapping[tuple[str, int], str],
) -> tuple[PedigreeIndividual, tuple[TransmissionRecord, TransmissionRecord]]:
    paternal_index = generator.choice((1, 2))
    maternal_index = generator.choice((1, 2))
    paternal_origin = (father.haplotype_1, father.haplotype_2)[paternal_index - 1]
    maternal_origin = (mother.haplotype_1, mother.haplotype_2)[maternal_index - 1]
    stratum_id, proportions = _realized_ancestry(
        (paternal_origin, maternal_origin), haplotype_strata
    )
    child = PedigreeIndividual(
        sample_id=sample_id,
        family_unit_id=family_unit_id,
        father_id=father.sample_id,
        mother_id=mother.sample_id,
        sex=sex,
        generation=generation,
        is_founder=False,
        stratum_id=stratum_id,
        mixture_proportions=proportions,
        haplotype_1=paternal_origin,
        haplotype_2=maternal_origin,
    )
    transmissions = (
        TransmissionRecord(
            family_unit_id,
            sample_id,
            father.sample_id,
            paternal_index,
            1,
        ),
        TransmissionRecord(
            family_unit_id,
            sample_id,
            mother.sample_id,
            maternal_index,
            2,
        ),
    )
    return child, transmissions


def _kinship_matrix(
    individuals: Sequence[PedigreeIndividual],
) -> dict[tuple[str, str], Fraction]:
    matrix: dict[tuple[str, str], Fraction] = {}
    seen: set[str] = set()
    for individual in individuals:
        if individual.sample_id in seen:
            _fail("pedigree_sample_duplicate")
        if individual.is_founder:
            if individual.father_id is not None or individual.mother_id is not None:
                _fail("pedigree_founder_parent_forbidden")
            for previous in seen:
                matrix[(individual.sample_id, previous)] = Fraction(0)
                matrix[(previous, individual.sample_id)] = Fraction(0)
            diagonal = Fraction(1, 2)
        else:
            if (
                individual.father_id not in seen
                or individual.mother_id not in seen
                or individual.father_id == individual.mother_id
            ):
                _fail("pedigree_parent_order_invalid")
            for previous in seen:
                coefficient = (
                    matrix[(individual.father_id, previous)]
                    + matrix[(individual.mother_id, previous)]
                ) / 2
                matrix[(individual.sample_id, previous)] = coefficient
                matrix[(previous, individual.sample_id)] = coefficient
            diagonal = (
                Fraction(1)
                + matrix[(individual.father_id, individual.mother_id)]
            ) / 2
        matrix[(individual.sample_id, individual.sample_id)] = diagonal
        seen.add(individual.sample_id)
    return matrix


def build_pedigree_panel(
    source_individuals: Sequence[SourceIndividual],
    family_specifications: Sequence[PedigreeFamilySpecification],
    *,
    unrelated_pair_count: int,
    random_seed: int,
    synthetic_identifier_prefix: str,
) -> PedigreePanel:
    """Construit des familles indépendantes couvrant les degrés un à trois."""
    if not source_individuals:
        _fail("pedigree_source_empty")
    if (
        type(random_seed) is not int
        or random_seed < 0
        or random_seed > 4_294_967_295
    ):
        _fail("pedigree_seed_invalid")
    if _STABLE_IDENTIFIER.fullmatch(f"{synthetic_identifier_prefix}1") is None:
        _fail("pedigree_output_prefix_invalid")
    if type(unrelated_pair_count) is not int or unrelated_pair_count < 1:
        _fail("pedigree_unrelated_pair_count_invalid")
    source_ids = [source.sample_id for source in source_individuals]
    if len(source_ids) != len({value.casefold() for value in source_ids}):
        _fail("pedigree_source_duplicate")
    grouped: defaultdict[str, list[SourceIndividual]] = defaultdict(list)
    haplotype_strata: dict[tuple[str, int], str] = {}
    for source in source_individuals:
        _identifier(source.sample_id, "d12_sample")
        _identifier(source.stratum_id, "d12_stratum")
        _identifier(source.haplotype_1_stratum_id, "d12_haplotype_stratum")
        _identifier(source.haplotype_2_stratum_id, "d12_haplotype_stratum")
        grouped[source.stratum_id].append(source)
        haplotype_strata[(source.sample_id, 1)] = source.haplotype_1_stratum_id
        haplotype_strata[(source.sample_id, 2)] = source.haplotype_2_stratum_id

    generator = random.Random(random_seed)
    for values in grouped.values():
        values.sort()
        generator.shuffle(values)
    family_requests: list[tuple[str, str]] = []
    used_group_ids: set[str] = set()
    used_strata: set[str] = set()
    for specification in family_specifications:
        _identifier(specification.family_group_id, "family_group")
        _identifier(specification.stratum_id, "family_stratum")
        if (
            type(specification.family_count) is not int
            or specification.family_count < 1
            or specification.family_group_id.casefold() in used_group_ids
            or specification.stratum_id in used_strata
        ):
            _fail("pedigree_family_specification_invalid")
        used_group_ids.add(specification.family_group_id.casefold())
        used_strata.add(specification.stratum_id)
        if 4 * specification.family_count > len(grouped[specification.stratum_id]):
            _fail("pedigree_family_source_capacity_insufficient")
        family_requests.extend(
            (
                f"{specification.family_group_id}_{index:04d}",
                specification.stratum_id,
            )
            for index in range(1, specification.family_count + 1)
        )

    width = max(6, len(str(len(source_individuals) + 4 * len(family_requests))))
    source_to_output = {
        source.sample_id: f"{synthetic_identifier_prefix}{index:0{width}d}"
        for index, source in enumerate(sorted(source_individuals), start=1)
    }
    individuals_by_id: dict[str, PedigreeIndividual] = {}
    for source in source_individuals:
        output_id = source_to_output[source.sample_id]
        individuals_by_id[output_id] = PedigreeIndividual(
            sample_id=output_id,
            family_unit_id=f"UNRELATED_{output_id}",
            father_id=None,
            mother_id=None,
            sex="UNKNOWN",
            generation=0,
            is_founder=True,
            stratum_id=source.stratum_id,
            mixture_proportions=source.mixture_proportions,
            haplotype_1=PedigreeHaplotypeOrigin(source.sample_id, 1),
            haplotype_2=PedigreeHaplotypeOrigin(source.sample_id, 2),
        )

    cursors: defaultdict[str, int] = defaultdict(int)
    descendants: list[PedigreeIndividual] = []
    transmissions: list[TransmissionRecord] = []
    family_members: dict[str, tuple[PedigreeIndividual, ...]] = {}
    next_index = len(source_individuals) + 1
    selected_founder_ids: set[str] = set()
    for family_unit_id, stratum_id in family_requests:
        start = cursors[stratum_id]
        selected_sources = grouped[stratum_id][start : start + 4]
        cursors[stratum_id] += 4
        selected = [
            individuals_by_id[source_to_output[value.sample_id]]
            for value in selected_sources
        ]
        proportions = {value.mixture_proportions for value in selected}
        if len(proportions) != 1:
            _fail("pedigree_family_ancestry_mismatch")
        role_sexes = ("MALE", "FEMALE", "FEMALE", "MALE")
        updated_founders: list[PedigreeIndividual] = []
        for founder, sex in zip(selected, role_sexes, strict=True):
            updated = replace(
                founder,
                family_unit_id=family_unit_id,
                sex=sex,
            )
            individuals_by_id[updated.sample_id] = updated
            updated_founders.append(updated)
            selected_founder_ids.add(updated.sample_id)
        grandfather, grandmother, spouse_1, spouse_2 = updated_founders

        parent_1_id = f"{synthetic_identifier_prefix}{next_index:0{width}d}"
        next_index += 1
        parent_1, records = _transmitted_child(
            sample_id=parent_1_id,
            family_unit_id=family_unit_id,
            father=grandfather,
            mother=grandmother,
            sex="MALE",
            generation=1,
            generator=generator,
            haplotype_strata=haplotype_strata,
        )
        descendants.append(parent_1)
        transmissions.extend(records)

        parent_2_id = f"{synthetic_identifier_prefix}{next_index:0{width}d}"
        next_index += 1
        parent_2, records = _transmitted_child(
            sample_id=parent_2_id,
            family_unit_id=family_unit_id,
            father=grandfather,
            mother=grandmother,
            sex="FEMALE",
            generation=1,
            generator=generator,
            haplotype_strata=haplotype_strata,
        )
        descendants.append(parent_2)
        transmissions.extend(records)

        cousin_1_id = f"{synthetic_identifier_prefix}{next_index:0{width}d}"
        next_index += 1
        cousin_1, records = _transmitted_child(
            sample_id=cousin_1_id,
            family_unit_id=family_unit_id,
            father=parent_1,
            mother=spouse_1,
            sex="FEMALE",
            generation=2,
            generator=generator,
            haplotype_strata=haplotype_strata,
        )
        descendants.append(cousin_1)
        transmissions.extend(records)

        cousin_2_id = f"{synthetic_identifier_prefix}{next_index:0{width}d}"
        next_index += 1
        cousin_2, records = _transmitted_child(
            sample_id=cousin_2_id,
            family_unit_id=family_unit_id,
            father=spouse_2,
            mother=parent_2,
            sex="MALE",
            generation=2,
            generator=generator,
            haplotype_strata=haplotype_strata,
        )
        descendants.append(cousin_2)
        transmissions.extend(records)
        family_members[family_unit_id] = (
            grandfather,
            grandmother,
            spouse_1,
            spouse_2,
            parent_1,
            parent_2,
            cousin_1,
            cousin_2,
        )

    founders = tuple(
        sorted(individuals_by_id.values(), key=lambda value: value.sample_id)
    )
    individuals = founders + tuple(descendants)
    kinship_records: list[KinshipRecord] = []
    for family_unit_id in sorted(family_members):
        members = family_members[family_unit_id]
        matrix = _kinship_matrix(members)
        member_ids = [member.sample_id for member in members]
        for left_index, sample_1 in enumerate(member_ids):
            for sample_2 in member_ids[left_index + 1 :]:
                coefficient = matrix[(sample_1, sample_2)]
                relationship_class = _RELATIONSHIP_BY_KINSHIP.get(coefficient)
                if relationship_class is None:
                    _fail("pedigree_kinship_class_unsupported")
                kinship_records.append(
                    KinshipRecord(
                        min(sample_1, sample_2),
                        max(sample_1, sample_2),
                        relationship_class,
                        float(coefficient),
                        family_unit_id,
                    )
                )

    unrelated_candidates = [
        individual
        for individual in founders
        if individual.sample_id not in selected_founder_ids
    ]
    generator.shuffle(unrelated_candidates)
    if 2 * unrelated_pair_count > len(unrelated_candidates):
        _fail("pedigree_unrelated_pair_capacity_insufficient")
    for index in range(unrelated_pair_count):
        left = unrelated_candidates[2 * index]
        right = unrelated_candidates[2 * index + 1]
        kinship_records.append(
            KinshipRecord(
                min(left.sample_id, right.sample_id),
                max(left.sample_id, right.sample_id),
                "UNRELATED",
                0.0,
                f"UNRELATED_PAIR_{index + 1:06d}",
            )
        )

    panel = PedigreePanel(
        individuals=individuals,
        transmissions=tuple(sorted(transmissions)),
        kinship_records=tuple(sorted(set(kinship_records))),
        source_individual_count=len(source_individuals),
        family_count=len(family_requests),
        trio_count=4 * len(family_requests),
        unrelated_pair_count=unrelated_pair_count,
    )
    validate_pedigree_panel(panel, source_individuals)
    return panel


def validate_pedigree_panel(
    panel: PedigreePanel,
    source_individuals: Sequence[SourceIndividual],
) -> None:
    """Vérifie les fondateurs, familles, transmissions et classes attendues."""
    if panel.source_individual_count != len(source_individuals):
        _fail("pedigree_source_count_mismatch")
    ids = [individual.sample_id for individual in panel.individuals]
    if len(ids) != len(set(ids)):
        _fail("pedigree_sample_duplicate")
    founders = [individual for individual in panel.individuals if individual.is_founder]
    descendants = [
        individual for individual in panel.individuals if not individual.is_founder
    ]
    expected_origins = Counter(
        (source.sample_id, index)
        for source in source_individuals
        for index in (1, 2)
    )
    observed_origins = Counter(
        (origin.d12_sample_id, origin.d12_haplotype_index)
        for founder in founders
        for origin in (founder.haplotype_1, founder.haplotype_2)
    )
    if len(founders) != len(source_individuals) or observed_origins != expected_origins:
        _fail("pedigree_founder_haplotype_bijection_failed")
    if len(descendants) != 4 * panel.family_count:
        _fail("pedigree_descendant_count_mismatch")
    if len(panel.transmissions) != 2 * len(descendants):
        _fail("pedigree_transmission_count_mismatch")
    classes = {record.relationship_class for record in panel.kinship_records}
    if classes != {
        "UNRELATED",
        "FIRST_DEGREE",
        "SECOND_DEGREE",
        "THIRD_DEGREE",
    }:
        _fail("pedigree_relationship_class_coverage_invalid")
    by_id = {individual.sample_id: individual for individual in panel.individuals}
    observed_transmissions: Counter[tuple[str, int]] = Counter()
    for record in panel.transmissions:
        child = by_id.get(record.child_id)
        parent = by_id.get(record.parent_id)
        if child is None or parent is None or child.is_founder:
            _fail("pedigree_transmission_reference_invalid")
        expected_parent = (
            child.father_id
            if record.child_haplotype_index == 1
            else child.mother_id
        )
        if expected_parent != parent.sample_id:
            _fail("pedigree_transmission_parent_mismatch")
        parent_origin = (parent.haplotype_1, parent.haplotype_2)[
            record.parental_haplotype_index - 1
        ]
        child_origin = (child.haplotype_1, child.haplotype_2)[
            record.child_haplotype_index - 1
        ]
        if parent_origin != child_origin:
            _fail("pedigree_transmission_origin_mismatch")
        observed_transmissions[(record.child_id, record.child_haplotype_index)] += 1
    if any(value != 1 for value in observed_transmissions.values()):
        _fail("pedigree_transmission_duplicate")


def canonical_pedigree_bytes(panel: PedigreePanel) -> bytes:
    """Encode la structure familiale privée D13."""
    output = io.StringIO(newline="")
    writer = csv.writer(output, delimiter="\t", lineterminator="\n")
    writer.writerow(
        (
            "sample_id",
            "family_unit_id",
            "father_id",
            "mother_id",
            "sex",
            "generation",
            "is_founder",
        )
    )
    for individual in panel.individuals:
        writer.writerow(
            (
                individual.sample_id,
                individual.family_unit_id,
                individual.father_id or "0",
                individual.mother_id or "0",
                individual.sex,
                individual.generation,
                str(individual.is_founder).lower(),
            )
        )
    return output.getvalue().encode("utf-8")


def canonical_kinship_bytes(panel: PedigreePanel) -> bytes:
    """Encode les coefficients de parenté théoriques D13."""
    output = io.StringIO(newline="")
    writer = csv.writer(output, delimiter="\t", lineterminator="\n")
    writer.writerow(
        (
            "sample_1_id",
            "sample_2_id",
            "relationship_class",
            "kinship_coefficient",
            "independence_unit_id",
        )
    )
    for record in panel.kinship_records:
        writer.writerow(
            (
                record.sample_1_id,
                record.sample_2_id,
                record.relationship_class,
                f"{record.kinship_coefficient:.6g}",
                record.independence_unit_id,
            )
        )
    return output.getvalue().encode("utf-8")


def canonical_transmission_bytes(panel: PedigreePanel) -> bytes:
    """Encode les choix de copies parentales transmis sans recombinaison."""
    output = io.StringIO(newline="")
    writer = csv.writer(output, delimiter="\t", lineterminator="\n")
    writer.writerow(
        (
            "family_unit_id",
            "child_id",
            "parent_id",
            "parental_haplotype_index",
            "child_haplotype_index",
            "recombination_count",
        )
    )
    for record in panel.transmissions:
        writer.writerow(
            (
                record.family_unit_id,
                record.child_id,
                record.parent_id,
                record.parental_haplotype_index,
                record.child_haplotype_index,
                0,
            )
        )
    return output.getvalue().encode("utf-8")


def canonical_provenance_bytes(panel: PedigreePanel) -> bytes:
    """Encode la provenance D13→D12 de chaque copie haplotypique."""
    output = io.StringIO(newline="")
    writer = csv.writer(output, delimiter="\t", lineterminator="\n")
    writer.writerow(
        (
            "sample_id",
            "d13_haplotype_index",
            "d12_sample_id",
            "d12_haplotype_index",
        )
    )
    for individual in panel.individuals:
        for index, origin in enumerate(
            (individual.haplotype_1, individual.haplotype_2), start=1
        ):
            writer.writerow(
                (
                    individual.sample_id,
                    index,
                    origin.d12_sample_id,
                    origin.d12_haplotype_index,
                )
            )
    return output.getvalue().encode("utf-8")


def canonical_ancestry_bytes(panel: PedigreePanel) -> bytes:
    """Propage la vérité d'ascendance aux descendants D13."""
    output = io.StringIO(newline="")
    writer = csv.writer(output, delimiter="\t", lineterminator="\n")
    writer.writerow(("sample_id", "stratum_id", "mixture_proportions"))
    for individual in panel.individuals:
        writer.writerow(
            (
                individual.sample_id,
                individual.stratum_id,
                json.dumps(
                    dict(individual.mixture_proportions),
                    separators=(",", ":"),
                    sort_keys=True,
                ),
            )
        )
    return output.getvalue().encode("utf-8")


def _parse_phased_gt(sample_field: str, gt_index: int) -> tuple[int, int]:
    values = sample_field.split(":")
    if gt_index >= len(values):
        _fail("pedigree_vcf_sample_incomplete")
    alleles = values[gt_index].split("|")
    if len(alleles) != 2 or any(value not in {"0", "1"} for value in alleles):
        _fail("pedigree_vcf_genotype_invalid")
    return int(alleles[0]), int(alleles[1])


def transform_pedigree_vcf(
    source_lines: Iterable[str],
    dense_output: TextIO,
    observed_output: TextIO,
    *,
    panel: PedigreePanel,
    observed_variant_keys: Mapping[tuple[str, int, str, str], str],
    chromosome: int,
) -> PedigreeVcfStats:
    """Matérialise les génotypes D13 et contrôle chaque transmission."""
    if not panel.individuals or not observed_variant_keys or not 1 <= chromosome <= 22:
        _fail("pedigree_vcf_arguments_invalid")
    output_ids = [individual.sample_id for individual in panel.individuals]
    source_name = f"{PEDIGREE_METHOD_ID}_{PEDIGREE_METHOD_VERSION}"
    header = (
        "##fileformat=VCFv4.2\n"
        f"##source={source_name}\n"
        f"##contig=<ID=chr{chromosome}>\n"
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
            *output_ids,
        )
    )
    dense_output.write(header + columns + "\n")
    observed_output.write(header + columns + "\n")

    source_index: dict[str, int] | None = None
    dense_count = 0
    observed_count = 0
    observed_seen: set[tuple[str, int, str, str]] = set()
    mendelian_errors = 0
    source_founder_alt_count = 0
    output_founder_alt_count = 0
    founders = [individual for individual in panel.individuals if individual.is_founder]
    expected_source_ids = {
        origin.d12_sample_id
        for founder in founders
        for origin in (founder.haplotype_1, founder.haplotype_2)
    }
    for raw_line in source_lines:
        if not raw_line.strip() or raw_line.startswith("##"):
            continue
        if raw_line.startswith("#CHROM\t"):
            source_ids = raw_line.rstrip("\n\r").split("\t")[9:]
            if (
                len(source_ids) != len(set(source_ids))
                or set(source_ids) != expected_source_ids
            ):
                _fail("pedigree_vcf_sample_universe_mismatch")
            source_index = {
                sample_id: index for index, sample_id in enumerate(source_ids)
            }
            continue
        if source_index is None:
            _fail("pedigree_vcf_header_missing")
        fields = raw_line.rstrip("\n\r").split("\t")
        if len(fields) != 9 + len(source_index):
            _fail("pedigree_vcf_column_count_invalid")
        try:
            position_bp = int(fields[1])
        except ValueError:
            _fail("pedigree_vcf_position_invalid")
        if fields[0] not in {str(chromosome), f"chr{chromosome}"}:
            _fail("pedigree_vcf_chromosome_mismatch")
        format_fields = fields[8].split(":")
        if "GT" not in format_fields:
            _fail("pedigree_vcf_gt_missing")
        gt_index = format_fields.index("GT")
        parsed = {
            sample_id: _parse_phased_gt(fields[9 + index], gt_index)
            for sample_id, index in source_index.items()
        }
        output_genotypes: dict[str, tuple[int, int]] = {}
        phased: list[str] = []
        unphased: list[str] = []
        for individual in panel.individuals:
            left = parsed[individual.haplotype_1.d12_sample_id][
                individual.haplotype_1.d12_haplotype_index - 1
            ]
            right = parsed[individual.haplotype_2.d12_sample_id][
                individual.haplotype_2.d12_haplotype_index - 1
            ]
            output_genotypes[individual.sample_id] = (left, right)
            phased.append(f"{left}|{right}")
            ordered = sorted((left, right))
            unphased.append(f"{ordered[0]}/{ordered[1]}")
            if individual.is_founder:
                output_founder_alt_count += left + right
        source_founder_alt_count += sum(sum(genotype) for genotype in parsed.values())
        for transmission in panel.transmissions:
            child = output_genotypes[transmission.child_id]
            parent = output_genotypes[transmission.parent_id]
            if (
                child[transmission.child_haplotype_index - 1]
                != parent[transmission.parental_haplotype_index - 1]
            ):
                mendelian_errors += 1
        fixed = (
            f"chr{chromosome}",
            fields[1],
            fields[2],
            fields[3],
            fields[4],
            ".",
            "PASS",
            ".",
            "GT",
        )
        dense_output.write("\t".join((*fixed, *phased)) + "\n")
        dense_count += 1
        key = (f"chr{chromosome}", position_bp, fields[3], fields[4])
        if key in observed_variant_keys:
            observed_fixed = (
                f"chr{chromosome}",
                fields[1],
                observed_variant_keys[key],
                fields[3],
                fields[4],
                ".",
                "PASS",
                ".",
                "GT",
            )
            observed_output.write("\t".join((*observed_fixed, *unphased)) + "\n")
            observed_seen.add(key)
            observed_count += 1
    if source_index is None:
        _fail("pedigree_vcf_header_missing")
    if observed_seen != set(observed_variant_keys):
        _fail("pedigree_observed_variant_mismatch")
    founder_origins = Counter(
        (origin.d12_sample_id, origin.d12_haplotype_index)
        for founder in founders
        for origin in (founder.haplotype_1, founder.haplotype_2)
    )
    expected_origins = Counter(
        (sample_id, index)
        for sample_id in source_index
        for index in (1, 2)
    )
    return PedigreeVcfStats(
        dense_variant_count=dense_count,
        observed_variant_count=observed_count,
        source_individual_count=len(source_index),
        output_individual_count=len(panel.individuals),
        descendant_count=len(panel.individuals) - len(source_index),
        founder_haplotype_copy_error_count=int(founder_origins != expected_origins),
        founder_allele_count_difference=(
            output_founder_alt_count - source_founder_alt_count
        ),
        mendelian_error_count=mendelian_errors,
    )
