"""Construction D12 de strates pures et de mélanges diploïdes connus.

D12 ne simule pas encore de recombinaison. Un mélange individuel peut donc
être représenté exactement uniquement en associant un haplotype entier de
chacune de deux strates, soit une proportion 50/50 sur le domaine. Les
proportions inégales et les mélanges à plus de deux composantes nécessitent des
mosaïques chromosomiques et sont différés à D14.
"""

from __future__ import annotations

from collections import Counter, defaultdict
from collections.abc import Iterable, Mapping, Sequence
import csv
from dataclasses import dataclass
import io
import json
import random
import re
from typing import TextIO


POPULATION_STRUCTURE_METHOD_ID = "reciprocal_whole_haplotype_mixture_v1"
POPULATION_STRUCTURE_METHOD_VERSION = "1.0.0"
_STABLE_IDENTIFIER = re.compile(r"^[A-Za-z0-9][A-Za-z0-9._-]*$")


class PopulationStructureError(RuntimeError):
    """Signale une violation D12 avec un code stable et non sensible."""

    def __init__(self, code: str, detail: str | None = None) -> None:
        self.code = code
        self.detail = detail
        super().__init__(code if detail is None else f"{code}: {detail}")


@dataclass(frozen=True, order=True)
class PopulationSourceDiploid:
    """Diploïde synthétique D11 et sa strate publique privée."""

    sample_id: str
    stratum_id: str


@dataclass(frozen=True, order=True)
class PopulationHaplotypeOrigin:
    """Copie haplotypique D11 utilisée par un diploïde D12."""

    sample_id: str
    haplotype_index: int
    stratum_id: str


@dataclass(frozen=True)
class MixtureSpecification:
    """Mélange réciproque exact entre deux strates, à parts égales."""

    mixture_id: str
    stratum_1: str
    stratum_2: str
    individual_count: int


@dataclass(frozen=True)
class PopulationDiploid:
    """Diploïde D12 avec vérité d'ascendance connue."""

    sample_id: str
    haplotype_1: PopulationHaplotypeOrigin
    haplotype_2: PopulationHaplotypeOrigin
    stratum_id: str
    mixture_proportions: tuple[tuple[str, float], ...]


@dataclass(frozen=True)
class PopulationVcfStats:
    """Contrôles agrégés de la matérialisation VCF D12."""

    dense_variant_count: int
    observed_variant_count: int
    individual_count: int
    source_haplotype_count: int
    output_haplotype_count: int
    haplotype_copy_error_count: int
    allele_count_difference: int
    genotype_haplotype_discordance_count: int


def _fail(code: str, detail: str | None = None) -> None:
    raise PopulationStructureError(code, detail)


def _validate_identifier(value: object, namespace: str) -> str:
    if not isinstance(value, str) or _STABLE_IDENTIFIER.fullmatch(value) is None:
        _fail("population_structure_identifier_invalid", namespace)
    return value


def parse_d11_population_provenance(
    lines: Iterable[str],
) -> tuple[PopulationSourceDiploid, ...]:
    """Lit la provenance D11 sans propager les identifiants donneurs.

    Les deux copies d'un échantillon D11 doivent appartenir à la même strate.
    Les colonnes donneur sont exigées par le contrat D11 mais volontairement
    ignorées dans l'objet retourné et dans tous les artefacts D12.
    """
    reader = csv.DictReader(lines, delimiter="\t")
    required = {
        "synthetic_identifier",
        "synthetic_haplotype_index",
        "source_identifier",
        "source_haplotype_index",
        "source_population",
    }
    if reader.fieldnames is None or set(reader.fieldnames) != required:
        _fail("population_structure_provenance_columns_invalid")
    rows: dict[str, dict[int, str]] = defaultdict(dict)
    for record in reader:
        sample_id = _validate_identifier(
            record["synthetic_identifier"], "d11_sample"
        )
        stratum_id = _validate_identifier(record["source_population"], "stratum")
        try:
            haplotype_index = int(record["synthetic_haplotype_index"])
            source_haplotype_index = int(record["source_haplotype_index"])
        except (TypeError, ValueError):
            _fail("population_structure_provenance_haplotype_invalid")
        if haplotype_index not in (1, 2) or source_haplotype_index not in (1, 2):
            _fail("population_structure_provenance_haplotype_invalid")
        if haplotype_index in rows[sample_id]:
            _fail("population_structure_provenance_haplotype_duplicate")
        rows[sample_id][haplotype_index] = stratum_id
    if not rows:
        _fail("population_structure_provenance_empty")
    diploids: list[PopulationSourceDiploid] = []
    for sample_id in sorted(rows, key=str.casefold):
        assignments = rows[sample_id]
        if set(assignments) != {1, 2}:
            _fail("population_structure_provenance_haplotype_missing")
        if assignments[1] != assignments[2]:
            _fail("population_structure_d11_stratum_mismatch")
        diploids.append(PopulationSourceDiploid(sample_id, assignments[1]))
    return tuple(diploids)


def parse_mixture_specifications(
    config: Mapping[str, object],
) -> tuple[int, tuple[MixtureSpecification, ...]]:
    """Valide la partie opérationnelle d'une configuration D12 finalisée."""
    if (
        set(config)
        != {
            "schema_version",
            "population_structure_config_id",
            "config_state",
            "method",
            "mixtures",
        }
        or config.get("schema_version") != "1.0.0"
        or config.get("config_state") != "FINALIZED"
    ):
        _fail("population_structure_config_invalid")
    _validate_identifier(
        config.get("population_structure_config_id"), "configuration"
    )
    method = config.get("method")
    mixtures = config.get("mixtures")
    if not isinstance(method, Mapping) or not isinstance(mixtures, list):
        _fail("population_structure_config_invalid")
    if (
        set(method)
        != {
            "method_id",
            "method_version",
            "random_seed",
            "pure_policy",
            "mixture_unit",
            "supported_proportions",
        }
        or method.get("pure_policy") != "RETAIN_UNUSED_D11_DIPLOIDS"
        or method.get("method_id") != POPULATION_STRUCTURE_METHOD_ID
        or method.get("method_version") != POPULATION_STRUCTURE_METHOD_VERSION
        or method.get("mixture_unit") != "RECIPROCAL_WHOLE_HAPLOTYPE_PAIR"
        or method.get("supported_proportions") != "PURE_OR_TWO_WAY_EQUAL"
    ):
        _fail("population_structure_method_invalid")
    random_seed = method.get("random_seed")
    if (
        type(random_seed) is not int
        or random_seed < 0
        or random_seed > 4_294_967_295
    ):
        _fail("population_structure_seed_invalid")

    specifications: list[MixtureSpecification] = []
    mixture_ids: set[str] = set()
    stratum_pairs: set[tuple[str, str]] = set()
    for entry in mixtures:
        if not isinstance(entry, Mapping) or set(entry) != {
            "mixture_id",
            "individual_count",
            "components",
        }:
            _fail("population_structure_mixture_invalid")
        mixture_id = _validate_identifier(entry.get("mixture_id"), "mixture")
        components = entry.get("components")
        individual_count = entry.get("individual_count")
        if (
            not isinstance(components, list)
            or len(components) != 2
            or type(individual_count) is not int
            or individual_count < 2
            or individual_count % 2
        ):
            _fail("population_structure_mixture_invalid")
        parsed_components: list[tuple[str, float]] = []
        for component in components:
            if not isinstance(component, Mapping) or set(component) != {
                "stratum_id",
                "proportion",
            }:
                _fail("population_structure_mixture_invalid")
            stratum_id = _validate_identifier(
                component.get("stratum_id"), "mixture_stratum"
            )
            proportion = component.get("proportion")
            if type(proportion) not in (int, float) or float(proportion) != 0.5:
                _fail("population_structure_mixture_proportion_unsupported")
            parsed_components.append((stratum_id, float(proportion)))
        if parsed_components[0][0] == parsed_components[1][0]:
            _fail("population_structure_mixture_strata_not_distinct")
        pair = tuple(sorted((parsed_components[0][0], parsed_components[1][0])))
        if mixture_id.casefold() in mixture_ids or pair in stratum_pairs:
            _fail("population_structure_mixture_duplicate")
        mixture_ids.add(mixture_id.casefold())
        stratum_pairs.add(pair)
        specifications.append(
            MixtureSpecification(
                mixture_id=mixture_id,
                stratum_1=parsed_components[0][0],
                stratum_2=parsed_components[1][0],
                individual_count=individual_count,
            )
        )
    return random_seed, tuple(specifications)


def build_population_panel(
    source_diploids: Sequence[PopulationSourceDiploid],
    mixture_specifications: Sequence[MixtureSpecification],
    *,
    random_seed: int,
    synthetic_identifier_prefix: str,
) -> tuple[PopulationDiploid, ...]:
    """Réapparie des diploïdes D11 sans modifier ni dupliquer un haplotype.

    Pour chaque paire de sources A/B, deux sorties réciproques sont produites :
    ``A1|B2`` et ``B1|A2``. Les sources non sélectionnées restent pures. Ainsi,
    chaque copie D11 est utilisée exactement une fois et les fréquences
    alléliques/LD haplotypiques du pool sont conservés par permutation.
    """
    if (
        type(random_seed) is not int
        or random_seed < 0
        or random_seed > 4_294_967_295
    ):
        _fail("population_structure_seed_invalid")
    if _STABLE_IDENTIFIER.fullmatch(f"{synthetic_identifier_prefix}1") is None:
        _fail("population_structure_output_prefix_invalid")
    if not source_diploids:
        _fail("population_structure_source_empty")
    source_ids = [diploid.sample_id for diploid in source_diploids]
    if len({value.casefold() for value in source_ids}) != len(source_ids):
        _fail("population_structure_source_duplicate")
    grouped: dict[str, list[PopulationSourceDiploid]] = defaultdict(list)
    for diploid in source_diploids:
        _validate_identifier(diploid.sample_id, "d11_sample")
        _validate_identifier(diploid.stratum_id, "stratum")
        grouped[diploid.stratum_id].append(diploid)

    mixture_ids: set[str] = set()
    mixture_pairs: set[tuple[str, str]] = set()
    for specification in mixture_specifications:
        _validate_identifier(specification.mixture_id, "mixture")
        _validate_identifier(specification.stratum_1, "mixture_stratum")
        _validate_identifier(specification.stratum_2, "mixture_stratum")
        pair = tuple(sorted((specification.stratum_1, specification.stratum_2)))
        if (
            specification.stratum_1 == specification.stratum_2
            or type(specification.individual_count) is not int
            or specification.individual_count < 2
            or specification.individual_count % 2
            or specification.mixture_id.casefold() in mixture_ids
            or pair in mixture_pairs
        ):
            _fail("population_structure_mixture_invalid")
        mixture_ids.add(specification.mixture_id.casefold())
        mixture_pairs.add(pair)

    generator = random.Random(random_seed)
    for members in grouped.values():
        members.sort(key=lambda value: value.sample_id.casefold())
        generator.shuffle(members)
    cursors = {stratum_id: 0 for stratum_id in grouped}
    unlabelled: list[
        tuple[
            PopulationHaplotypeOrigin,
            PopulationHaplotypeOrigin,
            str,
            tuple[tuple[str, float], ...],
        ]
    ] = []
    for specification in mixture_specifications:
        if (
            specification.stratum_1 not in grouped
            or specification.stratum_2 not in grouped
        ):
            _fail("population_structure_mixture_stratum_missing")
        source_pair_count = specification.individual_count // 2
        selected: list[list[PopulationSourceDiploid]] = []
        for stratum_id in (specification.stratum_1, specification.stratum_2):
            start = cursors[stratum_id]
            end = start + source_pair_count
            if end > len(grouped[stratum_id]):
                _fail("population_structure_mixture_capacity_insufficient")
            selected.append(grouped[stratum_id][start:end])
            cursors[stratum_id] = end
        proportions = tuple(
            sorted(
                ((specification.stratum_1, 0.5), (specification.stratum_2, 0.5))
            )
        )
        for left, right in zip(*selected, strict=True):
            unlabelled.extend(
                (
                    (
                        PopulationHaplotypeOrigin(left.sample_id, 1, left.stratum_id),
                        PopulationHaplotypeOrigin(right.sample_id, 2, right.stratum_id),
                        specification.mixture_id,
                        proportions,
                    ),
                    (
                        PopulationHaplotypeOrigin(right.sample_id, 1, right.stratum_id),
                        PopulationHaplotypeOrigin(left.sample_id, 2, left.stratum_id),
                        specification.mixture_id,
                        proportions,
                    ),
                )
            )

    for stratum_id in sorted(grouped, key=str.casefold):
        for source in grouped[stratum_id][cursors[stratum_id] :]:
            unlabelled.append(
                (
                    PopulationHaplotypeOrigin(source.sample_id, 1, stratum_id),
                    PopulationHaplotypeOrigin(source.sample_id, 2, stratum_id),
                    stratum_id,
                    ((stratum_id, 1.0),),
                )
            )
    generator.shuffle(unlabelled)
    width = max(6, len(str(len(unlabelled))))
    panel = tuple(
        PopulationDiploid(
            sample_id=f"{synthetic_identifier_prefix}{index:0{width}d}",
            haplotype_1=left,
            haplotype_2=right,
            stratum_id=stratum_id,
            mixture_proportions=proportions,
        )
        for index, (left, right, stratum_id, proportions) in enumerate(
            unlabelled, start=1
        )
    )
    validate_population_panel(panel, source_diploids)
    return panel


def validate_population_panel(
    panel: Sequence[PopulationDiploid],
    source_diploids: Sequence[PopulationSourceDiploid],
) -> None:
    """Contrôle bijection haplotypique, sommes d'ascendance et identifiants."""
    expected = Counter(
        (source.sample_id, haplotype_index)
        for source in source_diploids
        for haplotype_index in (1, 2)
    )
    observed: Counter[tuple[str, int]] = Counter()
    sample_ids: set[str] = set()
    for diploid in panel:
        _validate_identifier(diploid.sample_id, "d12_sample")
        _validate_identifier(diploid.stratum_id, "d12_stratum")
        if diploid.sample_id.casefold() in sample_ids:
            _fail("population_structure_output_identifier_duplicate")
        sample_ids.add(diploid.sample_id.casefold())
        total = sum(proportion for _, proportion in diploid.mixture_proportions)
        if abs(total - 1.0) > 1e-12 or any(
            proportion <= 0 or proportion > 1
            for _, proportion in diploid.mixture_proportions
        ):
            _fail("population_structure_proportion_sum_invalid")
        origin_strata = Counter(
            origin.stratum_id for origin in (diploid.haplotype_1, diploid.haplotype_2)
        )
        expected_proportions = tuple(
            sorted((stratum, count / 2) for stratum, count in origin_strata.items())
        )
        if diploid.mixture_proportions != expected_proportions:
            _fail("population_structure_proportion_origin_mismatch")
        for origin in (diploid.haplotype_1, diploid.haplotype_2):
            if origin.haplotype_index not in (1, 2):
                _fail("population_structure_haplotype_index_invalid")
            observed[(origin.sample_id, origin.haplotype_index)] += 1
    if len(panel) != len(source_diploids) or observed != expected:
        _fail("population_structure_haplotype_bijection_failed")


def canonical_ancestry_truth_bytes(panel: Sequence[PopulationDiploid]) -> bytes:
    """Encode la vérité d'ascendance D12 dans un TSV canonique."""
    output = io.StringIO(newline="")
    writer = csv.writer(output, delimiter="\t", lineterminator="\n")
    writer.writerow(("sample_id", "stratum_id", "mixture_proportions"))
    for diploid in sorted(panel, key=lambda value: value.sample_id):
        proportions = {key: value for key, value in diploid.mixture_proportions}
        writer.writerow(
            (
                diploid.sample_id,
                diploid.stratum_id,
                json.dumps(proportions, separators=(",", ":"), sort_keys=True),
            )
        )
    return output.getvalue().encode("utf-8")


def canonical_population_provenance_bytes(
    panel: Sequence[PopulationDiploid],
) -> bytes:
    """Encode la liaison privée D12→D11 sans identifiant donneur public."""
    output = io.StringIO(newline="")
    writer = csv.writer(output, delimiter="\t", lineterminator="\n")
    writer.writerow(
        (
            "sample_id",
            "output_haplotype_index",
            "d11_sample_id",
            "d11_haplotype_index",
            "source_stratum_id",
        )
    )
    for diploid in sorted(panel, key=lambda value: value.sample_id):
        for output_index, origin in enumerate(
            (diploid.haplotype_1, diploid.haplotype_2), start=1
        ):
            writer.writerow(
                (
                    diploid.sample_id,
                    output_index,
                    origin.sample_id,
                    origin.haplotype_index,
                    origin.stratum_id,
                )
            )
    return output.getvalue().encode("utf-8")


def _parse_phased_gt(sample_field: str, gt_index: int) -> tuple[int, int]:
    values = sample_field.split(":")
    if gt_index >= len(values):
        _fail("population_structure_vcf_sample_incomplete")
    genotype = values[gt_index]
    alleles = genotype.split("|")
    if len(alleles) != 2 or any(allele not in {"0", "1"} for allele in alleles):
        _fail("population_structure_vcf_genotype_invalid")
    return int(alleles[0]), int(alleles[1])


def transform_population_vcf(
    source_lines: Iterable[str],
    dense_phased_output: TextIO,
    observed_unphased_output: TextIO,
    *,
    panel: Sequence[PopulationDiploid],
    observed_variant_keys: Mapping[tuple[str, int, str, str], str],
    chromosome: int,
) -> PopulationVcfStats:
    """Réapparie le VCF phasé D11 et extrait l'univers observable D11b."""
    if not panel or not observed_variant_keys or not 1 <= chromosome <= 22:
        _fail("population_structure_vcf_arguments_invalid")
    output_ids = [diploid.sample_id for diploid in panel]
    source_name = (
        f"{POPULATION_STRUCTURE_METHOD_ID}_"
        f"{POPULATION_STRUCTURE_METHOD_VERSION}"
    )
    common_header = (
        "##fileformat=VCFv4.2\n"
        f"##source={source_name}\n"
        f"##contig=<ID=chr{chromosome}>\n"
        '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n'
    )
    column_header = "\t".join(
        (
            "#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER",
            "INFO", "FORMAT", *output_ids,
        )
    )
    dense_phased_output.write(common_header + column_header + "\n")
    observed_unphased_output.write(common_header + column_header + "\n")

    source_index: dict[str, int] | None = None
    dense_count = 0
    observed_count = 0
    source_alt_count = 0
    output_alt_count = 0
    genotype_discordance_count = 0
    observed_seen: set[tuple[str, int, str, str]] = set()
    for raw_line in source_lines:
        if not raw_line.strip() or raw_line.startswith("##"):
            continue
        if raw_line.startswith("#CHROM\t"):
            columns = raw_line.rstrip("\n\r").split("\t")
            source_ids = columns[9:]
            if len(source_ids) != len(set(source_ids)):
                _fail("population_structure_vcf_source_duplicate")
            source_index = {
                sample_id: index for index, sample_id in enumerate(source_ids)
            }
            required = {
                origin.sample_id
                for diploid in panel
                for origin in (diploid.haplotype_1, diploid.haplotype_2)
            }
            if set(source_ids) != required:
                _fail("population_structure_vcf_sample_universe_mismatch")
            continue
        if source_index is None:
            _fail("population_structure_vcf_header_missing")
        fields = raw_line.rstrip("\n\r").split("\t")
        if len(fields) != 9 + len(source_index):
            _fail("population_structure_vcf_column_count_invalid")
        chrom = fields[0]
        try:
            position_bp = int(fields[1])
        except ValueError:
            _fail("population_structure_vcf_position_invalid")
        if chrom not in {str(chromosome), f"chr{chromosome}"}:
            _fail("population_structure_vcf_chromosome_mismatch")
        format_fields = fields[8].split(":")
        if "GT" not in format_fields:
            _fail("population_structure_vcf_gt_missing")
        gt_index = format_fields.index("GT")
        parsed = {
            sample_id: _parse_phased_gt(fields[9 + index], gt_index)
            for sample_id, index in source_index.items()
        }
        phased_genotypes: list[str] = []
        unphased_genotypes: list[str] = []
        for diploid in panel:
            left = parsed[diploid.haplotype_1.sample_id][
                diploid.haplotype_1.haplotype_index - 1
            ]
            right = parsed[diploid.haplotype_2.sample_id][
                diploid.haplotype_2.haplotype_index - 1
            ]
            source_alt_count += left + right
            output_alt_count += left + right
            phased_genotypes.append(f"{left}|{right}")
            ordered = sorted((left, right))
            unphased_genotypes.append(f"{ordered[0]}/{ordered[1]}")
            if Counter((left, right)) != Counter(ordered):
                genotype_discordance_count += 1
        dense_fixed = (
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
        dense_phased_output.write(
            "\t".join((*dense_fixed, *phased_genotypes)) + "\n"
        )
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
            observed_unphased_output.write(
                "\t".join((*observed_fixed, *unphased_genotypes)) + "\n"
            )
            observed_seen.add(key)
            observed_count += 1
    if source_index is None:
        _fail("population_structure_vcf_header_missing")
    if observed_seen != set(observed_variant_keys):
        _fail("population_structure_observed_variant_mismatch")
    return PopulationVcfStats(
        dense_variant_count=dense_count,
        observed_variant_count=observed_count,
        individual_count=len(panel),
        source_haplotype_count=2 * len(panel),
        output_haplotype_count=2 * len(panel),
        haplotype_copy_error_count=0,
        allele_count_difference=output_alt_count - source_alt_count,
        genotype_haplotype_discordance_count=genotype_discordance_count,
    )


def parse_observed_variant_keys(
    lines: Iterable[str],
) -> dict[tuple[str, int, str, str], str]:
    """Extrait les clés D11b et conserve l'identifiant de sonde observable."""
    keys: dict[tuple[str, int, str, str], str] = {}
    header_seen = False
    for raw_line in lines:
        if raw_line.startswith("##") or not raw_line.strip():
            continue
        if raw_line.startswith("#CHROM\t"):
            header_seen = True
            continue
        if not header_seen:
            _fail("population_structure_observed_header_missing")
        fields = raw_line.rstrip("\n\r").split("\t")
        if len(fields) < 9:
            _fail("population_structure_observed_row_invalid")
        try:
            key = (fields[0], int(fields[1]), fields[3], fields[4])
        except ValueError:
            _fail("population_structure_observed_row_invalid")
        if key in keys:
            _fail("population_structure_observed_variant_duplicate")
        keys[key] = fields[2]
    if not header_seen or not keys:
        _fail("population_structure_observed_empty")
    return keys
