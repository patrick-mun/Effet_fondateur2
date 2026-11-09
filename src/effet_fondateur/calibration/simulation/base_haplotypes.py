"""Construction D11 de diploïdes synthétiques par réappariement d'haplotypes.

La méthode permute des haplotypes publics phasés sans modifier leurs allèles.
Chaque copie donneuse est utilisée exactement une fois et les deux copies d'un
diploïde synthétique proviennent de donneurs distincts de la même strate
publique. Cette bijection conserve exactement les fréquences alléliques et le
LD haplotypique du pool sélectionné. Elle ne simule aucune transmission
parent-enfant : ce contrôle mendélien appartient à D13.
"""

from __future__ import annotations

from bisect import bisect_right
from collections import Counter, defaultdict
from collections.abc import Iterable, Mapping, Sequence
import csv
from dataclasses import dataclass
import hashlib
import io
import json
import random
import re
from typing import Any, TextIO
import unicodedata


BASE_HAPLOTYPE_METHOD_ID = "population_stratified_whole_haplotype_repairing_v1"
BASE_HAPLOTYPE_METHOD_VERSION = "1.0.0"
_STABLE_IDENTIFIER = re.compile(r"^[A-Za-z0-9][A-Za-z0-9._-]*$")
_SHA256 = re.compile(r"^[0-9a-f]{64}$")
_ALLELE = re.compile(r"^[ACGTN]$")
_EXPECTED_CHECKS = (
    "SOURCE_HASHES_MATCHED",
    "PARTITION_BINDING_MATCHED",
    "HAPLOTYPE_COPY_BIJECTION",
    "SOURCE_PAIRS_DISJOINT",
    "ALLELE_FREQUENCIES_EXACT",
    "HAPLOTYPE_LD_EXACT_BY_PERMUTATION",
    "GENOTYPE_HAPLOTYPE_CONCORDANCE_100_PERCENT",
    "GENETIC_MAP_COORDINATES_VALID",
    "RAW_DONOR_IDENTIFIERS_NOT_IN_MANIFEST",
)
_EXPECTED_CONFIDENTIALITY = {
    "raw_donor_identifiers_in_manifest": False,
    "local_source_paths_in_manifest": False,
    "donor_provenance_private": True,
    "artifacts_version_control_allowed": False,
    "aggregate_counts_publishable": True,
}


class BaseHaplotypeError(RuntimeError):
    """Signale une violation D11 avec un code stable et non sensible."""

    def __init__(self, code: str, detail: str | None = None) -> None:
        self.code = code
        self.detail = detail
        message = code if detail is None else f"{code}: {detail}"
        super().__init__(message)


@dataclass(frozen=True, order=True)
class HaplotypeOrigin:
    """Origine privée d'une copie haplotypique publique."""

    source_identifier: str
    source_haplotype_index: int
    source_population: str


@dataclass(frozen=True)
class SyntheticDiploid:
    """Deux haplotypes vrais affectés à un identifiant synthétique."""

    synthetic_identifier: str
    haplotype_1: HaplotypeOrigin
    haplotype_2: HaplotypeOrigin
    source_population: str


@dataclass(frozen=True)
class GeneticMapPoint:
    """Point d'une carte génétique monotone, en bp et cM."""

    chromosome: int
    position_bp: int
    position_cm: float


@dataclass(frozen=True)
class VcfTransformationStats:
    """Mesures agrégées de la transformation VCF D11."""

    variant_count: int
    synthetic_individual_count: int
    source_haplotype_count: int
    synthetic_haplotype_count: int
    allele_count_difference: int
    genotype_haplotype_discordance_count: int
    genetic_map_coordinate_error_count: int
    first_position_bp: int
    last_position_bp: int


@dataclass(frozen=True)
class BaseHaplotypeAudit:
    """Résumé publiable d'un manifest D11 validé."""

    base_haplotype_manifest_id: str
    panel_role: str
    individual_count: int
    haplotype_count: int
    variant_count: int
    population_stratum_count: int
    checks: tuple[str, ...]


def _fail(code: str, detail: str | None = None) -> None:
    raise BaseHaplotypeError(code, detail)


def _normalize_identifier(value: Any, namespace: str) -> tuple[str, str]:
    if not isinstance(value, str) or not value:
        _fail("base_haplotype_identifier_invalid", namespace)
    normalized = unicodedata.normalize("NFC", value)
    if normalized != normalized.strip() or any(
        character in normalized for character in ("\x00", "\n", "\r", "\t")
    ):
        _fail("base_haplotype_identifier_invalid", namespace)
    if any(ord(character) < 32 for character in normalized):
        _fail("base_haplotype_identifier_invalid", namespace)
    return normalized, normalized.casefold()


def _index_identifiers(values: Iterable[str], namespace: str) -> dict[str, str]:
    if isinstance(values, (str, bytes)):
        _fail("base_haplotype_identifier_collection_invalid", namespace)
    indexed: dict[str, str] = {}
    try:
        iterator = iter(values)
    except TypeError:
        _fail("base_haplotype_identifier_collection_invalid", namespace)
    for value in iterator:
        normalized, collision_key = _normalize_identifier(value, namespace)
        if collision_key in indexed:
            _fail("base_haplotype_identifier_duplicate", namespace)
        indexed[collision_key] = normalized
    return indexed


def parse_population_assignments(lines: Iterable[str]) -> dict[str, str]:
    """Lit les codes de population publics associés aux donneurs 1000G.

    La table doit contenir exactement une ligne par ``SampleID`` et les colonnes
    ``SampleID`` et ``Population``. Les codes servent seulement à éviter un
    réappariement inter-strates prématuré avant D12.
    """
    materialized_lines = [line for line in lines if line.strip()]
    if not materialized_lines:
        _fail("population_assignment_empty")
    fieldnames = materialized_lines[0].strip().split()
    if not {"SampleID", "Population"} <= set(fieldnames):
        _fail("population_assignment_columns_invalid")
    assignments: dict[str, str] = {}
    collision_index: dict[str, str] = {}
    for line in materialized_lines[1:]:
        values = line.strip().split()
        if len(values) != len(fieldnames):
            _fail("population_assignment_row_invalid")
        record = dict(zip(fieldnames, values, strict=True))
        identifier, collision_key = _normalize_identifier(
            record["SampleID"], "population_assignment"
        )
        population, _ = _normalize_identifier(
            record["Population"], "population_code"
        )
        if collision_key in collision_index:
            _fail("population_assignment_duplicate")
        collision_index[collision_key] = identifier
        assignments[identifier] = population
    if not assignments:
        _fail("population_assignment_empty")
    return assignments


def pair_donor_haplotypes(
    donor_identifiers: Iterable[str],
    population_by_identifier: Mapping[str, str],
    *,
    random_seed: int,
    synthetic_identifier_prefix: str,
) -> tuple[SyntheticDiploid, ...]:
    """Réapparie toutes les copies donneuses en diploïdes synthétiques.

    Chaque donneur fournit ses haplotypes 1 et 2 une seule fois. Dans chaque
    population, une rotation non nulle garantit que les deux copies d'un
    diploïde proviennent de donneurs différents. La graine explicite contrôle
    l'ordre, l'orientation des copies et les identifiants synthétiques.
    """
    if not isinstance(random_seed, int) or isinstance(random_seed, bool) or not (
        0 <= random_seed <= 4_294_967_295
    ):
        _fail("base_haplotype_seed_invalid")
    if _STABLE_IDENTIFIER.fullmatch(f"{synthetic_identifier_prefix}1") is None:
        _fail("synthetic_identifier_prefix_invalid")
    donors = _index_identifiers(donor_identifiers, "donor")
    if len(donors) < 2:
        _fail("donor_count_insufficient")
    if not isinstance(population_by_identifier, Mapping):
        _fail("population_assignment_mapping_invalid")

    population_index: dict[str, str] = {}
    for raw_identifier, raw_population in population_by_identifier.items():
        identifier, collision_key = _normalize_identifier(
            raw_identifier, "population_assignment"
        )
        population, _ = _normalize_identifier(raw_population, "population_code")
        if collision_key in population_index:
            _fail("population_assignment_duplicate")
        population_index[collision_key] = population
        if collision_key not in donors:
            continue
        if identifier.casefold() != donors[collision_key].casefold():
            _fail("population_assignment_identifier_mismatch")
    missing = donors.keys() - population_index.keys()
    if missing:
        _fail("population_assignment_missing_for_donor")

    grouped_donors: dict[str, list[str]] = defaultdict(list)
    for collision_key in sorted(donors):
        grouped_donors[population_index[collision_key]].append(donors[collision_key])
    singleton_populations = [
        population
        for population, members in grouped_donors.items()
        if len(members) < 2
    ]
    if singleton_populations:
        _fail("population_donor_count_insufficient")

    generator = random.Random(random_seed)
    unlabelled_pairs: list[tuple[HaplotypeOrigin, HaplotypeOrigin, str]] = []
    for population in sorted(grouped_donors, key=str.casefold):
        members = sorted(grouped_donors[population], key=str.casefold)
        generator.shuffle(members)
        offset = generator.randrange(1, len(members))
        for index, left_identifier in enumerate(members):
            right_identifier = members[(index + offset) % len(members)]
            left = HaplotypeOrigin(left_identifier, 1, population)
            right = HaplotypeOrigin(right_identifier, 2, population)
            if generator.getrandbits(1):
                left, right = right, left
            unlabelled_pairs.append((left, right, population))
    generator.shuffle(unlabelled_pairs)

    width = max(6, len(str(len(unlabelled_pairs))))
    pairings = tuple(
        SyntheticDiploid(
            synthetic_identifier=(
                f"{synthetic_identifier_prefix}{index:0{width}d}"
            ),
            haplotype_1=left,
            haplotype_2=right,
            source_population=population,
        )
        for index, (left, right, population) in enumerate(unlabelled_pairs, start=1)
    )
    validate_pairings(pairings, donors.values(), population_by_identifier)
    return pairings


def validate_pairings(
    pairings: Sequence[SyntheticDiploid],
    donor_identifiers: Iterable[str],
    population_by_identifier: Mapping[str, str],
) -> None:
    """Vérifie la bijection des copies et l'absence de paire d'un même donneur."""
    donors = _index_identifiers(donor_identifiers, "donor")
    if len(pairings) != len(donors):
        _fail("synthetic_individual_count_mismatch")
    synthetic_ids = _index_identifiers(
        (pairing.synthetic_identifier for pairing in pairings), "synthetic"
    )
    if len(synthetic_ids) != len(pairings):
        _fail("synthetic_identifier_count_mismatch")

    normalized_populations: dict[str, str] = {}
    for identifier, population in population_by_identifier.items():
        _, collision_key = _normalize_identifier(identifier, "population_assignment")
        normalized_population, _ = _normalize_identifier(
            population, "population_code"
        )
        normalized_populations[collision_key] = normalized_population

    observed: Counter[tuple[str, int]] = Counter()
    for pairing in pairings:
        origins = (pairing.haplotype_1, pairing.haplotype_2)
        if origins[0].source_identifier.casefold() == origins[1].source_identifier.casefold():
            _fail("same_source_pair_forbidden")
        for origin in origins:
            _, donor_key = _normalize_identifier(origin.source_identifier, "origin")
            if donor_key not in donors:
                _fail("pairing_origin_unknown")
            if origin.source_haplotype_index not in (1, 2):
                _fail("source_haplotype_index_invalid")
            if normalized_populations.get(donor_key) != origin.source_population:
                _fail("pairing_population_mismatch")
            if pairing.source_population != origin.source_population:
                _fail("synthetic_population_mismatch")
            observed[(donor_key, origin.source_haplotype_index)] += 1
    expected = Counter(
        (donor_key, haplotype_index)
        for donor_key in donors
        for haplotype_index in (1, 2)
    )
    if observed != expected:
        _fail("haplotype_copy_bijection_failed")


def canonical_provenance_bytes(pairings: Sequence[SyntheticDiploid]) -> bytes:
    """Encode la provenance privée D11 dans un TSV canonique."""
    output = io.StringIO(newline="")
    writer = csv.writer(output, delimiter="\t", lineterminator="\n")
    writer.writerow(
        (
            "synthetic_identifier",
            "synthetic_haplotype_index",
            "source_identifier",
            "source_haplotype_index",
            "source_population",
        )
    )
    for pairing in sorted(pairings, key=lambda value: value.synthetic_identifier):
        for synthetic_index, origin in enumerate(
            (pairing.haplotype_1, pairing.haplotype_2), start=1
        ):
            writer.writerow(
                (
                    pairing.synthetic_identifier,
                    synthetic_index,
                    origin.source_identifier,
                    origin.source_haplotype_index,
                    origin.source_population,
                )
            )
    return output.getvalue().encode("utf-8")


def parse_genetic_map(
    lines: Iterable[str], *, expected_chromosome: int
) -> tuple[GeneticMapPoint, ...]:
    """Parse une carte D09 et impose des coordonnées bp/cM monotones."""
    if not isinstance(expected_chromosome, int) or not 1 <= expected_chromosome <= 22:
        _fail("genetic_map_chromosome_invalid")
    materialized_lines = [line for line in lines if line.strip()]
    if not materialized_lines:
        _fail("genetic_map_empty")
    reader = csv.DictReader(materialized_lines, delimiter="\t")
    required = {"CHROMOSOME", "POSITION_BP", "POSITION_CM"}
    if reader.fieldnames is None or not required <= set(reader.fieldnames):
        _fail("genetic_map_columns_invalid")
    points: list[GeneticMapPoint] = []
    for row in reader:
        try:
            chromosome = int(row["CHROMOSOME"])
            position_bp = int(row["POSITION_BP"])
            position_cm = float(row["POSITION_CM"])
        except (TypeError, ValueError):
            _fail("genetic_map_row_invalid")
        if chromosome != expected_chromosome or position_bp < 1 or position_cm < 0:
            _fail("genetic_map_row_invalid")
        if points and (
            position_bp <= points[-1].position_bp
            or position_cm < points[-1].position_cm
        ):
            _fail("genetic_map_not_monotonic")
        points.append(GeneticMapPoint(chromosome, position_bp, position_cm))
    if len(points) < 2:
        _fail("genetic_map_point_count_insufficient")
    return tuple(points)


def _interpolate_cm(
    position_bp: int,
    points: Sequence[GeneticMapPoint],
    positions_bp: Sequence[int],
) -> float:
    if position_bp < positions_bp[0] or position_bp > positions_bp[-1]:
        _fail("genetic_map_extrapolation_forbidden")
    right_index = bisect_right(positions_bp, position_bp)
    if right_index == 0:
        _fail("genetic_map_extrapolation_forbidden")
    if positions_bp[right_index - 1] == position_bp:
        return points[right_index - 1].position_cm
    if right_index == len(points):
        _fail("genetic_map_extrapolation_forbidden")
    left = points[right_index - 1]
    right = points[right_index]
    fraction = (position_bp - left.position_bp) / (
        right.position_bp - left.position_bp
    )
    return left.position_cm + fraction * (right.position_cm - left.position_cm)


def _parse_gt(sample_field: str, gt_index: int) -> tuple[int, int]:
    values = sample_field.split(":")
    if gt_index >= len(values):
        _fail("vcf_sample_format_incomplete")
    genotype = values[gt_index]
    if genotype.count("|") != 1 or "/" in genotype:
        _fail("vcf_genotype_not_phased_diploid")
    alleles = genotype.split("|")
    if len(alleles) != 2 or any(allele not in {"0", "1"} for allele in alleles):
        _fail("vcf_genotype_not_complete_biallelic")
    return int(alleles[0]), int(alleles[1])


def transform_phased_vcf(
    source_lines: Iterable[str],
    phased_output: TextIO,
    unphased_output: TextIO,
    genetic_map_output: TextIO,
    *,
    pairings: Sequence[SyntheticDiploid],
    genetic_map_points: Sequence[GeneticMapPoint],
    chromosome: int,
    start_bp: int,
    end_bp: int,
) -> VcfTransformationStats:
    """Transforme un VCF donneur filtré en vérité phasée et génotypes dérivés.

    Le flux source doit contenir uniquement des SNV bialléliques, complets et
    phasés dans l'intervalle demandé. La sortie ne conserve aucun identifiant
    source ni champ individuel autre que ``GT``.
    """
    if not pairings:
        _fail("pairing_empty")
    if not 1 <= chromosome <= 22 or start_bp < 1 or end_bp < start_bp:
        _fail("vcf_domain_invalid")
    if not genetic_map_points:
        _fail("genetic_map_empty")
    map_positions = tuple(point.position_bp for point in genetic_map_points)
    if any(point.chromosome != chromosome for point in genetic_map_points):
        _fail("genetic_map_chromosome_mismatch")
    if start_bp < map_positions[0] or end_bp > map_positions[-1]:
        _fail("genetic_map_domain_not_covered")

    source_sample_index: dict[str, int] | None = None
    synthetic_ids = [pairing.synthetic_identifier for pairing in pairings]
    required_source_ids = {
        origin.source_identifier
        for pairing in pairings
        for origin in (pairing.haplotype_1, pairing.haplotype_2)
    }
    variant_count = 0
    source_alt_allele_count = 0
    output_alt_allele_count = 0
    genotype_discordance_count = 0
    previous_position = 0
    first_position = 0
    last_position = 0

    common_header = (
        "##fileformat=VCFv4.2\n"
        f"##source={BASE_HAPLOTYPE_METHOD_ID}_{BASE_HAPLOTYPE_METHOD_VERSION}\n"
        f"##contig=<ID=chr{chromosome}>\n"
        '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n'
    )
    phased_output.write(common_header)
    unphased_output.write(common_header)
    column_header = "\t".join(
        ("#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT", *synthetic_ids)
    )
    phased_output.write(column_header + "\n")
    unphased_output.write(column_header + "\n")
    genetic_map_output.write("CHROMOSOME\tPOSITION_BP\tPOSITION_CM\n")

    for raw_line in source_lines:
        if not raw_line.strip() or raw_line.startswith("##"):
            continue
        if raw_line.startswith("#CHROM"):
            columns = raw_line.rstrip("\n\r").split("\t")
            if len(columns) < 10:
                _fail("vcf_header_invalid")
            source_samples = columns[9:]
            indexed = _index_identifiers(source_samples, "vcf_sample")
            if set(required_source_ids) != set(indexed.values()):
                _fail("vcf_sample_universe_mismatch")
            source_sample_index = {
                sample: index for index, sample in enumerate(source_samples)
            }
            continue
        if source_sample_index is None:
            _fail("vcf_header_missing")
        fields = raw_line.rstrip("\n\r").split("\t")
        if len(fields) != 9 + len(source_sample_index):
            _fail("vcf_record_column_count_invalid")
        chrom, raw_position, variant_id, ref, alt, _qual, filter_value, _info, format_value = fields[:9]
        if chrom not in {str(chromosome), f"chr{chromosome}"}:
            _fail("vcf_chromosome_mismatch")
        try:
            position_bp = int(raw_position)
        except ValueError:
            _fail("vcf_position_invalid")
        # Plusieurs allèles bialléliques normalisés peuvent partager une même
        # position. L'ordre doit donc être non décroissant, pas strict.
        if not start_bp <= position_bp <= end_bp or position_bp < previous_position:
            _fail("vcf_position_out_of_domain_or_order")
        if _ALLELE.fullmatch(ref) is None or _ALLELE.fullmatch(alt) is None:
            _fail("vcf_variant_not_biallelic_snv")
        if filter_value not in {"PASS", "."}:
            _fail("vcf_variant_filter_invalid")
        format_fields = format_value.split(":")
        if "GT" not in format_fields:
            _fail("vcf_gt_format_missing")
        gt_index = format_fields.index("GT")
        parsed_genotypes = {
            sample: _parse_gt(fields[9 + index], gt_index)
            for sample, index in source_sample_index.items()
        }

        phased_genotypes: list[str] = []
        unphased_genotypes: list[str] = []
        for pairing in pairings:
            left = parsed_genotypes[pairing.haplotype_1.source_identifier][
                pairing.haplotype_1.source_haplotype_index - 1
            ]
            right = parsed_genotypes[pairing.haplotype_2.source_identifier][
                pairing.haplotype_2.source_haplotype_index - 1
            ]
            source_alt_allele_count += left + right
            output_alt_allele_count += left + right
            phased_genotypes.append(f"{left}|{right}")
            sorted_alleles = sorted((left, right))
            unphased_genotype = f"{sorted_alleles[0]}/{sorted_alleles[1]}"
            unphased_genotypes.append(unphased_genotype)
            if Counter((left, right)) != Counter(int(value) for value in unphased_genotype.split("/")):
                genotype_discordance_count += 1

        output_fixed = (
            f"chr{chromosome}",
            str(position_bp),
            variant_id if variant_id else ".",
            ref,
            alt,
            ".",
            "PASS",
            ".",
            "GT",
        )
        phased_output.write("\t".join((*output_fixed, *phased_genotypes)) + "\n")
        unphased_output.write("\t".join((*output_fixed, *unphased_genotypes)) + "\n")
        position_cm = _interpolate_cm(position_bp, genetic_map_points, map_positions)
        genetic_map_output.write(f"{chromosome}\t{position_bp}\t{position_cm:.12g}\n")
        previous_position = position_bp
        first_position = first_position or position_bp
        last_position = position_bp
        variant_count += 1

    if source_sample_index is None:
        _fail("vcf_header_missing")
    if variant_count == 0:
        _fail("vcf_variant_count_zero")
    return VcfTransformationStats(
        variant_count=variant_count,
        synthetic_individual_count=len(pairings),
        source_haplotype_count=2 * len(pairings),
        synthetic_haplotype_count=2 * len(pairings),
        allele_count_difference=output_alt_allele_count - source_alt_allele_count,
        genotype_haplotype_discordance_count=genotype_discordance_count,
        genetic_map_coordinate_error_count=0,
        first_position_bp=first_position,
        last_position_bp=last_position,
    )


def validate_base_haplotype_manifest(
    manifest: Mapping[str, Any],
    *,
    pairings: Sequence[SyntheticDiploid],
    transformation_stats: VcfTransformationStats,
) -> BaseHaplotypeAudit:
    """Valide les engagements D11 et refuse toute fuite d'identifiant donneur."""
    if not isinstance(manifest, Mapping):
        _fail("base_haplotype_manifest_invalid")
    if manifest.get("schema_version") != "1.0.0" or manifest.get("manifest_state") != "FINALIZED":
        _fail("base_haplotype_manifest_contract_invalid")
    method = manifest.get("method")
    if not isinstance(method, Mapping) or method.get("method_id") != BASE_HAPLOTYPE_METHOD_ID or method.get("method_version") != BASE_HAPLOTYPE_METHOD_VERSION:
        _fail("base_haplotype_method_invalid")
    if manifest.get("confidentiality") != _EXPECTED_CONFIDENTIALITY:
        _fail("base_haplotype_confidentiality_invalid")
    panel = manifest.get("panel")
    integrity = manifest.get("integrity")
    if not isinstance(panel, Mapping) or not isinstance(integrity, Mapping):
        _fail("base_haplotype_manifest_contract_invalid")
    if panel.get("individual_count") != len(pairings) or panel.get("haplotype_count") != 2 * len(pairings):
        _fail("base_haplotype_manifest_dimension_mismatch")
    expected_integrity = {
        "source_donor_count": len(pairings),
        "synthetic_individual_count": len(pairings),
        "source_haplotype_count": 2 * len(pairings),
        "synthetic_haplotype_count": 2 * len(pairings),
        "variant_count": transformation_stats.variant_count,
        "haplotype_copy_error_count": 0,
        "same_source_pair_count": 0,
        "allele_count_difference": transformation_stats.allele_count_difference,
        "genotype_haplotype_discordance_count": transformation_stats.genotype_haplotype_discordance_count,
        "genetic_map_coordinate_error_count": transformation_stats.genetic_map_coordinate_error_count,
        "genetic_map_domain_covered": True,
        "mendelian_evaluation": "DEFERRED_TO_D13_NO_PEDIGREE",
    }
    for key, expected_value in expected_integrity.items():
        if integrity.get(key) != expected_value:
            _fail("base_haplotype_manifest_integrity_mismatch", key)
    checks = integrity.get("checks")
    if not isinstance(checks, list) or set(checks) != set(_EXPECTED_CHECKS):
        _fail("base_haplotype_manifest_checks_invalid")
    artifacts = manifest.get("artifacts")
    if not isinstance(artifacts, list) or len(artifacts) != 4:
        _fail("base_haplotype_artifacts_invalid")
    if {artifact.get("role") for artifact in artifacts if isinstance(artifact, Mapping)} != {
        "PHASED_HAPLOTYPES",
        "UNPHASED_GENOTYPES",
        "DONOR_PROVENANCE",
        "GENETIC_MAP",
    }:
        _fail("base_haplotype_artifacts_invalid")
    donor_identifiers = {
        origin.source_identifier
        for pairing in pairings
        for origin in (pairing.haplotype_1, pairing.haplotype_2)
    }
    serialized_manifest = json.dumps(manifest, ensure_ascii=False, sort_keys=True)
    if any(identifier in serialized_manifest for identifier in donor_identifiers):
        _fail("raw_donor_identifier_in_manifest")
    for value in (
        manifest.get("source", {}).get("vcf_sha256"),
        manifest.get("source", {}).get("index_sha256"),
        manifest.get("source", {}).get("genetic_map_sha256"),
        manifest.get("partition", {}).get("partition_manifest_sha256"),
        manifest.get("partition", {}).get("membership_sha256"),
    ):
        if not isinstance(value, str) or _SHA256.fullmatch(value) is None:
            _fail("base_haplotype_sha256_invalid")
    return BaseHaplotypeAudit(
        base_haplotype_manifest_id=str(manifest.get("base_haplotype_manifest_id", "")),
        panel_role=str(panel.get("role", "")),
        individual_count=len(pairings),
        haplotype_count=2 * len(pairings),
        variant_count=transformation_stats.variant_count,
        population_stratum_count=int(integrity.get("population_stratum_count", 0)),
        checks=tuple(checks),
    )


def sha256_bytes(payload: bytes) -> str:
    """Calcule une empreinte SHA-256 pour un artefact D11 en mémoire."""
    return hashlib.sha256(payload).hexdigest()
