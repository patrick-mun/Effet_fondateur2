"""Harmonisation D11b des marqueurs CytoScan avec la vérité dense 1000G.

La vérité haplotypique D11a reste dense et orientée selon la référence 1000G.
Ce module construit séparément l'univers observable par la puce à partir d'un
inventaire de marqueurs indépendant des échantillons. Un marqueur n'est retenu
que si son assemblage, sa coordonnée et sa paire d'allèles sont compatibles
sans ambiguïté avec un SNV biallélique de la vérité dense.
"""

from __future__ import annotations

from collections import Counter, defaultdict
from collections.abc import Iterable
import csv
from dataclasses import dataclass
import re
from typing import TextIO


CHIP_DESIGN_ID = "CytoScan750K_Accel"
CHIP_DESIGN_VERSION = "NA36.r1"
CHIP_ASSEMBLY = "GRCh38"
_BASES = frozenset("ACGT")
_COMPLEMENT = {"A": "T", "C": "G", "G": "C", "T": "A"}
_STABLE_IDENTIFIER = re.compile(r"^[A-Za-z0-9][A-Za-z0-9._:-]*$")
_REQUIRED_COLUMNS = (
    "MARKER_ID",
    "ARRAY_DESIGN",
    "DESIGN_VERSION",
    "ASSEMBLY",
    "CHROMOSOME",
    "POSITION_BP",
    "ALLELE_A",
    "ALLELE_B",
    "STRAND",
)


class ChipMarkerError(RuntimeError):
    """Signale une violation D11b avec un code stable et non sensible."""

    def __init__(self, code: str, detail: str | None = None) -> None:
        self.code = code
        self.detail = detail
        super().__init__(code if detail is None else f"{code}: {detail}")


@dataclass(frozen=True, order=True)
class ChipMarker:
    """Marqueur SNP indépendant de tout échantillon, en coordonnées GRCh38."""

    marker_id: str
    chromosome: int
    position_bp: int
    allele_a: str
    allele_b: str
    strand: str


@dataclass(frozen=True, order=True)
class ReferenceVariant:
    """SNV biallélique de la vérité dense, orienté comme le VCF 1000G."""

    chromosome: int
    position_bp: int
    variant_id: str
    reference_allele: str
    alternate_allele: str


@dataclass(frozen=True, order=True)
class MarkerMatch:
    """Correspondance non ambiguë entre une sonde SNP et un variant dense."""

    marker: ChipMarker
    reference_variant: ReferenceVariant
    allele_relation: str


@dataclass(frozen=True)
class MarkerHarmonizationAudit:
    """Comptes agrégés D11b, sans génotype ni identifiant d'individu."""

    chip_marker_count: int
    reference_variant_count: int
    matched_marker_count: int
    excluded_marker_count: int
    exclusion_counts: tuple[tuple[str, int], ...]
    relation_counts: tuple[tuple[str, int], ...]


def _fail(code: str, detail: str | None = None) -> None:
    raise ChipMarkerError(code, detail)


def _parse_chromosome(value: str) -> int:
    normalized = value.removeprefix("chr").removeprefix("CHR")
    try:
        chromosome = int(normalized)
    except ValueError:
        _fail("chip_marker_chromosome_invalid")
    if not 1 <= chromosome <= 22:
        _fail("chip_marker_chromosome_invalid")
    return chromosome


def parse_chip_marker_manifest(lines: Iterable[str]) -> tuple[ChipMarker, ...]:
    """Lit l'export SNP normalisé issu du paquet officiel NA36.r1.

    Le fichier doit être une table TSV sans colonne d'échantillon. Les champs
    de conception, version et assemblage sont répétés sur chaque ligne afin
    qu'une concaténation ou une mauvaise annotation ne passe pas silencieusement.
    """
    reader = csv.DictReader(lines, delimiter="\t")
    if tuple(reader.fieldnames or ()) != _REQUIRED_COLUMNS:
        _fail("chip_marker_columns_invalid")
    markers: list[ChipMarker] = []
    identifiers: set[str] = set()
    for row in reader:
        if row["ARRAY_DESIGN"] != CHIP_DESIGN_ID:
            _fail("chip_marker_design_mismatch")
        if row["DESIGN_VERSION"] != CHIP_DESIGN_VERSION:
            _fail("chip_marker_design_version_mismatch")
        if row["ASSEMBLY"] != CHIP_ASSEMBLY:
            _fail("chip_marker_assembly_mismatch")
        marker_id = row["MARKER_ID"]
        if _STABLE_IDENTIFIER.fullmatch(marker_id) is None:
            _fail("chip_marker_identifier_invalid")
        collision_key = marker_id.casefold()
        if collision_key in identifiers:
            _fail("chip_marker_identifier_duplicate")
        identifiers.add(collision_key)
        try:
            position_bp = int(row["POSITION_BP"])
        except ValueError:
            _fail("chip_marker_position_invalid")
        if position_bp < 1:
            _fail("chip_marker_position_invalid")
        allele_a = row["ALLELE_A"].upper()
        allele_b = row["ALLELE_B"].upper()
        if allele_a not in _BASES or allele_b not in _BASES or allele_a == allele_b:
            _fail("chip_marker_alleles_invalid")
        strand = row["STRAND"].upper()
        if strand not in {"+", "-", "UNKNOWN"}:
            _fail("chip_marker_strand_invalid")
        markers.append(
            ChipMarker(
                marker_id=marker_id,
                chromosome=_parse_chromosome(row["CHROMOSOME"]),
                position_bp=position_bp,
                allele_a=allele_a,
                allele_b=allele_b,
                strand=strand,
            )
        )
    if not markers:
        _fail("chip_marker_manifest_empty")
    return tuple(sorted(markers))


def parse_reference_vcf_variants(
    lines: Iterable[str],
    *,
    chromosome: int,
    start_bp: int,
    end_bp: int,
) -> tuple[ReferenceVariant, ...]:
    """Extrait les SNV bialléliques du VCF dense dans le domaine demandé."""
    if not 1 <= chromosome <= 22 or start_bp < 1 or end_bp < start_bp:
        _fail("chip_reference_domain_invalid")
    variants: list[ReferenceVariant] = []
    header_seen = False
    for line in lines:
        if line.startswith("##"):
            continue
        if line.startswith("#CHROM\t"):
            header_seen = True
            continue
        if line.startswith("#") or not line.strip():
            continue
        if not header_seen:
            _fail("chip_reference_vcf_header_missing")
        fields = line.rstrip("\n").split("\t")
        if len(fields) < 8:
            _fail("chip_reference_vcf_row_invalid")
        row_chromosome = _parse_chromosome(fields[0])
        try:
            position_bp = int(fields[1])
        except ValueError:
            _fail("chip_reference_vcf_position_invalid")
        if row_chromosome != chromosome or not start_bp <= position_bp <= end_bp:
            continue
        reference_allele = fields[3].upper()
        alternate_allele = fields[4].upper()
        if (
            reference_allele not in _BASES
            or alternate_allele not in _BASES
            or reference_allele == alternate_allele
        ):
            continue
        variants.append(
            ReferenceVariant(
                chromosome=row_chromosome,
                position_bp=position_bp,
                variant_id=fields[2],
                reference_allele=reference_allele,
                alternate_allele=alternate_allele,
            )
        )
    if not header_seen:
        _fail("chip_reference_vcf_header_missing")
    if not variants:
        _fail("chip_reference_variant_empty")
    return tuple(variants)


def _is_palindromic(marker: ChipMarker) -> bool:
    return {marker.allele_a, marker.allele_b} in ({"A", "T"}, {"C", "G"})


def _allele_relation(
    marker: ChipMarker,
    variant: ReferenceVariant,
) -> str | None:
    reference_pair = (variant.reference_allele, variant.alternate_allele)
    direct_pair = (marker.allele_a, marker.allele_b)
    complement_pair = tuple(_COMPLEMENT[allele] for allele in direct_pair)
    if marker.strand == "UNKNOWN" and _is_palindromic(marker):
        return None
    allowed_pairs: list[tuple[tuple[str, str], str]] = []
    if marker.strand in {"+", "UNKNOWN"}:
        allowed_pairs.extend(
            ((direct_pair, "DIRECT"), (direct_pair[::-1], "SWAPPED"))
        )
    if marker.strand in {"-", "UNKNOWN"}:
        allowed_pairs.extend(
            (
                (complement_pair, "COMPLEMENT"),
                (complement_pair[::-1], "COMPLEMENT_SWAPPED"),
            )
        )
    relations = {relation for pair, relation in allowed_pairs if pair == reference_pair}
    if len(relations) != 1:
        return None
    return relations.pop()


def harmonize_chip_markers(
    markers: Iterable[ChipMarker],
    variants: Iterable[ReferenceVariant],
    *,
    chromosome: int,
    start_bp: int,
    end_bp: int,
) -> tuple[tuple[MarkerMatch, ...], MarkerHarmonizationAudit]:
    """Construit l'intersection puce–1000G avec exclusions conservatrices."""
    marker_list = tuple(markers)
    variant_list = tuple(variants)
    variants_by_position: dict[int, list[ReferenceVariant]] = defaultdict(list)
    for variant in variant_list:
        if variant.chromosome == chromosome and start_bp <= variant.position_bp <= end_bp:
            variants_by_position[variant.position_bp].append(variant)
    markers_in_domain = tuple(
        marker
        for marker in marker_list
        if marker.chromosome == chromosome and start_bp <= marker.position_bp <= end_bp
    )
    marker_position_counts = Counter(marker.position_bp for marker in markers_in_domain)
    exclusions: Counter[str] = Counter()
    matches: list[MarkerMatch] = []
    for marker in markers_in_domain:
        if marker_position_counts[marker.position_bp] > 1:
            exclusions["DUPLICATE_CHIP_POSITION"] += 1
            continue
        candidates = variants_by_position.get(marker.position_bp, [])
        if not candidates:
            exclusions["NO_REFERENCE_VARIANT"] += 1
            continue
        if len(candidates) != 1:
            exclusions["REFERENCE_POSITION_AMBIGUOUS"] += 1
            continue
        if marker.strand == "UNKNOWN" and _is_palindromic(marker):
            exclusions["PALINDROMIC_STRAND_AMBIGUOUS"] += 1
            continue
        relation = _allele_relation(marker, candidates[0])
        if relation is None:
            exclusions["ALLELES_INCOMPATIBLE_OR_AMBIGUOUS"] += 1
            continue
        matches.append(MarkerMatch(marker, candidates[0], relation))
    relation_counts = Counter(match.allele_relation for match in matches)
    audit = MarkerHarmonizationAudit(
        chip_marker_count=len(markers_in_domain),
        reference_variant_count=sum(len(value) for value in variants_by_position.values()),
        matched_marker_count=len(matches),
        excluded_marker_count=sum(exclusions.values()),
        exclusion_counts=tuple(sorted(exclusions.items())),
        relation_counts=tuple(sorted(relation_counts.items())),
    )
    if audit.matched_marker_count + audit.excluded_marker_count != audit.chip_marker_count:
        _fail("chip_marker_harmonization_accounting_error")
    return tuple(sorted(matches)), audit


def write_marker_map(matches: Iterable[MarkerMatch], output: TextIO) -> int:
    """Écrit la table de liaison puce–1000G sans donnée individuelle."""
    writer = csv.writer(output, delimiter="\t", lineterminator="\n")
    writer.writerow(
        (
            "MARKER_ID",
            "CHROMOSOME",
            "POSITION_BP",
            "REFERENCE_ALLELE",
            "ALTERNATE_ALLELE",
            "ALLELE_RELATION",
        )
    )
    count = 0
    for match in sorted(matches):
        writer.writerow(
            (
                match.marker.marker_id,
                match.marker.chromosome,
                match.marker.position_bp,
                match.reference_variant.reference_allele,
                match.reference_variant.alternate_allele,
                match.allele_relation,
            )
        )
        count += 1
    return count


def filter_vcf_to_chip_markers(
    lines: Iterable[str],
    output: TextIO,
    matches: Iterable[MarkerMatch],
) -> int:
    """Filtre un VCF synthétique dense et remplace l'ID par celui de la puce.

    Les allèles et génotypes restent orientés comme 1000G. La relation de sonde
    est conservée dans la table de liaison et sera utilisée lors de l'export
    ACPA ; aucune réorientation silencieuse des génotypes n'est faite ici.
    """
    matches_by_key = {
        (match.marker.chromosome, match.marker.position_bp): match
        for match in matches
    }
    written = 0
    for line in lines:
        if line.startswith("#"):
            output.write(line)
            continue
        if not line.strip():
            continue
        fields = line.rstrip("\n").split("\t")
        if len(fields) < 8:
            _fail("chip_observed_vcf_row_invalid")
        try:
            key = (_parse_chromosome(fields[0]), int(fields[1]))
        except ValueError:
            _fail("chip_observed_vcf_position_invalid")
        match = matches_by_key.get(key)
        if match is None:
            continue
        if (
            fields[3].upper() != match.reference_variant.reference_allele
            or fields[4].upper() != match.reference_variant.alternate_allele
        ):
            _fail("chip_observed_vcf_allele_mismatch")
        fields[2] = match.marker.marker_id
        output.write("\t".join(fields) + "\n")
        written += 1
    if written != len(matches_by_key):
        _fail("chip_observed_vcf_marker_count_mismatch")
    return written
