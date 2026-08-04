"""Annotation de positions hg38 avec les identifiants RefSNP de NCBI."""

from __future__ import annotations

import re
import shutil
import subprocess
import tempfile
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable


DEFAULT_DBSNP_GRCH38_VCF = (
    "https://ftp.ncbi.nih.gov/snp/archive/b157/VCF/GCF_000001405.40.gz"
)

GRCH38_PRIMARY_ACCESSIONS = {
    "1": "NC_000001.11",
    "2": "NC_000002.12",
    "3": "NC_000003.12",
    "4": "NC_000004.12",
    "5": "NC_000005.10",
    "6": "NC_000006.12",
    "7": "NC_000007.14",
    "8": "NC_000008.11",
    "9": "NC_000009.12",
    "10": "NC_000010.11",
    "11": "NC_000011.10",
    "12": "NC_000012.12",
    "13": "NC_000013.11",
    "14": "NC_000014.9",
    "15": "NC_000015.10",
    "16": "NC_000016.10",
    "17": "NC_000017.11",
    "18": "NC_000018.10",
    "19": "NC_000019.10",
    "20": "NC_000020.11",
    "21": "NC_000021.9",
    "22": "NC_000022.11",
    "X": "NC_000023.11",
    "Y": "NC_000024.10",
    "MT": "NC_012920.1",
}

RSID_PATTERN = re.compile(r"^rs[0-9]+$")


class DbsnpAnnotationError(RuntimeError):
    """Erreur empêchant une annotation dbSNP reproductible."""


@dataclass(frozen=True)
class DbsnpVariant:
    """Variant SNV dbSNP placé sur l'assemblage GRCh38."""

    position_bp: int
    rsid: str
    reference: str
    alternates: tuple[str, ...]

    @property
    def alleles(self) -> frozenset[str]:
        return frozenset((self.reference, *self.alternates))


def normalize_dbsnp_chromosome(chromosome: str) -> str:
    """Normalise un chromosome pour retrouver son accession RefSeq GRCh38."""
    normalized = chromosome.strip().upper()
    if normalized.startswith("CHR"):
        normalized = normalized[3:]
    if normalized == "M":
        normalized = "MT"
    return normalized


def extract_rsids(identifier_field: str) -> tuple[str, ...]:
    """Extrait les identifiants rs valides de la colonne ID d'un VCF."""
    return tuple(
        identifier
        for identifier in identifier_field.split(";")
        if RSID_PATTERN.fullmatch(identifier)
    )


def parse_dbsnp_vcf_records(
    vcf_text: str, requested_positions: set[int]
) -> dict[int, list[DbsnpVariant]]:
    """Conserve les SNV dbSNP commençant exactement aux positions demandées."""
    variants_by_position: dict[int, list[DbsnpVariant]] = {}
    for line in vcf_text.splitlines():
        if not line or line.startswith("#"):
            continue
        columns = line.split("\t")
        if len(columns) < 5:
            continue
        try:
            position_bp = int(columns[1])
        except ValueError:
            continue
        if position_bp not in requested_positions:
            continue

        reference = columns[3].upper()
        alternates = tuple(allele.upper() for allele in columns[4].split(","))
        if len(reference) != 1 or not alternates:
            continue
        if any(len(allele) != 1 or allele == "*" for allele in alternates):
            continue

        for rsid in extract_rsids(columns[2]):
            variants_by_position.setdefault(position_bp, []).append(
                DbsnpVariant(
                    position_bp=position_bp,
                    rsid=rsid,
                    reference=reference,
                    alternates=alternates,
                )
            )
    return variants_by_position


def query_dbsnp_variants(
    positions: Iterable[int],
    chromosome: str,
    vcf_source: str = DEFAULT_DBSNP_GRCH38_VCF,
    bcftools_executable: str = "bcftools",
) -> dict[int, list[DbsnpVariant]]:
    """Interroge un VCF dbSNP indexé uniquement aux positions sélectionnées."""
    requested_positions = set(positions)
    if not requested_positions:
        return {}

    normalized_chromosome = normalize_dbsnp_chromosome(chromosome)
    accession = GRCH38_PRIMARY_ACCESSIONS.get(normalized_chromosome)
    if accession is None:
        raise DbsnpAnnotationError(
            f"Chromosome non pris en charge pour GRCh38: {chromosome}"
        )
    executable = shutil.which(bcftools_executable)
    if executable is None:
        raise DbsnpAnnotationError(
            f"Exécutable bcftools introuvable: {bcftools_executable}"
        )

    with tempfile.TemporaryDirectory(prefix="acpa_dbsnp_") as temporary_directory:
        regions_path = Path(temporary_directory) / "positions.tsv"
        regions_path.write_text(
            "".join(
                f"{accession}\t{position}\t{position}\n"
                for position in sorted(requested_positions)
            ),
            encoding="utf-8",
        )
        command = [
            executable,
            "view",
            "--no-version",
            "-H",
            "-R",
            str(regions_path),
            vcf_source,
        ]
        completed_process = subprocess.run(
            command,
            check=False,
            capture_output=True,
            text=True,
            cwd=temporary_directory,
        )
        if completed_process.returncode != 0:
            details = completed_process.stderr.strip() or "erreur bcftools non détaillée"
            raise DbsnpAnnotationError(f"Échec de l'interrogation dbSNP: {details}")

    return parse_dbsnp_vcf_records(
        completed_process.stdout,
        requested_positions,
    )


def select_compatible_rsid(
    variants: Iterable[DbsnpVariant],
    observed_alleles: frozenset[str],
    embedded_rsids: frozenset[str] = frozenset(),
) -> tuple[str, str]:
    """Choisit un rsID unique compatible avec les allèles forward-strand."""
    variants = tuple(variants)
    compatible = tuple(
        variant
        for variant in variants
        if not observed_alleles or observed_alleles.issubset(variant.alleles)
    )
    if len(observed_alleles) == 2:
        exact = tuple(
            variant for variant in compatible if variant.alleles == observed_alleles
        )
        exact_rsids = {variant.rsid for variant in exact}
        if len(exact_rsids) == 1:
            return next(iter(exact_rsids)), "dbsnp_exact"

    compatible_rsids = {variant.rsid for variant in compatible}
    confirmed_embedded = compatible_rsids.intersection(embedded_rsids)
    if len(confirmed_embedded) == 1:
        return next(iter(confirmed_embedded)), "dbsnp_embedded_confirmed"
    if len(compatible_rsids) == 1:
        return next(iter(compatible_rsids)), "dbsnp_compatible"
    if len(compatible_rsids) > 1:
        return "", "dbsnp_ambiguous"
    if variants:
        return "", "dbsnp_alleles_incompatible"
    return "", "dbsnp_not_found"
