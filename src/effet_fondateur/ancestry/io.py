"""Lecteurs stricts des génotypes nécessaires à l'étape d'ascendance."""

from __future__ import annotations

import csv
from pathlib import Path
from typing import Iterable, Sequence

import numpy as np

from effet_fondateur.ancestry.analysis import (
    AncestryAnalysisError,
    GenotypePanel,
    Variant,
    genotype_to_alt_dosage,
    phased_genotype_to_haplotypes,
)


def read_bim_variants(path: Path) -> tuple[Variant, ...]:
    """Lit un BIM en comptant A1 comme ALT, conformément au dosage PLINK ``A``."""

    variants: list[Variant] = []
    for line in path.read_text(encoding="utf-8").splitlines():
        if not line:
            continue
        fields = line.split()
        if len(fields) != 6 or fields[4] not in {"A", "C", "G", "T"} or fields[5] not in {"A", "C", "G", "T"}:
            raise AncestryAnalysisError("invalid_ancestry_bim_variant")
        try:
            position = int(fields[3])
        except ValueError as error:
            raise AncestryAnalysisError("invalid_ancestry_bim_position") from error
        variants.append(Variant(f"chr{fields[0].removeprefix('chr')}", position, fields[1], fields[5], fields[4]))
    if not variants:
        raise AncestryAnalysisError("empty_ancestry_bim")
    if (
        len({variant.locus for variant in variants}) != len(variants)
        or len({variant.variant_id for variant in variants}) != len(variants)
    ):
        raise AncestryAnalysisError("duplicate_ancestry_bim_variant")
    return tuple(variants)


def _dosage_column_indexes(
    header: Sequence[str], variants: tuple[Variant, ...],
) -> list[int]:
    """Résout les colonnes de dosage en temps linéaire dans la taille du panel."""

    column_by_name: dict[str, int] = {}
    duplicate_names: set[str] = set()
    for index, name in enumerate(header):
        if name in column_by_name:
            duplicate_names.add(name)
        else:
            column_by_name[name] = index

    dosage_indexes: list[int] = []
    for variant in variants:
        expected = f"{variant.variant_id}_{variant.alt}"
        if expected not in column_by_name or expected in duplicate_names:
            raise AncestryAnalysisError(
                "plink_ancestry_dosage_column_missing_or_ambiguous"
            )
        dosage_indexes.append(column_by_name[expected])
    return dosage_indexes


def read_plink_raw_panel(
    path: Path,
    variants: tuple[Variant, ...],
    sample_by_plink_id: dict[tuple[str, str], str],
) -> GenotypePanel:
    """Lit ``plink --recode A`` et restitue les dosages A1 dans l'ordre du BIM."""

    with path.open(encoding="utf-8", newline="") as handle:
        reader = csv.reader(handle, delimiter=" ", skipinitialspace=True)
        try:
            header = next(reader)
        except StopIteration as error:
            raise AncestryAnalysisError("empty_plink_ancestry_raw") from error
        if header[:6] != ["FID", "IID", "PAT", "MAT", "SEX", "PHENOTYPE"]:
            raise AncestryAnalysisError("invalid_plink_ancestry_raw_header")
        dosage_indexes = _dosage_column_indexes(header, variants)
        sample_ids: list[str] = []
        matrix: list[list[float]] = []
        for fields in reader:
            if not fields:
                continue
            if len(fields) != len(header):
                raise AncestryAnalysisError("invalid_plink_ancestry_raw_row")
            sample_id = sample_by_plink_id.get((fields[0], fields[1]))
            if sample_id is None:
                raise AncestryAnalysisError("plink_ancestry_sample_absent_from_registry")
            sample_ids.append(sample_id)
            row: list[float] = []
            for index in dosage_indexes:
                if fields[index] == "NA":
                    row.append(float("nan"))
                else:
                    try:
                        value = float(fields[index])
                    except ValueError as error:
                        raise AncestryAnalysisError("invalid_plink_ancestry_dosage") from error
                    if value not in {0.0, 1.0, 2.0}:
                        raise AncestryAnalysisError("invalid_plink_ancestry_dosage")
                    row.append(value)
            matrix.append(row)
    return GenotypePanel(tuple(sample_ids), variants, np.asarray(matrix, dtype=float), 2)


def parse_vcf_query_panel(
    lines: Iterable[str],
    sample_ids: tuple[str, ...],
    *,
    haplotypes: bool,
    phased_required: bool,
) -> GenotypePanel:
    """Transforme la sortie contrôlée de ``bcftools query`` en matrice ALT."""

    variants: list[Variant] = []
    dosage_columns: list[list[float]] = []
    for line in lines:
        if not line.strip():
            continue
        fields = line.rstrip("\n").split("\t")
        if len(fields) != 5 + len(sample_ids):
            raise AncestryAnalysisError("invalid_ancestry_vcf_query_record")
        try:
            position = int(fields[1])
        except ValueError as error:
            raise AncestryAnalysisError("invalid_ancestry_vcf_query_position") from error
        if "," in fields[4] or not fields[3] or not fields[4]:
            raise AncestryAnalysisError("ancestry_vcf_variant_not_biallelic")
        variants.append(Variant(fields[0], position, fields[2], fields[3], fields[4]))
        if haplotypes:
            haploid: list[float] = []
            for genotype in fields[5:]:
                haploid.extend(phased_genotype_to_haplotypes(genotype))
            dosage_columns.append(haploid)
        else:
            dosage_columns.append([
                genotype_to_alt_dosage(genotype, phased_required=phased_required)
                for genotype in fields[5:]
            ])
    if haplotypes:
        output_ids = tuple(f"{sample}:{haplotype}" for sample in sample_ids for haplotype in ("H1", "H2"))
    else:
        output_ids = sample_ids
    matrix = np.asarray(dosage_columns, dtype=float).T if dosage_columns else np.empty((len(output_ids), 0))
    return GenotypePanel(output_ids, tuple(variants), matrix, 1 if haplotypes else 2)
