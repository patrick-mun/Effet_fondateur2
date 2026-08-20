"""Harmonisation allélique et préparation des PCA globale et haplotypique locale."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Sequence

import numpy as np

from effet_fondateur.ancestry.pca import ReferencePca, fit_reference_pca, project_pca


class AncestryAnalysisError(ValueError):
    """Signale une identité de variant ou une matrice d'ascendance incohérente."""


@dataclass(frozen=True)
class Variant:
    """Identité canonique GRCh38 d'un variant biallélique."""

    chromosome: str
    position_bp: int
    variant_id: str
    ref: str
    alt: str

    @property
    def locus(self) -> tuple[str, int]:
        return (self.chromosome.removeprefix("chr"), self.position_bp)


@dataclass(frozen=True)
class GenotypePanel:
    """Matrice échantillons × variants comptant l'allèle ALT déclaré."""

    sample_ids: tuple[str, ...]
    variants: tuple[Variant, ...]
    alt_dosages: np.ndarray
    ploidy: int


@dataclass(frozen=True)
class HarmonizedPca:
    """PCA ajustée sur la référence et projection indépendante de l'étude."""

    model: ReferencePca
    variants: tuple[Variant, ...]
    reference_sample_ids: tuple[str, ...]
    study_sample_ids: tuple[str, ...]
    study_scores: np.ndarray
    candidate_variant_count: int


def _validate_panel(
    panel: GenotypePanel, name: str, *, allow_duplicate_loci: bool = False
) -> None:
    matrix = np.asarray(panel.alt_dosages, dtype=float)
    if panel.ploidy not in {1, 2}:
        raise AncestryAnalysisError(f"invalid_panel_ploidy:{name}")
    if (
        not panel.sample_ids
        or len(panel.sample_ids) != len(set(panel.sample_ids))
        or not panel.variants
        or matrix.shape != (len(panel.sample_ids), len(panel.variants))
        or np.isinf(matrix).any()
        or np.any((matrix[~np.isnan(matrix)] < 0))
        or np.any((matrix[~np.isnan(matrix)] > panel.ploidy))
    ):
        raise AncestryAnalysisError(f"invalid_genotype_panel:{name}")
    loci = [variant.locus for variant in panel.variants]
    if not allow_duplicate_loci and len(loci) != len(set(loci)):
        raise AncestryAnalysisError(f"duplicate_genotype_variant:{name}")


def harmonize_alt_dosages(
    reference: GenotypePanel,
    study: GenotypePanel,
    *,
    requested_components: int,
    minimum_variants: int,
    minimum_reference_call_rate: float = 0.0,
) -> HarmonizedPca:
    """Ajuste une PCA sur la référence après correspondance coordonnée/allèles.

    Les inversions REF/ALT de l'étude sont corrigées par ``ploidy - dosage``.
    Les compléments de brin et les correspondances par identifiant seul sont
    volontairement interdits afin d'éviter une harmonisation ambiguë.
    """

    # Un VCF public peut représenter plusieurs variants bialléliques à la même
    # position. Ils ne sont acceptés que si les allèles de l'étude désignent
    # ensuite un enregistrement de référence unique.
    _validate_panel(reference, "reference", allow_duplicate_loci=True)
    _validate_panel(study, "study")
    if reference.ploidy != study.ploidy:
        raise AncestryAnalysisError("ancestry_panel_ploidy_mismatch")
    if requested_components < 1 or minimum_variants < 1 or not 0 <= minimum_reference_call_rate <= 1:
        raise AncestryAnalysisError("invalid_ancestry_analysis_parameters")

    reference_by_locus: dict[tuple[str, int], list[int]] = {}
    for index, variant in enumerate(reference.variants):
        reference_by_locus.setdefault(variant.locus, []).append(index)
    reference_indexes: list[int] = []
    study_indexes: list[int] = []
    reverse_study: list[bool] = []
    selected_variants: list[Variant] = []
    for study_index, study_variant in enumerate(study.variants):
        locus_indexes = reference_by_locus.get(study_variant.locus, [])
        compatible_indexes = [
            index
            for index in locus_indexes
            if (study_variant.ref, study_variant.alt)
            in (
                (reference.variants[index].ref, reference.variants[index].alt),
                (reference.variants[index].alt, reference.variants[index].ref),
            )
        ]
        if not compatible_indexes:
            continue
        if len(compatible_indexes) > 1:
            raise AncestryAnalysisError("ambiguous_reference_variant")
        reference_index = compatible_indexes[0]
        reference_variant = reference.variants[reference_index]
        if (study_variant.ref, study_variant.alt) == (reference_variant.ref, reference_variant.alt):
            reverse = False
        elif (study_variant.ref, study_variant.alt) == (reference_variant.alt, reference_variant.ref):
            reverse = True
        else:
            continue
        reference_values = np.asarray(reference.alt_dosages, dtype=float)[:, reference_index]
        if float(np.mean(~np.isnan(reference_values))) < minimum_reference_call_rate:
            continue
        reference_indexes.append(reference_index)
        study_indexes.append(study_index)
        reverse_study.append(reverse)
        selected_variants.append(reference_variant)
    if len(selected_variants) < minimum_variants:
        raise AncestryAnalysisError("insufficient_harmonized_ancestry_variants")

    reference_matrix = np.asarray(reference.alt_dosages, dtype=float)[:, reference_indexes]
    study_matrix = np.asarray(study.alt_dosages, dtype=float)[:, study_indexes].copy()
    for column, reverse in enumerate(reverse_study):
        if reverse:
            observed = ~np.isnan(study_matrix[:, column])
            study_matrix[observed, column] = study.ploidy - study_matrix[observed, column]
    model = fit_reference_pca(
        reference_matrix,
        requested_components=requested_components,
        ploidy=reference.ploidy,
    )
    scores = project_pca(model, study_matrix, ploidy=study.ploidy)
    informative_variants = tuple(
        variant
        for variant, informative in zip(selected_variants, model.informative_variant_mask)
        if informative
    )
    return HarmonizedPca(
        model=model,
        variants=informative_variants,
        reference_sample_ids=reference.sample_ids,
        study_sample_ids=study.sample_ids,
        study_scores=scores,
        candidate_variant_count=len(study.variants),
    )


def genotype_to_alt_dosage(genotype: str, *, phased_required: bool) -> float:
    """Convertit un GT diploïde biallélique en dosage ALT, ou NaN s'il manque."""

    if phased_required and "|" not in genotype:
        raise AncestryAnalysisError("local_reference_genotype_not_phased")
    separator = "|" if "|" in genotype else "/"
    alleles = genotype.split(separator)
    if len(alleles) != 2 or any(allele not in {"0", "1", "."} for allele in alleles):
        raise AncestryAnalysisError("unsupported_ancestry_genotype")
    if "." in alleles:
        return float("nan")
    return float(sum(allele == "1" for allele in alleles))


def phased_genotype_to_haplotypes(genotype: str) -> tuple[float, float]:
    """Décompose un GT phasé en deux dosages haploïdes ALT."""

    if "|" not in genotype:
        raise AncestryAnalysisError("local_study_genotype_not_phased")
    alleles = genotype.split("|")
    if len(alleles) != 2 or any(allele not in {"0", "1", "."} for allele in alleles):
        raise AncestryAnalysisError("unsupported_ancestry_genotype")
    return tuple(float("nan") if allele == "." else float(allele) for allele in alleles)  # type: ignore[return-value]


def population_centroids(
    scores: np.ndarray,
    labels: Sequence[str],
) -> dict[str, tuple[int, np.ndarray]]:
    """Calcule des centroïdes descriptifs, sans transformer ceux-ci en attribution."""

    matrix = np.asarray(scores, dtype=float)
    if matrix.ndim != 2 or matrix.shape[0] != len(labels) or not labels:
        raise AncestryAnalysisError("invalid_population_centroid_input")
    result: dict[str, tuple[int, np.ndarray]] = {}
    for label in sorted(set(labels)):
        indexes = [index for index, value in enumerate(labels) if value == label]
        result[label] = (len(indexes), np.mean(matrix[indexes], axis=0))
    return result
