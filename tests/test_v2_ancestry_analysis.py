import numpy as np
import pytest

from effet_fondateur.ancestry import (
    AncestryAnalysisError,
    GenotypePanel,
    Variant,
    genotype_to_alt_dosage,
    harmonize_alt_dosages,
    phased_genotype_to_haplotypes,
    population_centroids,
)


def _variants() -> tuple[Variant, ...]:
    return tuple(
        Variant("chr2", 100 + index, f"v{index}", "A", "G")
        for index in range(4)
    )


def test_global_harmonization_fits_reference_and_projects_swapped_study() -> None:
    variants = _variants()
    reference = GenotypePanel(
        ("R1", "R2", "R3", "R4"),
        variants,
        np.array([[0, 0, 0, 1], [0, 1, 0, 0], [2, 2, 2, 1], [2, 1, 2, 2]], dtype=float),
        2,
    )
    swapped = tuple(
        Variant(value.chromosome, value.position_bp, value.variant_id, value.alt, value.ref)
        for value in variants
    )
    study = GenotypePanel(
        ("S1", "S2"),
        swapped,
        np.array([[2, 2, 2, 1], [0, 0, 0, 0]], dtype=float),
        2,
    )

    result = harmonize_alt_dosages(
        reference, study, requested_components=2, minimum_variants=4
    )

    assert result.study_scores.shape == (2, 2)
    assert result.model.reference_scores.shape == (4, 2)
    assert result.study_scores[0, 0] * result.study_scores[1, 0] < 0


def test_local_haplotype_analysis_is_generic_for_configured_region() -> None:
    variants = tuple(
        Variant("chr7", 55_000_000 + index, f"region_variant_{index}", "C", "T")
        for index in range(3)
    )
    reference = GenotypePanel(
        ("R1:H1", "R1:H2", "R2:H1", "R2:H2"),
        variants,
        np.array([[0, 0, 1], [0, 1, 0], [1, 1, 0], [1, 0, 1]], dtype=float),
        1,
    )
    study = GenotypePanel(
        ("S1:H1", "S1:H2"), variants, np.array([[0, 0, 1], [1, 1, 0]], dtype=float), 1
    )

    result = harmonize_alt_dosages(
        reference, study, requested_components=2, minimum_variants=3
    )

    assert {variant.chromosome for variant in result.variants} == {"chr7"}
    assert result.study_scores.shape == (2, 2)


def test_allele_mismatch_and_insufficient_overlap_block() -> None:
    reference = GenotypePanel(
        ("R1", "R2"),
        (Variant("chr3", 10, "r", "A", "G"),),
        np.array([[0], [2]], dtype=float),
        2,
    )
    study = GenotypePanel(
        ("S1",),
        (Variant("chr3", 10, "s", "A", "T"),),
        np.array([[1]], dtype=float),
        2,
    )
    with pytest.raises(AncestryAnalysisError, match="insufficient_harmonized"):
        harmonize_alt_dosages(reference, study, requested_components=1, minimum_variants=1)


def test_phased_local_genotype_and_centroids() -> None:
    assert genotype_to_alt_dosage("1/0", phased_required=False) == 1
    assert phased_genotype_to_haplotypes("1|0") == (1, 0)
    with pytest.raises(AncestryAnalysisError, match="not_phased"):
        phased_genotype_to_haplotypes("0/1")
    centroids = population_centroids(np.array([[0.0, 1.0], [2.0, 3.0]]), ["P", "P"])
    assert centroids["P"][0] == 2
    np.testing.assert_allclose(centroids["P"][1], [1.0, 2.0])
