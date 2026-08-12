import numpy as np
import pytest

from effet_fondateur.ancestry.pca import (
    AncestryPcaError,
    fit_reference_pca,
    project_pca,
)


def test_study_projection_cannot_change_reference_axes() -> None:
    reference = np.array(
        [[0, 0, 0, 1], [0, 0, 1, 0], [2, 2, 2, 1], [2, 2, 1, 2]],
        dtype=float,
    )
    study = np.array([[0, 0, 0, np.nan], [2, 2, 2, 2]], dtype=float)

    model = fit_reference_pca(reference, requested_components=2)
    before = model.loadings.copy()
    projected = project_pca(model, study)

    assert projected.shape == (2, 2)
    np.testing.assert_array_equal(model.loadings, before)
    assert projected[0, 0] * projected[1, 0] < 0


def test_haplotype_projection_uses_haploid_scaling() -> None:
    reference_haplotypes = np.array(
        [[0, 0, 1], [0, 1, 0], [1, 1, 0], [1, 0, 1]], dtype=float
    )
    model = fit_reference_pca(
        reference_haplotypes, requested_components=2, ploidy=1
    )

    scores = project_pca(model, np.array([[0, 0, 1]]), ploidy=1)

    assert scores.shape == (1, 2)


def test_projection_with_different_variant_set_blocks() -> None:
    model = fit_reference_pca(
        np.array([[0, 0], [0, 1], [2, 1], [2, 2]], dtype=float),
        requested_components=1,
    )

    with pytest.raises(AncestryPcaError, match="variant_mismatch"):
        project_pca(model, np.array([[0, 1, 2]], dtype=float))
