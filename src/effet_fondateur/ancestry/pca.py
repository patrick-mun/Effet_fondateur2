"""PCA de référence déterministe et projection sans réajustement sur l'étude."""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from scipy.sparse.linalg import svds


class AncestryPcaError(ValueError):
    """Signale une matrice insuffisante ou numériquement invalide."""


@dataclass(frozen=True)
class ReferencePca:
    """Modèle ajusté exclusivement sur les échantillons publics de référence."""

    allele_frequencies: np.ndarray
    scales: np.ndarray
    loadings: np.ndarray
    eigenvalues: np.ndarray
    total_variance: float
    reference_scores: np.ndarray
    informative_variant_mask: np.ndarray


def _matrix(values: np.ndarray, name: str) -> np.ndarray:
    matrix = np.asarray(values, dtype=float)
    if matrix.ndim != 2 or not matrix.size or np.isinf(matrix).any():
        raise AncestryPcaError(f"invalid_ancestry_matrix:{name}")
    return matrix


def fit_reference_pca(
    reference_dosages: np.ndarray,
    *,
    requested_components: int,
    ploidy: int = 2,
) -> ReferencePca:
    """Ajuste une PCA exclusivement sur les entités de référence.

    ``reference_dosages`` a la forme ``(entités, variants)``. Ses valeurs sont
    des nombres de copies de l'allèle alternatif entre zéro et ``ploidy`` ; les
    valeurs manquantes ``NaN`` sont imputées à la fréquence calculée sur la
    référence. Les variants non polymorphes sont exclus, puis les dosages sont
    centrés et réduits par ``sqrt(ploidy * p * (1-p))``.

    Le modèle retourne fréquences, échelles, axes, valeurs propres, scores de
    référence et masque des variants informatifs. Le signe de chaque axe est
    rendu déterministe en orientant positivement son chargement de plus grande
    valeur absolue. Aucun échantillon d'étude ne doit être fourni ici ; il doit
    être projeté ensuite avec :func:`project_pca`. Voir
    ``docs/modules/reference_ancestry.md``.
    """
    matrix = _matrix(reference_dosages, "reference")
    if ploidy not in {1, 2} or requested_components < 1:
        raise AncestryPcaError("invalid_ancestry_pca_parameters")
    frequencies = np.nanmean(matrix, axis=0) / ploidy
    informative = np.isfinite(frequencies) & (frequencies > 0) & (frequencies < 1)
    if informative.sum() < requested_components:
        raise AncestryPcaError("insufficient_informative_ancestry_variants")
    selected = matrix[:, informative]
    selected_frequencies = frequencies[informative]
    imputed = np.where(np.isnan(selected), ploidy * selected_frequencies, selected)
    scales = np.sqrt(ploidy * selected_frequencies * (1 - selected_frequencies))
    standardized = (imputed - ploidy * selected_frequencies) / scales
    component_count = min(
        requested_components, standardized.shape[0] - 1, standardized.shape[1]
    )
    if component_count < 1:
        raise AncestryPcaError("insufficient_reference_samples_for_pca")
    minimum_dimension = min(standardized.shape)
    if component_count < minimum_dimension:
        _, singular_values, loadings_transposed = svds(
            standardized,
            k=component_count,
            which="LM",
            random_state=0,
        )
        order = np.argsort(singular_values)[::-1]
        singular_values = singular_values[order]
        loadings_transposed = loadings_transposed[order]
    else:
        _, singular_values, loadings_transposed = np.linalg.svd(
            standardized, full_matrices=False
        )
        singular_values = singular_values[:component_count]
        loadings_transposed = loadings_transposed[:component_count]
    loadings = loadings_transposed[:component_count].T
    for component in range(component_count):
        pivot = int(np.argmax(np.abs(loadings[:, component])))
        if loadings[pivot, component] < 0:
            loadings[:, component] *= -1
    scores = standardized @ loadings
    eigenvalues = singular_values[:component_count] ** 2 / max(1, matrix.shape[0] - 1)
    total_variance = float(np.sum(standardized**2) / max(1, matrix.shape[0] - 1))
    return ReferencePca(
        frequencies,
        scales,
        loadings,
        eigenvalues,
        total_variance,
        scores,
        informative,
    )


def project_pca(model: ReferencePca, dosages: np.ndarray, *, ploidy: int = 2) -> np.ndarray:
    """Projette des entités sans réajuster le modèle de référence.

    ``dosages`` suit le même ordre de variants que la matrice ayant servi à
    l'ajustement et a la forme ``(entités, variants)``. L'imputation, le
    centrage et la réduction réutilisent strictement les fréquences et échelles
    du modèle. Une dimension ou une ploïdie incompatible lève
    :class:`AncestryPcaError`.
    """
    matrix = _matrix(dosages, "projection")
    if matrix.shape[1] != model.allele_frequencies.size or ploidy not in {1, 2}:
        raise AncestryPcaError("ancestry_projection_variant_mismatch")
    selected = matrix[:, model.informative_variant_mask]
    frequencies = model.allele_frequencies[model.informative_variant_mask]
    imputed = np.where(np.isnan(selected), ploidy * frequencies, selected)
    standardized = (imputed - ploidy * frequencies) / model.scales
    return standardized @ model.loadings
