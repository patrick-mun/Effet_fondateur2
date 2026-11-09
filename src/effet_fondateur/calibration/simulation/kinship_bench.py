"""D21 : gene-drop mendélien IID sur pedigrees artificiels indépendants.

Les loci sont indépendants, sauf copies LD déclarées. Il ne s'agit pas de
mosaïques chromosomiques ni d'une simulation de la puce ou des populations.
Les coefficients théoriques viennent d'une récurrence exacte de pedigree,
indépendante de KING ; les origines réalisent une vérité d'identité par copie.
"""

from __future__ import annotations

from dataclasses import dataclass
from fractions import Fraction
from typing import Any

import numpy as np
from numpy.typing import NDArray


RELATIONS = ("DUPLICATE", "PARENT_CHILD", "FULL_SIBLINGS", "HALF_SIBLINGS", "FIRST_COUSINS", "COUSINS_ONCE_REMOVED", "UNRELATED")
DEGREES = ("DUPLICATE_OR_MZ", "FIRST", "SECOND", "THIRD", "FOURTH", "UNRELATED")
EXPECTED_PHI = dict(zip(RELATIONS, (Fraction(1, 2), Fraction(1, 4), Fraction(1, 4), Fraction(1, 8), Fraction(1, 16), Fraction(1, 32), Fraction(0)), strict=True))


@dataclass(frozen=True)
class KinshipSimulation:
    """Allèles 0/1, origines fondatrices et unités indépendantes privées.

    Axes (locus, sujet observé, copie). Les unités comportent exactement une
    paire cible, afin de pouvoir rééchantillonner de vraies unités séparées.
    """

    alleles: NDArray[np.int8]
    origins: NDArray[np.int32]
    units: list[dict[str, Any]]
    mendelian_error_count: int


def pedigree_template(relation: str) -> tuple[list[tuple[int, int] | None], tuple[int, int], int | None]:
    """Définit un DAG topologique ; le clone est une copie, pas une filiation."""
    if relation in {"DUPLICATE", "UNRELATED"}:
        return [None, None], (0, 1), 0 if relation == "DUPLICATE" else None
    if relation == "PARENT_CHILD":
        return [None, None, (0, 1)], (0, 2), None
    if relation == "FULL_SIBLINGS":
        return [None, None, (0, 1), (0, 1)], (2, 3), None
    if relation == "HALF_SIBLINGS":
        return [None, None, None, (0, 1), (0, 2)], (3, 4), None
    common = [None, None, None, None, (0, 1), (0, 1), (4, 2), (5, 3)]
    if relation == "FIRST_COUSINS":
        return common, (6, 7), None
    if relation == "COUSINS_ONCE_REMOVED":
        return common + [None, (7, 8)], (6, 9), None
    raise ValueError("d21_unknown_relation")


def pedigree_phi(parents: list[tuple[int, int] | None], observed: tuple[int, int], clone_of: int | None = None) -> Fraction:
    """φ(i,j)=(φ(p,j)+φ(m,j))/2 ; φ(i,i)=(1+φ(p,m))/2.

    Fondateurs non consanguins indépendants. Refuse références non
    topologiques ; les coefficients restent rationnels jusqu'à la sortie.
    """
    size = len(parents)
    coefficients = [[Fraction(0) for _ in range(size)] for _ in range(size)]
    for child, pair in enumerate(parents):
        if clone_of is not None and child == observed[1]:
            if not 0 <= clone_of < child:
                raise ValueError("d21_clone_order_invalid")
            for previous in range(child):
                coefficients[child][previous] = coefficients[previous][child] = coefficients[clone_of][previous]
            coefficients[child][child] = coefficients[clone_of][clone_of]
        elif pair is None:
            coefficients[child][child] = Fraction(1, 2)
        else:
            father, mother = pair
            if not 0 <= father < child or not 0 <= mother < child or father == mother:
                raise ValueError("d21_pedigree_order_invalid")
            for previous in range(child):
                coefficients[child][previous] = coefficients[previous][child] = (coefficients[father][previous] + coefficients[mother][previous]) / 2
            coefficients[child][child] = (1 + coefficients[father][mother]) / 2
    return coefficients[observed[0]][observed[1]]


def _gene_drop(parents: list[tuple[int, int] | None], observed: tuple[int, int], clone_of: int | None,
               frequency: float, configuration: dict[str, Any], rng: np.random.Generator,
               origin_offset: int) -> tuple[np.ndarray, np.ndarray, int, int]:
    """Transmet les deux copies et leurs origines, contrôle Mendel par origine."""
    per_chromosome = configuration["markers_per_autosome"]
    marker_count = 22 * per_chromosome
    source_offset, destination_offset = configuration["ld_copy_offsets"]
    sources = np.arange(22) * per_chromosome + source_offset
    destinations = np.arange(22) * per_chromosome + destination_offset
    node_calls, node_origins = [], []
    errors = 0
    for index, pair in enumerate(parents):
        if clone_of is not None and index == observed[1]:
            calls, origins = node_calls[clone_of].copy(), node_origins[clone_of].copy()
        elif pair is None:
            calls = rng.binomial(1, frequency, (marker_count, 2)).astype(np.int8)
            calls[configuration["monomorphic_offset"]] = 0
            calls[destinations] = calls[sources]
            origins = np.broadcast_to(np.arange(origin_offset, origin_offset + 2, dtype=np.int32), (marker_count, 2)).copy()
            origin_offset += 2
        else:
            calls, origins = np.empty((marker_count, 2), dtype=np.int8), np.empty((marker_count, 2), dtype=np.int32)
            for copy, parent in enumerate(pair):
                choices = rng.integers(0, 2, marker_count)
                choices[destinations] = choices[sources]
                calls[:, copy] = node_calls[parent][np.arange(marker_count), choices]
                origins[:, copy] = node_origins[parent][np.arange(marker_count), choices]
                errors += int(np.count_nonzero(~np.any(origins[:, copy, None] == node_origins[parent], axis=1)))
        node_calls.append(calls)
        node_origins.append(origins)
    return np.stack([node_calls[index] for index in observed], axis=1), np.stack([node_origins[index] for index in observed], axis=1), origin_offset, errors


def generate_simulation(configuration: dict[str, Any]) -> KinshipSimulation:
    """Génère les unités indépendantes, sans lecture de données ou outil externe."""
    rng = np.random.default_rng(configuration["simulation_seed"])
    calls, origins, units = [], [], []
    origin_offset, errors = 0, 0
    for stratum in configuration["strata"]:
        for relation in RELATIONS:
            parents, observed, clone_of = pedigree_template(relation)
            theoretical = pedigree_phi(parents, observed, clone_of)
            if theoretical != EXPECTED_PHI[relation]:
                raise ValueError("d21_pedigree_phi_reference_mismatch")
            for replicate in range(configuration["families_per_relation_per_stratum"]):
                genotypes, ancestry, origin_offset, mendel = _gene_drop(parents, observed, clone_of, stratum["allele_frequency"], configuration, rng, origin_offset)
                index = len(units)
                calls.append(genotypes)
                origins.append(ancestry)
                errors += mendel
                degree = DEGREES[0] if theoretical == Fraction(1, 2) else DEGREES[1] if theoretical == Fraction(1, 4) else DEGREES[2] if theoretical == Fraction(1, 8) else DEGREES[3] if theoretical == Fraction(1, 16) else DEGREES[4] if theoretical == Fraction(1, 32) else DEGREES[5]
                units.append(dict(unit_index=index, stratum_id=stratum["stratum_id"], relation=relation, degree=degree,
                                  theoretical_phi=float(theoretical), parents=parents, observed_nodes=observed, clone_of=clone_of,
                                  sample_indices=[2 * index, 2 * index + 1], cryptic=relation == "DUPLICATE" and replicate == 0))
    if errors:
        raise ValueError("d21_gene_drop_mendelian_error")
    return KinshipSimulation(np.concatenate(calls, axis=1), np.concatenate(origins, axis=1), units, errors)
