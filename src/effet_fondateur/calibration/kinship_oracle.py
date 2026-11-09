"""Oracles D21 indépendants : métriques KING, graphe et incertitude familiale.

KING-robust : Manichaikul et al. (2010), équations 9 et 11,
https://www.cog-genomics.org/static/pdf/Manichaikuletal2010.pdf.
Les paires manquantes ne sont jamais imputées comme UNRELATED. Les seules
paires indépendantes pour l'incertitude sont les paires cibles d'unités
généalogiques séparées, une par unité ; les profils sont appariés.
"""

from __future__ import annotations

from collections import Counter
import math
from statistics import NormalDist
from typing import Any

import numpy as np

from effet_fondateur.calibration.simulation.kinship_bench import DEGREES
from effet_fondateur.calibration.technical_evaluation import GenotypeMatrix


def degree_and_relation(kinship: float, ibs0: float, parameters: dict[str, Any]) -> tuple[str, str]:
    """Reclasse les valeurs sérialisées, sans appeler le helper de production."""
    keys = ("duplicate_or_mz_min", "first_degree_min", "second_degree_min", "third_degree_min", "fourth_degree_min")
    for index, key in enumerate(keys):
        if kinship >= parameters[key]:
            if index == 1:
                return "FIRST", "PARENT_CHILD_COMPATIBLE" if ibs0 <= parameters["parent_child_ibs0_max"] else "FIRST_DEGREE_NON_PARENT_CHILD"
            return DEGREES[index], ("DUPLICATE_OR_MZ", "", "SECOND_DEGREE", "THIRD_DEGREE", "FOURTH_DEGREE")[index]
    return "UNRELATED", "UNRELATED"


def king_metrics(matrix: GenotypeMatrix, left: int, right: int) -> dict[str, int | float | None]:
    """Comptes appelés en commun et φ robuste intra/inter-FID, sans HWE ajusté.

    N_SNP est un compte ; HetHet et IBS0 sont des proportions. Hétérozygotie
    nulle => φ absent, jamais zéro fabriqué. Les allèles sont physiques.
    """
    first, second = matrix.alleles[:, left], matrix.alleles[:, right]
    called = (first[:, 0] >= 0) & (second[:, 0] >= 0)
    first, second = first[called], second[called]
    count = len(first)
    het_first = first[:, 0] != first[:, 1]
    het_second = second[:, 0] != second[:, 1]
    first_count, second_count = int(het_first.sum()), int(het_second.sum())
    both = int((het_first & het_second).sum())
    opposite = int((~het_first & ~het_second & (first[:, 0] != second[:, 0])).sum())
    denominator = first_count + second_count
    within = matrix.samples[left][0] == matrix.samples[right][0]
    if within:
        coefficient = (both - 2 * opposite) / denominator if denominator else None
    else:
        minimum = min(first_count, second_count)
        coefficient = .5 - (denominator - 2 * both + 4 * opposite) / (4 * minimum) if minimum else None
    return dict(N_SNP=count, HETHET=both / count if count else None, IBS0=opposite / count if count else None, KINSHIP=coefficient)


def independent_proposal(samples: list[str], edges: dict[tuple[str, str], set[str]], qc_rows: list[dict[str, str]]) -> tuple[dict[str, Any], dict[str, set[str]]]:
    """Oracle glouton explicite et contrôle de maximalité, pas d'optimisation max.

    Les priorités sont les entrées QC 05 fournies à 07 ; on ne réutilise
    aucune décision de proposition 07. Le graphe est construit par l'appelant.
    """
    neighbors = {sample: set() for sample in samples}
    for left, right in edges:
        neighbors[left].add(right)
        neighbors[right].add(left)
    quality = {row["SAMPLE_ID"]: (0 if row["QC_STATUS"] == "PASS" else 1, float(row["F_MISS"]),
                                  abs(float(row["HETEROZYGOSITY_Z"])) if row["HETEROZYGOSITY_Z"] else math.inf, row["SAMPLE_ID"]) for row in qc_rows}
    remaining = set(samples)
    selected, representatives = set(), {}
    # Le choix global par priorité est équivalent au choix composante par
    # composante : les composantes n'ont aucune arête commune.
    while remaining:
        winner = sorted(remaining, key=quality.__getitem__)[0]
        selected.add(winner)
        for blocked in remaining & neighbors[winner]:
            representatives[blocked] = winner
        remaining.difference_update(neighbors[winner] | {winner})
    return dict(selected=selected, representatives=representatives), neighbors


def wilson_interval(successes: int, count: int, confidence: float) -> list[float] | None:
    """Intervalle Wilson sur unités indépendantes ; dénominateur nul => None."""
    if count == 0:
        return None
    if not 0 <= successes <= count or not 0 < confidence < 1:
        raise ValueError("d21_interval_input_invalid")
    score = NormalDist().inv_cdf((1 + confidence) / 2)
    rate = successes / count
    divisor = 1 + score * score / count
    center = (rate + score * score / (2 * count)) / divisor
    radius = score * math.sqrt(rate * (1 - rate) / count + score * score / (4 * count * count)) / divisor
    return [max(0., center - radius), min(1., center + radius)]


def classification_metrics(records: list[dict[str, Any]], config: dict[str, Any]) -> dict[str, Any]:
    """Matrice complète, rappels Wilson et bootstrap par unité, sans seuil inventé."""
    columns = (*DEGREES, "MISSING")
    matrix = {truth: {observed: 0 for observed in columns} for truth in DEGREES}
    for record in records:
        matrix[record["true_degree"]][record["observed_degree"]] += 1
    recalls = {}
    for degree in DEGREES:
        count = sum(matrix[degree].values())
        correct = matrix[degree][degree]
        recalls[degree] = dict(unit_count=count, correct=correct, recall=correct / count if count else None,
                               wilson_interval=wilson_interval(correct, count, config["confidence_level"]))
    balanced = sum(value["recall"] for value in recalls.values()) / len(DEGREES) if all(value["recall"] is not None for value in recalls.values()) else None
    groups: dict[tuple[str, str], list[int]] = {}
    for record in records:
        groups.setdefault((record["stratum_id"], record["true_degree"]), []).append(int(record["true_degree"] == record["observed_degree"]))
    rng = np.random.default_rng(config["bootstrap_seed"])
    replicates = []
    for _ in range(config["bootstrap_replicates"]):
        correct_by_degree, total_by_degree = Counter(), Counter()
        for (_, degree), values in groups.items():
            drawn = rng.choice(values, len(values), replace=True)
            correct_by_degree[degree] += int(drawn.sum())
            total_by_degree[degree] += len(values)
        if len(total_by_degree) == len(DEGREES):
            replicates.append(sum(correct_by_degree[degree] / total_by_degree[degree] for degree in DEGREES) / len(DEGREES))
    tail = (1 - config["confidence_level"]) / 2
    return dict(confusion=matrix, per_degree=recalls, balanced_accuracy=balanced,
                family_bootstrap_interval=np.quantile(replicates, [tail, 1 - tail]).tolist() if replicates else None,
                bootstrap_seed=config["bootstrap_seed"], bootstrap_replicates=config["bootstrap_replicates"],
                independent_unit_count=len(records), missing_pair_count=sum(record["observed_degree"] == "MISSING" for record in records),
                qualification_minimum_per_class_met=all(value["unit_count"] >= 200 for value in recalls.values()))
