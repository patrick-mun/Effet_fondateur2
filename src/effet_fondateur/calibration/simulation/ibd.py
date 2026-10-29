"""D15 : copies haplotypiques partagées sur des intervalles préspecifiés.

La vérité décrit une origine commune injectée, pas les bornes maximales d'un
appel IBD. Le fond public peut prolonger une identité par état (IBS). Une paire
négative signifie absence d'injection, sans garantie d'absence d'IBD historique.
Les sujets choisis sont hors pedigree et hors copies porteuses D14.
"""

from __future__ import annotations

from collections.abc import Iterable, Mapping, Sequence
from dataclasses import dataclass
import math
import random
from typing import TextIO

from effet_fondateur.calibration.simulation.founder_recombination import (
    GeneticMapPoint, interpolate_bp, interpolate_cm,
)


METHOD_ID = "disjoint_unrelated_haplotype_copy_v1"


class IbdSimulationError(ValueError):
    """Violation d'un contrat D15, sans identifiant individuel dans l'erreur."""


@dataclass(frozen=True)
class IbdInjection:
    """Paire indépendante et intervalle vrai en bp/cM, ou bornes nulles."""

    scenario_id: str
    sample_1_id: str
    sample_2_id: str
    template_id: str | None
    is_positive: bool
    chromosome: int
    requested_length_cm: float
    start_bp: int | None
    end_bp: int | None
    start_cm: float | None
    end_cm: float | None
    minimum_markers: int


def build_ibd_injections(
    *, eligible_samples: Sequence[str], map_points: Sequence[GeneticMapPoint],
    chromosome: int, target_bp: int, scenarios: Sequence[Mapping[str, object]],
    random_seed: int,
) -> tuple[IbdInjection, ...]:
    """Tire cinq sujets distincts par longueur : deux positifs, deux négatifs,
    et une matrice privée. Tous les haplotypes injectés sont H1.

    La graine est explicitement dérivée du registre de recombinaison D05.
    Aucune sélection ne dépend des génotypes ni d'un appel IBD.
    """
    if type(random_seed) is not int or not 0 <= random_seed <= 4294967295:
        raise IbdSimulationError("d15_seed_invalid")
    if not scenarios or len(set(eligible_samples)) != len(eligible_samples):
        raise IbdSimulationError("d15_scenarios_or_samples_invalid")
    if len({scenario["scenario_id"] for scenario in scenarios}) != len(scenarios):
        raise IbdSimulationError("d15_duplicate_scenario")
    if type(chromosome) is not int or not 1 <= chromosome <= 22:
        raise IbdSimulationError("d15_chromosome_invalid")
    if len(eligible_samples) < 5 * len(scenarios):
        raise IbdSimulationError("d15_independent_samples_insufficient")
    samples = sorted(eligible_samples)
    random.Random(random_seed).shuffle(samples)
    target_cm = interpolate_cm(map_points, target_bp)
    injections: list[IbdInjection] = []
    for index, scenario in enumerate(scenarios):
        length = float(scenario["segment_length_cm"])
        if not math.isfinite(length) or length <= 0:
            raise IbdSimulationError("d15_length_invalid")
        start_cm, end_cm = target_cm - length / 2, target_cm + length / 2
        if start_cm <= map_points[0].position_cm or end_cm >= map_points[-1].position_cm:
            raise IbdSimulationError("d15_segment_outside_domain")
        start_bp, end_bp = interpolate_bp(map_points, start_cm), interpolate_bp(map_points, end_cm)
        if not start_bp < target_bp < end_bp:
            raise IbdSimulationError("d15_boundaries_invalid")
        realized_start, realized_end = interpolate_cm(map_points, start_bp), interpolate_cm(map_points, end_bp)
        first, second, negative_first, negative_second, template = samples[5 * index:5 * index + 5]
        common = (chromosome, length)
        injections.append(IbdInjection(
            str(scenario["scenario_id"]) + "_positive", first, second, template,
            True, *common, start_bp, end_bp, realized_start, realized_end,
            int(scenario["minimum_observed_markers"]),
        ))
        injections.append(IbdInjection(
            str(scenario["scenario_id"]) + "_negative", negative_first,
            negative_second, None, False, *common, None, None, None, None,
            int(scenario["minimum_observed_markers"]),
        ))
    return tuple(injections)


def transform_ibd_vcf(
    source: Iterable[str], dense: TextIO, observed: TextIO, *,
    injections: Sequence[IbdInjection],
    observed_keys: Mapping[tuple[str, int, str, str], str],
    transmissions: Sequence[tuple[str, str, int, int]],
) -> dict[str, int]:
    """Copie H1 matrice→H1 des deux sujets, puis contrôle l'identité et Mendel.

    Hors intervalle et sur toute paire négative, les GT restent identiques à
    l'entrée D14. La sortie dense conserve la phase ; l'observation conserve
    seulement l'univers puce et des GT non phasés. Erreur fatale si une sonde
    attendue manque, si un GT est invalide ou si Mendel est rompu.
    """
    positives = [item for item in injections if item.is_positive]
    participants = [sample for item in injections for sample in
                    (item.sample_1_id, item.sample_2_id, item.template_id) if sample]
    if len(set(participants)) != len(participants):
        raise IbdSimulationError("d15_pair_or_template_overlap")
    modified_samples = {sample for item in positives for sample in
                        (item.sample_1_id, item.sample_2_id)}
    if modified_samples & {sample for row in transmissions for sample in row[:2]}:
        raise IbdSimulationError("d15_pedigree_injection_forbidden")
    sample_index: dict[str, int] | None = None
    counts = {"dense_variant_count": 0, "observed_variant_count": 0,
              "changed_haplotype_allele_count": 0, "mendelian_error_count": 0,
              "positive_copy_error_count": 0, "negative_modified_count": 0}
    marker_counts = {item.scenario_id: 0 for item in positives}
    seen: set[tuple[str, int, str, str]] = set()
    for line in source:
        if line.startswith("##"):
            if line.startswith(("##fileformat=", "##contig=", "##FORMAT=")):
                dense.write(line)
                observed.write(line)
            continue
        if line.startswith("#CHROM\t"):
            fields = line.rstrip().split("\t")
            samples = fields[9:]
            if len(set(samples)) != len(samples):
                raise IbdSimulationError("d15_duplicate_sample")
            sample_index = {sample: index for index, sample in enumerate(samples)}
            required = {sample for item in injections for sample in
                        (item.sample_1_id, item.sample_2_id, item.template_id) if sample}
            if not required <= set(sample_index):
                raise IbdSimulationError("d15_sample_missing")
            dense.write(line)
            observed.write(line)
            continue
        if not line.strip():
            continue
        if sample_index is None:
            raise IbdSimulationError("d15_header_missing")
        fields = line.rstrip().split("\t")
        if len(fields) != 9 + len(sample_index) or fields[8] != "GT":
            raise IbdSimulationError("d15_vcf_columns_invalid")
        position = int(fields[1])
        genotypes = [value.split("|") for value in fields[9:]]
        if any(len(gt) != 2 or any(allele not in {"0", "1"} for allele in gt)
               for gt in genotypes):
            raise IbdSimulationError("d15_gt_invalid")
        original = [tuple(gt) for gt in genotypes]
        key = (fields[0], position, fields[3], fields[4])
        for item in positives:
            if item.start_bp <= position <= item.end_bp:
                allele = genotypes[sample_index[item.template_id]][0]
                for sample in (item.sample_1_id, item.sample_2_id):
                    genotype = genotypes[sample_index[sample]]
                    counts["changed_haplotype_allele_count"] += genotype[0] != allele
                    genotype[0] = allele
                if genotypes[sample_index[item.sample_1_id]][0] != genotypes[sample_index[item.sample_2_id]][0]:
                    raise IbdSimulationError("d15_positive_copy_error")
                if key in observed_keys:
                    marker_counts[item.scenario_id] += 1
        for item in injections:
            if not item.is_positive:
                for sample in (item.sample_1_id, item.sample_2_id):
                    index = sample_index[sample]
                    if tuple(genotypes[index]) != original[index]:
                        raise IbdSimulationError("d15_negative_modified")
        if any(fields[0] not in {str(item.chromosome), f"chr{item.chromosome}"} for item in injections):
            raise IbdSimulationError("d15_chromosome_mismatch")
        for child, parent, parent_hap, child_hap in transmissions:
            if genotypes[sample_index[child]][child_hap - 1] != genotypes[sample_index[parent]][parent_hap - 1]:
                raise IbdSimulationError("d15_mendelian_error")
        dense.write("\t".join(fields[:9] + ["|".join(gt) for gt in genotypes]) + "\n")
        counts["dense_variant_count"] += 1
        if key in observed_keys:
            if key in seen:
                raise IbdSimulationError("d15_duplicate_observed_variant")
            fixed = fields[:9]
            fixed[2] = observed_keys[key]
            observed.write("\t".join(fixed + ["/".join(sorted(gt)) for gt in genotypes]) + "\n")
            seen.add(key)
            counts["observed_variant_count"] += 1
    if seen != set(observed_keys) or counts["dense_variant_count"] == 0:
        raise IbdSimulationError("d15_variant_universe_mismatch")
    if any(marker_counts[item.scenario_id] < item.minimum_markers for item in positives):
        raise IbdSimulationError("d15_marker_density_insufficient")
    return counts
