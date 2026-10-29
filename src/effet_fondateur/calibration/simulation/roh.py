"""D16 : autozygotie contrôlée par copie d'une origine H1 sur H2.

L'intervalle d'origine commune connu est distinct des bornes maximales d'un
ROH appelé par PLINK. Les marqueurs du fond peuvent prolonger l'homozygotie.
Le témoin sans injection n'est pas une garantie d'absence de ROH historique.
"""

from __future__ import annotations

from collections.abc import Iterable, Mapping, Sequence
from dataclasses import dataclass
import random
from typing import TextIO

from effet_fondateur.calibration.simulation.founder_recombination import (
    GeneticMapPoint, interpolate_cm,
)


METHOD_ID = "within_sample_h1_to_h2_roh_copy_v1"


class RohSimulationError(ValueError):
    """Erreur D16 avec code stable, sans génotype ou identifiant individuel."""


@dataclass(frozen=True)
class RohInjection:
    """Vérité autozygote d'un sujet ; bornes nulles en absence d'injection."""

    scenario_id: str
    sample_id: str
    chromosome: int
    is_positive: bool
    relationship_to_target: str
    source_ibd_scenario_id: str | None
    start_bp: int | None
    end_bp: int | None
    start_cm: float | None
    end_cm: float | None
    length_bp: int | None
    length_cm: float | None
    source_haplotype_index: int | None
    copied_haplotype_index: int | None


def derive_roh_ibd_edges(
    injections: Sequence[RohInjection], ibd_truth: Sequence[Mapping[str, object]],
    map_points: Sequence[GeneticMapPoint], *, replication_id: str,
) -> list[dict[str, object]]:
    """Enregistre le partage H2–H1 induit dans chaque paire positive couplée.

    L'intervalle est l'intersection du ROH et du segment IBD source. La vérité
    D15 H1–H1 reste conservée séparément ; cette surcharge évite d'interpréter
    un appel H2 connu comme faux positif. La matrice D15 reste une origine de
    simulation, pas une paire supplémentaire du protocole d'évaluation.
    """
    source_by_id = {item["scenario_id"]: item for item in ibd_truth}
    edges: list[dict[str, object]] = []
    for item in injections:
        if item.source_ibd_scenario_id is None:
            continue
        source = source_by_id[item.source_ibd_scenario_id]
        start = max(item.start_bp, source["start_bp"])
        end = min(item.end_bp, source["end_bp"])
        if start >= end or item.sample_id != source["sample_1_id"]:
            raise RohSimulationError("d16_induced_ibd_interval_invalid")
        edges.append({
            "replication_id": replication_id, "roh_scenario_id": item.scenario_id,
            "source_ibd_scenario_id": item.source_ibd_scenario_id,
            "sample_1_id": item.sample_id, "sample_2_id": source["sample_2_id"],
            "haplotype_1_index": 2, "haplotype_2_index": 1,
            "chromosome": item.chromosome, "start_bp": start, "end_bp": end,
            "start_cm": interpolate_cm(map_points, start),
            "end_cm": interpolate_cm(map_points, end),
            "is_positive": True, "origin_policy": "ROH_COPY_OF_D15_SHARED_H1",
        })
    return edges


def build_roh_injections(
    *, scenarios: Sequence[Mapping[str, object]], eligible_samples: Sequence[str],
    ibd_truth: Sequence[Mapping[str, object]], map_points: Sequence[GeneticMapPoint],
    chromosome: int, target_bp: int, random_seed: int,
) -> tuple[RohInjection, ...]:
    """Sélectionne les sujets sans observer leurs génotypes ou un appel PLINK.

    Un scénario couplé reprend sample_1 d'une paire positive D15 explicite.
    Les autres utilisent des sujets hors pedigree et hors toute paire/matrice
    D15, fournis par le matérialisateur. Les bornes restent celles de la
    configuration ; aucune densité insuffisante ne provoque leur déplacement.
    """
    if type(random_seed) is not int or not 0 <= random_seed <= 4294967295:
        raise RohSimulationError("d16_seed_invalid")
    if type(chromosome) is not int or not 1 <= chromosome <= 22:
        raise RohSimulationError("d16_chromosome_invalid")
    if not scenarios or len(set(eligible_samples)) != len(eligible_samples):
        raise RohSimulationError("d16_scenarios_or_samples_invalid")
    if len({item["scenario_id"] for item in scenarios}) != len(scenarios):
        raise RohSimulationError("d16_duplicate_scenario")
    ibd_by_id = {item["scenario_id"]: item for item in ibd_truth}
    if len(ibd_by_id) != len(ibd_truth):
        raise RohSimulationError("d16_duplicate_ibd_scenario")
    protected = {sample for item in ibd_truth for sample in
                 (item["sample_1_id"], item["sample_2_id"], item["template_id"]) if sample}
    if protected & set(eligible_samples):
        raise RohSimulationError("d16_background_overlaps_ibd_participants")
    samples = sorted(eligible_samples)
    random.Random(random_seed).shuffle(samples)
    unused_count = sum(item["source_ibd_scenario_id"] is None for item in scenarios)
    if len(samples) < unused_count:
        raise RohSimulationError("d16_eligible_samples_insufficient")
    result: list[RohInjection] = []
    used: set[str] = set()
    for scenario in scenarios:
        source_id = scenario["source_ibd_scenario_id"]
        if source_id is None:
            sample_id = samples.pop()
        else:
            source = ibd_by_id.get(source_id)
            if source is None or source["is_positive"] is not True:
                raise RohSimulationError("d16_positive_ibd_context_required")
            if source["chromosome"] != chromosome:
                raise RohSimulationError("d16_ibd_chromosome_mismatch")
            sample_id = str(source["sample_1_id"])
        if sample_id in used:
            raise RohSimulationError("d16_sample_reused")
        used.add(sample_id)
        relationship = scenario["relationship_to_target"]
        start_bp, end_bp = scenario["start_bp"], scenario["end_bp"]
        if relationship == "NO_INJECTION":
            if start_bp is not None or end_bp is not None or source_id is not None:
                raise RohSimulationError("d16_negative_bounds_or_context_invalid")
            result.append(RohInjection(str(scenario["scenario_id"]), sample_id,
                                       chromosome, False, str(relationship), None,
                                       None, None, None, None, None, None, None, None))
            continue
        if type(start_bp) is not int or type(end_bp) is not int or start_bp >= end_bp:
            raise RohSimulationError("d16_boundaries_invalid")
        if start_bp <= map_points[0].position_bp or end_bp >= map_points[-1].position_bp:
            raise RohSimulationError("d16_interval_outside_domain")
        target_inside = start_bp <= target_bp <= end_bp
        if relationship not in {"TARGET_DISTINCT", "TARGET_OVERLAP"} or (
            target_inside != (relationship == "TARGET_OVERLAP")
        ):
            raise RohSimulationError("d16_target_relationship_mismatch")
        if source_id is not None:
            source = ibd_by_id[source_id]
            if max(start_bp, source["start_bp"]) > min(end_bp, source["end_bp"]):
                raise RohSimulationError("d16_ibd_roh_overlap_missing")
        start_cm = interpolate_cm(map_points, start_bp)
        end_cm = interpolate_cm(map_points, end_bp)
        result.append(RohInjection(
            str(scenario["scenario_id"]), sample_id, chromosome, True,
            str(relationship), source_id, start_bp, end_bp, start_cm, end_cm,
            end_bp - start_bp + 1, end_cm - start_cm, 1, 2,
        ))
    return tuple(result)


def validate_roh_density(
    injections: Sequence[RohInjection], observed_positions: Sequence[int], *,
    minimum_markers: int, maximum_density_kb: float, maximum_gap_kb: float,
    minimum_flanks: int,
) -> list[dict[str, object]]:
    """Contrôle la détectabilité géométrique sur l'univers puce, en bp/kb.

    La densité est la longueur fermée injectée divisée par le nombre de sondes.
    Elle ne prédit pas le résultat d'un scan PLINK ; aucun seuil n'est ajusté.
    """
    positions = sorted(observed_positions)
    audit: list[dict[str, object]] = []
    for item in injections:
        if not item.is_positive:
            continue
        inside = [position for position in positions if item.start_bp <= position <= item.end_bp]
        if not inside:
            raise RohSimulationError("d16_observable_interval_empty")
        density = item.length_bp / 1000 / len(inside)
        gap = max((right - left for left, right in zip(inside, inside[1:])), default=0) / 1000
        left_flanks = sum(position < item.start_bp for position in positions)
        right_flanks = sum(position > item.end_bp for position in positions)
        if len(inside) < minimum_markers or density > maximum_density_kb or gap > maximum_gap_kb or min(left_flanks, right_flanks) < minimum_flanks:
            raise RohSimulationError("d16_marker_density_insufficient")
        audit.append({"scenario_id": item.scenario_id, "observed_marker_count": len(inside),
                      "density_kb_per_marker": density, "maximum_internal_gap_kb": gap,
                      "left_flank_marker_count": left_flanks,
                      "right_flank_marker_count": right_flanks})
    return audit


def transform_roh_vcf(
    source: Iterable[str], dense: TextIO, observed: TextIO, *,
    injections: Sequence[RohInjection], observed_keys: Mapping[tuple[str, int, str, str], str],
    transmissions: Sequence[tuple[str, str, int, int]],
    ibd_truth: Sequence[Mapping[str, object]],
) -> dict[str, int]:
    """Copie H1→H2 et contrôle l'homozygotie, Mendel et l'IBD H1 existant.

    Toutes les copies H1 sont conservées. H2 ne change que sur les intervalles
    positifs. Les copies négatives et matrices D15 restent intouchées.
    Le VCF observable est non phasé et ne contient aucune étiquette de vérité.
    """
    positive = [item for item in injections if item.is_positive]
    modified = {item.sample_id for item in positive}
    if len({item.sample_id for item in injections}) != len(injections):
        raise RohSimulationError("d16_duplicate_injection_sample")
    if modified & {sample for row in transmissions for sample in row[:2]}:
        raise RohSimulationError("d16_pedigree_modification_forbidden")
    protected = {sample for item in ibd_truth for sample in (
        (item["template_id"],) if item["is_positive"] else
        (item["sample_1_id"], item["sample_2_id"])) if sample}
    if modified & protected:
        raise RohSimulationError("d16_negative_or_template_ibd_modification")
    sample_index: dict[str, int] | None = None
    seen: set[tuple[str, int, str, str]] = set()
    counts = {"dense_variant_count": 0, "observed_variant_count": 0,
              "changed_h2_allele_count": 0, "h1_modified_count": 0,
              "roh_copy_error_count": 0, "ibd_copy_error_count": 0,
              "mendelian_error_count": 0, "negative_modified_count": 0,
              "induced_ibd_copy_error_count": 0}
    ibd_by_id = {item["scenario_id"]: item for item in ibd_truth}
    for line in source:
        if line.startswith("##"):
            if line.startswith(("##fileformat=", "##contig=", "##FORMAT=")):
                dense.write(line)
                observed.write(line)
            continue
        if line.startswith("#CHROM\t"):
            samples = line.rstrip().split("\t")[9:]
            if len(samples) != len(set(samples)):
                raise RohSimulationError("d16_duplicate_vcf_sample")
            sample_index = {sample: index for index, sample in enumerate(samples)}
            required = {item.sample_id for item in injections} | {
                sample for item in ibd_truth for sample in
                (item["sample_1_id"], item["sample_2_id"], item["template_id"]) if sample}
            required |= {sample for row in transmissions for sample in row[:2]}
            if not required <= set(sample_index):
                raise RohSimulationError("d16_sample_missing")
            dense.write(line)
            observed.write(line)
            continue
        if not line.strip():
            continue
        if sample_index is None:
            raise RohSimulationError("d16_header_missing")
        fields = line.rstrip().split("\t")
        if len(fields) != 9 + len(sample_index) or fields[8] != "GT":
            raise RohSimulationError("d16_vcf_columns_invalid")
        position = int(fields[1])
        if any(fields[0] not in {str(item.chromosome), f"chr{item.chromosome}"} for item in injections):
            raise RohSimulationError("d16_chromosome_mismatch")
        genotypes = [value.split("|") for value in fields[9:]]
        if any(len(gt) != 2 or any(allele not in {"0", "1"} for allele in gt) for gt in genotypes):
            raise RohSimulationError("d16_gt_invalid")
        for item in positive:
            if item.start_bp <= position <= item.end_bp:
                genotype = genotypes[sample_index[item.sample_id]]
                counts["changed_h2_allele_count"] += genotype[1] != genotype[0]
                genotype[1] = genotype[0]
                if genotype[0] != genotype[1]:
                    raise RohSimulationError("d16_roh_copy_error")
        for item in ibd_truth:
            if item["is_positive"] and item["start_bp"] <= position <= item["end_bp"]:
                alleles = {genotypes[sample_index[sample]][0] for sample in
                           (item["sample_1_id"], item["sample_2_id"], item["template_id"])}
                if len(alleles) != 1:
                    raise RohSimulationError("d16_existing_ibd_copy_error")
        for item in positive:
            if item.source_ibd_scenario_id is not None:
                source_ibd = ibd_by_id[item.source_ibd_scenario_id]
                if max(item.start_bp, source_ibd["start_bp"]) <= position <= min(item.end_bp, source_ibd["end_bp"]):
                    if genotypes[sample_index[item.sample_id]][1] != genotypes[sample_index[source_ibd["sample_2_id"]]][0]:
                        raise RohSimulationError("d16_induced_ibd_copy_error")
        for child, parent, parent_hap, child_hap in transmissions:
            if genotypes[sample_index[child]][child_hap - 1] != genotypes[sample_index[parent]][parent_hap - 1]:
                raise RohSimulationError("d16_mendelian_error")
        dense.write("\t".join(fields[:9] + ["|".join(gt) for gt in genotypes]) + "\n")
        counts["dense_variant_count"] += 1
        key = (fields[0], position, fields[3], fields[4])
        if key in observed_keys:
            if key in seen:
                raise RohSimulationError("d16_duplicate_observable_variant")
            fixed = fields[:9]
            fixed[2] = observed_keys[key]
            observed.write("\t".join(fixed + ["/".join(sorted(gt)) for gt in genotypes]) + "\n")
            seen.add(key)
            counts["observed_variant_count"] += 1
    if seen != set(observed_keys) or counts["dense_variant_count"] == 0:
        raise RohSimulationError("d16_variant_universe_mismatch")
    return counts
