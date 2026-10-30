"""Oracles D23 indépendants de SHAPEIT5 et des résumés de production.

Unités : positions bp, taux sans unité. La permutation globale des copies
est libre ; l'ancrage de phase exclut toujours le variant cible.
"""
from __future__ import annotations

from dataclasses import dataclass
from typing import Any, Iterable

import numpy as np


@dataclass(frozen=True)
class PhasePanel:
    """GT bialléliques par locus/sujet/copie ; -1=manquant, PS vide=chromosome."""
    samples: tuple[str, ...]
    loci: tuple[tuple[int, int, str, str], ...]
    alleles: np.ndarray
    phased: np.ndarray
    phase_sets: np.ndarray
    variant_ids: tuple[str, ...] = ()


def parse_vcf(lines: Iterable[str]) -> PhasePanel:
    """Lit GT et PS sans confondre un génotype non phasé avec une phase vraie."""
    samples: tuple[str, ...] = ()
    loci, alleles, phased, phase_sets, variant_ids = [], [], [], [], []
    for line in lines:
        if line.startswith("##"):
            continue
        fields = line.rstrip().split("\t")
        if line.startswith("#CHROM"):
            samples = tuple(fields[9:])
            if not samples or len(set(samples)) != len(samples):
                raise ValueError("d23_vcf_samples_invalid")
            continue
        if not samples or len(fields) != 9 + len(samples):
            raise ValueError("d23_vcf_columns_invalid")
        chromosome, position = int(fields[0].removeprefix("chr")), int(fields[1])
        if not 1 <= chromosome <= 22 or position <= 0 or fields[3] not in ("A", "C", "G", "T") or fields[4] not in ("A", "C", "G", "T") or fields[3] == fields[4]:
            raise ValueError("d23_vcf_locus_invalid")
        formats = fields[8].split(":")
        if "GT" not in formats:
            raise ValueError("d23_vcf_gt_missing")
        gt_index = formats.index("GT")
        ps_index = formats.index("PS") if "PS" in formats else None
        calls, orientations, blocks = [], [], []
        for field in fields[9:]:
            values = field.split(":")
            genotype = values[gt_index]
            is_phased = "|" in genotype
            parts = genotype.replace("|", "/").split("/")
            if len(parts) != 2 or any(value not in {"0", "1", "."} for value in parts):
                raise ValueError("d23_vcf_gt_invalid")
            calls.append([int(value) if value != "." else -1 for value in parts])
            orientations.append(is_phased)
            block = values[ps_index] if ps_index is not None and ps_index < len(values) else ""
            blocks.append(f"UNKNOWN_AT_{chromosome}_{position}" if ps_index is not None and block in {"", "."} else block)
        loci.append((chromosome, position, fields[3], fields[4]))
        variant_ids.append(fields[2])
        alleles.append(calls)
        phased.append(orientations)
        phase_sets.append(blocks)
    if not loci or len({locus[:2] for locus in loci}) != len(loci):
        raise ValueError("d23_vcf_empty_or_duplicate_loci")
    return PhasePanel(samples, tuple(loci), np.asarray(alleles, dtype=np.int8), np.asarray(phased, dtype=bool), np.asarray(phase_sets, dtype=object), tuple(variant_ids))


def align_prediction(truth: PhasePanel, predicted: PhasePanel) -> tuple[np.ndarray, np.ndarray, np.ndarray, dict[str, int]]:
    """Aligne par identité, recode REF/ALT inversés et refuse univers inattendus."""
    if set(truth.samples) != set(predicted.samples):
        raise ValueError("d23_prediction_sample_universe_mismatch")
    sample_order = [predicted.samples.index(sample) for sample in truth.samples]
    truth_positions = {locus[:2]: (index, locus[2:]) for index, locus in enumerate(truth.loci)}
    calls = np.full_like(truth.alleles, -1)
    phased = np.zeros(truth.alleles.shape[:2], dtype=bool)
    blocks = np.full(truth.alleles.shape[:2], "", dtype=object)
    swapped = 0
    for source, locus in enumerate(predicted.loci):
        if locus[:2] not in truth_positions:
            raise ValueError("d23_unexpected_predicted_locus")
        destination, pair = truth_positions[locus[:2]]
        value = predicted.alleles[source, sample_order].copy()
        if pair == locus[2:][::-1]:
            value[value >= 0] = 1 - value[value >= 0]
            swapped += 1
        elif pair != locus[2:]:
            raise ValueError("d23_incompatible_predicted_alleles")
        calls[destination] = value
        phased[destination] = predicted.phased[source, sample_order]
        blocks[destination] = predicted.phase_sets[source, sample_order]
    return calls, phased, blocks, dict(absent_loci=len(truth.loci) - len(predicted.loci), swapped_loci=swapped,
                                     sample_order_changed=int(truth.samples != predicted.samples))


def bootstrap_ratio(numerators: np.ndarray, denominators: np.ndarray, *, seed: int, replicates: int) -> dict[str, Any]:
    """IC percentile familial descriptif ; jamais un IC IID par génotype."""
    if numerators.ndim != 1 or numerators.shape != denominators.shape or np.any(numerators < 0) or np.any(denominators < numerators) or replicates < 1:
        raise ValueError("d23_bootstrap_counts_invalid")
    numerator, denominator = int(numerators.sum()), int(denominators.sum())
    contributing = int(np.count_nonzero(denominators))
    result = dict(numerator=numerator, denominator=denominator, estimate=numerator / denominator if denominator else None,
                  status="EVALUABLE" if denominator else "NOT_EVALUABLE", contributing_units=contributing,
                  interval_95=None, bootstrap_valid_replicates=0, interval_method="FAMILY_CLUSTER_PERCENTILE_DESCRIPTIVE")
    if denominator == 0 or contributing < 2:
        result["interval_reason"] = "ZERO_DENOMINATOR" if denominator == 0 else "FEWER_THAN_TWO_CONTRIBUTING_UNITS"
        return result
    indices = np.random.default_rng(seed).integers(len(numerators), size=(replicates, len(numerators)))
    totals = denominators[indices].sum(axis=1)
    values = numerators[indices].sum(axis=1)[totals > 0] / totals[totals > 0]
    result.update(interval_95=np.quantile(values, [.025, .975]).tolist(), bootstrap_valid_replicates=len(values))
    return result


def score_phase(truth: PhasePanel, predicted: PhasePanel, *, founder_indices: list[int], unit_indices: list[int], target: tuple[int, int],
                rare_maf_upper_bound: float, density_gap_bp: int, bootstrap_seed: int, bootstrap_replicates: int) -> dict[str, Any]:
    """Mesure switches, Hamming, motifs isolés, phase rare et orientation cible.

    La liste de loci vrais est l'univers observé du profil considéré. Une
    exclusion de phasage interrompt une chaîne d'hétérozygotes ; les homozygotes
    ne contribuent jamais au dénominateur de phase. Les 24 unités familiales
    sont rééchantillonnées entières et de manière identique entre métriques.
    """
    if len(unit_indices) != len(truth.samples) or not founder_indices or min(unit_indices) < 0 or sorted(set(unit_indices)) != list(range(max(unit_indices) + 1)):
        raise ValueError("d23_independence_units_invalid")
    calls, phased, blocks, alignment = align_prediction(truth, predicted)
    if np.any(truth.alleles < 0) or not np.all(truth.phased):
        raise ValueError("d23_complete_phased_truth_required")
    unit_count = max(unit_indices) + 1
    frequencies = truth.alleles[:, founder_indices].mean(axis=(1, 2))
    maf = np.minimum(frequencies, 1 - frequencies)
    positions = np.array([locus[1] for locus in truth.loci])
    chromosomes = np.array([locus[0] for locus in truth.loci])
    target_mask = np.array([locus[:2] == target for locus in truth.loci])
    has_previous = np.zeros(len(positions), dtype=bool)
    has_previous[1:] = chromosomes[1:] == chromosomes[:-1]
    dense = np.zeros(len(positions), dtype=bool)
    dense[1:] = (np.diff(positions) <= density_gap_bp) & has_previous[1:]
    strata = dict(ALL=np.ones(len(positions), dtype=bool), MAF_LE_005=maf <= rare_maf_upper_bound,
                  MAF_GT_005=maf > rare_maf_upper_bound, GAP_LE_100KB=dense, GAP_GT_100KB=has_previous & ~dense)
    metrics: dict[tuple[str, str], np.ndarray] = {}
    for name in ("switch_error_rate", "haplotype_hamming_rate", "isolated_point_motif_rate", "rare_variant_phase_error_rate", "target_mutation_orientation_accuracy"):
        for stratum in strata:
            metrics[name, stratum] = np.zeros((unit_count, 2), dtype=np.int64)
    def count(name: str, unit: int, locus: int, error: bool) -> None:
        for label, selection in strata.items():
            if selection[locus]:
                metrics[name, label][unit] += (int(error), 1)
    complete = np.all(calls >= 0, axis=2)
    dosage_ok = complete & (calls.sum(axis=2) == truth.alleles.sum(axis=2))
    het = truth.alleles[:, :, 0] != truth.alleles[:, :, 1]
    valid = dosage_ok & phased & het
    runs = dict(length_1=0, length_2_to_5=0, length_6_plus=0)
    for sample, unit in enumerate(unit_indices):
        domains = dict.fromkeys((locus[0], blocks[index, sample]) for index, locus in enumerate(truth.loci))
        for chromosome, phase_set in domains:
            heterozygotes = np.flatnonzero(het[:, sample] & (chromosomes == chromosome) & (blocks[:, sample] == phase_set))
            eligible = heterozygotes[valid[heterozygotes, sample]]
            orientation = (calls[:, sample, 0] != truth.alleles[:, sample, 0]).astype(np.int8)
            anchors = eligible[~target_mask[eligible]]
            votes = int(orientation[anchors].sum())
            tied = 2 * votes == len(anchors)
            inversion = int(2 * votes > len(anchors))
            if not tied:
                for locus in eligible:
                    error = bool(orientation[locus] != inversion)
                    count("haplotype_hamming_rate", unit, locus, error)
                    if maf[locus] <= rare_maf_upper_bound:
                        count("rare_variant_phase_error_rate", unit, locus, error)
                    if target_mask[locus]:
                        count("target_mutation_orientation_accuracy", unit, locus, not error)
            continuous = []
            # Un hétérozygote d'un autre PS interposé ne doit pas être ponté.
            chromosome_hets = np.flatnonzero(het[:, sample] & (chromosomes == chromosome))
            ranks = {int(locus): rank for rank, locus in enumerate(chromosome_hets)}
            for left, right in zip(heterozygotes[:-1], heterozygotes[1:]):
                connected = bool(valid[left, sample] and valid[right, sample] and ranks[int(right)] == ranks[int(left)] + 1)
                continuous.append(connected)
                if connected:
                    count("switch_error_rate", unit, right, orientation[left] != orientation[right])
            for index in range(1, len(heterozygotes) - 1):
                if continuous[index - 1] and continuous[index]:
                    left, center, right = heterozygotes[index - 1:index + 2]
                    count("isolated_point_motif_rate", unit, center, orientation[left] == orientation[right] and orientation[center] != orientation[left])
            if not tied:
                length = 0
                for index, locus in enumerate(heterozygotes):
                    if index and not continuous[index - 1] and length:
                        runs["length_1" if length == 1 else "length_2_to_5" if length <= 5 else "length_6_plus"] += 1
                        length = 0
                    if valid[locus, sample] and orientation[locus] != inversion:
                        length += 1
                    elif length:
                        runs["length_1" if length == 1 else "length_2_to_5" if length <= 5 else "length_6_plus"] += 1
                        length = 0
                if length:
                    runs["length_1" if length == 1 else "length_2_to_5" if length <= 5 else "length_6_plus"] += 1
    results = {name: {label: bootstrap_ratio(values[:, 0], values[:, 1], seed=bootstrap_seed, replicates=bootstrap_replicates)
                      for (metric, label), values in metrics.items() if metric == name} for name in sorted({key[0] for key in metrics})}
    return dict(alignment=alignment, samples=len(truth.samples), observed_loci=len(truth.loci), independent_units=unit_count,
                complete_genotypes=int(complete.sum()), missing_genotypes=int((~complete).sum()), genotype_discordances=int((complete & ~dosage_ok).sum()),
                truth_heterozygotes=int(het.sum()), phased_heterozygotes=int(valid.sum()), unphased_heterozygotes=int((het & dosage_ok & ~phased).sum()),
                target_heterozygotes=int((het & target_mask[:, None]).sum()), metrics=results, error_run_lengths=runs,
                qualification_status="REVIEW_REQUIRED", calibrated=False, founder_effect_proven=False)
