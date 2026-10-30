"""Oracles de phase : exemples analytiques indépendants du logiciel de phasage."""
from dataclasses import replace

import numpy as np
import pytest

from effet_fondateur.calibration.phasing_oracle import PhasePanel, align_prediction, bootstrap_ratio, parse_vcf, score_phase


def fixture() -> PhasePanel:
    calls = np.tile(np.array([[[0, 1], [0, 1]]], dtype=np.int8), (7, 1, 1))
    return PhasePanel(("S1", "S2"), tuple((5, index + 1, "A", "G") for index in range(7)), calls,
                      np.ones((7, 2), dtype=bool), np.full((7, 2), "", dtype=object))


def score(predicted: PhasePanel, truth: PhasePanel | None = None) -> dict:
    return score_phase(truth or fixture(), predicted, founder_indices=[0, 1], unit_indices=[0, 1], target=(5, 4),
                       rare_maf_upper_bound=.05, density_gap_bp=100_000, bootstrap_seed=2305, bootstrap_replicates=200)


@pytest.mark.parametrize("flips,switches,points,hamming,target_correct", [([], 0, 0, 0, 2), (list(range(7)), 0, 0, 0, 2),
                          ([3, 4, 5, 6], 2, 0, 6, 0), ([3], 4, 2, 2, 0)])
def test_analytic_phase_examples(flips: list[int], switches: int, points: int, hamming: int, target_correct: int) -> None:
    original = fixture()
    calls = original.alleles.copy()
    calls[flips] = calls[flips, :, ::-1]
    result = score(replace(original, alleles=calls))
    metrics = result["metrics"]
    assert metrics["switch_error_rate"]["ALL"]["numerator"] == switches
    assert metrics["switch_error_rate"]["ALL"]["denominator"] == 12
    assert metrics["isolated_point_motif_rate"]["ALL"]["numerator"] == points
    # Dans l'exemple de bloc 3..6, les six ancrages hors cible sont à égalité.
    if flips == [3, 4, 5, 6]:
        assert metrics["haplotype_hamming_rate"]["ALL"]["denominator"] == 0
        assert metrics["target_mutation_orientation_accuracy"]["ALL"]["status"] == "NOT_EVALUABLE"
    else:
        assert metrics["haplotype_hamming_rate"]["ALL"]["numerator"] == hamming
        assert metrics["target_mutation_orientation_accuracy"]["ALL"]["numerator"] == target_correct


def test_missing_dosage_and_unphased_break_connections() -> None:
    original = fixture()
    calls = original.alleles.copy()
    calls[2] = -1
    calls[4] = 0
    phased = original.phased.copy()
    phased[6] = False
    result = score(replace(original, alleles=calls, phased=phased))
    assert result["missing_genotypes"] == 2
    assert result["genotype_discordances"] == 2
    assert result["unphased_heterozygotes"] == 2
    assert result["metrics"]["switch_error_rate"]["ALL"]["denominator"] == 2


def test_homozygotes_not_phase_observations_and_ps_free_orientation() -> None:
    original = fixture()
    calls = original.alleles.copy()
    calls[2] = 0
    truth = replace(original, alleles=calls)
    predicted = calls.copy()
    predicted[4:] = predicted[4:, :, ::-1]
    blocks = original.phase_sets.copy()
    blocks[4:] = "SECOND"
    result = score(replace(truth, alleles=predicted, phase_sets=blocks), truth)
    assert result["truth_heterozygotes"] == 12
    assert result["metrics"]["switch_error_rate"]["ALL"]["numerator"] == 0
    assert result["metrics"]["haplotype_hamming_rate"]["ALL"]["numerator"] == 0


def test_alignment_permutations_and_allele_swap() -> None:
    original = fixture()
    predicted = PhasePanel(original.samples[::-1], tuple((chrom, pos, alt, ref) for chrom, pos, ref, alt in original.loci[::-1]),
                           1-original.alleles[::-1, ::-1], original.phased[::-1, ::-1], original.phase_sets[::-1, ::-1])
    calls, _, _, info = align_prediction(original, predicted)
    np.testing.assert_array_equal(calls, original.alleles)
    assert info == dict(absent_loci=0, swapped_loci=7, sample_order_changed=1)
    with pytest.raises(ValueError, match="sample_universe"):
        align_prediction(original, replace(predicted, samples=("OTHER", "S2")))
    with pytest.raises(ValueError, match="incompatible"):
        align_prediction(original, replace(predicted, loci=((5, 7, "T", "C"),) + predicted.loci[1:]))


def test_bootstrap_units_reproducible_and_empty_denominators() -> None:
    first = bootstrap_ratio(np.array([1, 0, 1]), np.array([2, 0, 2]), seed=2305, replicates=100)
    assert first == bootstrap_ratio(np.array([1, 0, 1]), np.array([2, 0, 2]), seed=2305, replicates=100)
    assert first["estimate"] == .5
    assert first["interval_95"] == [.5, .5]
    assert bootstrap_ratio(np.array([0]), np.array([0]), seed=1, replicates=10)["status"] == "NOT_EVALUABLE"
    assert bootstrap_ratio(np.array([0]), np.array([2]), seed=1, replicates=10)["interval_95"] is None


def test_parser_preserves_gt_phase_and_rejects_duplicates() -> None:
    lines = ["#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\tS2\n",
             "chr5\t1\tm1\tA\tG\t.\tPASS\t.\tGT:PS\t0|1:1\t0/1:.\n"]
    parsed = parse_vcf(lines)
    assert parsed.phased.tolist() == [[True, False]]
    assert parsed.phase_sets[0, 1] == "UNKNOWN_AT_5_1"
    with pytest.raises(ValueError, match="duplicate"):
        parse_vcf(lines + lines[1:])


def test_target_without_independent_anchor_not_evaluable() -> None:
    original = fixture()
    one = PhasePanel(original.samples, (original.loci[3],), original.alleles[3:4], original.phased[3:4], original.phase_sets[3:4])
    assert score(one, one)["metrics"]["target_mutation_orientation_accuracy"]["ALL"]["denominator"] == 0


def test_rare_error_uses_founder_frequency_and_not_independent_gt_count() -> None:
    calls = np.tile(np.array([[[0, 1]]], dtype=np.int8), (7, 32, 1))
    calls[0] = 0
    calls[0, 0] = (0, 1)
    truth = PhasePanel(tuple(f"S{index}" for index in range(32)), tuple((5,index+1,"A","G") for index in range(7)), calls,
                       np.ones((7,32), dtype=bool), np.full((7,32), "", dtype=object))
    prediction = calls.copy()
    prediction[0,0] = (1,0)
    result = score_phase(truth, replace(truth, alleles=prediction), founder_indices=list(range(32)), unit_indices=list(range(32)), target=(5,4),
                         rare_maf_upper_bound=.05, density_gap_bp=100_000, bootstrap_seed=2305, bootstrap_replicates=100)
    rare = result["metrics"]["rare_variant_phase_error_rate"]["ALL"]
    assert rare["numerator"] == rare["denominator"] == 1
    assert rare["contributing_units"] == 1 and rare["interval_95"] is None
    assert result["metrics"]["rare_variant_phase_error_rate"]["GAP_GT_100KB"]["denominator"] == 0
