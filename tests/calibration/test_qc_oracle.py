"""Tests purs D20 : aucune invocation d'outil génétique réel."""

from copy import deepcopy
import itertools
import math
from pathlib import Path

import numpy as np
import pytest

from effet_fondateur.calibration.materialize_ibd import _json
from effet_fondateur.calibration.qc_bench import CASES, generate_cases, with_target
from effet_fondateur.calibration.qc_oracle import COHORTS, confusion, exact_hwe, final_expectation, preliminary_expectation
from effet_fondateur.calibration.technical_evaluation import GenotypeMatrix, TechnicalEvaluationError, observation_matrix


CONFIG = Path(__file__).resolve().parents[2] / "calibration/configurations/d20_qc_bench_v1.json"


@pytest.fixture(scope="module")
def bench():
    config = _json(CONFIG)
    cases = generate_cases(config)
    return config, {case: (panel, observation_matrix(panel, [dict(SAMPLE_ID=sample, FID=sample, IID=sample) for sample in panel.samples]), pairs)
                    for case, panel, pairs in cases}


@pytest.mark.parametrize("hets,hom_a,hom_b,expected", [(0, 0, 0, None), (0, 3, 0, 1.0), (3, 0, 0, .4), (0, 1, 1, 1 / 3),
                                                    (96, 0, 0, 2.373074178759543e-28)])
def test_hwe_known_exact_values(hets, hom_a, hom_b, expected):
    observed = exact_hwe(hets, hom_a, hom_b)
    assert observed is None if expected is None else observed == pytest.approx(expected, rel=1e-12, abs=0)


def test_hwe_matches_small_integer_enumeration():
    for counts in itertools.product(range(5), repeat=3):
        if not sum(counts):
            continue
        hets, hom_a, hom_b = counts
        total = sum(counts)
        copies = 2 * min(hom_a, hom_b) + hets
        weights = {heterozygotes: 2 ** heterozygotes / (math.factorial(heterozygotes)
                   * math.factorial((copies - heterozygotes) // 2)
                   * math.factorial(total - (copies + heterozygotes) // 2))
                   for heterozygotes in range(copies % 2, copies + 1, 2)}
        expected = sum(weight for weight in weights.values() if weight <= weights[hets] * (1 + 1e-12)) / sum(weights.values())
        assert exact_hwe(*counts) == pytest.approx(expected, abs=1e-12)


@pytest.mark.parametrize("counts", [(-1, 2, 3), (1.5, 0, 0), (True, 0, 0)])
def test_hwe_rejects_invalid_counts(counts):
    with pytest.raises(ValueError, match="counts_invalid"):
        exact_hwe(*counts)


def test_generator_is_paired_reproducible_and_autosomal(bench):
    config, cases = bench
    repeated = generate_cases(config)
    assert tuple(cases) == CASES
    all_ids = []
    for case, panel, _ in repeated:
        assert np.array_equal(panel.alleles, cases[case][0].alleles)
        assert {int(row[0]) for row in panel.variants} == set(range(1, 23))
        assert panel.alleles.shape == (1408, 192, 2)
        all_ids.extend(panel.samples)
    assert len(set(all_ids)) == len(all_ids)
    assert np.array_equal(cases["SAMPLE_MISSING"][0].alleles[:, 1:], cases["CLEAN"][0].alleles[:, 1:])


@pytest.mark.parametrize("field,value", [("seed", None), ("sample_count", 4), ("unexpected", 1), ("target_configuration", "../secret.json")])
def test_configuration_rejects_missing_unknown_or_unsafe(bench, field, value):
    config = deepcopy(bench[0])
    config[field] = value
    with pytest.raises(TechnicalEvaluationError):
        generate_cases(config)


@pytest.mark.parametrize("case,positive_field,positive_code", [
    ("SAMPLE_MISSING", "samples", "EXCLUDED"), ("VARIANT_MISSING", "variants", "EXCLUDED"),
    ("HET_OUTLIER", "samples", "ALERT"), ("BATCH_MISSING", "variants", "ALERT"), ("RARE_VARIANT", "variants", "ALERT")])
def test_preliminary_defects_have_known_positive_and_negative_units(bench, case, positive_field, positive_code):
    config, cases = bench
    matrix = cases[case][1]
    expected = preliminary_expectation(matrix, ["A", "B"] * 96, config["preliminary_parameters"])
    rows = getattr(expected, positive_field)
    assert any(row["QC_STATUS"] == positive_code for row in rows)
    assert any(row["QC_STATUS"] == "PASS" for row in rows)


def test_preliminary_exclusions_are_simultaneous_not_iterative(bench):
    config, cases = bench
    original = cases["CLEAN"][1]
    alleles = original.alleles.copy()
    alleles[:200, 0] = -1
    alleles[0, :10] = -1
    matrix = GenotypeMatrix(original.samples, original.variants, alleles)
    result = preliminary_expectation(matrix, ["A", "B"] * 96, config["preliminary_parameters"])
    assert result.samples[0]["QC_STATUS"] == "EXCLUDED"
    assert result.variants[0]["QC_STATUS"] == "EXCLUDED"  # 10/192, not 9/191 after subject exclusion


def test_control_specific_defect_passes_05_fails_10_not_carriers(bench):
    config, cases = bench
    matrix = cases["CONTROL_VARIANT_MISSING"][1]
    preliminary = preliminary_expectation(matrix, ["A", "B"] * 96, config["preliminary_parameters"])
    assert preliminary.variants[14]["QC_STATUS"] == "PASS"
    controls = GenotypeMatrix(matrix.samples[:96], matrix.variants, matrix.alleles[:, :96])
    carriers = GenotypeMatrix(matrix.samples[96:], matrix.variants, matrix.alleles[:, 96:])
    assert final_expectation(controls, ["A", "B"] * 48, COHORTS[0], config["final_parameters"], "target").variants[14]["FILTER_CODES"] == "variant_missingness"
    assert final_expectation(carriers, ["A", "B"] * 48, COHORTS[1], config["final_parameters"], "target").variants[14]["QC_STATUS"] == "PASS"


def test_hwe_only_controls_and_no_carrier_maf_filter(bench):
    config, cases = bench
    matrix = cases["CONTROL_HWE"][1]
    controls = GenotypeMatrix(matrix.samples[:96], matrix.variants, matrix.alleles[:, :96])
    result = final_expectation(controls, ["A", "B"] * 48, COHORTS[0], config["final_parameters"], "target")
    assert result.variants[13]["FILTER_CODES"] == "hwe_controls"
    result = final_expectation(controls, ["A", "B"] * 48, COHORTS[1], config["final_parameters"], "target")
    assert result.variants[13]["HWE_P"] is None
    assert result.variants[13]["QC_STATUS"] == "PASS"


def test_local_sample_and_target_exception_are_evaluable(bench):
    config, cases = bench
    target = _json(CONFIG.parent / config["target_configuration"])["target_metadata"]
    matrix = cases["LOCAL_SAMPLE_MISSING"][1]
    genotypes = ["T/T"] * 191 + ["T/A"]
    local = with_target(matrix, target, genotypes)
    result = final_expectation(local, ["A", "B"] * 96, COHORTS[2], config["final_parameters"], target["project_variant_id"])
    assert result.samples[0]["QC_STATUS"] == "EXCLUDED"
    assert result.variants[-1]["QC_STATUS"] == "RETAINED_EXCEPTION"
    assert result.variants[-1]["PRESENT_AFTER_QC"] == "true"
    assert len(result.retained.samples) == 191


def test_all_missing_hwe_is_not_fabricated():
    matrix = GenotypeMatrix((("A", "A"), ("B", "B")), ((1, 1, "v1"), (1, 2, "v2")),
                            np.array([[[-1, -1], [-1, -1]], [[0, 0], [0, 1]]], dtype=np.int8))
    config = _json(CONFIG)
    # Retain subjects to isolate variant-level all-missing semantics.
    config["final_parameters"]["sample_missingness_max"] = .9
    result = final_expectation(matrix, ["A", "B"], COHORTS[0], config["final_parameters"], "target")
    assert result.variants[0]["HWE_STATUS"] == "NOT_EVALUATED_DATA"
    assert "hwe_not_evaluable" in result.variants[0]["FILTER_CODES"]


def test_confusion_keeps_null_rates_and_rejects_unknown_units():
    assert confusion(set(), set(), {"A"}) == dict(tp=0, fp=0, fn=0, tn=1, sensitivity=None, specificity=1.0)
    assert confusion({"A"}, {"B"}, {"A", "B"}) == dict(tp=0, fp=1, fn=1, tn=0, sensitivity=0.0, specificity=0.0)
    with pytest.raises(ValueError, match="unknown_unit"):
        confusion(set(), {"C"}, {"A"})


def test_exact_thresholds_are_retained_and_strict_excess_is_excluded():
    config = _json(CONFIG)
    calls = np.zeros((10, 20, 2), dtype=np.int8)
    calls[0, 0] = -1  # 1/10 sujet et 1/20 variant : exactement les deux seuils.
    calls[1, 1:3] = -1
    matrix = GenotypeMatrix(tuple((str(index), str(index)) for index in range(20)),
                            tuple((1, index + 1, f"v{index}") for index in range(10)), calls)
    result = preliminary_expectation(matrix, ["A", "B"] * 10, config["preliminary_parameters"])
    assert result.samples[0]["QC_STATUS"] != "EXCLUDED"
    assert result.variants[0]["QC_STATUS"] != "EXCLUDED"
    assert result.variants[1]["QC_STATUS"] == "EXCLUDED"
    assert np.isnan(result.samples[0]["HETEROZYGOSITY_Z"])
    assert np.isnan(result.samples[0]["INBREEDING_F"])


def test_final_metrics_are_recomputed_after_removing_subjects():
    config = _json(CONFIG)
    calls = np.zeros((20, 100, 2), dtype=np.int8)
    calls[:3, 0] = -1
    matrix = GenotypeMatrix(tuple((str(index), str(index)) for index in range(100)),
                            tuple((1, index + 1, f"v{index}") for index in range(20)), calls)
    result = final_expectation(matrix, ["A", "B"] * 50, COHORTS[1], config["final_parameters"], "target")
    assert result.samples[0]["QC_STATUS"] == "EXCLUDED"
    assert result.variants[0]["N_GENO"] == 99
    assert result.variants[0]["N_MISS"] == 0
    assert result.variants[0]["F_MISS"] == 0
