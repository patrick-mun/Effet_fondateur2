"""D21 : références indépendantes, cas limites et absence de processus externe."""

from copy import deepcopy
from fractions import Fraction
from pathlib import Path

import numpy as np
import pytest

from effet_fondateur.calibration.kinship_bench import validate_configuration, profile_panel
from effet_fondateur.calibration.kinship_oracle import classification_metrics, degree_and_relation, independent_proposal, king_metrics, wilson_interval
from effet_fondateur.calibration.materialize_ibd import _json
from effet_fondateur.calibration.simulation.kinship_bench import DEGREES, EXPECTED_PHI, RELATIONS, generate_simulation, pedigree_phi, pedigree_template
from effet_fondateur.calibration.technical_evaluation import GenotypeMatrix


CONFIG = Path(__file__).resolve().parents[2] / "calibration/configurations/d21_kinship_bench_v1.json"


@pytest.fixture
def small_config():
    config = _json(CONFIG)
    config.update(markers_per_autosome=32, families_per_relation_per_stratum=1, bootstrap_replicates=100)
    return config


@pytest.mark.parametrize("relation", RELATIONS)
def test_exact_pedigree_coefficients_cover_all_degrees(relation):
    parents, observed, clone = pedigree_template(relation)
    assert pedigree_phi(parents, observed, clone) == EXPECTED_PHI[relation]
    assert EXPECTED_PHI["COUSINS_ONCE_REMOVED"] == Fraction(1, 32)


def test_invalid_pedigree_order_is_rejected():
    with pytest.raises(ValueError, match="order_invalid"):
        pedigree_phi([None, (0, 1)], (0, 1))


def test_gene_drop_is_deterministic_mendelian_and_separated(small_config):
    first, second = generate_simulation(small_config), generate_simulation(small_config)
    assert first.alleles.shape == (704, 28, 2)
    assert np.array_equal(first.alleles, second.alleles)
    assert np.array_equal(first.origins, second.origins)
    assert first.mendelian_error_count == 0
    assert set(unit["degree"] for unit in first.units) == set(DEGREES)
    assert len(first.units) == 14
    for unit in first.units:
        left, right = unit["sample_indices"]
        realized = np.mean(np.sum(first.origins[:, left, :, None] == first.origins[:, right, None, :], axis=(1, 2)) / 4)
        assert abs(realized - unit["theoretical_phi"]) < .05
        if unit["relation"] == "DUPLICATE":
            assert np.array_equal(first.alleles[:, left], first.alleles[:, right])
    for offset in range(22):
        assert np.array_equal(first.alleles[offset * 32], first.alleles[offset * 32 + 1])
    assert np.count_nonzero(first.alleles[2]) == 0
    first_unit = set(first.origins[:, :2].ravel())
    assert not first_unit & set(first.origins[:, 2:4].ravel())


def test_profiles_preserve_truth_and_have_distinct_pseudonyms(small_config):
    simulation = generate_simulation(small_config)
    complete, selected = profile_panel(simulation, small_config, 0)
    missing, _ = profile_panel(simulation, small_config, 1)
    sparse, subset = profile_panel(simulation, small_config, 2)
    assert len(complete.variants) == 704 and len(sparse.variants) == 176
    assert np.array_equal(complete.alleles[subset], sparse.alleles)
    assert not set(complete.samples) & set(missing.samples)
    called = missing.alleles[:, :, 0] >= 0
    assert np.array_equal(complete.alleles[called], missing.alleles[called])
    assert .03 < (~called).mean() < .07
    assert len(selected) == 704


@pytest.mark.parametrize("mutation", ["unknown_key", "seed_collision", "profile_collision", "overlap", "missing_autosome", "zero_spacing", "profile_semantics", "degree_scope"])
def test_configuration_refuses_uncontrolled_bench(small_config, mutation):
    config = deepcopy(small_config)
    if mutation == "unknown_key":
        config["unexpected"] = 1
    elif mutation == "seed_collision":
        config["bootstrap_seed"] = config["simulation_seed"]
    elif mutation == "profile_collision":
        config["profiles"][1]["profile_id"] = config["profiles"][0]["profile_id"]
    elif mutation == "overlap":
        config["monomorphic_offset"] = 0
    elif mutation == "missing_autosome":
        config["panel_parameters"]["required_autosomes"] = list(range(1, 22))
    elif mutation == "zero_spacing":
        config["marker_spacing_bp"] = 0
    elif mutation == "profile_semantics":
        config["profiles"][1]["missingness_fraction"] = 0
    else:
        config["kinship_parameters"]["relatedness_max_degree_for_independence"] = 3
    with pytest.raises(ValueError):
        validate_configuration(config)


@pytest.mark.parametrize("same_fid,expected", [(True, .2), (False, .125)])
def test_king_equations_have_independent_hand_calculated_oracles(same_fid, expected):
    calls = np.array([[[0, 0], [0, 2]], [[0, 2], [2, 2]], [[2, 2], [0, 2]], [[0, 2], [0, 2]]], dtype=np.int8)
    matrix = GenotypeMatrix((("F", "A"), ("F" if same_fid else "G", "B")), tuple((1, index + 1, str(index)) for index in range(4)), calls)
    metrics = king_metrics(matrix, 0, 1)
    assert metrics == dict(N_SNP=4, HETHET=.25, IBS0=0, KINSHIP=pytest.approx(expected))
    calls[0] = -1
    assert king_metrics(matrix, 0, 1)["N_SNP"] == 3


def test_zero_heterozygosity_is_unmeasured_not_unrelated():
    matrix = GenotypeMatrix((("F", "A"), ("F", "B")), ((1, 1, "v"),), np.zeros((1, 2, 2), dtype=np.int8))
    assert king_metrics(matrix, 0, 1)["KINSHIP"] is None


@pytest.mark.parametrize("value,degree", [(.354, "DUPLICATE_OR_MZ"), (.177, "FIRST"), (.0884, "SECOND"), (.0442, "THIRD"), (.0221, "FOURTH"), (-.1, "UNRELATED")])
def test_classification_boundaries_are_inclusive(value, degree):
    parameters = _json(CONFIG)["kinship_parameters"]
    assert degree_and_relation(value, 0, parameters)[0] == degree
    assert degree_and_relation(.25, .1, parameters)[1] == "FIRST_DEGREE_NON_PARENT_CHILD"


def test_proposal_is_maximal_not_one_person_per_component():
    rows = [dict(SAMPLE_ID=sample, QC_STATUS="PASS", F_MISS="0", HETEROZYGOSITY_Z="0") for sample in ("A", "B", "C")]
    result, _ = independent_proposal(["A", "B", "C"], {("A", "B"): {"KING"}, ("B", "C"): {"KING"}}, rows)
    assert result["selected"] == {"A", "C"}
    assert result["representatives"] == {"B": "A"}


def test_independent_priority_prefers_quality_not_identifier():
    rows = [dict(SAMPLE_ID="A", QC_STATUS="ALERT", F_MISS="0", HETEROZYGOSITY_Z="0"),
            dict(SAMPLE_ID="B", QC_STATUS="PASS", F_MISS=".01", HETEROZYGOSITY_Z=".1")]
    result, _ = independent_proposal(["A", "B"], {("A", "B"): {"KING"}}, rows)
    assert result["selected"] == {"B"}


def test_intervals_do_not_turn_zero_observed_errors_into_certainty(small_config):
    interval = wilson_interval(8, 8, .95)
    assert 0 < interval[0] < .8 and interval[1] == pytest.approx(1)
    assert wilson_interval(0, 0, .95) is None
    records = [dict(true_degree=degree, observed_degree=degree, stratum_id="A") for degree in DEGREES]
    records[2]["observed_degree"] = "MISSING"
    result = classification_metrics(records, small_config)
    assert result["confusion"]["SECOND"]["MISSING"] == 1
    assert result["balanced_accuracy"] == pytest.approx(5 / 6)
    assert result["missing_pair_count"] == 1
    assert result["qualification_minimum_per_class_met"] is False
