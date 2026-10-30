"""D22 : vérités, invariances, oracles de référence et décisions de cohortes."""

from copy import deepcopy
from pathlib import Path

import numpy as np
import pytest

from effet_fondateur.calibration.cohort_oracle import cohort_expectation
from effet_fondateur.calibration.materialize_ibd import _json
from effet_fondateur.calibration.structure_bench import CASES, fixture_tables, generate_cases, registry_rows, validate_configuration
from effet_fondateur.calibration.structure_oracle import compare_tables, fit_reference_pca, structure_metrics


CONFIG = Path(__file__).resolve().parents[2] / "calibration/configurations/d22_structure_bench_v1.json"


@pytest.fixture
def config():
    result = _json(CONFIG)
    result["markers_per_autosome"] = 32
    return result


def test_seeded_cases_preserve_reference_and_hide_strata(config):
    panels, truth = generate_cases(config)
    repeated, repeated_truth = generate_cases(config)
    assert truth == repeated_truth and truth["reference_count"] == 108
    assert panels[0].alleles.shape == (704, 128, 2)
    assert len(truth["injected_outlier_indices"]) == 8
    for panel, second in zip(panels, repeated, strict=True):
        assert np.array_equal(panel.alleles, second.alleles)
    assert np.array_equal(panels[0].alleles[:, :108], panels[1].alleles[:, :108])
    assert not np.array_equal(panels[0].alleles[:, 108:], panels[1].alleles[:, 108:])
    assert np.array_equal(panels[0].alleles, panels[2].alleles)
    assert np.all(panels[0].alleles[0, :108] == 0)
    assert np.all(panels[0].alleles[1, :108] == -1)
    registry = registry_rows(panels[0], CASES[0])
    assert not {"A", "B", "MIXED", "RARE"} & {value for row in registry for value in row.values()}
    assert registry != registry_rows(panels[2], CASES[2])


@pytest.mark.parametrize("mutation", ["seed_collision", "unknown_field", "third_outlier_axis", "unapproved", "real_review_role"])
def test_configuration_refuses_changed_prespecification(config, mutation):
    if mutation == "seed_collision":
        config["seed"] = config["missingness_seed"]
    elif mutation == "unknown_field":
        config["extra"] = 1
    elif mutation == "third_outlier_axis":
        config["population_parameters"]["outlier_components"] = 3
    elif mutation == "unapproved":
        config["manual_approval"]["approved"] = False
    else:
        config["fixture_review_policy"]["reviewer_role"] = "scientific_reviewer"
    with pytest.raises(ValueError):
        validate_configuration(config)


def test_gram_pca_agrees_with_separate_svd_and_ignores_projected_calls(config):
    rng = np.random.default_rng(7654)
    dosage = rng.binomial(2, .4, (20, 180)).astype(float)
    dosage[1, 8] = np.nan
    reference = np.arange(20) < 15
    parameters = dict(config["population_parameters"], min_informative_variants=10)
    oracle = fit_reference_pca(dosage, reference, parameters)
    frequency = np.nanmean(dosage[reference], axis=0) / 2
    standardized = np.nan_to_num((dosage - 2 * frequency) / np.sqrt(2 * frequency * (1 - frequency)), nan=0) / np.sqrt(180)
    _, singular, weights = np.linalg.svd(standardized[reference], full_matrices=False)
    weights = weights[:3].T
    for axis in range(3):
        if weights[np.argmax(np.abs(weights[:, axis])), axis] < 0:
            weights[:, axis] *= -1
    assert np.allclose(oracle["scores"], standardized @ weights, atol=1e-12)
    assert np.allclose(oracle["eigenvalues"], singular[:3] ** 2 / 14, atol=1e-12)
    assert np.allclose(oracle["p_values"], np.exp(-oracle["distances"] / 2))
    dosage[15:] = 2 - dosage[15:]
    perturbed = fit_reference_pca(dosage, reference, parameters)
    assert np.array_equal(oracle["weights"], perturbed["weights"])
    assert np.array_equal(oracle["scores"][:15], perturbed["scores"][:15])


def test_reference_specific_exclusions_and_no_rank_fabrication(config):
    dosage = np.random.default_rng(123).binomial(2, .4, (10, 30)).astype(float)
    reference = np.arange(10) < 8
    dosage[:8, 0] = 0
    dosage[:8, 1] = np.nan
    parameters = dict(config["population_parameters"], min_informative_variants=10)
    oracle = fit_reference_pca(dosage, reference, parameters)
    assert oracle["included"][:2].tolist() == [False, False]
    assert oracle["call_rate_eligible"][:2].tolist() == [True, False]
    with pytest.raises(ValueError, match="insufficient_informative"):
        fit_reference_pca(np.zeros((10, 30)), reference, parameters)


@pytest.mark.parametrize("field,value", [("VALUE", "0"), ("FLAG", "true"), ("KEY", "unknown")])
def test_numeric_or_decision_mutations_cannot_pass(config, field, value):
    expected = [dict(KEY="known", VALUE=1e-30, FLAG="false")]
    actual = [dict(KEY="known", VALUE="1e-30", FLAG="false")]
    strict = compare_tables(expected, actual, "KEY", config)
    assert all(count == 0 for key, count in strict.items() if key != "compared")
    actual[0][field] = value
    # P_VALUE interdit d'assimiler une probabilité minuscule à zéro.
    if field == "VALUE":
        expected[0]["P_VALUE"] = expected[0].pop("VALUE")
        actual[0]["P_VALUE"] = actual[0].pop("VALUE")
    changed = compare_tables(expected, actual, "KEY", config)
    assert any(count for key, count in changed.items() if key != "compared")


def test_cohort_roles_depend_on_explicit_genotype_and_reviews_not_stratum(config):
    panels, truth = generate_cases(config)
    registry = registry_rows(panels[0], CASES[0])
    target = dict(ref="T", alt="A")
    tables = fixture_tables(panels[0], registry, config, target)
    scores = [dict(SAMPLE_ID=sample, REFERENCE_INCLUDED=str(index < 108).lower()) for index, sample in enumerate(panels[0].samples)]
    outliers = [dict(SAMPLE_ID=sample, PROPOSED_EXCLUSION=str(index == 104).lower()) for index, sample in enumerate(panels[0].samples)]
    review = dict(review_id="fixture", reviewer_role="calibration_contract_fixture", reviewed_at="2026-10-30T07:00:00Z", reviewed_artifact_sha256="a" * 64)
    parameters = dict(control_group_labels=["CONTROL"], required_nonempty_cohorts=[], decision_reviews={key: review for key in ("kinship_exclusion_approval", "population_structure_exclusion_approval")})
    result = cohort_expectation(registry, tables["target_genotype_audit"], tables["qc_individual_metrics"], tables["independent_set_proposal"], scores, outliers, parameters, "A")
    assert len(result["rows"]) == 7 * 128
    carrier = registry[2]["SAMPLE_ID"]  # CONTROL mais A/A explicite : jamais témoin.
    row = next(row for row in result["rows"] if row["COHORT_ID"] == "controls_unrelated" and row["SAMPLE_ID"] == carrier)
    assert row["EXCLUSION_CODE"] == "target_variant_carrier"
    all_carriers = [row for row in result["rows"] if row["COHORT_ID"] == "target_carriers_all" and row["INCLUDED"] == "true"]
    assert len(all_carriers) == 85
    assert all("calibration_contract_fixture" == row["REVIEWER_ROLE"] for row in result["decisions"])
    altered = deepcopy(tables["target_genotype_audit"])
    altered[2]["GENOTYPE"] = "T/T"
    changed = cohort_expectation(registry, altered, tables["qc_individual_metrics"], tables["independent_set_proposal"], scores, outliers, parameters, "A")
    assert result["rows"] != changed["rows"]


def test_detection_intervals_exclude_projected_copy_pseudoreplication(config):
    panels, truth = generate_cases(config)
    dosage = np.where(panels[0].alleles[:, :, 0] < 0, np.nan, panels[0].alleles.sum(axis=2)).T
    reference = np.arange(128) < 108
    model = fit_reference_pca(dosage, reference, config["population_parameters"])
    metrics = structure_metrics(model["scores"], truth, model["outliers"], .95)
    independent = metrics["independent_reference_detection"]
    assert independent["tp"] + independent["fn"] == 4
    assert independent["sensitivity_wilson"][0] < 1
    assert "sensitivity_wilson" not in metrics["all_samples_detection_without_independent_interval"]
    assert metrics["unsupported_geographic_origin_claim_count"] == 0
