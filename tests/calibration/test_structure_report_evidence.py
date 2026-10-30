"""Descriptions graphiques agrégées, sans identifiant ou nouvelle analyse."""

import json

import pytest

from effet_fondateur.calibration.structure_bench import COHORTS
from effet_fondateur.calibration.structure_report_evidence import aggregate_cohorts, aggregate_population


def test_population_descriptions_preserve_counts_and_never_export_individuals():
    keys = ["PRIVATE_SAMPLE_A", "PRIVATE_SAMPLE_B", "PRIVATE_SAMPLE_C", "PRIVATE_SAMPLE_D"]
    scores = [dict(SAMPLE_ID=sample, PC1=str(index), PC2=str(2 * index), REFERENCE_INCLUDED=str(index < 3).lower()) for index, sample in enumerate(keys)]
    outliers = [dict(SAMPLE_ID=sample, DISTANCE_SQUARED=str(index * 10)) for index, sample in enumerate(keys)]
    result = aggregate_population(scores, outliers, dict(groups=["A", "A", "B", "RARE"], reference_count=3), .01)
    assert sum(row["count"] for row in result["centroids"]) == 4
    assert result["centroids"][0]["mean"] == [.5, 1.]
    assert result["centroids"][0]["standard_deviation"][0] == pytest.approx(2 ** -.5)
    histogram = result["distance_histogram"]
    assert sum(histogram["reference_counts"]) == 3 and sum(histogram["projected_counts"]) == 1
    assert histogram["threshold"] == pytest.approx(9.210340371976182)
    assert not any(key in json.dumps(result) for key in keys)


def test_visual_population_axis_errors_are_not_silently_aggregated():
    with pytest.raises(ValueError, match="sample_axis"):
        aggregate_population([dict(SAMPLE_ID="A")], [dict(SAMPLE_ID="B")], dict(groups=["A"]), .01)


def test_cohort_complement_is_all_subjects_not_only_candidates():
    rows = [dict(COHORT_ID=cohort, SAMPLE_ID=sample, INCLUDED=str(index == 0).lower(), EXCLUSION_CODE="" if index == 0 else "not_explicit_control_group")
            for cohort in COHORTS for index, sample in enumerate(("PRIVATE_A", "PRIVATE_B"))]
    result = aggregate_cohorts(rows, 2)
    assert result["controls_unrelated"] == dict(included=1, not_included=1, exclusion_counts=dict(not_explicit_control_group=1))
    assert "PRIVATE" not in json.dumps(result)
    rows[-1]["EXCLUSION_CODE"] = ""
    with pytest.raises(ValueError, match="partition_invalid"):
        aggregate_cohorts(rows, 2)


def test_missing_cohort_members_cannot_be_pictured_as_zero():
    with pytest.raises(ValueError, match="cohort_axis"):
        aggregate_cohorts([], 128)
