import csv
import json
from pathlib import Path

import pytest

from effet_fondateur.contracts import validate_tsv_table
from effet_fondateur.dating import (
    GammaAgeError,
    chance_sharing_correction_cm,
    estimate_correlated_gamma,
    estimate_independent_gamma,
    publish_variant_age,
)


def _write_tsv(path: Path, columns: list[str], rows: list[list[str]]) -> None:
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(columns)
        writer.writerows(rows)


def test_gamma_matches_official_six_individual_example() -> None:
    left = [1.2, 3.4, 1.6, 4.8, 4.8, 2.3]
    right = [2.0, 3.8, 2.7, 2.1, 3.8, 1.8]
    correction = chance_sharing_correction_cm(0.7, 68_378, 290)

    independent = estimate_independent_gamma(
        left, right, chance_sharing_correction_cm=correction
    )
    correlated = estimate_correlated_gamma(
        left, right, chance_sharing_correction_cm=correction
    )

    assert correction == pytest.approx(0.0358548778)
    assert independent.generations == pytest.approx(27.5, abs=0.1)
    assert independent.confidence_lower == pytest.approx(15.5, abs=0.1)
    assert independent.confidence_upper == pytest.approx(49.2, abs=0.1)
    assert correlated.generations == pytest.approx(22.3, abs=0.1)
    assert correlated.confidence_lower == pytest.approx(7.9, abs=0.1)
    assert correlated.confidence_upper == pytest.approx(67.5, abs=0.1)


def _founder_outputs(tmp_path: Path, unit_count: int) -> tuple[Path, Path, Path]:
    segment_path = tmp_path / "founder_segments.tsv"
    segment_columns = [
        "INDEPENDENT_UNIT_ID", "SAMPLE_ID", "FAMILY_ID", "CARRIER_HAPLOTYPE_ID",
        "TARGET_VARIANT_ID", "LEFT_BOUND_BP", "TARGET_BP", "RIGHT_BOUND_BP",
        "LEFT_LENGTH_CM", "RIGHT_LENGTH_CM", "PHASING_CONFIDENCE",
        "SEGMENT_METHOD", "SEGMENT_STATUS", "EXCLUSION_CODE",
    ]
    left = [1.2, 3.4, 1.6, 4.8, 4.8, 2.3]
    right = [2.0, 3.8, 2.7, 2.1, 3.8, 1.8]
    _write_tsv(segment_path, segment_columns, [
        [
            f"unit_{index}", f"sample_{index}", f"family_{index}", "H1", "target",
            "50", "100", "150", str(left[index - 1]), str(right[index - 1]),
            "0.99", "target_centered_pairwise_max_ibs_v1", "INCLUDED", "",
        ]
        for index in range(1, unit_count + 1)
    ])
    consensus_path = tmp_path / "founder_consensus.tsv"
    consensus_columns = [
        "ANALYSIS_ID", "TARGET_VARIANT_ID", "INTERPRETATION", "CONSENSUS_STATUS",
        "CARRIER_UNIT_COUNT", "LEFT_BOUND_BP", "TARGET_BP", "RIGHT_BOUND_BP",
        "LEFT_LENGTH_CM", "RIGHT_LENGTH_CM", "LEFT_MARKER_COUNT",
        "RIGHT_MARKER_COUNT", "SIGNATURE_MARKER_COUNT", "BACKGROUND_COHORT",
        "BACKGROUND_HAPLOTYPE_COUNT", "MATCHING_BACKGROUND_HAPLOTYPE_COUNT",
        "BACKGROUND_MATCH_FREQUENCY", "SEGMENT_METHOD", "MAP_VERSION",
    ]
    _write_tsv(consensus_path, consensus_columns, [[
        "primary_exact_ibs", "target", "IBS_SHARED_CANDIDATE",
        "SUPPORTED_IBS_CANDIDATE", str(unit_count), "80", "100", "120", "0.2",
        "0.2", "2", "2", "4", "controls_unrelated_plus_noncarriers", "20", "1",
        "0.05", "target_centered_exact_ibs_v1", "a" * 64,
    ]])
    summary_path = tmp_path / "founder_analysis_summary.json"
    summary_path.write_text(json.dumps({
        "schema_version": "1.0.0", "method_id": "target_centered_exact_ibs_v1",
        "interpretation": "IBS_SHARED_CANDIDATE", "status": "SUPPORTED_IBS_CANDIDATE",
        "selected_carrier_count": unit_count, "excluded_carrier_count": 0,
        "background_haplotype_count": 20, "matching_background_haplotype_count": 1,
        "minimum_independent_carriers": 3, "minimum_flank_markers": 2,
        "ibd_claimed": False,
    }), encoding="utf-8")
    return segment_path, consensus_path, summary_path


def test_variant_age_publishes_primary_and_sensitivity_scenarios(tmp_path: Path) -> None:
    segments, consensus, summary = _founder_outputs(tmp_path, 6)

    publication = publish_variant_age(
        founder_segments_path=segments, founder_consensus_path=consensus,
        founder_summary_path=summary, output_dir=tmp_path / "dating",
        confidence_level=0.95, minimum_units_for_estimate=3,
        minimum_units_for_primary=5, generation_years=[25, 28, 30],
        leave_one_family_out=True,
        chance_sharing_parameters={
            "median_allele_frequency": 0.7,
            "markers_on_chromosome": 68_378,
            "chromosome_length_cm": 290,
        },
    )

    assert publication.status == "PRIMARY_ESTIMATE"
    inputs = validate_tsv_table(
        publication.input_path, "variant_age_input.schema.json"
    )
    assert inputs.row_count == 6
    estimates = validate_tsv_table(
        publication.estimates_path, "variant_age_estimates.schema.json"
    ).rows
    correlated = next(row for row in estimates if row["MODEL"] == "CORRELATED")
    assert correlated["ANALYSIS_STATUS"] == "PRIMARY"
    assert float(correlated["ESTIMATE_GENERATIONS"]) == pytest.approx(22.3, abs=0.1)
    scenarios = validate_tsv_table(
        publication.scenarios_path, "variant_age_scenarios.schema.json"
    ).rows
    assert len(scenarios) == 11
    assert sum(row["SCENARIO_TYPE"] == "LEAVE_ONE_FAMILY_OUT" for row in scenarios) == 6
    years_28 = next(row for row in scenarios if row["SCENARIO_ID"] == "generation_28")
    assert float(years_28["ESTIMATE_YEARS"]) == pytest.approx(22.3396753271 * 28)


def test_variant_age_refuses_estimate_below_three_units(tmp_path: Path) -> None:
    segments, consensus, summary = _founder_outputs(tmp_path, 2)

    publication = publish_variant_age(
        founder_segments_path=segments, founder_consensus_path=consensus,
        founder_summary_path=summary, output_dir=tmp_path / "dating",
        confidence_level=0.95, minimum_units_for_estimate=3,
        minimum_units_for_primary=5, generation_years=[25, 28, 30],
        leave_one_family_out=True, chance_sharing_parameters=None,
    )

    assert publication.status == "NOT_ESTIMATED"
    estimates = validate_tsv_table(
        publication.estimates_path, "variant_age_estimates.schema.json"
    ).rows
    assert all(row["ANALYSIS_STATUS"] == "NOT_ESTIMATED" for row in estimates)
    assert all(row["ESTIMATE_GENERATIONS"] is None for row in estimates)


def test_variant_age_rejects_two_units_from_same_family(tmp_path: Path) -> None:
    segments, consensus, summary = _founder_outputs(tmp_path, 3)
    contents = segments.read_text(encoding="utf-8")
    segments.write_text(
        contents.replace("unit_2\tsample_2\tfamily_2", "unit_2\tsample_2\tfamily_1"),
        encoding="utf-8",
    )

    with pytest.raises(GammaAgeError, match="families_not_independent"):
        publish_variant_age(
            founder_segments_path=segments, founder_consensus_path=consensus,
            founder_summary_path=summary, output_dir=tmp_path / "dating",
            confidence_level=0.95, minimum_units_for_estimate=3,
            minimum_units_for_primary=5, generation_years=[25, 28, 30],
            leave_one_family_out=True, chance_sharing_parameters=None,
        )
