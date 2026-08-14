import pytest
from pathlib import Path

from effet_fondateur.founder_enrichment import (
    HaplotypeProfile,
    Marker,
    empirical_probability,
    enumerate_internal_null,
    evaluate_exact_sharing,
    sample_external_null,
    summarize_null,
)
from effet_fondateur.founder_enrichment.observed import (
    EnrichmentAnalysisError, distinct_background_count,
    validate_independent_family_profiles,
)
from effet_fondateur.stages.evaluate_founder_haplotype_enrichment import _classification
from effet_fondateur.contracts import validate_tsv_table
from effet_fondateur.founder_enrichment.publication import NULL_DRAW_COLUMNS, null_draw_rows, write_tsv
from itertools import chain


def _markers() -> tuple[Marker, ...]:
    return (
        Marker("l2", 80, 0.8), Marker("l1", 90, 0.9),
        Marker("target", 100, 1.0, True),
        Marker("r1", 110, 1.1), Marker("r2", 120, 1.25),
    )


def _profile(individual: str, haplotype: str, alleles: tuple[str | None, ...]) -> HaplotypeProfile:
    return HaplotypeProfile(individual, haplotype, alleles)


def test_observed_segment_is_asymmetric_and_target_is_excluded() -> None:
    profiles = (
        _profile("a", "H1", ("0", "1", "1", "0", "1")),
        _profile("b", "H1", ("0", "1", "0", "0", "0")),
        _profile("c", "H2", ("0", "1", None, "0", "1")),
    )
    statistic = evaluate_exact_sharing(_markers(), profiles)
    assert statistic.evaluation_status == "EVALUATED"
    assert statistic.left_marker_count == 2
    assert statistic.right_marker_count == 1
    assert statistic.left_shared_cm == pytest.approx(0.2)
    assert statistic.right_shared_cm == pytest.approx(0.1)


def test_exact_sharing_stops_at_first_missing_or_disagreement() -> None:
    profiles = (
        _profile("a", "H1", ("0", None, "1", "0", "1")),
        _profile("b", "H1", ("0", "1", "1", "0", "0")),
        _profile("c", "H1", ("0", "1", "1", "0", "1")),
    )
    statistic = evaluate_exact_sharing(_markers(), profiles)
    assert statistic.left_marker_count == 0
    assert statistic.right_marker_count == 1


def test_observed_segment_can_be_symmetric() -> None:
    profiles = tuple(_profile(individual, "H1", ("0", "0", "1", "0", "1")) for individual in ("a", "b", "c"))
    statistic = evaluate_exact_sharing(_markers(), profiles)
    assert statistic.left_marker_count == statistic.right_marker_count == 2
    assert statistic.left_shared_bp == statistic.right_shared_bp == 20


def test_internal_enumeration_excludes_two_haplotypes_from_same_individual() -> None:
    haplotypes = tuple(
        _profile(individual, haplotype, ("0", "0", "0", "0", "0"))
        for individual in ("a", "b", "c") for haplotype in ("H1", "H2")
    )
    draws = enumerate_internal_null(_markers(), haplotypes, unit_count=3)
    assert len(draws) == 8
    assert all(len(set(draw.individual_ids)) == 3 for draw in draws)


def test_external_monte_carlo_is_deterministic_and_keeps_non_evaluable_draws() -> None:
    haplotypes = (
        _profile("a", "H1", ("0", "0", "0", "0", "0")),
        _profile("b", "H1", ("0", None, "0", "0", "0")),
        _profile("c", "H1", ("0", "0", "0", "0", "0")),
        _profile("d", "H1", ("0", "0", "0", "0", "0")),
    )
    first = sample_external_null(_markers(), haplotypes, unit_count=3, evaluable_draws=5, max_attempts=50, random_seed=42, minimum_flank_markers=1)
    second = sample_external_null(_markers(), haplotypes, unit_count=3, evaluable_draws=5, max_attempts=50, random_seed=42, minimum_flank_markers=1)
    assert first == second
    assert sum(draw.statistic.evaluation_status == "EVALUATED" for draw in first) == 5
    assert any(draw.statistic.evaluation_status == "NOT_EVALUATED" for draw in first)


def test_empirical_plus_one_and_exact_internal_fraction_are_distinct() -> None:
    assert empirical_probability(0, 9) == 0.1
    haplotypes = tuple(
        _profile(individual, "H1", ("0", "0", "0", "0", "0"))
        for individual in ("a", "b", "c")
    )
    draws = enumerate_internal_null(_markers(), haplotypes, unit_count=3)
    summary = summarize_null(draws, 0.4, exhaustive=True)
    assert summary.exact_probability == 1.0
    assert summary.empirical_probability == 1.0


def test_family_units_are_grouped_once_and_multiple_backgrounds_are_visible() -> None:
    profiles = (
        HaplotypeProfile("a", "H1", ("0", "0", "1"), "family_1"),
        HaplotypeProfile("b", "H1", ("0", "0", "1"), "family_2"),
        HaplotypeProfile("c", "H1", ("0", "1", "1"), "family_3"),
    )
    validate_independent_family_profiles(profiles, 3)
    assert distinct_background_count(profiles, (0, 1)) == 2
    with pytest.raises(EnrichmentAnalysisError, match="not_independent"):
        validate_independent_family_profiles((profiles[0], HaplotypeProfile("d", "H2", ("0", "0", "1"), "family_1"), profiles[2]), 3)


def test_absent_threshold_stays_not_classified() -> None:
    assert _classification(0.0001, None, True, False) == "NOT_CLASSIFIED"
    assert _classification(0.0001, 0.001, True, False) == "UNUSUAL_TARGET_CENTERED_SHARING"
    assert _classification(0.5, 0.001, True, False) == "NO_UNUSUAL_SHARING_DETECTED"
    assert _classification(None, None, False, False) == "NOT_EVALUATED"


def test_synthetic_scientific_flow_runs_end_to_end_without_real_data(tmp_path: Path) -> None:
    markers = _markers()
    observed_profiles = tuple(
        HaplotypeProfile(individual, "H1", ("0", "0", "1", "0", "0"), f"family_{index}")
        for index, individual in enumerate(("carrier_a", "carrier_b", "carrier_c"), 1)
    )
    validate_independent_family_profiles(observed_profiles, 3)
    observed = evaluate_exact_sharing(markers, observed_profiles, minimum_flank_markers=1)
    background = tuple(
        _profile(individual, haplotype, ("0", "0", "0", "0", "0"))
        for individual in ("control_a", "control_b", "control_c", "control_d")
        for haplotype in ("H1", "H2")
    )
    internal = enumerate_internal_null(markers, background, unit_count=3, minimum_flank_markers=1)
    external = sample_external_null(markers, background, unit_count=3, evaluable_draws=10, max_attempts=20, random_seed=7, minimum_flank_markers=1)
    assert observed.total_shared_cm is not None
    assert summarize_null(internal, observed.total_shared_cm, exhaustive=True).evaluable_draws == 32
    assert summarize_null(external, observed.total_shared_cm, exhaustive=False).evaluable_draws == 10
    output = tmp_path / "null_draws.tsv.gz"
    write_tsv(output, NULL_DRAW_COLUMNS, null_draw_rows(chain(internal, external)))
    assert validate_tsv_table(output, "founder_haplotype_null_draws.schema.json").row_count == 42
