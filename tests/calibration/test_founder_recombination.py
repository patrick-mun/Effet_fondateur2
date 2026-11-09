"""Tests unitaires du fondateur et des recombinaisons D14."""

from __future__ import annotations

import io

import pytest

from effet_fondateur.calibration.simulation.founder_recombination import (
    FounderRecombinationError,
    FounderSimulationConfiguration,
    GeneticMapPoint,
    HaplotypeOrigin,
    PedigreeMember,
    build_founder_plan,
    interpolate_bp,
    interpolate_cm,
    parse_genetic_map,
    transform_founder_vcf,
)


CONFIGURATION = FounderSimulationConfiguration(
    "d14_test",
    1004,
    "D14_TEST_R001",
    150,
    "SYNTHETIC_TARGET_D14",
    19,
    200,
    "C",
    "T",
)
MAP = (
    GeneticMapPoint(100, 0.0),
    GeneticMapPoint(200, 50.0),
    GeneticMapPoint(300, 100.0),
)
ORIGINS = {
    "S1": (HaplotypeOrigin("D1", 1), HaplotypeOrigin("D1", 2)),
    "S2": (HaplotypeOrigin("D1", 1), HaplotypeOrigin("D2", 2)),
    "S3": (HaplotypeOrigin("D3", 1), HaplotypeOrigin("D3", 2)),
    "S4": (HaplotypeOrigin("D3", 1), HaplotypeOrigin("D4", 2)),
    "S5": (HaplotypeOrigin("D5", 1), HaplotypeOrigin("D5", 2)),
}
PEDIGREE = (
    PedigreeMember("S1", "FAM_A", 0, True),
    PedigreeMember("S2", "FAM_A", 2, False),
    PedigreeMember("S3", "FAM_B", 0, True),
    PedigreeMember("S4", "FAM_B", 2, False),
    PedigreeMember("S5", "UNRELATED_S5", 0, True),
)
STRATA = {
    origin: "POP_A" for copies in ORIGINS.values() for origin in copies
}


def _plan():
    return build_founder_plan(
        configuration=CONFIGURATION,
        genetic_map=MAP,
        pedigree=PEDIGREE,
        sample_origins=ORIGINS,
        background_strata=STRATA,
    )


def test_map_collapses_consistent_duplicate_positions() -> None:
    points = parse_genetic_map(
        io.StringIO(
            "CHROMOSOME\tPOSITION_BP\tPOSITION_CM\n"
            "19\t100\t0\n"
            "19\t100\t0\n"
            "19\t200\t1\n"
        ),
        19,
    )
    assert points == (GeneticMapPoint(100, 0.0), GeneticMapPoint(200, 1.0))


def test_map_interpolation_is_reversible_at_integer_bp() -> None:
    assert interpolate_cm(MAP, 150) == 25.0
    assert interpolate_bp(MAP, 25.0) == 150


def test_founder_plan_is_deterministic_and_has_internal_breakpoints() -> None:
    plan = _plan()
    assert plan == _plan()
    assert len(plan.segments) == 2
    assert len({segment.carrier_origin for segment in plan.segments}) == 2
    assert all(not segment.left_censored for segment in plan.segments)
    assert all(not segment.right_censored for segment in plan.segments)
    assert all(segment.recombination_event_count == 2 for segment in plan.segments)
    assert all(
        segment.left_length_cm + segment.right_length_cm
        == pytest.approx(segment.total_length_cm)
        for segment in plan.segments
    )


def test_short_map_records_edge_censoring_without_fake_event() -> None:
    short_map = (
        GeneticMapPoint(190, 49.9),
        GeneticMapPoint(200, 50.0),
        GeneticMapPoint(210, 50.1),
    )
    plan = build_founder_plan(
        configuration=CONFIGURATION,
        genetic_map=short_map,
        pedigree=PEDIGREE,
        sample_origins=ORIGINS,
        background_strata=STRATA,
    )
    assert all(segment.left_censored for segment in plan.segments)
    assert all(segment.right_censored for segment in plan.segments)
    assert all(segment.recombination_event_count == 0 for segment in plan.segments)


def test_vcf_mosaic_preserves_transmission_and_resets_target() -> None:
    source = io.StringIO(
        "##fileformat=VCFv4.2\n"
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t"
        "S1\tS2\tS3\tS4\tS5\n"
        "chr19\t100\tv1\tA\tG\t.\tPASS\t.\tGT\t"
        "0|1\t0|0\t1|0\t1|1\t1|0\n"
        "chr19\t200\tv2\tC\tT\t.\tPASS\t.\tGT\t"
        "0|0\t0|0\t0|0\t0|0\t0|0\n"
        "chr19\t300\tv3\tG\tA\t.\tPASS\t.\tGT\t"
        "0|1\t0|0\t1|0\t1|1\t1|0\n"
    )
    dense = io.StringIO()
    observed = io.StringIO()
    stats = transform_founder_vcf(
        source,
        dense,
        observed,
        plan=_plan(),
        observed_variant_keys={("chr19", 100, "A", "G"): "CHIP-1"},
        transmissions=(("S2", "S1", 1, 1), ("S4", "S3", 1, 1)),
    )
    assert stats.target_variant_count == 1
    assert stats.mendelian_error_count == 0
    assert stats.observed_variant_count == 1
    target = next(
        line for line in dense.getvalue().splitlines() if line.startswith("chr19\t200")
    )
    assert "\tSYNTHETIC_TARGET_D14\t" in target
    assert "SYNTHETIC_TARGET_D14" not in observed.getvalue()


def test_map_rejects_inconsistent_duplicate_position() -> None:
    with pytest.raises(
        FounderRecombinationError,
        match="founder_genetic_map_duplicate_inconsistent",
    ):
        parse_genetic_map(
            io.StringIO(
                "CHROMOSOME\tPOSITION_BP\tPOSITION_CM\n"
                "19\t100\t0\n"
                "19\t100\t1\n"
                "19\t200\t2\n"
            ),
            19,
        )
