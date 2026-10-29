"""Tests unitaires du dimensionnement génétique D14b."""

from __future__ import annotations

import io
import math

import pytest

from effet_fondateur.calibration.simulation.windows import (
    SimulationWindowError,
    WindowSpecification,
    WindowTarget,
    design_window,
    interpolate_window_bp,
    interpolate_window_cm,
    minimum_half_width_cm,
    parse_public_genetic_map,
)


def _map():
    return parse_public_genetic_map(
        io.StringIO(
            "MAP_ID\tASSEMBLY\tCHROMOSOME\tPOSITION_BP\tPOSITION_CM\n"
            "map\tGRCh38\t5\t100\t0\n"
            "map\tGRCh38\t5\t600\t5\n"
            "map\tGRCh38\t5\t1100\t10\n"
        ),
        expected_chromosome=5,
    )


def _spec(purpose: str, challenge_width: float | None = None):
    return WindowSpecification(
        window_id=f"test_{purpose.casefold()}",
        purpose=purpose,
        age_generations=100,
        chromosome=5,
        target=WindowTarget("SYNTHETIC_TARGET", 600, "C", "T"),
        challenge_half_width_cm=challenge_width,
    )


@pytest.mark.parametrize("age", [5, 10, 15, 20, 30, 40, 50, 75, 100, 150])
def test_minimum_half_width_matches_prespecified_formula(age: int) -> None:
    width = minimum_half_width_cm(age, 0.05)
    assert width == pytest.approx(299.5732273553991 / age)
    assert math.exp(-age * width / 100.0) == pytest.approx(0.05)


def test_primary_window_is_rounded_outward_and_within_probability() -> None:
    designed = design_window(
        _spec("PRIMARY"),
        _map(),
        maximum_arm_censor_probability=0.05,
    )
    assert designed.start_bp < 600 < designed.end_bp
    assert designed.realized_left_width_cm >= designed.requested_half_width_cm
    assert designed.realized_right_width_cm >= designed.requested_half_width_cm
    assert designed.theoretical_arm_censor_probability == pytest.approx(0.05)


def test_challenge_window_is_shorter_and_has_high_censor_probability() -> None:
    designed = design_window(
        _spec("CENSORING_CHALLENGE", 0.1),
        _map(),
        maximum_arm_censor_probability=0.05,
    )
    assert designed.requested_half_width_cm == 0.1
    assert designed.theoretical_arm_censor_probability == pytest.approx(
        math.exp(-0.1)
    )


def test_inverse_interpolation_uses_requested_outer_boundary() -> None:
    points = _map()
    lower = interpolate_window_bp(points, 2.995732273553991, boundary="LOWER")
    upper = interpolate_window_bp(points, 7.004267726446009, boundary="UPPER")
    assert interpolate_window_cm(points, lower) <= 2.995732273553991
    assert interpolate_window_cm(points, upper) >= 7.004267726446009


def test_infeasible_window_is_refused_without_clamping() -> None:
    with pytest.raises(
        SimulationWindowError, match="simulation_window_infeasible_window"
    ):
        design_window(
            WindowSpecification(
                "too_old_at_edge",
                "PRIMARY",
                10,
                5,
                WindowTarget("TARGET", 600, "C", "T"),
                None,
            ),
            _map(),
            maximum_arm_censor_probability=0.05,
        )


def test_non_monotonic_map_is_refused() -> None:
    with pytest.raises(
        SimulationWindowError, match="simulation_window_map_not_monotonic"
    ):
        parse_public_genetic_map(
            io.StringIO(
                "ASSEMBLY\tCHROMOSOME\tPOSITION_BP\tPOSITION_CM\n"
                "GRCh38\t5\t100\t1\n"
                "GRCh38\t5\t200\t0\n"
            ),
            expected_chromosome=5,
        )
