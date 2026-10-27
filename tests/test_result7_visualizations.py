from __future__ import annotations

import pytest

from scripts.plot_result7_roh_age import (
    Result7VisualizationError,
    _common_overlap,
    _quantiles,
)


def test_common_roh_overlap_must_contain_target() -> None:
    rows = [
        {"TARGET_IN_ROH": "true", "ROH_POS1_BP": "9000000", "ROH_POS2_BP": "12000000"},
        {"TARGET_IN_ROH": "true", "ROH_POS1_BP": "10500000", "ROH_POS2_BP": "13000000"},
    ]
    assert _common_overlap(rows) == (10_500_000, 12_000_000)


def test_common_roh_overlap_rejects_disjoint_target_segments() -> None:
    rows = [
        {"TARGET_IN_ROH": "true", "ROH_POS1_BP": "9000000", "ROH_POS2_BP": "11000000"},
        {"TARGET_IN_ROH": "true", "ROH_POS1_BP": "11500000", "ROH_POS2_BP": "13000000"},
    ]
    with pytest.raises(Result7VisualizationError, match="overlap_missing"):
        _common_overlap(rows)


def test_quantiles_are_stable_for_report_values() -> None:
    assert _quantiles([1, 2, 3, 4, 5]) == {
        "minimum": 1.0, "q25": 2.0, "median": 3.0, "q75": 4.0, "maximum": 5.0,
    }
