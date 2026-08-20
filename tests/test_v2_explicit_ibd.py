from pathlib import Path

import pytest

from effet_fondateur.explicit_ibd import (
    ExplicitIbdError,
    IbdSegment,
    Scenario,
    classify_explicit_ibd,
    evaluate_scenario,
    parse_hap_ibd,
    parse_refined_ibd,
)


FAMILIES = ("F1", "F2", "F3")
MUTANT = {family: frozenset({"H1"}) for family in FAMILIES}


def _segments(scenario="primary", end=300, tools=("HAP_IBD", "REFINED_IBD"), pairs=(("F1", "F2"), ("F1", "F3"), ("F2", "F3")), hap="H1"):
    rows = []
    for tool in tools:
        for index, (left, right) in enumerate(pairs):
            rows.append(IbdSegment(tool, scenario, f"S{left}", hap, left, f"S{right}", hap, right, "19", 100 + index, end + index, 1.0, 3.2, 120))
    return rows


def _evaluate(segments, scenario=None, **changes):
    values = dict(
        scenario=scenario or Scenario("primary", "PRIMARY", 2.0, 100, True),
        segments=segments, families=FAMILIES, mutant_haplotypes=MUTANT,
        target_position_bp=200, boundary_tolerance_bp=10,
        control_frequency=0.01, maximum_control_frequency=0.05,
        calibration_specificity_acceptable=True,
    )
    values.update(changes)
    return evaluate_scenario(**values)


def test_concordant_segment_containing_target_is_primary_support():
    primary = _evaluate(_segments())
    assert primary.all_pairs_by_both_tools
    assert classify_explicit_ibd(primary, ()) == "PRIMARY_CONCORDANT_IBD_SUPPORT"


def test_segment_not_containing_target_is_not_a_primary_call():
    assert classify_explicit_ibd(_evaluate(_segments(end=180)), ()) == "NO_PRIMARY_IBD_CALL"


def test_discordant_boundaries_are_reported():
    segments = _segments()
    segments[-1] = IbdSegment(**{**segments[-1].__dict__, "start_bp": 150})
    assert classify_explicit_ibd(_evaluate(segments), ()) == "METHOD_DISCORDANT"


@pytest.mark.parametrize("segments", [_segments(pairs=(("F1", "F2"), ("F1", "F3"))), _segments(tools=("HAP_IBD",))])
def test_missing_pair_or_single_tool_is_not_a_primary_call(segments):
    assert classify_explicit_ibd(_evaluate(segments), ()) == "NO_PRIMARY_IBD_CALL"


def test_sensitivity_signal_never_becomes_primary_support():
    primary = _evaluate([], scenario=Scenario("primary", "PRIMARY", 2.0, 100, True))
    sensitivity = _evaluate(_segments("s1"), scenario=Scenario("s1", "SENSITIVITY", 1.5, 50, True))
    assert classify_explicit_ibd(primary, (sensitivity,)) == "SENSITIVITY_ONLY_IBD_SUPPORT"


def test_unassignable_mutant_chromosome_blocks_call():
    result = _evaluate(_segments(), mutant_haplotypes={"F1": frozenset(), "F2": frozenset({"H1"}), "F3": frozenset({"H1"})})
    assert classify_explicit_ibd(result, ()) == "NO_PRIMARY_IBD_CALL"


def test_frequent_control_signal_blocks_positive_classification():
    result = _evaluate(_segments(), control_frequency=0.2)
    assert classify_explicit_ibd(result, ()) == "NO_PRIMARY_IBD_CALL"


def test_missing_calibration_is_not_evaluable():
    result = _evaluate(_segments(), calibration_specificity_acceptable=False)
    assert classify_explicit_ibd(result, ()) == "NOT_EVALUABLE"


def test_post_hoc_primary_threshold_is_refused():
    with pytest.raises(ExplicitIbdError, match="post_hoc"):
        _evaluate([], scenario=Scenario("primary", "PRIMARY", 1.5, 18, True))


def test_both_parsers_normalize_same_universe(tmp_path: Path):
    rows = "id1\thap1\tid2\thap2\tchrom\tstart\tend\nS1\t1\tS2\t1\t19\t100\t300\n"
    hap_path, refined_path = tmp_path / "hap.ibd", tmp_path / "refined.ibd"
    hap_path.write_text(rows, encoding="utf-8")
    refined_path.write_text(rows, encoding="utf-8")
    kwargs = dict(scenario_id="primary", family_by_sample={"S1": "F1", "S2": "F2"}, cm_at_bp={100: 1.0, 300: 3.2}, marker_positions=tuple(range(100, 301)))
    hap = parse_hap_ibd(hap_path, **kwargs)[0]
    refined = parse_refined_ibd(refined_path, **kwargs)[0]
    assert hap.tool == "HAP_IBD" and refined.tool == "REFINED_IBD"
    assert hap.marker_count == refined.marker_count == 201


def test_missing_marker_coordinate_is_rejected(tmp_path: Path):
    path = tmp_path / "hap.ibd"
    path.write_text("S1\t1\tS2\t1\t19\t100\t301\n", encoding="utf-8")
    with pytest.raises(ExplicitIbdError, match="coordinate"):
        parse_hap_ibd(path, "primary", {"S1": "F1", "S2": "F2"}, {100: 1.0, 300: 3.2}, (100, 300))


def test_within_family_calls_are_normalized_but_cannot_fill_interfamily_pairs(tmp_path: Path):
    path = tmp_path / "hap.ibd"
    path.write_text("S1\t1\tS2\t1\t19\t100\t300\n", encoding="utf-8")
    segments = parse_hap_ibd(path, "primary", {"S1": "F1", "S2": "F1"}, {100: 1.0, 300: 3.2}, (100, 300))
    assert segments[0].family_pair == ("F1", "F1")
