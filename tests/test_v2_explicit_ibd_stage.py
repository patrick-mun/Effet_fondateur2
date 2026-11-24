import subprocess

import pytest

from effet_fondateur.explicit_ibd.execution import ExplicitIbdExternalError
from effet_fondateur.external_tools import ExternalCommandError
from effet_fondateur.stages import call_explicit_ibd
from effet_fondateur.stages.call_explicit_ibd import (
    _audit_checks,
    _calibration_status,
    _empirical_control_status,
    _query_complete_variants,
    _scenario,
)


def test_missing_or_unphased_genotypes_are_excluded_from_hap_ibd(monkeypatch, tmp_path):
    output = "v1\t19\t100\t0|0\t0|1\nv2\t19\t200\t0|.\t0|0\nv3\t19\t300\t0/1\t0|0\n"
    monkeypatch.setattr(
        call_explicit_ibd,
        "default_command_runner",
        lambda command, timeout: subprocess.CompletedProcess(command, 0, output, ""),
    )
    rows = _query_complete_variants("bcftools", tmp_path / "input.bcf", 10)
    assert [row[3] for row in rows] == [True, False, False]


@pytest.mark.parametrize(
    ("technical_error", "expected_message"),
    [
        (ExternalCommandError("bcftools", timed_out=True), "bcftools_unavailable_or_timeout"),
        (ExternalCommandError("bcftools", 7), "bcftools_failed:7"),
    ],
)
def test_bcftools_errors_keep_the_stage_contract(
    monkeypatch, technical_error, expected_message
):
    def fail(*_args, **_kwargs):
        raise technical_error

    monkeypatch.setattr(call_explicit_ibd, "run_checked", fail)

    with pytest.raises(ExplicitIbdExternalError, match=expected_message):
        call_explicit_ibd._run(["bcftools", "query"], 10)


def _raw_scenario(
    *,
    calibrated=False,
    specificity_acceptable=False,
    control_frequency=None,
):
    return {
        "scenario_id": "primary",
        "role": "PRIMARY",
        "minimum_cm": 2.0,
        "minimum_markers": 100,
        "calibrated": calibrated,
        "calibration_specificity_acceptable": specificity_acceptable,
        "control_frequency": control_frequency,
    }


def test_scenario_preserves_disabled_prespecified_calibration():
    scenario = _scenario(_raw_scenario())

    assert not scenario.calibrated
    assert _calibration_status(scenario, 0.05) == "NOT_EVALUABLE"


def test_calibration_and_empirical_background_are_classified_separately():
    scenario = _scenario(
        _raw_scenario(
            calibrated=True,
            specificity_acceptable=True,
            control_frequency=0.01,
        )
    )

    checks = _audit_checks(
        (scenario,),
        {"primary": (105, 21, 0.2)},
        minimum_control_pairs=100,
        maximum_control_frequency=0.05,
    )
    statuses = {item["check"]: item["status"] for item in checks}

    assert statuses["prespecified_calibration_primary"] == "PASS"
    assert statuses["empirical_control_frequency_primary"] == "FAIL"


def test_empirical_background_requires_the_prespecified_pair_count():
    assert _empirical_control_status(91, 0.0, 100, 0.05) == "NOT_EVALUABLE"
    assert _empirical_control_status(105, 0.0, 100, 0.05) == "ACCEPTABLE"
