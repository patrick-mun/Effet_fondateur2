from pathlib import Path

import numpy as np

from effet_fondateur.calibration.contracts import STAGE_IDS
from effet_fondateur.calibration.pilot_panel import (
    CommandMeasure,
    _configuration,
    _regional_positions,
    _run_measured,
    _stage_statuses,
    _write_plink_panel,
    _write_regional_panel,
)


REPOSITORY_ROOT = Path(__file__).resolve().parents[2]


def test_r11_configuration_keeps_prespecified_pilot_dimensions() -> None:
    configuration = _configuration(REPOSITORY_ROOT)

    assert configuration["partition"] == "CALIBRATION"
    assert configuration["counts"]["study_identities"] == 500
    assert configuration["counts"]["phasing_reference"] == 500
    assert configuration["counts"]["ancestry_reference"] == 300
    assert configuration["counts"]["null_reference"] == 500
    assert configuration["variant_design"]["observed_densities"] == {
        "D_MAX": 600_000,
        "D_MID": 300_000,
        "D_LOW": 100_000,
    }
    assert configuration["counts"]["external_tool_seeds"] == 10
    assert configuration["null_draws_16B"] == 999
    assert configuration["qualification_flags"] == {
        "calibrated": False,
        "founder_effect_proven": False,
        "production_promotion": False,
    }


def test_plink_writer_emits_variant_major_bed_and_nested_lists(tmp_path: Path) -> None:
    target = {
        "chromosome": 7,
        "position_bp": 1_000,
        "variant_id": "SYN-TARGET",
    }
    panel = _write_plink_panel(
        tmp_path / "panel",
        sample_count=8,
        variant_count=44,
        seed=17,
        target=target,
    )

    prefix = panel["prefix"]
    assert prefix.with_suffix(".bed").read_bytes()[:3] == bytes((0x6C, 0x1B, 0x01))
    assert prefix.with_suffix(".bed").stat().st_size == 3 + 44 * 2
    assert len(prefix.with_suffix(".bim").read_text(encoding="utf-8").splitlines()) == 44
    assert "SYN-TARGET" in prefix.with_suffix(".bim").read_text(encoding="utf-8")
    assert len((prefix.parent / "dmid.extract").read_text(encoding="utf-8").splitlines()) == 22
    assert len((prefix.parent / "dlow.extract").read_text(encoding="utf-8").splitlines()) == 8


def test_regional_panel_keeps_truth_outside_blinded_inputs(tmp_path: Path) -> None:
    configuration = _configuration(REPOSITORY_ROOT)
    configuration["counts"] = {**configuration["counts"], "study_identities": 8, "phasing_reference": 6}
    configuration["variant_design"] = {
        **configuration["variant_design"],
        "regional_variants": 24,
        "regional_common_variants": 20,
        "regional_rare_variants": 4,
    }
    root = tmp_path / "release/blinded_inputs/CALIBRATION/regional"
    result = _write_regional_panel(root, configuration, seed=19)

    assert result["study_vcf"].is_file()
    assert result["reference_vcf"].is_file()
    assert result["truth_path"] == tmp_path / "release/private_truth/regional_haplotypes.npz"
    assert result["truth_path"].is_file()
    assert not list((tmp_path / "release/blinded_inputs").rglob("*.npz"))
    positions = _regional_positions(configuration)
    assert np.all(np.diff(positions) > 0)


def test_stage_statuses_fail_external_errors_after_readiness() -> None:
    measure = CommandMeasure(
        command_id="shapeit_common_seed_01",
        tool="shapeit5_phase_common",
        returncode=2,
        wall_seconds=0.1,
        maximum_rss_bytes=1,
        stdout_path="stdout.log",
        stderr_path="stderr.log",
        resource_path="resource.log",
        command=("shapeit",),
    )
    statuses = _stage_statuses(
        [measure],
        {"successful_shapeit_seed_count": 0, "ibd_outputs": {}},
    )

    assert tuple(statuses) == STAGE_IDS
    assert statuses["12"]["status"] == "FAIL"
    assert statuses["12"]["reason"].startswith("post_readiness_external_command_failed")
    assert statuses["12A"]["status"] == "NOT_EVALUABLE"


def test_measured_command_uses_macos_time_contract(tmp_path: Path) -> None:
    measure = _run_measured(
        command_id="echo",
        tool="shell_probe",
        command=["/bin/echo", "synthetic"],
        run_root=tmp_path,
        timeout_seconds=5,
    )

    assert measure.returncode == 0
    assert measure.maximum_rss_bytes is not None
    assert (tmp_path / measure.resource_path).is_file()
