import json
from pathlib import Path

import pytest

from simulation_genotype_famille.simulate_unrelated_controls import (
    ControlSimulationError,
    simulate_unrelated_controls,
)


def write_source_dataset(directory: Path) -> tuple[Path, Path]:
    source_prefix = directory / "source" / "genotype_data"
    source_prefix.parent.mkdir()
    source_prefix.with_suffix(".map").write_text(
        "19\trs1\t0\t100\n"
        "19\trs2\t0\t200\n",
        encoding="utf-8",
    )
    source_prefix.with_suffix(".ped").write_text(
        "F1\tCASE\t0\t0\t1\t2\tA\tA\t0\t0\n"
        "F2\tSOURCE_CONTROL\t0\t0\t2\t1\tA\tG\tC\tC\n",
        encoding="utf-8",
    )
    groups_path = source_prefix.parent / "groupes.txt"
    groups_path.write_text(
        "CASE\tATTEINT\nSOURCE_CONTROL\tSAINS\n",
        encoding="utf-8",
    )
    return source_prefix, groups_path


def test_simulation_creates_combined_reproducible_dataset(tmp_path):
    source_prefix, groups_path = write_source_dataset(tmp_path)
    original_ped = source_prefix.with_suffix(".ped").read_text()
    output_prefix = tmp_path / "output" / "genotype_data"

    report = simulate_unrelated_controls(
        source_prefix=source_prefix,
        output_prefix=output_prefix,
        control_count=3,
        random_seed=1234,
        source_groups_path=groups_path,
        simulate_missingness=False,
    )

    ped_rows = [line.split() for line in output_prefix.with_suffix(".ped").read_text().splitlines()]
    simulated_rows = ped_rows[2:]
    assert len(ped_rows) == 5
    assert all(len(row) == 10 for row in ped_rows)
    assert [row[0] for row in simulated_rows] == [
        "SIMCTRL_0001",
        "SIMCTRL_0002",
        "SIMCTRL_0003",
    ]
    assert all(row[0] == row[1] for row in simulated_rows)
    assert all(row[2:4] == ["0", "0"] for row in simulated_rows)
    assert all(row[5] == "1" for row in simulated_rows)
    assert output_prefix.with_suffix(".map").read_text() == (
        source_prefix.with_suffix(".map").read_text()
    )
    assert source_prefix.with_suffix(".ped").read_text() == original_ped
    assert (output_prefix.parent / "cas.txt").read_text().strip() == "F1\tCASE"
    assert len((output_prefix.parent / "temoins.txt").read_text().splitlines()) == 3
    assert (output_prefix.parent / "groupes.txt").read_text().splitlines() == [
        "CASE\tATTEINT",
        "SOURCE_CONTROL\tSAINS",
        "SIMCTRL_0001\tTEMOIN_SIMULE",
        "SIMCTRL_0002\tTEMOIN_SIMULE",
        "SIMCTRL_0003\tTEMOIN_SIMULE",
    ]
    saved_report = json.loads(
        (output_prefix.parent / "control_simulation_report.json").read_text()
    )
    assert saved_report["random_seed"] == 1234
    assert saved_report["simulated_control_count"] == 3
    assert saved_report["output_individual_count"] == 5
    assert report == saved_report


def test_same_seed_produces_same_controls(tmp_path):
    source_prefix, _groups_path = write_source_dataset(tmp_path)
    first_prefix = tmp_path / "first" / "genotype_data"
    second_prefix = tmp_path / "second" / "genotype_data"

    for output_prefix in (first_prefix, second_prefix):
        simulate_unrelated_controls(
            source_prefix=source_prefix,
            output_prefix=output_prefix,
            control_count=4,
            random_seed=9876,
            include_source=False,
        )

    assert first_prefix.with_suffix(".ped").read_text() == (
        second_prefix.with_suffix(".ped").read_text()
    )


def test_simulation_refuses_overwrite_and_invalid_count(tmp_path):
    source_prefix, _groups_path = write_source_dataset(tmp_path)
    output_prefix = tmp_path / "output" / "genotype_data"
    simulate_unrelated_controls(
        source_prefix=source_prefix,
        output_prefix=output_prefix,
        control_count=1,
    )

    with pytest.raises(ControlSimulationError, match="Sorties déjà présentes"):
        simulate_unrelated_controls(
            source_prefix=source_prefix,
            output_prefix=output_prefix,
            control_count=1,
        )
    with pytest.raises(ControlSimulationError, match="strictement positif"):
        simulate_unrelated_controls(
            source_prefix=source_prefix,
            output_prefix=tmp_path / "invalid" / "genotype_data",
            control_count=0,
        )
    with pytest.raises(ControlSimulationError, match="distinct du dossier source"):
        simulate_unrelated_controls(
            source_prefix=source_prefix,
            output_prefix=source_prefix.parent / "simulated",
            control_count=1,
            force=True,
        )
