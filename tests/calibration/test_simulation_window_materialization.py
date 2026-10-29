"""Test d'intégration locale du matérialisateur D14b."""

from __future__ import annotations

import json
from pathlib import Path

import jsonschema
import pytest

from effet_fondateur.calibration.materialize_simulation_windows import (
    materialize_simulation_windows,
)
from effet_fondateur.calibration.simulation.windows import SimulationWindowError


ROOT = Path(__file__).resolve().parents[2]


def _fixture(tmp_path: Path, *, minimum_markers: int = 5) -> dict[str, object]:
    source_vcf = (
        tmp_path
        / "1kGP_high_coverage_Illumina.chr5.filtered.SNV_INDEL_SV_phased_panel.vcf.gz"
    )
    source_vcf.write_bytes(b"public-reference-fixture")
    source_index = Path(f"{source_vcf}.tbi")
    source_index.write_bytes(b"public-index-fixture")
    map_path = tmp_path / "chr5.grch38.tsv"
    map_path.write_text(
        "MAP_ID\tASSEMBLY\tCHROMOSOME\tPOSITION_BP\tPOSITION_CM\n"
        "map\tGRCh38\t5\t100\t0\n"
        "map\tGRCh38\t5\t600\t5\n"
        "map\tGRCh38\t5\t1100\t10\n",
        encoding="utf-8",
    )
    positions = tuple(range(300, 901, 20))
    chip_path = tmp_path / "chip.tsv"
    chip_path.write_text(
        "MARKER_ID\tARRAY_DESIGN\tDESIGN_VERSION\tASSEMBLY\tCHROMOSOME\t"
        "POSITION_BP\tALLELE_A\tALLELE_B\tSTRAND\n"
        + "".join(
            f"M{position}\tCytoScan750K_Accel\tNA36.r1\tGRCh38\t5\t"
            f"{position}\tC\tT\t+\n"
            for position in positions
        ),
        encoding="utf-8",
    )
    vcf_payload = (
        "##fileformat=VCFv4.2\n"
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
        + "".join(
            f"chr5\t{position}\tv{position}\tC\tT\t.\tPASS\t.\n"
            for position in positions
        )
    )
    config = tmp_path / "d14b.json"
    config.write_text(
        json.dumps(
            {
                "schema_version": "1.0.0",
                "simulation_window_config_id": "d14b_test",
                "config_state": "FINALIZED",
                "panel_role": "CALIBRATION",
                "method": {
                    "method_id": "map_quantile_window_design_v1",
                    "method_version": "1.0.0",
                    "maximum_arm_censor_probability": 0.05,
                    "genetic_map_interpolation": (
                        "PIECEWISE_LINEAR_NO_EXTRAPOLATION"
                    ),
                    "candidate_selection": "PRESELECTED_PUBLIC_REFERENCE_TARGET",
                    "outcome_filtering": "FORBIDDEN",
                },
                "marker_requirements": {
                    "minimum_flank_markers_each_side": 2,
                    "ibd_scenarios": [
                        {
                            "scenario_id": "primary_2cm",
                            "segment_length_cm": 2,
                            "minimum_observed_markers": minimum_markers,
                        }
                    ],
                },
                "window_specs": [
                    {
                        "window_id": "primary",
                        "purpose": "PRIMARY",
                        "age_generations": 100,
                        "chromosome": 5,
                        "target": {
                            "variant_id": "SYNTHETIC_TARGET",
                            "position_bp": 600,
                            "reference_allele": "C",
                            "alternate_allele": "T",
                        },
                    },
                    {
                        "window_id": "challenge",
                        "purpose": "CENSORING_CHALLENGE",
                        "age_generations": 100,
                        "chromosome": 5,
                        "target": {
                            "variant_id": "SYNTHETIC_TARGET",
                            "position_bp": 600,
                            "reference_allele": "C",
                            "alternate_allele": "T",
                        },
                        "challenge_half_width_cm": 0.1,
                    },
                ],
            }
        ),
        encoding="utf-8",
    )

    def fake_runner(command: list[str]) -> None:
        output_path = Path(command[command.index("--output-file") + 1])
        output_path.write_text(vcf_payload, encoding="utf-8")

    return {
        "window_config_path": config,
        "source_vcf_path": source_vcf,
        "source_index_path": source_index,
        "genetic_map_path": map_path,
        "normalized_chip_markers_path": chip_path,
        "workspace_root": tmp_path / "calibration_workspace",
        "run_id": "d14b_test",
        "created_at": "2026-10-29T10:00:00Z",
        "command_runner": fake_runner,
    }


def test_materializer_is_deterministic_and_manifest_valid(tmp_path: Path) -> None:
    arguments = _fixture(tmp_path)
    manifest = materialize_simulation_windows(**arguments)
    schema = json.loads(
        (ROOT / "calibration/schemas/simulation_window_manifest.schema.json").read_text(
            encoding="utf-8"
        )
    )
    jsonschema.Draft202012Validator(
        schema, format_checker=jsonschema.FormatChecker()
    ).validate(manifest)
    assert manifest["integrity"]["ibd_not_evaluable_count"] == 0
    assert manifest["windows"][0]["ibd_scenarios"][0]["status"] == "EVALUABLE"
    manifest_text = (
        arguments["workspace_root"]
        / "d14b_test/private_truth/d14b/simulation_window_manifest.json"
    ).read_text(encoding="utf-8")

    second_root = tmp_path / "second/calibration_workspace"
    second_arguments = dict(arguments)
    second_arguments["workspace_root"] = second_root
    second_manifest = materialize_simulation_windows(**second_arguments)
    second_text = (
        second_root
        / "d14b_test/private_truth/d14b/simulation_window_manifest.json"
    ).read_text(encoding="utf-8")
    assert second_manifest == manifest
    assert second_text == manifest_text


def test_materializer_refuses_insufficient_marker_density(tmp_path: Path) -> None:
    arguments = _fixture(tmp_path, minimum_markers=100)
    with pytest.raises(
        SimulationWindowError,
        match="simulation_window_ibd_marker_threshold_not_met",
    ):
        materialize_simulation_windows(**arguments)
    assert not (arguments["workspace_root"] / "d14b_test").exists()
