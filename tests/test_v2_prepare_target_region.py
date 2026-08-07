import csv
import json
import sys
from pathlib import Path

import pytest
import yaml

from effet_fondateur.orchestrator import StageExecutionError
from effet_fondateur.orchestrator.pipeline import run_pipeline
from test_v2_qc_final import prepare_final_qc_inputs


GENETIC_MAP_HEADER = "MAP_ID\tASSEMBLY\tCHROMOSOME\tPOSITION_BP\tPOSITION_CM\n"


def write_target_region_plink(
    path: Path,
    delegated_plink: Path,
    *,
    fail_region: bool = False,
) -> None:
    script = f'''#!{sys.executable}
import os
import pathlib
import sys

if sys.argv[1:] == ["--version"]:
    print("PLINK v1.90 synthetic")
    raise SystemExit(0)
output_prefix = pathlib.Path(sys.argv[sys.argv.index("--out") + 1]) if "--out" in sys.argv else None
if output_prefix is None or output_prefix.name != "target_region":
    os.execv({str(delegated_plink)!r}, [{str(delegated_plink)!r}, *sys.argv[1:]])
if {fail_region!r}:
    raise SystemExit(8)
base_prefix = pathlib.Path(sys.argv[sys.argv.index("--bfile") + 1])
chromosome = sys.argv[sys.argv.index("--chr") + 1]
start_bp = int(sys.argv[sys.argv.index("--from-bp") + 1])
end_bp = int(sys.argv[sys.argv.index("--to-bp") + 1])
fam_rows = [line for line in base_prefix.with_suffix(".fam").read_text().splitlines() if line]
bim_rows = [line.split() for line in base_prefix.with_suffix(".bim").read_text().splitlines() if line]
bim_rows = [
    row for row in bim_rows
    if row[0] == chromosome and start_bp <= int(row[3]) <= end_bp
]
output_prefix.with_suffix(".fam").write_text("".join(line + "\\n" for line in fam_rows))
output_prefix.with_suffix(".bim").write_text("".join(" ".join(row) + "\\n" for row in bim_rows))
output_prefix.with_suffix(".bed").write_bytes(b"\\x6c\\x1b\\x01")
'''
    path.write_text(script, encoding="utf-8")
    path.chmod(0o755)


def write_genetic_map(
    path: Path,
    *,
    assembly: str = "GRCh38",
    points: list[tuple[int, str]] | None = None,
) -> None:
    map_points = points or [(1, "0"), (100_000, "1"), (200_000, "2")]
    path.write_text(
        GENETIC_MAP_HEADER
        + "".join(
            f"synthetic_grch38_map\t{assembly}\t19\t{position_bp}\t{position_cm}\n"
            for position_bp, position_cm in map_points
        ),
        encoding="utf-8",
    )


def prepare_target_region_inputs(
    tmp_path: Path,
    *,
    map_assembly: str = "GRCh38",
    map_points: list[tuple[int, str]] | None = None,
    fail_region: bool = False,
    parameter_overrides: dict[str, object] | None = None,
) -> tuple[Path, Path]:
    config_path, runs_dir = prepare_final_qc_inputs(tmp_path)
    config = yaml.safe_load(config_path.read_text(encoding="utf-8"))
    config["stages"]["qc_final"]["parameters"][
        "target_region_variant_missingness_max"
    ] = 0.6
    config["stages"]["qc_final"]["parameters"]["target_region_maf_min"] = 0
    delegated_plink = Path(config["tools"]["plink"])
    region_plink = tmp_path / "target_region_plink"
    write_target_region_plink(region_plink, delegated_plink, fail_region=fail_region)
    genetic_map_path = tmp_path / "genetic_map.tsv"
    write_genetic_map(
        genetic_map_path,
        assembly=map_assembly,
        points=map_points,
    )
    config["tools"]["plink"] = str(region_plink)
    config["inputs"]["genetic_map"] = str(genetic_map_path)
    parameters: dict[str, object] = {
        "window_left_bp": 100_000,
        "window_right_bp": 100_000,
        "max_interpolation_gap_bp": 150_000,
        "min_region_variants": 2,
        "plink_timeout_seconds": 30,
    }
    parameters.update(parameter_overrides or {})
    config["stages"]["prepare_target_region"] = {
        "enabled": True,
        "parameters": parameters,
    }
    config_path.write_text(
        yaml.safe_dump(config, allow_unicode=True, sort_keys=False),
        encoding="utf-8",
    )
    return config_path, runs_dir


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8") as input_file:
        return list(csv.DictReader(input_file, delimiter="\t"))


def test_target_region_interpolates_map_and_preserves_variant_order(
    tmp_path: Path,
) -> None:
    config_path, runs_dir = prepare_target_region_inputs(tmp_path)

    run_dir = run_pipeline(config_path, runs_dir)

    stage_dir = run_dir / "stages" / "11_prepare_target_region"
    map_rows = read_tsv(stage_dir / "target_genetic_map.tsv")
    assert [row["VARIANT_ID"] for row in map_rows] == ["probe_19", "target_GRCh38_1_100000_A_G"]
    assert [row["VARIANT_ORDER"] for row in map_rows] == ["1", "2"]
    assert map_rows[0]["MAP_STATUS"] == "INTERPOLATED"
    assert map_rows[0]["LEFT_MAP_BP"] == "1"
    assert map_rows[0]["RIGHT_MAP_BP"] == "100000"
    assert map_rows[1]["MAP_STATUS"] == "EXACT"
    assert map_rows[1]["POSITION_CM"] == "1"
    bim_rows = [
        line.split()
        for line in (stage_dir / "target_region.bim").read_text(encoding="utf-8").splitlines()
    ]
    assert [row[1] for row in bim_rows] == [row["VARIANT_ID"] for row in map_rows]
    assert [row[2] for row in bim_rows] == [row["POSITION_CM"] for row in map_rows]


def test_target_region_publishes_generic_phasing_manifest(tmp_path: Path) -> None:
    config_path, runs_dir = prepare_target_region_inputs(tmp_path)

    run_dir = run_pipeline(config_path, runs_dir)

    stage_dir = run_dir / "stages" / "11_prepare_target_region"
    manifest = json.loads(
        (stage_dir / "phasing_input_manifest.json").read_text(encoding="utf-8")
    )
    descriptor = json.loads(
        (stage_dir / "target_region.dataset.json").read_text(encoding="utf-8")
    )
    assert manifest["format"] == "PLINK_BED_BIM_FAM_WITH_CM"
    assert manifest["target_variant_order"] == 2
    assert manifest["variant_set_id"] == descriptor["variant_set_id"]
    assert descriptor["source_format"] == "PLINK_TARGET_REGION"
    assert descriptor["sample_count"] == 3
    assert descriptor["variant_count"] == 2


def test_stage_16_publishes_non_evaluation_on_sparse_synthetic_chain(
    tmp_path: Path,
) -> None:
    config_path, runs_dir = prepare_target_region_inputs(tmp_path)
    config = yaml.safe_load(config_path.read_text(encoding="utf-8"))
    config["stages"]["analyze_roh"] = {
        "enabled": True,
        "parameters": _parameters_for_sparse_roh_test(),
    }
    config_path.write_text(
        yaml.safe_dump(config, allow_unicode=True, sort_keys=False), encoding="utf-8"
    )

    run_dir = run_pipeline(config_path, runs_dir)

    stage_dir = run_dir / "stages" / "16_analyze_roh"
    summary = json.loads(
        (stage_dir / "roh" / "roh_analysis_summary.json").read_text(encoding="utf-8")
    )
    assert set(summary["scope_statuses"].values()) == {
        "NOT_EVALUATED_DENSITY_OR_COVERAGE"
    }
    assert summary["feeds_founder_haplotype"] is False
    assert summary["feeds_variant_age"] is False
    assert summary["consumes_local_ld"] is False
    target_rows = read_tsv(stage_dir / "roh" / "target_in_roh.tsv")
    assert len(target_rows) == 3
    assert {row["INTERPRETATION_STATUS"] for row in target_rows} == {
        "NOT_EVALUATED_SCOPE"
    }


def _parameters_for_sparse_roh_test() -> dict[str, object]:
    return {
        "method": "plink19_array_roh_secondary_v1",
        "minimum_roh_kb": 1500,
        "minimum_roh_snps": 50,
        "maximum_density_kb_per_snp": 50,
        "maximum_gap_kb": 1000,
        "window_snps": 50,
        "window_max_heterozygotes": 1,
        "window_max_missing": 5,
        "window_threshold": 0.05,
        "minimum_genomewide_variants": 10000,
        "minimum_genomewide_autosomes": 22,
        "minimum_target_chromosome_variants": 100,
        "minimum_group_samples": 5,
        "minimum_primary_samples": 20,
        "autosomal_denominator_kb": None,
        "autosomal_denominator_source": None,
        "plink_timeout_seconds": 30,
    }


def test_target_region_blocks_map_extrapolation(tmp_path: Path) -> None:
    config_path, runs_dir = prepare_target_region_inputs(
        tmp_path,
        map_points=[(2_000, "0.1"), (200_000, "2")],
    )

    with pytest.raises(StageExecutionError) as error:
        run_pipeline(config_path, runs_dir)

    assert error.value.return_code == 4


def test_target_region_blocks_nonmonotonic_genetic_map(tmp_path: Path) -> None:
    config_path, runs_dir = prepare_target_region_inputs(
        tmp_path,
        map_points=[(1, "0"), (100_000, "2"), (200_000, "1")],
    )

    with pytest.raises(StageExecutionError) as error:
        run_pipeline(config_path, runs_dir)

    assert error.value.return_code == 4


def test_target_region_blocks_large_interpolation_gap(tmp_path: Path) -> None:
    config_path, runs_dir = prepare_target_region_inputs(
        tmp_path,
        map_points=[(1, "0"), (100_000, "1"), (1_000_000, "3")],
        parameter_overrides={"max_interpolation_gap_bp": 50_000},
    )

    with pytest.raises(StageExecutionError) as error:
        run_pipeline(config_path, runs_dir)

    assert error.value.return_code == 4


def test_target_region_rejects_map_assembly_mismatch(tmp_path: Path) -> None:
    config_path, runs_dir = prepare_target_region_inputs(
        tmp_path,
        map_assembly="GRCh37",
    )

    with pytest.raises(StageExecutionError) as error:
        run_pipeline(config_path, runs_dir)

    assert error.value.return_code == 2


def test_target_region_plink_failure_uses_external_tool_code(tmp_path: Path) -> None:
    config_path, runs_dir = prepare_target_region_inputs(tmp_path, fail_region=True)

    with pytest.raises(StageExecutionError) as error:
        run_pipeline(config_path, runs_dir)

    assert error.value.return_code == 3


def test_target_region_invalid_parameter_uses_input_code(tmp_path: Path) -> None:
    config_path, runs_dir = prepare_target_region_inputs(
        tmp_path,
        parameter_overrides={"window_left_bp": 0},
    )

    with pytest.raises(StageExecutionError) as error:
        run_pipeline(config_path, runs_dir)

    assert error.value.return_code == 2
