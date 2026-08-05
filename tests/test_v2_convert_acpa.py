import csv
import json
import sys
from pathlib import Path

import pytest
import yaml

from effet_fondateur.orchestrator import StageExecutionError
from effet_fondateur.orchestrator.pipeline import run_pipeline


REPOSITORY_ROOT = Path(__file__).resolve().parents[1]
SAMPLES_HEADER = (
    "SAMPLE_ID\tSOURCE_FILE\tFID\tIID\tPID\tMID\tSEX\tCLINICAL_STATUS\t"
    "GROUP_LABEL\tTARGET_GENOTYPE\tTARGET_GENOTYPE_SOURCE\tARRAY_BATCH\t"
    "INCLUDE_GENOMEWIDE\tINCLUDE_TARGET_CHROMOSOME\tNOTES_CODE\n"
)
ACPA_HEADER = (
    "# Array Type Name: CytoScan 750K Accel Array\n"
    "# UCSC Genomic Version: hg38\n"
    "Probe Set ID\tForward Strand Base Calls\tdbSNP RS ID\tChromosome\t"
    "Chromosomal Position\n"
)


def write_fake_plink(path: Path, *, fail: bool = False) -> None:
    script = f"""#!{sys.executable}
import pathlib
import sys

if '--version' in sys.argv:
    print('PLINK fake 1.0')
    raise SystemExit(0)
if {fail!r}:
    raise SystemExit(8)
input_prefix = pathlib.Path(sys.argv[sys.argv.index('--file') + 1])
output_prefix = pathlib.Path(sys.argv[sys.argv.index('--out') + 1])
ped_rows = [
    line.split()
    for line in input_prefix.with_suffix('.ped').read_text().splitlines()
    if line
]
map_rows = [
    line.split()
    for line in input_prefix.with_suffix('.map').read_text().splitlines()
    if line
]
output_prefix.with_suffix('.fam').write_text(
    ''.join(' '.join(row[:6]) + '\\n' for row in ped_rows)
)
output_prefix.with_suffix('.bim').write_text(
    ''.join(f'{{row[0]}} {{row[1]}} {{row[2]}} {{row[3]}} A G\\n' for row in map_rows)
)
output_prefix.with_suffix('.bed').write_bytes(b'\\x6c\\x1b\\x01')
"""
    path.write_text(script, encoding="utf-8")
    path.chmod(0o755)


def acpa_rows(
    *,
    extra_target: bool = True,
    extra_position: int = 1950,
    invalid_chromosome_1_genotype: bool = False,
) -> list[str]:
    rows = []
    for chromosome in range(1, 23):
        genotype = "NN" if chromosome == 1 and invalid_chromosome_1_genotype else "AA"
        rows.append(
            f"probe_{chromosome}\t{genotype}\trs{chromosome}\t{chromosome}\t"
            f"{chromosome * 100}"
        )
    if extra_target:
        rows.append(f"probe_target_extra\tAG\trs1950\t19\t{extra_position}")
    return rows


def write_acpa_source(path: Path, rows: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(ACPA_HEADER + "\n".join(rows) + "\n", encoding="utf-8")


def write_samples_metadata(path: Path) -> None:
    path.write_text(
        SAMPLES_HEADER
        + "sample_1\tsample_1.txt\tFAM1\tI1\t0\t0\tMALE\tAFFECTED\tFAMILY\t"
        "A/G\tlaboratory_report\tbatch_1\ttrue\ttrue\t\n"
        + "sample_2\tsample_2.txt\tFAM2\tI2\t0\t0\tFEMALE\tUNAFFECTED\tFAMILY\t"
        "\t\tbatch_1\ttrue\ttrue\t\n",
        encoding="utf-8",
    )


def write_config(
    path: Path,
    source_dir: Path,
    metadata_path: Path,
    plink_path: Path,
    *,
    marker_mode: str = "intersection",
    approved: bool = True,
) -> None:
    config = yaml.safe_load(
        (REPOSITORY_ROOT / "config" / "pipeline.example.yaml").read_text(
            encoding="utf-8"
        )
    )
    config["inputs"]["acpa_samples_dir"] = str(source_dir)
    config["inputs"]["sample_metadata"] = str(metadata_path)
    config["tools"]["plink"] = str(plink_path)
    config["target"]["chromosome"] = 19
    config["stages"]["validate_sources"]["enabled"] = True
    registry_parameters = {}
    if approved:
        registry_parameters["manual_approval"] = {
            "approved": True,
            "decision_id": "synthetic_registry_review",
            "reviewer_role": "data_steward",
            "approved_at": "2026-08-05T10:00:00Z",
        }
    config["stages"]["build_sample_registry"] = {
        "enabled": True,
        "parameters": registry_parameters,
    }
    config["stages"]["convert_acpa"] = {
        "enabled": True,
        "parameters": {"marker_mode": marker_mode},
    }
    path.write_text(
        yaml.safe_dump(config, allow_unicode=True, sort_keys=False),
        encoding="utf-8",
    )


def prepare_inputs(
    tmp_path: Path,
    *,
    sample_1_rows: list[str] | None = None,
    sample_2_rows: list[str] | None = None,
    marker_mode: str = "intersection",
    approved: bool = True,
    plink_fails: bool = False,
) -> tuple[Path, Path]:
    source_dir = tmp_path / "sources"
    write_acpa_source(source_dir / "sample_1.txt", sample_1_rows or acpa_rows())
    write_acpa_source(source_dir / "sample_2.txt", sample_2_rows or acpa_rows())
    metadata_path = tmp_path / "samples.master.input.tsv"
    write_samples_metadata(metadata_path)
    plink_path = tmp_path / "fake_plink"
    write_fake_plink(plink_path, fail=plink_fails)
    config_path = tmp_path / "config.yaml"
    write_config(
        config_path,
        source_dir,
        metadata_path,
        plink_path,
        marker_mode=marker_mode,
        approved=approved,
    )
    return config_path, tmp_path / "runs"


def read_manifest(run_dir: Path) -> dict[str, object]:
    return json.loads((run_dir / "manifest.json").read_text(encoding="utf-8"))


def test_convert_acpa_publishes_two_aligned_plink_datasets(tmp_path: Path) -> None:
    config_path, runs_dir = prepare_inputs(tmp_path)

    run_dir = run_pipeline(config_path, runs_dir)

    stage_dir = run_dir / "stages" / "03_convert_acpa"
    assert (stage_dir / "genomewide_base.bed").read_bytes()[:3] == b"\x6c\x1b\x01"
    assert (stage_dir / "target_chromosome_base.bed").is_file()
    assert len((stage_dir / "genomewide_base.bim").read_text().splitlines()) == 23
    assert len((stage_dir / "target_chromosome_base.bim").read_text().splitlines()) == 2
    assert len((stage_dir / "genomewide_base.fam").read_text().splitlines()) == 2
    assert (stage_dir / "genomewide_base.fam").read_text() == (
        stage_dir / "target_chromosome_base.fam"
    ).read_text()
    report = json.loads(
        (stage_dir / "acpa_conversion_report.json").read_text(encoding="utf-8")
    )
    assert report["source_read_counts"] == {"sample_1.txt": 1, "sample_2.txt": 1}
    assert read_manifest(run_dir)["global_status"] == "TECHNICALLY_VALID"


def test_union_keeps_marker_missing_from_one_sample(tmp_path: Path) -> None:
    config_path, runs_dir = prepare_inputs(
        tmp_path,
        sample_2_rows=acpa_rows(extra_target=False),
        marker_mode="union",
    )

    run_dir = run_pipeline(config_path, runs_dir)

    audit_path = run_dir / "stages" / "03_convert_acpa" / "acpa_variant_audit.tsv"
    with audit_path.open(encoding="utf-8") as input_file:
        rows = list(csv.DictReader(input_file, delimiter="\t"))
    extra_marker = next(row for row in rows if row["PROBE_ID"] == "probe_target_extra")
    assert extra_marker["STATUS"] == "ACCEPTED"
    assert extra_marker["PRESENCE_COUNT"] == "1"


def test_intersection_excludes_marker_missing_from_one_sample(tmp_path: Path) -> None:
    config_path, runs_dir = prepare_inputs(
        tmp_path,
        sample_2_rows=acpa_rows(extra_target=False),
        marker_mode="intersection",
    )

    run_dir = run_pipeline(config_path, runs_dir)

    stage_dir = run_dir / "stages" / "03_convert_acpa"
    assert len((stage_dir / "target_chromosome_base.bim").read_text().splitlines()) == 1
    with (stage_dir / "acpa_variant_audit.tsv").open(encoding="utf-8") as input_file:
        rows = list(csv.DictReader(input_file, delimiter="\t"))
    extra_marker = next(row for row in rows if row["PROBE_ID"] == "probe_target_extra")
    assert extra_marker["STATUS"] == "EXCLUDED"
    assert extra_marker["EXCLUSION_CODE"] == "missing_in_samples"


def test_coordinate_conflict_is_excluded_from_both_datasets(tmp_path: Path) -> None:
    config_path, runs_dir = prepare_inputs(
        tmp_path,
        sample_2_rows=acpa_rows(extra_position=1960),
    )

    run_dir = run_pipeline(config_path, runs_dir)

    stage_dir = run_dir / "stages" / "03_convert_acpa"
    assert len((stage_dir / "target_chromosome_base.bim").read_text().splitlines()) == 1
    with (stage_dir / "acpa_variant_audit.tsv").open(encoding="utf-8") as input_file:
        rows = list(csv.DictReader(input_file, delimiter="\t"))
    extra_marker = next(row for row in rows if row["PROBE_ID"] == "probe_target_extra")
    assert extra_marker["STATUS"] == "EXCLUDED"
    assert extra_marker["EXCLUSION_CODE"] == "coordinate_conflict"


def test_invalid_genotype_is_converted_to_missing_and_audited(tmp_path: Path) -> None:
    config_path, runs_dir = prepare_inputs(
        tmp_path,
        sample_1_rows=acpa_rows(invalid_chromosome_1_genotype=True),
    )

    run_dir = run_pipeline(config_path, runs_dir)

    report = json.loads(
        (
            run_dir
            / "stages"
            / "03_convert_acpa"
            / "acpa_conversion_report.json"
        ).read_text(encoding="utf-8")
    )
    assert report["invalid_genotype_count"] == 1
    qc_path = (
        run_dir / "stages" / "03_convert_acpa" / "acpa_sample_chromosome_qc.tsv"
    )
    with qc_path.open(encoding="utf-8") as input_file:
        qc_rows = list(csv.DictReader(input_file, delimiter="\t"))
    chromosome_1_qc = next(
        row
        for row in qc_rows
        if row["SAMPLE_ID"] == "sample_1" and row["CHROMOSOME"] == "1"
    )
    assert chromosome_1_qc["CALLED_GENOTYPE_COUNT"] == "0"
    assert chromosome_1_qc["GENOTYPING_RATE"] == "0.000000"


def test_conversion_is_blocked_without_registry_approval(tmp_path: Path) -> None:
    config_path, runs_dir = prepare_inputs(tmp_path, approved=False)

    with pytest.raises(StageExecutionError) as error:
        run_pipeline(config_path, runs_dir)

    assert error.value.return_code == 4
    run_dir = next(path for path in runs_dir.iterdir() if not path.name.startswith("."))
    manifest = read_manifest(run_dir)
    assert "sample_registry_approval" in manifest["manual_decisions_required"]
    assert manifest["stages"][-1]["state"] == "FAILED"


def test_plink_failure_uses_external_tool_return_code(tmp_path: Path) -> None:
    config_path, runs_dir = prepare_inputs(tmp_path, plink_fails=True)

    with pytest.raises(StageExecutionError) as error:
        run_pipeline(config_path, runs_dir)

    assert error.value.return_code == 3
