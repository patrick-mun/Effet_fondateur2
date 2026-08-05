import csv
import json
from pathlib import Path

import pytest
import yaml

from effet_fondateur.orchestrator import IntegrityError, StageExecutionError
from effet_fondateur.orchestrator.pipeline import resume_pipeline, run_pipeline


REPOSITORY_ROOT = Path(__file__).resolve().parents[1]
ACPA_HEADER = (
    "# Array Type Name: CytoScan 750K Accel Array\n"
    "# UCSC Genomic Version: hg38\n"
    "Probe Set ID\tForward Strand Base Calls\tdbSNP RS ID\tChromosome\t"
    "Chromosomal Position\n"
)


def write_acpa_source(path: Path, *, missing_column: bool = False) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    header = ACPA_HEADER
    if missing_column:
        header = header.replace("\tdbSNP RS ID", "")
    rows = [
        f"probe_{chromosome}\tAA\trs{chromosome}\t{chromosome}\t{chromosome * 100}"
        for chromosome in range(1, 23)
    ]
    path.write_text(header + "\n".join(rows) + "\n", encoding="utf-8")


def write_config(path: Path, source_dir: Path, **parameters: object) -> None:
    config = yaml.safe_load(
        (REPOSITORY_ROOT / "config" / "pipeline.example.yaml").read_text(
            encoding="utf-8"
        )
    )
    config["inputs"]["acpa_samples_dir"] = str(source_dir)
    config["stages"]["validate_sources"] = {
        "enabled": True,
        "parameters": parameters,
    }
    path.write_text(
        yaml.safe_dump(config, allow_unicode=True, sort_keys=False),
        encoding="utf-8",
    )


def read_manifest(run_dir: Path) -> dict[str, object]:
    return json.loads((run_dir / "manifest.json").read_text(encoding="utf-8"))


def test_validate_sources_publishes_inventory_and_qc(tmp_path: Path) -> None:
    source_dir = tmp_path / "sources"
    write_acpa_source(source_dir / "sample_1.txt")
    write_acpa_source(source_dir / "sample_2.txt")
    config_path = tmp_path / "config.yaml"
    write_config(config_path, source_dir, expected_file_count=2)

    run_dir = run_pipeline(config_path, tmp_path / "runs")

    stage_dir = run_dir / "stages" / "01_validate_sources"
    manifest = read_manifest(run_dir)
    assert manifest["stages"][-1]["state"] == "SUCCEEDED"
    assert manifest["global_status"] == "TECHNICALLY_VALID"
    with (stage_dir / "source_inventory.tsv").open(encoding="utf-8") as input_file:
        inventory = list(csv.DictReader(input_file, delimiter="\t"))
    assert len(inventory) == 2
    assert {row["STATUS"] for row in inventory} == {"PASS"}
    assert {row["AUTOSOMES_PRESENT"] for row in inventory} == {
        ",".join(str(chromosome) for chromosome in range(1, 23))
    }
    audit = json.loads((stage_dir / "audit.json").read_text(encoding="utf-8"))
    assert audit["counts"] == {
        "source_files": 2,
        "passed_files": 2,
        "failed_checks": 0,
    }


def test_validate_sources_blocks_missing_required_column(tmp_path: Path) -> None:
    source_dir = tmp_path / "sources"
    write_acpa_source(source_dir / "invalid.txt", missing_column=True)
    config_path = tmp_path / "config.yaml"
    write_config(config_path, source_dir)
    runs_dir = tmp_path / "runs"

    with pytest.raises(StageExecutionError) as error:
        run_pipeline(config_path, runs_dir)

    assert error.value.return_code == 4
    run_dir = next(path for path in runs_dir.iterdir() if not path.name.startswith("."))
    failed_attempt = next((run_dir / "stages" / "attempts").iterdir())
    qc_text = (failed_attempt / "source_qc.tsv").read_text(encoding="utf-8")
    assert "acpa_structure\tFAIL\tmissing_required_columns" in qc_text
    assert read_manifest(run_dir)["global_status"] == "BLOCKED"


def test_validate_sources_rejects_duplicate_basenames(tmp_path: Path) -> None:
    source_dir = tmp_path / "sources"
    write_acpa_source(source_dir / "batch_1" / "sample.txt")
    write_acpa_source(source_dir / "batch_2" / "sample.txt")
    config_path = tmp_path / "config.yaml"
    write_config(config_path, source_dir)

    with pytest.raises(StageExecutionError) as error:
        run_pipeline(config_path, tmp_path / "runs")

    assert error.value.return_code == 4


def test_resume_blocks_source_modified_after_publication(tmp_path: Path) -> None:
    source_dir = tmp_path / "sources"
    source_path = source_dir / "sample.txt"
    write_acpa_source(source_path)
    config_path = tmp_path / "config.yaml"
    write_config(config_path, source_dir)
    run_dir = run_pipeline(config_path, tmp_path / "runs")
    source_path.write_text(
        source_path.read_text(encoding="utf-8") + "probe_extra\tAA\trs99\t19\t9999\n",
        encoding="utf-8",
    )

    with pytest.raises(IntegrityError, match="Entrées modifiées"):
        resume_pipeline(run_dir)

    assert read_manifest(run_dir)["global_status"] == "BLOCKED"


def test_validate_sources_rejects_symlinked_input(tmp_path: Path) -> None:
    source_dir = tmp_path / "sources"
    external_source = tmp_path / "external.txt"
    write_acpa_source(external_source)
    source_dir.mkdir()
    (source_dir / "linked.txt").symlink_to(external_source)
    config_path = tmp_path / "config.yaml"
    write_config(config_path, source_dir)
    runs_dir = tmp_path / "runs"

    with pytest.raises(StageExecutionError) as error:
        run_pipeline(config_path, runs_dir)

    assert error.value.return_code == 2
    run_dir = next(path for path in runs_dir.iterdir() if not path.name.startswith("."))
    manifest = read_manifest(run_dir)
    assert manifest["global_status"] == "BLOCKED"
    assert manifest["stages"][-1]["state"] == "FAILED"
