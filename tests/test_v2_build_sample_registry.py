import csv
import json
from pathlib import Path

import pytest
import yaml

from effet_fondateur.orchestrator import IntegrityError, StageExecutionError
from effet_fondateur.orchestrator.pipeline import resume_pipeline, run_pipeline


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


def write_acpa_source(path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    rows = [
        f"probe_{chromosome}\tAA\trs{chromosome}\t{chromosome}\t{chromosome * 100}"
        for chromosome in range(1, 23)
    ]
    path.write_text(ACPA_HEADER + "\n".join(rows) + "\n", encoding="utf-8")


def write_samples_metadata(
    path: Path,
    *,
    duplicate_source: bool = False,
    genotype_source: str = "laboratory_report",
) -> None:
    rows = (
        "sample_1\tsample_1.txt\tFAM1\tI1\t0\t0\tMALE\tAFFECTED\tFAMILY\t"
        f"A/G\t{genotype_source}\tbatch_1\ttrue\ttrue\t\n"
    )
    if duplicate_source:
        rows += (
            "sample_2\tsample_1.txt\tFAM2\tI2\t0\t0\tFEMALE\tUNAFFECTED\t"
            "FAMILY\t\t\tbatch_1\ttrue\ttrue\t\n"
        )
    path.write_text(SAMPLES_HEADER + rows, encoding="utf-8")


def write_config(
    path: Path,
    source_dir: Path,
    metadata_path: Path,
    approval: dict[str, object] | None = None,
) -> None:
    config = yaml.safe_load(
        (REPOSITORY_ROOT / "config" / "pipeline.example.yaml").read_text(
            encoding="utf-8"
        )
    )
    config["inputs"]["acpa_samples_dir"] = str(source_dir)
    config["inputs"]["sample_metadata"] = str(metadata_path)
    config["stages"]["validate_sources"]["enabled"] = True
    registry_parameters = {} if approval is None else {"manual_approval": approval}
    config["stages"]["build_sample_registry"] = {
        "enabled": True,
        "parameters": registry_parameters,
    }
    path.write_text(
        yaml.safe_dump(config, allow_unicode=True, sort_keys=False),
        encoding="utf-8",
    )


def read_manifest(run_dir: Path) -> dict[str, object]:
    return json.loads((run_dir / "manifest.json").read_text(encoding="utf-8"))


def test_registry_is_published_with_manual_review_pending(tmp_path: Path) -> None:
    source_dir = tmp_path / "sources"
    write_acpa_source(source_dir / "sample_1.txt")
    metadata_path = tmp_path / "samples.master.input.tsv"
    write_samples_metadata(metadata_path)
    config_path = tmp_path / "config.yaml"
    write_config(config_path, source_dir, metadata_path)

    run_dir = run_pipeline(config_path, tmp_path / "runs")

    stage_dir = run_dir / "stages" / "02_build_sample_registry"
    manifest = read_manifest(run_dir)
    assert manifest["stages"][-1]["state"] == "SUCCEEDED"
    assert manifest["global_status"] == "INCOMPLETE"
    assert "sample_registry_approval" in manifest["manual_decisions_required"]
    with (stage_dir / "sample_registry_review.tsv").open(
        encoding="utf-8"
    ) as input_file:
        review_rows = list(csv.DictReader(input_file, delimiter="\t"))
    assert review_rows[0]["REVIEW_STATUS"] == "PENDING"
    assert review_rows[0]["PSEUDONYM"].startswith("P_")


def test_explicit_approval_makes_registry_technically_valid(tmp_path: Path) -> None:
    source_dir = tmp_path / "sources"
    write_acpa_source(source_dir / "sample_1.txt")
    metadata_path = tmp_path / "samples.master.input.tsv"
    write_samples_metadata(metadata_path)
    config_path = tmp_path / "config.yaml"
    write_config(
        config_path,
        source_dir,
        metadata_path,
        approval={
            "approved": True,
            "decision_id": "synthetic_review_001",
            "reviewer_role": "data_steward",
            "approved_at": "2026-08-05T10:00:00Z",
        },
    )

    run_dir = run_pipeline(config_path, tmp_path / "runs")

    manifest = read_manifest(run_dir)
    assert manifest["global_status"] == "TECHNICALLY_VALID"
    assert "sample_registry_approval" not in manifest["manual_decisions_required"]
    audit = json.loads(
        (
            run_dir / "stages" / "02_build_sample_registry" / "audit.json"
        ).read_text(encoding="utf-8")
    )
    assert audit["manual_validation_required"] is False


def test_registry_rejects_source_assigned_twice(tmp_path: Path) -> None:
    source_dir = tmp_path / "sources"
    write_acpa_source(source_dir / "sample_1.txt")
    metadata_path = tmp_path / "samples.master.input.tsv"
    write_samples_metadata(metadata_path, duplicate_source=True)
    config_path = tmp_path / "config.yaml"
    write_config(config_path, source_dir, metadata_path)
    runs_dir = tmp_path / "runs"

    with pytest.raises(StageExecutionError) as error:
        run_pipeline(config_path, runs_dir)

    assert error.value.return_code == 2
    run_dir = next(path for path in runs_dir.iterdir() if not path.name.startswith("."))
    assert read_manifest(run_dir)["stages"][-1]["state"] == "FAILED"


def test_registry_rejects_clinical_genotype_source(tmp_path: Path) -> None:
    source_dir = tmp_path / "sources"
    write_acpa_source(source_dir / "sample_1.txt")
    metadata_path = tmp_path / "samples.master.input.tsv"
    write_samples_metadata(metadata_path, genotype_source="clinical_status")
    config_path = tmp_path / "config.yaml"
    write_config(config_path, source_dir, metadata_path)

    with pytest.raises(StageExecutionError) as error:
        run_pipeline(config_path, tmp_path / "runs")

    assert error.value.return_code == 2


def test_resume_blocks_modified_sample_metadata(tmp_path: Path) -> None:
    source_dir = tmp_path / "sources"
    write_acpa_source(source_dir / "sample_1.txt")
    metadata_path = tmp_path / "samples.master.input.tsv"
    write_samples_metadata(metadata_path)
    config_path = tmp_path / "config.yaml"
    write_config(config_path, source_dir, metadata_path)
    run_dir = run_pipeline(config_path, tmp_path / "runs")
    metadata_path.write_text(
        metadata_path.read_text(encoding="utf-8").replace("batch_1", "batch_2"),
        encoding="utf-8",
    )

    with pytest.raises(IntegrityError, match="Entrées modifiées"):
        resume_pipeline(run_dir)

    assert read_manifest(run_dir)["global_status"] == "BLOCKED"


def test_registry_is_audited_when_source_validation_is_disabled(
    tmp_path: Path,
) -> None:
    source_dir = tmp_path / "sources"
    write_acpa_source(source_dir / "sample_1.txt")
    metadata_path = tmp_path / "samples.master.input.tsv"
    write_samples_metadata(metadata_path)
    config_path = tmp_path / "config.yaml"
    write_config(config_path, source_dir, metadata_path)
    config = yaml.safe_load(config_path.read_text(encoding="utf-8"))
    config["stages"]["validate_sources"]["enabled"] = False
    config_path.write_text(
        yaml.safe_dump(config, allow_unicode=True, sort_keys=False),
        encoding="utf-8",
    )
    runs_dir = tmp_path / "runs"

    with pytest.raises(StageExecutionError) as error:
        run_pipeline(config_path, runs_dir)

    assert error.value.return_code == 4
    run_dir = next(path for path in runs_dir.iterdir() if not path.name.startswith("."))
    manifest = read_manifest(run_dir)
    assert manifest["global_status"] == "BLOCKED"
    assert manifest["stages"][-1]["state"] == "FAILED"
