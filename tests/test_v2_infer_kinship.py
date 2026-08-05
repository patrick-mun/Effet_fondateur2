import csv
import json
import sys
from pathlib import Path

import pytest
import yaml

from effet_fondateur.orchestrator import StageExecutionError
from effet_fondateur.orchestrator.pipeline import run_pipeline
from effet_fondateur.stages.infer_kinship import _classify_relationship
from test_v2_qc_preliminary import prepare_inputs


def write_fake_king(
    path: Path,
    *,
    between_rows: list[str] | None = None,
    within_rows: list[str] | None = None,
    fail: bool = False,
) -> None:
    script = f'''#!{sys.executable}
import pathlib
import sys

if "--version" in sys.argv:
    print("KING fake 2.3.2")
    raise SystemExit(0)
if {fail!r}:
    raise SystemExit(8)
prefix = pathlib.Path(sys.argv[sys.argv.index("--prefix") + 1])
between_rows = {between_rows or []!r}
within_rows = {within_rows or []!r}
prefix.with_suffix(".kin0").write_text(
    "FID1 ID1 FID2 ID2 N_SNP HetHet IBS0 Kinship\\n"
    + "".join(row + "\\n" for row in between_rows)
)
if within_rows:
    prefix.with_suffix(".kin").write_text(
        "FID ID1 ID2 N_SNP Z0 Phi HetHet IBS0 Kinship Error\\n"
        + "".join(row + "\\n" for row in within_rows)
    )
'''
    path.write_text(script, encoding="utf-8")
    path.chmod(0o755)


def update_metadata_for_parent(
    config_path: Path, *, external_parent: bool = False
) -> None:
    config = yaml.safe_load(config_path.read_text(encoding="utf-8"))
    metadata_path = Path(config["inputs"]["sample_metadata"])
    with metadata_path.open(encoding="utf-8") as input_file:
        rows = list(csv.DictReader(input_file, delimiter="\t"))
        fieldnames = list(rows[0])
    if external_parent:
        rows[1]["PID"] = "EXTERNAL_PARENT"
    else:
        rows[1]["FID"] = rows[0]["FID"]
        rows[1]["PID"] = rows[0]["IID"]
    with metadata_path.open("w", encoding="utf-8", newline="") as output_file:
        writer = csv.DictWriter(
            output_file,
            fieldnames=fieldnames,
            delimiter="\t",
            lineterminator="\n",
        )
        writer.writeheader()
        writer.writerows(rows)


def prepare_kinship_inputs(
    tmp_path: Path,
    *,
    between_rows: list[str] | None = None,
    within_rows: list[str] | None = None,
    fail_king: bool = False,
    declared_parent: bool = False,
    external_parent: bool = False,
) -> tuple[Path, Path]:
    config_path, runs_dir, _ = prepare_inputs(
        tmp_path,
        enable_kinship_panel=True,
    )
    if declared_parent or external_parent:
        update_metadata_for_parent(config_path, external_parent=external_parent)
    king_path = tmp_path / "fake_king"
    write_fake_king(
        king_path,
        between_rows=between_rows,
        within_rows=within_rows,
        fail=fail_king,
    )
    config = yaml.safe_load(config_path.read_text(encoding="utf-8"))
    config["tools"]["king"] = str(king_path)
    config["stages"]["infer_kinship"] = {
        "enabled": True,
        "parameters": {
            "king_degree": 4,
            "duplicate_or_mz_min": 0.354,
            "first_degree_min": 0.177,
            "second_degree_min": 0.0884,
            "third_degree_min": 0.0442,
            "fourth_degree_min": 0.0221,
            "parent_child_ibs0_max": 0.01,
            "relatedness_max_degree_for_independence": 4,
            "threads": 2,
            "king_timeout_seconds": 30,
        },
    }
    config_path.write_text(
        yaml.safe_dump(config, allow_unicode=True, sort_keys=False),
        encoding="utf-8",
    )
    return config_path, runs_dir


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8") as input_file:
        return list(csv.DictReader(input_file, delimiter="\t"))


def test_infer_kinship_accepts_empty_related_pair_set(tmp_path: Path) -> None:
    config_path, runs_dir = prepare_kinship_inputs(tmp_path)

    run_dir = run_pipeline(config_path, runs_dir)

    stage_dir = run_dir / "stages" / "07_infer_kinship"
    assert read_tsv(stage_dir / "kinship_pairs.tsv") == []
    proposal = read_tsv(stage_dir / "independent_set_proposal.tsv")
    assert all(row["PROPOSED_INCLUDED"] == "true" for row in proposal)
    manifest = json.loads((run_dir / "manifest.json").read_text(encoding="utf-8"))
    assert manifest["global_status"] == "TECHNICALLY_VALID"
    assert manifest["manual_decisions_required"] == []


def test_cryptic_first_degree_pair_requires_manual_exclusion_review(
    tmp_path: Path,
) -> None:
    config_path, runs_dir = prepare_kinship_inputs(
        tmp_path,
        between_rows=["F1 I1 F2 I2 22 0.25 0.10 0.25"],
    )

    run_dir = run_pipeline(config_path, runs_dir)

    stage_dir = run_dir / "stages" / "07_infer_kinship"
    pairs = read_tsv(stage_dir / "kinship_pairs.tsv")
    assert pairs[0]["INFERRED_DEGREE"] == "FIRST"
    assert pairs[0]["INFERRED_RELATION"] == "FIRST_DEGREE_NON_PARENT_CHILD"
    assert pairs[0]["CONCORDANCE_STATUS"] == "CRYPTIC_RELATEDNESS"
    proposal = read_tsv(stage_dir / "independent_set_proposal.tsv")
    assert sum(row["PROPOSED_INCLUDED"] == "false" for row in proposal) == 1
    excluded = next(row for row in proposal if row["PROPOSED_INCLUDED"] == "false")
    assert excluded["SAMPLE_ID"] == "sample_1"
    assert excluded["REPRESENTATIVE_SAMPLE_ID"] == "sample_2"
    assert excluded["RELATION_BASIS"] == "KING_RELATEDNESS"
    manifest = json.loads((run_dir / "manifest.json").read_text(encoding="utf-8"))
    assert manifest["manual_decisions_required"] == ["kinship_exclusion_approval"]


def test_declared_parent_child_is_compared_with_king(tmp_path: Path) -> None:
    config_path, runs_dir = prepare_kinship_inputs(
        tmp_path,
        declared_parent=True,
        within_rows=["F1 I1 I2 22 0 0.25 0.25 0 0.25 0"],
    )

    run_dir = run_pipeline(config_path, runs_dir)

    stage_dir = run_dir / "stages" / "07_infer_kinship"
    pedigree = read_tsv(stage_dir / "pedigree_concordance.tsv")
    assert pedigree[0]["PARENT_ROLE"] == "FATHER"
    assert pedigree[0]["INFERRED_RELATION"] == "PARENT_CHILD_COMPATIBLE"
    assert pedigree[0]["CONCORDANCE_STATUS"] == "CONCORDANT"
    proposal = read_tsv(stage_dir / "independent_set_proposal.tsv")
    excluded = next(row for row in proposal if row["PROPOSED_INCLUDED"] == "false")
    assert excluded["RELATION_BASIS"] == "KING_AND_DECLARED"


def test_declared_parent_child_degree_discordance_is_audited(tmp_path: Path) -> None:
    config_path, runs_dir = prepare_kinship_inputs(
        tmp_path,
        declared_parent=True,
        within_rows=["F1 I1 I2 22 0 0.25 0.25 0.10 0.05 1"],
    )

    run_dir = run_pipeline(config_path, runs_dir)

    stage_dir = run_dir / "stages" / "07_infer_kinship"
    pedigree = read_tsv(stage_dir / "pedigree_concordance.tsv")
    assert pedigree[0]["INFERRED_DEGREE"] == "THIRD"
    assert pedigree[0]["CONCORDANCE_STATUS"] == "DISCORDANT_DEGREE"
    audit = json.loads((stage_dir / "audit.json").read_text(encoding="utf-8"))
    assert audit["manual_validation_required"] is True


def test_external_declared_parent_is_not_silently_ignored(tmp_path: Path) -> None:
    config_path, runs_dir = prepare_kinship_inputs(
        tmp_path,
        external_parent=True,
    )

    run_dir = run_pipeline(config_path, runs_dir)

    pedigree = read_tsv(
        run_dir / "stages" / "07_infer_kinship" / "pedigree_concordance.tsv"
    )
    assert pedigree[0]["PARENT_SAMPLE_ID"] == ""
    assert pedigree[0]["CONCORDANCE_STATUS"] == "NOT_EVALUATED_PARENT_ABSENT"


def test_king_failure_uses_external_tool_return_code(tmp_path: Path) -> None:
    config_path, runs_dir = prepare_kinship_inputs(tmp_path, fail_king=True)

    with pytest.raises(StageExecutionError) as error:
        run_pipeline(config_path, runs_dir)

    assert error.value.return_code == 3


def test_invalid_kinship_threshold_order_uses_input_error_code(tmp_path: Path) -> None:
    config_path, runs_dir = prepare_kinship_inputs(tmp_path)
    config = yaml.safe_load(config_path.read_text(encoding="utf-8"))
    config["stages"]["infer_kinship"]["parameters"]["second_degree_min"] = 0.2
    config_path.write_text(
        yaml.safe_dump(config, allow_unicode=True, sort_keys=False),
        encoding="utf-8",
    )

    with pytest.raises(StageExecutionError) as error:
        run_pipeline(config_path, runs_dir)

    assert error.value.return_code == 2


def test_king_degree_boundaries_are_inclusive() -> None:
    parameters = {
        "duplicate_or_mz_min": 0.354,
        "first_degree_min": 0.177,
        "second_degree_min": 0.0884,
        "third_degree_min": 0.0442,
        "fourth_degree_min": 0.0221,
        "parent_child_ibs0_max": 0.01,
    }

    assert _classify_relationship(0.354, 0, parameters)[0] == "DUPLICATE_OR_MZ"
    assert _classify_relationship(0.177, 0, parameters) == (
        "FIRST",
        "PARENT_CHILD_COMPATIBLE",
    )
    assert _classify_relationship(0.0884, 0.1, parameters)[0] == "SECOND"
    assert _classify_relationship(0.0442, 0.1, parameters)[0] == "THIRD"
    assert _classify_relationship(0.0221, 0.1, parameters)[0] == "FOURTH"
    assert _classify_relationship(0.02, 0.1, parameters)[0] == "UNRELATED"
