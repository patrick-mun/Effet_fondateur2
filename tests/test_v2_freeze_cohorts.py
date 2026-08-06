import csv
import json
import sys
from pathlib import Path

import pytest
import yaml

from effet_fondateur.orchestrator import StageExecutionError
from effet_fondateur.orchestrator.pipeline import run_pipeline
from test_v2_population_structure import prepare_population_inputs


GENOTYPE_HEADER = "SAMPLE_ID\tGENOTYPE\tGENOTYPE_SOURCE\tMEASUREMENT_METHOD\n"


def write_cohort_plink(path: Path, delegated_plink: Path) -> None:
    script = f'''#!{sys.executable}
import os
import pathlib
import shutil
import sys

if "--bmerge" in sys.argv:
    output_prefix = pathlib.Path(sys.argv[sys.argv.index("--out") + 1])
    base_prefix = pathlib.Path(sys.argv[sys.argv.index("--bfile") + 1])
    merge_bim = pathlib.Path(sys.argv[sys.argv.index("--bmerge") + 2])
    shutil.copyfile(base_prefix.with_suffix(".fam"), output_prefix.with_suffix(".fam"))
    output_prefix.with_suffix(".bim").write_text(base_prefix.with_suffix(".bim").read_text() + merge_bim.read_text())
    output_prefix.with_suffix(".bed").write_bytes(b"\\x6c\\x1b\\x01")
    raise SystemExit(0)
if "--mendel" in sys.argv:
    output_prefix = pathlib.Path(sys.argv[sys.argv.index("--out") + 1])
    output_prefix.with_suffix(".mendel").write_text("FID KID CHR SNP CODE ERROR\\n")
    raise SystemExit(0)
os.execv({str(delegated_plink)!r}, [{str(delegated_plink)!r}, *sys.argv[1:]])
'''
    path.write_text(script, encoding="utf-8")
    path.chmod(0o755)


def write_target_inputs(
    config: dict[str, object],
    root: Path,
    *,
    carrier_sample_id: str,
) -> None:
    target = config["target"]
    variant_path = root / "target_variant.yaml"
    variant_path.write_text(
        yaml.safe_dump(
            {
                "schema_version": "1.0.0",
                "assembly": config["project"]["assembly"],
                "gene": target["gene"],
                "chromosome": target["chromosome"],
                "position_bp": target["position_bp"],
                "ref": target["ref"],
                "alt": target["alt"],
                "transcript": target["transcript"],
                "hgvs_c": "c.100A>G",
                "hgvs_p": "p.Lys34Arg",
                "project_variant_id": target["project_variant_id"],
                "rsid": None,
                "reference_validation": {
                    "status": "CONFIRMED",
                    "source_id": "synthetic_reference",
                    "confirmed_at": "2026-08-05T10:00:00Z",
                },
                "annotation_validation": {
                    "status": "CONFIRMED",
                    "source_id": "synthetic_annotation",
                    "confirmed_at": "2026-08-05T10:00:00Z",
                },
            },
            sort_keys=False,
        ),
        encoding="utf-8",
    )
    genotype_path = root / "target_genotypes.tsv"
    genotype_path.write_text(
        GENOTYPE_HEADER
        + "".join(
            f"sample_{sample_index}\t{'A/G' if f'sample_{sample_index}' == carrier_sample_id else 'A/A'}\t"
            "laboratory_report\tvalidated_assay\n"
            for sample_index in range(1, 4)
        ),
        encoding="utf-8",
    )
    config["inputs"]["target_variant_metadata"] = str(variant_path)
    config["inputs"]["target_variant_genotypes"] = str(genotype_path)


def prepare_cohort_inputs(
    tmp_path: Path,
    *,
    related_pair: bool = False,
    population_outliers: bool = False,
    carrier_sample_id: str = "sample_1",
    enable_freeze: bool = True,
    decision_reviews: dict[str, object] | None = None,
    required_cohorts: list[str] | None = None,
) -> tuple[Path, Path]:
    config_path, runs_dir = prepare_population_inputs(
        tmp_path,
        related_pair=related_pair,
        parameter_overrides={"outlier_alpha": 1.0 if population_outliers else 1e-12},
    )
    config = yaml.safe_load(config_path.read_text(encoding="utf-8"))
    write_target_inputs(config, tmp_path, carrier_sample_id=carrier_sample_id)
    delegated_plink = Path(config["tools"]["plink"])
    cohort_plink = tmp_path / "cohort_plink"
    write_cohort_plink(cohort_plink, delegated_plink)
    config["tools"]["plink"] = str(cohort_plink)
    config["stages"]["prepare_target_variant_dataset"] = {
        "enabled": True,
        "parameters": {"plink_timeout_seconds": 30},
    }
    config["stages"]["freeze_cohorts"] = {
        "enabled": enable_freeze,
        "parameters": {
            "control_group_labels": ["FAMILY"],
            "required_nonempty_cohorts": required_cohorts
            or [
                "controls_unrelated",
                "target_carriers_all",
                "target_carriers_independent",
                "target_chromosome_all_qc",
                "genomewide_structure_reference",
            ],
            "decision_reviews": decision_reviews or {},
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


def artifact_sha256(run_dir: Path, stage_directory: str, artifact_id: str) -> str:
    stage_outputs = json.loads(
        (run_dir / "stages" / stage_directory / "stage_outputs.json").read_text(
            encoding="utf-8"
        )
    )
    return next(
        artifact["sha256"]
        for artifact in stage_outputs["artifacts"]
        if artifact["artifact_id"] == artifact_id
    )


def review(artifact_sha: str, sample_ids: list[str], review_id: str) -> dict[str, object]:
    return {
        "review_id": review_id,
        "reviewer_role": "scientific_reviewer",
        "reviewed_at": "2026-08-05T11:00:00Z",
        "reviewed_artifact_sha256": artifact_sha,
        "approved_exclusion_sample_ids": sample_ids,
    }


def test_freeze_cohorts_publishes_explicit_memberships_and_keep_files(
    tmp_path: Path,
) -> None:
    config_path, runs_dir = prepare_cohort_inputs(tmp_path)

    run_dir = run_pipeline(config_path, runs_dir)

    stage_dir = run_dir / "stages" / "09_freeze_cohorts"
    cohorts = read_tsv(stage_dir / "cohorts.frozen.tsv")
    assert len(cohorts) == 21
    carrier = next(
        row
        for row in cohorts
        if row["COHORT_ID"] == "target_carriers_all"
        and row["SAMPLE_ID"] == "sample_1"
    )
    assert carrier["INCLUDED"] == "true"
    assert carrier["ROLE"] == "CARRIER"
    affected_noncarrier = next(
        row
        for row in cohorts
        if row["COHORT_ID"] == "affected_noncarriers"
        and row["SAMPLE_ID"] == "sample_1"
    )
    assert affected_noncarrier["INCLUDED"] == "false"
    controls_keep = (stage_dir / "keep" / "controls_unrelated.keep").read_text(
        encoding="utf-8"
    )
    assert "F2\tI2" in controls_keep
    assert "F3\tI3" in controls_keep
    assert "F1\tI1" not in controls_keep
    assert read_tsv(stage_dir / "cohort_decisions.tsv") == []
    manifest = json.loads((run_dir / "manifest.json").read_text(encoding="utf-8"))
    assert manifest["global_status"] == "TECHNICALLY_VALID"


def test_unreviewed_kinship_exclusion_blocks_freeze(tmp_path: Path) -> None:
    config_path, runs_dir = prepare_cohort_inputs(
        tmp_path,
        related_pair=True,
        carrier_sample_id="sample_2",
    )

    with pytest.raises(StageExecutionError) as error:
        run_pipeline(config_path, runs_dir)

    assert error.value.return_code == 4
    manifests = list(runs_dir.glob("*/manifest.json"))
    assert len(manifests) == 1
    manifest = json.loads(manifests[0].read_text(encoding="utf-8"))
    assert "kinship_exclusion_approval" in manifest["manual_decisions_required"]


def test_reviewed_kinship_exclusion_is_applied_and_resolved(tmp_path: Path) -> None:
    first_config, runs_dir = prepare_cohort_inputs(
        tmp_path,
        related_pair=True,
        carrier_sample_id="sample_2",
        enable_freeze=False,
    )
    first_run = run_pipeline(first_config, runs_dir)
    proposal_sha = artifact_sha256(
        first_run, "07_infer_kinship", "independent_set_proposal"
    )
    second_config = yaml.safe_load(first_config.read_text(encoding="utf-8"))
    second_config["stages"]["freeze_cohorts"] = {
        "enabled": True,
        "parameters": {
            "control_group_labels": ["FAMILY"],
            "required_nonempty_cohorts": [
                "controls_unrelated",
                "target_carriers_all",
                "target_carriers_independent",
                "target_chromosome_all_qc",
                "genomewide_structure_reference",
            ],
            "decision_reviews": {
                "kinship_exclusion_approval": review(
                    proposal_sha, ["sample_1"], "kinship_review_001"
                )
            },
        },
    }
    reviewed_config_path = tmp_path / "reviewed_config.yaml"
    reviewed_config_path.write_text(
        yaml.safe_dump(second_config, allow_unicode=True, sort_keys=False),
        encoding="utf-8",
    )

    second_run = run_pipeline(reviewed_config_path, runs_dir)

    decisions = read_tsv(
        second_run / "stages" / "09_freeze_cohorts" / "cohort_decisions.tsv"
    )
    assert len(decisions) == 3
    assert {row["SAMPLE_ID"] for row in decisions} == {"sample_1"}
    manifest = json.loads((second_run / "manifest.json").read_text(encoding="utf-8"))
    assert "kinship_exclusion_approval" not in manifest["manual_decisions_required"]
    assert manifest["global_status"] == "TECHNICALLY_VALID"
    audit = json.loads(
        (second_run / "stages" / "09_freeze_cohorts" / "audit.json").read_text(
            encoding="utf-8"
        )
    )
    assert audit["parameters"]["decision_reviews"]["kinship_exclusion_approval"][
        "approved_exclusion_count"
    ] == 1
    assert "approved_exclusion_sample_ids" not in audit["parameters"]["decision_reviews"][
        "kinship_exclusion_approval"
    ]


def test_review_with_wrong_artifact_hash_blocks_freeze(tmp_path: Path) -> None:
    config_path, runs_dir = prepare_cohort_inputs(
        tmp_path,
        related_pair=True,
        carrier_sample_id="sample_2",
        decision_reviews={
            "kinship_exclusion_approval": review(
                "0" * 64, ["sample_1"], "kinship_review_wrong_hash"
            )
        },
    )

    with pytest.raises(StageExecutionError) as error:
        run_pipeline(config_path, runs_dir)

    assert error.value.return_code == 4


def test_reviewed_population_outliers_are_resolved(tmp_path: Path) -> None:
    first_config, runs_dir = prepare_cohort_inputs(
        tmp_path,
        population_outliers=True,
        enable_freeze=False,
    )
    first_run = run_pipeline(first_config, runs_dir)
    outliers = read_tsv(
        first_run
        / "stages"
        / "08_analyze_population_structure"
        / "population_outliers.tsv"
    )
    proposed_ids = sorted(
        row["SAMPLE_ID"]
        for row in outliers
        if row["PROPOSED_EXCLUSION"] == "true"
    )
    assert proposed_ids
    outlier_sha = artifact_sha256(
        first_run, "08_analyze_population_structure", "population_outliers"
    )
    second_config = yaml.safe_load(first_config.read_text(encoding="utf-8"))
    second_config["stages"]["freeze_cohorts"] = {
        "enabled": True,
        "parameters": {
            "control_group_labels": ["FAMILY"],
            "required_nonempty_cohorts": [
                "target_carriers_all",
                "target_chromosome_all_qc",
            ],
            "decision_reviews": {
                "population_structure_exclusion_approval": review(
                    outlier_sha, proposed_ids, "population_review_001"
                )
            },
        },
    }
    reviewed_config_path = tmp_path / "population_reviewed_config.yaml"
    reviewed_config_path.write_text(
        yaml.safe_dump(second_config, allow_unicode=True, sort_keys=False),
        encoding="utf-8",
    )

    second_run = run_pipeline(reviewed_config_path, runs_dir)

    manifest = json.loads((second_run / "manifest.json").read_text(encoding="utf-8"))
    assert "population_structure_exclusion_approval" not in manifest["manual_decisions_required"]
    assert manifest["global_status"] == "TECHNICALLY_VALID"


def test_required_empty_cohort_is_scientific_block(tmp_path: Path) -> None:
    config_path, runs_dir = prepare_cohort_inputs(
        tmp_path,
        required_cohorts=["family_noncarriers"],
    )

    with pytest.raises(StageExecutionError) as error:
        run_pipeline(config_path, runs_dir)

    assert error.value.return_code == 4


def test_invalid_control_group_configuration_uses_input_code(tmp_path: Path) -> None:
    config_path, runs_dir = prepare_cohort_inputs(tmp_path)
    config = yaml.safe_load(config_path.read_text(encoding="utf-8"))
    config["stages"]["freeze_cohorts"]["parameters"]["control_group_labels"] = []
    config_path.write_text(
        yaml.safe_dump(config, allow_unicode=True, sort_keys=False),
        encoding="utf-8",
    )

    with pytest.raises(StageExecutionError) as error:
        run_pipeline(config_path, runs_dir)

    assert error.value.return_code == 2
