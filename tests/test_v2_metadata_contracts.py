from pathlib import Path

import pytest

from effet_fondateur.cli import main
from effet_fondateur.contracts import (
    TableValidationError,
    validate_cohorts_frozen,
    validate_samples_master,
)


SAMPLES_HEADER = (
    "SAMPLE_ID\tSOURCE_FILE\tFID\tIID\tPID\tMID\tSEX\tCLINICAL_STATUS\t"
    "GROUP_LABEL\tTARGET_GENOTYPE\tTARGET_GENOTYPE_SOURCE\tARRAY_BATCH\t"
    "INCLUDE_GENOMEWIDE\tINCLUDE_TARGET_CHROMOSOME\tNOTES_CODE\n"
)
COHORTS_HEADER = (
    "COHORT_ID\tSAMPLE_ID\tROLE\tINCLUDED\tEXCLUSION_CODE\tDECISION_SOURCE\t"
    "INDEPENDENT_UNIT_ID\tFAMILY_REPRESENTATIVE\n"
)


def write_valid_samples(table_path: Path, source_root: Path) -> None:
    source_root.mkdir()
    (source_root / "parent.txt").write_text("source technique\n", encoding="utf-8")
    (source_root / "child.txt").write_text("source technique\n", encoding="utf-8")
    table_path.write_text(
        SAMPLES_HEADER
        + "sample_parent\tparent.txt\tFAM1\tPARENT\t0\t0\tMALE\tUNAFFECTED\t"
        "FAMILY\tA/G\tlaboratory_report\tbatch_1\ttrue\ttrue\t\n"
        + "sample_child\tchild.txt\tFAM1\tCHILD\tPARENT\t0\tUNKNOWN\tAFFECTED\t"
        "FAMILY\t\t\tbatch_1\ttrue\ttrue\t\n",
        encoding="utf-8",
    )


def test_samples_master_validates_sources_and_genotype_provenance(tmp_path: Path) -> None:
    table_path = tmp_path / "samples.master.tsv"
    source_root = tmp_path / "sources"
    write_valid_samples(table_path, source_root)

    validated_table = validate_samples_master(table_path, source_root)

    assert validated_table.row_count == 2
    assert validated_table.rows[0]["INCLUDE_GENOMEWIDE"] is True
    assert validated_table.rows[1]["TARGET_GENOTYPE"] is None
    assert len(validated_table.sha256) == 64


def test_samples_master_rejects_clinical_genotype_inference(tmp_path: Path) -> None:
    table_path = tmp_path / "samples.master.tsv"
    source_root = tmp_path / "sources"
    write_valid_samples(table_path, source_root)
    table_path.write_text(
        table_path.read_text(encoding="utf-8").replace(
            "A/G\tlaboratory_report", "A/G\tclinical_status"
        ),
        encoding="utf-8",
    )

    with pytest.raises(TableValidationError, match="statut clinique"):
        validate_samples_master(table_path, source_root)


def test_samples_master_rejects_duplicate_plink_identifier(tmp_path: Path) -> None:
    table_path = tmp_path / "samples.master.tsv"
    source_root = tmp_path / "sources"
    write_valid_samples(table_path, source_root)
    table_path.write_text(
        table_path.read_text(encoding="utf-8").replace(
            "FAM1\tCHILD\tPARENT", "FAM1\tPARENT\tPARENT"
        ),
        encoding="utf-8",
    )

    with pytest.raises(TableValidationError, match=r"FID \+ IID"):
        validate_samples_master(table_path, source_root)


def test_samples_master_rejects_source_assigned_to_multiple_samples(
    tmp_path: Path,
) -> None:
    table_path = tmp_path / "samples.master.tsv"
    source_root = tmp_path / "sources"
    write_valid_samples(table_path, source_root)
    table_path.write_text(
        table_path.read_text(encoding="utf-8").replace("child.txt", "parent.txt"),
        encoding="utf-8",
    )

    with pytest.raises(TableValidationError, match="plusieurs individus"):
        validate_samples_master(table_path, source_root)


def test_samples_master_rejects_missing_source(tmp_path: Path) -> None:
    table_path = tmp_path / "samples.master.tsv"
    source_root = tmp_path / "sources"
    write_valid_samples(table_path, source_root)
    (source_root / "child.txt").unlink()

    with pytest.raises(TableValidationError, match="SOURCE_FILE"):
        validate_samples_master(table_path, source_root)


def test_samples_master_rejects_noncanonical_boolean(tmp_path: Path) -> None:
    table_path = tmp_path / "samples.master.tsv"
    source_root = tmp_path / "sources"
    write_valid_samples(table_path, source_root)
    table_path.write_text(
        table_path.read_text(encoding="utf-8").replace("\ttrue\ttrue\t", "\t1\ttrue\t", 1),
        encoding="utf-8",
    )

    with pytest.raises(TableValidationError, match="INCLUDE_GENOMEWIDE"):
        validate_samples_master(table_path, source_root)


def test_samples_master_rejects_extra_field(tmp_path: Path) -> None:
    table_path = tmp_path / "samples.master.tsv"
    source_root = tmp_path / "sources"
    write_valid_samples(table_path, source_root)
    lines = table_path.read_text(encoding="utf-8").splitlines()
    lines[1] += "\textra"
    table_path.write_text("\n".join(lines) + "\n", encoding="utf-8")

    with pytest.raises(TableValidationError, match="nombre de champs"):
        validate_samples_master(table_path, source_root)


def test_cohorts_validate_membership_and_exclusion_rules(tmp_path: Path) -> None:
    table_path = tmp_path / "cohorts.frozen.tsv"
    table_path.write_text(
        COHORTS_HEADER
        + "controls_v1\tsample_parent\tCONTROL\ttrue\t\tstage_09\tcontrol_1\ttrue\n"
        + "controls_v1\tsample_child\tRELATED\tfalse\tkinship_degree4\tstage_09\t\tfalse\n",
        encoding="utf-8",
    )

    validated_table = validate_cohorts_frozen(
        table_path, {"sample_parent", "sample_child"}
    )

    assert validated_table.row_count == 2
    assert validated_table.rows[1]["INCLUDED"] is False


def test_cohorts_reject_exclusion_without_reason(tmp_path: Path) -> None:
    table_path = tmp_path / "cohorts.frozen.tsv"
    table_path.write_text(
        COHORTS_HEADER
        + "controls_v1\tsample_parent\tCONTROL\tfalse\t\tstage_09\t\tfalse\n",
        encoding="utf-8",
    )

    with pytest.raises(TableValidationError, match="EXCLUSION_CODE"):
        validate_cohorts_frozen(table_path, {"sample_parent"})


def test_cohorts_reject_unknown_sample(tmp_path: Path) -> None:
    table_path = tmp_path / "cohorts.frozen.tsv"
    table_path.write_text(
        COHORTS_HEADER
        + "controls_v1\tunknown_sample\tCONTROL\ttrue\t\tstage_09\tunit_1\ttrue\n",
        encoding="utf-8",
    )

    with pytest.raises(TableValidationError, match="table maître"):
        validate_cohorts_frozen(table_path, {"sample_parent"})


def test_cohorts_reject_duplicate_membership(tmp_path: Path) -> None:
    table_path = tmp_path / "cohorts.frozen.tsv"
    membership = "controls_v1\tsample_parent\tCONTROL\ttrue\t\tstage_09\tunit_1\ttrue\n"
    table_path.write_text(COHORTS_HEADER + membership + membership, encoding="utf-8")

    with pytest.raises(TableValidationError, match="clé primaire"):
        validate_cohorts_frozen(table_path, {"sample_parent"})


def test_validate_samples_cli_reports_only_row_count(
    tmp_path: Path, capsys: pytest.CaptureFixture[str]
) -> None:
    table_path = tmp_path / "samples.master.tsv"
    source_root = tmp_path / "sources"
    write_valid_samples(table_path, source_root)

    exit_code = main(
        [
            "validate-samples",
            "--table",
            str(table_path),
            "--source-root",
            str(source_root),
        ]
    )

    assert exit_code == 0
    assert capsys.readouterr().out == "Table maître valide : 2 lignes\n"
