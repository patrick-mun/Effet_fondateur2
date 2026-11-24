import pytest

from effet_fondateur.phasing.inputs import (
    Shapeit5InputBlockError,
    _mendel_exclusion_rows,
    _sample_rows,
    _variant_selection,
)


def test_shapeit5_pedigree_uses_master_sample_ids_for_trios_and_duos() -> None:
    fam_rows = [
        ["F1", "father", "0", "0", "1", "1"],
        ["F1", "mother", "0", "0", "2", "1"],
        ["F1", "child", "father", "mother", "1", "1"],
        ["F1", "duo", "father", "absent", "2", "1"],
    ]
    master_rows = [
        {"SAMPLE_ID": "sample_father", "FID": "F1", "IID": "father", "PID": "0", "MID": "0"},
        {"SAMPLE_ID": "sample_mother", "FID": "F1", "IID": "mother", "PID": "0", "MID": "0"},
        {"SAMPLE_ID": "sample_child", "FID": "F1", "IID": "child", "PID": "father", "MID": "mother"},
        {"SAMPLE_ID": "sample_duo", "FID": "F1", "IID": "duo", "PID": "father", "MID": "absent"},
    ]

    mapping_rows, pedigree_rows = _sample_rows(fam_rows, master_rows)

    assert pedigree_rows == [
        ("sample_child", "sample_father", "sample_mother"),
        ("sample_duo", "sample_father", "NA"),
    ]
    assert [row["PEDIGREE_INCLUDED"] for row in mapping_rows] == [
        False,
        False,
        True,
        True,
    ]


def test_mendel_policy_excludes_only_incompatible_non_target_variant() -> None:
    samples = ["father", "mother", "child"]
    pedigree = [("child", "father", "mother")]
    records = [
        ("chr19", 100, "probe_bad", ("0/0", "0/0", "1/1")),
        ("chr19", 200, "target", ("0/0", "0/1", "0/1")),
    ]

    rows = _mendel_exclusion_rows(
        records,
        samples,
        pedigree,
        "target",
        "exclude_non_target_variants",
    )

    assert [row["VARIANT_ID"] for row in rows] == ["probe_bad"]
    assert rows[0]["AFFECTED_PEDIGREE_RECORD_COUNT"] == 1
    assert rows[0]["IS_TARGET_VARIANT"] is False


def test_mendel_policy_blocks_without_explicit_exclusion_approval() -> None:
    records = [("chr19", 100, "probe_bad", ("0/0", "0/0", "1/1"))]

    with pytest.raises(Shapeit5InputBlockError, match="mendel_errors_before_phasing"):
        _mendel_exclusion_rows(
            records,
            ["father", "mother", "child"],
            [("child", "father", "mother")],
            "target",
            "block",
        )


def test_mendel_policy_never_excludes_target_variant() -> None:
    records = [("chr19", 100, "target", ("0/0", "0/0", "1/1"))]

    with pytest.raises(Shapeit5InputBlockError, match="target_variant_mendel_error"):
        _mendel_exclusion_rows(
            records,
            ["father", "mother", "child"],
            [("child", "father", "mother")],
            "target",
            "exclude_non_target_variants",
        )


def test_variant_selection_exposes_included_common_and_target_rows() -> None:
    target = {
        "VARIANT_ID": "target",
        "SHAPEIT_ROLE": "RARE_TARGET",
        "INCLUDED": True,
        "IS_TARGET_VARIANT": True,
    }
    common = {
        "VARIANT_ID": "common",
        "SHAPEIT_ROLE": "COMMON",
        "INCLUDED": True,
        "IS_TARGET_VARIANT": False,
    }
    excluded = {
        "VARIANT_ID": "excluded",
        "SHAPEIT_ROLE": "EXCLUDED",
        "INCLUDED": False,
        "IS_TARGET_VARIANT": False,
    }

    selection = _variant_selection([target, common, excluded])

    assert selection.included_rows == [target, common]
    assert selection.common_rows == [common]
    assert selection.target_row is target


def test_variant_selection_uses_post_mendel_error_code() -> None:
    rows = [
        {
            "SHAPEIT_ROLE": "RARE_TARGET",
            "INCLUDED": True,
            "IS_TARGET_VARIANT": True,
        },
        {
            "SHAPEIT_ROLE": "EXCLUDED",
            "INCLUDED": False,
            "IS_TARGET_VARIANT": False,
        },
    ]

    with pytest.raises(
        Shapeit5InputBlockError,
        match="insufficient_variants_after_mendel_exclusions",
    ):
        _variant_selection(
            rows,
            insufficient_variants_error=(
                "insufficient_variants_after_mendel_exclusions"
            ),
        )
