from effet_fondateur.phasing.inputs import _sample_rows


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
