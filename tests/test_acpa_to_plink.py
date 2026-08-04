import json
from pathlib import Path

import pytest

from simulation_genotype_famille.acpa_to_plink import (
    ConversionError,
    convert_acpa_to_plink,
    create_metadata_template,
)


ACPA_HEADER = """# Array Type Name: CytoScan 750K Accel Array
# UCSC Genomic Version: hg38
Probe Set ID\tCall Codes\tConfidence\tSignal A\tSignal B\tForward Strand Base Calls\tdbSNP RS ID\tChromosome\tChromosomal Position\t
"""


def write_acpa_sample(path: Path, rows: list[str]) -> None:
    path.write_text(ACPA_HEADER + "\n".join(rows) + "\n", encoding="utf-8")


def test_convert_acpa_to_plink_generates_expected_files(tmp_path):
    input_dir = tmp_path / "acpa_samples"
    input_dir.mkdir()
    write_acpa_sample(
        input_dir / "patient.txt",
        [
            "probe_1\tAB\t0\t1\t1\tAG\t\t19\t100",
            "probe_2\tAA\t0\t1\t1\tCC\trs2\t19\t200",
            "probe_chr1\tAA\t0\t1\t1\tAA\trs1\t1\t50",
        ],
    )
    write_acpa_sample(
        input_dir / "temoin.txt",
        [
            "probe_1\tAA\t0\t1\t1\tAA\t\t19\t100",
            "probe_2\tAB\t0\t1\t1\tCT\trs2\t19\t200",
        ],
    )
    metadata_path = tmp_path / "samples.tsv"
    metadata_path.write_text(
        "FILE\tFID\tIID\tPID\tMID\tSEX\tPHENOTYPE\tGROUP\n"
        "patient.txt\tF1\tPATIENT\t0\t0\t1\t2\tATTEINT\n"
        "temoin.txt\tCTRL\tTEMOIN\t0\t0\t2\t1\tTEMOIN\n",
        encoding="utf-8",
    )
    output_prefix = tmp_path / "plink" / "genotype_data"

    report = convert_acpa_to_plink(input_dir, metadata_path, output_prefix)

    assert report["sample_count"] == 2
    assert report["output_marker_count"] == 2
    assert output_prefix.with_suffix(".map").read_text().splitlines() == [
        "19\tprobe_1\t0\t100",
        "19\tprobe_2\t0\t200",
    ]
    assert output_prefix.with_suffix(".ped").read_text().splitlines() == [
        "F1\tPATIENT\t0\t0\t1\t2\tA\tG\tC\tC",
        "CTRL\tTEMOIN\t0\t0\t2\t1\tA\tA\tC\tT",
    ]
    assert (output_prefix.parent / "groupes.txt").read_text().splitlines() == [
        "PATIENT\tATTEINT",
        "TEMOIN\tTEMOIN",
    ]
    assert (output_prefix.parent / "cas.txt").read_text().strip() == "F1\tPATIENT"
    assert (output_prefix.parent / "temoins.txt").read_text().strip() == "CTRL\tTEMOIN"
    saved_report = json.loads(
        (output_prefix.parent / "acpa_conversion_report.json").read_text()
    )
    assert saved_report["genomic_builds"] == ["hg38"]


def test_converter_excludes_multiallelic_marker(tmp_path):
    input_dir = tmp_path / "samples"
    input_dir.mkdir()
    write_acpa_sample(
        input_dir / "sample_1.txt",
        [
            "probe_bad\tAB\t0\t1\t1\tAG\trs_bad\t19\t100",
            "probe_ok\tAA\t0\t1\t1\tCC\trs_ok\t19\t200",
        ],
    )
    write_acpa_sample(
        input_dir / "sample_2.txt",
        [
            "probe_bad\tAB\t0\t1\t1\tCT\trs_bad\t19\t100",
            "probe_ok\tAA\t0\t1\t1\tCC\trs_ok\t19\t200",
        ],
    )
    metadata_path = tmp_path / "metadata.tsv"
    metadata_path.write_text(
        "FILE\tFID\tIID\tPID\tMID\tSEX\tPHENOTYPE\tGROUP\n"
        "sample_1.txt\tF1\tI1\t0\t0\t1\t2\tATTEINT\n"
        "sample_2.txt\tF2\tI2\t0\t0\t2\t1\tTEMOIN\n",
        encoding="utf-8",
    )

    report = convert_acpa_to_plink(
        input_dir,
        metadata_path,
        tmp_path / "output" / "genotype_data",
    )

    assert report["output_marker_count"] == 1
    assert report["excluded_marker_count"] == 1
    excluded_report = (tmp_path / "output" / "acpa_excluded_markers.tsv").read_text()
    assert "probe_bad\tplus_de_deux_alleles\tA,C,G,T" in excluded_report


def test_create_metadata_template_and_prevent_overwrite(tmp_path):
    input_dir = tmp_path / "samples"
    input_dir.mkdir()
    write_acpa_sample(
        input_dir / "PATIENT_01_(CytoScan750K_Array).txt",
        ["probe_1\tAA\t0\t1\t1\tAA\trs1\t19\t100"],
    )
    template_path = tmp_path / "samples.tsv"

    create_metadata_template(input_dir, template_path)

    template_lines = template_path.read_text().splitlines()
    assert template_lines[1].startswith(
        "PATIENT_01_(CytoScan750K_Array).txt\tPATIENT_01\tPATIENT_01"
    )
    with pytest.raises(ConversionError, match="Sorties déjà présentes"):
        create_metadata_template(input_dir, template_path)


def test_converter_rejects_plink_identifier_with_whitespace(tmp_path):
    input_dir = tmp_path / "samples"
    input_dir.mkdir()
    write_acpa_sample(
        input_dir / "sample.txt",
        ["probe_1\tAA\t0\t1\t1\tAA\trs1\t19\t100"],
    )
    metadata_path = tmp_path / "metadata.tsv"
    metadata_path.write_text(
        "FILE\tFID\tIID\tPID\tMID\tSEX\tPHENOTYPE\tGROUP\n"
        "sample.txt\tF 1\tI1\t0\t0\t1\t2\tATTEINT\n",
        encoding="utf-8",
    )

    with pytest.raises(ConversionError, match="Identifiant PLINK contenant un espace"):
        convert_acpa_to_plink(
            input_dir,
            metadata_path,
            tmp_path / "output" / "genotype_data",
        )
