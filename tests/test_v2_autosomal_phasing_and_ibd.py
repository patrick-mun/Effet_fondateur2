from pathlib import Path

import pytest

from effet_fondateur.explicit_ibd.analysis import IbdSegment
from effet_fondateur.stages import call_explicit_ibd
from effet_fondateur.stages.call_explicit_ibd import (
    _align_mutant_haplotypes,
    _control_frequency,
    _map,
)
from effet_fondateur.stages.phase_autosomal_panel import (
    AutosomalPhasingInputError,
    _interpolate,
    _fill_shapeit5_info_tags_command,
    _map_points,
    _md5_file,
    _plink_bgz_vcf_path,
    _shapeit5_thread_arguments,
    _validated_local_reference_source,
)


def test_local_autosomal_reference_requires_matching_vcf_and_index(tmp_path: Path):
    vcf = tmp_path / "panel.vcf.gz"
    index = tmp_path / "panel.vcf.gz.tbi"
    vcf.write_bytes(b"reference-vcf")
    index.write_bytes(b"reference-index")

    assert _validated_local_reference_source(
        tmp_path, "panel.vcf.gz", _md5_file(vcf), _md5_file(index)
    ) == vcf.resolve()
    with pytest.raises(AutosomalPhasingInputError, match="checksum_mismatch"):
        _validated_local_reference_source(
            tmp_path, "panel.vcf.gz", "0" * 32, _md5_file(index)
        )


def test_plink_bgz_vcf_path_preserves_dotted_prefix(tmp_path: Path):
    assert _plink_bgz_vcf_path(tmp_path / "chr1.study") == tmp_path / "chr1.study.vcf.gz"


def test_shapeit5_input_command_fills_ac_and_an_tags(tmp_path: Path):
    command = _fill_shapeit5_info_tags_command(
        "bcftools", tmp_path / "renamed.vcf.gz", tmp_path / "study.vcf.gz"
    )
    assert command[:2] == ["bcftools", "+fill-tags"]
    assert command[-2:] == ["-t", "AC,AN"]


def test_shapeit5_single_thread_does_not_enable_multithread_mode():
    assert _shapeit5_thread_arguments(1) == []
    assert _shapeit5_thread_arguments(4) == ["--thread", "4"]


def test_ibd_map_requires_vcf_compatible_chr_label_and_is_monotonic(tmp_path: Path):
    path = tmp_path / "chr19.map"
    path.write_text("chr19\tv1\t1.25\t100\nchr19\tv2\t1.75\t200\n", encoding="utf-8")
    cm_at_bp, positions, rows = _map(path, 19)
    assert cm_at_bp == {100: 1.25, 200: 1.75}
    assert positions == (100, 200)
    assert [row["CHROMOSOME"] for row in rows] == ["19", "19"]


def test_interpolated_map_uses_bp_fraction_and_flat_extrapolation(tmp_path: Path):
    path = tmp_path / "map.tsv"
    path.write_text(
        "MAP_ID\tASSEMBLY\tCHROMOSOME\tPOSITION_BP\tPOSITION_CM\n"
        "m\tGRCh38\t2\t100\t1\n"
        "m\tGRCh38\t2\t300\t3\n",
        encoding="utf-8",
    )
    positions, cms = _map_points(path, 2)
    assert _interpolate(200, positions, cms) == pytest.approx(2.0)
    assert _interpolate(50, positions, cms) == 1.0
    assert _interpolate(400, positions, cms) == 3.0


def test_regional_h1_is_swapped_when_whole_chromosome_labels_are_reversed(monkeypatch, tmp_path: Path):
    samples = ["C1"]
    regional = {
        ("19", position, "A", "G"): ("0|1",)
        for position in range(100, 106)
    }
    whole = {
        ("19", position, "A", "G"): ("1|0",)
        for position in range(100, 106)
    }
    calls = iter(((samples, regional), (samples, whole)))
    monkeypatch.setattr(call_explicit_ibd, "_query_phased", lambda *args: next(calls))
    carriers = [{"SAMPLE_ID": "C1", "RELIABILITY_STATUS": "PASS", "ALT_COPY_COUNT": "1", "CARRIER_HAPLOTYPE": "H1"}]
    with pytest.raises(call_explicit_ibd.ExplicitIbdInputError, match="insufficient_aligned"):
        _align_mutant_haplotypes("bcftools", tmp_path / "regional.bcf", tmp_path / "whole.bcf", carriers, {"C1": "F1"}, 10)
    # L'orientation a été calculée avant le garde-fou des trois familles.
    calls = iter(((samples, regional), (samples, whole)))
    monkeypatch.setattr(call_explicit_ibd, "_query_phased", lambda *args: next(calls))
    tripled = [
        {**carriers[0], "SAMPLE_ID": sample}
        for sample in ("C1", "C2", "C3")
    ]
    samples3 = ["C1", "C2", "C3"]
    regional3 = {key: ("0|1", "0|1", "0|1") for key in regional}
    whole3 = {key: ("1|0", "1|0", "1|0") for key in whole}
    calls = iter(((samples3, regional3), (samples3, whole3)))
    monkeypatch.setattr(call_explicit_ibd, "_query_phased", lambda *args: next(calls))
    mutant, audit = _align_mutant_haplotypes("bcftools", tmp_path / "regional.bcf", tmp_path / "whole.bcf", tripled, {"C1": "F1", "C2": "F2", "C3": "F3"}, 10)
    assert mutant == {"F1": {"H2"}, "F2": {"H2"}, "F3": {"H2"}}
    assert {row["ORIENTATION"] for row in audit} == {"SWAPPED"}


def test_homozygous_alternate_carriers_do_not_require_phase_orientation(monkeypatch, tmp_path: Path):
    samples = ["C1", "C2", "C3"]
    regional = {("19", 100, "A", "G"): ("0|0", "0|0", "0|0")}
    whole = {("19", 100, "A", "G"): ("0|0", "0|0", "0|0")}
    calls = iter(((samples, regional), (samples, whole)))
    monkeypatch.setattr(call_explicit_ibd, "_query_phased", lambda *args: next(calls))
    carriers = [
        {
            "SAMPLE_ID": sample,
            "RELIABILITY_STATUS": "PASS",
            "ALT_COPY_COUNT": "2",
            "CARRIER_HAPLOTYPE": "BOTH",
        }
        for sample in samples
    ]

    mutant, audit = _align_mutant_haplotypes(
        "bcftools",
        tmp_path / "regional.bcf",
        tmp_path / "whole.bcf",
        carriers,
        {"C1": "F1", "C2": "F2", "C3": "F3"},
        10,
    )

    assert mutant == {"F1": {"H1", "H2"}, "F2": {"H1", "H2"}, "F3": {"H1", "H2"}}
    assert {row["ORIENTATION"] for row in audit} == {"NOT_REQUIRED_HOMOZYGOUS"}
    assert {row["CONCORDANCE"] for row in audit} == {None}


def test_control_frequency_is_computed_from_both_tools_at_target():
    segments = [
        IbdSegment(tool, "primary", "C1", "H1", "C1", "C2", "H2", "C2", "19", 100, 300, 1.0, 3.0, 120)
        for tool in ("HAP_IBD", "REFINED_IBD")
    ]
    evaluable, positive, frequency = _control_frequency(segments, "primary", {"C1", "C2", "C3"}, 200)
    assert (evaluable, positive) == (3, 1)
    assert frequency == pytest.approx(1 / 3)
