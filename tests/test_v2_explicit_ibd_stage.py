import subprocess

from effet_fondateur.stages.call_explicit_ibd import _query_complete_variants


def test_missing_or_unphased_genotypes_are_excluded_from_hap_ibd(monkeypatch, tmp_path):
    output = "v1\t19\t100\t0|0\t0|1\nv2\t19\t200\t0|.\t0|0\nv3\t19\t300\t0/1\t0|0\n"
    monkeypatch.setattr(subprocess, "run", lambda *args, **kwargs: subprocess.CompletedProcess(args[0], 0, output, ""))
    rows = _query_complete_variants("bcftools", tmp_path / "input.bcf", 10)
    assert [row[3] for row in rows] == [True, False, False]
