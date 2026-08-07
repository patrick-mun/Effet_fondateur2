from __future__ import annotations

from pathlib import Path
import shutil
import subprocess

import pytest

from effet_fondateur.contracts import validate_tsv_table
from effet_fondateur.orchestrator.pipeline import ANALYZE_ROH_STAGE
from effet_fondateur.roh import publish_roh
from effet_fondateur.stages.analyze_roh import AnalyzeRohInputError, _parameters


def _pairs(start: int, count: int) -> list[tuple[str, str]]:
    return [(f"family_{index}", f"sample_{index}") for index in range(start, start + count)]


def _write_native(prefix: Path, pairs: list[tuple[str, str]], segments: list[str]) -> None:
    prefix.with_suffix(".hom").write_text(
        "FID IID PHE CHR SNP1 SNP2 POS1 POS2 KB NSNP DENSITY PHOM PHET\n"
        + "".join(segments),
        encoding="utf-8",
    )
    prefix.with_suffix(".hom.indiv").write_text(
        "FID IID PHE NSEG KB KBAVG\n"
        + "".join(f"{fid} {iid} -9 0 0 0\n" for fid, iid in pairs),
        encoding="utf-8",
    )
    prefix.with_suffix(".hom.summary").write_text(
        "CHR SNP BP AFF UNAFF\n1 marker 100 0 0\n",
        encoding="utf-8",
    )


def test_roh_separates_genomewide_burden_and_individual_target_intersection(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch,
) -> None:
    controls = _pairs(1, 5)
    carriers = _pairs(6, 5)
    target_pairs = [controls[0], carriers[0], ("family_11", "sample_11")]

    def fake_run(_executable: str, arguments: list[str], _timeout: int) -> None:
        prefix = Path(arguments[arguments.index("--out") + 1])
        if prefix.name == "controls_unrelated_genomewide":
            _write_native(prefix, controls, [
                "family_1 sample_1 -9 1 v1 v2 100 2000000 1999.901 50 39.998 1 0\n"
            ])
        elif prefix.name == "target_carriers_independent_genomewide":
            _write_native(prefix, carriers, [])
        else:
            _write_native(prefix, target_pairs, [
                "family_1 sample_1 -9 19 t1 t3 50 250 0.201 3 0.067 1 0\n",
                "family_6 sample_6 -9 19 t1 t3 50 250 0.201 3 0.067 0.9 0.1\n",
            ])

    monkeypatch.setattr("effet_fondateur.roh.analysis._run_plink", fake_run)
    all_pairs = controls + carriers + [("family_11", "sample_11")]
    sample_by_pair = {pair: pair[1] for pair in all_pairs}
    parameters = {
        **_parameters({}),
        "minimum_roh_kb": 0.1, "minimum_roh_snps": 2,
        "window_snps": 2, "minimum_genomewide_variants": 3,
        "minimum_genomewide_autosomes": 1,
        "minimum_target_chromosome_variants": 3,
    }

    publication = publish_roh(
        plink_executable="plink",
        genomewide_datasets={
            "controls_unrelated": {"dataset_id": "controls_unrelated_qc", "prefix": tmp_path / "controls", "fam_pairs": controls},
            "target_carriers_independent": {"dataset_id": "target_carriers_independent_qc", "prefix": tmp_path / "carriers", "fam_pairs": carriers},
        },
        target_dataset={"dataset_id": "target_chromosome_all_qc", "prefix": tmp_path / "target", "fam_pairs": target_pairs},
        common_variant_ids=["v1", "v2", "v3"], common_autosome_count=1,
        target_variant_count=3, sample_by_pair=sample_by_pair,
        target_genotypes={"sample_1": "A/A", "sample_6": "A/G", "sample_11": "A/A"},
        roles_by_sample={
            "sample_1": {"controls_unrelated"},
            "sample_6": {"target_carriers_independent"},
            "sample_11": {"family_noncarriers"},
        },
        target_variant_id="target", target_chromosome="19", target_bp=100,
        output_dir=tmp_path / "output", parameters=parameters,
    )

    segments = validate_tsv_table(publication.segments_path, "roh_segments.schema.json").rows
    assert len(segments) == 3
    assert sum(row["INVOLVES_TARGET"] for row in segments) == 2
    burden = validate_tsv_table(publication.burden_path, "roh_burden.schema.json").rows
    assert len(burden) == 10
    assert next(row for row in burden if row["SAMPLE_ID"] == "sample_1")["N_ROH"] == "1"
    assert next(row for row in burden if row["SAMPLE_ID"] == "sample_6")["N_ROH"] == "0"
    assert all(row["F_ROH"] is None for row in burden)
    target = validate_tsv_table(publication.target_path, "target_in_roh.schema.json").rows
    assert next(row for row in target if row["SAMPLE_ID"] == "sample_1")["INTERPRETATION_STATUS"] == "IN_ROH_HOMOZYGOUS"
    assert next(row for row in target if row["SAMPLE_ID"] == "sample_6")["INTERPRETATION_STATUS"] == "IN_ROH_HETEROZYGOUS_WARNING"
    assert next(row for row in target if row["SAMPLE_ID"] == "sample_11")["INTERPRETATION_STATUS"] == "NOT_IN_ROH"
    assert publication.target_in_roh_count == 2
    assert len(publication.native_paths) == 9


def test_roh_publishes_not_evaluated_without_relaxing_density_thresholds(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch,
) -> None:
    def fail_if_called(*_args, **_kwargs) -> None:
        raise AssertionError("PLINK must not run when both scopes fail preflight")

    monkeypatch.setattr("effet_fondateur.roh.analysis._run_plink", fail_if_called)
    control = [("family_1", "sample_1")]
    carrier = [("family_2", "sample_2")]
    parameters = _parameters({})
    publication = publish_roh(
        plink_executable="plink",
        genomewide_datasets={
            "controls_unrelated": {"dataset_id": "controls_unrelated_qc", "prefix": tmp_path / "controls", "fam_pairs": control},
            "target_carriers_independent": {"dataset_id": "target_carriers_independent_qc", "prefix": tmp_path / "carriers", "fam_pairs": carrier},
        },
        target_dataset={"dataset_id": "target_chromosome_all_qc", "prefix": tmp_path / "target", "fam_pairs": control + carrier},
        common_variant_ids=["v1"], common_autosome_count=1,
        target_variant_count=1,
        sample_by_pair={control[0]: "sample_1", carrier[0]: "sample_2"},
        target_genotypes={"sample_1": "A/A", "sample_2": "A/G"},
        roles_by_sample={"sample_1": {"controls_unrelated"}, "sample_2": {"target_carriers_independent"}},
        target_variant_id="target", target_chromosome="19", target_bp=100,
        output_dir=tmp_path / "output", parameters=parameters,
    )

    assert set(publication.scope_statuses.values()) == {"NOT_EVALUATED_DENSITY_OR_COVERAGE"}
    burden = validate_tsv_table(publication.burden_path, "roh_burden.schema.json").rows
    assert all(row["ANALYSIS_STATUS"] == "NOT_EVALUATED_SCOPE" for row in burden)
    target = validate_tsv_table(publication.target_path, "target_in_roh.schema.json").rows
    assert all(row["INTERPRETATION_STATUS"] == "NOT_EVALUATED_SCOPE" for row in target)


def test_roh_requires_a_sourced_denominator_for_froh() -> None:
    with pytest.raises(AnalyzeRohInputError, match="requires_value_and_source"):
        _parameters({"autosomal_denominator_kb": 2_800_000})


def test_roh_stage_is_secondary_and_does_not_depend_on_13_15() -> None:
    assert ANALYZE_ROH_STAGE.critical is False
    assert "infer_founder_haplotype" not in ANALYZE_ROH_STAGE.dependencies
    assert "estimate_variant_age" not in ANALYZE_ROH_STAGE.dependencies
    assert "analyze_local_ld" not in ANALYZE_ROH_STAGE.dependencies


@pytest.mark.skipif(shutil.which("plink") is None, reason="PLINK 1.9 absent")
def test_roh_smoke_with_real_plink_on_temporary_synthetic_data(tmp_path: Path) -> None:
    marker_ids = [f"v{index}" for index in range(1, 11)]
    (tmp_path / "base.map").write_text(
        "".join(f"1 {variant_id} 0 {index * 1000}\n" for index, variant_id in enumerate(marker_ids, start=1)),
        encoding="utf-8",
    )
    all_pairs = _pairs(1, 10)
    genotype_calls = " ".join("A A" for _ in marker_ids)
    (tmp_path / "base.ped").write_text(
        "".join(f"{fid} {iid} 0 0 0 -9 {genotype_calls}\n" for fid, iid in all_pairs),
        encoding="utf-8",
    )
    base_prefix = tmp_path / "target_all"
    completed = subprocess.run(
        ["plink", "--file", str(tmp_path / "base"), "--make-bed", "--out", str(base_prefix)],
        capture_output=True, check=False, text=True,
    )
    assert completed.returncode == 0, completed.stderr
    dataset_prefixes = {}
    for cohort_id, pairs in (
        ("controls_unrelated", all_pairs[:5]),
        ("target_carriers_independent", all_pairs[5:]),
    ):
        keep = tmp_path / f"{cohort_id}.keep"
        keep.write_text("".join(f"{fid}\t{iid}\n" for fid, iid in pairs), encoding="utf-8")
        prefix = tmp_path / cohort_id
        completed = subprocess.run(
            ["plink", "--bfile", str(base_prefix), "--keep", str(keep), "--make-bed", "--out", str(prefix)],
            capture_output=True, check=False, text=True,
        )
        assert completed.returncode == 0, completed.stderr
        dataset_prefixes[cohort_id] = prefix
    parameters = {
        **_parameters({}),
        "minimum_roh_kb": 1.0, "minimum_roh_snps": 5,
        "maximum_density_kb_per_snp": 1000.0, "window_snps": 5,
        "minimum_genomewide_variants": 10, "minimum_genomewide_autosomes": 1,
        "minimum_target_chromosome_variants": 10,
    }

    publication = publish_roh(
        plink_executable="plink",
        genomewide_datasets={
            cohort_id: {"dataset_id": f"{cohort_id}_qc", "prefix": dataset_prefixes[cohort_id], "fam_pairs": pairs}
            for cohort_id, pairs in (
                ("controls_unrelated", all_pairs[:5]),
                ("target_carriers_independent", all_pairs[5:]),
            )
        },
        target_dataset={"dataset_id": "target_chromosome_all_qc", "prefix": base_prefix, "fam_pairs": all_pairs},
        common_variant_ids=marker_ids, common_autosome_count=1,
        target_variant_count=10, sample_by_pair={pair: pair[1] for pair in all_pairs},
        target_genotypes={pair[1]: "A/A" for pair in all_pairs},
        roles_by_sample={
            pair[1]: {"controls_unrelated" if index < 5 else "target_carriers_independent"}
            for index, pair in enumerate(all_pairs)
        },
        target_variant_id="v5", target_chromosome="1", target_bp=5000,
        output_dir=tmp_path / "roh", parameters=parameters,
    )

    assert publication.scope_statuses == {
        "GENOMEWIDE_BURDEN": "EVALUATED",
        "TARGET_CHROMOSOME": "EVALUATED",
    }
    assert validate_tsv_table(publication.burden_path, "roh_burden.schema.json").row_count == 10
    assert validate_tsv_table(publication.target_path, "target_in_roh.schema.json").row_count == 10
