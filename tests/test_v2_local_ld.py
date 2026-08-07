from __future__ import annotations

from pathlib import Path
import shutil
import subprocess

import pytest

from effet_fondateur.contracts import validate_tsv_table
from effet_fondateur.ld import publish_local_ld
from effet_fondateur.orchestrator.pipeline import ANALYZE_LOCAL_LD_STAGE
from effet_fondateur.stages.analyze_local_ld import (
    AnalyzeLocalLdInputError,
    _parameters,
)


def _write_keep(path: Path, count: int) -> None:
    path.write_text(
        "".join(f"family_{index}\tsample_{index}\n" for index in range(1, count + 1)),
        encoding="utf-8",
    )


def _fake_plink(arguments: list[str]) -> None:
    output_prefix = Path(arguments[arguments.index("--out") + 1])
    if "--freq" in arguments:
        output_prefix.with_suffix(".frq").write_text(
            "CHR SNP A1 A2 MAF NCHROBS\n"
            "1 v1 G A 0.2 10\n"
            "1 target T C 0.3 10\n"
            "1 v3 C G 0.4 10\n",
            encoding="utf-8",
        )
        output_prefix.with_suffix(".traw").write_text(
            "CHR SNP (C)M POS COUNTED ALT F1_I1 F2_I2 F3_I3 F4_I4 F5_I5\n"
            "1 v1 0 100 G A 0 1 2 0 1\n"
            "1 target 0.1 200 T C 0 1 1 2 0\n"
            "1 v3 0.3 400 C G 2 1 0 1 2\n",
            encoding="utf-8",
        )
        return
    if "dprime" in arguments:
        output_prefix.with_suffix(".ld").write_text(
            "CHR_A BP_A SNP_A MAF_A CHR_B BP_B SNP_B MAF_B R2 DP\n"
            "1 100 v1 0.2 1 200 target 0.3 0.5 0.8\n"
            "1 100 v1 0.2 1 400 v3 0.4 0.6 0.7\n"
            "1 200 target 0.3 1 400 v3 0.4 0.7 0.9\n",
            encoding="utf-8",
        )
    else:
        output_prefix.with_suffix(".ld").write_text(
            "CHR_A BP_A SNP_A MAF_A CHR_B BP_B SNP_B MAF_B R2\n"
            "1 100 v1 0.2 1 200 target 0.3 0.4\n"
            "1 100 v1 0.2 1 400 v3 0.4 0.45\n"
            "1 200 target 0.3 1 400 v3 0.4 0.55\n",
            encoding="utf-8",
        )


def test_local_ld_keeps_genotype_r2_separate_from_haplotype_dprime(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch,
) -> None:
    def fake_run(_executable: str, arguments: list[str], _timeout: int) -> None:
        _fake_plink(arguments)

    monkeypatch.setattr("effet_fondateur.ld.local._run_plink", fake_run)
    keep_paths = {}
    for cohort_id, count in (
        ("controls_unrelated", 5),
        ("target_carriers_independent", 2),
        ("family_noncarriers", 5),
    ):
        keep_path = tmp_path / f"{cohort_id}.keep"
        _write_keep(keep_path, count)
        keep_paths[cohort_id] = keep_path
    variants = [
        {"VARIANT_ORDER": 1, "VARIANT_ID": "v1", "CHROMOSOME": "1", "POSITION_BP": 100, "POSITION_CM": 0.0, "A1": "G", "A2": "A", "IS_TARGET_VARIANT": False},
        {"VARIANT_ORDER": 2, "VARIANT_ID": "target", "CHROMOSOME": "1", "POSITION_BP": 200, "POSITION_CM": 0.1, "A1": "T", "A2": "C", "IS_TARGET_VARIANT": True},
        {"VARIANT_ORDER": 3, "VARIANT_ID": "v3", "CHROMOSOME": "1", "POSITION_BP": 400, "POSITION_CM": 0.3, "A1": "C", "A2": "G", "IS_TARGET_VARIANT": False},
    ]

    publication = publish_local_ld(
        plink_executable="plink", dataset_prefix=tmp_path / "region",
        variants=variants, cohort_keep_paths=keep_paths,
        cohort_sample_counts={"controls_unrelated": 5, "target_carriers_independent": 2, "family_noncarriers": 5},
        cohort_independent={"controls_unrelated": True, "target_carriers_independent": True, "family_noncarriers": False},
        output_dir=tmp_path / "output", minimum_samples=5, primary_samples=20,
        minimum_called_samples_per_pair=5, minimum_variant_maf=0.05,
        minimum_pairs_per_bin=1,
        max_pair_distance_bp=1_000, max_pair_distance_cm=0.3,
        distance_bins_cm=[0.1, 0.3], max_pair_count=10,
        timeout_seconds=30,
    )

    pairs = validate_tsv_table(publication.pairs_path, "local_ld_pairs.schema.json").rows
    target_pair = next(row for row in pairs if row["VARIANT_A"] == "v1" and row["VARIANT_B"] == "target")
    assert float(target_pair["R2_GENOTYPE"]) == pytest.approx(0.4)
    assert float(target_pair["R2_HAPLOTYPE_ML"]) == pytest.approx(0.5)
    assert float(target_pair["D_PRIME_ABS"]) == pytest.approx(0.8)
    assert target_pair["INVOLVES_TARGET"] is True
    assert publication.cohort_statuses == {
        "controls_unrelated": "EXPLORATORY_SMALL_N",
        "target_carriers_independent": "NOT_EVALUATED",
        "family_noncarriers": "NOT_EVALUATED",
    }
    variants_table = validate_tsv_table(publication.variants_path, "local_ld_variants.schema.json").rows
    assert {row["VARIANT_STATUS"] for row in variants_table if row["COHORT_ID"] == "target_carriers_independent"} == {"COHORT_TOO_SMALL"}
    assert {row["VARIANT_STATUS"] for row in variants_table if row["COHORT_ID"] == "family_noncarriers"} == {"COHORT_NOT_INDEPENDENT"}
    summaries = validate_tsv_table(publication.summary_path, "local_ld_summary.schema.json").rows
    assert any(row["SUMMARY_STATUS"] == "EVALUATED" for row in summaries if row["COHORT_ID"] == "controls_unrelated")
    assert len(publication.native_paths) == 2


def test_local_ld_contract_is_secondary_and_independent_from_stages_13_14() -> None:
    assert ANALYZE_LOCAL_LD_STAGE.critical is False
    assert ANALYZE_LOCAL_LD_STAGE.dependencies == (
        "freeze_cohorts", "qc_final", "prepare_target_region",
    )
    assert "infer_founder_haplotype" not in ANALYZE_LOCAL_LD_STAGE.dependencies
    assert "estimate_variant_age" not in ANALYZE_LOCAL_LD_STAGE.dependencies


def test_local_ld_parameters_require_bins_to_cover_configured_window() -> None:
    with pytest.raises(AnalyzeLocalLdInputError, match="distance_bins_cm"):
        _parameters({"max_pair_distance_cm": 1.0, "distance_bins_cm": [0.1, 0.5]})


@pytest.mark.skipif(shutil.which("plink") is None, reason="PLINK 1.9 absent")
def test_local_ld_smoke_with_real_plink_on_temporary_synthetic_data(tmp_path: Path) -> None:
    map_path = tmp_path / "region.map"
    ped_path = tmp_path / "region.ped"
    map_path.write_text(
        "1 v1 0 100\n1 target 0.1 200\n1 v3 0.3 400\n",
        encoding="utf-8",
    )
    genotypes = [
        "A A C C G G", "A G C T G C", "G G T T C C",
        "A G C T G C", "A A C C C C",
    ]
    ped_path.write_text(
        "".join(
            f"family_{index} sample_{index} 0 0 0 -9 {genotype}\n"
            for index, genotype in enumerate(genotypes, start=1)
        ),
        encoding="utf-8",
    )
    dataset_prefix = tmp_path / "region_binary"
    completed = subprocess.run(
        ["plink", "--file", str(tmp_path / "region"), "--make-bed", "--out", str(dataset_prefix)],
        capture_output=True, check=False, text=True,
    )
    assert completed.returncode == 0, completed.stderr
    keep_paths = {}
    for cohort_id, count in (
        ("controls_unrelated", 5),
        ("target_carriers_independent", 0),
        ("family_noncarriers", 0),
    ):
        path = tmp_path / f"{cohort_id}.keep"
        _write_keep(path, count)
        keep_paths[cohort_id] = path
    variants = [
        {"VARIANT_ORDER": 1, "VARIANT_ID": "v1", "CHROMOSOME": "1", "POSITION_BP": 100, "POSITION_CM": 0.0, "A1": "G", "A2": "A", "IS_TARGET_VARIANT": False},
        {"VARIANT_ORDER": 2, "VARIANT_ID": "target", "CHROMOSOME": "1", "POSITION_BP": 200, "POSITION_CM": 0.1, "A1": "T", "A2": "C", "IS_TARGET_VARIANT": True},
        {"VARIANT_ORDER": 3, "VARIANT_ID": "v3", "CHROMOSOME": "1", "POSITION_BP": 400, "POSITION_CM": 0.3, "A1": "C", "A2": "G", "IS_TARGET_VARIANT": False},
    ]

    publication = publish_local_ld(
        plink_executable="plink", dataset_prefix=dataset_prefix,
        variants=variants, cohort_keep_paths=keep_paths,
        cohort_sample_counts={"controls_unrelated": 5, "target_carriers_independent": 0, "family_noncarriers": 0},
        cohort_independent={cohort_id: True for cohort_id in keep_paths},
        output_dir=tmp_path / "ld", minimum_samples=5, primary_samples=20,
        minimum_called_samples_per_pair=5, minimum_variant_maf=0.05,
        minimum_pairs_per_bin=1,
        max_pair_distance_bp=1_000, max_pair_distance_cm=0.3,
        distance_bins_cm=[0.1, 0.3], max_pair_count=10,
        timeout_seconds=30,
    )

    assert publication.evaluated_pair_count == 3
    assert validate_tsv_table(publication.pairs_path, "local_ld_pairs.schema.json").row_count == 3
