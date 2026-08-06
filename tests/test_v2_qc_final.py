import csv
import json
import sys
from pathlib import Path

import pytest
import yaml

from effet_fondateur.orchestrator import StageExecutionError
from effet_fondateur.orchestrator.pipeline import run_pipeline
from test_v2_freeze_cohorts import prepare_cohort_inputs


def write_final_qc_plink(
    path: Path,
    delegated_plink: Path,
    *,
    exclude_sample_iid: str | None = None,
    exclude_all_samples: bool = False,
    fail_metrics: bool = False,
) -> None:
    script = f'''#!{sys.executable}
import os
import pathlib
import sys

dataset_names = (
    "controls_unrelated_qc",
    "target_carriers_independent_qc",
    "target_chromosome_all_qc",
)
output_prefix = pathlib.Path(sys.argv[sys.argv.index("--out") + 1]) if "--out" in sys.argv else None
is_final_qc = output_prefix is not None and any(name in output_prefix.name for name in dataset_names)
if not is_final_qc:
    os.execv({str(delegated_plink)!r}, [{str(delegated_plink)!r}, *sys.argv[1:]])
if {fail_metrics!r} and "--missing" in sys.argv:
    raise SystemExit(8)
base_prefix = pathlib.Path(sys.argv[sys.argv.index("--bfile") + 1])
fam_rows = [line.split() for line in base_prefix.with_suffix(".fam").read_text().splitlines() if line]
bim_rows = [line.split() for line in base_prefix.with_suffix(".bim").read_text().splitlines() if line]
if "--make-bed" in sys.argv:
    if "--keep" in sys.argv:
        keep_path = pathlib.Path(sys.argv[sys.argv.index("--keep") + 1])
        kept = {{tuple(line.split()[:2]) for line in keep_path.read_text().splitlines() if line}}
        fam_rows = [row for row in fam_rows if tuple(row[:2]) in kept]
    if "--remove" in sys.argv:
        remove_path = pathlib.Path(sys.argv[sys.argv.index("--remove") + 1])
        removed = {{tuple(line.split()[:2]) for line in remove_path.read_text().splitlines() if line}}
        fam_rows = [row for row in fam_rows if tuple(row[:2]) not in removed]
    if "--exclude" in sys.argv:
        exclude_path = pathlib.Path(sys.argv[sys.argv.index("--exclude") + 1])
        excluded = {{line.strip() for line in exclude_path.read_text().splitlines() if line}}
        bim_rows = [row for row in bim_rows if row[1] not in excluded]
    output_prefix.with_suffix(".fam").write_text("".join(" ".join(row) + "\\n" for row in fam_rows))
    output_prefix.with_suffix(".bim").write_text("".join(" ".join(row) + "\\n" for row in bim_rows))
    output_prefix.with_suffix(".bed").write_bytes(b"\\x6c\\x1b\\x01")
    raise SystemExit(0)
if "--missing" in sys.argv:
    imiss = ["FID IID MISS_PHENO N_MISS N_GENO F_MISS\\n"]
    for fam_row in fam_rows:
        f_miss = 1.0 if {exclude_all_samples!r} or fam_row[1] == {exclude_sample_iid!r} else 0.0
        n_miss = len(bim_rows) if f_miss else 0
        imiss.append(f"{{fam_row[0]}} {{fam_row[1]}} N {{n_miss}} {{len(bim_rows)}} {{f_miss}}\\n")
    output_prefix.with_suffix(".imiss").write_text("".join(imiss))
    lmiss = ["CHR SNP N_MISS N_GENO F_MISS\\n"]
    frequencies = ["CHR SNP A1 A2 MAF NCHROBS\\n"]
    hwe = ["CHR SNP TEST A1 A2 GENO O_HET E_HET P\\n"]
    for variant in bim_rows:
        is_target = variant[1].startswith("target_")
        high_missingness = variant[1] == "probe_2" or variant[1] == "probe_19" or is_target
        f_miss = 0.5 if high_missingness else 0.0
        n_miss = max(1, int(len(fam_rows) * f_miss)) if high_missingness else 0
        maf = 0.001 if is_target else 0.01 if variant[1] == "probe_3" else 0.2
        hwe_p = 0.00000001 if variant[1] == "probe_4" else 0.5
        lmiss.append(f"{{variant[0]}} {{variant[1]}} {{n_miss}} {{len(fam_rows)}} {{f_miss}}\\n")
        frequencies.append(f"{{variant[0]}} {{variant[1]}} A G {{maf}} {{2 * len(fam_rows)}}\\n")
        hwe.append(f"{{variant[0]}} {{variant[1]}} ALL A G 0/0/0 0 0 {{hwe_p}}\\n")
    output_prefix.with_suffix(".lmiss").write_text("".join(lmiss))
    if "--freq" in sys.argv:
        output_prefix.with_suffix(".frq").write_text("".join(frequencies))
    if "--hardy" in sys.argv:
        output_prefix.with_suffix(".hwe").write_text("".join(hwe))
    raise SystemExit(0)
raise SystemExit(9)
'''
    path.write_text(script, encoding="utf-8")
    path.chmod(0o755)


def prepare_final_qc_inputs(
    tmp_path: Path,
    *,
    exclude_sample_iid: str | None = None,
    exclude_all_samples: bool = False,
    fail_metrics: bool = False,
    parameter_overrides: dict[str, object] | None = None,
) -> tuple[Path, Path]:
    config_path, runs_dir = prepare_cohort_inputs(tmp_path)
    config = yaml.safe_load(config_path.read_text(encoding="utf-8"))
    delegated_plink = Path(config["tools"]["plink"])
    final_qc_plink = tmp_path / "final_qc_plink"
    write_final_qc_plink(
        final_qc_plink,
        delegated_plink,
        exclude_sample_iid=exclude_sample_iid,
        exclude_all_samples=exclude_all_samples,
        fail_metrics=fail_metrics,
    )
    config["tools"]["plink"] = str(final_qc_plink)
    parameters: dict[str, object] = {
        "sample_missingness_max": 0.1,
        "controls_variant_missingness_max": 0.05,
        "controls_maf_min": 0.05,
        "controls_hwe_p_min": 0.000001,
        "carriers_variant_missingness_max": 0.05,
        "target_region_variant_missingness_max": 0.05,
        "target_region_maf_min": 0.01,
        "batch_missingness_delta_alert": 0.02,
        "plink_timeout_seconds": 30,
    }
    parameters.update(parameter_overrides or {})
    config["stages"]["qc_final"] = {
        "enabled": True,
        "parameters": parameters,
    }
    config_path.write_text(
        yaml.safe_dump(config, allow_unicode=True, sort_keys=False),
        encoding="utf-8",
    )
    return config_path, runs_dir


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8") as input_file:
        return list(csv.DictReader(input_file, delimiter="\t"))


def test_final_qc_applies_distinct_cohort_policies(tmp_path: Path) -> None:
    config_path, runs_dir = prepare_final_qc_inputs(tmp_path)

    run_dir = run_pipeline(config_path, runs_dir)

    stage_dir = run_dir / "stages" / "10_qc_final"
    variants = read_tsv(stage_dir / "variant_qc_final.tsv")
    controls = {
        row["VARIANT_ID"]: row
        for row in variants
        if row["DATASET_ID"] == "controls_unrelated_qc"
    }
    carriers = {
        row["VARIANT_ID"]: row
        for row in variants
        if row["DATASET_ID"] == "target_carriers_independent_qc"
    }
    target_rows = [
        row
        for row in variants
        if row["DATASET_ID"] == "target_chromosome_all_qc"
    ]
    assert controls["probe_2"]["FILTER_CODES"] == "variant_missingness"
    assert controls["probe_3"]["FILTER_CODES"] == "low_maf_controls"
    assert controls["probe_4"]["FILTER_CODES"] == "hwe_controls"
    assert all(controls[variant_id]["QC_STATUS"] == "EXCLUDED" for variant_id in ("probe_2", "probe_3", "probe_4"))
    assert carriers["probe_2"]["QC_STATUS"] == "EXCLUDED"
    assert carriers["probe_3"]["QC_STATUS"] == "PASS"
    assert carriers["probe_4"]["HWE_STATUS"] == "NOT_EVALUATED_POLICY"
    target = next(row for row in target_rows if row["IS_TARGET_VARIANT"] == "true")
    assert target["QC_STATUS"] == "RETAINED_EXCEPTION"
    assert target["FILTER_EXCEPTION"] == "true"
    assert target["PRESENT_AFTER_QC"] == "true"
    target_bim = (stage_dir / "target_chromosome_all_qc.bim").read_text(encoding="utf-8")
    assert target["VARIANT_ID"] in target_bim
    assert "probe_19" not in target_bim


def test_final_qc_publishes_versioned_datasets_and_audit(tmp_path: Path) -> None:
    config_path, runs_dir = prepare_final_qc_inputs(tmp_path)

    run_dir = run_pipeline(config_path, runs_dir)

    stage_dir = run_dir / "stages" / "10_qc_final"
    report = json.loads((stage_dir / "qc_final_report.json").read_text(encoding="utf-8"))
    assert set(report["dataset_counts"]) == {
        "controls_unrelated_qc",
        "target_carriers_independent_qc",
        "target_chromosome_all_qc",
    }
    assert report["target_variant_exception_count"] == 1
    for dataset_id, counts in report["dataset_counts"].items():
        descriptor = json.loads(
            (stage_dir / f"{dataset_id}.dataset.json").read_text(encoding="utf-8")
        )
        assert descriptor["source_format"] == "PLINK_FINAL_QC"
        assert descriptor["sample_count"] == counts["samples"]
        assert descriptor["variant_count"] == counts["variants"]
        assert descriptor["sample_set_id"].startswith(f"{dataset_id}_samples_")
        assert descriptor["variant_set_id"].startswith(f"{dataset_id}_variants_")
    audit = json.loads((stage_dir / "audit.json").read_text(encoding="utf-8"))
    assert audit["warnings"] == [
        {"code": "target_variant_forced_retention", "count": 1}
    ]
    assert (stage_dir / "batch_qc_final.tsv").is_file()


def test_final_qc_excludes_high_missingness_sample_and_reports_batch_alert(
    tmp_path: Path,
) -> None:
    config_path, runs_dir = prepare_final_qc_inputs(
        tmp_path,
        exclude_sample_iid="I3",
    )

    run_dir = run_pipeline(config_path, runs_dir)

    stage_dir = run_dir / "stages" / "10_qc_final"
    samples = read_tsv(stage_dir / "sample_qc_final.tsv")
    sample_3_rows = [row for row in samples if row["IID"] == "I3"]
    assert len(sample_3_rows) == 2
    assert all(row["QC_STATUS"] == "EXCLUDED" for row in sample_3_rows)
    assert all(row["EXCLUSION_CODE"] == "sample_missingness" for row in sample_3_rows)
    target_fam = (stage_dir / "target_chromosome_all_qc.fam").read_text(
        encoding="utf-8"
    )
    assert " I3 " not in f" {target_fam} "
    batches = read_tsv(stage_dir / "batch_qc_final.tsv")
    batch_2_rows = [row for row in batches if row["ARRAY_BATCH"] == "batch_2"]
    assert len(batch_2_rows) == 2
    assert all(row["QC_STATUS"] == "ALERT" for row in batch_2_rows)
    assert all(row["ALERT_CODE"] == "batch_missingness_delta" for row in batch_2_rows)


def test_final_qc_blocks_if_every_sample_fails(tmp_path: Path) -> None:
    config_path, runs_dir = prepare_final_qc_inputs(
        tmp_path,
        exclude_all_samples=True,
    )

    with pytest.raises(StageExecutionError) as error:
        run_pipeline(config_path, runs_dir)

    assert error.value.return_code == 4


def test_final_qc_plink_failure_uses_external_tool_code(tmp_path: Path) -> None:
    config_path, runs_dir = prepare_final_qc_inputs(tmp_path, fail_metrics=True)

    with pytest.raises(StageExecutionError) as error:
        run_pipeline(config_path, runs_dir)

    assert error.value.return_code == 3


def test_final_qc_invalid_parameter_uses_input_code(tmp_path: Path) -> None:
    config_path, runs_dir = prepare_final_qc_inputs(
        tmp_path,
        parameter_overrides={"controls_maf_min": 1.1},
    )

    with pytest.raises(StageExecutionError) as error:
        run_pipeline(config_path, runs_dir)

    assert error.value.return_code == 2
