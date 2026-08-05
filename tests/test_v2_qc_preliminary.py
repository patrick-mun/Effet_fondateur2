import csv
import json
import sys
from pathlib import Path

import pytest
import yaml

from effet_fondateur.orchestrator import StageExecutionError
from effet_fondateur.orchestrator.pipeline import run_pipeline


REPOSITORY_ROOT = Path(__file__).resolve().parents[1]
SAMPLES_HEADER = (
    "SAMPLE_ID\tSOURCE_FILE\tFID\tIID\tPID\tMID\tSEX\tCLINICAL_STATUS\t"
    "GROUP_LABEL\tTARGET_GENOTYPE\tTARGET_GENOTYPE_SOURCE\tARRAY_BATCH\t"
    "INCLUDE_GENOMEWIDE\tINCLUDE_TARGET_CHROMOSOME\tNOTES_CODE\n"
)
ACPA_HEADER = (
    "# Array Type Name: CytoScan 750K Accel Array\n"
    "# UCSC Genomic Version: hg38\n"
    "Probe Set ID\tForward Strand Base Calls\tdbSNP RS ID\tChromosome\t"
    "Chromosomal Position\n"
)


def write_fake_plink(
    path: Path,
    command_audit_path: Path,
    *,
    exclude_sample: bool = False,
    exclude_all_samples: bool = False,
    exclude_variant: bool = False,
    duplicate_pair: bool = False,
    fail_qc: bool = False,
) -> None:
    script = f'''#!{sys.executable}
import pathlib
import shutil
import sys

audit_path = pathlib.Path({str(command_audit_path)!r})
with audit_path.open("a") as audit:
    audit.write(" ".join(sys.argv[1:]) + "\\n")
if "--version" in sys.argv:
    print("PLINK fake 1.0")
    raise SystemExit(0)
output_prefix = pathlib.Path(sys.argv[sys.argv.index("--out") + 1])
if "--file" in sys.argv:
    input_prefix = pathlib.Path(sys.argv[sys.argv.index("--file") + 1])
    ped_rows = [line.split() for line in input_prefix.with_suffix(".ped").read_text().splitlines() if line]
    map_rows = [line.split() for line in input_prefix.with_suffix(".map").read_text().splitlines() if line]
    output_prefix.with_suffix(".fam").write_text("".join(" ".join(row[:6]) + "\\n" for row in ped_rows))
    bim_lines = []
    for marker_index, map_row in enumerate(map_rows):
        alleles = sorted({{row[6 + marker_index * 2] for row in ped_rows if row[6 + marker_index * 2] != "0"}} | {{row[7 + marker_index * 2] for row in ped_rows if row[7 + marker_index * 2] != "0"}})
        allele_1 = alleles[0] if alleles else "0"
        allele_2 = alleles[1] if len(alleles) > 1 else "0"
        bim_lines.append(f"{{map_row[0]}} {{map_row[1]}} {{map_row[2]}} {{map_row[3]}} {{allele_1}} {{allele_2}}\\n")
    output_prefix.with_suffix(".bim").write_text("".join(bim_lines))
    output_prefix.with_suffix(".bed").write_bytes(b"\\x6c\\x1b\\x01")
elif "--genome" in sys.argv:
    content = "FID1 IID1 FID2 IID2 PI_HAT\\n"
    if {duplicate_pair!r}:
        content += "F1 I1 F3 I3 0.995\\n"
    output_prefix.with_suffix(".genome").write_text(content)
elif "--make-bed" in sys.argv and "--bfile" in sys.argv:
    base_prefix = pathlib.Path(sys.argv[sys.argv.index("--bfile") + 1])
    fam_rows = [line for line in base_prefix.with_suffix(".fam").read_text().splitlines() if line]
    bim_rows = [line for line in base_prefix.with_suffix(".bim").read_text().splitlines() if line]
    if "--remove" in sys.argv:
        remove_path = pathlib.Path(sys.argv[sys.argv.index("--remove") + 1])
        removed = {{tuple(line.split()[:2]) for line in remove_path.read_text().splitlines() if line}}
        fam_rows = [line for line in fam_rows if tuple(line.split()[:2]) not in removed]
    if "--exclude" in sys.argv:
        exclude_path = pathlib.Path(sys.argv[sys.argv.index("--exclude") + 1])
        excluded = {{line.strip() for line in exclude_path.read_text().splitlines() if line}}
        bim_rows = [line for line in bim_rows if line.split()[1] not in excluded]
    output_prefix.with_suffix(".fam").write_text("".join(line + "\\n" for line in fam_rows))
    output_prefix.with_suffix(".bim").write_text("".join(line + "\\n" for line in bim_rows))
    output_prefix.with_suffix(".bed").write_bytes(b"\\x6c\\x1b\\x01")
elif "--missing" in sys.argv and "--within" in sys.argv:
    if {fail_qc!r}:
        raise SystemExit(8)
    base_prefix = pathlib.Path(sys.argv[sys.argv.index("--bfile") + 1])
    cluster_path = pathlib.Path(sys.argv[sys.argv.index("--within") + 1])
    variants = [line.split() for line in base_prefix.with_suffix(".bim").read_text().splitlines() if line]
    clusters = sorted({{line.split()[2] for line in cluster_path.read_text().splitlines() if line}})
    lines = ["CHR SNP CLST N_MISS N_CLST N_GENO F_MISS\\n"]
    for variant in variants:
        for cluster in clusters:
            missingness = 0.5 if variant[1] == "probe_3" and cluster == "batch_1" else 0
            sample_count = 2 if cluster == "batch_1" else 1
            lines.append(f"{{variant[0]}} {{variant[1]}} {{cluster}} {{int(missingness * sample_count)}} {{sample_count}} {{sample_count}} {{missingness}}\\n")
    output_prefix.with_suffix(".lmiss").write_text("".join(lines))
    output_prefix.with_suffix(".imiss").write_text("FID IID MISS_PHENO N_MISS N_GENO F_MISS\\n")
elif "--missing" in sys.argv and "--freq" in sys.argv and "--het" in sys.argv:
    if {fail_qc!r}:
        raise SystemExit(8)
    base_prefix = pathlib.Path(sys.argv[sys.argv.index("--bfile") + 1])
    fam_rows = [line.split() for line in base_prefix.with_suffix(".fam").read_text().splitlines() if line]
    variants = [line.split() for line in base_prefix.with_suffix(".bim").read_text().splitlines() if line]
    imiss = ["FID IID MISS_PHENO N_MISS N_GENO F_MISS\\n"]
    het = ["FID IID O(HOM) E(HOM) N(NM) F\\n"]
    for index, row in enumerate(fam_rows):
        if {exclude_all_samples!r} or ({exclude_sample!r} and row[1] == "I2"):
            n_miss = len(variants)
            f_miss = 1
        else:
            n_miss = 0
            f_miss = 0
        imiss.append(f"{{row[0]}} {{row[1]}} N {{n_miss}} {{len(variants)}} {{f_miss}}\\n")
        het.append(f"{{row[0]}} {{row[1]}} {{10 + index}} 10.5 20 {{[-0.1, 0, 0.1][index]}}\\n")
    lmiss = ["CHR SNP N_MISS N_GENO F_MISS\\n"]
    frq = ["CHR SNP A1 A2 MAF NCHROBS\\n"]
    for variant in variants:
        missing = 1 if {exclude_variant!r} and variant[1] == "probe_2" else 0
        f_miss = missing / len(fam_rows)
        lmiss.append(f"{{variant[0]}} {{variant[1]}} {{missing}} {{len(fam_rows)}} {{f_miss}}\\n")
        maf = 0.01 if variant[1] == "probe_1" else 0.2
        frq.append(f"{{variant[0]}} {{variant[1]}} A G {{maf}} {{2 * (len(fam_rows) - missing)}}\\n")
    output_prefix.with_suffix(".imiss").write_text("".join(imiss))
    output_prefix.with_suffix(".lmiss").write_text("".join(lmiss))
    output_prefix.with_suffix(".frq").write_text("".join(frq))
    output_prefix.with_suffix(".het").write_text("".join(het))
else:
    raise SystemExit(9)
'''
    path.write_text(script, encoding="utf-8")
    path.chmod(0o755)


def write_acpa_source(path: Path, sample_index: int) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    rows = []
    for chromosome in range(1, 23):
        if sample_index == 1:
            genotype = "AA"
        elif sample_index == 2:
            genotype = "AG" if chromosome % 2 else "AA"
        else:
            genotype = "GG" if chromosome % 2 else "AA"
        rows.append(
            f"probe_{chromosome}\t{genotype}\trs{chromosome}\t{chromosome}\t"
            f"{chromosome * 100}"
        )
    path.write_text(ACPA_HEADER + "\n".join(rows) + "\n", encoding="utf-8")


def write_samples_metadata(path: Path) -> None:
    path.write_text(
        SAMPLES_HEADER
        + "sample_1\tsample_1.txt\tF1\tI1\t0\t0\tMALE\tAFFECTED\tFAMILY\t"
        "\t\tbatch_1\ttrue\ttrue\t\n"
        + "sample_2\tsample_2.txt\tF2\tI2\t0\t0\tFEMALE\tUNAFFECTED\tFAMILY\t"
        "\t\tbatch_1\ttrue\ttrue\t\n"
        + "sample_3\tsample_3.txt\tF3\tI3\t0\t0\tMALE\tUNAFFECTED\tFAMILY\t"
        "\t\tbatch_2\ttrue\ttrue\t\n",
        encoding="utf-8",
    )


def prepare_inputs(
    tmp_path: Path,
    *,
    exclude_sample: bool = False,
    exclude_all_samples: bool = False,
    exclude_variant: bool = False,
    duplicate_pair: bool = False,
    fail_qc: bool = False,
    duplicate_scan_minimum: int = 100,
) -> tuple[Path, Path, Path]:
    source_dir = tmp_path / "sources"
    for sample_index in range(1, 4):
        write_acpa_source(source_dir / f"sample_{sample_index}.txt", sample_index)
    samples_path = tmp_path / "samples.tsv"
    write_samples_metadata(samples_path)
    command_audit_path = tmp_path / "plink_commands.txt"
    plink_path = tmp_path / "fake_plink"
    write_fake_plink(
        plink_path,
        command_audit_path,
        exclude_sample=exclude_sample,
        exclude_all_samples=exclude_all_samples,
        exclude_variant=exclude_variant,
        duplicate_pair=duplicate_pair,
        fail_qc=fail_qc,
    )
    config = yaml.safe_load(
        (REPOSITORY_ROOT / "config" / "pipeline.example.yaml").read_text(
            encoding="utf-8"
        )
    )
    config["inputs"]["acpa_samples_dir"] = str(source_dir)
    config["inputs"]["sample_metadata"] = str(samples_path)
    config["tools"]["plink"] = str(plink_path)
    config["target"]["chromosome"] = 19
    config["stages"]["validate_sources"]["enabled"] = True
    config["stages"]["build_sample_registry"] = {
        "enabled": True,
        "parameters": {
            "manual_approval": {
                "approved": True,
                "decision_id": "synthetic_registry_review",
                "reviewer_role": "data_steward",
                "approved_at": "2026-08-05T10:00:00Z",
            }
        },
    }
    config["stages"]["convert_acpa"] = {
        "enabled": True,
        "parameters": {"marker_mode": "intersection"},
    }
    config["stages"]["qc_preliminary"] = {
        "enabled": True,
        "parameters": {
            "sample_missingness_max": 0.1,
            "variant_missingness_max": 0.05,
            "kinship_maf_alert_threshold": 0.05,
            "heterozygosity_z_alert": 3.0,
            "batch_missingness_delta_alert": 0.02,
            "duplicate_pi_hat_alert": 0.98,
            "duplicate_scan_min_polymorphic_variants": duplicate_scan_minimum,
            "plink_timeout_seconds": 30,
        },
    }
    config_path = tmp_path / "config.yaml"
    config_path.write_text(
        yaml.safe_dump(config, allow_unicode=True, sort_keys=False),
        encoding="utf-8",
    )
    return config_path, tmp_path / "runs", command_audit_path


def test_qc_preliminary_filters_only_missingness_failures(tmp_path: Path) -> None:
    config_path, runs_dir, command_audit_path = prepare_inputs(
        tmp_path, exclude_sample=True, exclude_variant=True
    )

    run_dir = run_pipeline(config_path, runs_dir)

    stage_dir = run_dir / "stages" / "05_qc_preliminary"
    assert len((stage_dir / "genomewide_pre_qc.fam").read_text().splitlines()) == 2
    assert len((stage_dir / "genomewide_pre_qc.bim").read_text().splitlines()) == 21
    report = json.loads(
        (stage_dir / "qc_preliminary_report.json").read_text(encoding="utf-8")
    )
    assert report["excluded_sample_count"] == 1
    assert report["excluded_variant_count"] == 1
    assert report["hwe_filter_applied"] is False
    assert "--hwe" not in command_audit_path.read_text(encoding="utf-8")


def test_qc_preliminary_records_maf_batch_and_non_evaluable_checks(
    tmp_path: Path,
) -> None:
    config_path, runs_dir, _ = prepare_inputs(tmp_path)

    run_dir = run_pipeline(config_path, runs_dir)

    stage_dir = run_dir / "stages" / "05_qc_preliminary"
    with (stage_dir / "qc_variant_metrics.tsv").open(encoding="utf-8") as file:
        variants = list(csv.DictReader(file, delimiter="\t"))
    assert "low_maf_for_kinship" in variants[0]["ALERT_CODES"]
    probe_3 = next(row for row in variants if row["VARIANT_ID"] == "probe_3")
    assert "batch_missingness_delta" in probe_3["ALERT_CODES"]
    alerts = (stage_dir / "qc_alerts.tsv").read_text(encoding="utf-8")
    assert "autosomal_dataset_without_sex_chromosomes" in alerts
    assert "insufficient_polymorphic_variants" in alerts
    assert "cohorts_not_frozen" in alerts


def test_duplicate_scan_publishes_only_an_alert(tmp_path: Path) -> None:
    config_path, runs_dir, _ = prepare_inputs(
        tmp_path, duplicate_pair=True, duplicate_scan_minimum=1
    )

    run_dir = run_pipeline(config_path, runs_dir)

    stage_dir = run_dir / "stages" / "05_qc_preliminary"
    report = json.loads(
        (stage_dir / "qc_preliminary_report.json").read_text(encoding="utf-8")
    )
    assert report["duplicate_scan_status"] == "EVALUATED"
    assert report["potential_duplicate_pair_count"] == 1
    assert len((stage_dir / "genomewide_pre_qc.fam").read_text().splitlines()) == 3
    assert "potential_technical_duplicate" in (
        stage_dir / "qc_alerts.tsv"
    ).read_text(encoding="utf-8")


def test_all_samples_excluded_is_a_scientific_block(tmp_path: Path) -> None:
    config_path, runs_dir, _ = prepare_inputs(tmp_path, exclude_all_samples=True)

    with pytest.raises(StageExecutionError) as error:
        run_pipeline(config_path, runs_dir)

    assert error.value.return_code == 4


def test_invalid_qc_threshold_uses_input_error_code(tmp_path: Path) -> None:
    config_path, runs_dir, _ = prepare_inputs(tmp_path)
    config = yaml.safe_load(config_path.read_text(encoding="utf-8"))
    config["stages"]["qc_preliminary"]["parameters"][
        "sample_missingness_max"
    ] = 1.5
    config_path.write_text(yaml.safe_dump(config, sort_keys=False), encoding="utf-8")

    with pytest.raises(StageExecutionError) as error:
        run_pipeline(config_path, runs_dir)

    assert error.value.return_code == 2


def test_plink_qc_failure_uses_external_tool_return_code(tmp_path: Path) -> None:
    config_path, runs_dir, _ = prepare_inputs(tmp_path, fail_qc=True)

    with pytest.raises(StageExecutionError) as error:
        run_pipeline(config_path, runs_dir)

    assert error.value.return_code == 3
