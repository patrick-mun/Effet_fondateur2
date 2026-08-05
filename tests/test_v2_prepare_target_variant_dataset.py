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
GENOTYPE_HEADER = "SAMPLE_ID\tGENOTYPE\tGENOTYPE_SOURCE\tMEASUREMENT_METHOD\n"


def write_fake_plink(
    path: Path, *, fail_merge: bool = False, mendel_error: bool = False
) -> None:
    script = f'''#!{sys.executable}
import pathlib
import shutil
import sys

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
        alleles = sorted({{row[6 + marker_index * 2] for row in ped_rows}} | {{row[7 + marker_index * 2] for row in ped_rows}})
        allele_1 = alleles[0]
        allele_2 = alleles[1] if len(alleles) > 1 else "0"
        bim_lines.append(f"{{map_row[0]}} {{map_row[1]}} {{map_row[2]}} {{map_row[3]}} {{allele_1}} {{allele_2}}\\n")
    output_prefix.with_suffix(".bim").write_text("".join(bim_lines))
    output_prefix.with_suffix(".bed").write_bytes(b"\\x6c\\x1b\\x01")
elif "--bmerge" in sys.argv:
    if {fail_merge!r}:
        raise SystemExit(8)
    base_prefix = pathlib.Path(sys.argv[sys.argv.index("--bfile") + 1])
    merge_bed = pathlib.Path(sys.argv[sys.argv.index("--bmerge") + 1])
    merge_bim = pathlib.Path(sys.argv[sys.argv.index("--bmerge") + 2])
    shutil.copyfile(base_prefix.with_suffix(".fam"), output_prefix.with_suffix(".fam"))
    output_prefix.with_suffix(".bim").write_text(base_prefix.with_suffix(".bim").read_text() + merge_bim.read_text())
    output_prefix.with_suffix(".bed").write_bytes(b"\\x6c\\x1b\\x01")
elif "--mendel" in sys.argv:
    input_prefix = pathlib.Path(sys.argv[sys.argv.index("--bfile") + 1])
    target_id = input_prefix.with_suffix(".bim").read_text().splitlines()[-1].split()[1]
    content = "FID KID CHR SNP CODE ERROR\\n"
    if {mendel_error!r}:
        content += f"FAM1 I1 19 {{target_id}} 1 1\\n"
    output_prefix.with_suffix(".mendel").write_text(content)
else:
    raise SystemExit(9)
'''
    path.write_text(script, encoding="utf-8")
    path.chmod(0o755)


def write_acpa_source(path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    rows = [
        f"probe_{chromosome}\tAA\trs{chromosome}\t{chromosome}\t{chromosome * 100}"
        for chromosome in range(1, 23)
    ]
    path.write_text(ACPA_HEADER + "\n".join(rows) + "\n", encoding="utf-8")


def write_samples_metadata(path: Path) -> None:
    path.write_text(
        SAMPLES_HEADER
        + "sample_1\tsample_1.txt\tFAM1\tI1\t0\t0\tMALE\tAFFECTED\tFAMILY\t"
        "A/G\tlaboratory_report\tbatch_1\ttrue\ttrue\t\n"
        + "sample_2\tsample_2.txt\tFAM2\tI2\t0\t0\tFEMALE\tUNAFFECTED\tFAMILY\t"
        "\t\tbatch_1\ttrue\ttrue\t\n",
        encoding="utf-8",
    )


def write_variant_metadata(path: Path, *, position_bp: int = 100000) -> None:
    path.write_text(
        yaml.safe_dump(
            {
                "schema_version": "1.0.0",
                "assembly": "GRCh38",
                "gene": "GENE_EXAMPLE",
                "chromosome": 19,
                "position_bp": position_bp,
                "ref": "A",
                "alt": "G",
                "transcript": "TRANSCRIPT_EXAMPLE",
                "hgvs_c": "c.100A>G",
                "hgvs_p": "p.Lys34Arg",
                "project_variant_id": "target_GRCh38_19_100000_A_G",
                "rsid": None,
                "reference_validation": {
                    "status": "CONFIRMED",
                    "source_id": "synthetic_reference",
                    "confirmed_at": "2026-08-05T10:00:00Z",
                },
                "annotation_validation": {
                    "status": "CONFIRMED",
                    "source_id": "synthetic_transcript_annotation",
                    "confirmed_at": "2026-08-05T10:00:00Z",
                },
            },
            sort_keys=False,
        ),
        encoding="utf-8",
    )


def write_genotypes(
    path: Path,
    *,
    include_sample_2: bool = True,
    sample_2_source: str = "laboratory_report",
) -> None:
    rows = "sample_1\tA/G\tlaboratory_report\tvalidated_assay\n"
    if include_sample_2:
        rows += f"sample_2\tA/A\t{sample_2_source}\tvalidated_assay\n"
    path.write_text(GENOTYPE_HEADER + rows, encoding="utf-8")


def prepare_inputs(
    tmp_path: Path,
    *,
    include_sample_2: bool = True,
    sample_2_source: str = "laboratory_report",
    metadata_position_bp: int = 100000,
    fail_merge: bool = False,
    mendel_error: bool = False,
) -> tuple[Path, Path]:
    source_dir = tmp_path / "sources"
    write_acpa_source(source_dir / "sample_1.txt")
    write_acpa_source(source_dir / "sample_2.txt")
    samples_path = tmp_path / "samples.tsv"
    write_samples_metadata(samples_path)
    variant_path = tmp_path / "variant.yaml"
    write_variant_metadata(variant_path, position_bp=metadata_position_bp)
    genotypes_path = tmp_path / "genotypes.tsv"
    write_genotypes(
        genotypes_path,
        include_sample_2=include_sample_2,
        sample_2_source=sample_2_source,
    )
    plink_path = tmp_path / "fake_plink"
    write_fake_plink(
        plink_path, fail_merge=fail_merge, mendel_error=mendel_error
    )

    config = yaml.safe_load(
        (REPOSITORY_ROOT / "config" / "pipeline.example.yaml").read_text(
            encoding="utf-8"
        )
    )
    config["inputs"].update(
        {
            "acpa_samples_dir": str(source_dir),
            "sample_metadata": str(samples_path),
            "target_variant_metadata": str(variant_path),
            "target_variant_genotypes": str(genotypes_path),
        }
    )
    config["tools"]["plink"] = str(plink_path)
    config["target"].update(
        {
            "chromosome": 19,
            "position_bp": 100000,
            "project_variant_id": "target_GRCh38_19_100000_A_G",
        }
    )
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
    config["stages"]["prepare_target_variant_dataset"] = {
        "enabled": True,
        "parameters": {"plink_timeout_seconds": 30},
    }
    config_path = tmp_path / "config.yaml"
    config_path.write_text(
        yaml.safe_dump(config, allow_unicode=True, sort_keys=False),
        encoding="utf-8",
    )
    return config_path, tmp_path / "runs"


def test_target_variant_dataset_is_published_without_modifying_base(
    tmp_path: Path,
) -> None:
    config_path, runs_dir = prepare_inputs(tmp_path)

    run_dir = run_pipeline(config_path, runs_dir)

    base_dir = run_dir / "stages" / "03_convert_acpa"
    stage_dir = run_dir / "stages" / "04_prepare_target_variant_dataset"
    report = json.loads(
        (stage_dir / "target_variant_injection_report.json").read_text(
            encoding="utf-8"
        )
    )
    assert report["input_checksums"]["bed"] == json.loads(
        (base_dir / "target_chromosome_base.dataset.json").read_text(
            encoding="utf-8"
        )
    )["files"]["bed"]["sha256"]
    assert len((stage_dir / "target_variant.bim").read_text().splitlines()) == 2
    target_row = (stage_dir / "target_variant.bim").read_text().splitlines()[-1]
    assert set(target_row.split()[4:6]) == {"A", "G"}
    with (stage_dir / "target_genotype_audit.tsv").open(encoding="utf-8") as file:
        audit_rows = list(csv.DictReader(file, delimiter="\t"))
    assert [row["SAMPLE_ID"] for row in audit_rows] == ["sample_1", "sample_2"]
    assert report["mendel_status"] == "NOT_APPLICABLE"


def test_missing_individual_genotype_blocks_stage(tmp_path: Path) -> None:
    config_path, runs_dir = prepare_inputs(tmp_path, include_sample_2=False)

    with pytest.raises(StageExecutionError) as error:
        run_pipeline(config_path, runs_dir)

    assert error.value.return_code == 4


def test_clinical_genotype_source_blocks_stage(tmp_path: Path) -> None:
    config_path, runs_dir = prepare_inputs(
        tmp_path, sample_2_source="clinical_status"
    )

    with pytest.raises(StageExecutionError) as error:
        run_pipeline(config_path, runs_dir)

    assert error.value.return_code == 4


def test_metadata_config_mismatch_blocks_stage(tmp_path: Path) -> None:
    config_path, runs_dir = prepare_inputs(tmp_path, metadata_position_bp=100001)

    with pytest.raises(StageExecutionError) as error:
        run_pipeline(config_path, runs_dir)

    assert error.value.return_code == 4


def test_mendel_error_blocks_publication_and_keeps_failed_attempt(
    tmp_path: Path,
) -> None:
    config_path, runs_dir = prepare_inputs(tmp_path, mendel_error=True)

    with pytest.raises(StageExecutionError) as error:
        run_pipeline(config_path, runs_dir)

    assert error.value.return_code == 4
    run_dir = next(path for path in runs_dir.iterdir() if not path.name.startswith("."))
    assert not (run_dir / "stages" / "04_prepare_target_variant_dataset").exists()
    failed_attempts = list(
        (run_dir / "stages" / "attempts").glob(
            "04_prepare_target_variant_dataset.*.failed"
        )
    )
    assert len(failed_attempts) == 1
    assert (failed_attempts[0] / "target_variant_mendel.tsv").is_file()
    assert (failed_attempts[0] / "target_variant_injection_report.json").is_file()


def test_plink_merge_failure_uses_external_tool_return_code(tmp_path: Path) -> None:
    config_path, runs_dir = prepare_inputs(tmp_path, fail_merge=True)

    with pytest.raises(StageExecutionError) as error:
        run_pipeline(config_path, runs_dir)

    assert error.value.return_code == 3
