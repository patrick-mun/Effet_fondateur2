import csv
import hashlib
import json
import sys
from pathlib import Path

import pytest
import yaml

from effet_fondateur.orchestrator import IntegrityError
from effet_fondateur.orchestrator.pipeline import resume_pipeline, run_pipeline
from effet_fondateur.references import cache_reference_panel, resolve_reference_panel
from effet_fondateur.stages import phase_target_region
from effet_fondateur.references import (
    ReferenceHarmonizationIntegrityError,
    ReferenceWindowError,
)
from test_v2_prepare_target_region import prepare_target_region_inputs


def _md5(content: bytes) -> str:
    return hashlib.md5(content, usedforsecurity=False).hexdigest()


def _sha256(content: bytes) -> str:
    return hashlib.sha256(content).hexdigest()


def _write_synthetic_catalog(path: Path, contents: dict[str, bytes]) -> None:
    base_url = "https://ftp.1000genomes.ebi.ac.uk/synthetic_stage_12"
    vcf_content = contents["synthetic.chr19.panel.vcf.gz"]
    index_content = contents["synthetic.chr19.panel.vcf.gz.tbi"]
    catalog = {
        "schema_version": "1.0.0",
        "catalog_id": "synthetic_stage_12_catalog_v1",
        "panels": [
            {
                "panel_id": "synthetic_stage_12_panel",
                "provider": "IGSR_1000_GENOMES",
                "release_id": "synthetic_release_v1",
                "assembly": "GRCh38",
                "sample_count": 3_202,
                "sample_scope": "ALL_RELEASE_SAMPLES",
                "population_scope": "ALL_AVAILABLE_1KGP_POPULATIONS",
                "phased": True,
                "variant_classes": ["SNV", "INDEL"],
                "base_url": base_url,
                "readme_url": f"{base_url}/README",
                "readme_sha256": _sha256(contents["README"]),
                "manifest_url": f"{base_url}/manifest.txt",
                "manifest_sha256": _sha256(contents["manifest.txt"]),
                "vcf_filename_template": "synthetic.chr{chromosome}.panel.vcf.gz",
                "index_filename_template": "synthetic.chr{chromosome}.panel.vcf.gz.tbi",
                "chromosomes": [
                    {
                        "chromosome": chromosome,
                        "vcf_md5": _md5(vcf_content),
                        "index_md5": _md5(index_content),
                    }
                    for chromosome in range(1, 23)
                ],
            }
        ],
    }
    path.write_text(json.dumps(catalog), encoding="utf-8")


def _populate_cache(catalog_path: Path, cache_dir: Path, contents: dict[str, bytes]) -> None:
    resolved = resolve_reference_panel(
        panel_id="synthetic_stage_12_panel",
        assembly="GRCh38",
        chromosome=19,
        catalog_path=catalog_path,
    )

    def downloader(url: str, destination: Path, timeout: int, chunk_size: int) -> None:
        destination.write_bytes(contents[Path(url).name])

    cache_reference_panel(resolved, cache_dir, downloader=downloader)


def _write_fake_bcftools(path: Path, target_confidence: float = 0.95) -> None:
    script = f'''#!{sys.executable}
import pathlib
import shutil
import sys

arguments = sys.argv[1:]
if arguments == ["--version"]:
    print("bcftools 1.21")
    raise SystemExit(0)
if arguments[:2] == ["view", "--header-only"]:
    print("##fileformat=VCFv4.2")
    print("##contig=<ID=chr19,length=58617616>")
    print('##FORMAT=<ID=GT,Number=1,Type=String,Description="Phased genotypes">')
    raise SystemExit(0)
if arguments[:2] == ["query", "--list-samples"]:
    sample_sidecar = pathlib.Path(arguments[-1] + ".samples")
    if sample_sidecar.is_file():
        print(sample_sidecar.read_text(), end="")
        raise SystemExit(0)
    if pathlib.Path(arguments[-1]).name in {
        "study.shapeit5.vcf.gz", "common.phased.bcf", "target.phased.bcf"
    }:
        print("sample_1\\nsample_2\\nsample_3")
        raise SystemExit(0)
    for sample_index in range(1, 3203):
        print(f"HG{{sample_index:05d}}")
    raise SystemExit(0)
if arguments[0] == "view" and "--regions" in arguments:
    pathlib.Path(arguments[arguments.index("--output") + 1]).write_bytes(b"window vcf")
    raise SystemExit(0)
if arguments[:2] == ["index", "--tbi"]:
    pathlib.Path(arguments[-1] + ".tbi").write_bytes(b"tabix index")
    raise SystemExit(0)
if arguments[:2] == ["index", "--nrecords"]:
    print("1")
    raise SystemExit(0)
if arguments[0] == "+fill-tags":
    source = pathlib.Path(arguments[1])
    destination = pathlib.Path(arguments[arguments.index("--output") + 1])
    shutil.copyfile(source, destination)
    for suffix in (".samples", ".variants"):
        source_sidecar = pathlib.Path(str(source) + suffix)
        if source_sidecar.is_file():
            shutil.move(source_sidecar, pathlib.Path(str(destination) + suffix))
    raise SystemExit(0)
if arguments[0] == "query" and "--include" in arguments:
    raise SystemExit(0)
if arguments[0] == "norm":
    pathlib.Path(arguments[arguments.index("--output") + 1]).write_bytes(b"decomposed bcf")
    raise SystemExit(0)
if arguments[0] == "view" and "--types" in arguments:
    pathlib.Path(arguments[arguments.index("--output") + 1]).write_bytes(b"harmonized vcf")
    raise SystemExit(0)
if arguments[0] == "query" and "--format" in arguments:
    if arguments[arguments.index("--format") + 1] == "%INFO/AC\\\\t%INFO/AN\\\\n":
        print("1\\t6")
        print("2\\t6")
        raise SystemExit(0)
    format_text = arguments[arguments.index("--format") + 1]
    input_name = pathlib.Path(arguments[-1]).name
    if "%GT" in format_text and input_name == "study.shapeit5.vcf.gz":
        print("chr19\\t1900\\trs19\\tA\\tG\\t0/1\\t0/0\\t0/0")
        print("chr19\\t100000\\ttarget_GRCh38_1_100000_A_G\\tA\\tG\\t0/1\\t0/0\\t0/0")
        raise SystemExit(0)
    if "%GT" in format_text and input_name == "common.phased.bcf":
        print("chr19\\t1900\\trs19\\tA\\tG\\t1|0\\t0|0\\t0|0")
        raise SystemExit(0)
    if "%PP" in format_text and input_name == "target.phased.bcf":
        print("chr19\\t1900\\trs19\\tA\\tG\\t1|0\\t.\\t0|0\\t.\\t0|0\\t.")
        print("chr19\\t100000\\ttarget_GRCh38_1_100000_A_G\\tA\\tG\\t1|0\\t{target_confidence}\\t0|0\\t.\\t0|0\\t.")
        raise SystemExit(0)
    variant_sidecar = pathlib.Path(arguments[-1] + ".variants")
    if variant_sidecar.is_file():
        print(variant_sidecar.read_text(), end="")
        pathlib.Path(arguments[-1] + ".samples").unlink()
        variant_sidecar.unlink()
        raise SystemExit(0)
    print("chr19\\t1900\\trs19\\tA\\tG")
    raise SystemExit(0)
raise SystemExit(9)
'''
    path.write_text(script, encoding="utf-8")
    path.chmod(0o755)


def _write_shapeit5_plink(path: Path, delegated_plink: Path) -> None:
    script = f'''#!{sys.executable}
import os
import pathlib
import sys

arguments = sys.argv[1:]
if arguments == ["--version"]:
    print("PLINK v1.90 synthetic")
    raise SystemExit(0)
output_prefix = pathlib.Path(arguments[arguments.index("--out") + 1]) if "--out" in arguments else None
if output_prefix is None or output_prefix.name != "study.plink_export":
    os.execv({str(delegated_plink)!r}, [{str(delegated_plink)!r}, *arguments])
base_prefix = pathlib.Path(arguments[arguments.index("--bfile") + 1])
selected_ids = set(pathlib.Path(arguments[arguments.index("--extract") + 1]).read_text().split())
update_rows = [line.split() for line in pathlib.Path(arguments[arguments.index("--update-ids") + 1]).read_text().splitlines()]
reference_rows = [line.split() for line in pathlib.Path(arguments[arguments.index("--a2-allele") + 1]).read_text().splitlines()]
reference_by_id = {{row[1]: row[0] for row in reference_rows}}
variants = [
    row for row in (line.split() for line in base_prefix.with_suffix(".bim").read_text().splitlines())
    if row[1] in selected_ids
]
output_vcf = pathlib.Path(str(output_prefix) + ".vcf.gz")
output_vcf.write_bytes(b"synthetic study vcf")
pathlib.Path(str(output_vcf) + ".samples").write_text(
    "".join(row[2] + "\\n" for row in update_rows)
)
pathlib.Path(str(output_vcf) + ".variants").write_text(
    "".join(
        f"chr{{row[0]}}\\t{{row[3]}}\\t{{row[1]}}\\t{{reference_by_id[row[1]]}}\\t{{row[4] if row[5] == reference_by_id[row[1]] else row[5]}}\\n"
        for row in variants
    )
)
'''
    path.write_text(script, encoding="utf-8")
    path.chmod(0o755)


def _write_fake_shapeit5(path: Path, component: str) -> None:
    script = f'''#!{sys.executable}
import pathlib
import sys

arguments = sys.argv[1:]
if arguments == ["--help"]:
    print("[SHAPEIT5] {component}")
    print("  * Version       : 5.1.1 / commit = test")
    raise SystemExit(0)
output_path = pathlib.Path(arguments[arguments.index("--output") + 1])
output_path.write_bytes(b"synthetic phased bcf")
pathlib.Path(str(output_path) + ".csi").write_bytes(b"synthetic csi")
if "--log" in arguments:
    pathlib.Path(arguments[arguments.index("--log") + 1]).write_text("synthetic common log\\n")
else:
    print("synthetic rare log")
'''
    path.write_text(script, encoding="utf-8")
    path.chmod(0o755)


def _prepare_phase_inputs(
    tmp_path: Path, *, target_confidence: float = 0.95
) -> tuple[Path, Path]:
    config_path, runs_dir = prepare_target_region_inputs(tmp_path)
    config = yaml.safe_load(config_path.read_text(encoding="utf-8"))
    contents = {
        "synthetic.chr19.panel.vcf.gz": b"cached source vcf",
        "synthetic.chr19.panel.vcf.gz.tbi": b"cached source index",
        "README": b"synthetic readme",
        "manifest.txt": b"synthetic manifest",
    }
    catalog_path = tmp_path / "reference_catalog.json"
    cache_dir = tmp_path / "reference_cache"
    bcftools_path = tmp_path / "bcftools"
    _write_synthetic_catalog(catalog_path, contents)
    _populate_cache(catalog_path, cache_dir, contents)
    _write_fake_bcftools(bcftools_path, target_confidence)
    shapeit5_plink_path = tmp_path / "shapeit5_plink"
    _write_shapeit5_plink(shapeit5_plink_path, Path(config["tools"]["plink"]))
    common_path = tmp_path / "SHAPEIT5_phase_common"
    rare_path = tmp_path / "SHAPEIT5_phase_rare"
    _write_fake_shapeit5(common_path, "phase_common")
    _write_fake_shapeit5(rare_path, "phase_rare")
    config["inputs"]["reference_panel_catalog"] = str(catalog_path)
    config["tools"]["bcftools"] = str(bcftools_path)
    config["tools"]["plink"] = str(shapeit5_plink_path)
    config["tools"]["phasing_adapter"] = {
        "adapter_id": "shapeit5_phase_common_rare_v1",
        "phase_common_command": str(common_path),
        "phase_rare_command": str(rare_path),
        "expected_version": "5.1.1",
    }
    config["stages"]["phase_target_region"] = {
        "enabled": True,
        "parameters": {
            "reference_panel_id": "synthetic_stage_12_panel",
            "reference_cache_dir": str(cache_dir),
            "reference_cache_offline": True,
            "reference_lock_timeout_seconds": 10,
            "reference_download_timeout_seconds": 30,
            "reference_download_chunk_size": 1_024,
            "reference_extract_timeout_seconds": 30,
            "reference_extract_threads": 1,
            "minimum_common_variants": 1,
            "shapeit5_input_timeout_seconds": 30,
            "shapeit5_execution_timeout_seconds": 30,
            "shapeit5_threads": 1,
            "shapeit5_seed": 15052011,
            "shapeit5_effective_size": 15000,
            "minimum_phase_confidence": 0.9,
        },
    }
    config_path.write_text(
        yaml.safe_dump(config, allow_unicode=True, sort_keys=False), encoding="utf-8"
    )
    return config_path, runs_dir


def _read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def test_stage_12_publishes_reference_and_harmonization_then_reuses(
    tmp_path: Path,
) -> None:
    config_path, runs_dir = _prepare_phase_inputs(tmp_path)

    run_dir = run_pipeline(config_path, runs_dir)

    stage_dir = run_dir / "stages" / "12_phase_target_region"
    outputs = json.loads((stage_dir / "stage_outputs.json").read_text(encoding="utf-8"))
    assert {artifact["artifact_id"] for artifact in outputs["artifacts"]} == {
        "reference_window_vcf",
        "reference_window_index",
        "reference_window_manifest",
        "harmonized_reference_vcf",
        "harmonized_reference_index",
        "variant_harmonization",
        "reference_harmonization_manifest",
        "shapeit5_study_vcf",
        "shapeit5_study_index",
        "shapeit5_genetic_map",
        "shapeit5_pedigree",
        "shapeit5_variant_selection",
        "shapeit5_sample_mapping",
        "shapeit5_inputs_manifest",
        "shapeit5_common_bcf",
        "shapeit5_common_index",
        "shapeit5_final_bcf",
        "shapeit5_final_index",
        "shapeit5_common_log",
        "shapeit5_rare_log",
        "carrier_haplotypes",
        "phasing_transmissions",
        "phasing_unreliable_regions",
        "shapeit5_phasing_manifest",
        "phasing_qc",
        "carrier_haplotype_summary",
        "phase_target_region_summary",
    }
    rows = _read_tsv(stage_dir / "reference_harmonization" / "variant_harmonization.tsv")
    assert rows[0]["HARMONIZATION_STATUS"] == "MATCHED_DIRECT"
    assert rows[1]["HARMONIZATION_STATUS"] == "STUDY_ONLY"
    audit = json.loads((stage_dir / "audit.json").read_text(encoding="utf-8"))
    assert audit["metrics"]["cache_status"] == "HIT"
    assert audit["metrics"]["target_variant_status"] == "STUDY_ONLY"
    assert audit["counts"]["harmonization_statuses"] == {
        "MATCHED_DIRECT": 1,
        "STUDY_ONLY": 1,
    }
    assert audit["exclusions"] == [
        {"scope": "common_phase", "reason": "STUDY_ONLY", "count": 1}
    ]
    assert audit["warnings"] == [
        {"code": "target_variant_absent_from_reference", "count": 1}
    ]
    shapeit5_manifest = json.loads(
        (stage_dir / "shapeit5_inputs" / "shapeit5_inputs_manifest.json").read_text(
            encoding="utf-8"
        )
    )
    assert shapeit5_manifest["sample_count"] == 3
    assert shapeit5_manifest["study_variant_count"] == 2
    assert shapeit5_manifest["common_variant_count"] == 1
    assert shapeit5_manifest["pedigree_record_count"] == 0
    assert shapeit5_manifest["target_variant_role"] == "RARE_TARGET"
    sample_rows = _read_tsv(
        stage_dir / "shapeit5_inputs" / "shapeit5_sample_mapping.tsv"
    )
    assert [row["VCF_SAMPLE_ID"] for row in sample_rows] == [
        "sample_1",
        "sample_2",
        "sample_3",
    ]
    carrier_rows = _read_tsv(
        stage_dir / "shapeit5_phasing" / "carrier_haplotypes.tsv"
    )
    assert carrier_rows[0]["SAMPLE_ID"] == "sample_1"
    assert carrier_rows[0]["CARRIER_HAPLOTYPE"] == "H1"
    assert carrier_rows[0]["CONFIDENCE_STATUS"] == "SCORED_PASS"
    assert carrier_rows[0]["RELIABILITY_STATUS"] == "PASS"
    phasing_manifest = json.loads(
        (stage_dir / "shapeit5_phasing" / "shapeit5_phasing_manifest.json").read_text(
            encoding="utf-8"
        )
    )
    assert phasing_manifest["carrier_count"] == 1
    assert phasing_manifest["reliable_carrier_count"] == 1
    qc_rows = _read_tsv(stage_dir / "phasing_qc" / "phasing_qc.tsv")
    assert len(qc_rows) == 12
    assert qc_rows[-1]["CHECK_ID"] == "carrier_phase_confidence"
    assert qc_rows[-1]["STATUS"] == "PASS"
    summary = json.loads(
        (stage_dir / "phasing_qc" / "phase_target_region_summary.json").read_text(
            encoding="utf-8"
        )
    )
    assert summary["haplotype_counts"] == {
        "NONE": 2,
        "H1": 1,
        "H2": 0,
        "BOTH": 0,
    }
    assert summary["manual_validation_required"] is False
    assert audit["expected_visualizations"] == [
        "plot_phasing_haplotypes",
        "plot_phasing_transmissions",
        "plot_phasing_confidence",
    ]

    resume_pipeline(run_dir)

    manifest = json.loads((run_dir / "manifest.json").read_text(encoding="utf-8"))
    stage_record = next(
        stage for stage in manifest["stages"] if stage["stage_name"] == "phase_target_region"
    )
    assert stage_record["state"] == "SUCCEEDED"
    assert stage_record["attempt_count"] == 1

    study_bim = run_dir / "stages" / "11_prepare_target_region" / "target_region.bim"
    study_bim.write_text(
        study_bim.read_text(encoding="utf-8") + "19 extra 0 200000 A G\n",
        encoding="utf-8",
    )
    with pytest.raises(IntegrityError):
        resume_pipeline(run_dir)


def test_stage_12_error_mapping_uses_standard_codes(monkeypatch, tmp_path: Path) -> None:
    arguments = [
        "--stage-inputs",
        str(tmp_path / "inputs.json"),
        "--output-dir",
        str(tmp_path / "output"),
    ]

    def raise_external(*args, **kwargs):
        raise ReferenceWindowError("bcftools_norm_failed:8")

    monkeypatch.setattr(phase_target_region, "execute", raise_external)
    assert phase_target_region.main(arguments) == 3

    def raise_shapeit_external(*args, **kwargs):
        raise phase_target_region.Shapeit5ExecutionExternalError(
            "shapeit5_phase_rare_failed:8"
        )

    monkeypatch.setattr(phase_target_region, "execute", raise_shapeit_external)
    assert phase_target_region.main(arguments) == 3

    def raise_block(*args, **kwargs):
        raise ReferenceHarmonizationIntegrityError(
            "target_variant_harmonization_ambiguous"
        )

    monkeypatch.setattr(phase_target_region, "execute", raise_block)
    assert phase_target_region.main(arguments) == 4

    def raise_shapeit_block(*args, **kwargs):
        raise phase_target_region.Shapeit5ExecutionBlockError(
            "phased_target_genotype_discordant"
        )

    monkeypatch.setattr(phase_target_region, "execute", raise_shapeit_block)
    assert phase_target_region.main(arguments) == 4

    def raise_input(*args, **kwargs):
        raise ValueError("invalid_parameter:minimum_common_variants")

    monkeypatch.setattr(phase_target_region, "execute", raise_input)
    assert phase_target_region.main(arguments) == 2


def test_stage_12_7_requires_manual_review_for_low_carrier_confidence(
    tmp_path: Path,
) -> None:
    config_path, runs_dir = _prepare_phase_inputs(
        tmp_path, target_confidence=0.7
    )

    run_dir = run_pipeline(config_path, runs_dir)

    stage_dir = run_dir / "stages" / "12_phase_target_region"
    audit = json.loads((stage_dir / "audit.json").read_text(encoding="utf-8"))
    summary = json.loads(
        (stage_dir / "phasing_qc" / "phase_target_region_summary.json").read_text(
            encoding="utf-8"
        )
    )
    qc_rows = _read_tsv(stage_dir / "phasing_qc" / "phasing_qc.tsv")
    confidence_qc = next(
        row for row in qc_rows if row["CHECK_ID"] == "carrier_phase_confidence"
    )
    assert confidence_qc["STATUS"] == "WARN"
    assert confidence_qc["DETAIL_CODE"] == "carrier_phase_unreliable"
    assert summary["unreliable_carrier_count"] == 1
    assert summary["manual_validation_required"] is True
    assert audit["manual_validation_required"] is True
    assert {warning["code"] for warning in audit["warnings"]} == {
        "target_variant_absent_from_reference",
        "carrier_phase_unreliable",
    }
