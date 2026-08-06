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


def _write_fake_bcftools(path: Path) -> None:
    script = f'''#!{sys.executable}
import pathlib
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
if arguments[0] == "query" and "--include" in arguments:
    raise SystemExit(0)
if arguments[0] == "norm":
    pathlib.Path(arguments[arguments.index("--output") + 1]).write_bytes(b"decomposed bcf")
    raise SystemExit(0)
if arguments[0] == "view" and "--types" in arguments:
    pathlib.Path(arguments[arguments.index("--output") + 1]).write_bytes(b"harmonized vcf")
    raise SystemExit(0)
if arguments[0] == "query" and "--format" in arguments:
    print("chr19\\t1900\\trs19\\tA\\tG")
    raise SystemExit(0)
raise SystemExit(9)
'''
    path.write_text(script, encoding="utf-8")
    path.chmod(0o755)


def _prepare_phase_inputs(tmp_path: Path) -> tuple[Path, Path]:
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
    _write_fake_bcftools(bcftools_path)
    config["inputs"]["reference_panel_catalog"] = str(catalog_path)
    config["tools"]["bcftools"] = str(bcftools_path)
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

    def raise_block(*args, **kwargs):
        raise ReferenceHarmonizationIntegrityError(
            "target_variant_harmonization_ambiguous"
        )

    monkeypatch.setattr(phase_target_region, "execute", raise_block)
    assert phase_target_region.main(arguments) == 4

    def raise_input(*args, **kwargs):
        raise ValueError("invalid_parameter:minimum_common_variants")

    monkeypatch.setattr(phase_target_region, "execute", raise_input)
    assert phase_target_region.main(arguments) == 2
