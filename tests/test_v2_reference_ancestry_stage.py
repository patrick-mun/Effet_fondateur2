import json
import csv
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pytest
import yaml

from effet_fondateur.ancestry import GenotypePanel, ReferenceSample, Variant, harmonize_alt_dosages
from effet_fondateur.audit import atomic_write_json, sha256_file
from effet_fondateur.contracts import build_file_artifact, validate_json_document, validate_tsv_table
from effet_fondateur.orchestrator.pipeline import DEFAULT_STAGE_DEFINITIONS, build_stage_catalog
from effet_fondateur.stages.analyze_reference_ancestry import (
    AUDIT_COLUMNS,
    CENTROID_COLUMNS,
    EIGENVALUE_COLUMNS,
    LOADING_COLUMNS,
    SCORE_COLUMNS,
    _centroid_rows,
    _eigenvalue_rows,
    _loading_rows,
    _parameters,
    _require_complete_local_target,
    _score_rows,
    _variant_audit_rows,
    _write_tsv,
    execute,
)
import effet_fondateur.stages.analyze_reference_ancestry as ancestry_stage


def _result():
    variants = tuple(Variant("chr5", 100 + index, ".", "A", "G") for index in range(3))
    reference = GenotypePanel(
        ("R1:H1", "R1:H2", "R2:H1", "R2:H2"),
        variants,
        np.array([[0, 0, 1], [0, 1, 0], [1, 1, 0], [1, 0, 1]], dtype=float),
        1,
    )
    study = GenotypePanel(
        ("S1:H1", "S1:H2"), variants, np.array([[0, 0, 1], [1, 1, 0]], dtype=float), 1
    )
    return reference, study, harmonize_alt_dosages(
        reference, study, requested_components=2, minimum_variants=3
    )


def test_stage_16a_is_registered_after_roh() -> None:
    catalog = build_stage_catalog(DEFAULT_STAGE_DEFINITIONS)
    stage = catalog["analyze_reference_ancestry"]
    assert stage.stage_id == "16A"
    assert "analyze_roh" in stage.dependencies
    assert stage.config_input_files == (
        "target_variant_metadata",
        "reference_panel_catalog",
        "ancestry_reference_catalog",
    )


def test_stage_16a_parameters_cover_global_local_and_offline_cache() -> None:
    parameters = _parameters({"ancestry_cache_offline": True})
    assert parameters["global_requested_components"] == 10
    assert parameters["local_requested_components"] == 10
    assert parameters["ancestry_cache_offline"] is True


def test_stage_16a_versioned_tables_validate(tmp_path: Path) -> None:
    reference, study, result = _result()
    metadata = {
        "R1": ReferenceSample("R1", "P1", "SP"),
        "R2": ReferenceSample("R2", "P2", "SP"),
    }
    status = {"S1:H1": "CARRIER_COPY", "S1:H2": "NON_CARRIER_COPY"}
    specifications = (
        ("scores.tsv", SCORE_COLUMNS, _score_rows("LOCAL", result, metadata, status), "ancestry_scores.schema.json"),
        ("eigen.tsv", EIGENVALUE_COLUMNS, _eigenvalue_rows("LOCAL", result), "ancestry_eigenvalues.schema.json"),
        ("loadings.tsv", LOADING_COLUMNS, _loading_rows("LOCAL", result), "ancestry_variant_loadings.schema.json"),
        ("audit.tsv", AUDIT_COLUMNS, _variant_audit_rows("LOCAL", reference, study, result, 0.95), "ancestry_variant_audit.schema.json"),
        ("centroids.tsv", CENTROID_COLUMNS, _centroid_rows("LOCAL", result, metadata), "ancestry_population_centroids.schema.json"),
    )
    for filename, columns, rows, schema in specifications:
        path = tmp_path / filename
        _write_tsv(path, columns, rows)
        assert validate_tsv_table(path, schema).row_count > 0


def test_stage_16a_summary_for_arbitrary_target_validates() -> None:
    summary = {
        "schema_version": "1.0.0",
        "method_id": "reference_only_global_local_pca_v1",
        "assembly": "GRCh38",
        "target": {"variant_id": "configured_variant", "chromosome": 5, "position_bp": 102, "ref": "A", "alt": "G"},
        "global": {"reference_entity_count": 2504, "study_entity_count": 4, "candidate_variant_count": 100, "informative_variant_count": 90, "component_count": 10},
        "local": {"reference_entity_count": 5008, "study_entity_count": 8, "candidate_variant_count": 20, "informative_variant_count": 18, "component_count": 10, "region_start_bp": 90, "region_end_bp": 120},
        "cache": {"metadata_status": "HIT", "extract_hits": 22, "extract_populated": 0, "offline": True},
        "interpretation": {"policy": "RELATIVE_REFERENCE_POSITIONING_ONLY", "ethnic_identity_assigned": False, "genealogical_ancestor_identified": False, "local_ancestry_proven": False, "ibd_proven": False},
        "checks": {"reference_unrelated_only": "PASS", "study_not_used_for_axes": "PASS", "global_local_separated": "PASS", "target_and_region_resolved": "PASS", "cache_integrity": "PASS", "variant_harmonization": "PASS"},
    }
    validate_json_document(summary, "reference_ancestry_summary.schema.json")


def test_local_target_must_be_present_and_complete() -> None:
    target = {
        "project_variant_id": "configured_variant", "chromosome": 5,
        "position_bp": 102, "ref": "A", "alt": "G",
    }
    variant = Variant("chr5", 102, "configured_variant", "A", "G")
    complete = GenotypePanel(
        ("S1:H1", "S1:H2"), (variant,), np.array([[1.0], [0.0]]), 1
    )
    assert _require_complete_local_target(complete, target) == variant
    absent = GenotypePanel(
        ("S1:H1", "S1:H2"), (Variant("chr5", 103, "other", "A", "G"),),
        np.array([[1.0], [0.0]]), 1,
    )
    with pytest.raises(ValueError, match="local_target_missing"):
        _require_complete_local_target(absent, target)
    missing = GenotypePanel(
        ("S1:H1", "S1:H2"), (variant,), np.array([[np.nan], [0.0]]), 1
    )
    with pytest.raises(ValueError, match="local_target_genotype_missing"):
        _require_complete_local_target(missing, target)


def test_all_16a_json_schemas_are_valid_json() -> None:
    for path in Path("schemas").glob("*ancestry*.schema.json"):
        json.loads(path.read_text(encoding="utf-8"))


def test_stage_16a_executes_both_scopes_with_synthetic_external_tools(
    tmp_path: Path, monkeypatch
) -> None:
    run_dir = tmp_path / "run"
    output_dir = run_dir / "stages" / ".16A_attempt"
    inputs_dir = run_dir / "synthetic_inputs"
    output_dir.mkdir(parents=True)
    inputs_dir.mkdir()
    target = {
        "schema_version": "1.0.0", "assembly": "GRCh38", "gene": "GENE5",
        "chromosome": 5, "position_bp": 102, "ref": "A", "alt": "G",
        "transcript": "NM_SYNTHETIC.1", "hgvs_c": "c.1A>G", "hgvs_p": None,
        "project_variant_id": "configured_variant", "rsid": None,
        "reference_validation": {"status": "CONFIRMED", "source_id": "synthetic", "confirmed_at": "2026-08-12T00:00:00Z"},
        "annotation_validation": {"status": "CONFIRMED", "source_id": "synthetic", "confirmed_at": "2026-08-12T00:00:00Z"},
    }
    target_path = inputs_dir / "target.yaml"
    target_path.write_text(yaml.safe_dump(target), encoding="utf-8")
    reference_catalog = Path("config/references/1000g_high_coverage_grch38_phased.json").resolve()
    ancestry_catalog = Path("config/references/1000g_ancestry_grch38.json").resolve()

    bed = inputs_dir / "kinship_panel.bed"; bed.write_bytes(b"\x6c\x1b\x01synthetic")
    bim = inputs_dir / "kinship_panel.bim"; bim.write_text("5 v1 0 100 G A\n5 v2 0 101 G A\n5 v3 0 103 G A\n", encoding="utf-8")
    fam = inputs_dir / "kinship_panel.fam"; fam.write_text("F I 0 0 1 -9\n", encoding="utf-8")
    descriptor = {
        "schema_version": "1.0.0", "dataset_id": "kinship_panel", "assembly": "GRCh38",
        "scope": "autosomal_genomewide", "target_chromosome": None,
        "sample_set_id": "study", "variant_set_id": "global", "sample_count": 1,
        "variant_count": 3, "marker_mode": "intersection", "source_format": "PLINK_KINSHIP_PANEL",
        "allele_orientation": "forward_strand_calls",
        "files": {suffix: {"name": path.name, "sha256": sha256_file(path)} for suffix, path in (("bed", bed), ("bim", bim), ("fam", fam))},
    }
    dataset = inputs_dir / "kinship_panel.dataset.json"; atomic_write_json(dataset, descriptor)
    samples = inputs_dir / "samples.master.tsv"
    sample_columns = ["SAMPLE_ID", "SOURCE_FILE", "FID", "IID", "PID", "MID", "SEX", "CLINICAL_STATUS", "GROUP_LABEL", "TARGET_GENOTYPE", "TARGET_GENOTYPE_SOURCE", "ARRAY_BATCH", "INCLUDE_GENOMEWIDE", "INCLUDE_TARGET_CHROMOSOME", "NOTES_CODE"]
    with samples.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n"); writer.writerow(sample_columns)
        writer.writerow(["S1", "synthetic.txt", "F", "I", "0", "0", "MALE", "UNKNOWN", "STUDY", "A/G", "MOLECULAR", "B1", "true", "true", ""])
    validate_tsv_table(samples, "samples_master.schema.json")

    final_bcf = inputs_dir / "target.phased.bcf"; final_bcf.write_bytes(b"synthetic bcf")
    final_index = inputs_dir / "target.phased.bcf.csi"; final_index.write_bytes(b"synthetic index")
    reference_vcf = inputs_dir / "reference.vcf.gz"; reference_vcf.write_bytes(b"synthetic reference")
    reference_index = inputs_dir / "reference.vcf.gz.tbi"; reference_index.write_bytes(b"synthetic index")
    carriers = inputs_dir / "carrier_haplotypes.tsv"
    carrier_columns = ["SAMPLE_ORDER", "SAMPLE_ID", "TARGET_VARIANT_ID", "EXPLICIT_GENOTYPE", "PHASED_GT", "ALT_COPY_COUNT", "CARRIER_HAPLOTYPE", "PHASE_CONFIDENCE", "CONFIDENCE_STATUS", "RELIABILITY_STATUS", "UNRELIABLE_REASON"]
    with carriers.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n"); writer.writerow(carrier_columns)
        writer.writerow([1, "S1", "configured_variant", "A/G", "1|0", 1, "H1", 0.99, "SCORED_PASS", "PASS", ""])
    harmonization = inputs_dir / "reference_harmonization_manifest.json"
    file_record = lambda name: {"filename": name, "sha256": "aa" * 32, "size_bytes": 1}
    atomic_write_json(harmonization, {
        "schema_version": "1.0.0", "created_at": "2026-08-12T00:00:00Z", "method_id": "canonical_coordinate_allele_harmonization_v1",
        "assembly": "GRCh38", "chromosome": 5, "region": {"start_bp": 90, "end_bp": 110}, "target_variant_id": "configured_variant",
        "sample_count": 3202, "reference_variant_count": 3, "study_variant_count": 4, "common_variant_count": 3,
        "source": {"window_manifest_sha256": "aa" * 32, "vcf_sha256": "bb" * 32, "index_sha256": "cc" * 32},
        "study": {"phasing_manifest_sha256": "dd" * 32, "bim_sha256": "ee" * 32, "variant_set_id": "local"},
        "tool": {"name": "bcftools", "version": "synthetic"},
        "files": {"vcf": file_record("reference.vcf.gz"), "index": file_record("reference.vcf.gz.tbi"), "harmonization": file_record("audit.tsv")},
        "checks": {"input_integrity": "PASS", "assembly_and_region": "PASS", "sample_set_and_order": "PASS", "phased_genotypes": "PASS", "biallelic_sequence_variants": "PASS", "canonical_variant_identity": "PASS", "target_variant": "STUDY_ONLY", "minimum_common_variants": "PASS", "tabix_index": "PASS"},
    })
    roh = inputs_dir / "roh_analysis_summary.json"
    atomic_write_json(roh, {"schema_version": "1.0.0", "method_id": "plink19_array_roh_secondary_v1", "analysis_role": "SECONDARY_DESCRIPTIVE_AUTOZYGOSITY", "scope_statuses": {"GLOBAL": "EVALUATED"}, "genomewide_common_variant_count": 3, "target_chromosome_variant_count": 4, "target_variant_id": "configured_variant", "target_in_roh_count": 0, "f_roh_calculated": False, "feeds_founder_haplotype": False, "feeds_variant_age": False, "consumes_local_ld": False})

    config = yaml.safe_load(Path("config/pipeline.example.yaml").read_text(encoding="utf-8"))
    config["target"].update({"gene": "GENE5", "chromosome": 5, "position_bp": 102, "ref": "A", "alt": "G", "transcript": "NM_SYNTHETIC.1", "project_variant_id": "configured_variant", "rsid": None})
    (run_dir / "config.resolved.yaml").write_text(yaml.safe_dump(config, sort_keys=False), encoding="utf-8")

    source_paths = {
        "config_input_target_variant_metadata": target_path,
        "config_input_reference_panel_catalog": reference_catalog,
        "config_input_ancestry_reference_catalog": ancestry_catalog,
        "samples_master": samples, "kinship_panel_bed": bed, "kinship_panel_bim": bim,
        "kinship_panel_fam": fam, "kinship_panel_dataset": dataset,
        "shapeit5_final_bcf": final_bcf, "shapeit5_final_index": final_index,
        "carrier_haplotypes": carriers, "harmonized_reference_vcf": reference_vcf,
        "harmonized_reference_index": reference_index,
        "reference_harmonization_manifest": harmonization, "roh_analysis_summary": roh,
    }
    artifacts = [build_file_artifact(
        physical_path=path, published_path=str(path), artifact_id=artifact_id,
        artifact_type=artifact_id, media_type="application/octet-stream",
        producer_stage="synthetic", producer_signature="ab" * 32,
        schema_name=None, schema_version=None, assembly="GRCh38",
    ) for artifact_id, path in source_paths.items()]
    stage_inputs = {"schema_version": "1.0.0", "run_id": "synthetic_16a", "stage_id": "16A", "stage_name": "analyze_reference_ancestry", "signature": "16" * 32, "attempt_number": 1, "published_output_dir": "stages/16A_analyze_reference_ancestry", "parameters": {"global_requested_components": 2, "local_requested_components": 2, "minimum_global_variants": 3, "minimum_local_variants": 3}, "artifacts": artifacts}
    stage_inputs_path = output_dir / "stage_inputs.json"; atomic_write_json(stage_inputs_path, stage_inputs)

    reference_ids = tuple(f"R{index:04d}" for index in range(2504))
    reference_metadata = tuple(ReferenceSample(sample, "POP", "SUPER") for sample in reference_ids)
    global_variants = tuple(Variant("chr5", position, f"v{position}", "A", "G") for position in (100, 101, 103))
    global_values = np.asarray([[index % 2, (index // 2) % 2, index % 2] for index in range(2504)], dtype=float)
    global_reference = GenotypePanel(reference_ids, global_variants, global_values, 2)
    global_study = GenotypePanel(("S1",), global_variants, np.array([[1, 0, 1]], dtype=float), 2)
    local_reference_ids = tuple(f"{sample}:{hap}" for sample in reference_ids for hap in ("H1", "H2"))
    local_variants = tuple(Variant("chr5", position, f"local{position}", "A", "G") for position in (100, 101, 103))
    local_values = np.asarray([[index % 2, (index // 2) % 2, index % 2] for index in range(5008)], dtype=float)
    local_reference = GenotypePanel(local_reference_ids, local_variants, local_values, 1)
    local_study_variants = (Variant("chr5", 102, "configured_variant", "A", "G"), *local_variants)
    local_study = GenotypePanel(("S1:H1", "S1:H2"), local_study_variants, np.array([[1, 0, 1, 0], [0, 1, 0, 1]], dtype=float), 1)
    cached_vcf = inputs_dir / "cached.vcf.gz"; cached_vcf.write_bytes(b"cache")

    monkeypatch.setattr(ancestry_stage, "_resolve_tool", lambda command, name: name)
    monkeypatch.setattr(ancestry_stage, "_tool_version", lambda executable: "synthetic")
    monkeypatch.setattr(ancestry_stage, "_run", lambda *args, **kwargs: None)
    monkeypatch.setattr(ancestry_stage, "cache_ancestry_metadata", lambda **kwargs: SimpleNamespace(status="HIT"))
    monkeypatch.setattr(ancestry_stage, "load_reference_samples", lambda cached: reference_metadata)
    monkeypatch.setattr(ancestry_stage, "cache_reference_extract", lambda **kwargs: SimpleNamespace(status="HIT", vcf_path=cached_vcf))
    monkeypatch.setattr(ancestry_stage, "read_plink_raw_panel", lambda *args, **kwargs: global_study)
    monkeypatch.setattr(ancestry_stage, "_samples", lambda bcftools, path, timeout: reference_ids if path in {cached_vcf, reference_vcf} else ("S1",))
    monkeypatch.setattr(ancestry_stage, "_query_panel", lambda bcftools, path, sample_ids, output_path, timeout, **kwargs: global_reference if path == cached_vcf else local_reference if path == reference_vcf else local_study)

    assert execute(stage_inputs_path, output_dir) == 0
    summary = json.loads((output_dir / "ancestry" / "reference_ancestry_summary.json").read_text(encoding="utf-8"))
    assert summary["target"]["chromosome"] == 5
    assert summary["global"]["reference_entity_count"] == 2504
    assert summary["local"]["reference_entity_count"] == 5008
    assert summary["checks"]["study_not_used_for_axes"] == "PASS"
    assert summary["interpretation"]["ethnic_identity_assigned"] is False
    outputs = json.loads((output_dir / "stage_outputs.json").read_text(encoding="utf-8"))
    assert len(outputs["artifacts"]) == 6
