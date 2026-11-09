"""D18 : formats réels, aveuglement, manque préservé, chemins et publication."""

from __future__ import annotations

import gzip
import hashlib
import json
from pathlib import Path

import jsonschema
import numpy as np
import pytest

from effet_fondateur.calibration.export_blinded_inputs import (
    BlindedExportError, audit_bundle, blind_panel, export_blinded_inputs,
    load_checked_export, validate_inventory, validate_schema,
)
from effet_fondateur.calibration.report_blinded_export import build_report
from effet_fondateur.calibration.materialize_measurement_errors import materialize_measurement_errors
from effet_fondateur.calibration.simulation.measurement_errors import read_panel
from effet_fondateur.contracts import load_pipeline_config
from effet_fondateur.orchestrator.pipeline import run_pipeline
from effet_fondateur.stages.convert_acpa import _read_source_to_spool, _load_spooled_calls

# Même minuscule contrat synthétique que D17, jamais un artefact réel local.
from test_measurement_errors import _fixture as d17_fixture, _panel


ROOT = Path(__file__).resolve().parents[2]


def _config() -> dict:
    document = json.loads((ROOT / "calibration/configurations/d18_micro_chr5_calibration_v1.json").read_text())
    document["target"] = dict(chromosome=5, position_bp=150000, ref="A", alt="G")
    return document


def _fixture(tmp_path: Path) -> dict:
    arguments = d17_fixture(tmp_path)
    materialize_measurement_errors(**arguments)
    config = tmp_path / "d18_config.json"
    config.write_text(json.dumps(_config()))
    return dict(d17_root=arguments["workspace_root"] / "d17_test", config_path=config,
                workspace_root=arguments["workspace_root"], run_id="d18_test", created_at="2026-10-30T00:00:00Z")


def _manifest_root(arguments: dict) -> Path:
    return arguments["workspace_root"] / "d18_test"


def test_schemas_configuration_and_pipeline_draft_are_valid(tmp_path):
    validate_schema(_config(), "blinded_export_config")
    for name in ("blinded_export_config", "blinded_export_manifest", "blinded_input_inventory", "blinding_map"):
        jsonschema.Draft202012Validator.check_schema(json.loads((ROOT / f"calibration/schemas/{name}.schema.json").read_text()))
    config = _config()
    config["private_truth"] = "forbidden"
    with pytest.raises(BlindedExportError):
        validate_schema(config, "blinded_export_config")


def test_blinding_permutation_is_invertible_and_preserves_missing_calls():
    panel = _panel()
    panel.alleles[0, 0] = (-1, -1)
    before = panel.alleles.copy()
    blinded, mapping = blind_panel(panel, _config(), "B" + "a" * 16)
    for index, item in enumerate(mapping):
        original = panel.samples.index(item["source_sample_id"])
        assert item["blinded_sample_id"] == blinded.samples[index]
        assert np.array_equal(blinded.alleles[:, index], panel.alleles[:, original])
    assert not set(blinded.samples) & set(panel.samples)
    assert np.array_equal(panel.alleles, before)
    other, _ = blind_panel(panel, _config(), "B" + "b" * 16)
    assert not set(other.samples) & set(blinded.samples)


def test_export_roundtrip_inventory_permissions_and_determinism(tmp_path, monkeypatch):
    arguments = _fixture(tmp_path)
    monkeypatch.setattr(gzip.time, "time", lambda: 1000000000)
    manifest = export_blinded_inputs(**arguments)
    root = _manifest_root(arguments)
    assert load_checked_export(root) == manifest
    assert manifest["bundle_count"] == 8 and manifest["file_count"] == 56
    assert all(item["genotype_mismatch_count"] == item["truth_leak_count"] == 0 for item in manifest["audits"])
    for path in root.rglob("*"):
        assert path.stat().st_mode & 0o777 == (0o700 if path.is_dir() else 0o600)
    forbidden = (*_panel().samples, "private_truth", "phase_point_02", "combined", "profile_id", "measurement_error", "blinding_seed")
    for path in (root / "blinded_inputs").rglob("*"):
        if path.is_file():
            content = gzip.open(path, "rt").read() if path.suffix == ".gz" else path.read_text()
            assert not any(token in content for token in forbidden)
    with pytest.raises(BlindedExportError, match="run_exists"):
        export_blinded_inputs(**arguments)
    monkeypatch.setattr(gzip.time, "time", lambda: 1000000002)
    assert export_blinded_inputs(**_fixture(tmp_path / "second")) == manifest


def test_actual_acpa_converter_preserves_all_observed_alleles(tmp_path):
    arguments = _fixture(tmp_path)
    export_blinded_inputs(**arguments)
    root = _manifest_root(arguments)
    mapping = json.loads((root / "private_truth/d18/blinding_map.json").read_text())
    combined = next(item for item in mapping if item["profile_id"] == "combined")
    bundle = root / "blinded_inputs/bundles" / combined["bundle_id"]
    with gzip.open(bundle / "observed.vcf.gz", "rt") as handle:
        observed = read_panel(handle, phased=False)
    total_missing = 0
    for sample_index, sample in enumerate(observed.samples):
        spool = tmp_path / f"spool_{sample_index}.tsv"
        markers = {}
        invalid, skipped, counts, called = _read_source_to_spool(bundle / "acpa_samples" / f"{sample}.txt", spool, "GRCh38", markers)
        assert skipped == 0
        calls = _load_spooled_calls(spool)
        for index, variant in enumerate(observed.variants):
            pair = observed.alleles[index, sample_index]
            expected = ("0", "0") if pair[0] < 0 else tuple((variant[3], variant[4])[int(allele)] for allele in pair)
            assert calls[variant[2]] == expected
        assert invalid == np.count_nonzero(observed.alleles[:, sample_index, 0] < 0)
        total_missing += invalid
    assert total_missing > 0


def test_generated_configuration_runs_only_00_to_02_without_autoapproval(tmp_path, monkeypatch):
    arguments = _fixture(tmp_path)
    export_blinded_inputs(**arguments)
    root = _manifest_root(arguments)
    inventory = json.loads((root / "blinded_inputs/inventory.json").read_text())
    bundle = root / "blinded_inputs/bundles" / inventory["bundles"][0]["bundle_id"]
    config = load_pipeline_config(bundle / "pipeline.yaml")
    assert set(config["stages"]) == {"initialize_run", "validate_sources", "build_sample_registry"}
    monkeypatch.chdir(bundle)
    run = run_pipeline(bundle / "pipeline.yaml", tmp_path / "pipeline_runs")
    manifest = json.loads((run / "manifest.json").read_text())
    assert {stage["state"] for stage in manifest["stages"]} == {"SUCCEEDED"}
    assert "sample_registry_approval" in manifest["manual_decisions_required"]
    assert not config["policies"]["allow_clinical_genotype_inference"]


@pytest.mark.parametrize("mutation", ["role", "seed", "target", "alleles", "path", "source_hash", "symlink", "hardlink"])
def test_input_errors_are_refused_before_publication(tmp_path, mutation):
    arguments = _fixture(tmp_path)
    config = _config()
    if mutation == "role": config["panel_role"] = "VALIDATION"
    elif mutation == "seed": config["blinding_seed"] = 1005
    elif mutation == "target": config["target"]["position_bp"] += 1
    elif mutation == "alleles": config["target"]["alt"] = "T"
    elif mutation == "path": arguments["run_id"] = "../escape"
    else:
        source = arguments["d17_root"] / "private_truth/d17/baseline.phased.vcf.gz"
        if mutation == "source_hash": source.write_bytes(b"tampered")
        elif mutation == "symlink":
            source.rename(source.with_suffix(".backup"))
            source.symlink_to(source.with_suffix(".backup").name)
        else: (source.parent / "alias.gz").hardlink_to(source)
    arguments["config_path"].write_text(json.dumps(config))
    with pytest.raises(BlindedExportError): export_blinded_inputs(**arguments)
    assert not _manifest_root(arguments).exists()


@pytest.mark.parametrize("mutation", ["extra", "empty_directory", "hash", "inventory_truth", "inventory_role", "private_map"])
def test_final_export_rejects_tampering_or_uninventoried_files(tmp_path, mutation):
    arguments = _fixture(tmp_path)
    export_blinded_inputs(**arguments)
    root = _manifest_root(arguments)
    inventory_path = root / "blinded_inputs/inventory.json"
    inventory = json.loads(inventory_path.read_text())
    if mutation == "extra": (root / "blinded_inputs/private_truth.json").write_text("leak")
    elif mutation == "empty_directory": (root / "blinded_inputs/private_truth").mkdir()
    elif mutation == "hash": (root / "blinded_inputs" / inventory["artifacts"][0]["path"]).write_text("tampered")
    elif mutation == "private_map": (root / "private_truth/d18/blinding_map.json").write_text("[]")
    else:
        if mutation == "inventory_truth": inventory["seed"] = 1006
        else: inventory["artifacts"][0]["path"] = "../private_truth/leak.json"
        with pytest.raises(BlindedExportError): validate_inventory(inventory)
        return
    with pytest.raises(BlindedExportError): load_checked_export(root)


def test_export_never_decodes_candidate_phase_or_error_registry(tmp_path, monkeypatch):
    arguments = _fixture(tmp_path)
    original = gzip.open
    def guarded(path, mode="rb", *args, **kwargs):
        if "t" in mode and ("candidate.phased" in str(path) or "measurement_errors.jsonl" in str(path) or "baseline.phased" in str(path)):
            raise AssertionError("private_genotype_truth_decoded")
        return original(path, mode, *args, **kwargs)
    monkeypatch.setattr(gzip, "open", guarded)
    export_blinded_inputs(**arguments)


def test_failed_audit_removes_temporary_run(tmp_path, monkeypatch):
    arguments = _fixture(tmp_path)
    def fail(*args, **kwargs): raise BlindedExportError("synthetic_audit_failure")
    monkeypatch.setattr("effet_fondateur.calibration.export_blinded_inputs.audit_bundle", fail)
    with pytest.raises(BlindedExportError, match="synthetic_audit_failure"):
        export_blinded_inputs(**arguments)
    assert not _manifest_root(arguments).exists()
    assert not list(arguments["workspace_root"].glob(".d18_test.*"))


@pytest.mark.parametrize("mutation", ["vcf", "metadata", "acpa", "config"])
def test_bundle_audit_detects_mismatch_and_truth_injection(tmp_path, mutation):
    arguments = _fixture(tmp_path)
    export_blinded_inputs(**arguments)
    root = _manifest_root(arguments)
    inventory = json.loads((root / "blinded_inputs/inventory.json").read_text())
    bundle = root / "blinded_inputs/bundles" / inventory["bundles"][0]["bundle_id"]
    with gzip.open(bundle / "observed.vcf.gz", "rt") as handle:
        expected = read_panel(handle, phased=False)
    if mutation == "vcf":
        with gzip.open(bundle / "observed.vcf.gz", "wt") as handle:
            handle.write("invalid")
    elif mutation == "metadata":
        path = bundle / "metadata/samples.tsv"
        path.write_text(path.read_text().replace("UNKNOWN", "MALE", 1))
    elif mutation == "acpa":
        path = next((bundle / "acpa_samples").glob("*.txt"))
        lines = path.read_text().splitlines()
        fields = lines[3].split("\t")
        fields[1] = "AT"
        lines[3] = "\t".join(fields)
        path.write_text("\n".join(lines) + "\n")
    else:
        path = bundle / "pipeline.yaml"
        path.write_text(path.read_text().replace("acpa_samples_dir: acpa_samples", "acpa_samples_dir: private_truth"))
    with pytest.raises((BlindedExportError, ValueError)):
        audit_bundle(bundle, expected, _config())


def test_uninventoried_truth_is_blocked_before_atomic_publication(tmp_path, monkeypatch):
    from effet_fondateur.calibration import export_blinded_inputs as module
    arguments = _fixture(tmp_path)
    original = module._write_bundle
    def inject(root, panel, config):
        original(root, panel, config)
        (root / "truth.json").write_text("leak")
    monkeypatch.setattr(module, "_write_bundle", inject)
    with pytest.raises(BlindedExportError, match="uninventoried_file"):
        export_blinded_inputs(**arguments)
    assert not _manifest_root(arguments).exists()


def test_report_is_aggregate_and_independent_of_private_sample_mapping(tmp_path):
    arguments = _fixture(tmp_path)
    export_blinded_inputs(**arguments)
    root = _manifest_root(arguments)
    template = tmp_path / "template.html"
    template.write_text('<html><head><title>Template</title></head><body></body></html>')
    output = build_report(run_root=root, template_path=template, output_root=arguments["workspace_root"] / "report")
    html = output.read_text()
    assert html.count('<section id="s') == 14 and html.count('class="toc-entry"') == 14
    assert "NOT_EVALUABLE" in html and "calibrated=false" in html
    assert "Figure 1" in html and "Figure 2" in html and "Tableau 1" in html
    data = (output.parent / "report_data.json").read_text()
    mapping = json.loads((root / "private_truth/d18/blinding_map.json").read_text())
    for group in mapping:
        for sample in group["samples"]:
            assert sample["source_sample_id"] not in html + data
            assert sample["blinded_sample_id"] not in html + data
    with pytest.raises(BlindedExportError, match="report_output_invalid"):
        build_report(run_root=root, template_path=template, output_root=output.parent)


@pytest.mark.parametrize("mutation", ["duplicate", "sample_count", "sample_order", "seed"])
def test_private_correspondence_and_configuration_are_cross_checked(tmp_path, mutation):
    arguments = _fixture(tmp_path)
    export_blinded_inputs(**arguments)
    root = _manifest_root(arguments)
    private = root / "private_truth/d18"
    mapping = json.loads((private / "blinding_map.json").read_text())
    manifest_path = private / "blinded_export_manifest.json"
    manifest = json.loads(manifest_path.read_text())
    if mutation == "duplicate":
        mapping[0]["samples"][0]["blinded_sample_id"] = mapping[1]["samples"][0]["blinded_sample_id"]
    elif mutation == "sample_count":
        mapping[0]["samples"].pop()
    elif mutation == "sample_order":
        mapping[0]["samples"].reverse()
    else:
        manifest["blinding_seed"] += 1
    path = private / "blinding_map.json"
    path.write_text(json.dumps(mapping))
    manifest["blinding_map_sha256"] = hashlib.sha256(path.read_bytes()).hexdigest()
    manifest_path.write_text(json.dumps(manifest))
    with pytest.raises(BlindedExportError):
        load_checked_export(root)


def test_source_header_and_info_truth_annotations_are_not_copied(tmp_path):
    arguments = _fixture(tmp_path)
    manifest_path = arguments["d17_root"] / "private_truth/d17/measurement_error_manifest.json"
    manifest = json.loads(manifest_path.read_text())
    artifact = next(item for item in manifest["artifacts"] if item["role"] == "OBSERVED_MEASUREMENT_PANEL")
    path = arguments["d17_root"] / artifact["path"]
    with gzip.open(path, "rt") as handle:
        text = handle.read()
    text = text.replace("##fileformat=VCFv4.2", "##fileformat=VCFv4.2\n##private_truth=secret_profile")
    text = text.replace("\tPASS\t.\tGT", "\tPASS\ttrue_age=100\tGT")
    with gzip.open(path, "wt") as handle:
        handle.write(text)
    artifact["sha256"] = hashlib.sha256(path.read_bytes()).hexdigest()
    artifact["size_bytes"] = path.stat().st_size
    manifest_path.write_text(json.dumps(manifest))
    export_blinded_inputs(**arguments)
    for path in (_manifest_root(arguments) / "blinded_inputs").rglob("*.gz"):
        with gzip.open(path, "rt") as handle:
            exported = handle.read()
        assert "secret_profile" not in exported and "true_age" not in exported
