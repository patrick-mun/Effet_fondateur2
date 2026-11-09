"""Sélection signée de référence : outils simulés, aucun VCF réel."""

from dataclasses import replace
import hashlib
import json
from pathlib import Path
import subprocess
import yaml

import pytest

from effet_fondateur.audit import sha256_file
from effet_fondateur.contracts import DocumentValidationError, validate_json_document
from effet_fondateur.orchestrator.pipeline import PHASE_AUTOSOMAL_PANEL_STAGE, PHASE_TARGET_REGION_STAGE
from effet_fondateur.references.selection import load_reference_selection, stage_reference_selection
from effet_fondateur.references.window import ReferenceWindowIntegrityError, extract_reference_window
from test_v2_reference_window import FakeBcftools, _cached_reference


def selection_file(tmp_path: Path, samples: tuple[str, ...] = ("HG00001",)) -> Path:
    document = dict(schema_version="1.0.0", method_id="explicit_partition_reference_selection_v1", assembly="GRCh38",
                    source_panel_id="synthetic_panel_v1", partition_manifest_sha256="a" * 64,
                    membership_sha256=hashlib.sha256("".join(sample + "\n" for sample in samples).encode()).hexdigest(),
                    sample_count=len(samples), sample_ids=list(samples))
    path = tmp_path / "selection.json"
    path.write_text(json.dumps(document))
    return path


def test_loads_canonical_selection_without_exporting_identifiers(tmp_path: Path) -> None:
    path = selection_file(tmp_path)
    selection = load_reference_selection(path, expected_sha256=sha256_file(path), panel_id="synthetic_panel_v1")
    assert selection.checked_samples(["HG00002", "HG00001"]) == ["HG00001"]
    assert "HG00001" not in json.dumps(selection.audit_summary())
    assert selection.audit_summary()["sample_count"] == 1


@pytest.mark.parametrize("field,value", [("sample_count", 2), ("membership_sha256", "0" * 64),
                                         ("sample_ids", ["HG00001", "HG00001"]), ("extra", True)])
def test_refuses_changed_content_or_unknown_keys(tmp_path: Path, field: str, value: object) -> None:
    path = selection_file(tmp_path)
    document = json.loads(path.read_text())
    document[field] = value
    path.write_text(json.dumps(document))
    with pytest.raises((ValueError, DocumentValidationError)):
        load_reference_selection(path, expected_sha256=sha256_file(path), panel_id="synthetic_panel_v1")


def test_refuses_wrong_panel_unordered_list_and_modified_file(tmp_path: Path) -> None:
    path = selection_file(tmp_path, ("HG00002", "HG00001"))
    with pytest.raises(ValueError, match="content_invalid"):
        load_reference_selection(path, expected_sha256=sha256_file(path), panel_id="synthetic_panel_v1")
    path = selection_file(tmp_path)
    with pytest.raises(ValueError, match="content_invalid"):
        load_reference_selection(path, expected_sha256=sha256_file(path), panel_id="other")
    with pytest.raises(ValueError, match="integrity_invalid"):
        load_reference_selection(path, expected_sha256="0" * 64, panel_id="synthetic_panel_v1")


def test_unknown_membership_is_not_silently_intersected(tmp_path: Path) -> None:
    path = selection_file(tmp_path)
    selection = load_reference_selection(path, expected_sha256=sha256_file(path), panel_id="synthetic_panel_v1")
    for available in (["HG00002"], ["HG00001", "HG00001"]):
        with pytest.raises(ValueError, match="not_in_source"):
            selection.checked_samples(available)
    with pytest.raises(ValueError, match="not_in_source"):
        replace(selection, sample_ids=()).checked_samples(["HG00001"])


def test_optional_input_is_declared_by_both_stage_definitions() -> None:
    for definition in (PHASE_AUTOSOMAL_PANEL_STAGE, PHASE_TARGET_REGION_STAGE):
        assert definition.optional_config_input_files == ("phasing_reference_selection",)


def test_optional_file_is_hashed_by_orchestrator_and_cannot_disappear(tmp_path: Path) -> None:
    from effet_fondateur.orchestrator.inputs import resolve_stage_input_artifacts
    from effet_fondateur.orchestrator.models import StageDefinition
    from effet_fondateur.orchestrator.errors import StageExecutionError
    definition = StageDefinition("test", "test", "synthetic", optional_config_input_files=("phasing_reference_selection",))
    config = dict(inputs={}, project=dict(assembly="GRCh38"))
    assert resolve_stage_input_artifacts(config, definition, "a" * 64, tmp_path) == []
    path = selection_file(tmp_path)
    config["inputs"]["phasing_reference_selection"] = str(path)
    first = resolve_stage_input_artifacts(config, definition, "a" * 64, tmp_path)
    path.write_text(path.read_text() + "\n")
    second = resolve_stage_input_artifacts(config, definition, "a" * 64, tmp_path)
    assert first[0]["artifact_id"] == "config_input_phasing_reference_selection"
    assert first[0]["sha256"] != second[0]["sha256"]
    config["inputs"]["phasing_reference_selection"] = str(tmp_path / "absent.json")
    with pytest.raises(StageExecutionError, match="source invalide"):
        resolve_stage_input_artifacts(config, definition, "a" * 64, tmp_path)


def test_stage_requires_config_and_signed_artifact_without_silent_fallback(tmp_path: Path) -> None:
    inputs = dict(parameters={}, artifacts=[])
    config = dict(inputs={})
    assert stage_reference_selection(inputs, config, tmp_path, "synthetic_panel_v1") is None
    inputs["parameters"]["require_reference_sample_selection"] = True
    with pytest.raises(ValueError, match="missing_or_ambiguous"):
        stage_reference_selection(inputs, config, tmp_path, "synthetic_panel_v1")
    path = selection_file(tmp_path)
    config["inputs"]["phasing_reference_selection"] = str(path)
    inputs["artifacts"] = [dict(artifact_id="config_input_phasing_reference_selection", path=str(path), sha256=sha256_file(path))]
    assert stage_reference_selection(inputs, config, tmp_path, "synthetic_panel_v1").sample_ids == ("HG00001",)
    config["inputs"]["phasing_reference_selection"] = str(tmp_path / "different.json")
    with pytest.raises(ValueError, match="config_binding_invalid"):
        stage_reference_selection(inputs, config, tmp_path, "synthetic_panel_v1")


def test_extracts_selected_reference_and_records_distinct_method(tmp_path: Path, monkeypatch) -> None:
    from effet_fondateur.references import window
    monkeypatch.setattr(window, "_resolve_executable", lambda command: command)
    resolved, cached = _cached_reference(tmp_path)
    path = selection_file(tmp_path)
    selection = load_reference_selection(path, expected_sha256=sha256_file(path), panel_id=resolved.panel_id)
    runner = FakeBcftools(output_samples=("HG00001",))
    result = extract_reference_window(resolved, cached, tmp_path / "window", start_bp=1000, end_bp=2000,
                                      sample_selection=selection, command_runner=runner)
    document = json.loads(result.manifest_path.read_text())
    assert result.sample_count == 1
    assert document["sample_scope"] == "EXPLICIT_PARTITION_SELECTION"
    assert document["method_id"] == "bcftools_selected_samples_region_extract_v1"
    assert document["sample_selection"]["selection_sha256"] == sha256_file(path)
    assert not (result.output_dir / "reference.samples.txt").exists()
    assert any("--samples-file" in command for command in runner.commands)
    validate_json_document(document, "reference_window_manifest.schema.json")


def test_blocks_reference_that_returns_unselected_donors(tmp_path: Path, monkeypatch) -> None:
    from effet_fondateur.references import window
    monkeypatch.setattr(window, "_resolve_executable", lambda command: command)
    resolved, cached = _cached_reference(tmp_path)
    path = selection_file(tmp_path)
    selection = load_reference_selection(path, expected_sha256=sha256_file(path), panel_id=resolved.panel_id)
    with pytest.raises(ReferenceWindowIntegrityError, match="sample_set_or_order"):
        extract_reference_window(resolved, cached, tmp_path / "window", start_bp=1000, end_bp=2000,
                                 sample_selection=selection, command_runner=FakeBcftools())
    assert not (tmp_path / "window").exists()


@pytest.mark.parametrize("observed", [["HG00001", "HG00002"], [], ["HG00002"]])
def test_autosomal_extracted_membership_is_checked_not_just_cache_key(monkeypatch, observed: list[str]) -> None:
    from effet_fondateur.stages import phase_autosomal_panel as stage
    monkeypatch.setattr(stage, "_run", lambda command, timeout, code: subprocess.CompletedProcess(command, 0, "\n".join(observed), ""))
    with pytest.raises(stage.AutosomalPhasingBlockError, match="sample_set_or_order"):
        stage._validate_selected_reference("mock-bcftools", Path("synthetic.vcf"), ["HG00001"], 1)


@pytest.mark.parametrize("stage_name", ["phase_autosomal_panel", "phase_target_region"])
def test_both_stage_entrypoints_block_missing_required_selection_before_tools(tmp_path: Path, monkeypatch, stage_name: str) -> None:
    from importlib import import_module
    stage = import_module("effet_fondateur.stages." + stage_name)
    inputs = dict(parameters=dict(require_reference_sample_selection=True,
                                  reference_panel_id="synthetic_panel_v1", reference_cache_dir="unused"), artifacts=[])
    # L'enveloppe est testée par les contrats séparés ; ce test isole l'ordre
    # des gardes du point d'entrée, avant toute résolution d'outil/référence.
    monkeypatch.setattr(stage, "read_json", lambda path: inputs)
    monkeypatch.setattr(stage, "validate_json_document", lambda *args: None)
    monkeypatch.setattr(stage, "load_pipeline_config", lambda path: dict(inputs={}))
    output = tmp_path / "stages/12_synthetic"
    with pytest.raises(ValueError, match="selection_missing_or_ambiguous"):
        stage.execute(tmp_path / "inputs.json", output)
    assert not output.exists()


def test_stage_12_publishes_selected_reference_with_signed_input_and_distinct_universe(tmp_path: Path) -> None:
    from test_v2_phase_target_region import _prepare_phase_inputs
    from effet_fondateur.orchestrator.pipeline import run_pipeline
    config_path, runs = _prepare_phase_inputs(tmp_path)
    samples = tuple(f"HG{index:05d}" for index in range(1, 3203))
    path = selection_file(tmp_path, samples)
    document = json.loads(path.read_text())
    document["source_panel_id"] = "synthetic_stage_12_panel"
    path.write_text(json.dumps(document))
    config = yaml.safe_load(config_path.read_text())
    config["inputs"]["phasing_reference_selection"] = str(path)
    config["stages"]["phase_target_region"]["parameters"]["require_reference_sample_selection"] = True
    config_path.write_text(yaml.safe_dump(config, sort_keys=False))
    # Les exécutables sont les petites fixtures synthétiques du test de 12,
    # pas PLINK/bcftools/SHAPEIT réels. Le contrôle de rejet d'un donneur
    # surnuméraire dans une extraction est testé séparément ci-dessus.
    run = run_pipeline(config_path, runs)
    stage = run / "stages/12_phase_target_region"
    inputs = json.loads((stage / "stage_inputs.json").read_text())
    assert any(row["artifact_id"] == "config_input_phasing_reference_selection" and row["sha256"] == sha256_file(path) for row in inputs["artifacts"])
    window = json.loads((stage / "reference_window/reference_window_manifest.json").read_text())
    assert window["sample_scope"] == "EXPLICIT_PARTITION_SELECTION"
    assert window["sample_selection"]["selection_sha256"] == sha256_file(path)
    artifacts = json.loads((stage / "stage_outputs.json").read_text())["artifacts"]
    assert next(row["sample_set_id"] for row in artifacts if row["artifact_id"] == "reference_window_vcf").endswith("_selected_" + document["membership_sha256"][:12])


def test_autosomal_selection_summary_schema_is_self_contained_offline() -> None:
    artifact = dict(path="synthetic.bcf", sha256="a" * 64, size_bytes=1)
    document = dict(schema_version="1.0.0", method_id="shapeit5_reference_common_autosomes_v1", assembly="GRCh38",
                    sample_set_id="synthetic_samples", variant_set_id="synthetic_markers", sample_count=2, variant_count=44,
                    sample_selection=dict(selection_sha256="b" * 64, membership_sha256="c" * 64,
                                          partition_manifest_sha256="d" * 64, sample_count=2),
                    chromosomes=[dict(chromosome=chromosome, study_input_variants=2, reference_matched_variants=2,
                                      phased_variants=2, reference_cache_status="HIT", bcf=artifact, index=artifact,
                                      ibd_map=artifact, log=artifact) for chromosome in range(1,23)],
                    checks={name: "PASS" for name in ("same_samples_all_chromosomes", "complete_genotypes",
                                                     "all_genotypes_phased", "reference_alleles", "maps_monotonic", "all_22_autosomes")})
    validate_json_document(document, "autosomal_phasing_manifest.schema.json")
