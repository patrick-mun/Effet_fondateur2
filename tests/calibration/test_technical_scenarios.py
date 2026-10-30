"""D19 03–04 : préparation isolée et contrats, sans outil génétique réel."""

import json
from pathlib import Path
import shutil
import sys
from types import SimpleNamespace

import pytest

from effet_fondateur.calibration.evaluate_technical import evaluate_technical, load_checked_evaluation, validate_schema
from effet_fondateur.calibration.export_blinded_inputs import export_blinded_inputs, load_checked_export
from effet_fondateur.calibration.materialize_ibd import _sha
from effet_fondateur.calibration.run_technical_scenarios import run_technical_scenarios
from effet_fondateur.calibration.report_technical_evaluation import build_reports
from effet_fondateur.calibration.technical_evaluation import TechnicalEvaluationError
from effet_fondateur.calibration.technical_rejection import evaluate_target_rejection
from effet_fondateur.calibration.technical_scenarios import (
    MANIFEST_PATH, SCENARIOS, load_checked_panel, prepare_technical_panel,
    target_expectation, technical_configuration, validate_derivation, validate_scenario_config,
)
from test_blinded_export import _fixture as export_fixture
from test_technical_evaluation import _matrix, write_plink


ROOT = Path(__file__).resolve().parents[2]
CREATED_AT = "2026-10-30T00:00:00Z"


def config():
    return json.loads((ROOT / "calibration/configurations/d19_technical_target_v1.json").read_text())


@pytest.fixture(scope="module")
def prepared(tmp_path_factory):
    arguments = export_fixture(tmp_path_factory.mktemp("technical_panels"))
    export_blinded_inputs(**arguments)
    workspace = arguments["workspace_root"]
    d18 = workspace / "d18_test"
    config_path = ROOT / "calibration/configurations/d19_technical_target_v1.json"
    for scenario, run_id in zip(SCENARIOS, ("positive", "negative"), strict=True):
        prepare_technical_panel(d18_root=d18, configuration_path=config_path, workspace_root=workspace,
                                run_id=run_id, created_at=CREATED_AT, scenario=scenario)
    return workspace, d18


def test_target_rule_is_explicit_and_independent_of_clinical_labels():
    rows = [dict(SAMPLE_ID=f"S{index}", GROUP_LABEL="ANY", CLINICAL_STATUS="UNKNOWN") for index in range(7)]
    first = target_expectation(rows, config(), "B" + "a" * 16)
    for row in rows:
        row.update(GROUP_LABEL="OTHER", CLINICAL_STATUS="AFFECTED")
    assert first == target_expectation(rows, config(), "B" + "a" * 16)
    assert [row["genotype"] for row in first["genotypes"]] == ["T/T", "T/A", "A/A", "T/T", "T/A", "A/A", "T/T"]


@pytest.mark.parametrize("defect", ["extra", "approval", "alleles", "source", "rule"])
def test_scenario_contract_refuses_unsafe_or_incomplete_parameters(defect):
    document = config()
    if defect == "extra": document["unknown"] = True
    elif defect == "approval": document["manual_approval"]["approved"] = False
    elif defect == "alleles": document["target_metadata"]["alt"] = document["target_metadata"]["ref"]
    elif defect == "source": document["target_metadata"]["annotation_validation"]["source_id"] = "INVENTED"
    else: document["genotype_rule"] = "FROM_CLINICAL_STATUS"
    with pytest.raises(ValueError): validate_scenario_config(document)


def test_configuration_is_bounded_and_retains_supplied_approval():
    document = config()
    result = technical_configuration(document)
    assert len(result["stages"]) == 5
    assert result["stages"]["build_sample_registry"]["parameters"]["manual_approval"] == document["manual_approval"]
    assert all(value is None for key, value in result["tools"].items() if key != "plink")
    assert result["policies"]["allow_clinical_genotype_inference"] is False


@pytest.mark.parametrize("scenario", ["positive", "negative"])
def test_prepared_inputs_are_sealed_derived_and_private(prepared, scenario):
    workspace, d18 = prepared
    panel = load_checked_panel(workspace / scenario)
    validate_derivation(workspace / scenario, d18)
    assert len(panel["inventory"]["bundles"]) == (8 if scenario == "positive" else 1)
    assert panel["calibrated"] is False
    assert all(path.stat().st_mode & 0o077 == 0 for path in (workspace / scenario).rglob("*"))
    load_checked_export(d18)
    with pytest.raises(TechnicalEvaluationError, match="run_exists"):
        prepare_technical_panel(d18_root=d18, configuration_path=ROOT / "calibration/configurations/d19_technical_target_v1.json",
                                workspace_root=workspace, run_id=scenario, created_at=CREATED_AT)


@pytest.mark.parametrize("defect", ["unlisted", "modified", "symlink"])
def test_panel_reader_rejects_tree_or_checksum_changes(prepared, tmp_path, defect):
    workspace, _ = prepared
    panel = tmp_path / "panel"
    shutil.copytree(workspace / "negative", panel)
    if defect == "unlisted": (panel / "extra").mkdir()
    else:
        path = next((panel / "blinded_inputs").glob("bundles/*/target_genotypes.tsv"))
        if defect == "modified": path.write_text("wrong")
        else:
            target = tmp_path / "other"
            path.rename(target)
            path.symlink_to(target)
    with pytest.raises(ValueError): load_checked_panel(panel)


@pytest.mark.parametrize("defect", ["acpa", "target_input", "target_truth"])
def test_derivation_catches_substitutions_even_if_panel_checksums_are_resealed(prepared, tmp_path, defect):
    workspace, d18 = prepared
    panel = tmp_path / "panel"
    shutil.copytree(workspace / "positive", panel)
    if defect == "acpa":
        path = next((panel / "blinded_inputs").glob("bundles/*/acpa_samples/*.txt"))
        path.write_text(path.read_text() + "injected\tAA\t\t5\t999999\n")
    elif defect == "target_input":
        path = next((panel / "blinded_inputs").glob("bundles/*/target_genotypes.tsv"))
        path.write_text(path.read_text().replace("T/T", "A/A", 1))
    else:
        path = next((panel / "private_truth/d19_inputs").glob("*.target.json"))
        document = json.loads(path.read_text())
        document["genotypes"][0]["genotype"] = "A/A"
        path.write_text(json.dumps(document))
    manifest_path = panel / MANIFEST_PATH
    manifest = json.loads(manifest_path.read_text())
    for artifact in manifest["artifacts"]:
        artifact["sha256"] = _sha(panel / artifact["path"])
    for artifact in manifest["inventory"]["artifacts"]:
        actual = panel / "blinded_inputs" / artifact["path"]
        artifact.update(sha256=_sha(actual), size_bytes=actual.stat().st_size)
    manifest_path.write_text(json.dumps(manifest))
    load_checked_panel(panel)
    with pytest.raises(TechnicalEvaluationError, match="derivation"):
        validate_derivation(panel, d18)


def _mock_pipeline_run(command, *, cwd, stdout, stderr, timeout, check):
    """Simule la clôture de l'orchestrateur ; aucun binaire externe n'est lancé."""
    output = Path(command[-1]) / "synthetic_run"
    output.mkdir()
    negative = "negative" in str(cwd)
    stages = []
    for stage, name in zip(("00", "01", "02", "03", "04"), ("initialize_run", "validate_sources", "build_sample_registry", "convert_acpa", "prepare_target_variant_dataset"), strict=True):
        stages.append(dict(stage_id=stage, stage_name=name, state="FAILED" if negative and stage == "04" else "SUCCEEDED", critical=True,
                           signature="a" * 64, started_at=CREATED_AT, completed_at=CREATED_AT, duration_seconds=0,
                           audit_path=None, audit_sha256=None, stage_outputs_sha256=None, attempt_count=1, last_error_code=4 if negative and stage == "04" else None))
    manifest = dict(schema_version="1.0.0", run_id="synthetic_run", project_id="technical_03_04_only", created_at=CREATED_AT,
                    updated_at=CREATED_AT, global_status="BLOCKED" if negative else "TECHNICALLY_VALID", config_path="config.resolved.yaml",
                    config_sha256="b" * 64, environment_path="environment.json", stages=stages, manual_decisions_required=[])
    (output / "manifest.json").write_text(json.dumps(manifest))
    (output / "environment.json").write_text("{}")
    if negative:
        attempt = output / "stages/attempts/04_prepare_target_variant_dataset.abc.failed"
        attempt.mkdir(parents=True)
        (attempt / "stderr.log").write_text("target_variant_already_present\n")
    return SimpleNamespace(returncode=2 if negative else 0)


@pytest.mark.parametrize("scenario", ["positive", "negative"])
def test_runner_publishes_only_closed_bounded_receipts(prepared, monkeypatch, scenario):
    workspace, _ = prepared
    monkeypatch.setattr("effet_fondateur.calibration.run_technical_scenarios.subprocess.run", _mock_pipeline_run)
    index = run_technical_scenarios(panel_root=workspace / scenario, workspace_root=workspace, run_id=f"mock_{scenario}")
    validate_schema(index, "technical_run_index")
    assert index["closure_state"] == "CLOSED"
    assert index["technical_panel_manifest_sha256"] == _sha(workspace / scenario / MANIFEST_PATH)
    assert len(index["runs"]) == (8 if scenario == "positive" else 1)
    assert all(("rejection_evidence" in entry) == (scenario == "negative") for entry in index["runs"])
    assert all("private_truth" in entry["target_truth_path"] for entry in index["runs"])


def test_unexpected_runner_failure_preserves_diagnostics_without_closed_receipt(prepared, monkeypatch):
    workspace, _ = prepared
    def unexpected(*args, **kwargs):
        _mock_pipeline_run(*args, **kwargs)
        return SimpleNamespace(returncode=9)
    monkeypatch.setattr("effet_fondateur.calibration.run_technical_scenarios.subprocess.run", unexpected)
    with pytest.raises(TechnicalEvaluationError, match="process_or_manifest"):
        run_technical_scenarios(panel_root=workspace / "positive", workspace_root=workspace, run_id="unexpected")
    assert not (workspace / "unexpected").exists()
    assert (workspace / "unexpected.failed").is_dir()
    assert not (workspace / "unexpected.failed/run_index.json").exists()
    with pytest.raises(TechnicalEvaluationError, match="diagnostics_exist"):
        run_technical_scenarios(panel_root=workspace / "positive", workspace_root=workspace, run_id="unexpected")


@pytest.mark.parametrize("relative", ["data/input/calibration_workspace", "dock6_real/calibration_workspace"])
def test_preparation_refuses_real_data_roots_before_loading(tmp_path, relative):
    with pytest.raises(TechnicalEvaluationError, match="workspace"):
        prepare_technical_panel(d18_root=tmp_path / "missing", configuration_path=tmp_path / "missing", workspace_root=tmp_path / relative,
                                run_id="bad", created_at=CREATED_AT)


def _rejection_fixture(tmp_path, monkeypatch):
    matrix = _matrix()
    prefix = tmp_path / "base"
    write_plink(prefix, matrix)
    metadata = [dict(FID=fid, IID=iid, PID="0", MID="0", SEX="UNKNOWN", CLINICAL_STATUS="UNKNOWN") for fid, iid in matrix.samples]
    attempt = tmp_path / "stages/attempts/04_prepare_target_variant_dataset.abc.failed"
    attempt.mkdir(parents=True)
    path = attempt / "stderr.log"
    path.write_text("target_variant_already_present\n")
    inputs = dict(schema_version="1.0.0", run_id="run", stage_id="04", stage_name="prepare_target_variant_dataset", signature="a" * 64,
                  attempt_number=1, published_output_dir="stages/04_prepare_target_variant_dataset", parameters={}, artifacts=[])
    (attempt / "stage_inputs.json").write_text(json.dumps(inputs))
    manifest = dict(run_id="run", stages=[dict(stage_id="04", stage_name=inputs["stage_name"], state="FAILED", last_error_code=4, attempt_count=1, signature="a" * 64)])
    receipt = dict(rejection_evidence=dict(path=path.relative_to(tmp_path).as_posix(), sha256=_sha(path)))
    monkeypatch.setattr("effet_fondateur.calibration.technical_rejection._stage_artifacts", lambda *_: dict(target_chromosome_base_fam=prefix.with_suffix(".fam")))
    monkeypatch.setattr("effet_fondateur.calibration.technical_rejection._dataset", lambda *_: matrix)
    return dict(run_root=tmp_path, manifest=manifest, receipt=receipt, expected=matrix, metadata=metadata,
                target=dict(chromosome=5, position_bp=100, project_variant_id="chip1"))


def test_duplicate_control_accepts_exact_refusal_but_never_claims_insertion(tmp_path, monkeypatch):
    arguments = _rejection_fixture(tmp_path, monkeypatch)
    verdict = evaluate_target_rejection(**arguments)
    assert verdict["status"] == "PASS"
    assert verdict["metrics"]["expected_duplicate_target_rejection_count"] == 1
    assert verdict["metrics"]["target_insertion_performed_count"] == 0
    assert verdict["metrics"]["target_genotype_discordance_rate"] is None


@pytest.mark.parametrize("defect", ["reason", "checksum", "state", "signature", "locus", "published"])
def test_duplicate_control_never_accepts_unrelated_failure(tmp_path, monkeypatch, defect):
    arguments = _rejection_fixture(tmp_path, monkeypatch)
    evidence = arguments["receipt"]["rejection_evidence"]
    path = tmp_path / evidence["path"]
    if defect == "reason":
        path.write_text("unrelated_failure\n")
        evidence["sha256"] = _sha(path)
    elif defect == "checksum": evidence["sha256"] = "0" * 64
    elif defect == "state": arguments["manifest"]["stages"][0]["state"] = "SUCCEEDED"
    elif defect == "signature": arguments["manifest"]["stages"][0]["signature"] = "b" * 64
    elif defect == "locus": arguments["target"]["position_bp"] = 999
    else: (tmp_path / "stages/04_prepare_target_variant_dataset").mkdir()
    try:
        assert evaluate_target_rejection(**arguments)["status"] == "FAIL"
    except TechnicalEvaluationError:
        assert defect in {"checksum", "signature"}


@pytest.mark.parametrize("scenario", SCENARIOS)
def test_full_python_pipeline_with_simulated_plink_is_evaluable_00_04(prepared, tmp_path, scenario):
    workspace, d18 = prepared
    executable = tmp_path / "fake_plink"
    source = ROOT / "tests/calibration/fixtures/fake_plink_technical.py"
    executable.write_text(f"#!{sys.executable}\n" + source.read_text().split("\n", 1)[1])
    executable.chmod(0o700)
    document = config()
    document["plink_executable"] = str(executable)
    configuration = tmp_path / "config.json"
    configuration.write_text(json.dumps(document))
    label = "fake_positive" if scenario == SCENARIOS[0] else "fake_negative"
    prepare_technical_panel(d18_root=d18, configuration_path=configuration, workspace_root=workspace,
                            run_id=label, created_at=CREATED_AT, scenario=scenario)
    index = run_technical_scenarios(panel_root=workspace / label, workspace_root=workspace, run_id=label + "_runs")
    manifest = evaluate_technical(d18_root=d18, technical_panel_root=workspace / label,
                                  run_index_path=workspace / (label + "_runs") / "run_index.json", workspace_root=workspace,
                                  run_id=label + "_evaluation", created_at=CREATED_AT)
    assert all(item["status"] == "PASS" for item in manifest["stages"]), manifest["stages"]
    assert manifest["scenario"] == scenario
    assert manifest["bundle_count"] == (8 if scenario == SCENARIOS[0] else 1)
    assert load_checked_evaluation(workspace / (label + "_evaluation")) == manifest
    template = tmp_path / "template.html"
    template.write_text("<!doctype html><html><head><title>Template</title></head><body></body></html>")
    paths = build_reports(run_root=workspace / (label + "_evaluation"), template_path=template, output_root=workspace / (label + "_reports"))
    html = paths[-1].read_text()
    assert ("valide le refus de doublon" in html) == (scenario == SCENARIOS[1])
    assert "mesure expérimentale" in html and "elle n’est pas disponible" not in html
    if scenario == SCENARIOS[0]:
        evaluations = json.loads((workspace / (label + "_evaluation") / "private_truth/d19/evaluation_04.json").read_text())["evaluations"]
        assert all(item["metrics"]["target_genotype_discordance_rate"] == 0 for item in evaluations)
    with pytest.raises(TechnicalEvaluationError, match="technical_panel_required"):
        evaluate_technical(d18_root=d18, run_index_path=workspace / (label + "_runs") / "run_index.json", workspace_root=workspace,
                           run_id=label + "_wrong_evaluation", created_at=CREATED_AT)
