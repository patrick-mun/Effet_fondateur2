"""Contrats, mutations et gardes D20 ; aucun PLINK réel dans ces tests."""

from copy import deepcopy
from pathlib import Path

import pytest

from effet_fondateur.calibration.evaluate_qc import aggregate_rules, compare_metrics, compare_preliminary_alerts, verify_final_input_binding
from effet_fondateur.calibration.materialize_ibd import _json
from effet_fondateur.calibration import qc_bench
from effet_fondateur.calibration.qc_bench import checked_inventory, file_inventory, prepare_qc_bench, write_json
from effet_fondateur.calibration.qc_oracle import confusion
from effet_fondateur.calibration.run_qc_bench import run_qc_bench
from effet_fondateur.calibration.technical_evaluation import TechnicalEvaluationError


CONFIG = Path(__file__).resolve().parents[2] / "calibration/configurations/d20_qc_bench_v1.json"


@pytest.fixture
def prepared(tmp_path, monkeypatch):
    workspace = tmp_path / "calibration_workspace"
    workspace.mkdir()
    original = qc_bench.generate_cases
    # Une fixture de composant suffit à tester publication et confinement ;
    # l'évaluateur de campagne exige séparément les quatorze cas réels.
    monkeypatch.setattr(qc_bench, "generate_cases", lambda config: original(config)[:1])
    prepare_qc_bench(configuration_path=CONFIG, workspace_root=workspace, run_id="inputs")
    return workspace, workspace / "inputs"


def test_public_private_separation_and_no_overwrite(prepared):
    workspace, panel = prepared
    document = checked_inventory(panel, "panel_manifest.json", "manifest_state", "FINALIZED")
    assert len(document["bundles"]) == 1
    source = panel / "blinded_inputs/bundles/B00"
    text = (source / "pipeline.yaml").read_text()
    assert "private_truth" not in text and "dock6_real" not in text
    assert "qc_preliminary:" in text and "qc_final:" not in text
    cohort_text = (panel / "blinded_inputs/cohort_fixtures/B00/cohorts.tsv").read_text()
    assert "D20_COHORT_CONTRACT_FIXTURE" in cohort_text
    assert "09_freeze" not in cohort_text
    with pytest.raises(TechnicalEvaluationError, match="run_exists"):
        prepare_qc_bench(configuration_path=CONFIG, workspace_root=workspace, run_id="inputs")


@pytest.mark.parametrize("mutation", ["changed", "extra", "symlink", "traversal", "wrong_state", "duplicate_inventory"])
def test_inventory_rejects_tampering(prepared, mutation):
    _, panel = prepared
    manifest = _json(panel / "panel_manifest.json")
    path = panel / "blinded_inputs/qc_execution.json"
    if mutation == "changed":
        path.write_text("{}\n")
    elif mutation == "extra":
        (panel / "extra.txt").write_text("unlisted\n")
    elif mutation == "symlink":
        path.unlink()
        path.symlink_to(panel / "private_truth/configuration.json")
    elif mutation == "traversal":
        manifest["artifacts"][0]["path"] = "../outside.json"
    elif mutation == "wrong_state":
        manifest["manifest_state"] = "OPEN"
    else:
        manifest["artifacts"].append(manifest["artifacts"][0])
    if mutation in {"traversal", "wrong_state", "duplicate_inventory"}:
        write_json(panel / "panel_manifest.json", manifest)
    with pytest.raises(TechnicalEvaluationError):
        checked_inventory(panel, "panel_manifest.json", "manifest_state", "FINALIZED")


def test_runner_rejects_unapproved_stage_scope_before_external_call(prepared, monkeypatch):
    import yaml
    workspace, panel = prepared
    path = panel / "blinded_inputs/bundles/B00/pipeline.yaml"
    config = yaml.safe_load(path.read_text())
    config["stages"]["qc_final"] = dict(enabled=True, parameters={})
    path.write_text(yaml.safe_dump(config))
    manifest = _json(panel / "panel_manifest.json")
    manifest["artifacts"] = file_inventory(panel, ("panel_manifest.json",))
    write_json(panel / "panel_manifest.json", manifest)
    monkeypatch.setattr("effet_fondateur.calibration.run_qc_bench.subprocess.run", lambda *args, **kwargs: pytest.fail("external process forbidden"))
    with pytest.raises(TechnicalEvaluationError, match="scope_forbidden"):
        run_qc_bench(panel_root=panel, workspace_root=workspace, run_id="attempt")
    assert not (workspace / "attempt/run_index.json").exists()
    assert (workspace / "attempt.failed").is_dir()


def test_runner_rejects_missing_closure(prepared):
    workspace, panel = prepared
    (panel / "panel_manifest.json").unlink()
    with pytest.raises(TechnicalEvaluationError, match="closure_document_missing"):
        run_qc_bench(panel_root=panel, workspace_root=workspace, run_id="attempt")
    assert not (workspace / "attempt").exists()


@pytest.mark.parametrize("mutation,metric", [("delete", "lost"), ("duplicate", "duplicate"), ("add", "extra"), ("swap", "order"),
                                             ("status", "decision"), ("reason", "decision"), ("count", "decision"), ("numeric", "numeric"), ("nan", "numeric")])
def test_metric_comparator_catches_semantic_mutations(mutation, metric):
    expected = [dict(ID="A", STATUS="PASS", REASONS="", N=2, VALUE=.25), dict(ID="B", STATUS="EXCLUDED", REASONS="missingness", N=1, VALUE=.5)]
    actual = [{key: str(value) for key, value in row.items()} for row in expected]
    if mutation == "delete":
        actual.pop()
    elif mutation == "duplicate":
        actual.append(deepcopy(actual[0]))
    elif mutation == "add":
        actual.append(dict(ID="C", STATUS="PASS", REASONS="", N="0", VALUE="0"))
    elif mutation == "swap":
        actual.reverse()
    else:
        field, value = {"status": ("STATUS", "EXCLUDED"), "reason": ("REASONS", "arbitrary"), "count": ("N", "2.0"),
                        "numeric": ("VALUE", "0.4"), "nan": ("VALUE", "nan")}[mutation]
        actual[0][field] = value
    assert compare_metrics(expected, actual, "ID")[metric] > 0


def test_numeric_rounding_tolerance_is_not_decision_tolerance():
    expected = [dict(ID="A", VALUE=.33333333333, COUNT=1)]
    assert compare_metrics(expected, [dict(ID="A", VALUE=".3333", COUNT="1")], "ID")["numeric"] == 0
    assert compare_metrics(expected, [dict(ID="A", VALUE=".3333", COUNT="1.0000")], "ID")["decision"] == 1


def test_tiny_hwe_p_values_are_not_approximated_to_zero_or_across_threshold():
    expected = [dict(ID="A", HWE_P=1.234567e-28)]
    assert compare_metrics(expected, [dict(ID="A", HWE_P="1.235e-28")], "ID")["numeric"] == 0
    for value in ("0", "4e-6", "1e-20"):
        assert compare_metrics(expected, [dict(ID="A", HWE_P=value)], "ID")["numeric"] == 1


def test_confusion_aggregation_uses_counts_and_detects_empty_classes():
    rows = [dict(rules=dict(rule=dict(unit="variant", **confusion({"A"}, {"A"}, {"A", "B"})))),
            dict(rules=dict(rule=dict(unit="variant", **confusion({"A"}, set(), {"A", "B", "C"}))))]
    aggregate = aggregate_rules(rows, "05")["rule"]
    assert aggregate["sensitivity"] == .5
    assert aggregate["specificity"] == 1.0
    assert aggregate["coverage"] == "POSITIVE_AND_NEGATIVE"
    assert aggregate_rules([dict(rules=dict(rule=dict(unit="variant", **confusion(set(), set(), {"A"}))))], "05")["rule"]["coverage"] == "MISSING_CLASS"


def test_final_input_binding_is_order_independent_and_rejects_changes():
    parameters = dict(sample_missingness_max=.1)
    first = dict(artifact_id="samples_master", sha256="a" * 64)
    second = dict(artifact_id="genomewide_pre_qc_bed", sha256="b" * 64)
    inputs = dict(artifacts=[first, second], parameters=parameters, signature="c" * 64)
    audit = dict(inputs=[second, first], signature="c" * 64)
    verify_final_input_binding(inputs, audit, parameters)
    audit["inputs"] = [dict(second, sha256="d" * 64), first]
    with pytest.raises(TechnicalEvaluationError, match="binding_invalid"):
        verify_final_input_binding(inputs, audit, parameters)


def test_only_known_unused_v1_reports_are_allowed():
    first = dict(artifact_id="samples_master", sha256="a" * 64)
    inputs = dict(artifacts=[first, dict(artifact_id="target_variant_injection_report", sha256="b" * 64)], parameters={}, signature="c" * 64)
    audit = dict(inputs=[first], signature="c" * 64)
    verify_final_input_binding(inputs, audit, {})
    inputs["artifacts"][-1]["artifact_id"] = "unknown_input"
    with pytest.raises(TechnicalEvaluationError, match="binding_invalid"):
        verify_final_input_binding(inputs, audit, {})


@pytest.mark.parametrize("mutation", ["duplicate", "missing_consumed", "signature", "parameters"])
def test_final_input_binding_rejects_ambiguous_or_unbound_contract(mutation):
    first = dict(artifact_id="samples_master", sha256="a" * 64)
    inputs = dict(artifacts=[first], parameters={}, signature="c" * 64)
    audit = dict(inputs=[first], signature="c" * 64)
    if mutation == "duplicate":
        inputs["artifacts"].append(first)
    elif mutation == "missing_consumed":
        inputs["artifacts"] = []
    elif mutation == "signature":
        audit["signature"] = "d" * 64
    else:
        inputs["parameters"] = dict(sample_missingness_max=.9)
    with pytest.raises(TechnicalEvaluationError, match="binding_invalid"):
        verify_final_input_binding(inputs, audit, {})


@pytest.mark.parametrize("field,value", [("THRESHOLD", "0.2"), ("OBSERVED_VALUE", "0.5"), ("STATUS", "PASS"), ("SCOPE", "PAIR"), ("REASON_CODE", "unknown")])
def test_alert_comparator_checks_values_thresholds_scope_and_reasons(field, value):
    from effet_fondateur.calibration.qc_oracle import QcExpectation
    expected = QcExpectation([dict(IID="A", EXCLUSION_CODE="sample_missingness", QC_STATUS="EXCLUDED", F_MISS=.3)], [], [], None)
    actual = [dict(SCOPE="SAMPLE", ENTITY_ID="A", METRIC="F_MISS", OBSERVED_VALUE=".3", THRESHOLD=".1", STATUS="EXCLUDED", REASON_CODE="sample_missingness"),
              dict(SCOPE="GLOBAL", ENTITY_ID="genomewide_pre_qc", METRIC="SEX_DISCORDANCE", OBSERVED_VALUE="", THRESHOLD="", STATUS="NOT_EVALUATED", REASON_CODE="autosomal_dataset_without_sex_chromosomes"),
              dict(SCOPE="GLOBAL", ENTITY_ID="genomewide_pre_qc", METRIC="HWE_FILTER", OBSERVED_VALUE="not_applied", THRESHOLD="", STATUS="INFO", REASON_CODE="cohorts_not_frozen")]
    parameters = _json(CONFIG)["preliminary_parameters"]
    baseline = compare_preliminary_alerts(expected, actual, parameters)
    assert not any(count for key, count in baseline.items() if key != "compared")
    actual[0][field] = value
    comparison = compare_preliminary_alerts(expected, actual, parameters)
    assert any(count for key, count in comparison.items() if key != "compared")
