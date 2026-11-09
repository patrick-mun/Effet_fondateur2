"""D22 : préparation confinée, inventaire et refus avant processus externe."""

from pathlib import Path

import pytest
import yaml

from effet_fondateur.calibration.materialize_ibd import _json
from effet_fondateur.calibration.qc_bench import file_inventory, write_json
from effet_fondateur.calibration.run_structure_bench import run_structure_bench
from effet_fondateur.calibration.structure_bench import checked_document, prepare_structure_bench


@pytest.fixture
def prepared(tmp_path):
    repository = Path(__file__).resolve().parents[2]
    config = _json(repository / "calibration/configurations/d22_structure_bench_v1.json")
    config["markers_per_autosome"] = 32
    config_path = tmp_path / "config.json"
    write_json(config_path, config)
    (tmp_path / config["target_configuration"]).write_bytes((repository / "calibration/configurations" / config["target_configuration"]).read_bytes())
    workspace = tmp_path / "calibration_workspace"
    workspace.mkdir()
    document = prepare_structure_bench(configuration_path=config_path, workspace_root=workspace, run_id="inputs")
    return workspace, workspace / "inputs", document, config_path


def test_closed_preparation_refuses_overwrite_and_changed_truth(prepared):
    workspace, root, document, configuration = prepared
    assert checked_document(root, "panel_manifest.json", "manifest_state", "FINALIZED") == document
    assert len(document["bundles"]) == 4
    for path in (root / "blinded_inputs").rglob("*.yaml"):
        assert "private_truth" not in path.read_text()
    with pytest.raises(ValueError, match="run_exists"):
        prepare_structure_bench(configuration_path=configuration, workspace_root=workspace, run_id="inputs")
    (root / "private_truth/structure_truth.json").write_text("{}")
    with pytest.raises(ValueError, match="artifact_changed"):
        checked_document(root, "panel_manifest.json", "manifest_state", "FINALIZED")


@pytest.mark.parametrize("mutation", ["approval", "scope", "external_tool", "review_role"])
def test_scope_or_review_changes_fail_before_any_external_call(prepared, monkeypatch, mutation):
    workspace, root, document, _ = prepared
    pipeline = root / "blinded_inputs/bundles/B00/pipeline.yaml"
    config = yaml.safe_load(pipeline.read_text())
    if mutation == "approval":
        config["stages"]["build_sample_registry"]["parameters"]["manual_approval"]["approved"] = False
    elif mutation == "scope":
        config["stages"]["analyze_population_structure"] = dict(enabled=True, parameters={})
    elif mutation == "external_tool":
        config["tools"]["king"] = "king"
    else:
        path = root / "blinded_inputs/execution.json"
        execution = _json(path)
        execution["fixture_review_policy"]["reviewer_role"] = "scientific_reviewer"
        write_json(path, execution)
    pipeline.write_text(yaml.safe_dump(config))
    document["artifacts"] = file_inventory(root, ("panel_manifest.json",))
    write_json(root / "panel_manifest.json", document)

    def forbidden_process(*args, **kwargs):
        pytest.fail("Le garde D22 doit refuser avant tout processus externe")

    monkeypatch.setattr("effet_fondateur.calibration.run_structure_bench.subprocess.run", forbidden_process)
    with pytest.raises(ValueError, match="scope"):
        run_structure_bench(panel_root=root, workspace_root=workspace, run_id="forbidden")
