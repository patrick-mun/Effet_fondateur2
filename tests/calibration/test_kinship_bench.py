"""D21 : confinement, empreintes, oracles LD et défauts réellement détectés."""

from copy import deepcopy
from pathlib import Path

import numpy as np
import pytest
import yaml

from effet_fondateur.calibration.kinship_bench import checked_document, prepare_kinship_bench
from effet_fondateur.calibration.kinship_panel_oracle import allele_frequencies, dosage_r2, founder_mask, residual_ld
from effet_fondateur.calibration.kinship_stage_oracle import check_pairs, evaluate_kinship_stage
from effet_fondateur.calibration.materialize_ibd import _json
from effet_fondateur.calibration.qc_bench import file_inventory, write_json, write_table
from effet_fondateur.calibration.run_kinship_bench import run_kinship_bench
from effet_fondateur.calibration.simulation.kinship_bench import DEGREES
from effet_fondateur.calibration.technical_evaluation import GenotypeMatrix


CONFIG = Path(__file__).resolve().parents[2] / "calibration/configurations/d21_kinship_bench_v1.json"


def test_preparation_is_neutral_confined_and_immutable(tmp_path):
    config = _json(CONFIG)
    config.update(markers_per_autosome=32, families_per_relation_per_stratum=1, bootstrap_replicates=100)
    source = tmp_path / "configuration.json"
    write_json(source, config)
    workspace = tmp_path / "calibration_workspace"
    workspace.mkdir()
    document = prepare_kinship_bench(configuration_path=source, workspace_root=workspace, run_id="inputs")
    root = workspace / "inputs"
    assert checked_document(root, "panel_manifest.json", "manifest_state", "FINALIZED") == document
    assert [entry["sample_count"] for entry in document["bundles"]] == [28] * 3
    for path in (root / "blinded_inputs").rglob("*.yaml"):
        text = path.read_text()
        assert "private_truth" not in text and "ARTIFICIAL_A" not in text
    with pytest.raises(ValueError, match="run_exists"):
        prepare_kinship_bench(configuration_path=source, workspace_root=workspace, run_id="inputs")
    (root / "private_truth/pedigrees.json").write_text("{}")
    with pytest.raises(ValueError, match="artifact_changed"):
        checked_document(root, "panel_manifest.json", "manifest_state", "FINALIZED")


def test_frequency_and_ld_use_declared_founders_and_joint_calls():
    matrix = GenotypeMatrix((("F", "A"), ("F", "B"), ("F", "C")), ((1, 1, "v1"), (1, 2, "v2")),
                            np.array([[[0, 0], [0, 2], [2, 2]], [[0, 0], [0, 2], [2, 2]]], dtype=np.int8))
    registry = [dict(FID="F", IID=sample, PID="A" if sample == "C" else "0", MID="0") for sample in ("A", "B", "C")]
    founders = founder_mask(matrix, registry)
    assert founders.tolist() == [True, True, False]
    assert allele_frequencies(matrix, founders).tolist() == [.25, .25]
    first = np.array([[0, 1, 2, -1], [0, 0, 0, 0], [0, 1, 2, 0]])
    second = np.array([[2, 1, 0, 1], [0, 1, 2, 1], [0, 1, 2, -1]])
    result = dosage_r2(first, second)
    assert result[0] == 1 and np.isnan(result[1]) and result[2] == 1
    parameters = dict(residual_ld_window_variants=2, residual_ld_window_kb=1)
    assert residual_ld(matrix, founders, parameters).tolist() == [1]
    parameters["residual_ld_window_variants"] = 1
    assert residual_ld(matrix, founders, parameters).size == 0


@pytest.mark.parametrize("mutation", ["approval", "stage_scope", "external_tool"])
def test_runner_refuses_changed_scope_before_external_process(tmp_path, monkeypatch, mutation):
    config = _json(CONFIG)
    config.update(markers_per_autosome=32, families_per_relation_per_stratum=1, bootstrap_replicates=100)
    source = tmp_path / "config.json"
    write_json(source, config)
    workspace = tmp_path / "calibration_workspace"
    workspace.mkdir()
    document = prepare_kinship_bench(configuration_path=source, workspace_root=workspace, run_id="inputs")
    root = workspace / "inputs"
    pipeline = root / "blinded_inputs/bundles/B00/pipeline.yaml"
    public = yaml.safe_load(pipeline.read_text())
    if mutation == "approval":
        public["stages"]["build_sample_registry"]["parameters"]["manual_approval"]["approved"] = False
    elif mutation == "stage_scope":
        public["stages"]["infer_kinship"]["enabled"] = False
    else:
        public["tools"]["bcftools"] = "bcftools"
    pipeline.write_text(yaml.safe_dump(public, sort_keys=False))
    document["artifacts"] = file_inventory(root, ("panel_manifest.json",))
    write_json(root / "panel_manifest.json", document)

    def forbidden_process(*args, **kwargs):
        pytest.fail("Le garde doit refuser avant tout processus externe")

    monkeypatch.setattr("effet_fondateur.calibration.run_kinship_bench.subprocess.run", forbidden_process)
    with pytest.raises(ValueError, match="runner_scope_forbidden"):
        run_kinship_bench(panel_root=root, workspace_root=workspace, run_id="forbidden")
    assert (workspace / "forbidden.failed").is_dir()


@pytest.fixture
def unrelated_pair():
    calls = np.array([[[0, 2], [0, 0]], [[0, 0], [0, 2]]] * 2, dtype=np.int8)
    matrix = GenotypeMatrix((("F", "A"), ("F", "B")), tuple((1, index, str(index)) for index in range(4)), calls)
    registry = [dict(SAMPLE_ID=sample, FID="F", IID=sample, PID="0", MID="0") for sample in ("A", "B")]
    pair = dict(PAIR_ID="kinship_pair_000001", SAMPLE_ID_1="A", SAMPLE_ID_2="B", FID_1="F", IID_1="A", FID_2="F", IID_2="B",
                N_SNP="4", HETHET="0", IBS0="0", KINSHIP="0", INFERRED_DEGREE="UNRELATED", INFERRED_RELATION="UNRELATED",
                DECLARED_RELATION="SAME_FAMILY_UNSPECIFIED", PAIR_ORIGIN="WITHIN_FAMILY", CONCORDANCE_STATUS="NO_DECLARED_RELATION")
    return matrix, registry, pair


@pytest.mark.parametrize("mutation", ["numeric", "classification", "duplicate", "noncanonical"])
def test_pair_oracle_detects_changed_outputs(unrelated_pair, mutation):
    matrix, registry, pair = unrelated_pair
    config = _json(CONFIG)
    good, _, _ = check_pairs(matrix, [pair], registry, config["kinship_parameters"], config["king_absolute_tolerance"])
    assert all(value == 0 for key, value in good.items() if key != "compared")
    rows = [deepcopy(pair)]
    if mutation == "numeric":
        rows[0]["KINSHIP"] = ".01"
    elif mutation == "classification":
        rows[0]["INFERRED_DEGREE"] = "FOURTH"
    elif mutation == "duplicate":
        rows.append(deepcopy(pair))
    else:
        rows[0].update(SAMPLE_ID_1="B", SAMPLE_ID_2="A")
    result, _, _ = check_pairs(matrix, rows, registry, config["kinship_parameters"], config["king_absolute_tolerance"])
    assert any(value for key, value in result.items() if key != "compared")


def test_conformant_code_can_fail_true_independence(unrelated_pair, tmp_path):
    matrix, registry, pair = unrelated_pair
    config = _json(CONFIG)
    config["bootstrap_replicates"] = 100
    qc = [dict(SAMPLE_ID=sample, QC_STATUS="PASS", F_MISS="0", HETEROZYGOSITY_Z="0") for sample in ("A", "B")]
    proposal = [dict(SAMPLE_ID=sample, FID="F", IID=sample, COMPONENT_ID=f"kinship_component_{index:06d}", PROPOSED_INCLUDED="true",
                     PROPOSAL_STATUS="RETAIN", REPRESENTATIVE_SAMPLE_ID="", RELATION_BASIS="", F_MISS="0", ABS_HETEROZYGOSITY_Z="0",
                     DECISION_SOURCE="stage_07_quality_aware_greedy_maximal", MANUAL_VALIDATION_REQUIRED="false") for index, sample in enumerate(("A", "B"), 1)]
    summary = [dict(INFERRED_DEGREE=degree, PAIR_COUNT=int(degree == "UNRELATED"), DECLARED_PAIR_COUNT=int(degree == "UNRELATED"),
                    CRYPTIC_PAIR_COUNT=0, CONCORDANT_PAIR_COUNT=0, DISCORDANT_PAIR_COUNT=0) for degree in DEGREES]
    artifacts = {key: tmp_path / (key + ".tsv") for key in ("kinship_pairs", "independent_set_proposal", "kinship_degree_summary", "pedigree_concordance")}
    for key, rows in (("kinship_pairs", [pair]), ("independent_set_proposal", proposal), ("kinship_degree_summary", summary)):
        write_table(artifacts[key], tuple(rows[0]), rows)
    write_table(artifacts["pedigree_concordance"], ("RELATION_ID",), [])
    artifacts["kinship_report"] = tmp_path / "report.json"
    write_json(artifacts["kinship_report"], dict(proposed_inclusion_count=2, reported_pair_count=1, manual_validation_required=True))
    unit = dict(sample_indices=[0, 1], degree="FOURTH", stratum_id="ARTIFICIAL_A", theoretical_phi=1 / 32, cryptic=False)
    result = evaluate_kinship_stage(matrix, registry, qc, artifacts, [unit], ("A", "B"), config)
    assert result["functional_concordance"] == "PASS"
    assert result["status"] == "FAIL" and result["reasons"] == ["true_related_pair_in_independent_proposal"]
    assert result["performance"]["confusion"]["FOURTH"]["UNRELATED"] == 1
    assert result["frozen_cohort_09_evaluated"] is False
