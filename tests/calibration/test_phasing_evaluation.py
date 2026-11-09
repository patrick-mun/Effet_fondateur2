"""Contrats, reprise et clôture D23 sur fixtures minuscules, sans outil réel."""
from copy import deepcopy
import json
from pathlib import Path

import numpy as np
import pytest

from effet_fondateur.audit import sha256_file
from effet_fondateur.calibration.evaluate_phasing import acceptance_results, checked_closure, evaluate_autosomal_input, evaluate_region
from effet_fondateur.calibration.phasing_bench import METHOD_ID
from effet_fondateur.calibration.qc_bench import file_inventory, write_json, write_table
from effet_fondateur.calibration.run_phasing_regional import bind_dataset_artifacts, canonical_genotype_audit
from effet_fondateur.calibration.technical_evaluation import GenotypeMatrix
from effet_fondateur.contracts import validate_tsv_table


def test_binding_fixture_set_identifiers_without_mutating_inputs(tmp_path: Path) -> None:
    artifacts = []
    for prefix in ("genomewide_pre_qc", "target_chromosome_all_qc"):
        records = {}
        for suffix in ("bed", "bim", "fam"):
            path = tmp_path / (prefix + "." + suffix)
            path.write_text("synthetic fixture")
            records[suffix] = dict(name=path.name, sha256=sha256_file(path))
            artifacts.append(dict(artifact_id=prefix + "_" + suffix, path=str(path), sha256=sha256_file(path), sample_set_id=None, variant_set_id=None))
        descriptor = dict(schema_version="1.0.0", dataset_id=prefix, assembly="GRCh38", scope="autosomal_genomewide", target_chromosome=None,
                          sample_set_id="samples_fixed", variant_set_id=prefix + "_variants", sample_count=2, variant_count=3, marker_mode="intersection",
                          source_format="PLINK_PRELIMINARY_QC", allele_orientation="forward_strand_calls", files=records)
        path = tmp_path / (prefix + ".json")
        write_json(path, descriptor)
        artifacts.append(dict(artifact_id=prefix + "_dataset", path=str(path), sha256=sha256_file(path), sample_set_id=None, variant_set_id=None))
    original = deepcopy(artifacts)
    bound = bind_dataset_artifacts(artifacts)
    assert all(item["sample_set_id"] == "samples_fixed" and item["variant_set_id"] for item in bound)
    assert artifacts == original
    Path(artifacts[0]["path"]).write_text("changed")
    altered = deepcopy(artifacts)
    altered[0]["sha256"] = sha256_file(Path(altered[0]["path"]))
    with pytest.raises(ValueError, match="descriptor_artifact_mismatch"):
        bind_dataset_artifacts(altered)


def test_canonical_audit_reorders_only_columns_and_validates(tmp_path: Path) -> None:
    row = dict(SAMPLE_ID="SYNTHETIC1", GENOTYPE="A/G", GENOTYPE_SOURCE="SIMULATION", MEASUREMENT_METHOD="NOT_MEASURED", FID="F1", IID="SYNTHETIC1", STATUS="ACCEPTED")
    path = tmp_path / "initial.tsv"
    write_table(path, tuple(row), [row])
    with pytest.raises(ValueError):
        validate_tsv_table(path, "target_genotype_audit.schema.json")
    artifacts = [dict(artifact_id="target_genotype_audit", path=str(path), sha256=sha256_file(path))]
    result = canonical_genotype_audit(artifacts, tmp_path / "canonical")
    actual = Path(result[0]["path"])
    assert actual.read_text().splitlines()[0] == "SAMPLE_ID\tFID\tIID\tGENOTYPE\tGENOTYPE_SOURCE\tMEASUREMENT_METHOD\tSTATUS"
    assert artifacts[0]["path"] == str(path)
    validate_tsv_table(actual, "target_genotype_audit.schema.json")


def test_closed_inventory_rejects_addition_tamper_and_symlink(tmp_path: Path) -> None:
    data = tmp_path / "aggregate.json"
    data.write_text('{}')
    document = dict(schema_version="1.0.0", method_id=METHOD_ID, closure_state="CLOSED", calibrated=False, founder_effect_proven=False,
                    artifacts=file_inventory(tmp_path))
    write_json(tmp_path / "run_index.json", document)
    assert checked_closure(tmp_path) == document
    extra = tmp_path / "extra"
    extra.write_text("added")
    with pytest.raises(ValueError, match="uninventoried"):
        checked_closure(tmp_path)
    extra.unlink()
    data.write_text('{"changed":true}')
    with pytest.raises(ValueError, match="artifact_changed"):
        checked_closure(tmp_path)
    data.unlink()
    data.symlink_to(tmp_path / "run_index.json")
    with pytest.raises(ValueError):
        checked_closure(tmp_path)


def test_10a_oracle_complete_polymorphic_filter_and_mismatch(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    samples = (("F1", "S1"), ("F2", "S2"))
    variants = tuple((1, index + 1, "m" + str(index)) for index in range(4))
    calls = np.array([[[0, 0], [0, 2]], [[0, 0], [0, 0]], [[-1, -1], [0, 2]], [[2, 2], [0, 2]]], dtype=np.int8)
    source = GenotypeMatrix(samples, variants, calls)
    expected = GenotypeMatrix(samples, (variants[0], variants[3]), calls[[0, 3]])
    (tmp_path / "source.fam").write_text("F1 S1 0 0 0 -9\nF2 S2 0 0 0 -9\n")
    monkeypatch.setattr("effet_fondateur.calibration.evaluate_phasing.read_plink", lambda path: source if path.name == "source" else expected)
    result = evaluate_autosomal_input(tmp_path / "source", tmp_path / "actual")
    assert result["status"] == "PASS"
    assert result["expected_markers"] == 2
    assert result["excluded_monomorphic"] == result["excluded_incomplete"] == 1
    wrong = GenotypeMatrix(samples, expected.variants, expected.alleles.copy())
    wrong.alleles[0, 0] = 2
    monkeypatch.setattr("effet_fondateur.calibration.evaluate_phasing.read_plink", lambda path: source if path.name == "source" else wrong)
    assert evaluate_autosomal_input(tmp_path / "source", tmp_path / "actual")["status"] == "FAIL"


def test_11_oracle_known_interpolation_and_target_loss(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    samples = (("F1", "S1"),)
    variants = ((5, 100, "m1"), (5, 150, "target"), (5, 200, "m3"))
    matrix = GenotypeMatrix(samples, variants, np.array([[[0, 2]]] * 3, dtype=np.int8))
    monkeypatch.setattr("effet_fondateur.calibration.evaluate_phasing.read_plink", lambda path: matrix)
    source_map, output_map = tmp_path / "source.tsv", tmp_path / "actual.tsv"
    write_table(source_map, ("POSITION_BP", "POSITION_CM"), [dict(POSITION_BP=100, POSITION_CM=1), dict(POSITION_BP=200, POSITION_CM=3)])
    rows = [dict(POSITION_BP=100, POSITION_CM=1), dict(POSITION_BP=150, POSITION_CM=2), dict(POSITION_BP=200, POSITION_CM=3)]
    write_table(output_map, tuple(rows[0]), rows)
    config = dict(target_flank_bp=50, map_absolute_tolerance_cm=1e-8)
    result = evaluate_region(Path("source"), Path("actual"), output_map, source_map, dict(chromosome=5, position_bp=150), config)
    assert result["status"] == "PASS" and result["map_max_absolute_error_cm"] == 0
    rows[1]["POSITION_CM"] = 2.1
    write_table(output_map, tuple(rows[0]), rows)
    result = evaluate_region(Path("source"), Path("actual"), output_map, source_map, dict(chromosome=5, position_bp=150), config)
    assert result["status"] == "FAIL" and result["map_coordinate_error_count"] == 1


def test_acceptance_keeps_measurability_separate_from_unsealed_qualification() -> None:
    result = dict(target_retained=True, metrics={name: dict(ALL=dict(estimate=.1,denominator=100)) for name in
                  ("switch_error_rate", "rare_variant_phase_error_rate", "target_mutation_orientation_accuracy")})
    criteria = acceptance_results("12",result)
    assert [item["criterion_id"] for item in criteria] == ["CAL-PHASE-003", "CAL-PHASE-004", "CAL-PHASE-005", "CAL-PHASE-006"]
    assert criteria[0]["status"] == "PASS"
    assert all(item["status"] == "REVIEW_REQUIRED" and item["threshold"] is None and item["denominator"] == 100 for item in criteria[1:])
    result["metrics"]["rare_variant_phase_error_rate"]["ALL"] = dict(estimate=None,denominator=0)
    assert acceptance_results("12",result)[2]["status"] == "NOT_EVALUABLE"
