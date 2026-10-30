"""D19 : oracles indépendants, matrices BED, défauts injectés et sorties 00–02."""

from __future__ import annotations

from copy import deepcopy
import hashlib
import json
from pathlib import Path
import shutil

import numpy as np
import pytest

from effet_fondateur.calibration.evaluate_technical import (
    _dataset, compare_fam, contained_file, evaluate_pipeline_run, evaluate_technical,
    load_checked_evaluation, validate_schema,
)
from effet_fondateur.calibration.export_blinded_inputs import export_blinded_inputs, load_checked_export
from effet_fondateur.calibration.run_technical_inputs import run_technical_inputs
from effet_fondateur.calibration.report_technical_evaluation import build_reports
from effet_fondateur.calibration.technical_evaluation import (
    GenotypeMatrix, TechnicalEvaluationError, compare_genotypes, compare_rows,
    evaluate_target, read_plink, row_result,
)
from test_blinded_export import _fixture as export_fixture


def _matrix() -> GenotypeMatrix:
    return GenotypeMatrix((("F1", "I1"), ("F2", "I2"), ("F3", "I3")),
                          ((5, 100, "chip1"), (5, 200, "chip2")),
                          np.array([[[0, 0], [0, 2], [-1, -1]], [[1, 3], [3, 3], [1, 1]]], dtype=np.int8))


def write_plink(prefix: Path, matrix: GenotypeMatrix, *, reverse_alleles: bool = False) -> None:
    """Encodeur de fixture ; le décodeur a aussi un oracle binaire officiel séparé."""
    prefix.parent.mkdir(parents=True, exist_ok=True)
    prefix.with_suffix(".fam").write_text("".join(f"{fid} {iid} 0 0 0 -9\n" for fid, iid in matrix.samples))
    bim = []
    payload = bytearray(b"\x6c\x1b\x01")
    for index, (chromosome, position, name) in enumerate(matrix.variants):
        values = sorted({int(value) for value in matrix.alleles[index].ravel() if value >= 0})
        alleles = ["ACGT"[value] for value in values]
        alleles = (alleles + ["0", "0"])[:2]
        if reverse_alleles: alleles.reverse()
        bim.append(f"{chromosome} {name} 0 {position} {alleles[0]} {alleles[1]}\n")
        block = bytearray((len(matrix.samples) + 3) // 4)
        for sample, pair in enumerate(matrix.alleles[index]):
            if pair[0] < 0: code = 1
            else:
                bases = tuple("ACGT"[int(value)] for value in pair)
                code = 0 if bases == (alleles[0], alleles[0]) else 3 if bases == (alleles[1], alleles[1]) else 2
            block[sample // 4] |= code << (2 * (sample % 4))
        payload.extend(block)
    prefix.with_suffix(".bim").write_text("".join(bim))
    prefix.with_suffix(".bed").write_bytes(payload)


def test_official_plink_bed_blocks_are_an_independent_oracle(tmp_path):
    # Premier et troisième blocs du petit exemple officiel PLINK 1.9 BED.
    prefix = tmp_path / "official"
    prefix.with_suffix(".fam").write_text("".join(f"F I{index} 0 0 0 -9\n" for index in range(6)))
    prefix.with_suffix(".bim").write_text("1 snp1 0 1 G A\n1 snp3 0 3 A C\n")
    prefix.with_suffix(".bed").write_bytes(bytes.fromhex("6c1b01dc0f6b01"))
    matrix = read_plink(prefix)
    assert matrix.alleles.tolist() == [[[2, 2], [0, 0], [-1, -1], [0, 0], [0, 0], [0, 0]],
                                      [[1, 1], [0, 1], [0, 1], [-1, -1], [-1, -1], [0, 0]]]


@pytest.mark.parametrize("reverse", [False, True])
def test_bed_roundtrip_compares_bases_not_a1_dosage(tmp_path, reverse):
    prefix = tmp_path / "panel"
    write_plink(prefix, _matrix(), reverse_alleles=reverse)
    assert compare_genotypes(_matrix(), read_plink(prefix))["status"] == "PASS"


@pytest.mark.parametrize("mutation", ["header", "truncated", "extra", "padding", "duplicate_sample", "duplicate_variant", "called_zero"])
def test_malformed_plink_is_never_accepted(tmp_path, mutation):
    prefix = tmp_path / "panel"
    write_plink(prefix, _matrix())
    payload = bytearray(prefix.with_suffix(".bed").read_bytes())
    if mutation == "header": payload[2] = 0
    elif mutation == "truncated": payload = payload[:-1]
    elif mutation == "extra": payload += b"\x00"
    elif mutation == "padding": payload[-1] |= 192
    elif mutation == "duplicate_sample":
        rows = prefix.with_suffix(".fam").read_text().splitlines()
        rows[1] = rows[0]
        prefix.with_suffix(".fam").write_text("\n".join(rows) + "\n")
    elif mutation == "duplicate_variant":
        rows = prefix.with_suffix(".bim").read_text().splitlines()
        rows[1] = rows[0]
        prefix.with_suffix(".bim").write_text("\n".join(rows) + "\n")
    else:
        text = prefix.with_suffix(".bim").read_text().replace("A G", "0 G")
        prefix.with_suffix(".bim").write_text(text)
    prefix.with_suffix(".bed").write_bytes(payload)
    with pytest.raises(TechnicalEvaluationError): read_plink(prefix)


@pytest.mark.parametrize("mutation", ["lost", "extra", "duplicate", "order", "value"])
def test_table_detector_catches_every_injected_defect(mutation):
    expected = [{"id": "a", "value": "one"}, {"id": "b", "value": "two"}]
    actual = deepcopy(expected)
    if mutation == "lost": actual.pop()
    elif mutation == "extra": actual.append({"id": "c", "value": "three"})
    elif mutation == "duplicate": actual.append(dict(actual[0]))
    elif mutation == "order": actual.reverse()
    else: actual[0]["value"] = "wrong"
    assert row_result("02", expected, actual, keys=("id",), fields=("id", "value"))["status"] == "FAIL"
    assert not any(token in json.dumps(row_result("02", expected, actual, keys=("id",), fields=("value",))) for token in ("one", "two", "wrong"))


def test_genotype_substitution_and_missing_imputation_are_detected():
    expected = _matrix()
    for marker, sample, pair in ((0, 0, (2, 2)), (0, 2, (0, 0))):
        changed = expected.alleles.copy()
        changed[marker, sample] = pair
        verdict = compare_genotypes(expected, GenotypeMatrix(expected.samples, expected.variants, changed))
        assert verdict["status"] == "FAIL" and verdict["metrics"]["genotype_mismatch_count"] == 1
        assert verdict["metrics"]["genotype_discordance_rate"] == 1 / 6


def test_reordered_axes_are_detected_even_with_correct_realigned_genotypes():
    expected = _matrix()
    changed = GenotypeMatrix(expected.samples[::-1], expected.variants, expected.alleles[:, ::-1])
    verdict = compare_genotypes(expected, changed)
    assert verdict["status"] == "FAIL"
    assert verdict["metrics"]["sample_order_mismatch_count"] > 0
    assert verdict["metrics"]["genotype_mismatch_count"] == 0


@pytest.mark.parametrize("axis", ["sample", "variant"])
def test_missing_axis_never_reports_zero_discordance_rate(axis):
    expected = _matrix()
    actual = GenotypeMatrix(expected.samples[:-1], expected.variants, expected.alleles[:, :-1]) if axis == "sample" else GenotypeMatrix(expected.samples, expected.variants[:-1], expected.alleles[:-1])
    verdict = compare_genotypes(expected, actual)
    assert verdict["status"] == "FAIL"
    assert verdict["metrics"]["genotype_discordance_rate"] is None


def test_target_evaluator_checks_target_and_unchanged_background():
    base = _matrix()
    target = GenotypeMatrix(base.samples, ((5, 300, "target"),), np.array([[[0, 2], [2, 2], [0, 0]]], dtype=np.int8))
    enriched = GenotypeMatrix(base.samples, base.variants + target.variants, np.concatenate((base.alleles, target.alleles)))
    assert evaluate_target(base, target, enriched)["status"] == "PASS"
    for marker in (0, 2):
        changed = enriched.alleles.copy()
        changed[marker, 0] = (3, 3)
        verdict = evaluate_target(base, target, GenotypeMatrix(base.samples, enriched.variants, changed))
        assert verdict["status"] == "FAIL"
        assert verdict["metrics"]["target_genotype_discordance_rate"] == (1 / 3 if marker == 2 else 0.0)
    assert evaluate_target(base, None, enriched)["status"] == "NOT_EVALUABLE"
    assert evaluate_target(base, target, None)["status"] == "NOT_EVALUABLE"


@pytest.fixture(scope="module")
def closed_inputs(tmp_path_factory):
    root = tmp_path_factory.mktemp("d19")
    arguments = export_fixture(root)
    export_blinded_inputs(**arguments)
    d18 = arguments["workspace_root"] / "d18_test"
    index = run_technical_inputs(d18_root=d18, workspace_root=arguments["workspace_root"], run_id="pipeline_d19")
    return root, arguments["workspace_root"], d18, index


def test_real_partial_pipeline_gives_pass_00_02_and_not_evaluable_03_04(closed_inputs):
    _, workspace, d18, index = closed_inputs
    manifest = evaluate_technical(d18_root=d18, run_index_path=workspace / "pipeline_d19/run_index.json", workspace_root=workspace,
                                  run_id="evaluation_d19", created_at="2026-10-30T00:00:00Z")
    assert [item["status"] for item in manifest["stages"]] == ["PASS", "PASS", "PASS", "NOT_EVALUABLE", "NOT_EVALUABLE"]
    assert all(item["pass_count"] == 8 for item in manifest["stages"][:3])
    assert load_checked_evaluation(workspace / "evaluation_d19") == manifest
    with pytest.raises(TechnicalEvaluationError, match="run_exists"):
        evaluate_technical(d18_root=d18, run_index_path=workspace / "pipeline_d19/run_index.json", workspace_root=workspace,
                           run_id="evaluation_d19", created_at="2026-10-30T00:00:00Z")
    load_checked_export(d18)


def _refresh_run(run: Path, stage: str, receipt: dict) -> None:
    manifest = json.loads((run / "manifest.json").read_text())
    entry = next(item for item in manifest["stages"] if item["stage_id"] == stage)
    audit_path = run / entry["audit_path"]
    directory = audit_path.parent
    outputs = json.loads((directory / "stage_outputs.json").read_text())
    for artifact in outputs["artifacts"]:
        path = run / artifact["path"]
        artifact["sha256"] = hashlib.sha256(path.read_bytes()).hexdigest()
    audit = json.loads(audit_path.read_text())
    audit["outputs"] = outputs["artifacts"]
    audit_path.write_text(json.dumps(audit))
    (directory / "stage_outputs.json").write_text(json.dumps(outputs))
    entry["audit_sha256"] = hashlib.sha256(audit_path.read_bytes()).hexdigest()
    entry["stage_outputs_sha256"] = hashlib.sha256((directory / "stage_outputs.json").read_bytes()).hexdigest()
    (run / "manifest.json").write_text(json.dumps(manifest))
    receipt["manifest_sha256"] = hashlib.sha256((run / "manifest.json").read_bytes()).hexdigest()


@pytest.mark.parametrize("mutation", ["lost_source", "duplicate_source", "substituted_hash", "registry_order", "registry_value", "summary"])
def test_pipeline_output_mutations_detected_even_when_local_hashes_are_resealed(closed_inputs, tmp_path, mutation):
    _, workspace, d18, index = closed_inputs
    receipt = deepcopy(index["runs"][0])
    run = tmp_path / "run"
    shutil.copytree(workspace / receipt["run_path"], run)
    bundle = d18 / "blinded_inputs/bundles" / receipt["bundle_id"]
    inventory = json.loads((d18 / "blinded_inputs/inventory.json").read_text())
    stage = "00" if mutation == "summary" else "02" if mutation.startswith("registry") else "01"
    if stage == "00":
        path = run / "stages/00_initialize_run/initialization_summary.json"
        summary = json.loads(path.read_text())
        summary["enabled_stage_count"] = 999
        path.write_text(json.dumps(summary))
    else:
        path = run / ("stages/01_validate_sources/source_inventory.tsv" if stage == "01" else "stages/02_build_sample_registry/samples.master.tsv")
        lines = path.read_text().splitlines()
        if mutation == "lost_source": lines.pop()
        elif mutation == "duplicate_source": lines.append(lines[1])
        elif mutation == "substituted_hash":
            fields = lines[1].split("\t")
            fields[2] = "0" * 64
            lines[1] = "\t".join(fields)
        elif mutation == "registry_order": lines = lines[:1] + lines[:0:-1]
        else:
            fields = lines[1].split("\t")
            fields[8] = "WRONG_GROUP"
            lines[1] = "\t".join(fields)
        path.write_text("\n".join(lines) + "\n")
    _refresh_run(run, stage, receipt)
    verdicts = evaluate_pipeline_run(bundle_root=bundle, inventory=inventory, bundle_id=receipt["bundle_id"], run_root=run, receipt=receipt)
    assert next(item for item in verdicts if item["stage_id"] == stage)["status"] == "FAIL"


def test_missing_execution_is_not_a_pass(closed_inputs):
    _, _, d18, index = closed_inputs
    receipt = index["runs"][0]
    inventory = json.loads((d18 / "blinded_inputs/inventory.json").read_text())
    verdicts = evaluate_pipeline_run(bundle_root=d18 / "blinded_inputs/bundles" / receipt["bundle_id"], inventory=inventory,
                                     bundle_id=receipt["bundle_id"], run_root=None, receipt=None)
    assert {item["status"] for item in verdicts} == {"NOT_EVALUABLE"}
    assert all(item["metrics"]["evaluated_count"] is None for item in verdicts)


@pytest.mark.parametrize("relative", ["../outside", "/absolute", "folder/../../escape", "back\\slash"])
def test_unsafe_relative_path_is_refused_before_read(tmp_path, relative):
    with pytest.raises(TechnicalEvaluationError, match="unsafe_path"):
        contained_file(tmp_path, relative)


def test_symlinks_and_hardlinks_are_refused(tmp_path):
    source = tmp_path / "source"
    source.write_text("fixture")
    (tmp_path / "link").symlink_to(source)
    with pytest.raises(TechnicalEvaluationError, match="unsafe_path"):
        contained_file(tmp_path, "link")
    (tmp_path / "hard").hardlink_to(source)
    with pytest.raises(TechnicalEvaluationError, match="unsafe_path"):
        contained_file(tmp_path, "hard")


def test_report_renders_six_independent_documents(closed_inputs, tmp_path):
    _, workspace, d18, _ = closed_inputs
    evaluate_technical(d18_root=d18, run_index_path=workspace / "pipeline_d19/run_index.json", workspace_root=workspace,
                       run_id="report_source", created_at="2026-10-30T00:00:00Z")
    template = tmp_path / "template.html"
    template.write_text("<!doctype html><html><head><title>Template</title></head><body></body></html>")
    paths = build_reports(run_root=workspace / "report_source", template_path=template, output_root=workspace / "reports")
    assert len(paths) == 6
    for path in paths:
        html = path.read_text()
        assert html.count('<section id="s') == 14
        assert html.count('<figure>') == 2
        assert html.count('<table>') == 1
        assert 'calibrated=false' in html and '[To complete:' not in html
        assert not any(token in html for token in ('sample_id', 'genotype":', 'SAMPLE_ID'))
    with pytest.raises(TechnicalEvaluationError, match="report_output_invalid"):
        build_reports(run_root=workspace / "report_source", template_path=template, output_root=workspace / "reports")


def test_evaluation_reader_rejects_resealed_wrong_summary(closed_inputs, tmp_path):
    _, workspace, d18, _ = closed_inputs
    evaluate_technical(d18_root=d18, run_index_path=workspace / "pipeline_d19/run_index.json", workspace_root=workspace,
                       run_id="reader_source", created_at="2026-10-30T00:00:00Z")
    run = tmp_path / "evaluation"
    shutil.copytree(workspace / "reader_source", run)
    path = run / "private_truth/d19/technical_evaluation_manifest.json"
    manifest = json.loads(path.read_text())
    manifest['stages'][0]['pass_count'] = 0
    path.write_text(json.dumps(manifest))
    with pytest.raises(TechnicalEvaluationError, match="summary_mismatch"):
        load_checked_evaluation(run)


@pytest.mark.parametrize("corruption", ["manifest", "resealed_invalid_manifest", "artifact", "missing_artifact"])
def test_successful_stage_with_corrupted_or_missing_output_fails(closed_inputs, tmp_path, corruption):
    _, workspace, d18, index = closed_inputs
    receipt = deepcopy(index["runs"][0])
    run = tmp_path / "run"
    shutil.copytree(workspace / receipt["run_path"], run)
    if corruption == "manifest":
        (run / "manifest.json").write_text("{}")
    elif corruption == "resealed_invalid_manifest":
        (run / "manifest.json").write_text("not JSON")
        receipt["manifest_sha256"] = hashlib.sha256((run / "manifest.json").read_bytes()).hexdigest()
    else:
        path = run / "stages/02_build_sample_registry/samples.master.tsv"
        if corruption == "artifact": path.write_text("wrong")
        else: path.unlink()
    verdicts = evaluate_pipeline_run(bundle_root=d18 / "blinded_inputs/bundles" / receipt["bundle_id"],
                                     inventory=json.loads((d18 / "blinded_inputs/inventory.json").read_text()),
                                     bundle_id=receipt["bundle_id"], run_root=run, receipt=receipt)
    assert verdicts[2]["status"] == "FAIL"


def test_dataset_adapter_checks_descriptor_and_decoded_genotypes(tmp_path):
    prefix = tmp_path / "base"
    write_plink(prefix, _matrix())
    descriptor = dict(schema_version="1.0.0", dataset_id="base", assembly="GRCh38", scope="target_chromosome",
                      target_chromosome=5, sample_set_id="samples", variant_set_id="markers", sample_count=3, variant_count=2,
                      marker_mode="union", source_format="ChAS_ACPA", allele_orientation="forward_strand_calls",
                      files={suffix:dict(name=prefix.with_suffix('.' + suffix).name, sha256=hashlib.sha256(prefix.with_suffix('.' + suffix).read_bytes()).hexdigest()) for suffix in ("bed", "bim", "fam")})
    descriptor_path = tmp_path / "dataset.json"
    descriptor_path.write_text(json.dumps(descriptor))
    artifacts = {"base_" + suffix: prefix.with_suffix('.' + suffix) for suffix in ("bed", "bim", "fam")}
    artifacts["base_dataset"] = descriptor_path
    assert compare_genotypes(_matrix(), _dataset(artifacts, "base"))["status"] == "PASS"
    descriptor['sample_count'] = 99
    descriptor_path.write_text(json.dumps(descriptor))
    with pytest.raises(TechnicalEvaluationError, match="descriptor_count_mismatch"):
        _dataset(artifacts, "base")


def test_running_process_is_rejected_before_private_export_is_opened(closed_inputs, monkeypatch):
    _, workspace, d18, source_index = closed_inputs
    run = workspace / "not_closed"
    receipt = deepcopy(source_index["runs"][0])
    shutil.copytree(workspace / receipt["run_path"], run)
    path = run / "manifest.json"
    document = json.loads(path.read_text())
    document["stages"][0]["state"] = "RUNNING"
    path.write_text(json.dumps(document))
    receipt["run_path"] = "not_closed"
    receipt["manifest_sha256"] = hashlib.sha256(path.read_bytes()).hexdigest()
    index = deepcopy(source_index)
    index["runs"] = [receipt]
    index_path = workspace / "active_index.json"
    index_path.write_text(json.dumps(index))
    def forbidden_private_access(_):
        pytest.fail("private export opened before closure check")
    monkeypatch.setattr("effet_fondateur.calibration.evaluate_technical.load_checked_export", forbidden_private_access)
    with pytest.raises(TechnicalEvaluationError, match="run_not_closed"):
        evaluate_technical(d18_root=d18, run_index_path=index_path, workspace_root=workspace,
                           run_id="active_evaluation", created_at="2026-10-30T00:00:00Z")
    assert not (workspace / "active_evaluation").exists()


@pytest.mark.parametrize("field", ["PID", "MID", "SEX", "CLINICAL_STATUS"])
def test_fam_metadata_substitutions_are_detected(tmp_path, field):
    path = tmp_path / "samples.fam"
    path.write_text("F I 0 0 0 -9\n")
    metadata = [dict(FID="F", IID="I", PID="0", MID="0", SEX="UNKNOWN", CLINICAL_STATUS="UNKNOWN")]
    assert compare_fam(path, metadata, stage="03")["status"] == "PASS"
    metadata[0][field] = {"PID":"father", "MID":"mother", "SEX":"MALE", "CLINICAL_STATUS":"AFFECTED"}[field]
    assert compare_fam(path, metadata, stage="03")["status"] == "FAIL"


@pytest.mark.parametrize("missing", [False, True])
def test_monomorphic_or_all_missing_variant_preserves_absence(tmp_path, missing):
    expected = _matrix()
    expected.alleles[0] = -1 if missing else 0
    prefix = tmp_path / "monomorphic"
    write_plink(prefix, expected)
    assert compare_genotypes(expected, read_plink(prefix))["status"] == "PASS"
