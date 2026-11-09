"""D15 : identité injectée, fond négatif inchangé et protection de Mendel."""

from __future__ import annotations

from dataclasses import asdict, replace
import gzip
import hashlib
import io
import json
from pathlib import Path

import jsonschema
import pytest

from effet_fondateur.calibration.simulation.founder_recombination import GeneticMapPoint
from effet_fondateur.calibration.simulation.ibd import (
    IbdSimulationError, build_ibd_injections, transform_ibd_vcf,
)
from effet_fondateur.calibration.materialize_ibd import materialize_ibd
from effet_fondateur.calibration.report_ibd_injection import build_report


SAMPLES = [f"SYN_CAL_{index:03d}" for index in range(30)]
POINTS = (GeneticMapPoint(100, 0), GeneticMapPoint(600, 5), GeneticMapPoint(1100, 10))
SCENARIOS = [{"scenario_id": f"length_{index}", "segment_length_cm": length,
              "minimum_observed_markers": 1} for index, length in enumerate((1, 1.5, 2, 3, 5))]


def _plan(seed: int = 1004):
    return build_ibd_injections(eligible_samples=SAMPLES, map_points=POINTS,
                               chromosome=5, target_bp=600, scenarios=SCENARIOS,
                               random_seed=seed)


def _vcf():
    rows = ["##fileformat=VCFv4.2", "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t" + "\t".join(SAMPLES)]
    for position in range(100, 1101, 10):
        rows.append(f"chr5\t{position}\tv{position}\tA\tG\t.\tPASS\t.\tGT\t" + "\t".join(
            f"{(index + position // 10) % 2}|{index % 2}" for index in range(30)))
    return "\n".join(rows) + "\n"


def _keys():
    return {("chr5", position, "A", "G"): f"CHIP_{position}"
            for position in range(100, 1101, 10)}


def _rows(text: str):
    return {int(line.split("\t")[1]): line.split("\t")[9:]
            for line in text.splitlines() if not line.startswith("#")}


def test_plan_is_deterministic_disjoint_and_negative_bounds_null():
    plan = _plan()
    assert plan == _plan()
    assert plan != _plan(1005)
    participants = [sample for item in plan for sample in
                    (item.sample_1_id, item.sample_2_id, item.template_id) if sample]
    assert len(participants) == len(set(participants)) == 25
    assert sum(item.is_positive for item in plan) == 5
    assert all(item.start_bp is None and item.end_bp is None
               for item in plan if not item.is_positive)


def test_truth_validates_schema():
    schema = json.loads((Path(__file__).resolve().parents[2] /
                         "calibration/schemas/ibd_injection_truth.schema.json").read_text())
    jsonschema.Draft202012Validator(schema).validate([
        {"replication_id": "test", **asdict(item)} for item in _plan()])
    bad = [{"replication_id": "test", **asdict(item)} for item in _plan()]
    bad[1]["start_bp"] = 100
    with pytest.raises(jsonschema.ValidationError):
        jsonschema.validate(bad, schema)


def test_positive_copy_negative_unchanged_and_outside_unchanged():
    plan = _plan()
    dense, observed = io.StringIO(), io.StringIO()
    result = transform_ibd_vcf(io.StringIO(_vcf()), dense, observed,
                               injections=plan, observed_keys=_keys(), transmissions=())
    original, output = _rows(_vcf()), _rows(dense.getvalue())
    assert result["dense_variant_count"] == 101
    assert result["observed_variant_count"] == 101
    for item in plan:
        first, second = SAMPLES.index(item.sample_1_id), SAMPLES.index(item.sample_2_id)
        for position, genotypes in output.items():
            if item.is_positive and item.start_bp <= position <= item.end_bp:
                template = SAMPLES.index(item.template_id)
                assert genotypes[first][0] == genotypes[second][0] == original[position][template][0]
                assert genotypes[first][-1] == original[position][first][-1]
            else:
                assert genotypes[first] == original[position][first]
                assert genotypes[second] == original[position][second]
    assert "_positive" not in observed.getvalue()
    assert "_negative" not in observed.getvalue()
    assert "|" not in observed.getvalue()


def test_mendel_transmission_remains_exact():
    # Two unused synthetic subjects with identical GT at every site.
    used = {sample for item in _plan() for sample in
            (item.sample_1_id, item.sample_2_id, item.template_id) if sample}
    remaining = [sample for sample in SAMPLES if sample not in used]
    first, second = next((a, b) for a in remaining for b in remaining
                         if a != b and SAMPLES.index(a) % 2 == SAMPLES.index(b) % 2)
    result = transform_ibd_vcf(io.StringIO(_vcf()), io.StringIO(), io.StringIO(),
                               injections=_plan(), observed_keys=_keys(),
                               transmissions=((first, second, 1, 1),))
    assert result["mendelian_error_count"] == 0


def test_pedigree_member_cannot_be_modified():
    with pytest.raises(IbdSimulationError, match="pedigree_injection_forbidden"):
        transform_ibd_vcf(io.StringIO(_vcf()), io.StringIO(), io.StringIO(),
                          injections=_plan(), observed_keys=_keys(),
                          transmissions=((_plan()[0].sample_1_id, SAMPLES[-1], 1, 1),))


def test_pair_overlap_is_rejected():
    plan = list(_plan())
    plan[1] = replace(plan[1], sample_1_id=plan[0].sample_1_id)
    with pytest.raises(IbdSimulationError, match="pair_or_template_overlap"):
        transform_ibd_vcf(io.StringIO(_vcf()), io.StringIO(), io.StringIO(),
                          injections=plan, observed_keys=_keys(), transmissions=())


def test_insufficient_observation_is_not_accepted():
    with pytest.raises(IbdSimulationError, match="marker_density_insufficient"):
        transform_ibd_vcf(io.StringIO(_vcf()), io.StringIO(), io.StringIO(),
                          injections=_plan(), observed_keys={}, transmissions=())


def test_window_cannot_be_enlarged_silently():
    with pytest.raises(IbdSimulationError, match="segment_outside_domain"):
        build_ibd_injections(eligible_samples=SAMPLES, map_points=POINTS,
                             chromosome=5, target_bp=150, scenarios=SCENARIOS,
                             random_seed=1004)


def test_missing_and_unphased_genotypes_are_rejected():
    for replacement in ("./.", "0/1"):
        with pytest.raises(IbdSimulationError, match="gt_invalid"):
            transform_ibd_vcf(io.StringIO(_vcf().replace("0|0", replacement)),
                              io.StringIO(), io.StringIO(), injections=_plan(),
                              observed_keys=_keys(), transmissions=())


def _materialization_fixture(tmp_path: Path):
    workspace = tmp_path / "calibration_workspace"
    roots = {key: workspace / key for key in ("d14", "d11", "design", "replay")}
    def write(root, relative, text):
        path = root / relative
        path.parent.mkdir(parents=True, exist_ok=True)
        if path.suffix == ".gz":
            # Les entrées du test de rejeu ne doivent pas dépendre de l'heure.
            with path.open("wb") as raw:
                with gzip.GzipFile(filename="", fileobj=raw, mode="wb", mtime=0) as handle:
                    handle.write(text.encode("utf-8"))
        else:
            path.write_text(text)
        return path
    def digest(path):
        return hashlib.sha256(path.read_bytes()).hexdigest()
    domain = dict(chromosome=5, start_bp=100, end_bp=1100, start_cm=0,
                  end_cm=10, target_position_bp=600, target_position_cm=5,
                  position_unit="BP_1_BASED_CLOSED", genetic_distance_unit="cM")
    artifacts = []
    payloads = {
        "DENSE_PHASED_FOUNDER_PANEL": _vcf(),
        "BLINDED_OBSERVED_FOUNDER_PANEL": _vcf(),
        "PEDIGREE_TRUTH": "sample_id\tfamily_unit_id\n" + "".join(f"{sample}\tUNRELATED_{sample}\n" for sample in SAMPLES),
        "FOUNDER_CARRIER_TRUTH": "sample_id\n",
        "MENDELIAN_TRANSMISSION_TRUTH": "child_id\tparent_id\tparental_haplotype_index\tchild_haplotype_index\n",
    }
    for index, (role, payload) in enumerate(payloads.items()):
        path = write(roots["d14"], f"private_truth/d14/artifact_{index}.gz", payload)
        artifacts.append(dict(role=role, path=path.relative_to(roots["d14"]).as_posix(), sha256=digest(path)))
    d14 = dict(manifest_state="FINALIZED", panel={"role": "CALIBRATION"}, domain=domain,
               configuration={"random_seed": 1004}, artifacts=artifacts)
    d14_path = write(roots["d14"], "private_truth/d14/founder_simulation_manifest.json", json.dumps(d14))
    map_path = write(roots["d11"], "private_truth/d11/map.tsv", "CHROMOSOME\tPOSITION_BP\tPOSITION_CM\n5\t100\t0\n5\t600\t5\n5\t1100\t10\n")
    d11 = dict(manifest_state="FINALIZED", panel={"role": "CALIBRATION"}, artifacts=[dict(role="GENETIC_MAP", path="private_truth/d11/map.tsv", sha256=digest(map_path))])
    d11_path = write(roots["d11"], "private_truth/d11/base_haplotype_manifest.json", json.dumps(d11))
    scenarios = [{**item, "status": "EVALUABLE"} for item in SCENARIOS]
    design = dict(manifest_state="FINALIZED", panel={"role": "CALIBRATION"}, windows=[dict(purpose="PRIMARY", target={"position_bp": 600}, ibd_scenarios=scenarios)])
    design_path = write(roots["design"], "private_truth/d14b/simulation_window_manifest.json", json.dumps(design))
    replay = dict(receipt_state="FINALIZED", stages=[dict(stage="D14", manifest_sha256=digest(d14_path)), dict(stage="D11A", manifest_sha256=digest(d11_path))], design=dict(manifest_sha256=digest(design_path)))
    write(roots["replay"], "private_truth/d14b_replay/simulation_window_replay.json", json.dumps(replay))
    return dict(d14_root=roots["d14"], d11a_root=roots["d11"], design_root=roots["design"], replay_root=roots["replay"], workspace_root=workspace,
                run_id="d15_test", created_at="2026-10-29T10:00:00Z")


def test_materialization_is_deterministic_and_refuses_overwrite(tmp_path, monkeypatch):
    monkeypatch.setattr(gzip.time, "time", lambda: 1000000000)
    arguments = _materialization_fixture(tmp_path)
    manifest = materialize_ibd(**arguments)
    assert manifest["positive_pair_count"] == manifest["negative_pair_count"] == 5
    assert manifest["calibrated"] is False
    with pytest.raises(IbdSimulationError, match="run_exists"):
        materialize_ibd(**arguments)
    monkeypatch.setattr(gzip.time, "time", lambda: 1000000002)
    second = _materialization_fixture(tmp_path / "second")
    assert materialize_ibd(**second) == manifest


def test_tampering_blocks_materialization(tmp_path):
    arguments = _materialization_fixture(tmp_path)
    source = arguments["d14_root"] / "private_truth/d14/artifact_0.gz"
    source.write_bytes(b"tampered")
    with pytest.raises(IbdSimulationError, match="artifact_hash_or_path_invalid"):
        materialize_ibd(**arguments)
    assert not (arguments["workspace_root"] / arguments["run_id"]).exists()


def test_real_input_path_is_refused_before_reading(tmp_path):
    arguments = _materialization_fixture(tmp_path / "data/input")
    with pytest.raises(IbdSimulationError, match="real_or_traversing_path_forbidden"):
        materialize_ibd(**arguments)


def test_report_retains_limits_and_contains_no_sample_truth(tmp_path):
    arguments = _materialization_fixture(tmp_path)
    materialize_ibd(**arguments)
    template = tmp_path / "template.html"
    template.write_text('<!doctype html><html lang="fr"><head><title>Modèle</title></head><body></body></html>')
    report = build_report(run_root=arguments["workspace_root"] / "d15_test",
                          template_path=template,
                          output_root=arguments["workspace_root"] / "d15_report")
    text = report.read_text()
    assert "calibrated=false" in text
    assert "NOT_EVALUABLE" in text
    assert "Tableau 1" in text
    assert all(sample not in text for sample in SAMPLES)
    assert text.count('<section id="s') == 14
    assert (report.parent / "checksums.sha256").is_file()


def test_duplicate_scenario_is_rejected():
    with pytest.raises(IbdSimulationError, match="duplicate_scenario"):
        build_ibd_injections(eligible_samples=SAMPLES, map_points=POINTS,
                             chromosome=5, target_bp=600,
                             scenarios=[SCENARIOS[0], SCENARIOS[0]], random_seed=1004)


def test_incomplete_grid_cannot_be_reported_as_full_d15(tmp_path):
    arguments = _materialization_fixture(tmp_path)
    path = arguments["design_root"] / "private_truth/d14b/simulation_window_manifest.json"
    document = json.loads(path.read_text())
    document["windows"][0]["ibd_scenarios"].pop()
    path.write_text(json.dumps(document))
    replay_path = arguments["replay_root"] / "private_truth/d14b_replay/simulation_window_replay.json"
    replay = json.loads(replay_path.read_text())
    replay["design"]["manifest_sha256"] = hashlib.sha256(path.read_bytes()).hexdigest()
    replay_path.write_text(json.dumps(replay))
    with pytest.raises(IbdSimulationError, match="scientific_grid_incomplete"):
        materialize_ibd(**arguments)
