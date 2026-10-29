"""D16 : autozygotie connue, conservation H1/IBD et isolation des résultats."""

from __future__ import annotations

from dataclasses import asdict, replace
import gzip
import hashlib
import io
import json
from pathlib import Path

import jsonschema
import pytest

from effet_fondateur.calibration.materialize_roh import materialize_roh
from effet_fondateur.calibration.report_roh_injection import build_report
from effet_fondateur.calibration.simulation.founder_recombination import GeneticMapPoint
from effet_fondateur.calibration.simulation.ibd import build_ibd_injections
from effet_fondateur.calibration.simulation.roh import (
    RohSimulationError, build_roh_injections, derive_roh_ibd_edges,
    transform_roh_vcf, validate_roh_density,
)


ROOT = Path(__file__).resolve().parents[2]
SAMPLES = [f"SYN_CAL_{index:03d}" for index in range(35)]
POINTS = (GeneticMapPoint(1000000, 0), GeneticMapPoint(6000000, 5), GeneticMapPoint(11000000, 10))
POSITIONS = tuple(range(1000000, 11000001, 10000))
IBD_SCENARIOS = [{"scenario_id": name, "segment_length_cm": length,
                  "minimum_observed_markers": 1} for name, length in (
                      ("one", 1), ("one_half", 1.5), ("primary_2cm", 2),
                      ("three", 3), ("five", 5))]


def _ibd_truth() -> list[dict[str, object]]:
    injections = build_ibd_injections(eligible_samples=SAMPLES[:25], map_points=POINTS,
                                     chromosome=5, target_bp=6000000,
                                     scenarios=IBD_SCENARIOS, random_seed=1004)
    return [{"replication_id": "d15_test", **asdict(item)} for item in injections]


def _config() -> dict[str, object]:
    document = json.loads((ROOT / "calibration/configurations/d16_micro_chr5_calibration_v1.json").read_text())
    for item in document["scenarios"]:
        if item["relationship_to_target"] == "TARGET_OVERLAP":
            item["start_bp"], item["end_bp"] = 5000000, 6999999
        elif item["relationship_to_target"] == "TARGET_DISTINCT":
            item["start_bp"], item["end_bp"] = 8000000, 9999999
    return document


def _plan(seed: int = 1006):
    return build_roh_injections(scenarios=_config()["scenarios"],
                                eligible_samples=SAMPLES[25:33], ibd_truth=_ibd_truth(),
                                map_points=POINTS, chromosome=5, target_bp=6000000,
                                random_seed=seed)


def _vcf() -> str:
    header = "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t" + "\t".join(SAMPLES) + "\n"
    return header + "".join(
        f"chr5\t{position}\tv{position}\tA\tG\t.\tPASS\t.\tGT\t" +
        "\t".join(f"0|{(index + position // 10000) % 2}" for index in range(35)) + "\n"
        for position in POSITIONS
    )


def _keys():
    return {("chr5", position, "A", "G"): f"CHIP_{position}" for position in POSITIONS}


def _rows(text: str):
    return {int(line.split("\t")[1]): line.split("\t")[9:]
            for line in text.splitlines() if not line.startswith("#")}


def test_plan_deterministic_four_distinct_contexts():
    plan = _plan()
    assert plan == _plan()
    assert plan != _plan(1007)
    assert len({item.sample_id for item in plan}) == 4
    assert sum(item.is_positive for item in plan) == 3
    negative = next(item for item in plan if not item.is_positive)
    assert negative.start_bp is None and negative.end_bp is None
    assert negative.source_haplotype_index is None
    coupled = next(item for item in plan if item.source_ibd_scenario_id)
    source = next(item for item in _ibd_truth() if item["scenario_id"] == coupled.source_ibd_scenario_id)
    assert coupled.sample_id == source["sample_1_id"]


def test_self_copy_preserves_h1_ibd_and_every_unmodified_genotype():
    dense, observed = io.StringIO(), io.StringIO()
    plan = _plan()
    stats = transform_roh_vcf(io.StringIO(_vcf()), dense, observed, injections=plan,
                              observed_keys=_keys(), transmissions=(), ibd_truth=_ibd_truth())
    before, after = _rows(_vcf()), _rows(dense.getvalue())
    by_sample = {item.sample_id: item for item in plan if item.is_positive}
    for position, genotypes in after.items():
        for index, sample in enumerate(SAMPLES):
            original, result = before[position][index], genotypes[index]
            assert original[0] == result[0]
            item = by_sample.get(sample)
            if item and item.start_bp <= position <= item.end_bp:
                assert result[0] == result[2]
            else:
                assert result == original
    assert stats["dense_variant_count"] == stats["observed_variant_count"] == 1001
    assert stats["changed_h2_allele_count"] > 0
    assert "roh_" not in observed.getvalue()
    assert "|" not in observed.getvalue()


def test_density_is_a_gate_not_an_adjustment():
    audit = validate_roh_density(_plan(), POSITIONS, minimum_markers=50,
                                maximum_density_kb=50, maximum_gap_kb=1000, minimum_flanks=2)
    assert all(item["observed_marker_count"] == 200 for item in audit)
    with pytest.raises(RohSimulationError, match="marker_density_insufficient"):
        validate_roh_density(_plan(), POSITIONS[::100], minimum_markers=50,
                             maximum_density_kb=50, maximum_gap_kb=1000, minimum_flanks=2)


def test_coupled_roh_records_additional_ibd_copy_on_h2():
    plan, truth = _plan(), _ibd_truth()
    edges = derive_roh_ibd_edges(plan, truth, POINTS, replication_id="test")
    assert len(edges) == 1
    edge = edges[0]
    coupled = next(item for item in plan if item.source_ibd_scenario_id)
    source = next(item for item in truth if item["scenario_id"] == coupled.source_ibd_scenario_id)
    assert edge["start_bp"] == max(coupled.start_bp, source["start_bp"])
    assert edge["end_bp"] == min(coupled.end_bp, source["end_bp"])
    assert (edge["haplotype_1_index"], edge["haplotype_2_index"]) == (2, 1)
    dense = io.StringIO()
    transform_roh_vcf(io.StringIO(_vcf()), dense, io.StringIO(), injections=plan,
                      observed_keys=_keys(), transmissions=(), ibd_truth=truth)
    rows = _rows(dense.getvalue())
    first, second = SAMPLES.index(edge["sample_1_id"]), SAMPLES.index(edge["sample_2_id"])
    for position, genotypes in rows.items():
        if edge["start_bp"] <= position <= edge["end_bp"]:
            assert genotypes[first][2] == genotypes[second][0]


def test_incorrect_target_relationship_is_refused():
    config = _config()
    config["scenarios"][1]["relationship_to_target"] = "TARGET_OVERLAP"
    with pytest.raises(RohSimulationError, match="target_relationship_mismatch"):
        build_roh_injections(scenarios=config["scenarios"], eligible_samples=SAMPLES[25:33],
                             ibd_truth=_ibd_truth(), map_points=POINTS,
                             chromosome=5, target_bp=6000000, random_seed=1006)


def test_negative_ibd_cannot_be_claimed_as_positive_context():
    config = _config()
    config["scenarios"][-1]["source_ibd_scenario_id"] = "primary_2cm_negative"
    with pytest.raises(RohSimulationError, match="positive_ibd_context_required"):
        build_roh_injections(scenarios=config["scenarios"], eligible_samples=SAMPLES[25:33],
                             ibd_truth=_ibd_truth(), map_points=POINTS,
                             chromosome=5, target_bp=6000000, random_seed=1006)


def test_transmission_subject_cannot_be_modified():
    with pytest.raises(RohSimulationError, match="pedigree_modification_forbidden"):
        transform_roh_vcf(io.StringIO(_vcf()), io.StringIO(), io.StringIO(),
                          injections=_plan(), observed_keys=_keys(),
                          transmissions=((_plan()[1].sample_id, SAMPLES[-1], 1, 1),),
                          ibd_truth=_ibd_truth())


def test_d15_negative_subject_is_protected():
    plan = list(_plan())
    negative = next(item for item in _ibd_truth() if not item["is_positive"])
    plan[1] = replace(plan[1], sample_id=negative["sample_1_id"])
    with pytest.raises(RohSimulationError, match="negative_or_template_ibd_modification"):
        transform_roh_vcf(io.StringIO(_vcf()), io.StringIO(), io.StringIO(),
                          injections=plan, observed_keys=_keys(), transmissions=(),
                          ibd_truth=_ibd_truth())


def test_corrupted_ibd_is_not_silently_accepted():
    rows = _vcf().splitlines()
    source = next(item for item in _ibd_truth() if item["is_positive"])
    index = SAMPLES.index(source["sample_1_id"]) + 9
    for row_index, line in enumerate(rows):
        if line.startswith("chr5\t6000000\t"):
            fields = line.split("\t")
            fields[index] = "1|0"
            rows[row_index] = "\t".join(fields)
    with pytest.raises(RohSimulationError, match="existing_ibd_copy_error"):
        transform_roh_vcf(io.StringIO("\n".join(rows)), io.StringIO(), io.StringIO(),
                          injections=_plan(), observed_keys=_keys(), transmissions=(),
                          ibd_truth=_ibd_truth())


def _fixture(tmp_path: Path) -> dict[str, object]:
    workspace = tmp_path / "calibration_workspace"
    d15, d14, d11 = (workspace / name for name in ("d15", "d14", "d11"))
    def sha(path):
        return hashlib.sha256(path.read_bytes()).hexdigest()
    def write(root, relative, payload):
        path = root / relative
        path.parent.mkdir(parents=True, exist_ok=True)
        if path.suffix == ".gz":
            # Les empreintes des entrées doivent rester indépendantes de
            # l'heure et du chemin des deux jeux comparés pour le rejeu.
            with path.open("wb") as raw:
                with gzip.GzipFile(filename="", fileobj=raw, mode="wb", mtime=0) as handle:
                    handle.write(payload.encode("utf-8"))
        else:
            path.write_text(payload)
        return path
    def artifact(root, relative, role, payload, count):
        path = write(root, relative, payload)
        return dict(role=role, path=relative, sha256=sha(path), size_bytes=path.stat().st_size, record_count=count)
    domain = dict(chromosome=5, start_bp=1000000, end_bp=11000000, start_cm=0,
                  end_cm=10, target_position_bp=6000000, target_position_cm=5,
                  position_unit="BP_1_BASED_CLOSED", genetic_distance_unit="cM")
    pedigree = "sample_id\tfamily_unit_id\n" + "".join(f"{sample}\tUNRELATED_{sample}\n" for sample in SAMPLES)
    d14_artifacts = [artifact(d14, "private_truth/d14/pedigree.gz", "PEDIGREE_TRUTH", pedigree, 35),
                     artifact(d14, "private_truth/d14/carriers.gz", "FOUNDER_CARRIER_TRUTH", "sample_id\n", 0),
                     artifact(d14, "private_truth/d14/transmissions.gz", "MENDELIAN_TRANSMISSION_TRUTH", "child_id\tparent_id\tparental_haplotype_index\tchild_haplotype_index\n", 0)]
    d14_doc = dict(manifest_state="FINALIZED", domain=domain,
                   panel={"role": "CALIBRATION", "individual_count": 35}, artifacts=d14_artifacts)
    d14_path = write(d14, "private_truth/d14/founder_simulation_manifest.json", json.dumps(d14_doc))
    map_artifact = artifact(d11, "private_truth/d11/map.tsv", "GENETIC_MAP",
                            "CHROMOSOME\tPOSITION_BP\tPOSITION_CM\n5\t1000000\t0\n5\t6000000\t5\n5\t11000000\t10\n", 3)
    d11_path = write(d11, "private_truth/d11/base_haplotype_manifest.json", json.dumps(dict(
        manifest_state="FINALIZED", domain=domain, panel={"role": "CALIBRATION"}, artifacts=[map_artifact])))
    source_manifest = json.loads((ROOT / "calibration/schemas/ibd_simulation_manifest.schema.json").read_text())
    d15_doc = dict(schema_version="1.0.0", manifest_state="FINALIZED", ibd_simulation_manifest_id="fixture.d15",
                   created_at="2026-10-29T10:00:00Z", method_id="disjoint_unrelated_haplotype_copy_v1",
                   panel_role="CALIBRATION", scope="REGIONAL_TECHNICAL_INJECTION", domain=domain,
                   inputs={"D14": sha(d14_path), "D11A": sha(d11_path), "DESIGN": "a" * 64, "REPLAY": "b" * 64},
                   seed={"parent_recombination_seed": 1004, "derived_seed": 1005,
                         "derivation": "SHA256(method_id:D15:parent_seed:panel_role)[0:4]_BIG_ENDIAN"},
                   truth_policy="INJECTED_ORIGIN_INTERVAL_BACKGROUND_IBD_UNKNOWN", positive_pair_count=5,
                   negative_pair_count=5, calibrated=False, founder_effect_proven=False,
                   artifacts=[artifact(d15, "private_truth/d15/dense.gz", "DENSE_PHASED_IBD_PANEL", _vcf(), 1001),
                              artifact(d15, "private_truth/d15/observed.gz", "OBSERVED_IBD_PANEL", _vcf(), 1001),
                              artifact(d15, "private_truth/d15/truth.json", "IBD_INJECTION_TRUTH", json.dumps(_ibd_truth()), 10)],
                   integrity=dict(dense_variant_count=1001, observed_variant_count=1001,
                                  changed_haplotype_allele_count=0, mendelian_error_count=0,
                                  positive_copy_error_count=0, negative_modified_count=0))
    jsonschema.validate(d15_doc, source_manifest)
    write(d15, "private_truth/d15/ibd_simulation_manifest.json", json.dumps(d15_doc))
    config = tmp_path / "d16.json"
    config.write_text(json.dumps(_config()))
    return dict(d15_root=d15, d14_root=d14, d11a_root=d11, config_path=config,
                workspace_root=workspace, run_id="d16_test", created_at="2026-10-29T10:00:00Z")


def test_materialization_validates_and_preserves_ibd_truth(tmp_path, monkeypatch):
    monkeypatch.setattr(gzip.time, "time", lambda: 1000000000)
    arguments = _fixture(tmp_path)
    manifest = materialize_roh(**arguments)
    assert manifest["positive_roh_count"] == 3
    assert manifest["negative_roh_count"] == manifest["ibd_coupled_roh_count"] == 1
    assert manifest["induced_ibd_edge_count"] == 1
    assert manifest["integrity"]["mendelian_error_count"] == 0
    truth_item = next(item for item in manifest["artifacts"] if item["role"] == "PRESERVED_IBD_TRUTH")
    before = arguments["d15_root"] / "private_truth/d15/truth.json"
    after = arguments["workspace_root"] / "d16_test" / truth_item["path"]
    assert before.read_bytes() == after.read_bytes()
    assert after.stat().st_mode & 0o777 == 0o600
    with pytest.raises(RohSimulationError, match="run_exists"):
        materialize_roh(**arguments)
    monkeypatch.setattr(gzip.time, "time", lambda: 1000000002)
    second = _fixture(tmp_path / "second")
    assert materialize_roh(**second) == manifest


def test_input_hash_mismatch_blocks_output(tmp_path):
    arguments = _fixture(tmp_path)
    path = arguments["d14_root"] / "private_truth/d14/founder_simulation_manifest.json"
    path.write_text(path.read_text() + " ")
    with pytest.raises(RohSimulationError, match="lineage_binding_mismatch"):
        materialize_roh(**arguments)
    assert not (arguments["workspace_root"] / "d16_test").exists()


@pytest.mark.parametrize("genotype", ["0/1", "0|.", ".|."])
def test_invalid_or_missing_gt_blocks_injection(genotype):
    with pytest.raises(RohSimulationError, match="gt_invalid"):
        transform_roh_vcf(io.StringIO(_vcf().replace("0|1", genotype)),
                          io.StringIO(), io.StringIO(), injections=_plan(),
                          observed_keys=_keys(), transmissions=(), ibd_truth=_ibd_truth())


def test_no_extrapolation_or_silent_window_expansion():
    config = _config()
    config["scenarios"][1]["end_bp"] = 12000000
    with pytest.raises(RohSimulationError, match="interval_outside_domain"):
        build_roh_injections(scenarios=config["scenarios"], eligible_samples=SAMPLES[25:33],
                             ibd_truth=_ibd_truth(), map_points=POINTS, chromosome=5,
                             target_bp=6000000, random_seed=1006)


def test_negative_truth_cannot_have_roh_bounds():
    schema = json.loads((ROOT / "calibration/schemas/roh_injection_truth.schema.json").read_text())
    truth = [{"replication_id": "test", **asdict(item)} for item in _plan()]
    jsonschema.validate(truth, schema)
    truth[0]["start_bp"] = 1000000
    with pytest.raises(jsonschema.ValidationError):
        jsonschema.validate(truth, schema)


def test_real_data_path_is_refused_before_reading(tmp_path):
    arguments = _fixture(tmp_path / "data/input")
    with pytest.raises(RohSimulationError, match="input_integrity_error"):
        materialize_roh(**arguments)


def test_corrupted_source_artifact_is_rejected(tmp_path):
    arguments = _fixture(tmp_path)
    path = arguments["d15_root"] / "private_truth/d15/dense.gz"
    path.write_bytes(b"corrupted")
    with pytest.raises(RohSimulationError, match="input_integrity_error"):
        materialize_roh(**arguments)
    assert not (arguments["workspace_root"] / "d16_test").exists()


def test_configuration_rejects_unknown_fields(tmp_path):
    arguments = _fixture(tmp_path)
    config = _config()
    config["unexpected"] = True
    arguments["config_path"].write_text(json.dumps(config))
    with pytest.raises(RohSimulationError, match="schema_invalid:roh_simulation_config"):
        materialize_roh(**arguments)


def test_independent_report_contains_limits_but_no_sample_identifiers(tmp_path):
    arguments = _fixture(tmp_path)
    materialize_roh(**arguments)
    template = tmp_path / "template.html"
    template.write_text('<!doctype html><html lang="fr"><head><title>Modèle</title></head><body></body></html>')
    report = build_report(run_root=arguments["workspace_root"] / "d16_test",
                          template_path=template,
                          output_root=arguments["workspace_root"] / "report")
    html = report.read_text()
    assert all(sample not in html for sample in SAMPLES)
    assert html.count('<section id="s') == 14
    assert "NOT_EVALUABLE" in html
    assert "gamma_inputs_modified=false" in html
    assert "Tableau 1" in html and "Figure 2" in html
    assert html.index('class="page-footer"') < html.index('<main>')
