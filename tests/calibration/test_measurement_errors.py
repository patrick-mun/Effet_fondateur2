"""D17 : quotas, phase, rejeu de chaque altération, contrats et isolation."""

from __future__ import annotations

from copy import deepcopy
import gzip
import hashlib
import io
import json
from pathlib import Path

import jsonschema
import numpy as np
import pytest

from effet_fondateur.calibration.materialize_measurement_errors import (
    deterministic_gzip, materialize_measurement_errors, validate_configuration,
    validate_manifest, validate_schema, load_checked_run, _audit_exports, _write_profile,
)
from effet_fondateur.calibration.report_measurement_errors import build_report
from effet_fondateur.calibration.simulation.measurement_errors import (
    MeasurementError, MeasurementPanel, derive_seed, inject_profile,
    project_dense_truth, quota, read_panel, replay_events, write_panel,
)


ROOT = Path(__file__).resolve().parents[2]


def _config() -> dict:
    return json.loads((ROOT / "calibration/configurations/d17_micro_chr5_calibration_v1.json").read_text())


def _panel() -> MeasurementPanel:
    variants = tuple(("chr5", str(10000 * index), f"CHIP_{index}", "A", "G", ".", "PASS", ".", "GT") for index in range(1, 31))
    alleles = np.array([[(0, 1), (1, 0), (index % 2, index % 2), (0, 1)]
                        for index in range(30)], dtype=np.int8)
    return MeasurementPanel(("SYN_A", "SYN_B", "SYN_C", "SYN_D"), variants, alleles)


def _inject(profile: dict, panel: MeasurementPanel | None = None, role: str = "CALIBRATION"):
    return inject_profile(panel or _panel(), profile, parent_seed=1005, panel_role=role,
                          target_chromosome=5, target_bp=150000)


def _profile(profile_id: str) -> dict:
    return next(item for item in _config()["profiles"] if item["profile_id"] == profile_id)


def test_configuration_and_schemas_are_strict():
    validate_configuration(_config())
    for path in (ROOT / "calibration/schemas").glob("measurement_error_*.schema.json"):
        jsonschema.Draft202012Validator.check_schema(json.loads(path.read_text()))
    invalid = _config()
    invalid["profiles"][0]["truth"] = "leak"
    with pytest.raises(MeasurementError, match="schema_invalid"):
        validate_configuration(invalid)


@pytest.mark.parametrize("mutation", ["seed", "control", "duplicate", "rate"])
def test_bad_configuration_is_refused(mutation):
    config = _config()
    if mutation == "seed":
        del config["measurement_error_seed"]
    elif mutation == "control":
        config["profiles"][0]["missing_rate"] = 0.5
    elif mutation == "duplicate":
        config["profiles"][1]["profile_id"] = "control"
    else:
        config["profiles"][1]["missing_rate"] = 1.1
    with pytest.raises(MeasurementError):
        validate_configuration(config)


def test_quota_half_up_and_seed_domain_separation():
    assert quota(0.5, 3) == 2
    assert quota(0.01, 120) == 1
    assert len({derive_seed(1005, role, profile, operation) for role in ("CALIBRATION", "VALIDATION")
                for profile in ("a", "b") for operation in ("missing", "density")}) == 8
    for invalid in (float("nan"), -0.1, 1.1):
        with pytest.raises(MeasurementError):
            quota(invalid, 20)


@pytest.mark.parametrize("profile", _config()["profiles"], ids=lambda item: item["profile_id"])
def test_every_profile_is_exactly_replayable_and_preserves_biological_truth(profile):
    panel = _panel()
    before = panel.alleles.copy()
    result = _inject(profile, panel)
    reconstructed, retained = replay_events(panel, result.events, profile_id=profile["profile_id"])
    assert np.array_equal(reconstructed, result.alleles)
    assert np.array_equal(retained, result.retained)
    assert np.array_equal(before, panel.alleles)
    for event in result.events:
        validate_schema(event, "measurement_error_event")
    repeated = _inject(profile, panel)
    assert repeated.events == result.events
    assert repeated.summary == result.summary
    if profile["profile_id"] == "control":
        assert not result.events
        assert np.array_equal(result.alleles, panel.alleles)


def test_exact_rates_and_density_target_policy():
    panel = _panel()
    profile = {**_profile("combined"), "missing_rate": 0.2, "genotyping_error_rate": 0.1}
    result = _inject(profile, panel)
    summary = result.summary
    assert summary["retained_marker_count"] == 15
    assert result.retained[14]
    assert summary["event_counts"]["GENOTYPE"] == 6
    assert summary["event_counts"]["MISSING_GT"] == 12
    assert np.count_nonzero(result.alleles[result.retained, :, 0] < 0) == 12
    assert all(event["position_bp"] != 150000 for event in result.events if event["error_type"] == "MARKER_DROP")
    all_missing = _inject({**_profile("missing_01"), "missing_rate": 1}, panel)
    assert np.all(all_missing.alleles == -1)


def test_profile_order_is_irrelevant_and_event_ids_are_distinct():
    profiles = _config()["profiles"]
    first = {item["profile_id"]: _inject(item) for item in profiles}
    second = {item["profile_id"]: _inject(item) for item in reversed(profiles)}
    for profile_id, result in first.items():
        assert result.events == second[profile_id].events
        assert np.array_equal(result.alleles, second[profile_id].alleles)
    event_ids = [event["error_event_id"] for result in first.values() for event in result.events]
    assert len(set(event_ids)) == len(event_ids)


def test_absent_target_is_not_added_and_empty_retained_panel_is_refused():
    panel = _panel()
    panel = MeasurementPanel(panel.samples, tuple((row[0], str(int(row[1]) + 1), *row[2:]) for row in panel.variants), panel.alleles)
    result = _inject(_profile("density_50"), panel)
    assert not result.summary["target_observed"]
    assert result.retained.sum() == 15
    with pytest.raises(MeasurementError, match="empty_retained_panel"):
        _inject({**_profile("density_50"), "marker_retention_fraction": 0.001}, panel)


@pytest.mark.parametrize("profile_id", ["phase_point_02", "phase_switch_02"])
def test_phase_never_changes_unphased_dosage_or_homozygotes(profile_id):
    panel = _panel()
    result = _inject({**_profile(profile_id), "phase_point_rate": 0.2 if profile_id == "phase_point_02" else 0,
                      "phase_switch_rate": 0.2 if profile_id == "phase_switch_02" else 0}, panel)
    assert result.events
    assert np.array_equal(np.sort(result.alleles, axis=2), np.sort(panel.alleles, axis=2))
    assert np.array_equal(result.alleles[:, 2], panel.alleles[:, 2])
    source, exported = io.StringIO(), io.StringIO()
    write_panel(source, panel, panel.alleles, result.retained, phased=False)
    write_panel(exported, panel, result.alleles, result.retained, phased=False)
    assert exported.getvalue() == source.getvalue()


def test_genotyping_error_changes_dosage_not_only_orientation():
    panel = _panel()
    result = _inject({**_profile("genotype_01"), "genotyping_error_rate": 1}, panel)
    assert np.all(result.alleles.sum(axis=2) != panel.alleles.sum(axis=2))


def test_switch_semantics_are_cumulative_and_reset_by_chromosome():
    panel = _panel()
    variants = tuple(("chr6", *row[1:]) if index >= 15 else row for index, row in enumerate(panel.variants))
    panel = MeasurementPanel(panel.samples, variants, panel.alleles)
    result = _inject({**_profile("phase_switch_02"), "phase_switch_rate": 1}, panel)
    assert result.summary["eligible_switch_boundary_count"] == 3 * (14 + 14)
    assert np.array_equal(result.alleles[0], panel.alleles[0])
    assert np.array_equal(result.alleles[15], panel.alleles[15])
    assert np.array_equal(result.alleles[1, 0], panel.alleles[1, 0, ::-1])
    assert np.array_equal(result.alleles[2], panel.alleles[2])
    assert np.array_equal(replay_events(panel, result.events, profile_id=result.summary["profile_id"])[0], result.alleles)


@pytest.mark.parametrize("mutation", ["before", "dosage", "id", "duplicate", "missing_switch"])
def test_replay_refuses_corrupt_or_incomplete_registry(mutation):
    panel = _panel()
    profile = {**_profile("phase_switch_02"), "phase_switch_rate": 0.4}
    events = deepcopy(list(_inject(profile, panel).events))
    genotype_event = next(item for item in events if item["error_type"] == "PHASE_SWITCH_GT")
    if mutation == "before":
        genotype_event["original_value"] = "0|0"
    elif mutation == "dosage":
        genotype_event["injected_value"] = "1|1"
    elif mutation == "id":
        events[0]["error_profile_id"] = "other"
    elif mutation == "duplicate":
        events.insert(1, dict(events[0]))
    else:
        events.remove(genotype_event)
    for index, event in enumerate(events, 1):
        event["error_event_id"] = f'{profile["profile_id"]}.D17_ERROR_{index:09d}'
    with pytest.raises(MeasurementError):
        replay_events(panel, events, profile_id=profile["profile_id"])


def test_no_phase_eligible_sites_and_no_truth_missingness():
    panel = _panel()
    panel.alleles[:] = 0
    result = _inject(_profile("phase_switch_02"), panel)
    assert not result.events
    panel.alleles[0, 0] = (-1, -1)
    with pytest.raises(MeasurementError, match="complete_truth_required"):
        _inject(_profile("control"), panel)


def test_projection_checks_sample_order_genotypes_and_complete_marker_set():
    panel = _panel()
    dense, observed = io.StringIO(), io.StringIO()
    retained = np.ones(30, dtype=bool)
    write_panel(dense, panel, panel.alleles, retained, phased=True)
    write_panel(observed, panel, panel.alleles, retained, phased=False)
    unphased = read_panel(observed.getvalue().splitlines(), phased=False)
    projected = project_dense_truth(dense.getvalue().splitlines(), unphased)
    assert np.array_equal(projected.alleles, panel.alleles)
    with pytest.raises(MeasurementError, match="sample_order"):
        project_dense_truth(dense.getvalue().replace("SYN_A\tSYN_B", "SYN_B\tSYN_A").splitlines(), unphased)
    with pytest.raises(MeasurementError, match="genotype_mismatch"):
        project_dense_truth(dense.getvalue().replace("0|1", "0|0", 1).splitlines(), unphased)
    with pytest.raises(MeasurementError, match="projection_incomplete"):
        project_dense_truth(dense.getvalue().splitlines()[:-1], unphased)


@pytest.mark.parametrize("bad", ["0|.", "0/1", "0|2", "0|1:50"])
def test_invalid_phased_genotypes_are_refused(bad):
    panel = _panel()
    output = io.StringIO()
    write_panel(output, panel, panel.alleles, np.ones(30, dtype=bool), phased=True)
    with pytest.raises(MeasurementError):
        read_panel(output.getvalue().replace("0|1", bad, 1).splitlines(), phased=True)


def _fixture(tmp_path: Path) -> dict:
    """Contrat D16 synthétique minimal : aucune donnée réelle ou sortie existante."""
    workspace = tmp_path / "calibration_workspace"
    source = workspace / "d16"
    private = source / "private_truth/d16"
    private.mkdir(parents=True)
    panel = _panel()
    roles = ("DENSE_PHASED_ROH_PANEL", "OBSERVED_ROH_PANEL", "ROH_INJECTION_TRUTH",
             "ROH_OBSERVABILITY_AUDIT", "ROH_INDUCED_IBD_TRUTH", "PRESERVED_IBD_TRUTH",
             "PRESERVED_PEDIGREE_TRUTH", "PRESERVED_TRANSMISSION_TRUTH")
    artifacts = []
    for index, role in enumerate(roles):
        path = private / f"artifact_{index}.gz"
        with deterministic_gzip(path) as handle:
            if index < 2:
                write_panel(handle, panel, panel.alleles, np.ones(30, dtype=bool), phased=index == 0)
            else:
                handle.write("[]\n")
        artifacts.append(dict(role=role, path=path.relative_to(source).as_posix(),
                              sha256=hashlib.sha256(path.read_bytes()).hexdigest(),
                              size_bytes=path.stat().st_size, record_count=30 if index < 2 else 0))
    document = dict(schema_version="1.0.0", manifest_state="FINALIZED", roh_simulation_manifest_id="synthetic.d16",
                    created_at="2026-10-30T00:00:00Z", method_id="within_sample_h1_to_h2_roh_copy_v1",
                    panel_role="CALIBRATION", scope="REGIONAL_TECHNICAL_INJECTION",
                    domain=dict(chromosome=5, start_bp=10000, end_bp=300000, start_cm=0, end_cm=3,
                                target_position_bp=150000, target_position_cm=1.5,
                                position_unit="BP_1_BASED_CLOSED", genetic_distance_unit="cM"),
                    configuration=dict(configuration_id="synthetic", configuration_sha256="a" * 64),
                    inputs={name: "b" * 64 for name in ("D15", "D14", "D11A")},
                    seed=dict(parent_recombination_seed=1004, derived_seed=1006,
                              derivation="SHA256(method_id:D16:parent_seed:panel_role)[0:4]_BIG_ENDIAN"),
                    truth_policy="INJECTED_AUTOZYGOUS_ORIGIN_BACKGROUND_ROH_UNKNOWN",
                    positive_roh_count=3, negative_roh_count=1, ibd_coupled_roh_count=1, induced_ibd_edge_count=1,
                    observability=[dict(scenario_id=f"synthetic_{index}", observed_marker_count=50,
                                        density_kb_per_marker=20, maximum_internal_gap_kb=30,
                                        left_flank_marker_count=2, right_flank_marker_count=2) for index in range(3)],
                    artifacts=artifacts,
                    integrity=dict(dense_variant_count=30, observed_variant_count=30, changed_h2_allele_count=0,
                                   h1_modified_count=0, roh_copy_error_count=0, ibd_copy_error_count=0,
                                   mendelian_error_count=0, negative_modified_count=0, induced_ibd_copy_error_count=0),
                    calibrated=False, gamma_inputs_modified=False, founder_effect_proven=False)
    (private / "roh_simulation_manifest.json").write_text(json.dumps(document))
    config = tmp_path / "configuration.json"
    config.write_text(json.dumps(_config()))
    return dict(d16_root=source, config_path=config, workspace_root=workspace,
                run_id="d17_test", created_at="2026-10-30T00:00:00Z")


def test_materialization_is_private_deterministic_and_fully_audited(tmp_path, monkeypatch):
    monkeypatch.setattr(gzip.time, "time", lambda: 1000000000)
    arguments = _fixture(tmp_path)
    manifest = materialize_measurement_errors(**arguments)
    assert len(manifest["profiles"]) == 8
    assert len(manifest["artifacts"]) == 26
    validate_manifest(manifest)
    assert all(item["audit"]["status"] == "PASS" for item in manifest["profiles"])
    assert not manifest["biological_truth_modified"]
    for item in manifest["artifacts"]:
        path = arguments["workspace_root"] / "d17_test" / item["path"]
        assert path.stat().st_mode & 0o777 == 0o600
        assert hashlib.sha256(path.read_bytes()).hexdigest() == item["sha256"]
    with pytest.raises(MeasurementError, match="run_exists"):
        materialize_measurement_errors(**arguments)
    monkeypatch.setattr(gzip.time, "time", lambda: 1000000002)
    assert materialize_measurement_errors(**_fixture(tmp_path / "second")) == manifest


@pytest.mark.parametrize("mutation", ["hash", "role", "symlink", "hardlink", "real", "traversal"])
def test_input_isolation_and_tampering(tmp_path, mutation):
    arguments = _fixture(tmp_path / "data/input" if mutation == "real" else tmp_path)
    source = arguments["d16_root"] / "private_truth/d16/artifact_0.gz"
    if mutation == "hash":
        source.write_bytes(b"tampered")
    elif mutation == "role":
        config = _config()
        config["panel_role"] = "VALIDATION"
        arguments["config_path"].write_text(json.dumps(config))
    elif mutation == "symlink":
        source.rename(source.with_suffix(".backup"))
        source.symlink_to(source.with_suffix(".backup").name)
    elif mutation == "hardlink":
        (source.parent / "alias.gz").hardlink_to(source)
    elif mutation == "traversal":
        arguments["run_id"] = "../escape"
    with pytest.raises(MeasurementError):
        materialize_measurement_errors(**arguments)
    assert not (arguments["workspace_root"] / "d17_test").exists()


def test_failed_export_never_publishes_partial_run(tmp_path, monkeypatch):
    arguments = _fixture(tmp_path)
    def fail(*args, **kwargs):
        raise MeasurementError("synthetic_export_failure")
    monkeypatch.setattr("effet_fondateur.calibration.materialize_measurement_errors._audit_exports", fail)
    with pytest.raises(MeasurementError, match="synthetic_export_failure"):
        materialize_measurement_errors(**arguments)
    assert not (arguments["workspace_root"] / "d17_test").exists()
    assert not list(arguments["workspace_root"].glob(".d17_test.*"))


@pytest.mark.parametrize("filename", ["observed.vcf.gz", "candidate.phased.vcf.gz", "measurement_errors.jsonl.gz"])
def test_independent_audit_detects_unregistered_export_changes(tmp_path, filename):
    panel = _panel()
    result = _inject({**_profile("missing_01"), "missing_rate": 0.2}, panel)
    root = tmp_path / "profile"
    _write_profile(root, panel, result)
    path = root / filename
    with gzip.open(path, "rt") as handle:
        payload = handle.read()
    if filename.endswith("jsonl.gz"):
        payload = "\n".join(payload.splitlines()[1:]) + "\n"
    elif filename == "observed.vcf.gz":
        payload = payload.replace("0/1", "1/1", 1)
    else:
        payload = payload.replace("0|1", "1|1", 1)
    with deterministic_gzip(path) as handle:
        handle.write(payload)
    with pytest.raises(MeasurementError):
        _audit_exports(panel, result, root)


def test_report_is_aggregate_complete_and_rejects_artifact_tampering(tmp_path):
    arguments = _fixture(tmp_path)
    manifest = materialize_measurement_errors(**arguments)
    root = arguments["workspace_root"] / "d17_test"
    template = tmp_path / "template.html"
    template.write_text('<!doctype html><html><head><title>Template</title></head><body></body></html>')
    output = build_report(run_root=root, template_path=template,
                          output_root=arguments["workspace_root"] / "report")
    html = output.read_text()
    assert html.count('<section id="s') == 14
    assert html.count('<a class="toc-entry"') == 14
    assert "NOT_EVALUABLE" in html and "calibrated=false" in html
    assert "Figure 1" in html and "Figure 2" in html and "Tableau 2" in html
    assert all(sample not in html for sample in _panel().samples)
    assert "SYN_A" not in (output.parent / "report_data.json").read_text()
    with pytest.raises(MeasurementError, match="report_output_invalid"):
        build_report(run_root=root, template_path=template, output_root=output.parent)
    artifact = root / manifest["artifacts"][0]["path"]
    artifact.write_text("tampered")
    with pytest.raises(MeasurementError, match="artifact_hash"):
        load_checked_run(root)


@pytest.mark.parametrize("mutation", ["seed", "quota", "artifact", "density", "denominator", "record_count", "control"])
def test_manifest_rejects_inconsistent_aggregates(tmp_path, mutation):
    manifest = materialize_measurement_errors(**_fixture(tmp_path))
    item = manifest["profiles"][-1]
    if mutation == "seed":
        item["seeds"]["genotype"] += 1
    elif mutation == "quota":
        item["event_counts"]["MISSING_GT"] += 1
    elif mutation == "artifact":
        manifest["artifacts"].pop()
    elif mutation == "density":
        item["requested"]["marker_retention_fraction"] = 0.9
    elif mutation == "denominator":
        item["retained_cell_count"] += 1
    elif mutation == "record_count":
        manifest["artifacts"][-1]["record_count"] += 1
    else:
        manifest["profiles"][0]["requested"]["profile_id"] = "not_control"
    with pytest.raises(MeasurementError):
        validate_manifest(manifest)
