"""Contrats D23a : aucun outil génétique réel et seulement tmp_path."""

import json
from pathlib import Path

import numpy as np
import pytest
import yaml

from effet_fondateur.calibration.phasing_preparation import (
    _public_path, _source_rows, build_chromosome, choose_founders, evenly_spaced,
    inventory, load_configuration, transmission_track, verify_preparation,
    finalize_preparation,
)


CONFIG_PATH = Path(__file__).resolve().parents[2] / "calibration/configurations/d23a_phasing_panel_v2.json"


def test_prespecified_configuration_is_strict(tmp_path: Path) -> None:
    config = load_configuration(CONFIG_PATH)
    assert 3 * config["trio_count"] + config["unrelated_count"] == 40
    for change in ({"unrelated_count": True}, {"founder_seed": -1}, {"rare_maf_upper_bound": float("nan")}, {"unexpected": 1}, {"phasing_seed": config["meiosis_seed"]}):
        path = tmp_path / "invalid.json"
        path.write_text(json.dumps(config | change))
        with pytest.raises(ValueError, match="configuration_invalid"):
            load_configuration(path)


def test_meioses_are_deterministic_and_domain_separated() -> None:
    positions = np.linspace(0, 300, 1000)
    first, events = transmission_track(positions, 2302, 1, 0, "father")
    repeated, repeated_events = transmission_track(positions, 2302, 1, 0, "father")
    other, _ = transmission_track(positions, 2302, 1, 0, "mother")
    assert np.array_equal(first, repeated) and events == repeated_events
    assert not np.array_equal(first, other)
    assert set(np.unique(first)) <= {0, 1}
    assert events == sorted(events) and all(0 <= value <= 300 for value in events)


@pytest.mark.parametrize("positions", [np.array([1.0]), np.array([2.0, 1.0]), np.array([0.0, float("nan")])])
def test_invalid_map_is_blocking(positions: np.ndarray) -> None:
    with pytest.raises(ValueError, match="map_invalid"):
        transmission_track(positions, 1, 1, 0, "father")


def test_flat_map_has_no_spurious_recombination() -> None:
    track, events = transmission_track(np.zeros(20), 2302, 1, 0, "mother")
    assert events == [] and len(set(track)) == 1


def test_children_copy_the_two_actual_parental_tracks_without_mutating_founders() -> None:
    config = load_configuration(CONFIG_PATH) | dict(trio_count=1, unrelated_count=1)
    positions = np.linspace(0, 200, 200)
    founders = np.zeros((200, 3, 2), dtype=np.int8)
    founders[:, :, 1] = 1
    original = founders.copy()
    panel, records = build_chromosome(founders, positions, 1, config)
    assert np.array_equal(founders, original) and np.array_equal(panel[:, :3], original)
    for side, name in enumerate(("father", "mother")):
        track, _ = transmission_track(positions, config["meiosis_seed"], 1, 0, name)
        assert np.array_equal(panel[:, 3, side], founders[np.arange(200), side, track])
    assert len(records) == 2


def test_founders_never_reuse_a_donor_independence_group() -> None:
    config = load_configuration(CONFIG_PATH) | dict(trio_count=1, unrelated_count=1)
    members = [f"D{index:03d}" for index in range(24)]
    groups = {sample: f"G{index // 2}" for index, sample in enumerate(members)}
    populations = {sample: "STRATUM" for sample in members}
    first = choose_founders(members, groups, populations, config)
    assert first == choose_founders(list(reversed(members)), groups, populations, config)
    assert len(first) == 3
    assert len({groups[origin[0]] for pair in first for origin in pair}) == 6
    with pytest.raises(ValueError, match="insufficient"):
        choose_founders(members[:4], groups, populations, config)


def test_marker_selection_is_spaced_unique_and_keeps_ends() -> None:
    selected = evenly_spaced(list(range(100)), 10)
    assert len(selected) == len(set(selected)) == 10
    assert selected[0] == 0 and selected[-1] == 99
    assert evenly_spaced([1, 2], 3) == [1, 2]


def test_source_parser_refuses_partial_unphased_and_ambiguous_records() -> None:
    records = "chr1\t10\tv1\tA\tG\t0|1\t1|1\nchr1\t20\tv2\tA\tG\t0/1\t1|1\n"
    variants, calls = _source_rows(records, 1, 2)
    assert len(variants) == 2 and set(calls) == {10}
    variants, calls = _source_rows(records + "chr1\t10\tv3\tG\tA\t1|0\t0|0\n", 1, 2)
    assert 10 not in calls


def test_public_path_forbids_real_inputs_caches_and_symlinks(tmp_path: Path) -> None:
    for suffix in ("data/input/source.vcf", "data/cache/source.vcf", "dock6_real/source.vcf"):
        path = tmp_path / suffix
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text("synthetic")
        with pytest.raises(ValueError, match="forbidden"):
            _public_path(path)
    target = tmp_path / "public.vcf"
    target.write_text("synthetic")
    link = tmp_path / "linked.vcf"
    link.symlink_to(target)
    with pytest.raises(ValueError, match="forbidden"):
        _public_path(link)


def test_closed_inventory_detects_any_added_or_modified_file(tmp_path: Path, monkeypatch) -> None:
    from effet_fondateur.calibration import phasing_preparation
    monkeypatch.setattr(phasing_preparation, "_validate_document", lambda *args: None)
    path = tmp_path / "result.txt"
    path.write_text("synthetic")
    document = dict(preparation_state="FINALIZED", phasing_executed=False, artifacts=inventory(tmp_path))
    (tmp_path / "preparation_manifest.json").write_text(json.dumps(document))
    assert verify_preparation(tmp_path, require_blinded_seeds=False) == document
    path.write_text("modified")
    with pytest.raises(ValueError, match="inventory_invalid"):
        verify_preparation(tmp_path, require_blinded_seeds=False)


def _draft_preparation(tmp_path: Path, monkeypatch) -> tuple[Path, dict]:
    from effet_fondateur.calibration import phasing_preparation
    from effet_fondateur.calibration.export_blinded_inputs import pipeline_configuration
    workspace = tmp_path / "calibration_workspace"
    source = workspace / "draft"
    (source / "private_truth").mkdir(parents=True)
    (source / "blinded_inputs").mkdir()
    config = load_configuration(CONFIG_PATH)
    config.pop("phasing_seed")
    config["configuration_id"] = "draft"
    (source / "private_truth/configuration.json").write_text(json.dumps(config))
    (source / "unchanged.txt").write_text("synthetic immutable genotype placeholder")
    target = dict(chromosome=5, position_bp=100, ref="A", alt="G")
    (source / "blinded_inputs/pipeline.yaml").write_text(yaml.safe_dump(pipeline_configuration(target)))
    (source / "blinded_inputs/phasing_parameters.yaml").write_text(yaml.safe_dump({
        name: dict(shapeit5_seed=config["meiosis_seed"], require_reference_sample_selection=True)
        for name in ("phase_target_region", "phase_autosomal_panel")}))
    document = dict(configuration=config, audit=dict(genotype_mismatch_count=0), artifacts=inventory(source))
    (source / "preparation_manifest.json").write_text(json.dumps(document))
    monkeypatch.setattr(phasing_preparation, "verify_preparation", lambda root, **kwargs: json.loads((root / "preparation_manifest.json").read_text()))
    monkeypatch.setattr(phasing_preparation, "audit_exports", lambda root, configuration: document["audit"])
    return source, document


def test_seed_revision_preserves_source_and_every_genetic_artifact(tmp_path: Path, monkeypatch) -> None:
    source, document = _draft_preparation(tmp_path, monkeypatch)
    source_bytes = {path.relative_to(source).as_posix(): path.read_bytes() for path in source.rglob("*") if path.is_file()}
    result = finalize_preparation(source_root=source, configuration=CONFIG_PATH, workspace_root=source.parent, run_id="final")
    target = source.parent / "final"
    parameters = yaml.safe_load((target / "blinded_inputs/phasing_parameters.yaml").read_text())
    assert all(row["shapeit5_seed"] == 2303 and row["shapeit5_seed"] != document["configuration"]["meiosis_seed"] for row in parameters.values())
    assert (source / "unchanged.txt").read_bytes() == (target / "unchanged.txt").read_bytes()
    assert result["source_preparation_sha256"] and result["implementation_sha256"]
    assert source_bytes == {path.relative_to(source).as_posix(): path.read_bytes() for path in source.rglob("*") if path.is_file()}


def test_seed_revision_refuses_any_scientific_parameter_change(tmp_path: Path, monkeypatch) -> None:
    source, _ = _draft_preparation(tmp_path, monkeypatch)
    configuration = load_configuration(CONFIG_PATH) | dict(founder_seed=999)
    path = tmp_path / "changed.json"
    path.write_text(json.dumps(configuration))
    with pytest.raises(ValueError, match="scientific_configuration_change_forbidden"):
        finalize_preparation(source_root=source, configuration=path, workspace_root=source.parent, run_id="final")
    assert not (source.parent / "final").exists()


def test_seed_revision_refuses_a_changed_genetic_file(tmp_path: Path, monkeypatch) -> None:
    from effet_fondateur.calibration import phasing_preparation
    source, _ = _draft_preparation(tmp_path, monkeypatch)
    copytree = phasing_preparation.shutil.copytree
    def altered_copy(origin: Path, target: Path, *args, **kwargs):
        result = copytree(origin, target, *args, **kwargs)
        if Path(origin) == source:
            (Path(target) / "unchanged.txt").write_text("altered")
        return result
    monkeypatch.setattr(phasing_preparation.shutil, "copytree", altered_copy)
    with pytest.raises(ValueError, match="genetic_artifact_changed"):
        finalize_preparation(source_root=source, configuration=CONFIG_PATH, workspace_root=source.parent, run_id="final")
    assert not (source.parent / "final").exists()
