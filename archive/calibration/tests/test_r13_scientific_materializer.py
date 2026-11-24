from __future__ import annotations

import gzip
from pathlib import Path
import shutil

from effet_fondateur.calibration.contracts import load_json
from effet_fondateur.calibration.r13_scientific_materializer import (
    SmokeDimensions,
    dense_variant_allocation,
    load_materializer_configuration,
    materialize_ibd_world,
    materialize_shapeit_shard,
    materializer_contract_audit,
)
from effet_fondateur.calibration.r13_scientific_runner import (
    build_external_call_schedule,
)


REPOSITORY_ROOT = Path(__file__).resolve().parents[2]
RELEASE_ROOT = (
    REPOSITORY_ROOT
    / "calibration_workspace/releases/panel_master_generic_v1_r12c_scientific_validation"
)


def _fake_compress_and_index(source: Path, output: Path, _bcftools: Path) -> None:
    shutil.copyfile(source, output)
    Path(f"{output}.tbi").write_text("synthetic-index\n", encoding="utf-8")


def _variant_lines(path: Path) -> list[str]:
    with gzip.open(path, "rt", encoding="utf-8") as handle:
        return [line for line in handle if not line.startswith("#")]


def test_materializer_contract_covers_all_calls_and_nested_density_totals() -> None:
    configuration = load_materializer_configuration(REPOSITORY_ROOT)
    allocation = dense_variant_allocation(configuration)
    audit = materializer_contract_audit(REPOSITORY_ROOT)

    assert sum(allocation.values()) == 1_200_000
    assert audit["state"] == "FULL_CALL_MATERIALIZATION_READY"
    assert audit["supported_call_count"] == 1_180
    assert audit["observed_variant_totals"] == {
        "D_MAX": 600_000,
        "D_MID": 300_000,
        "D_LOW": 100_000,
    }


def test_shapeit_materializer_separates_blind_input_and_phase_truth(
    tmp_path: Path,
    monkeypatch,
) -> None:
    monkeypatch.setattr(
        "effet_fondateur.calibration.r13_scientific_materializer._compress_and_index_vcf",
        _fake_compress_and_index,
    )
    run_root = tmp_path / "run"
    truth_root = tmp_path / "truth"
    run_root.mkdir()
    call = build_external_call_schedule()["shapeit"][0]

    manifest = materialize_shapeit_shard(
        REPOSITORY_ROOT,
        RELEASE_ROOT,
        run_root,
        truth_root,
        call,
        test_dimensions=SmokeDimensions(12, 9, 2, 120),
    )
    repeated = materialize_shapeit_shard(
        REPOSITORY_ROOT,
        RELEASE_ROOT,
        run_root,
        truth_root,
        call,
        test_dimensions=SmokeDimensions(12, 9, 2, 120),
    )

    assert manifest == repeated
    assert manifest["variant_count"] == 60
    assert manifest["study_count"] == 12
    assert manifest["reference_count"] == 9
    assert manifest["complete_trios"] == 2
    assert manifest["truth_available_to_worker"] is False
    assert "truth" not in str(manifest).lower().replace(
        "truth_available_to_worker", ""
    ).replace("source_truth_binding_sha256", "")
    blind_vcf = run_root / manifest["artifacts"][0]["relative_path"]
    private_manifest = load_json(
        truth_root / "shapeit/chr01/D_MAX/manifest.json"
    )
    truth_vcf = truth_root / private_manifest["artifacts"][0]["relative_path"]
    assert len(_variant_lines(blind_vcf)) == 60
    assert "/" in _variant_lines(blind_vcf)[0].split("\t")[9]
    assert "|" in _variant_lines(truth_vcf)[0].split("\t")[9]


def test_ibd_materializer_creates_100_pair_world_without_truth_leak(
    tmp_path: Path,
    monkeypatch,
) -> None:
    monkeypatch.setattr(
        "effet_fondateur.calibration.r13_scientific_materializer._compress_and_index_vcf",
        _fake_compress_and_index,
    )
    run_root = tmp_path / "run"
    truth_root = tmp_path / "truth"
    run_root.mkdir()
    call = build_external_call_schedule()["ibd"][0]

    manifest = materialize_ibd_world(
        REPOSITORY_ROOT,
        RELEASE_ROOT,
        run_root,
        truth_root,
        call,
        test_dense_variant_count=120,
    )
    private = load_json(
        truth_root / f"ibd/{call['world_id']}/{call['density']}/manifest.json"
    )

    assert manifest["variant_count"] == 60
    assert manifest["sample_count"] == 200
    assert manifest["pair_count"] == 100
    assert manifest["truth_available_to_worker"] is False
    assert len(private["segments"]) == 100
    assert all("private_truth" not in item["relative_path"] for item in manifest["artifacts"])
