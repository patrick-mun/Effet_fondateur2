from __future__ import annotations

import json
from pathlib import Path

import pytest

from effet_fondateur.ld.matched import (
    PhasedVariant,
    centered_pair_metrics,
    count_neighbors_by_bin,
    density_is_matched,
    local_background_metrics,
    parse_phased_genotype,
    r_squared,
)
from scripts.analyze_phase15_ld_v2 import (
    Phase15LdV2Error,
    _load_verified_stage_artifacts,
    _sha256,
    _validated_artifact_path,
)


def _variant(
    variant_id: str, position_cm: float, haplotypes: tuple[int | None, ...],
    *, position_bp: int | None = None,
) -> PhasedVariant:
    return PhasedVariant(
        chromosome="19",
        variant_id=variant_id,
        position_bp=position_bp or round(position_cm * 1_000_000),
        position_cm=position_cm,
        haplotypes=haplotypes,
    )


def test_phased_r2_is_pairwise_and_refuses_monomorphic_marker() -> None:
    first = (0, 0, 1, 1, None, 1)
    second = (0, 0, 1, 1, 0, None)

    value, called = r_squared(first, second, minimum_called=4)

    assert called == 4
    assert value == pytest.approx(1.0)
    assert r_squared(first, (0, 0, 0, 0, 0, 0), minimum_called=4)[0] is None


def test_genotype_parser_requires_phased_biallelic_calls() -> None:
    assert parse_phased_genotype("0|1:0.99") == (0, 1)
    assert parse_phased_genotype(".|.") == (None, None)
    with pytest.raises(ValueError, match="unphased"):
        parse_phased_genotype("0/1")
    with pytest.raises(ValueError, match="non_biallelic"):
        parse_phased_genotype("0|2")


def test_centered_ld_keeps_rare_anchor_and_filters_neighbor_maf() -> None:
    anchor = _variant("target", 10.0, (0, 0, 0, 0, 0, 1))
    common = _variant("common", 10.08, (0, 0, 0, 1, 0, 1))
    rare = _variant("rare", 10.12, (0, 0, 0, 0, 0, 1))
    distant = _variant("distant", 11.2, (0, 1, 0, 1, 0, 1))

    rows = centered_pair_metrics(
        anchor, [anchor, common, rare, distant], distance_bins_cm=(0.1, 0.25, 1.0),
        minimum_neighbor_maf=0.2, minimum_called_haplotypes=6,
    )

    assert [row.neighbor_id for row in rows] == ["common"]
    assert rows[0].distance_bin == "0-0.1_cM"


def test_local_background_excludes_target_and_uses_common_common_pairs() -> None:
    variants = [
        _variant("target", 1.0, (0, 0, 0, 0, 0, 1)),
        _variant("a", 1.05, (0, 0, 1, 1, 0, 1)),
        _variant("b", 1.09, (0, 0, 1, 1, 0, 1)),
        _variant("c", 1.4, (0, 1, 0, 1, 0, 1)),
    ]

    result = local_background_metrics(
        variants, excluded_variant_id="target", distance_bins_cm=(0.1, 0.5),
        minimum_maf=0.2, minimum_called_haplotypes=6,
    )

    assert result["0-0.1_cM"] == pytest.approx([1.0])
    assert len(result["0.1-0.5_cM"]) == 2


def test_density_matching_is_prespecified_per_distance_bin() -> None:
    target = {"0-0.1_cM": 10, "0.1-0.25_cM": 20}

    assert density_is_matched(
        {"0-0.1_cM": 5, "0.1-0.25_cM": 40}, target, ratio=2.0,
    )
    assert not density_is_matched(
        {"0-0.1_cM": 4, "0.1-0.25_cM": 40}, target, ratio=2.0,
    )
    assert density_is_matched(
        {"0-0.1_cM": 8, "0.1-0.25_cM": 999},
        {"0-0.1_cM": 10, "0.1-0.25_cM": 0}, ratio=2.0,
    )


def test_density_count_can_require_the_same_evaluable_r2_as_analysis() -> None:
    anchor = _variant("target", 1.0, (0, 0, 0, 1, 1, 1))
    complete = _variant("complete", 1.05, (0, 0, 1, 1, 0, 1))
    too_missing = _variant("too_missing", 1.08, (0, 1, 1, None, None, None))
    variants = [anchor, complete, too_missing]

    eligible = count_neighbors_by_bin(
        anchor, variants, distance_bins_cm=(0.1,), minimum_neighbor_maf=0.2,
    )
    evaluable = count_neighbors_by_bin(
        anchor, variants, distance_bins_cm=(0.1,), minimum_neighbor_maf=0.2,
        minimum_called_haplotypes=4,
    )
    metrics = centered_pair_metrics(
        anchor, variants, distance_bins_cm=(0.1,), minimum_neighbor_maf=0.2,
        minimum_called_haplotypes=4,
    )

    assert eligible == {"0-0.1_cM": 2}
    assert evaluable == {"0-0.1_cM": 1}
    assert [metric.neighbor_id for metric in metrics] == ["complete"]


def _stage_outputs_document(run_id: str, signature: str) -> dict[str, object]:
    return {
        "schema_version": "1.0.0",
        "run_id": run_id,
        "stage_id": "12A",
        "stage_name": "phase_autosomal_panel",
        "signature": signature,
        "artifacts": [{
            "artifact_id": "autosomal_phased_chr1_bcf",
            "artifact_type": "autosomal_phased_bcf",
            "path": "stages/12A_phase_autosomal_panel/chr1/study.phased.bcf",
            "media_type": "application/octet-stream",
            "schema_name": None,
            "schema_version": None,
            "sha256": "b" * 64,
            "producer_stage": "12A_phase_autosomal_panel",
            "producer_signature": signature,
            "assembly": "GRCh38",
            "sample_set_id": "synthetic_samples",
            "variant_set_id": "synthetic_variants",
            "sensitivity": "sensitive_genetic",
        }],
    }


def test_stage_outputs_must_match_checksum_anchored_in_run_manifest(
    tmp_path: Path,
) -> None:
    run_id = "synthetic_run"
    signature = "a" * 64
    outputs_path = tmp_path / "stages/12A_phase_autosomal_panel/stage_outputs.json"
    outputs_path.parent.mkdir(parents=True)
    outputs_path.write_text(
        json.dumps(_stage_outputs_document(run_id, signature)), encoding="utf-8",
    )
    manifest = {
        "run_id": run_id,
        "stages": [{
            "stage_id": "12A",
            "stage_name": "phase_autosomal_panel",
            "state": "SUCCEEDED",
            "signature": signature,
            "audit_path": "stages/12A_phase_autosomal_panel/audit.json",
            "stage_outputs_sha256": _sha256(outputs_path),
        }],
    }

    artifacts = _load_verified_stage_artifacts(
        tmp_path, manifest, "phase_autosomal_panel", outputs_path,
    )
    assert set(artifacts) == {"autosomal_phased_chr1_bcf"}

    outputs_path.write_text("{}\n", encoding="utf-8")
    with pytest.raises(Phase15LdV2Error, match="outputs_checksum_mismatch"):
        _load_verified_stage_artifacts(
            tmp_path, manifest, "phase_autosomal_panel", outputs_path,
        )


def test_declared_artifact_checksum_is_verified_before_use(tmp_path: Path) -> None:
    relative_path = "stages/12A_phase_autosomal_panel/chr1/study.phased.bcf"
    artifact_path = tmp_path / relative_path
    artifact_path.parent.mkdir(parents=True)
    artifact_path.write_bytes(b"verified-bcf")
    artifacts = {
        "autosomal_phased_chr1_bcf": {
            "artifact_id": "autosomal_phased_chr1_bcf",
            "artifact_type": "autosomal_phased_bcf",
            "path": relative_path,
            "sha256": _sha256(artifact_path),
        }
    }

    assert _validated_artifact_path(
        tmp_path, artifacts, "autosomal_phased_chr1_bcf", relative_path,
        "autosomal_phased_bcf",
    ) == artifact_path

    artifact_path.write_bytes(b"modified-bcf")
    with pytest.raises(Phase15LdV2Error, match="artifact_checksum_mismatch"):
        _validated_artifact_path(
            tmp_path, artifacts, "autosomal_phased_chr1_bcf", relative_path,
            "autosomal_phased_bcf",
        )
