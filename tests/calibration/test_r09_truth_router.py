"""Accès minimal à la vérité, exclusivement après clôture."""

from copy import deepcopy

import pytest

from effet_fondateur.calibration.contracts import TRUTH_CONSUMERS
from effet_fondateur.calibration.truth_router import (
    TruthRoutingError,
    allowed_truth_for_stage,
    authorize_truth_access,
)


def _truth_manifest() -> dict[str, object]:
    return {
        "release_id": "release-generic-v1",
        "disclosure_state": "SEALED_PRIVATE",
        "truth_artifacts": [
            {
                "truth_id": truth_id,
                "relative_path": f"private_truth/resource/{truth_id.lower()}.tsv",
                "sha256": f"{index:x}" * 64,
            }
            for index, truth_id in enumerate(TRUTH_CONSUMERS, start=1)
        ],
    }


def _closed_run() -> dict[str, object]:
    return {
        "run_id": "run-resource-001",
        "release_id": "release-generic-v1",
        "run_state": "CLOSED",
        "closed_at": "2026-11-12T12:00:00Z",
        "manifest_sha256": "f" * 64,
    }


def test_stage_receives_only_its_declared_truth() -> None:
    grant = authorize_truth_access(
        evaluation_id="EVAL-16B-001",
        stage_id="16B",
        requested_truth_ids=["T-FOUNDER"],
        truth_manifest=_truth_manifest(),
        run_manifest=_closed_run(),
    )
    assert grant.stage_id == "16B"
    assert [artifact[0] for artifact in grant.artifacts] == ["T-FOUNDER"]


def test_truth_cannot_be_opened_before_run_closure() -> None:
    run = _closed_run()
    run["run_state"] = "RUNNING"
    with pytest.raises(TruthRoutingError, match="run_not_closed"):
        authorize_truth_access(
            evaluation_id="EVAL-13-001",
            stage_id="13",
            requested_truth_ids=["T-FOUNDER"],
            truth_manifest=_truth_manifest(),
            run_manifest=run,
        )


def test_cross_stage_truth_request_is_refused() -> None:
    with pytest.raises(TruthRoutingError, match="truth_scope_exceeded"):
        authorize_truth_access(
            evaluation_id="EVAL-16C-001",
            stage_id="16C",
            requested_truth_ids=["T-FOUNDER"],
            truth_manifest=_truth_manifest(),
            run_manifest=_closed_run(),
        )


def test_truth_release_must_match_closed_run() -> None:
    truth = deepcopy(_truth_manifest())
    truth["release_id"] = "another-release"
    with pytest.raises(TruthRoutingError, match="truth_release_mismatch"):
        authorize_truth_access(
            evaluation_id="EVAL-03-001",
            stage_id="03",
            requested_truth_ids=["T-GENOTYPE"],
            truth_manifest=truth,
            run_manifest=_closed_run(),
        )


def test_stage_00_needs_no_private_truth_family() -> None:
    assert allowed_truth_for_stage("00") == frozenset()
