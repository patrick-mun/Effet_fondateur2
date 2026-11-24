"""Exemples minimaux valides pour chacun des neuf schémas R09."""

from pathlib import Path

import pytest

from effet_fondateur.calibration.contracts import (
    STAGE_IDS,
    TRUTH_CONSUMERS,
    VIEW_CONSUMERS,
    load_json,
    validate_document,
)


ROOT = Path(__file__).resolve().parents[2]
SCHEMAS = ROOT / "calibration/schemas"
SHA = "a" * 64


def _documents() -> dict[str, dict[str, object]]:
    truth_artifacts = [
        {
            "truth_id": truth_id,
            "relative_path": f"private_truth/calibration/{truth_id.lower()}.tsv",
            "sha256": SHA,
            "authorized_evaluators": [f"EVAL-{stage}" for stage in consumers],
        }
        for truth_id, consumers in TRUTH_CONSUMERS.items()
    ]
    views = [
        {
            "view_type": view_id,
            "opaque_view_id": f"VIEW-{index:016x}",
            "parent_sha256": SHA,
            "transformation_id": "derive_view_v1",
            "consumer_stages": list(consumers),
            "inventory_path": f"blinded_inputs/resource/view-{index}/inventory.json",
        }
        for index, (view_id, consumers) in enumerate(VIEW_CONSUMERS.items(), start=1)
    ]
    return {
        "panel_release_manifest.schema.json": {
            "schema_version": "1.0.0",
            "release_id": "panel-master-generic-v1-resource",
            "release_design": "panel_master_generic_v1",
            "assembly": "GRCh38",
            "stage_ids": list(STAGE_IDS),
            "partitions": ["CALIBRATION", "VALIDATION_BLIND", "RESOURCE_ONLY"],
            "scenario_design": {"scale": "micro"},
            "target": {"selection": "PROFILE_CONFIGURED_AUTOSOMAL", "hard_coded": False},
            "qualification_flags": {"calibrated": False, "founder_effect_proven": False, "production_promotion": False},
        },
        "subject_registry.schema.json": {
            "schema_version": "1.0.0",
            "release_id": "release-resource-v1",
            "relative_path": "governance/subject_registry.tsv",
            "sha256": SHA,
            "row_count": 60,
            "columns": ["master_id", "partition_id", "isolation_unit_id", "family_id", "role", "batch", "activation_rank", "pseudonym"],
            "partition_values": ["CALIBRATION", "VALIDATION_BLIND", "RESOURCE_ONLY"],
        },
        "master_truth_manifest.schema.json": {
            "schema_version": "1.0.0",
            "release_id": "release-resource-v1",
            "disclosure_state": "SEALED_PRIVATE",
            "truth_artifacts": truth_artifacts,
        },
        "panel_view_manifest.schema.json": {
            "schema_version": "1.0.0",
            "release_id": "release-resource-v1",
            "partition_id": "RESOURCE_ONLY",
            "views": views,
        },
        "blinded_input_inventory.schema.json": {
            "schema_version": "1.0.0",
            "inventory_id": "inventory-resource-001",
            "release_id": "release-resource-v1",
            "partition_id": "RESOURCE_ONLY",
            "opaque_view_id": "VIEW-0000000000000001",
            "artifacts": [{"artifact_id": "input", "relative_path": "inputs/panel.vcf.gz", "sha256": SHA, "size_bytes": 10, "media_type": "application/gzip"}],
        },
        "pipeline_run_manifest.schema.json": {
            "schema_version": "1.0.0",
            "run_id": "run-resource-001",
            "release_id": "release-resource-v1",
            "partition_id": "RESOURCE_ONLY",
            "run_state": "CLOSED",
            "closed_at": "2026-11-12T12:00:00Z",
            "manifest_sha256": SHA,
            "consumed_inventory_sha256": SHA,
            "stage_outputs": [{"stage_id": "00", "artifact_id": "summary", "relative_path": "runs/run-resource-001/00/summary.json", "sha256": SHA}],
            "dependencies": {},
        },
        "stage_evaluation_manifest.schema.json": {
            "schema_version": "1.0.0",
            "evaluation_id": "evaluation-00-resource-001",
            "stage_id": "00",
            "release_id": "release-resource-v1",
            "closed_run_id": "run-resource-001",
            "truth_grants": [],
            "metrics_path": "evaluations/00/ETAPE_00.metrics.tsv",
            "metrics_sha256": SHA,
            "status": "PASS",
            "decision_scope": "TECHNICAL_INTEGRITY",
            "scientific_performance_evaluated": False,
            "technical_checks": ["OUTPUTS_DECLARED_AND_HASHED"],
            "counts": {"N_EXPECTED": 8, "N_ELIGIBLE": 8, "N_EVALUABLE": 8, "N_EXCLUDED": 0},
            "criteria_version": "1.0.0",
            "closed": True,
        },
        "readiness_receipt.schema.json": load_json(ROOT / "calibration/readiness/R09_READY_TO_RUN.json"),
        "calibration_decision.schema.json": {
            "schema_version": "1.0.0",
            "decision_id": "decision-validation-001",
            "release_id": "release-resource-v1",
            "stage_decisions": [{"stage_id": stage, "status": "PASS"} for stage in STAGE_IDS],
            "qualification": "NOT_QUALIFIED",
            "restrictions": ["validation_not_run"],
            "compatibility": {"release_design": "panel_master_generic_v1"},
            "human_signature": {"reviewer_id": "reviewer-opaque", "signed_sha256": SHA},
            "decided_at": "2026-11-12T12:00:00Z",
        },
    }


@pytest.mark.parametrize("schema_name", tuple(_documents()))
def test_minimal_document_validates_against_schema(schema_name: str) -> None:
    validate_document(_documents()[schema_name], SCHEMAS / schema_name)
