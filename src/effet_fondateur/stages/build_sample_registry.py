"""Étape 02 : validation et publication du registre maître des échantillons."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import re
import shutil
from datetime import datetime
from pathlib import Path, PurePosixPath
from time import monotonic
from typing import Any, Sequence

from effet_fondateur.audit import atomic_write_json, read_json, sha256_file
from effet_fondateur.contracts import (
    DocumentValidationError,
    TableValidationError,
    build_file_artifact,
    load_pipeline_config,
    validate_json_document,
    validate_samples_master,
    validate_tsv_table,
)
from effet_fondateur.orchestrator.state import utc_now


REVIEW_COLUMNS = (
    "SAMPLE_ID",
    "PSEUDONYM",
    "SOURCE_INVENTORY_STATUS",
    "TARGET_GENOTYPE_STATUS",
    "PROVENANCE_STATUS",
    "REVIEW_STATUS",
)
STABLE_IDENTIFIER_PATTERN = re.compile(r"^[A-Za-z0-9][A-Za-z0-9._-]*$")


class SampleRegistryError(ValueError):
    """Signale une incohérence entre métadonnées, inventaire et approbation."""


def _build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(prog="build-sample-registry")
    parser.add_argument("--stage-inputs", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    return parser


def _physical_artifact_path(artifact: dict[str, Any], run_dir: Path) -> Path:
    artifact_path = Path(artifact["path"])
    if artifact_path.is_absolute():
        return artifact_path
    if PurePosixPath(artifact["path"]).parts[:1] == ("stages",):
        return run_dir / artifact_path
    return Path.cwd() / artifact_path


def _input_artifact(
    stage_inputs: dict[str, Any],
    artifact_id: str,
    run_dir: Path,
) -> tuple[dict[str, Any], Path]:
    matches = [
        artifact
        for artifact in stage_inputs["artifacts"]
        if artifact["artifact_id"] == artifact_id
    ]
    if len(matches) != 1:
        raise SampleRegistryError(f"input_artifact_{artifact_id}_missing_or_ambiguous")
    artifact = matches[0]
    physical_path = _physical_artifact_path(artifact, run_dir)
    if (
        not physical_path.is_file()
        or sha256_file(physical_path) != artifact["sha256"]
    ):
        raise SampleRegistryError(f"input_artifact_{artifact_id}_modified")
    return artifact, physical_path


def _validate_manual_approval(parameters: dict[str, Any]) -> dict[str, Any] | None:
    approval = parameters.get("manual_approval")
    if approval is None:
        return None
    if not isinstance(approval, dict) or set(approval) != {
        "approved",
        "decision_id",
        "reviewer_role",
        "approved_at",
    }:
        raise SampleRegistryError("invalid_manual_approval_contract")
    if approval["approved"] is not True:
        raise SampleRegistryError("manual_approval_must_be_true_or_omitted")
    for field in ("decision_id", "reviewer_role"):
        value = approval[field]
        if (
            not isinstance(value, str)
            or STABLE_IDENTIFIER_PATTERN.fullmatch(value) is None
        ):
            raise SampleRegistryError(f"invalid_manual_approval_{field}")
    try:
        approved_at = datetime.fromisoformat(
            str(approval["approved_at"]).replace("Z", "+00:00")
        )
    except ValueError as error:
        raise SampleRegistryError("invalid_manual_approval_approved_at") from error
    if approved_at.tzinfo is None:
        raise SampleRegistryError("manual_approval_approved_at_requires_timezone")
    return approval


def _pseudonym(project_id: str, sample_id: str) -> str:
    digest = hashlib.sha256(f"{project_id}:{sample_id}".encode("utf-8")).hexdigest()
    return f"P_{digest[:12].upper()}"


def _write_review_table(
    review_path: Path,
    rows: tuple[dict[str, Any], ...],
    project_id: str,
    approved: bool,
) -> None:
    review_rows = []
    pseudonyms: set[str] = set()
    for row in rows:
        pseudonym = _pseudonym(project_id, row["SAMPLE_ID"])
        if pseudonym in pseudonyms:
            raise SampleRegistryError("pseudonym_collision")
        pseudonyms.add(pseudonym)
        genotype_recorded = row["TARGET_GENOTYPE"] is not None
        review_rows.append(
            {
                "SAMPLE_ID": row["SAMPLE_ID"],
                "PSEUDONYM": pseudonym,
                "SOURCE_INVENTORY_STATUS": "MATCHED",
                "TARGET_GENOTYPE_STATUS": (
                    "RECORDED" if genotype_recorded else "MISSING"
                ),
                "PROVENANCE_STATUS": (
                    "RECORDED" if genotype_recorded else "NOT_APPLICABLE"
                ),
                "REVIEW_STATUS": "APPROVED" if approved else "PENDING",
            }
        )
    with review_path.open("w", encoding="utf-8", newline="") as output_file:
        writer = csv.DictWriter(
            output_file,
            fieldnames=REVIEW_COLUMNS,
            delimiter="\t",
            lineterminator="\n",
        )
        writer.writeheader()
        writer.writerows(review_rows)


def execute(stage_inputs_path: Path, output_dir: Path) -> int:
    """Valide, recoupe et publie le registre sans créer de génotype implicite."""
    started_at = utc_now()
    started_clock = monotonic()
    stage_inputs = read_json(stage_inputs_path)
    validate_json_document(stage_inputs, "stage_inputs.schema.json")
    run_dir = output_dir.parent.parent
    config = load_pipeline_config(run_dir / "config.resolved.yaml")
    configured_source_root = config["inputs"]["acpa_samples_dir"]
    if configured_source_root is None:
        raise SampleRegistryError("missing_acpa_samples_dir")
    source_root = Path(configured_source_root)
    if not source_root.is_absolute():
        source_root = Path.cwd() / source_root

    metadata_artifact, metadata_path = _input_artifact(
        stage_inputs,
        "config_input_sample_metadata",
        run_dir,
    )
    inventory_artifact, inventory_path = _input_artifact(
        stage_inputs,
        "source_inventory",
        run_dir,
    )
    qc_artifact, qc_path = _input_artifact(stage_inputs, "source_qc", run_dir)
    inventory = validate_tsv_table(inventory_path, "source_inventory.schema.json")
    source_qc = validate_tsv_table(qc_path, "source_qc.schema.json")
    if any(row["STATUS"] == "FAIL" for row in source_qc.rows):
        raise SampleRegistryError("source_qc_contains_failure")

    validated_metadata = validate_samples_master(metadata_path, source_root)
    inventory_source_files = {row["SOURCE_FILE"] for row in inventory.rows}
    metadata_source_files = {
        row["SOURCE_FILE"] for row in validated_metadata.rows
    }
    missing_from_inventory = metadata_source_files - inventory_source_files
    if missing_from_inventory:
        raise SampleRegistryError("metadata_source_missing_from_inventory")

    approval = _validate_manual_approval(stage_inputs["parameters"])
    approved = approval is not None
    output_dir.mkdir(parents=True, exist_ok=True)
    master_path = output_dir / "samples.master.tsv"
    shutil.copyfile(metadata_path, master_path)
    published_master = validate_samples_master(master_path, source_root)
    review_path = output_dir / "sample_registry_review.tsv"
    _write_review_table(
        review_path,
        published_master.rows,
        config["project"]["project_id"],
        approved,
    )
    validate_tsv_table(review_path, "sample_registry_review.schema.json")

    sample_set_id = f"samples_{published_master.sha256[:12]}"
    producer_stage = f"{stage_inputs['stage_id']}_{stage_inputs['stage_name']}"
    published_dir = PurePosixPath(stage_inputs["published_output_dir"])
    master_artifact = build_file_artifact(
        physical_path=master_path,
        published_path=str(published_dir / master_path.name),
        artifact_id="samples_master",
        artifact_type="sample_table",
        media_type="text/tab-separated-values",
        producer_stage=producer_stage,
        producer_signature=stage_inputs["signature"],
        schema_name="samples_master.schema.json",
        schema_version="1.0.0",
        assembly=config["project"]["assembly"],
        sample_set_id=sample_set_id,
        sensitivity="sensitive_genetic",
    )
    review_artifact = build_file_artifact(
        physical_path=review_path,
        published_path=str(published_dir / review_path.name),
        artifact_id="sample_registry_review",
        artifact_type="sample_registry_review",
        media_type="text/tab-separated-values",
        producer_stage=producer_stage,
        producer_signature=stage_inputs["signature"],
        schema_name="sample_registry_review.schema.json",
        schema_version="1.0.0",
        assembly=config["project"]["assembly"],
        sample_set_id=sample_set_id,
        sensitivity="sensitive_genetic",
    )
    artifacts = [master_artifact, review_artifact]
    stage_outputs = {
        "schema_version": "1.0.0",
        "run_id": stage_inputs["run_id"],
        "stage_id": stage_inputs["stage_id"],
        "stage_name": stage_inputs["stage_name"],
        "signature": stage_inputs["signature"],
        "artifacts": artifacts,
    }
    validate_json_document(stage_outputs, "stage_outputs.schema.json")
    atomic_write_json(output_dir / "stage_outputs.json", stage_outputs)

    recorded_genotypes = sum(
        row["TARGET_GENOTYPE"] is not None for row in published_master.rows
    )
    unassigned_sources = inventory_source_files - metadata_source_files
    audit = {
        "schema_version": "1.0.0",
        "run_id": stage_inputs["run_id"],
        "stage_id": stage_inputs["stage_id"],
        "stage_name": stage_inputs["stage_name"],
        "method_id": "build_sample_registry_v1",
        "signature": stage_inputs["signature"],
        "started_at": started_at,
        "completed_at": utc_now(),
        "duration_seconds": monotonic() - started_clock,
        "inputs": [metadata_artifact, inventory_artifact, qc_artifact],
        "outputs": artifacts,
        "parameters": stage_inputs["parameters"],
        "tools": [],
        "counts": {
            "samples": published_master.row_count,
            "target_genotypes_recorded": recorded_genotypes,
            "target_genotypes_missing": published_master.row_count - recorded_genotypes,
            "unassigned_inventory_sources": len(unassigned_sources),
        },
        "metrics": {},
        "exclusions": [],
        "warnings": (
            [{"code": "unassigned_inventory_sources", "count": len(unassigned_sources)}]
            if unassigned_sources
            else []
        ),
        "checks": [
            {"check": "samples_master_contract", "status": "PASS"},
            {"check": "source_inventory_join", "status": "PASS"},
            {
                "check": "manual_approval",
                "status": "PASS" if approved else "PENDING",
            },
        ],
        "known_limits": [
            "Les génotypes cibles manquants restent manquants et ne sont jamais inférés."
        ],
        "expected_visualizations": [],
        "manual_validation_required": not approved,
    }
    validate_json_document(audit, "stage_audit.schema.json")
    atomic_write_json(output_dir / "audit.json", audit)
    (output_dir / "checksums.sha256").write_text(
        "".join(
            f"{artifact['sha256']}  {Path(str(artifact['path'])).name}\n"
            for artifact in artifacts
        ),
        encoding="utf-8",
    )
    return 0


def main(arguments: Sequence[str] | None = None) -> int:
    """Point d'entrée autonome respectant les codes de retour V2."""
    parser = _build_parser()
    parsed_arguments = parser.parse_args(arguments)
    try:
        return execute(parsed_arguments.stage_inputs, parsed_arguments.output_dir)
    except (
        DocumentValidationError,
        TableValidationError,
        SampleRegistryError,
        OSError,
        ValueError,
        json.JSONDecodeError,
    ) as error:
        parser.error(str(error))
    return 2


if __name__ == "__main__":
    raise SystemExit(main())
