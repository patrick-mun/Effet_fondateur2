"""Étape 18 : figures consolidées et index audité du run courant."""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path
from time import monotonic
from typing import Any, Sequence

from effet_fondateur.audit import atomic_write_json, read_json, sha256_file
from effet_fondateur.contracts import build_file_artifact, validate_json_document
from effet_fondateur.orchestrator.state import utc_now
from effet_fondateur.visualization import build_consolidated_figures


def _build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(prog="build-visualizations")
    parser.add_argument("--stage-inputs", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    return parser


def _parameters(parameters: dict[str, Any]) -> dict[str, Any]:
    if parameters.get("method", "validated_current_run_figures_v1") != "validated_current_run_figures_v1":
        raise ValueError("invalid_parameter:method")
    return {"method": "validated_current_run_figures_v1"}


def _validate_producer_controls(run_dir: Path, input_artifacts: list[dict[str, Any]]) -> None:
    """Lie chaque entrée au manifest, à l'audit et aux sorties de son producteur."""
    manifest = read_json(run_dir / "manifest.json")
    records = {record["stage_name"]: record for record in manifest.get("stages", [])}
    for producer in sorted({artifact["producer_stage"] for artifact in input_artifacts}):
        producer_artifacts = [artifact for artifact in input_artifacts if artifact["producer_stage"] == producer]
        record = records.get(producer)
        if record is None or record.get("state") not in {"SUCCEEDED", "CACHED"}:
            raise ValueError("producer_stage_not_validated_in_current_run")
        expected_signatures = {artifact["producer_signature"] for artifact in producer_artifacts}
        if len(expected_signatures) != 1 or record.get("signature") not in expected_signatures:
            raise ValueError("producer_signature_mismatch")
        audit_path = run_dir / record["audit_path"]
        outputs_path = audit_path.parent / "stage_outputs.json"
        if sha256_file(audit_path) != record.get("audit_sha256"):
            raise ValueError("producer_audit_checksum_mismatch")
        if sha256_file(outputs_path) != record.get("stage_outputs_sha256"):
            raise ValueError("producer_outputs_checksum_mismatch")
        audit = read_json(audit_path); validate_json_document(audit, "stage_audit.schema.json")
        outputs = read_json(outputs_path); validate_json_document(outputs, "stage_outputs.schema.json")
        if audit["signature"] != record["signature"] or outputs["signature"] != record["signature"]:
            raise ValueError("producer_control_signature_mismatch")
        published = {artifact["artifact_id"]: artifact for artifact in outputs["artifacts"]}
        for artifact in producer_artifacts:
            if published.get(artifact["artifact_id"]) != artifact:
                raise ValueError("producer_artifact_descriptor_mismatch")


def execute(stage_inputs_path: Path, output_dir: Path) -> int:
    """Valide les artefacts 13–17, rend les domaines et publie leur index."""
    started_at, started_clock = utc_now(), monotonic()
    stage_inputs = read_json(stage_inputs_path)
    validate_json_document(stage_inputs, "stage_inputs.schema.json")
    parameters = _parameters(stage_inputs["parameters"])
    run_dir = output_dir.parent.parent
    _validate_producer_controls(run_dir, stage_inputs["artifacts"])
    figures_dir = output_dir / "figures"
    results = build_consolidated_figures(run_dir=run_dir, output_dir=figures_dir, stage_inputs=stage_inputs)
    index_rows = []
    for result in results:
        index_rows.append({
            "figure_id": result.figure_id, "domain": result.domain, "status": result.status,
            "figure_path": None if result.figure_path is None else f"figures/{result.figure_path.name}",
            "provenance_path": f"figures/{result.provenance_path.name}",
            "sensitivity": "sensitive_genetic", "blocked_reason": result.blocked_reason,
        })
    index = {"schema_version": "1.0.0", "run_id": stage_inputs["run_id"], "method_id": parameters["method"], "composite_founder_score_calculated": False, "figures": index_rows}
    validate_json_document(index, "figure_index.schema.json")
    index_path = output_dir / "figure_index.json"; atomic_write_json(index_path, index)
    rendered_count = sum(item.status == "RENDERED" for item in results)
    not_evaluated_count = sum(item.status == "NOT_EVALUATED" for item in results)
    blocked_count = sum(item.status == "BLOCKED" for item in results)
    completeness = {"schema_version": "1.0.0", "run_id": stage_inputs["run_id"], "expected_domain_count": 6, "rendered_count": rendered_count, "not_evaluated_count": not_evaluated_count, "blocked_count": blocked_count, "complete_for_scientific_report": blocked_count == 0}
    validate_json_document(completeness, "visualization_completeness.schema.json")
    completeness_path = output_dir / "visualization_completeness.json"; atomic_write_json(completeness_path, completeness)
    specs: list[tuple[str, Path, str, str]] = [
        ("figure_index", index_path, "figure_index.schema.json", "internal"),
        ("visualization_completeness", completeness_path, "visualization_completeness.schema.json", "internal"),
    ]
    for result in results:
        specs.append((f"figure_provenance_{result.figure_id}", result.provenance_path, "figure_provenance.schema.json", "sensitive_genetic"))
        if result.figure_path is not None:
            specs.append((f"figure_{result.figure_id}", result.figure_path, "", "sensitive_genetic"))
    artifacts = [build_file_artifact(
        physical_path=path,
        published_path=f"{stage_inputs['published_output_dir']}/{path.relative_to(output_dir).as_posix()}",
        artifact_id=artifact_id, artifact_type=artifact_id,
        media_type="image/svg+xml" if path.suffix == ".svg" else "application/json",
        producer_stage=stage_inputs["stage_name"], producer_signature=stage_inputs["signature"],
        schema_name=schema_name or None, schema_version="1.0.0" if schema_name else None,
        assembly=None, sample_set_id=None, variant_set_id=None, sensitivity=sensitivity,
    ) for artifact_id, path, schema_name, sensitivity in specs]
    stage_outputs = {"schema_version": "1.0.0", "run_id": stage_inputs["run_id"], "stage_id": stage_inputs["stage_id"], "stage_name": stage_inputs["stage_name"], "signature": stage_inputs["signature"], "artifacts": artifacts}
    validate_json_document(stage_outputs, "stage_outputs.schema.json"); atomic_write_json(output_dir / "stage_outputs.json", stage_outputs)
    audit = {
        "schema_version": "1.0.0", "run_id": stage_inputs["run_id"], "stage_id": stage_inputs["stage_id"], "stage_name": stage_inputs["stage_name"],
        "method_id": parameters["method"], "signature": stage_inputs["signature"], "started_at": started_at, "completed_at": utc_now(), "duration_seconds": monotonic() - started_clock,
        "inputs": stage_inputs["artifacts"], "outputs": artifacts, "parameters": parameters,
        "tools": [{"tool": "python_svg_renderer", "configured": Path(sys.executable).name, "version": sys.version.split()[0]}],
        "counts": {"expected_domains": 6, "rendered": rendered_count, "not_evaluated": not_evaluated_count, "blocked": blocked_count},
        "metrics": {"complete_for_scientific_report": blocked_count == 0, "pseudonymized": True, "sensitivity": "sensitive_genetic", "composite_founder_score_calculated": False},
        "exclusions": [], "warnings": [{"code": "figure_blocked", "count": blocked_count}] if blocked_count else [],
        "checks": [{"check": "current_run_only", "status": "PASS"}, {"check": "source_checksums_and_signatures", "status": "PASS" if blocked_count == 0 else "WARN"}, {"check": "domain_separation", "status": "PASS"}, {"check": "primary_exploratory_separation", "status": "PASS"}, {"check": "no_scientific_recalculation", "status": "PASS"}, {"check": "no_composite_founder_score", "status": "PASS"}],
        "known_limits": ["Les figures reprennent les limites et non-évaluations des producteurs.", "Une figure ne constitue pas une preuve causale, IBD ou d'origine fondatrice unique.", "Les SVG et provenances restent classés sensitive_genetic malgré la pseudonymisation."],
        "expected_visualizations": [item.figure_id for item in results], "manual_validation_required": True,
    }
    validate_json_document(audit, "stage_audit.schema.json"); atomic_write_json(output_dir / "audit.json", audit)
    (output_dir / "checksums.sha256").write_text("".join(f"{item['sha256']}  {item['path'].removeprefix(stage_inputs['published_output_dir'] + '/')}\n" for item in artifacts), encoding="utf-8")
    return 0


def main(arguments: Sequence[str] | None = None) -> int:
    parsed = _build_parser().parse_args(arguments)
    try:
        return execute(parsed.stage_inputs, parsed.output_dir)
    except Exception as error:
        sys.stderr.write(f"{error}\n")
        return 2 if isinstance(error, (ValueError, OSError, json.JSONDecodeError)) else 5


if __name__ == "__main__":
    raise SystemExit(main())
