"""Étape 19 : rapport HTML révisable, prompt prudent et validation humaine."""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path
from time import monotonic
from typing import Any, Sequence

from effet_fondateur.audit import atomic_write_json, read_json
from effet_fondateur.contracts import build_file_artifact, validate_json_document
from effet_fondateur.orchestrator.state import utc_now
from effet_fondateur.reporting import build_report_draft
from effet_fondateur.stages.build_visualizations import _validate_producer_controls


METHOD_ID = "reviewable_scientific_report_v1"


def _build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(prog="build-report")
    parser.add_argument("--stage-inputs", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    return parser


def _parameters(parameters: dict[str, Any]) -> dict[str, Any]:
    if parameters.get("method", METHOD_ID) != METHOD_ID:
        raise ValueError("invalid_parameter:method")
    if parameters.get("ai_provider", "disabled") != "disabled":
        raise ValueError("external_ai_provider_not_implemented_or_approved")
    return {"method": METHOD_ID, "ai_provider": "disabled"}


def execute(stage_inputs_path: Path, output_dir: Path) -> int:
    """Assemble un brouillon local depuis les seuls artefacts signés du run."""
    started_at, started_clock = utc_now(), monotonic()
    stage_inputs = read_json(stage_inputs_path); validate_json_document(stage_inputs, "stage_inputs.schema.json")
    parameters = _parameters(stage_inputs["parameters"])
    run_dir = output_dir.parent.parent
    _validate_producer_controls(run_dir, stage_inputs["artifacts"])
    publication = build_report_draft(run_dir=run_dir, output_dir=output_dir, stage_inputs=stage_inputs)
    specs = [
        ("interpretation_facts", publication.facts_path, "interpretation_facts.schema.json", "application/json", "sensitive_genetic"),
        ("interpretation_prompt", publication.prompt_path, None, "text/plain", "sensitive_genetic"),
        ("interpretation_draft", publication.draft_path, "interpretation_draft.schema.json", "application/json", "sensitive_genetic"),
        ("report_draft_html", publication.html_path, None, "text/html", "sensitive_genetic"),
        ("report_editor_script", publication.editor_script_path, None, "text/javascript", "internal"),
        ("kinship_network", publication.kinship_network_path, None, "image/svg+xml", "sensitive_genetic"),
        ("report_validation", publication.validation_path, "report_validation.schema.json", "application/json", "internal"),
    ]
    artifacts = [build_file_artifact(
        physical_path=path, published_path=f"{stage_inputs['published_output_dir']}/{path.name}",
        artifact_id=artifact_id, artifact_type=artifact_id, media_type=media_type,
        producer_stage=stage_inputs["stage_name"], producer_signature=stage_inputs["signature"],
        schema_name=schema_name, schema_version="1.0.0" if schema_name else None,
        assembly=None, sample_set_id=None, variant_set_id=None, sensitivity=sensitivity,
    ) for artifact_id, path, schema_name, media_type, sensitivity in specs]
    outputs = {"schema_version": "1.0.0", "run_id": stage_inputs["run_id"], "stage_id": stage_inputs["stage_id"], "stage_name": stage_inputs["stage_name"], "signature": stage_inputs["signature"], "artifacts": artifacts}
    validate_json_document(outputs, "stage_outputs.schema.json"); atomic_write_json(output_dir / "stage_outputs.json", outputs)
    validation = read_json(publication.validation_path)
    audit = {
        "schema_version": "1.0.0", "run_id": stage_inputs["run_id"], "stage_id": stage_inputs["stage_id"], "stage_name": stage_inputs["stage_name"],
        "method_id": METHOD_ID, "signature": stage_inputs["signature"], "started_at": started_at, "completed_at": utc_now(), "duration_seconds": monotonic() - started_clock,
        "inputs": stage_inputs["artifacts"], "outputs": artifacts, "parameters": parameters,
        "tools": [{"tool": "deterministic_report_draft", "configured": "built-in", "version": "1.0.0"}],
        "counts": {"report_sections": 7},
        "metrics": {"report_status": validation["status"], "human_validation_required": True, "external_ai_called": False, "scientific_recalculation_performed": False, "composite_founder_score_calculated": False, "pseudonymized": True},
        "exclusions": [], "warnings": [{"code": "human_report_review_required"}],
        "checks": validation["checks"],
        "known_limits": ["Le brouillon n'est pas un rapport final.", "Une proposition rédactionnelle, même approuvée, ne démontre pas un effet fondateur.", "Le PDF n'est produit qu'après finalisation humaine contrôlée."],
        "expected_visualizations": ["kinship_network", "validated_stage18_figures"], "manual_validation_required": True,
    }
    validate_json_document(audit, "stage_audit.schema.json"); atomic_write_json(output_dir / "audit.json", audit)
    (output_dir / "checksums.sha256").write_text("".join(f"{item['sha256']}  {Path(item['path']).name}\n" for item in artifacts), encoding="utf-8")
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
