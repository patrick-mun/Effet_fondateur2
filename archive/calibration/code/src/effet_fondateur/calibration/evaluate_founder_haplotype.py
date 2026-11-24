"""D24 : évaluation close de l’haplotype fondateur contre une vérité privée."""

from __future__ import annotations

import argparse
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
import shutil
import tempfile
from typing import Any

from effet_fondateur.audit import sha256_file
from effet_fondateur.calibration.founder_haplotype_oracle import (
    METHOD_ID,
    acceptance_results,
    aggregate_results,
    evaluate_case,
    validate_case_balance,
)
from effet_fondateur.calibration.qc_bench import write_json
from effet_fondateur.calibration.technical_scenarios import new_workspace_run
from effet_fondateur.contracts import validate_json_document, validate_tsv_table


_REPOSITORY = Path(__file__).resolve().parents[3]
_CALIBRATION_SCHEMAS = _REPOSITORY / "calibration/schemas"
_PROTOCOL = _REPOSITORY / "calibration/protocoles/D24_FOUNDER_HAPLOTYPE_PRESPECIFICATION.md"


def load_configuration(path: Path) -> dict[str, Any]:
    """Charge la configuration D24 fermée et refuse les variantes implicites."""
    if not path.is_file() or path.is_symlink():
        raise ValueError("d24_configuration_invalid")
    document = json.loads(path.read_text(encoding="utf-8"))
    validate_json_document(
        document,
        "founder_haplotype_bench_config.schema.json",
        schemas_dir=_CALIBRATION_SCHEMAS,
    )
    if document["method_id"] != METHOD_ID:
        raise ValueError("d24_method_changed")
    scenario_ids = [row["scenario_id"] for row in document["scenarios"]]
    scenario_seeds = [row["random_seed"] for row in document["scenarios"]]
    if len(scenario_ids) != len(set(scenario_ids)):
        raise ValueError("d24_scenario_definition_duplicate")
    if len(scenario_seeds) != len(set(scenario_seeds)):
        raise ValueError("d24_scenario_seed_duplicate")
    roles = {role: 0 for role in ("POSITIVE", "NEGATIVE")}
    for scenario in document["scenarios"]:
        roles[scenario["truth_role"]] += 1
    if roles["POSITIVE"] < document["minimum_positive_scenarios"]:
        raise ValueError("d24_positive_scenario_count_insufficient")
    if roles["NEGATIVE"] < document["minimum_negative_scenarios"]:
        raise ValueError("d24_negative_scenario_count_insufficient")
    return document


def _safe_case_root(path: Path) -> Path:
    root = path.absolute()
    folded = [part.casefold() for part in root.parts]
    if (
        not root.is_dir()
        or root.is_symlink()
        or "calibration_workspace" not in folded
        or "dock6_real" in folded
        or "cache" in folded
        or any(left == "data" and right == "input" for left, right in zip(folded, folded[1:]))
        or any(parent.is_symlink() for parent in (root, *root.parents))
    ):
        raise ValueError("d24_case_root_outside_workspace")
    return root


def _case_digest(paths: list[Path]) -> str:
    payload = [sha256_file(path) for path in paths]
    return hashlib.sha256(json.dumps(payload, separators=(",", ":")).encode()).hexdigest()


def load_case(case_root: Path, config: dict[str, Any]) -> tuple[dict[str, Any], str]:
    """Ouvre la vérité seulement après présence des sorties de production 13."""
    root = _safe_case_root(case_root)
    truth_path = root / "private_truth/founder_haplotype_case.json"
    output_root = root / "stage13"
    summary_path = output_root / "founder_analysis_summary.json"
    consensus_path = output_root / "founder_consensus.tsv"
    segments_path = output_root / "founder_segments.tsv"
    paths = [truth_path, summary_path, consensus_path, segments_path]
    if any(not path.is_file() or path.is_symlink() for path in paths):
        raise ValueError("d24_case_artifact_missing")
    truth = json.loads(truth_path.read_text(encoding="utf-8"))
    validate_json_document(
        truth,
        "founder_haplotype_case.schema.json",
        schemas_dir=_CALIBRATION_SCHEMAS,
    )
    summary = json.loads(summary_path.read_text(encoding="utf-8"))
    validate_json_document(summary, "founder_analysis_summary.schema.json")
    consensus_table = validate_tsv_table(consensus_path, "founder_consensus.schema.json")
    segment_table = validate_tsv_table(segments_path, "founder_segments.schema.json")
    if consensus_table.row_count != 1:
        raise ValueError("d24_consensus_row_count_invalid")
    result = evaluate_case(
        truth,
        summary,
        consensus_table.rows[0],
        list(segment_table.rows),
        minimum_independent_units=config["minimum_independent_units_per_scenario"],
        minimum_flank_markers=config["minimum_flank_markers"],
    )
    return result, _case_digest(paths)


def evaluate_founder_haplotype(
    *,
    case_roots: list[Path],
    configuration_path: Path,
    workspace_root: Path,
    run_id: str,
    created_at: str | None = None,
) -> dict[str, Any]:
    """Publie uniquement les agrégats D24, sans outil génétique ni identifiant."""
    destination = new_workspace_run(workspace_root, run_id)
    config = load_configuration(configuration_path)
    if len(case_roots) != len({path.absolute() for path in case_roots}):
        raise ValueError("d24_case_root_duplicate")
    loaded = [load_case(path, config) for path in case_roots]
    results = [item[0] for item in loaded]
    validate_case_balance(results, config)
    aggregate = aggregate_results(results, primary_profile=config["primary_profile"])
    criteria = acceptance_results(aggregate["metrics"])
    technical_status = "FAIL" if any(
        row["status"] == "FAIL" for row in criteria
    ) else "PASS"
    document = {
        "schema_version": "1.0.0",
        "method_id": METHOD_ID,
        "report_id": "CAL-REP-13",
        "evaluation_state": "CLOSED",
        "created_at": created_at or datetime.now(timezone.utc).isoformat(),
        "configuration_sha256": sha256_file(configuration_path),
        "protocol_sha256": sha256_file(_PROTOCOL),
        "source_case_sha256s": sorted(item[1] for item in loaded),
        "technical_status": technical_status,
        "performance_qualification": "REVIEW_REQUIRED",
        **aggregate,
        "criteria": criteria,
        "calibrated": False,
        "founder_effect_proven": False,
    }
    validate_json_document(
        document,
        "founder_haplotype_evaluation.schema.json",
        schemas_dir=_CALIBRATION_SCHEMAS,
    )
    temporary = Path(tempfile.mkdtemp(prefix=f".{run_id}.", dir=destination.parent))
    try:
        result_path = temporary / "CAL-REP-13.json"
        write_json(result_path, document)
        (temporary / "checksums.sha256").write_text(
            f"{sha256_file(result_path)}  CAL-REP-13.json\n",
            encoding="utf-8",
        )
        temporary.rename(destination)
    except Exception:
        shutil.rmtree(temporary)
        raise
    return document


def main() -> int:
    """CLI de l’évaluateur D24 ; aucune exécution de l’étape 13."""
    parser = argparse.ArgumentParser(description="Évaluation D24 de l’étape 13")
    parser.add_argument("--case-root", dest="case_roots", type=Path, action="append", required=True)
    parser.add_argument("--configuration-path", type=Path, required=True)
    parser.add_argument("--workspace-root", type=Path, required=True)
    parser.add_argument("--run-id", required=True)
    document = evaluate_founder_haplotype(**vars(parser.parse_args()))
    print(document["technical_status"], document["performance_qualification"])
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
