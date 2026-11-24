"""Préflight du nouveau candidat scientifique R12c.

Ce module ferme le contrat de construction à partir de R07, R08 et R11b et
audite les composants requis avant l'émission du reçu G5 séparé.
"""

from __future__ import annotations

import argparse
from collections.abc import Mapping
import json
import os
from pathlib import Path
import sys
from typing import Any

from effet_fondateur.calibration.contracts import (
    STAGE_IDS,
    canonical_sha256,
    file_sha256,
    load_json,
    load_yaml,
)
from effet_fondateur.calibration.r11b_closure import (
    RECEIPT_RELATIVE as R11B_RECEIPT_PATH,
    receipt_document as r11b_receipt_document,
)
from effet_fondateur.calibration.r07_revision import (
    CONFIGURATION_PATH as R07_CONFIGURATION_PATH,
)
from effet_fondateur.calibration.r13_scientific_runner import (
    RUNNER_IMPLEMENTATION_STATE,
)
from effet_fondateur.calibration.r13_scientific_evaluator import (
    METRIC_RULE_MANIFEST_PATH,
    SCIENTIFIC_EVALUATOR_STATE,
    load_metric_rule_bundle,
    validate_metric_rule_manifest,
)
from effet_fondateur.calibration.r13_scientific_call_worker import (
    WORKER_STATE,
    worker_contract_audit,
)
from effet_fondateur.calibration.r13_scientific_materializer import (
    MATERIALIZER_STATE,
    materializer_contract_audit,
)
from effet_fondateur.calibration.scientific_metric_producers import (
    audit_metric_producer_methods,
)
from effet_fondateur.calibration.scientific_stage_producers import (
    stage_producer_audit,
)


CONFIGURATION_PATH = Path("calibration/configurations/scientific_validation_v1.yaml")
CRITERIA_PATH = Path("calibration/configurations/stage_criteria_v1.yaml")
PREFLIGHT_PATH = Path("calibration/readiness/R12C_SCIENTIFIC_PREFLIGHT.json")
EXPECTED_FLAGS = {
    "calibrated": False,
    "founder_effect_proven": False,
    "production_promotion": False,
}


class R12cScientificError(RuntimeError):
    """Violation bloquante du contrat scientifique R12c."""

    def __init__(self, code: str, detail: str | None = None) -> None:
        self.code = code
        self.detail = detail
        super().__init__(code if detail is None else f"{code}: {detail}")


def _fail(code: str, detail: str | None = None) -> None:
    raise R12cScientificError(code, detail)


def preflight_document(repository_root: Path) -> dict[str, Any]:
    """Valide le dessin R12c et décrit explicitement le travail restant."""
    repository_root = repository_root.resolve(strict=True)
    configuration = load_yaml(repository_root / CONFIGURATION_PATH)
    dimensioning = load_yaml(repository_root / R07_CONFIGURATION_PATH)
    criteria = load_yaml(repository_root / CRITERIA_PATH)
    tracked_r11b = load_json(repository_root / R11B_RECEIPT_PATH)
    if tracked_r11b != r11b_receipt_document(repository_root):
        _fail("r12c_r11b_receipt_drift")
    if not tracked_r11b.get("new_scientific_release_construction_authorized"):
        _fail("r12c_construction_not_authorized")
    if configuration.get("partition") != "VALIDATION_BLIND":
        _fail("r12c_partition_mismatch")
    if configuration.get("scope") != "R13_SCIENTIFIC_BLIND_VALIDATION":
        _fail("r12c_scope_mismatch")
    if configuration.get("qualification_flags") != EXPECTED_FLAGS:
        _fail("r12c_qualification_flags_mismatch")
    if tuple(criteria.get("stages", {})) != STAGE_IDS:
        _fail("r12c_stage_criteria_mismatch")

    minima = configuration.get("scientific_minima", {})
    r07_minima = dimensioning.get("scientific_minima", {})
    expected_minima = {
        "study_identities": r07_minima["study_identities_per_partition"],
        "phasing_reference": r07_minima["phasing_reference_per_partition"],
        "ancestry_reference": r07_minima["ancestry_reference_per_partition"],
        "null_reference": r07_minima["null_reference_per_partition"],
        "primary_independent_units": r07_minima["primary_independent_units"],
        "primary_negative_units": r07_minima["primary_negative_units"],
        "founder_events_per_age": r07_minima["founder_events_per_age"],
        "external_tool_seeds": r07_minima["external_tool_seed_minimum"],
        "autosomes": r07_minima["autosomes"],
    }
    if minima != expected_minima:
        _fail("r12c_scientific_minima_mismatch")
    if configuration.get("observed_densities") != r07_minima["observed_densities"]:
        _fail("r12c_density_mismatch")
    if configuration.get("phasing_schedule") != dimensioning["phasing_schedule"][
        "expected_calls"
    ]:
        _fail("r12c_phasing_schedule_mismatch")
    ibd = configuration.get("ibd_schedule", {})
    r07_ibd = dimensioning["ibd_worlds"]
    if (
        ibd.get("independent_worlds") != r07_ibd["expected_independent_worlds"]["total"]
        or ibd.get("tool_density_calls") != r07_ibd["expected_tool_density_calls"]
        or ibd.get("tools") != r07_ibd["tools"]
        or ibd.get("densities") != r07_ibd["densities"]
    ):
        _fail("r12c_ibd_schedule_mismatch")
    execution_contract = configuration.get("execution_contract", {})
    schedule_check = next(
        (
            check
            for check in tracked_r11b.get("checks", [])
            if check.get("check") == "QUALIFICATION_SCHEDULE_DAYS_AT_24_WORKERS"
        ),
        None,
    )
    if not isinstance(schedule_check, Mapping):
        _fail("r12c_r11b_schedule_projection_missing")
    projected_days_at_eight = float(schedule_check["projected"]) * 24 / 8
    if (
        execution_contract.get("required_workers") != 8
        or execution_contract.get("scheduling")
        != "BATCHED_RESUMABLE_NO_REDUCED_SCIENTIFIC_VOLUME"
        or execution_contract.get("local_eight_cpu_host_authorized_for_full_campaign")
        is not True
        or abs(
            float(execution_contract.get("projected_calendar_days", -1))
            - projected_days_at_eight
        )
        > 1e-12
        or projected_days_at_eight
        > float(dimensioning["resources"]["qualification_calendar_days_max"])
    ):
        _fail("r12c_eight_worker_schedule_invalid")

    metric_rule_manifest = repository_root / METRIC_RULE_MANIFEST_PATH
    metric_rule_audit = (
        validate_metric_rule_manifest(repository_root)
        if metric_rule_manifest.exists()
        else None
    )
    producer_method_audit = (
        audit_metric_producer_methods(load_metric_rule_bundle(repository_root))
        if metric_rule_manifest.exists()
        else None
    )
    required_components = {
        "scientific_runner": Path(
            "src/effet_fondateur/calibration/r13_scientific_runner.py"
        ),
        "scientific_evaluator": Path(
            "src/effet_fondateur/calibration/r13_scientific_evaluator.py"
        ),
        "scientific_materializer": Path(
            "src/effet_fondateur/calibration/r13_scientific_materializer.py"
        ),
        "single_call_worker": Path(
            "src/effet_fondateur/calibration/r13_scientific_call_worker.py"
        ),
        "materialized_release": Path(
            "calibration_workspace/releases/"
            "panel_master_generic_v1_r12c_scientific_validation"
        ),
        "scientific_metric_rules": METRIC_RULE_MANIFEST_PATH,
    }
    missing = []
    for name, relative in required_components.items():
        if not (repository_root / relative).exists():
            missing.append(name)
        elif (
            name == "scientific_runner"
            and RUNNER_IMPLEMENTATION_STATE != "FULL_EXECUTION_READY"
        ):
            missing.append(name)
        elif (
            name == "scientific_evaluator"
            and SCIENTIFIC_EVALUATOR_STATE != "FULL_EVALUATION_READY"
        ):
            missing.append(name)
        elif (
            name == "scientific_materializer"
            and MATERIALIZER_STATE != "FULL_CALL_MATERIALIZATION_READY"
        ):
            missing.append(name)
        elif (
            name == "single_call_worker"
            and WORKER_STATE != "FULL_SINGLE_CALL_EXECUTION_READY"
        ):
            missing.append(name)
    return {
        "schema_version": "1.0.0",
        "receipt_id": "R12C-SCIENTIFIC-PREFLIGHT-001",
        "release_id": configuration["release_id"],
        "scope": configuration["scope"],
        "partition_id": configuration["partition"],
        "status": "DESIGN_READY_COMPONENTS_REQUIRED" if missing else "READY_FOR_G5_AUDIT",
        "design_checks": [
            "R07_MINIMA_PRESERVED",
            "R08_ALL_25_STAGES_BOUND",
            "R11B_RESOURCE_READY_BOUND",
            "R11B_PROJECTION_REBATCHED_TO_8_WORKERS_UNDER_120_DAYS",
            "VALIDATION_BLIND_NAMESPACE_INDEPENDENT",
            "TARGET_PROFILE_CONFIGURED_NOT_HARD_CODED",
            "ALL_1180_EXTERNAL_CALLS_MATERIALIZABLE",
            "SINGLE_CALL_WORKER_TOOLCHAIN_BOUND",
            "PRIVATE_TRUTH_SEPARATED_FROM_WORKER",
        ],
        "missing_before_g5": missing,
        "gate_g5_passed": False,
        "scientific_validation_authorized": False,
        "current_host_logical_cpus": os.cpu_count() or 1,
        "required_workers": configuration["execution_contract"]["required_workers"],
        "projected_calendar_days": projected_days_at_eight,
        "qualification_flags": EXPECTED_FLAGS,
        "metric_rule_bundle": metric_rule_audit,
        "metric_producer_methods": producer_method_audit,
        "stage_metric_producers": stage_producer_audit(),
        "scientific_materializer": materializer_contract_audit(repository_root),
        "single_call_worker": worker_contract_audit(repository_root),
        "evidence_sha256": {
            CONFIGURATION_PATH.as_posix(): file_sha256(
                repository_root / CONFIGURATION_PATH
            ),
            R07_CONFIGURATION_PATH.as_posix(): file_sha256(
                repository_root / R07_CONFIGURATION_PATH
            ),
            CRITERIA_PATH.as_posix(): file_sha256(repository_root / CRITERIA_PATH),
            R11B_RECEIPT_PATH.as_posix(): file_sha256(
                repository_root / R11B_RECEIPT_PATH
            ),
            **(
                {
                    METRIC_RULE_MANIFEST_PATH.as_posix(): file_sha256(
                        metric_rule_manifest
                    )
                }
                if metric_rule_manifest.exists()
                else {}
            ),
        },
        "configuration_sha256": canonical_sha256(configuration),
        "next_phase": "BUILD_AND_FREEZE_EXECUTABLE_R12C_RELEASE",
    }


def write_preflight(repository_root: Path) -> dict[str, Any]:
    """Écrit le préflight suivi, sans créer de release ni de résultat R13."""
    repository_root = repository_root.resolve(strict=True)
    document = preflight_document(repository_root)
    path = repository_root / PREFLIGHT_PATH
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(
        json.dumps(document, ensure_ascii=False, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    return document


def _parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--repository-root", type=Path, default=Path.cwd())
    parser.add_argument("--write", action="store_true")
    return parser


def main(argv: list[str] | None = None) -> int:
    arguments = _parser().parse_args(argv)
    document = (
        write_preflight(arguments.repository_root)
        if arguments.write
        else preflight_document(arguments.repository_root)
    )
    print(json.dumps(document, ensure_ascii=False, sort_keys=True))
    return 0


if __name__ == "__main__":
    sys.exit(main())
