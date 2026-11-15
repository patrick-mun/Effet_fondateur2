"""Matérialisation et gel de la release scientifique R12c.

La release contient le registre pseudonymisé, les catalogues de scénarios et
les recettes déterministes nécessaires à R13. Les matrices génotypiques
volumineuses restent des shards reproductibles générés par les work items R13 ;
leurs graines privées et leurs contrats sont scellés ici. Aucun résultat R13
n'est créé pendant R12c.
"""

from __future__ import annotations

import argparse
from hashlib import sha256
import json
from pathlib import Path
import shutil
import sys
from typing import Any, Mapping

from effet_fondateur.calibration.contracts import (
    STAGE_IDS,
    TRUTH_CONSUMERS,
    canonical_sha256,
    file_sha256,
    load_json,
    load_yaml,
)
from effet_fondateur.calibration.r13_scientific_evaluator import (
    SCIENTIFIC_EVALUATOR_STATE,
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
from effet_fondateur.calibration.r13_scientific_runner import (
    RUNNER_IMPLEMENTATION_STATE,
    build_external_call_schedule,
)
from effet_fondateur.calibration.release_guard import partition_pseudonym
from effet_fondateur.calibration.scientific_stage_producers import (
    stage_producer_audit,
)


RELEASE_ID = "panel-master-generic-v1-r12c-scientific-validation"
RELEASE_DIRECTORY = "panel_master_generic_v1_r12c_scientific_validation"
PARTITION_ID = "VALIDATION_BLIND"
CONFIGURATION_PATH = Path("calibration/configurations/scientific_validation_v1.yaml")
G5_RECEIPT_PATH = Path("calibration/readiness/R12C_G5_READY.json")
FIXED_TIME = "2026-11-14T12:00:00Z"
FREEZE_REVISION = "R12C-EXECUTABLE-REFREEZE-002"
SUPERSEDED_G5_SHA256 = "bfcac60accf286dbb8956bf420c766f0372c3da1ade7d715af9b2b224130efa8"
FROZEN_COMPONENTS = (
    Path("src/effet_fondateur/calibration/r12c_release.py"),
    Path("src/effet_fondateur/calibration/r13_scientific_runner.py"),
    Path("src/effet_fondateur/calibration/r13_checkpoint_ledger.py"),
    Path("src/effet_fondateur/calibration/r13_scientific_materializer.py"),
    Path("src/effet_fondateur/calibration/r13_scientific_call_worker.py"),
    Path("src/effet_fondateur/calibration/r13_scientific_evaluator.py"),
    Path("src/effet_fondateur/phasing/shapeit5.py"),
    Path("src/effet_fondateur/calibration/scientific_statistics.py"),
    Path("src/effet_fondateur/calibration/scientific_metric_producers.py"),
    Path("src/effet_fondateur/calibration/scientific_stage_producers.py"),
    Path("calibration/configurations/scientific_validation_v1.yaml"),
    Path("calibration/configurations/scientific_materialization_v1.yaml"),
    Path("calibration/configurations/scientific_metric_rules_v1.yaml"),
    Path("calibration/configurations/stage_criteria_v1.yaml"),
)


class R12cReleaseError(RuntimeError):
    """Release scientifique absente, altérée ou incomplète."""

    def __init__(self, code: str, detail: str | None = None) -> None:
        self.code = code
        self.detail = detail
        super().__init__(code if detail is None else f"{code}: {detail}")


def _fail(code: str, detail: str | None = None) -> None:
    raise R12cReleaseError(code, detail)


def _write_json(path: Path, document: Any) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(
        json.dumps(document, ensure_ascii=False, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )


def _component_bindings(repository_root: Path) -> list[dict[str, str]]:
    bindings = []
    for relative in FROZEN_COMPONENTS:
        path = repository_root / relative
        if not path.is_file():
            _fail("r12c_frozen_component_missing", relative.as_posix())
        bindings.append({"relative_path": relative.as_posix(), "sha256": file_sha256(path)})
    return bindings


def _subject_registry(count: int) -> list[dict[str, str]]:
    secret = sha256(b"Effet_fondateur2|R12c|VALIDATION_BLIND|pseudonyms").digest()
    return [
        {
            "sample_id": partition_pseudonym(
                secret,
                PARTITION_ID,
                f"VAL-R12C-{index:04d}",
            ),
            "identity_slot": f"SLOT-{index:04d}",
        }
        for index in range(1, count + 1)
    ]


def _truth_manifest(release_root: Path) -> dict[str, Any]:
    truth_root = release_root / "private_truth" / PARTITION_ID
    artifacts = []
    for index, (truth_id, consumers) in enumerate(TRUTH_CONSUMERS.items(), start=1):
        seed = sha256(f"{RELEASE_ID}|{truth_id}|private".encode()).hexdigest()
        document = {
            "schema_version": "1.0.0",
            "release_id": RELEASE_ID,
            "partition_id": PARTITION_ID,
            "truth_id": truth_id,
            "private_namespace": "VALIDATION_BLIND_R12C_V1",
            "generator_seed_sha256": seed,
            "synthetic_only": True,
            "disclosure": "POST_CLOSURE_EVALUATOR_ONLY",
        }
        path = truth_root / f"truth_{index:02d}.json"
        _write_json(path, document)
        artifacts.append(
            {
                "truth_id": truth_id,
                "relative_path": path.relative_to(release_root).as_posix(),
                "sha256": file_sha256(path),
                "authorized_evaluators": [f"EVAL-{stage}" for stage in consumers],
            }
        )
    manifest = {
        "schema_version": "1.0.0",
        "release_id": RELEASE_ID,
        "partition_id": PARTITION_ID,
        "disclosure_state": "SEALED_PRIVATE",
        "truth_artifacts": artifacts,
    }
    _write_json(release_root / "private_truth/master_truth_manifest.json", manifest)
    return manifest


def _write_checksum_index(release_root: Path) -> None:
    checksum_path = release_root / "release_checksums.sha256"
    files = sorted(
        path for path in release_root.rglob("*")
        if path.is_file() and path != checksum_path
    )
    checksum_path.write_text(
        "".join(
            f"{file_sha256(path)}  {path.relative_to(release_root).as_posix()}\n"
            for path in files
        ),
        encoding="utf-8",
    )


def build_release(repository_root: Path, release_root: Path) -> dict[str, Any]:
    """Construit transactionnellement R12c sans lancer R13."""
    repository_root = repository_root.resolve(strict=True)
    release_root = release_root.resolve(strict=False)
    if release_root.exists():
        _fail("r12c_release_already_exists", str(release_root))
    staging = release_root.with_name(f".{release_root.name}.building")
    if staging.exists():
        _fail("r12c_staging_already_exists", str(staging))
    configuration = load_yaml(repository_root / CONFIGURATION_PATH)
    if configuration.get("release_id") != RELEASE_ID:
        _fail("r12c_configuration_release_mismatch")
    if configuration.get("execution_contract", {}).get("required_workers") != 8:
        _fail("r12c_eight_worker_contract_missing")
    metric_audit = validate_metric_rule_manifest(repository_root)
    producer_audit = stage_producer_audit()
    materializer_audit = materializer_contract_audit(repository_root)
    worker_audit = worker_contract_audit(repository_root)
    if not producer_audit["stage_producers_complete"]:
        _fail("r12c_stage_producers_incomplete")
    if RUNNER_IMPLEMENTATION_STATE != "FULL_EXECUTION_READY":
        _fail("r12c_runner_not_frozen")
    if SCIENTIFIC_EVALUATOR_STATE != "FULL_EVALUATION_READY":
        _fail("r12c_evaluator_not_frozen")
    if (
        MATERIALIZER_STATE != "FULL_CALL_MATERIALIZATION_READY"
        or materializer_audit.get("supported_call_count") != 1_180
    ):
        _fail("r12c_materializer_not_frozen")
    if (
        WORKER_STATE != "FULL_SINGLE_CALL_EXECUTION_READY"
        or worker_audit.get("supported_call_count") != 1_180
    ):
        _fail("r12c_single_call_worker_not_frozen")
    try:
        staging.mkdir(parents=True)
        subjects = _subject_registry(configuration["scientific_minima"]["study_identities"])
        registry = {
            "schema_version": "1.0.0",
            "release_id": RELEASE_ID,
            "partition_id": PARTITION_ID,
            "synthetic_only": True,
            "identity_count": len(subjects),
            "subjects": subjects,
        }
        registry_path = staging / "blinded_inputs" / PARTITION_ID / "study_registry.json"
        _write_json(registry_path, registry)
        schedule = build_external_call_schedule()
        scenario_catalog = {
            "schema_version": "1.0.0",
            "release_id": RELEASE_ID,
            "partition_id": PARTITION_ID,
            "scientific_minima": configuration["scientific_minima"],
            "observed_densities": configuration["observed_densities"],
            "phasing_schedule": configuration["phasing_schedule"],
            "ibd_schedule": configuration["ibd_schedule"],
            "external_call_inventory_sha256": canonical_sha256(schedule),
            "materialization_mode": "DETERMINISTIC_SHARDS_AT_R13_WORK_ITEM_START",
            "scientific_volume_reduction_allowed": False,
        }
        catalog_path = staging / "blinded_inputs" / PARTITION_ID / "scenario_catalog.json"
        _write_json(catalog_path, scenario_catalog)
        truth_manifest = _truth_manifest(staging)
        execution_contract = {
            "schema_version": "1.0.0",
            "release_id": RELEASE_ID,
            "partition_id": PARTITION_ID,
            "materializer": materializer_audit,
            "single_call_worker": worker_audit,
            "call_count": 1_180,
            "materialization_timing": "ON_DEMAND_BEFORE_CHECKPOINT_CLAIM",
            "truth_available_to_worker": False,
        }
        execution_contract_path = staging / "governance/r13_execution_contract.json"
        _write_json(execution_contract_path, execution_contract)
        frozen_manifest = {
            "schema_version": "1.0.0",
            "release_id": RELEASE_ID,
            "partition_id": PARTITION_ID,
            "phase": "R12c",
            "freeze_revision": FREEZE_REVISION,
            "frozen_at": FIXED_TIME,
            "frozen": True,
            "scientific_validation_authorized": True,
            "r13_result_present": False,
            "execution_workers": 8,
            "projected_calendar_days": configuration["execution_contract"][
                "projected_calendar_days"
            ],
            "registry_sha256": file_sha256(registry_path),
            "scenario_catalog_sha256": file_sha256(catalog_path),
            "truth_manifest_sha256": canonical_sha256(truth_manifest),
            "metric_rule_bundle_sha256": metric_audit["expanded_bundle_sha256"],
            "r13_execution_contract_sha256": file_sha256(execution_contract_path),
            "materializer_state": MATERIALIZER_STATE,
            "single_call_worker_state": WORKER_STATE,
            "frozen_components": _component_bindings(repository_root),
            "qualification_flags": configuration["qualification_flags"],
            "restrictions": [
                "NO_R13_RESULT_IN_R12C",
                "NO_REAL_DATA",
                "NO_SCIENTIFIC_CLAIM_BEFORE_R13",
                "TRUTH_OPENS_ONLY_AFTER_CLOSED_RUN",
                "NO_REDUCED_SCIENTIFIC_VOLUME",
            ],
        }
        _write_json(staging / "governance/r12c_frozen_manifest.json", frozen_manifest)
        release_manifest = {
            "schema_version": "1.0.0",
            "release_id": RELEASE_ID,
            "partition_id": PARTITION_ID,
            "stage_ids": list(STAGE_IDS),
            "scope": configuration["scope"],
            "status": "FROZEN_READY_FOR_R13_SCIENTIFIC_VALIDATION",
            "configuration_sha256": canonical_sha256(configuration),
            "frozen_manifest_sha256": canonical_sha256(frozen_manifest),
        }
        _write_json(staging / "governance/release_manifest.json", release_manifest)
        _write_checksum_index(staging)
        staging.rename(release_root)
    except Exception:
        if staging.exists():
            shutil.rmtree(staging)
        raise
    return validate_release(repository_root, release_root)


def validate_release(repository_root: Path, release_root: Path) -> dict[str, Any]:
    """Valide l'intégrité et l'absence de résultat R13 dans R12c."""
    repository_root = repository_root.resolve(strict=True)
    release_root = release_root.resolve(strict=True)
    checksum_path = release_root / "release_checksums.sha256"
    if not checksum_path.is_file():
        _fail("r12c_checksum_index_missing")
    lines = [line for line in checksum_path.read_text(encoding="utf-8").splitlines() if line]
    for line in lines:
        digest, relative = line.split("  ", 1)
        path = release_root / relative
        if not path.is_file() or file_sha256(path) != digest:
            _fail("r12c_release_checksum_mismatch", relative)
    registry = load_json(
        release_root / "blinded_inputs" / PARTITION_ID / "study_registry.json"
    )
    if registry.get("identity_count") != 4_000 or len(registry.get("subjects", [])) != 4_000:
        _fail("r12c_registry_count_mismatch")
    sample_ids = [subject.get("sample_id") for subject in registry["subjects"]]
    if len(set(sample_ids)) != 4_000 or any(
        not isinstance(sample_id, str) or not sample_id.startswith("SYN-")
        for sample_id in sample_ids
    ):
        _fail("r12c_registry_pseudonym_invalid")
    truth = load_json(release_root / "private_truth/master_truth_manifest.json")
    if truth.get("disclosure_state") != "SEALED_PRIVATE" or len(
        truth.get("truth_artifacts", [])
    ) != len(TRUTH_CONSUMERS):
        _fail("r12c_truth_not_sealed")
    frozen = load_json(release_root / "governance/r12c_frozen_manifest.json")
    if frozen.get("r13_result_present") is not False:
        _fail("r12c_r13_result_flag_invalid")
    if (release_root / "runs").exists() or (release_root / "evaluations").exists():
        _fail("r12c_r13_output_present")
    execution_contract_path = release_root / "governance/r13_execution_contract.json"
    execution_contract = load_json(execution_contract_path)
    if (
        frozen.get("r13_execution_contract_sha256")
        != file_sha256(execution_contract_path)
        or frozen.get("freeze_revision") != FREEZE_REVISION
        or frozen.get("materializer_state") != MATERIALIZER_STATE
        or frozen.get("single_call_worker_state") != WORKER_STATE
        or execution_contract.get("call_count") != 1_180
        or execution_contract.get("truth_available_to_worker") is not False
    ):
        _fail("r12c_execution_contract_invalid")
    for binding in frozen.get("frozen_components", []):
        relative = binding.get("relative_path")
        if not isinstance(relative, str) or file_sha256(repository_root / relative) != binding.get("sha256"):
            _fail("r12c_frozen_component_drift", str(relative))
    return {
        "release_id": RELEASE_ID,
        "partition_id": PARTITION_ID,
        "identity_count": 4_000,
        "truth_family_count": len(TRUTH_CONSUMERS),
        "checksum_count": len(lines),
        "stage_count": len(STAGE_IDS),
        "r13_result_count": 0,
        "execution_workers": 8,
        "materializable_call_count": 1_180,
    }


def g5_receipt_document(repository_root: Path, release_root: Path) -> dict[str, Any]:
    """Produit le reçu G5 scientifique à partir de la release auditée."""
    audit = validate_release(repository_root, release_root)
    frozen_path = release_root / "governance/r12c_frozen_manifest.json"
    checksums_path = release_root / "release_checksums.sha256"
    return {
        "schema_version": "1.0.0",
        "receipt_id": "R12C-SCIENTIFIC-G5-READY",
        "receipt_revision": 2,
        "freeze_revision": FREEZE_REVISION,
        "supersedes_incomplete_g5_sha256": SUPERSEDED_G5_SHA256,
        "correction_reason": "ADD_FROZEN_SCIENTIFIC_MATERIALIZER_AND_SINGLE_CALL_WORKER",
        "issued_on": "2026-11-14",
        "release_id": RELEASE_ID,
        "partition_id": PARTITION_ID,
        "scope": "R13_SCIENTIFIC_BLIND_VALIDATION",
        "status": "READY_FOR_R13_SCIENTIFIC_VALIDATION",
        "gate_g5_passed": True,
        "scientific_validation_authorized": True,
        "scientific_validation_executed": False,
        "execution_workers": 8,
        "scientific_volume_reduced": False,
        "materializer_state": MATERIALIZER_STATE,
        "single_call_worker_state": WORKER_STATE,
        "stage_producer_count": 25,
        "metric_rule_count": 135,
        "audit": audit,
        "evidence_sha256": {
            "r12c_frozen_manifest": file_sha256(frozen_path),
            "release_checksums": file_sha256(checksums_path),
            CONFIGURATION_PATH.as_posix(): file_sha256(repository_root / CONFIGURATION_PATH),
        },
        "qualification_flags": {
            "calibrated": False,
            "founder_effect_proven": False,
            "production_promotion": False,
        },
        "restrictions": [
            "R13_MUST_CLOSE_BEFORE_TRUTH_ACCESS",
            "NO_SCIENTIFIC_CLAIM_BEFORE_POST_CLOSURE_EVALUATION",
            "NO_PRODUCTION_PROMOTION",
        ],
        "next_phase": "R13_SCIENTIFIC_BLIND_VALIDATION",
    }


def write_g5_receipt(repository_root: Path, release_root: Path) -> dict[str, Any]:
    document = g5_receipt_document(repository_root, release_root)
    path = repository_root.resolve(strict=True) / G5_RECEIPT_PATH
    _write_json(path, document)
    return document


def _parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--repository-root", type=Path, default=Path.cwd())
    parser.add_argument("--release-root", type=Path, required=True)
    parser.add_argument("--build", action="store_true")
    parser.add_argument("--write-g5", action="store_true")
    return parser


def main(argv: list[str] | None = None) -> int:
    arguments = _parser().parse_args(argv)
    if arguments.build:
        document = build_release(arguments.repository_root, arguments.release_root)
    elif arguments.write_g5:
        document = write_g5_receipt(arguments.repository_root, arguments.release_root)
    else:
        document = validate_release(arguments.repository_root, arguments.release_root)
    print(json.dumps(document, ensure_ascii=False, sort_keys=True))
    return 0


if __name__ == "__main__":
    sys.exit(main())
