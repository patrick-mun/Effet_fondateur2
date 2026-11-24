"""Construction et audit de la release aveugle technique R12b.

R12b est indépendante du candidat R12 rejeté. Elle matérialise seulement une
partition synthétique ``VALIDATION_BLIND`` et sa vérité privée scellée. Elle ne
produit aucun run ni aucune évaluation R13. Le reçu G5 associé autorise donc
uniquement la validation d'intégrité technique, jamais une qualification
scientifique ni une promotion en production.
"""

from __future__ import annotations

import argparse
from dataclasses import dataclass
from hashlib import sha256
import json
from pathlib import Path
import shutil
import sys
from typing import Any, Mapping

from effet_fondateur.calibration.contracts import (
    STAGE_IDS,
    TRUTH_CONSUMERS,
    VIEW_CONSUMERS,
    canonical_sha256,
    file_sha256,
    load_json,
    load_yaml,
    validate_document,
)
from effet_fondateur.calibration.release_guard import (
    partition_pseudonym,
    scan_blinded_document,
    validate_blinded_inventory,
)


RELEASE_ID = "panel-master-generic-v1-r12b-validation-blind"
RELEASE_DIRECTORY = "panel_master_generic_v1_r12b_validation_blind"
PARTITION_ID = "VALIDATION_BLIND"
FIXED_TIME = "2026-11-13T08:00:00Z"
CONFIGURATION_PATH = Path("calibration/configurations/validation_blind_v1.yaml")
SOURCE_R12_RECEIPT = Path("calibration/readiness/R12_FROZEN_CANDIDATE.json")
FROZEN_COMPONENTS = (
    Path("src/effet_fondateur/calibration/r12b_validation.py"),
    Path("src/effet_fondateur/calibration/r13_runner.py"),
    Path("src/effet_fondateur/calibration/r13_evaluator.py"),
    Path("src/effet_fondateur/calibration/micro_stage_worker.py"),
    Path("src/effet_fondateur/calibration/micro_contract.py"),
    Path("calibration/configurations/validation_blind_v1.yaml"),
    Path("calibration/configurations/stage_criteria_v1.yaml"),
    Path("calibration/configurations/compatibility_matrix_v1.yaml"),
    Path("calibration/configurations/isolation_policy_v1.yaml"),
)


class R12bValidationError(RuntimeError):
    """Violation bloquante du contrat R12b/G5."""

    def __init__(self, code: str, detail: str | None = None) -> None:
        self.code = code
        self.detail = detail
        super().__init__(code if detail is None else f"{code}: {detail}")


@dataclass(frozen=True)
class R12bAudit:
    """Résumé non sensible de la release R12b auditée."""

    release_id: str
    partition_id: str
    blind_artifact_count: int
    truth_family_count: int
    view_count: int
    frozen_component_count: int
    run_count: int
    evaluation_count: int
    r07_status: str
    r13_scope: str


def _fail(code: str, detail: str | None = None) -> None:
    raise R12bValidationError(code, detail)


def _write_json(path: Path, document: Any) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(
        json.dumps(document, ensure_ascii=False, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )


def _synthetic_subjects(count: int) -> list[dict[str, str]]:
    secret = sha256(b"Effet_fondateur2|R12b|VALIDATION_BLIND|pseudonyms").digest()
    subjects: list[dict[str, str]] = []
    for index in range(1, count + 1):
        private_id = f"VAL-R12B-{index:04d}"
        subjects.append(
            {
                "sample_id": partition_pseudonym(secret, PARTITION_ID, private_id),
                "cohort_group": f"GROUP-{((index - 1) % 3) + 1}",
                "batch_id": f"VB-{((index - 1) % 4) + 1}",
            }
        )
    return subjects


def _synthetic_variants(count: int) -> list[dict[str, Any]]:
    variants: list[dict[str, Any]] = []
    for index in range(count):
        chromosome = (index % 22) + 1
        position = 1_000_000 + (index // 22) * 1_000_000 + chromosome * 10_000
        variants.append(
            {
                "chromosome": chromosome,
                "position_bp": position,
                "variant_id": f"SYN-R12B-{index + 1:04d}",
                "ref": "A" if index % 2 == 0 else "C",
                "alt": "G" if index % 2 == 0 else "T",
            }
        )
    return variants


def _view_payload(
    view_type: str,
    subjects: list[dict[str, str]],
    variants: list[dict[str, Any]],
    configuration: Mapping[str, Any],
) -> dict[str, Any]:
    selected = variants
    if view_type == "V-REGIONAL":
        selected = [item for item in variants if item["chromosome"] == 7]
    elif view_type == "V-TARGET":
        selected = [dict(configuration["target"])]
    elif view_type == "V-DENSITY":
        selected = variants[::2]
    payload: dict[str, Any] = {
        "schema_version": "1.0.0",
        "release_id": RELEASE_ID,
        "partition_id": PARTITION_ID,
        "view_type": view_type,
        "scope": "R13_TECHNICAL_BLIND_VALIDATION",
        "subjects": subjects,
        "variants": selected,
        "scientific_performance_decision_allowed": False,
    }
    if view_type == "V-SENSITIVITY":
        payload["signed_factors"] = ["BASELINE", "CHALLENGE"]
        payload["expected_direction_exposed"] = False
    scan_blinded_document(payload)
    return payload


def _materialize_blind_inputs(
    release_root: Path,
    configuration: Mapping[str, Any],
) -> tuple[dict[str, Any], dict[str, Any]]:
    blind_root = release_root / "blinded_inputs" / PARTITION_ID
    views_root = blind_root / "views"
    views_root.mkdir(parents=True)
    subjects = _synthetic_subjects(int(configuration["counts"]["study_identities"]))
    variants = _synthetic_variants(int(configuration["counts"]["variants"]))
    artifacts: list[dict[str, Any]] = []
    view_records: list[dict[str, Any]] = []
    for index, (view_type, consumers) in enumerate(VIEW_CONSUMERS.items(), start=1):
        payload = _view_payload(view_type, subjects, variants, configuration)
        path = views_root / f"{view_type.lower().replace('-', '_')}.json"
        _write_json(path, payload)
        relative = path.relative_to(blind_root).as_posix()
        digest = file_sha256(path)
        opaque_id = f"VIEW-{sha256(f'{RELEASE_ID}|{view_type}'.encode()).hexdigest()[:24]}"
        artifacts.append(
            {
                "artifact_id": f"BLIND-R12B-{index:02d}",
                "relative_path": relative,
                "sha256": digest,
                "size_bytes": path.stat().st_size,
                "media_type": "application/json",
            }
        )
        view_records.append(
            {
                "view_type": view_type,
                "opaque_view_id": opaque_id,
                "parent_sha256": digest,
                "transformation_id": f"R12B-{view_type}-V1",
                "consumer_stages": list(consumers),
                "inventory_path": (
                    f"blinded_inputs/{PARTITION_ID}/blinded_input_inventory.json"
                ),
            }
        )
    inventory = {
        "schema_version": "1.0.0",
        "inventory_id": "INVENTORY-R12B-VALIDATION-BLIND-001",
        "release_id": RELEASE_ID,
        "partition_id": PARTITION_ID,
        "opaque_view_id": f"VIEW-{sha256(RELEASE_ID.encode()).hexdigest()[:24]}",
        "artifacts": artifacts,
    }
    view_manifest = {
        "schema_version": "1.0.0",
        "release_id": RELEASE_ID,
        "partition_id": PARTITION_ID,
        "views": view_records,
    }
    _write_json(blind_root / "blinded_input_inventory.json", inventory)
    _write_json(release_root / "views/panel_view_manifest.json", view_manifest)
    return inventory, view_manifest


def _materialize_truth(release_root: Path) -> dict[str, Any]:
    truth_root = release_root / "private_truth" / PARTITION_ID
    truth_root.mkdir(parents=True)
    artifacts: list[dict[str, Any]] = []
    for index, (truth_id, consumers) in enumerate(TRUTH_CONSUMERS.items(), start=1):
        path = truth_root / f"{truth_id.lower().replace('-', '_')}.json"
        _write_json(
            path,
            {
                "schema_version": "1.0.0",
                "release_id": RELEASE_ID,
                "partition_id": PARTITION_ID,
                "truth_id": truth_id,
                "private_namespace": "VALIDATION_BLIND_R12B_V1",
                "master_identity_prefix": "VAL-R12B",
                "synthetic_only": True,
                "expected_technical_verdict": "PASS",
                "scientific_performance_truth": None,
            },
        )
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
        "disclosure_state": "SEALED_PRIVATE",
        "truth_artifacts": artifacts,
    }
    _write_json(release_root / "private_truth/master_truth_manifest.json", manifest)
    return manifest


def _component_bindings(repository_root: Path) -> list[dict[str, str]]:
    bindings: list[dict[str, str]] = []
    for relative in FROZEN_COMPONENTS:
        path = repository_root / relative
        if not path.is_file():
            _fail("frozen_component_missing", relative.as_posix())
        bindings.append({"relative_path": relative.as_posix(), "sha256": file_sha256(path)})
    return bindings


def _write_checksum_index(release_root: Path) -> None:
    checksum_path = release_root / "release_checksums.sha256"
    paths = sorted(
        path for path in release_root.rglob("*")
        if path.is_file() and path != checksum_path
    )
    checksum_path.write_text(
        "".join(
            f"{file_sha256(path)}  {path.relative_to(release_root).as_posix()}\n"
            for path in paths
        ),
        encoding="utf-8",
    )


def build_frozen_candidate(repository_root: Path, release_root: Path) -> R12bAudit:
    """Construit R12b de manière transactionnelle, sans lancer R13."""
    repository_root = repository_root.resolve(strict=True)
    release_root = release_root.resolve(strict=False)
    if release_root.exists():
        _fail("release_already_exists", str(release_root))
    staging = release_root.with_name(f".{release_root.name}.building")
    if staging.exists():
        _fail("staging_already_exists", str(staging))
    configuration = load_yaml(repository_root / CONFIGURATION_PATH)
    if configuration.get("partition") != PARTITION_ID:
        _fail("configuration_partition_mismatch")
    try:
        staging.mkdir(parents=True)
        inventory, view_manifest = _materialize_blind_inputs(staging, configuration)
        truth_manifest = _materialize_truth(staging)
        release_manifest = {
            "schema_version": "1.0.0",
            "release_id": RELEASE_ID,
            "release_design": "panel_master_generic_v1",
            "assembly": "GRCh38",
            "stage_ids": list(STAGE_IDS),
            "partitions": ["CALIBRATION", "VALIDATION_BLIND", "RESOURCE_ONLY"],
            "scenario_design": {
                "materialized_partition": PARTITION_ID,
                "scope": "R13_TECHNICAL_BLIND_VALIDATION",
                "configuration_sha256": canonical_sha256(configuration),
                "blind_inventory_sha256": canonical_sha256(inventory),
                "truth_manifest_sha256": canonical_sha256(truth_manifest),
                "scientific_performance_decision_allowed": False,
                "r07_status": "REVISION_REQUIRED",
            },
            "target": {
                "selection": "PROFILE_CONFIGURED_AUTOSOMAL",
                "hard_coded": False,
            },
            "qualification_flags": {
                "calibrated": False,
                "founder_effect_proven": False,
                "production_promotion": False,
            },
        }
        _write_json(staging / "governance/panel_release_manifest.json", release_manifest)
        source_receipt = load_json(repository_root / SOURCE_R12_RECEIPT)
        frozen_manifest = {
            "schema_version": "1.0.0",
            "release_id": RELEASE_ID,
            "candidate_phase": "R12b",
            "frozen_at": FIXED_TIME,
            "frozen": True,
            "partition_id": PARTITION_ID,
            "r13_authorization": "TECHNICAL_BLIND_VALIDATION_ONLY",
            "scientific_validation_authorized": False,
            "r07_status": "REVISION_REQUIRED",
            "r07_revision_topics": list(configuration["r07_revision_topics"]),
            "source_r12": {
                "release_id": source_receipt["release_id"],
                "candidate_state": source_receipt["candidate_state"],
                "receipt_path": SOURCE_R12_RECEIPT.as_posix(),
                "receipt_sha256": file_sha256(repository_root / SOURCE_R12_RECEIPT),
                "inputs_reused": False,
            },
            "frozen_components": _component_bindings(repository_root),
            "blind_inventory_sha256": canonical_sha256(inventory),
            "view_manifest_sha256": canonical_sha256(view_manifest),
            "truth_manifest_sha256": canonical_sha256(truth_manifest),
            "restrictions": [
                "NO_R13_RESULT_IN_R12B",
                "NO_SCIENTIFIC_PERFORMANCE_CLAIM",
                "NO_FOUNDER_EFFECT_CLAIM",
                "NO_PRODUCTION_PROMOTION",
                "NO_REAL_DATA",
                "R07_REVISION_REMAINS_REQUIRED",
                "TRUTH_OPENS_ONLY_AFTER_CLOSED_R13_RUN",
            ],
        }
        _write_json(staging / "governance/r12b_frozen_manifest.json", frozen_manifest)
        validate_document(
            release_manifest,
            repository_root / "calibration/schemas/panel_release_manifest.schema.json",
        )
        validate_document(
            view_manifest,
            repository_root / "calibration/schemas/panel_view_manifest.schema.json",
        )
        validate_document(
            inventory,
            repository_root / "calibration/schemas/blinded_input_inventory.schema.json",
        )
        validate_document(
            truth_manifest,
            repository_root / "calibration/schemas/master_truth_manifest.schema.json",
        )
        _write_checksum_index(staging)
        staging.rename(release_root)
    except Exception:
        if staging.exists():
            shutil.rmtree(staging)
        raise
    return validate_frozen_candidate(repository_root, release_root)


def _validate_checksum_index(release_root: Path) -> None:
    checksum_path = release_root / "release_checksums.sha256"
    if not checksum_path.is_file():
        _fail("release_checksum_index_missing")
    declared: dict[str, str] = {}
    for line in checksum_path.read_text(encoding="utf-8").splitlines():
        digest, separator, relative = line.partition("  ")
        if not separator or relative in declared:
            _fail("release_checksum_index_invalid", line)
        declared[relative] = digest
    actual_paths = {
        path.relative_to(release_root).as_posix(): path
        for path in release_root.rglob("*")
        if path.is_file() and path != checksum_path
    }
    if set(declared) != set(actual_paths):
        _fail("release_checksum_file_set_mismatch")
    for relative, path in actual_paths.items():
        if path.is_symlink() or path.stat(follow_symlinks=False).st_nlink != 1:
            _fail("release_link_forbidden", relative)
        if file_sha256(path) != declared[relative]:
            _fail("release_checksum_mismatch", relative)


def validate_frozen_candidate(repository_root: Path, release_root: Path) -> R12bAudit:
    """Vérifie la release, les empreintes et l'absence de toute exécution R13."""
    repository_root = repository_root.resolve(strict=True)
    release_root = release_root.resolve(strict=True)
    if release_root.is_symlink():
        _fail("release_root_symlink_forbidden")
    frozen = load_json(release_root / "governance/r12b_frozen_manifest.json")
    release_manifest = load_json(release_root / "governance/panel_release_manifest.json")
    inventory = load_json(
        release_root / f"blinded_inputs/{PARTITION_ID}/blinded_input_inventory.json"
    )
    view_manifest = load_json(release_root / "views/panel_view_manifest.json")
    truth_manifest = load_json(release_root / "private_truth/master_truth_manifest.json")
    if frozen.get("release_id") != RELEASE_ID or frozen.get("frozen") is not True:
        _fail("frozen_manifest_identity_mismatch")
    if frozen.get("r13_authorization") != "TECHNICAL_BLIND_VALIDATION_ONLY":
        _fail("r13_scope_not_technical_only")
    if frozen.get("scientific_validation_authorized") is not False:
        _fail("scientific_validation_must_remain_forbidden")
    if frozen.get("r07_status") != "REVISION_REQUIRED":
        _fail("r07_revision_not_recorded")
    validate_document(
        release_manifest,
        repository_root / "calibration/schemas/panel_release_manifest.schema.json",
    )
    validate_document(
        inventory,
        repository_root / "calibration/schemas/blinded_input_inventory.schema.json",
    )
    validate_document(
        view_manifest,
        repository_root / "calibration/schemas/panel_view_manifest.schema.json",
    )
    validate_document(
        truth_manifest,
        repository_root / "calibration/schemas/master_truth_manifest.schema.json",
    )
    if inventory.get("partition_id") != PARTITION_ID:
        _fail("blind_partition_mismatch")
    validate_blinded_inventory(
        inventory,
        view_root=release_root / f"blinded_inputs/{PARTITION_ID}",
        repository_root=repository_root,
        verify_files=True,
    )
    if canonical_sha256(inventory) != frozen.get("blind_inventory_sha256"):
        _fail("frozen_inventory_hash_mismatch")
    if canonical_sha256(view_manifest) != frozen.get("view_manifest_sha256"):
        _fail("frozen_view_hash_mismatch")
    if canonical_sha256(truth_manifest) != frozen.get("truth_manifest_sha256"):
        _fail("frozen_truth_hash_mismatch")
    truth_ids = {artifact["truth_id"] for artifact in truth_manifest["truth_artifacts"]}
    if truth_ids != set(TRUTH_CONSUMERS):
        _fail("truth_family_set_mismatch")
    for artifact in truth_manifest["truth_artifacts"]:
        path = release_root / artifact["relative_path"]
        if not path.is_file() or file_sha256(path) != artifact["sha256"]:
            _fail("private_truth_hash_mismatch", artifact["truth_id"])
    bindings = frozen.get("frozen_components")
    if not isinstance(bindings, list) or len(bindings) != len(FROZEN_COMPONENTS):
        _fail("frozen_component_set_mismatch")
    expected_paths = {path.as_posix() for path in FROZEN_COMPONENTS}
    if {item.get("relative_path") for item in bindings} != expected_paths:
        _fail("frozen_component_path_mismatch")
    for binding in bindings:
        path = repository_root / binding["relative_path"]
        if not path.is_file() or file_sha256(path) != binding.get("sha256"):
            _fail("frozen_component_hash_mismatch", binding["relative_path"])
    _validate_checksum_index(release_root)
    runs = list(release_root.glob("runs/*"))
    evaluations = list(release_root.glob("evaluations/*"))
    if runs or evaluations:
        _fail("r13_output_present_in_r12b")
    return R12bAudit(
        release_id=RELEASE_ID,
        partition_id=PARTITION_ID,
        blind_artifact_count=len(inventory["artifacts"]),
        truth_family_count=len(truth_manifest["truth_artifacts"]),
        view_count=len(view_manifest["views"]),
        frozen_component_count=len(bindings),
        run_count=0,
        evaluation_count=0,
        r07_status="REVISION_REQUIRED",
        r13_scope="TECHNICAL_BLIND_VALIDATION_ONLY",
    )


def g5_receipt_document(repository_root: Path, release_root: Path) -> dict[str, Any]:
    """Retourne le reçu G5 positif après un audit physique complet de R12b."""
    audit = validate_frozen_candidate(repository_root, release_root)
    frozen_path = release_root / "governance/r12b_frozen_manifest.json"
    checksum_path = release_root / "release_checksums.sha256"
    return {
        "schema_version": "1.0.0",
        "receipt_id": "R12B-G5-READY-panel_master_generic_v1",
        "release_design": "panel_master_generic_v1",
        "release_id": RELEASE_ID,
        "scope": "R13_TECHNICAL_BLIND_VALIDATION",
        "ready_to_run": True,
        "frozen": True,
        "gate_g5_passed": True,
        "candidate_state": "READY_TECHNICAL_VALIDATION",
        "issued_on": "2026-11-13",
        "partition_id": audit.partition_id,
        "r07_status": audit.r07_status,
        "scientific_validation_authorized": False,
        "checks": [
            "VALIDATION_BLIND_MATERIALIZED",
            "VALIDATION_IDENTITIES_USE_INDEPENDENT_NAMESPACE",
            "FOURTEEN_TRUTH_FAMILIES_SEALED",
            "ELEVEN_BLINDED_VIEWS_CLOSED",
            "BLINDED_INPUT_HASHES_VERIFIED",
            "R13_RUNNER_FROZEN",
            "R13_EVALUATOR_FROZEN",
            "TRUTH_ROUTING_AFTER_RUN_CLOSURE_ONLY",
            "NO_R13_RESULT_PRESENT",
            "NO_REAL_DATA",
            "R07_REVISION_RECORDED",
        ],
        "evidence_sha256": {
            "r12b_frozen_manifest": file_sha256(frozen_path),
            "release_checksums": file_sha256(checksum_path),
            "validation_configuration": file_sha256(repository_root / CONFIGURATION_PATH),
            "r13_runner": file_sha256(repository_root / FROZEN_COMPONENTS[1]),
            "r13_evaluator": file_sha256(repository_root / FROZEN_COMPONENTS[2]),
        },
        "qualification_flags": {
            "calibrated": False,
            "founder_effect_proven": False,
            "production_promotion": False,
        },
        "restrictions": [
            "R13_TECHNICAL_BLIND_VALIDATION_ONLY",
            "NO_SCIENTIFIC_QUALIFICATION_CLAIM",
            "NO_PRODUCTION_PROMOTION",
            "R07_MUST_BE_REVISED_BEFORE_FULL_QUALIFICATION",
            "VALIDATION_TRUTH_MUST_REMAIN_SEALED_UNTIL_RUN_CLOSURE",
            "NO_REAL_DATA",
            "NO_HTML_OR_PDF_BEFORE_R14",
        ],
    }


def _parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--repository-root", type=Path, default=Path.cwd())
    parser.add_argument("--release-root", type=Path, required=True)
    parser.add_argument("--validate-only", action="store_true")
    parser.add_argument("--print-g5-receipt", action="store_true")
    return parser


def main(argv: list[str] | None = None) -> int:
    arguments = _parser().parse_args(argv)
    if arguments.validate_only:
        audit = validate_frozen_candidate(arguments.repository_root, arguments.release_root)
    else:
        audit = build_frozen_candidate(arguments.repository_root, arguments.release_root)
    print(json.dumps(audit.__dict__, ensure_ascii=False, sort_keys=True))
    if arguments.print_g5_receipt:
        print(
            json.dumps(
                g5_receipt_document(arguments.repository_root, arguments.release_root),
                ensure_ascii=False,
                indent=2,
                sort_keys=True,
            )
        )
    return 0


if __name__ == "__main__":
    sys.exit(main())
