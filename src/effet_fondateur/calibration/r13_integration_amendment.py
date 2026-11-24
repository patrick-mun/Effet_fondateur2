"""Audit et clôture de l'amendement R13 limité à l'intégration des outils.

L'amendement conserve le plan R12c et son registre de 1 180 appels comme
historique. Il retient un checkpoint SHAPEIT5 grandeur réelle et les contrôles
positifs IBD R12 comme preuves locales d'intégration. Les publications servent
uniquement de références à la validité intrinsèque des méthodes : aucune
performance R07/R08 ni conclusion d'effet fondateur n'est déduite ici.
"""

from __future__ import annotations

import argparse
from datetime import datetime, timezone
import json
from pathlib import Path
import re
import subprocess
import sys
from typing import Any, Callable, Mapping, Sequence

from effet_fondateur.calibration.contracts import (
    canonical_sha256,
    file_sha256,
    load_json,
    load_yaml,
)
from effet_fondateur.calibration.r13_checkpoint_ledger import (
    AMENDMENT_RELATIVE,
    LEDGER_RELATIVE,
    checkpoint_status,
    load_checkpoint_ledger,
)


CONFIGURATION_PATH = Path("calibration/configurations/r13_integration_amendment_v1.yaml")
G5_RECEIPT_PATH = Path("calibration/readiness/R12C_G5_READY.json")
READINESS_RECEIPT_PATH = Path("calibration/readiness/R13_INTEGRATION_AMENDMENT.json")
R12C_MANIFEST_RELATIVE = Path("governance/r12c_frozen_manifest.json")
R12C_CONTRACT_RELATIVE = Path("governance/r13_execution_contract.json")
DISPOSITION_RELATIVE = Path("governance/r13_integration_amendment.md")
AMENDMENT_ID = "R13-LIGHTWEIGHT-INTEGRATION-001"
RELEASE_ID = "panel-master-generic-v1-r12c-scientific-validation"
RUN_ID = "RUN-R13-SCIENTIFIC-001"
RETAINED_CALL_ID = "SHAPEIT-DMAX-CHR01-SEED01"


class R13IntegrationAmendmentError(RuntimeError):
    """Preuve absente, altérée ou incompatible avec l'amendement."""

    def __init__(self, code: str, detail: str | None = None) -> None:
        self.code = code
        self.detail = detail
        super().__init__(code if detail is None else f"{code}: {detail}")


def _fail(code: str, detail: str | None = None) -> None:
    raise R13IntegrationAmendmentError(code, detail)


def _timestamp() -> str:
    return datetime.now(timezone.utc).replace(microsecond=0).isoformat().replace("+00:00", "Z")


def _write_json(path: Path, document: Any) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(
        json.dumps(document, ensure_ascii=False, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )


def _path_within(root: Path, relative: str, *, code: str) -> Path:
    if not isinstance(relative, str):
        _fail(code, "relative_path_missing")
    try:
        path = (root / relative).resolve(strict=True)
        path.relative_to(root)
    except (FileNotFoundError, ValueError) as error:
        raise R13IntegrationAmendmentError(code, relative) from error
    if not path.is_file():
        _fail(code, relative)
    return path


def _require_equal(actual: Any, expected: Any, code: str) -> None:
    if actual != expected:
        _fail(code, f"expected={expected!r} actual={actual!r}")


def _parse_shapeit_log(log_path: Path) -> dict[str, Any]:
    text = log_path.read_text(encoding="utf-8")
    patterns = {
        "shapeit_version": r"\* Version\s+: ([^ /]+)",
        "variant_count": r"Variants \[#sites=(\d+) / region=1\]",
        "study_samples": r"Samples \[#target=(\d+) / #reference=(\d+)\]",
        "trios": r"#trios = (\d+)",
        "mendel_errors": r"%mendel_errors = [0-9.]+% \(n=(\d+)\)",
        "threads": r"\* Threads : (\d+) threads",
        "seed": r"\* Seed\s+: (\d+)",
        "output_dimensions": r"VCF/BCF writing \[N=(\d+) / L=(\d+)\]",
        "total_seconds": r"Total running time = (\d+) seconds",
    }
    matches: dict[str, re.Match[str]] = {}
    for name, pattern in patterns.items():
        match = re.search(pattern, text)
        if match is None:
            _fail("r13_amendment_shapeit_log_field_missing", name)
        matches[name] = match
    return {
        "shapeit_version": matches["shapeit_version"].group(1),
        "variant_count": int(matches["variant_count"].group(1)),
        "study_samples": int(matches["study_samples"].group(1)),
        "reference_samples": int(matches["study_samples"].group(2)),
        "trios": int(matches["trios"].group(1)),
        "mendel_errors": int(matches["mendel_errors"].group(1)),
        "threads": int(matches["threads"].group(1)),
        "seed": int(matches["seed"].group(1)),
        "output_samples": int(matches["output_dimensions"].group(1)),
        "output_variants": int(matches["output_dimensions"].group(2)),
        "shapeit_reported_seconds": int(matches["total_seconds"].group(1)),
    }


CommandRunner = Callable[..., subprocess.CompletedProcess[str]]


def _bcf_dimensions(
    bcf_path: Path,
    *,
    bcftools: str = "bcftools",
    command_runner: CommandRunner = subprocess.run,
) -> dict[str, int]:
    """Contrôle le BCF avec bcftools sans lire ni exposer les génotypes."""
    commands: Sequence[Sequence[str]] = (
        (bcftools, "query", "-l", str(bcf_path)),
        (bcftools, "index", "-n", str(bcf_path)),
    )
    outputs: list[str] = []
    for command in commands:
        try:
            result = command_runner(
                list(command),
                check=True,
                capture_output=True,
                text=True,
            )
        except (OSError, subprocess.CalledProcessError) as error:
            raise R13IntegrationAmendmentError(
                "r13_amendment_bcf_inspection_failed", command[1]
            ) from error
        outputs.append(result.stdout)
    try:
        variants = int(outputs[1].strip())
    except ValueError as error:
        raise R13IntegrationAmendmentError(
            "r13_amendment_bcf_variant_count_invalid", outputs[1].strip()
        ) from error
    return {
        "samples": len([line for line in outputs[0].splitlines() if line]),
        "variants": variants,
    }


def _validate_identity(configuration: Mapping[str, Any]) -> None:
    _require_equal(configuration.get("schema_version"), "1.0.0", "r13_amendment_schema_mismatch")
    _require_equal(configuration.get("amendment_id"), AMENDMENT_ID, "r13_amendment_id_mismatch")
    _require_equal(configuration.get("release_id"), RELEASE_ID, "r13_amendment_release_mismatch")
    _require_equal(configuration.get("run_id"), RUN_ID, "r13_amendment_run_mismatch")
    _require_equal(configuration.get("scope"), "R13_TOOL_INTEGRATION_VERIFICATION", "r13_amendment_scope_mismatch")
    campaign = configuration.get("original_campaign", {})
    _require_equal(campaign.get("call_count"), 1_180, "r13_amendment_original_count_mismatch")
    decision = configuration.get("decision", {})
    _require_equal(decision.get("retained_completed_calls"), 1, "r13_amendment_retained_count_mismatch")
    _require_equal(decision.get("withdrawn_unexecuted_calls"), 1_179, "r13_amendment_withdrawn_count_mismatch")
    _require_equal(decision.get("newly_authorized_calls"), 0, "r13_amendment_authorized_count_mismatch")
    _require_equal(decision.get("truth_opening_authorized"), False, "r13_amendment_truth_rule_invalid")
    if any(decision.get("qualification_flags", {}).values()):
        _fail("r13_amendment_qualification_flag_true")


def _validate_original_governance(
    repository_root: Path,
    release_root: Path,
) -> dict[str, str]:
    g5_path = repository_root / G5_RECEIPT_PATH
    manifest_path = release_root / R12C_MANIFEST_RELATIVE
    contract_path = release_root / R12C_CONTRACT_RELATIVE
    g5 = load_json(g5_path)
    manifest = load_json(manifest_path)
    contract = load_json(contract_path)
    _require_equal(g5.get("gate_g5_passed"), True, "r13_amendment_g5_not_passed")
    _require_equal(g5.get("receipt_revision"), 2, "r13_amendment_g5_revision_mismatch")
    _require_equal(g5.get("release_id"), RELEASE_ID, "r13_amendment_g5_release_mismatch")
    _require_equal(manifest.get("freeze_revision"), "R12C-EXECUTABLE-REFREEZE-002", "r13_amendment_freeze_mismatch")
    _require_equal(contract.get("call_count"), 1_180, "r13_amendment_contract_count_mismatch")
    _require_equal(contract.get("truth_available_to_worker"), False, "r13_amendment_worker_truth_exposed")
    return {
        f"repository/{G5_RECEIPT_PATH.as_posix()}": file_sha256(g5_path),
        f"r12c_release/{R12C_MANIFEST_RELATIVE.as_posix()}": file_sha256(manifest_path),
        f"r12c_release/{R12C_CONTRACT_RELATIVE.as_posix()}": file_sha256(contract_path),
    }


def _validate_checkpoint(
    run_root: Path,
    configuration: Mapping[str, Any],
    *,
    bcftools: str,
    command_runner: CommandRunner,
) -> tuple[dict[str, Any], dict[str, str]]:
    status = checkpoint_status(run_root, verify_output_hashes=True)
    expected_status = {
        "completed": 1,
        "pending": 1_179,
        "running": 0,
        "failed": 0,
        "invalid_completed_outputs": 0,
        "resume_safe": True,
    }
    for key, expected in expected_status.items():
        _require_equal(status.get(key), expected, f"r13_amendment_ledger_{key}_mismatch")
    ledger = load_checkpoint_ledger(run_root)
    completed = [row for row in ledger["calls"] if row["status"] == "COMPLETED"]
    _require_equal(completed[0].get("call_id"), RETAINED_CALL_ID, "r13_amendment_retained_call_mismatch")

    evidence = configuration.get("local_evidence", {})
    receipt_path = _path_within(run_root, evidence.get("shapeit_receipt"), code="r13_amendment_shapeit_receipt_missing")
    bcf_path = _path_within(run_root, evidence.get("shapeit_bcf"), code="r13_amendment_shapeit_bcf_missing")
    csi_path = _path_within(run_root, evidence.get("shapeit_csi"), code="r13_amendment_shapeit_csi_missing")
    log_path = _path_within(run_root, evidence.get("shapeit_log"), code="r13_amendment_shapeit_log_missing")
    receipt = load_json(receipt_path)
    retained = configuration.get("retained_checkpoint", {})
    checks = {
        "call_id": RETAINED_CALL_ID,
        "attempt": retained.get("expected_attempt"),
        "release_id": RELEASE_ID,
        "state": "COMPLETED",
        "shell_used": False,
        "truth_available_to_worker": False,
    }
    for key, expected in checks.items():
        _require_equal(receipt.get(key), expected, f"r13_amendment_receipt_{key}_mismatch")
    _require_equal(receipt.get("measurement", {}).get("returncode"), 0, "r13_amendment_shapeit_exit_nonzero")
    receipt_outputs = {item.get("relative_path"): item for item in receipt.get("outputs", [])}
    for path in (bcf_path, csi_path, log_path):
        relative = path.relative_to(run_root).as_posix()
        item = receipt_outputs.get(relative)
        if not isinstance(item, Mapping):
            _fail("r13_amendment_receipt_output_missing", relative)
        _require_equal(item.get("size_bytes"), path.stat().st_size, "r13_amendment_output_size_mismatch")
        _require_equal(item.get("sha256"), file_sha256(path), "r13_amendment_output_hash_mismatch")

    log_audit = _parse_shapeit_log(log_path)
    expected_log = {
        "variant_count": retained.get("expected_variants"),
        "study_samples": retained.get("expected_study_samples"),
        "reference_samples": retained.get("expected_reference_samples"),
        "trios": retained.get("expected_trios"),
        "mendel_errors": retained.get("expected_mendel_errors"),
        "threads": retained.get("expected_threads"),
        "seed": retained.get("expected_seed"),
        "output_samples": retained.get("expected_study_samples"),
        "output_variants": retained.get("expected_variants"),
    }
    for key, expected in expected_log.items():
        _require_equal(log_audit.get(key), expected, f"r13_amendment_shapeit_{key}_mismatch")
    dimensions = _bcf_dimensions(bcf_path, bcftools=bcftools, command_runner=command_runner)
    _require_equal(dimensions["samples"], retained.get("expected_study_samples"), "r13_amendment_bcf_samples_mismatch")
    _require_equal(dimensions["variants"], retained.get("expected_variants"), "r13_amendment_bcf_variants_mismatch")
    return (
        {
            "call_id": RETAINED_CALL_ID,
            "attempt": receipt["attempt"],
            "tool_version": log_audit["shapeit_version"],
            "study_samples": log_audit["study_samples"],
            "reference_samples": log_audit["reference_samples"],
            "trios": log_audit["trios"],
            "mendel_errors": log_audit["mendel_errors"],
            "variants": log_audit["variant_count"],
            "output_samples": dimensions["samples"],
            "output_variants": dimensions["variants"],
            "wall_seconds": receipt["measurement"]["wall_seconds"],
            "maximum_rss_native": receipt["measurement"]["maximum_rss_native"],
            "shell_used": receipt["shell_used"],
            "truth_available_to_worker": receipt["truth_available_to_worker"],
            "verdict": "INTEGRATION_VERIFIED",
        },
        {
            f"r13_run/{LEDGER_RELATIVE.as_posix()}": file_sha256(run_root / LEDGER_RELATIVE),
            f"r13_run/{receipt_path.relative_to(run_root).as_posix()}": file_sha256(receipt_path),
            f"r13_run/{bcf_path.relative_to(run_root).as_posix()}": file_sha256(bcf_path),
            f"r13_run/{csi_path.relative_to(run_root).as_posix()}": file_sha256(csi_path),
            f"r13_run/{log_path.relative_to(run_root).as_posix()}": file_sha256(log_path),
        },
    )


def _verify_checksum_index(release_root: Path) -> dict[str, str]:
    checksum_path = release_root / "governance/release_checksums.sha256"
    if not checksum_path.is_file():
        checksum_path = release_root / "release_checksums.sha256"
    if not checksum_path.is_file():
        _fail("r13_amendment_ibd_checksum_index_missing", str(release_root))
    verified: dict[str, str] = {}
    for line in checksum_path.read_text(encoding="utf-8").splitlines():
        try:
            expected, relative = line.split("  ", 1)
        except ValueError as error:
            raise R13IntegrationAmendmentError(
                "r13_amendment_ibd_checksum_line_invalid", line
            ) from error
        path = _path_within(
            release_root,
            relative,
            code="r13_amendment_ibd_checksum_target_missing",
        )
        _require_equal(
            file_sha256(path),
            expected,
            "r13_amendment_ibd_checksum_mismatch",
        )
        verified[relative] = expected
    return {
        checksum_path.relative_to(release_root).as_posix(): file_sha256(checksum_path),
    }


def _validate_ibd_evidence(
    ibd_release_root: Path,
    configuration: Mapping[str, Any],
) -> tuple[dict[str, Any], dict[str, str]]:
    local_evidence = configuration.get("local_evidence", {})
    manifest_path = ibd_release_root / "governance/frozen_candidate_manifest.json"
    manifest = load_json(manifest_path)
    _require_equal(
        manifest.get("release_id"),
        local_evidence.get("ibd_source_release_id"),
        "r13_amendment_ibd_release_mismatch",
    )
    checksum_audit = _verify_checksum_index(ibd_release_root)
    relative = configuration.get("local_evidence", {}).get("ibd_evaluation")
    path = _path_within(ibd_release_root, relative, code="r13_amendment_ibd_evidence_missing")
    evidence = load_json(path)
    tools = evidence.get("tools", {})
    result: dict[str, Any] = {}
    for tool in ("hap_ibd", "refined_ibd"):
        metrics = tools.get(tool, {})
        expected = metrics.get("designated_duplicate_expected")
        detected = metrics.get("designated_duplicate_detected")
        if not isinstance(expected, int) or expected <= 0 or not isinstance(detected, int) or detected <= 0:
            _fail("r13_amendment_ibd_positive_control_missing", tool)
        if metrics.get("segment_count", 0) <= 0 or metrics.get("pair_count", 0) <= 0:
            _fail("r13_amendment_ibd_output_empty", tool)
        result[tool] = {
            "designated_duplicate_detected": detected,
            "designated_duplicate_expected": expected,
            "pair_count": metrics["pair_count"],
            "segment_count": metrics["segment_count"],
            "verdict": "LOCAL_INTEGRATION_SUPPORTED",
        }
    evidence_hashes = {
        "r12_ibd_release/governance/frozen_candidate_manifest.json": file_sha256(manifest_path),
        f"r12_ibd_release/{relative}": file_sha256(path),
        **{
            f"r12_ibd_release/{key}": value
            for key, value in checksum_audit.items()
        },
    }
    outputs = local_evidence.get("ibd_outputs")
    if not isinstance(outputs, list) or len(outputs) != 6:
        _fail("r13_amendment_ibd_output_inventory_invalid")
    for output_relative in outputs:
        output_path = _path_within(
            ibd_release_root,
            output_relative,
            code="r13_amendment_ibd_output_missing",
        )
        if output_path.stat().st_size <= 0:
            _fail("r13_amendment_ibd_output_empty", output_relative)
        evidence_hashes[f"r12_ibd_release/{output_relative}"] = file_sha256(output_path)
    return result, evidence_hashes


def audit_amendment(
    repository_root: Path,
    release_root: Path,
    run_root: Path,
    *,
    ibd_release_root: Path | None = None,
    configuration_path: Path | None = None,
    bcftools: str = "bcftools",
    command_runner: CommandRunner = subprocess.run,
) -> dict[str, Any]:
    """Audite toutes les preuves sans modifier le registre historique."""
    repository_root = repository_root.resolve(strict=True)
    release_root = release_root.resolve(strict=True)
    run_root = run_root.resolve(strict=True)
    ibd_release_root = (ibd_release_root or release_root).resolve(strict=True)
    config_path = configuration_path or repository_root / CONFIGURATION_PATH
    config_path = config_path.resolve(strict=True)
    configuration = load_yaml(config_path)
    _validate_identity(configuration)
    governance_hashes = _validate_original_governance(repository_root, release_root)
    shapeit, checkpoint_hashes = _validate_checkpoint(
        run_root,
        configuration,
        bcftools=bcftools,
        command_runner=command_runner,
    )
    ibd, ibd_hashes = _validate_ibd_evidence(ibd_release_root, configuration)
    literature = configuration.get("literature")
    if not isinstance(literature, list) or {entry.get("tool") for entry in literature} != {"SHAPEIT5", "Hap-IBD", "Refined IBD"}:
        _fail("r13_amendment_literature_incomplete")
    decision = configuration["decision"]
    return {
        "schema_version": "1.0.0",
        "receipt_id": "R13-INTEGRATION-AMENDMENT-PASS",
        "amendment_id": AMENDMENT_ID,
        "scope": configuration["scope"],
        "release_id": RELEASE_ID,
        "run_id": RUN_ID,
        "partition_id": configuration["partition_id"],
        "status": "CLOSED_LIGHTWEIGHT_INTEGRATION_VERIFIED",
        "shapeit5_integration": shapeit,
        "explicit_ibd_integration": ibd,
        "literature": literature,
        "original_campaign": {
            "historical_call_count": 1_180,
            "retained_completed_calls": decision["retained_completed_calls"],
            "withdrawn_unexecuted_calls": decision["withdrawn_unexecuted_calls"],
            "newly_authorized_calls": decision["newly_authorized_calls"],
            "disposition": "SUPERSEDED_NOT_EXECUTED",
        },
        "evidence_sha256": {
            f"repository/{config_path.relative_to(repository_root).as_posix()}": file_sha256(config_path),
            **governance_hashes,
            **checkpoint_hashes,
            **ibd_hashes,
        },
        "amendment_implementation_sha256": {
            "r13_integration_amendment.py": file_sha256(Path(__file__)),
            "r13_checkpoint_ledger.py": file_sha256(
                Path(__file__).with_name("r13_checkpoint_ledger.py")
            ),
        },
        "truth_opened": False,
        "full_performance_qualification_evaluated": False,
        "scientific_method_modified": False,
        "qualification_flags": decision["qualification_flags"],
        "restrictions": decision["restrictions"],
    }


def _execution_overlay(audit: Mapping[str, Any]) -> dict[str, Any]:
    return {
        "schema_version": "1.0.0",
        "amendment_id": AMENDMENT_ID,
        "release_id": RELEASE_ID,
        "run_id": RUN_ID,
        "state": "EXECUTION_CLOSED_BY_AMENDMENT",
        "historical_call_count": 1_180,
        "retained_call_ids": [RETAINED_CALL_ID],
        "retained_completed_calls": 1,
        "withdrawn_unexecuted_calls": 1_179,
        "authorized_call_ids": [],
        "original_ledger_sha256": audit["evidence_sha256"][
            f"r13_run/{LEDGER_RELATIVE.as_posix()}"
        ],
        "reason": "INTRINSIC_TOOL_VALIDITY_REFERENCED_FROM_PUBLICATIONS_LOCAL_INTEGRATION_VERIFIED",
    }


def _disposition_markdown(audit: Mapping[str, Any], overlay: Mapping[str, Any]) -> str:
    shapeit = audit["shapeit5_integration"]
    return "\n".join(
        (
            "# Amendement allégé R13 — suivi des checkpoints",
            "",
            f"- Amendement : `{AMENDMENT_ID}`",
            "- Portée : vérification locale de l'intégration des outils",
            f"- Appel retenu et terminé : {overlay['retained_completed_calls']}",
            f"- Appels non exécutés retirés : {overlay['withdrawn_unexecuted_calls']}",
            "- Nouveaux appels autorisés : 0",
            "- Qualification scientifique complète : NON ÉVALUÉE",
            "- Effet fondateur : NON DÉMONTRÉ",
            "",
            "| Appel | Décision | Preuve principale |",
            "|---|---|---|",
            f"| `{RETAINED_CALL_ID}` | RETENU — INTÉGRATION VÉRIFIÉE | SHAPEIT5 {shapeit['tool_version']}, {shapeit['study_samples']} sujets, {shapeit['reference_samples']} références, {shapeit['variants']} variants, {shapeit['trios']} trios, {shapeit['mendel_errors']} erreur mendélienne |",
            "| 1 179 autres appels du registre | RETIRÉS PAR AMENDEMENT — NON EXÉCUTÉS | répétitions non requises pour la seule preuve d'intégration |",
            "",
            "Le registre JSON initial reste inchangé pour préserver l'audit. Cette vue de",
            "disposition et le verrou d'exécution font autorité après l'amendement.",
            "",
        )
    )


def apply_amendment(
    repository_root: Path,
    release_root: Path,
    run_root: Path,
    *,
    ibd_release_root: Path | None = None,
    configuration_path: Path | None = None,
    bcftools: str = "bcftools",
    command_runner: CommandRunner = subprocess.run,
) -> dict[str, Any]:
    """Audite, ferme l'exécution et publie le reçu suivi de l'amendement."""
    repository_root = repository_root.resolve(strict=True)
    run_root = run_root.resolve(strict=True)
    audit = audit_amendment(
        repository_root,
        release_root,
        run_root,
        ibd_release_root=ibd_release_root,
        configuration_path=configuration_path,
        bcftools=bcftools,
        command_runner=command_runner,
    )
    overlay = _execution_overlay(audit)
    overlay_path = run_root / AMENDMENT_RELATIVE
    _write_json(overlay_path, overlay)
    disposition_path = run_root / DISPOSITION_RELATIVE
    disposition_path.write_text(_disposition_markdown(audit, overlay), encoding="utf-8")
    audit["closed_at"] = _timestamp()
    audit["execution_overlay_sha256"] = file_sha256(overlay_path)
    audit["checkpoint_disposition_sha256"] = file_sha256(disposition_path)
    audit["amended_execution_state"] = overlay["state"]
    receipt_path = repository_root / READINESS_RECEIPT_PATH
    _write_json(receipt_path, audit)
    return audit


def _parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("action", choices=("audit", "apply"))
    parser.add_argument("--repository-root", type=Path, default=Path.cwd())
    parser.add_argument("--release-root", type=Path, required=True)
    parser.add_argument("--run-root", type=Path, required=True)
    parser.add_argument("--ibd-release-root", type=Path)
    parser.add_argument("--configuration", type=Path)
    parser.add_argument("--bcftools", default="bcftools")
    return parser


def main(argv: list[str] | None = None) -> int:
    arguments = _parser().parse_args(argv)
    function = apply_amendment if arguments.action == "apply" else audit_amendment
    result = function(
        arguments.repository_root,
        arguments.release_root,
        arguments.run_root,
        ibd_release_root=arguments.ibd_release_root,
        configuration_path=arguments.configuration,
        bcftools=arguments.bcftools,
    )
    print(json.dumps(result, ensure_ascii=False, sort_keys=True))
    return 0


if __name__ == "__main__":
    sys.exit(main())
