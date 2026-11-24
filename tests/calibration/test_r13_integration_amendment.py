from __future__ import annotations

import json
from pathlib import Path
import subprocess

import pytest
import yaml

from effet_fondateur.calibration.contracts import file_sha256, load_json
from effet_fondateur.calibration.r13_checkpoint_ledger import (
    AMENDMENT_RELATIVE,
    R13CheckpointError,
    claim_next_batch,
    complete_call,
    initialize_checkpoint_ledger,
)
from effet_fondateur.calibration.r13_integration_amendment import (
    R13IntegrationAmendmentError,
    apply_amendment,
    audit_amendment,
)


def _write_json(path: Path, value: object) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(value, indent=2) + "\n", encoding="utf-8")


def _configuration() -> dict[str, object]:
    return {
        "schema_version": "1.0.0",
        "amendment_id": "R13-LIGHTWEIGHT-INTEGRATION-001",
        "scope": "R13_TOOL_INTEGRATION_VERIFICATION",
        "release_id": "panel-master-generic-v1-r12c-scientific-validation",
        "run_id": "RUN-R13-SCIENTIFIC-001",
        "partition_id": "VALIDATION_BLIND",
        "original_campaign": {"call_count": 1180},
        "retained_checkpoint": {
            "expected_attempt": 1,
            "expected_seed": 15052001,
            "expected_threads": 1,
            "expected_study_samples": 4000,
            "expected_reference_samples": 3000,
            "expected_trios": 200,
            "expected_variants": 46524,
            "expected_mendel_errors": 0,
        },
        "local_evidence": {
            "shapeit_receipt": "outputs/attempt/receipt.json",
            "shapeit_bcf": "outputs/attempt/phased.bcf",
            "shapeit_csi": "outputs/attempt/phased.bcf.csi",
            "shapeit_log": "outputs/attempt/shapeit.log",
            "ibd_source_release_id": "panel_master_generic_v1_r12_frozen_candidate",
            "ibd_evaluation": "evaluations/explicit_ibd_r12.json",
            "ibd_outputs": [
                "outputs/explicit_ibd/hap_ibd.ibd.gz",
                "outputs/explicit_ibd/hap_ibd.hbd.gz",
                "outputs/explicit_ibd/hap_ibd.log",
                "outputs/explicit_ibd/refined_ibd.ibd.gz",
                "outputs/explicit_ibd/refined_ibd.hbd.gz",
                "outputs/explicit_ibd/refined_ibd.log",
            ],
        },
        "literature": [
            {"tool": "SHAPEIT5", "doi": "10.1038/s41588-023-01415-w"},
            {"tool": "Hap-IBD", "doi": "10.1016/j.ajhg.2020.02.010"},
            {"tool": "Refined IBD", "doi": "10.1534/genetics.113.150029"},
        ],
        "decision": {
            "retained_completed_calls": 1,
            "withdrawn_unexecuted_calls": 1179,
            "newly_authorized_calls": 0,
            "truth_opening_authorized": False,
            "qualification_flags": {
                "calibrated": False,
                "founder_effect_proven": False,
                "production_promotion": False,
            },
            "restrictions": ["INTEGRATION_EVIDENCE_ONLY"],
        },
    }


def _fixture(tmp_path: Path) -> tuple[Path, Path, Path, Path]:
    repository_root = tmp_path / "repository"
    release_root = tmp_path / "release"
    run_root = tmp_path / "run"
    repository_root.mkdir()
    release_root.mkdir()
    config_path = repository_root / "amendment.yaml"
    config_path.write_text(yaml.safe_dump(_configuration()), encoding="utf-8")
    _write_json(
        repository_root / "calibration/readiness/R12C_G5_READY.json",
        {
            "gate_g5_passed": True,
            "receipt_revision": 2,
            "release_id": "panel-master-generic-v1-r12c-scientific-validation",
        },
    )
    _write_json(
        release_root / "governance/r12c_frozen_manifest.json",
        {"freeze_revision": "R12C-EXECUTABLE-REFREEZE-002"},
    )
    _write_json(
        release_root / "governance/r13_execution_contract.json",
        {"call_count": 1180, "truth_available_to_worker": False},
    )
    _write_json(
        release_root / "governance/frozen_candidate_manifest.json",
        {"release_id": "panel_master_generic_v1_r12_frozen_candidate"},
    )
    _write_json(
        release_root / "evaluations/explicit_ibd_r12.json",
        {
            "tools": {
                "hap_ibd": {
                    "designated_duplicate_expected": 25,
                    "designated_duplicate_detected": 23,
                    "segment_count": 178,
                    "pair_count": 101,
                },
                "refined_ibd": {
                    "designated_duplicate_expected": 25,
                    "designated_duplicate_detected": 22,
                    "segment_count": 122,
                    "pair_count": 44,
                },
            }
        },
    )
    ibd_outputs = [
        "outputs/explicit_ibd/hap_ibd.ibd.gz",
        "outputs/explicit_ibd/hap_ibd.hbd.gz",
        "outputs/explicit_ibd/hap_ibd.log",
        "outputs/explicit_ibd/refined_ibd.ibd.gz",
        "outputs/explicit_ibd/refined_ibd.hbd.gz",
        "outputs/explicit_ibd/refined_ibd.log",
    ]
    for relative in ibd_outputs:
        path = release_root / relative
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_bytes(f"evidence:{relative}".encode())
    checksum_targets = sorted(
        path
        for path in release_root.rglob("*")
        if path.is_file() and path.name != "release_checksums.sha256"
    )
    checksum_path = release_root / "governance/release_checksums.sha256"
    checksum_path.write_text(
        "".join(
            f"{file_sha256(path)}  {path.relative_to(release_root).as_posix()}\n"
            for path in checksum_targets
        ),
        encoding="utf-8",
    )
    initialize_checkpoint_ledger(run_root)
    call = claim_next_batch(run_root, batch_size=1, worker_id="TEST")[0]
    attempt = run_root / "outputs/attempt"
    attempt.mkdir(parents=True)
    bcf = attempt / "phased.bcf"
    csi = attempt / "phased.bcf.csi"
    log = attempt / "shapeit.log"
    bcf.write_bytes(b"synthetic-bcf")
    csi.write_bytes(b"synthetic-index")
    log.write_text(
        """* Version       : 5.1.1 / commit = 5.1.1
* Seed    : 15052001
* Threads : 1 threads
Variants [#sites=46524 / region=1]
Samples [#target=4000 / #reference=3000]
#trios = 200
%mendel_errors = 0.00% (n=0)
VCF/BCF writing [N=4000 / L=46524]
Total running time = 21547 seconds
""",
        encoding="utf-8",
    )
    receipt = {
        "attempt": 1,
        "call_id": call["call_id"],
        "release_id": "panel-master-generic-v1-r12c-scientific-validation",
        "state": "COMPLETED",
        "shell_used": False,
        "truth_available_to_worker": False,
        "measurement": {
            "returncode": 0,
            "wall_seconds": 21548.0,
            "maximum_rss_native": 2248622080,
        },
        "outputs": [],
    }
    for path in (bcf, log, csi):
        receipt["outputs"].append(
            {
                "relative_path": path.relative_to(run_root).as_posix(),
                "sha256": file_sha256(path),
                "size_bytes": path.stat().st_size,
            }
        )
    receipt_path = attempt / "receipt.json"
    _write_json(receipt_path, receipt)
    complete_call(
        run_root,
        call_id=call["call_id"],
        output_paths=[bcf, log, csi, receipt_path],
    )
    return repository_root, release_root, run_root, config_path


def _fake_bcftools(command: list[str], **_: object) -> subprocess.CompletedProcess[str]:
    output = "\n".join(f"SAMPLE-{index}" for index in range(4000)) + "\n"
    if command[1] == "index":
        output = "46524\n"
    return subprocess.CompletedProcess(command, 0, stdout=output, stderr="")


def test_amendment_audits_evidence_and_closes_new_claims(tmp_path: Path) -> None:
    repository_root, release_root, run_root, config_path = _fixture(tmp_path)

    receipt = apply_amendment(
        repository_root,
        release_root,
        run_root,
        configuration_path=config_path,
        command_runner=_fake_bcftools,
    )

    assert receipt["status"] == "CLOSED_LIGHTWEIGHT_INTEGRATION_VERIFIED"
    assert receipt["shapeit5_integration"]["verdict"] == "INTEGRATION_VERIFIED"
    assert receipt["original_campaign"]["withdrawn_unexecuted_calls"] == 1179
    assert receipt["full_performance_qualification_evaluated"] is False
    assert not any(receipt["qualification_flags"].values())
    assert (run_root / AMENDMENT_RELATIVE).is_file()
    with pytest.raises(R13CheckpointError, match="execution_closed_by_amendment"):
        claim_next_batch(run_root, batch_size=1, worker_id="TEST")


def test_amendment_rejects_shapeit_dimension_drift(tmp_path: Path) -> None:
    repository_root, release_root, run_root, config_path = _fixture(tmp_path)

    def wrong_dimensions(command: list[str], **kwargs: object) -> subprocess.CompletedProcess[str]:
        result = _fake_bcftools(command, **kwargs)
        if command[1] == "index":
            result.stdout = "46523\n"
        return result

    with pytest.raises(R13IntegrationAmendmentError, match="bcf_variants_mismatch"):
        audit_amendment(
            repository_root,
            release_root,
            run_root,
            configuration_path=config_path,
            command_runner=wrong_dimensions,
        )


def test_tracked_receipt_binds_current_amendment_implementation() -> None:
    repository_root = Path(__file__).resolve().parents[2]
    receipt = load_json(
        repository_root / "calibration/readiness/R13_INTEGRATION_AMENDMENT.json"
    )

    assert receipt["amendment_implementation_sha256"] == {
        "r13_checkpoint_ledger.py": file_sha256(
            repository_root
            / "src/effet_fondateur/calibration/r13_checkpoint_ledger.py"
        ),
        "r13_integration_amendment.py": file_sha256(
            repository_root
            / "src/effet_fondateur/calibration/r13_integration_amendment.py"
        ),
    }
    assert receipt["evidence_sha256"][
        "repository/calibration/configurations/r13_integration_amendment_v1.yaml"
    ] == file_sha256(
        repository_root
        / "calibration/configurations/r13_integration_amendment_v1.yaml"
    )
