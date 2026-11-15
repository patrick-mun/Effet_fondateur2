"""Contrat de planification du runner scientifique R13.

Le module valide le plan complet avant tout lancement. Son état reste
Le plan conserve le volume scientifique intégral et autorise une exécution
reprise par lots de huit workers. Le parallélisme ne modifie ni les scénarios,
ni les graines, ni les 1 180 appels externes.
"""

from __future__ import annotations

import argparse
from concurrent.futures import ThreadPoolExecutor, as_completed
import json
import os
from pathlib import Path, PurePosixPath
import subprocess
import sys
from typing import Any, Mapping

from effet_fondateur.calibration.contracts import STAGE_IDS, canonical_sha256, load_json
from effet_fondateur.calibration.release_guard import scan_blinded_document


RELEASE_ID = "panel-master-generic-v1-r12c-scientific-validation"
PARTITION_ID = "VALIDATION_BLIND"
SCOPE = "R13_SCIENTIFIC_BLIND_VALIDATION"
RUNNER_IMPLEMENTATION_STATE = "FULL_EXECUTION_READY"
REQUIRED_WORKERS = 8
EXPECTED_SCHEDULE = {
    "shapeit": {"D_MAX": 340, "D_MID": 150, "D_LOW": 150, "total": 640},
    "ibd": {"independent_worlds": 90, "tool_density_calls": 540},
}
SENTINEL_AUTOSOMES = (
    "SENTINEL_SHORTEST",
    "SENTINEL_MEDIAN_LOWER",
    "SENTINEL_LONGEST",
)
IBD_TOOLS = ("HAP_IBD", "REFINED_IBD")
DENSITIES = ("D_MAX", "D_MID", "D_LOW")
IBD_LENGTHS_CM = (1.0, 1.5, 2.0, 3.0, 5.0)


class R13ScientificRunnerError(RuntimeError):
    """Plan R13 invalide ou environnement d'exécution insuffisant."""

    def __init__(self, code: str, detail: str | None = None) -> None:
        self.code = code
        self.detail = detail
        super().__init__(code if detail is None else f"{code}: {detail}")


def _fail(code: str, detail: str | None = None) -> None:
    raise R13ScientificRunnerError(code, detail)


def _safe_relative(value: Any, detail: str) -> str:
    if not isinstance(value, str) or not value or "\\" in value:
        _fail("r13_plan_path_invalid", detail)
    path = PurePosixPath(value)
    if path.is_absolute() or any(part in {"", ".", ".."} for part in path.parts):
        _fail("r13_plan_path_invalid", detail)
    folded = {part.casefold() for part in path.parts}
    if folded.intersection({"private_truth", "governance"}):
        _fail("r13_private_path_forbidden", detail)
    if tuple(part.casefold() for part in path.parts[:2]) == ("data", "input"):
        _fail("r13_real_input_path_forbidden", detail)
    return path.as_posix()


def build_external_call_schedule() -> dict[str, list[dict[str, Any]]]:
    """Construit l'inventaire R07 exact sans lancer d'outil externe."""
    shapeit_calls: list[dict[str, Any]] = []
    for chromosome in range(1, 23):
        for seed_index in range(1, 11):
            shapeit_calls.append(
                {
                    "call_id": f"SHAPEIT-DMAX-CHR{chromosome:02d}-SEED{seed_index:02d}",
                    "density": "D_MAX",
                    "autosome": chromosome,
                    "seed_index": seed_index,
                    "schedule_class": "ALL_AUTOSOMES_BASE",
                }
            )
    for sentinel in SENTINEL_AUTOSOMES:
        for seed_index in range(11, 51):
            shapeit_calls.append(
                {
                    "call_id": f"SHAPEIT-DMAX-{sentinel}-SEED{seed_index:02d}",
                    "density": "D_MAX",
                    "autosome": sentinel,
                    "seed_index": seed_index,
                    "schedule_class": "SENTINEL_DENSITY_EXTENSION",
                }
            )
        for density in ("D_MID", "D_LOW"):
            for seed_index in range(1, 51):
                shapeit_calls.append(
                    {
                        "call_id": f"SHAPEIT-{density}-{sentinel}-SEED{seed_index:02d}",
                        "density": density,
                        "autosome": sentinel,
                        "seed_index": seed_index,
                        "schedule_class": "SENTINEL_DENSITY_SENSITIVITY",
                    }
                )

    worlds: list[dict[str, Any]] = []
    for length_cm in IBD_LENGTHS_CM:
        for world_index in range(1, 6):
            worlds.append(
                {
                    "world_id": f"IBD-POS-{str(length_cm).replace('.', '_')}-{world_index:02d}",
                    "world_class": "PRIMARY_POSITIVE",
                    "length_cm": length_cm,
                    "pair_count": 100,
                }
            )
    for world_index in range(1, 31):
        worlds.append(
            {
                "world_id": f"IBD-NEG-{world_index:02d}",
                "world_class": "PRIMARY_NEGATIVE",
                "pair_count": 100,
            }
        )
    locations = ("TARGET", "OFF_TARGET_SAME_AUTOSOME", "OTHER_AUTOSOME")
    for length_cm in IBD_LENGTHS_CM:
        for location in locations:
            for world_index in range(1, 3):
                worlds.append(
                    {
                        "world_id": (
                            f"IBD-SECONDARY-{str(length_cm).replace('.', '_')}-"
                            f"{location}-{world_index:02d}"
                        ),
                        "world_class": "SECONDARY_LOCATION",
                        "length_cm": length_cm,
                        "location": location,
                        "pair_count": 100,
                    }
                )
    for length_cm in IBD_LENGTHS_CM:
        worlds.append(
            {
                "world_id": f"IBD-COMPOUND-{str(length_cm).replace('.', '_')}",
                "world_class": "COMPOUND",
                "length_cm": length_cm,
                "pair_count": 100,
            }
        )
    ibd_calls = [
        {
            "call_id": f"{world['world_id']}-{tool}-{density}",
            "world_id": world["world_id"],
            "world_class": world["world_class"],
            "tool": tool,
            "density": density,
            "pair_count": world["pair_count"],
            **(
                {"length_cm": world["length_cm"]}
                if "length_cm" in world
                else {}
            ),
            **({"location": world["location"]} if "location" in world else {}),
        }
        for world in worlds
        for tool in IBD_TOOLS
        for density in DENSITIES
    ]
    return {"shapeit": shapeit_calls, "ibd": ibd_calls}


def _validate_external_call_schedule(external_calls: Any) -> None:
    if not isinstance(external_calls, Mapping) or set(external_calls) != {
        "shapeit",
        "ibd",
    }:
        _fail("r13_external_call_schedule_missing")
    expected = build_external_call_schedule()
    for family in ("shapeit", "ibd"):
        calls = external_calls.get(family)
        if not isinstance(calls, list) or calls != expected[family]:
            _fail("r13_external_call_schedule_mismatch", family)
        identifiers = [call.get("call_id") for call in calls if isinstance(call, Mapping)]
        if len(identifiers) != len(set(identifiers)):
            _fail("r13_external_call_id_duplicate", family)


def validate_execution_plan(document: Mapping[str, Any]) -> dict[str, Any]:
    """Valide un plan fermé sans exécuter de commande ni ouvrir la vérité."""
    scan_blinded_document(document)
    if document.get("schema_version") != "1.0.0":
        _fail("r13_plan_schema_mismatch")
    if document.get("release_id") != RELEASE_ID:
        _fail("r13_plan_release_mismatch")
    if document.get("partition_id") != PARTITION_ID:
        _fail("r13_plan_partition_mismatch")
    if document.get("scope") != SCOPE:
        _fail("r13_plan_scope_mismatch")
    if tuple(document.get("stage_ids", ())) != STAGE_IDS:
        _fail("r13_plan_stage_set_mismatch")
    if document.get("required_workers") != REQUIRED_WORKERS:
        _fail("r13_plan_worker_count_mismatch")
    if document.get("schedule") != EXPECTED_SCHEDULE:
        _fail("r13_plan_schedule_mismatch")
    _validate_external_call_schedule(document.get("external_calls"))
    if document.get("truth_available_to_runner") is not False:
        _fail("r13_plan_truth_access_forbidden")

    work_items = document.get("work_items")
    if not isinstance(work_items, list) or not work_items:
        _fail("r13_plan_work_items_missing")
    identifiers: set[str] = set()
    stages: set[str] = set()
    for item in work_items:
        if not isinstance(item, Mapping):
            _fail("r13_plan_work_item_invalid")
        identifier = item.get("work_item_id")
        stage_id = item.get("stage_id")
        if not isinstance(identifier, str) or not identifier or identifier in identifiers:
            _fail("r13_plan_work_item_id_invalid", str(identifier))
        if stage_id not in STAGE_IDS:
            _fail("r13_plan_work_item_stage_invalid", str(stage_id))
        identifiers.add(identifier)
        stages.add(str(stage_id))
        command = item.get("command")
        if (
            not isinstance(command, list)
            or not command
            or any(not isinstance(token, str) or not token for token in command)
        ):
            _fail("r13_plan_command_invalid", identifier)
        for output in item.get("expected_outputs", []):
            _safe_relative(output, identifier)
        dependencies = item.get("dependencies", [])
        if not isinstance(dependencies, list) or any(
            not isinstance(value, str) for value in dependencies
        ):
            _fail("r13_plan_dependencies_invalid", identifier)
    if stages != set(STAGE_IDS):
        _fail("r13_plan_work_item_stage_coverage_mismatch")
    for item in work_items:
        unknown = set(item.get("dependencies", ())) - identifiers
        if unknown:
            _fail("r13_plan_dependency_unknown", sorted(unknown)[0])

    return {
        "release_id": RELEASE_ID,
        "partition_id": PARTITION_ID,
        "stage_count": len(stages),
        "work_item_count": len(work_items),
        "shapeit_call_count": len(document["external_calls"]["shapeit"]),
        "ibd_call_count": len(document["external_calls"]["ibd"]),
        "required_workers": REQUIRED_WORKERS,
        "plan_sha256": canonical_sha256(document),
        "runner_implementation_state": RUNNER_IMPLEMENTATION_STATE,
        "execution_authorized": False,
    }


def validate_plan_file(path: Path) -> dict[str, Any]:
    """Charge et valide un plan JSON fermé."""
    return validate_execution_plan(load_json(path.resolve(strict=True)))


def assert_full_campaign_capacity(logical_cpus: int | None = None) -> None:
    """Refuse la campagne si les huit workers gelés ne sont pas disponibles."""
    available = logical_cpus if logical_cpus is not None else (os.cpu_count() or 1)
    if available < REQUIRED_WORKERS:
        _fail(
            "r13_execution_host_underprovisioned",
            f"available={available}:required={REQUIRED_WORKERS}",
        )
    if RUNNER_IMPLEMENTATION_STATE != "FULL_EXECUTION_READY":
        _fail("r13_runner_not_implementation_complete")


def _write_json(path: Path, document: Any) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(
        json.dumps(document, ensure_ascii=False, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )


def execute_frozen_plan(
    repository_root: Path,
    plan: Mapping[str, Any],
    run_root: Path,
    g5_receipt_path: Path,
) -> dict[str, Any]:
    """Exécute un plan G5 par lots de huit avec checkpoints reprenables.

    Les commandes sont lancées sans shell et reçoivent uniquement le chemin du
    workspace via ``EFFET_FONDATEUR_R13_RUN_ROOT``. Une reprise ne relance pas
    un work item dont les sorties attendues existent et dont le checkpoint est
    ``COMPLETED``.
    """
    repository_root = repository_root.resolve(strict=True)
    run_root = run_root.resolve(strict=False)
    audit = validate_execution_plan(plan)
    assert_full_campaign_capacity()
    receipt = load_json(g5_receipt_path.resolve(strict=True))
    if (
        receipt.get("receipt_id") != "R12C-SCIENTIFIC-G5-READY"
        or receipt.get("release_id") != RELEASE_ID
        or receipt.get("gate_g5_passed") is not True
        or receipt.get("scientific_validation_authorized") is not True
    ):
        _fail("r13_g5_receipt_invalid")
    run_root.mkdir(parents=True, exist_ok=True)
    checkpoint_path = run_root / "governance/execution_checkpoint.json"
    checkpoint: dict[str, Any] = (
        load_json(checkpoint_path)
        if checkpoint_path.is_file()
        else {
            "schema_version": "1.0.0",
            "release_id": RELEASE_ID,
            "plan_sha256": audit["plan_sha256"],
            "state": "RUNNING",
            "work_items": {},
        }
    )
    if checkpoint.get("plan_sha256") != audit["plan_sha256"]:
        _fail("r13_resume_plan_mismatch")
    items = {item["work_item_id"]: item for item in plan["work_items"]}

    def outputs_exist(item: Mapping[str, Any]) -> bool:
        return all((run_root / output).is_file() for output in item.get("expected_outputs", []))

    completed = {
        identifier
        for identifier, state in checkpoint["work_items"].items()
        if state.get("state") == "COMPLETED" and outputs_exist(items[identifier])
    }
    environment = dict(os.environ)
    environment["EFFET_FONDATEUR_R13_RUN_ROOT"] = str(run_root)

    def run_item(identifier: str) -> tuple[str, int, str, str]:
        item = items[identifier]
        completed_process = subprocess.run(
            item["command"],
            cwd=repository_root,
            env=environment,
            capture_output=True,
            text=True,
            check=False,
            timeout=int(plan.get("per_work_item_timeout_seconds", 604_800)),
        )
        return (
            identifier,
            completed_process.returncode,
            completed_process.stdout,
            completed_process.stderr,
        )

    while len(completed) < len(items):
        ready = sorted(
            identifier
            for identifier, item in items.items()
            if identifier not in completed
            and set(item.get("dependencies", ())).issubset(completed)
        )
        if not ready:
            _fail("r13_execution_dependency_deadlock")
        batch = ready[:REQUIRED_WORKERS]
        with ThreadPoolExecutor(max_workers=REQUIRED_WORKERS) as executor:
            futures = {executor.submit(run_item, identifier): identifier for identifier in batch}
            for future in as_completed(futures):
                identifier, returncode, stdout, stderr = future.result()
                log_root = run_root / "logs" / identifier
                log_root.mkdir(parents=True, exist_ok=True)
                (log_root / "stdout.log").write_text(stdout, encoding="utf-8")
                (log_root / "stderr.log").write_text(stderr, encoding="utf-8")
                item = items[identifier]
                state = "COMPLETED" if returncode == 0 and outputs_exist(item) else "FAILED"
                checkpoint["work_items"][identifier] = {
                    "state": state,
                    "returncode": returncode,
                    "expected_outputs_present": outputs_exist(item),
                }
                _write_json(checkpoint_path, checkpoint)
                if state != "COMPLETED":
                    checkpoint["state"] = "FAILED"
                    _write_json(checkpoint_path, checkpoint)
                    _fail("r13_work_item_failed", identifier)
                completed.add(identifier)
    checkpoint["state"] = "CLOSED_EXECUTION_COMPLETE"
    checkpoint["completed_work_item_count"] = len(completed)
    _write_json(checkpoint_path, checkpoint)
    return checkpoint


def _parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--plan", type=Path, required=True)
    return parser


def main(argv: list[str] | None = None) -> int:
    arguments = _parser().parse_args(argv)
    audit = validate_plan_file(arguments.plan)
    print(json.dumps(audit, ensure_ascii=False, sort_keys=True))
    return 0


if __name__ == "__main__":
    sys.exit(main())
