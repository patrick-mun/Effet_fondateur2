"""Validation de la révision prospective R07 pour la qualification scientifique.

Le module ferme le dimensionnement et le calendrier de calcul sans exécuter de
simulation ni ouvrir de vérité. Il démontre que le blocage SHAPEIT5 connu est
traité par un plan exact, puis maintient la porte scientifique fermée jusqu'aux
mesures R11b encore requises pour l'IBD, le stockage et l'échelle cible.
"""

from __future__ import annotations

import argparse
import json
from math import ceil
from pathlib import Path
import sys
from typing import Any

from effet_fondateur.calibration.contracts import file_sha256, load_json, load_yaml


CONFIGURATION_PATH = Path("calibration/configurations/scientific_dimensioning_v2.yaml")
PROJECTION_PATH = Path(
    "calibration_workspace/resource_probes/"
    "r11_shapeit_sample_scaling_v1/projection.json"
)
RECEIPT_PATH = Path("calibration/readiness/R07_SCIENTIFIC_DIMENSIONING_V2.json")


class R07RevisionError(ValueError):
    """Erreur de contrat stable de la révision R07."""

    def __init__(self, code: str, detail: str | None = None) -> None:
        self.code = code
        self.detail = detail
        super().__init__(code if detail is None else f"{code}: {detail}")


def _expect(condition: bool, code: str, detail: str | None = None) -> None:
    if not condition:
        raise R07RevisionError(code, detail)


def _qualification_flags() -> dict[str, bool]:
    return {
        "calibrated": False,
        "founder_effect_proven": False,
        "production_promotion": False,
    }


def assess_r07_revision(repository_root: Path) -> dict[str, Any]:
    """Valide R07 v2 contre la projection R11 et retourne son reçu déterministe."""
    repository_root = repository_root.resolve(strict=True)
    configuration_path = repository_root / CONFIGURATION_PATH
    projection_path = repository_root / PROJECTION_PATH
    configuration = load_yaml(configuration_path)
    projection = load_json(projection_path)

    _expect(configuration.get("configuration_version") == "2.0.0", "r07_version_mismatch")
    _expect(
        configuration.get("status") == "REVISED_PROSPECTIVE_DIMENSIONING_COMPLETE",
        "r07_status_mismatch",
    )
    _expect(
        projection.get("decision", {}).get("status") == "RESOURCE_BLOCKED",
        "r11_source_decision_mismatch",
    )
    minima = configuration.get("scientific_minima", {})
    expected_minima = {
        "study_identities_per_partition": 4_000,
        "phasing_reference_per_partition": 3_000,
        "ancestry_reference_per_partition": 1_800,
        "null_reference_per_partition": 3_000,
        "primary_independent_units": 500,
        "primary_negative_units": 3_000,
        "founder_events_per_age": 500,
        "external_tool_seed_minimum": 50,
        "autosomes": 22,
    }
    for field, expected in expected_minima.items():
        _expect(minima.get(field) == expected, "r07_scientific_minimum_mismatch", field)

    phasing = configuration.get("phasing_schedule", {})
    sentinels = phasing.get("sentinel_autosomes", {})
    all_autosome_calls = int(phasing.get("all_autosome_seed_count", 0)) * 22
    sentinel_extra_calls = (
        int(sentinels.get("total_seed_count", 0))
        - int(phasing.get("all_autosome_seed_count", 0))
    ) * int(sentinels.get("count", 0))
    expected_calls = phasing.get("expected_calls", {})
    _expect(all_autosome_calls + sentinel_extra_calls == 340, "r07_dmax_call_count_mismatch")
    _expect(expected_calls.get("D_MAX") == 340, "r07_dmax_schedule_mismatch")
    _expect(expected_calls.get("D_MID") == 150, "r07_dmid_schedule_mismatch")
    _expect(expected_calls.get("D_LOW") == 150, "r07_dlow_schedule_mismatch")
    _expect(expected_calls.get("total") == 640, "r07_total_schedule_mismatch")

    ibd = configuration.get("ibd_worlds", {})
    worlds = ibd.get("expected_independent_worlds", {})
    _expect(worlds.get("total") == 90, "r07_ibd_world_count_mismatch")
    expected_ibd_calls = (
        int(worlds.get("total", 0))
        * len(ibd.get("tools", ()))
        * len(ibd.get("densities", ()))
    )
    _expect(expected_ibd_calls == 540, "r07_ibd_call_count_mismatch")
    _expect(
        ibd.get("expected_tool_density_calls") == expected_ibd_calls,
        "r07_ibd_schedule_mismatch",
    )

    resources = configuration.get("resources", {})
    allocations = resources.get("core_hour_allocations", {})
    _expect(
        sum(int(value) for value in allocations.values())
        == resources.get("aggregate_core_hours_max")
        == 120_000,
        "r07_core_hour_allocation_mismatch",
    )
    safe_rss_gib = projection["fits"]["maximum_rss_bytes"]["safe_projection_gibibytes"]
    workers = int(resources.get("r2_workers_simultaneous_max", 0))
    aggregate_rss = float(resources.get("aggregate_rss_gib_max", 0))
    _expect(workers * safe_rss_gib <= aggregate_rss, "r07_rss_concurrency_exceeded")
    _expect(
        resources.get("r2_autosomal_call_timeout_hours", 0)
        > projection["fits"]["wall_seconds"]["safe_projection_hours"],
        "r07_shapeit_timeout_insufficient",
    )

    dmax_core_hours = projection["fits"]["cpu_seconds"][
        "safe_projection_core_hours_per_call"
    ]
    dmid_ratio = projection["density_control"]["cpu_seconds_ratio"]
    density_exponent = projection["density_control"]["cpu_density_exponent"]
    dlow_ratio = (1 / 6) ** density_exponent
    projected_phasing_core_hours = (
        expected_calls["D_MAX"] * dmax_core_hours
        + expected_calls["D_MID"] * dmax_core_hours * dmid_ratio
        + expected_calls["D_LOW"] * dmax_core_hours * dlow_ratio
    )
    reserved_phasing_core_hours = projected_phasing_core_hours * (
        1 + resources["phasing_projection"]["safety_reserve_fraction"]
    )
    _expect(
        reserved_phasing_core_hours <= allocations.get("phasing", 0),
        "r07_phasing_allocation_exceeded",
    )
    _expect(
        resources["phasing_projection"]["projected_core_hours_max"]
        >= ceil(reserved_phasing_core_hours),
        "r07_phasing_projection_cap_insufficient",
    )

    required_evidence = configuration.get("required_resource_evidence_before_g5_scientific")
    _expect(isinstance(required_evidence, list) and len(required_evidence) == 7, "r07_evidence_set_mismatch")
    gate = configuration.get("gate", {})
    _expect(gate.get("r07_revision_complete") is True, "r07_revision_gate_mismatch")
    _expect(gate.get("resource_confirmation_complete") is False, "r07_resource_gate_mismatch")
    _expect(gate.get("scientific_g5_authorized") is False, "r07_scientific_gate_mismatch")
    _expect(configuration.get("qualification_flags") == _qualification_flags(), "r07_flags_mismatch")

    return {
        "schema_version": "1.0.0",
        "receipt_id": "R07-SCIENTIFIC-DIMENSIONING-V2",
        "revision_id": configuration["revision_id"],
        "status": "REVISED_READY_FOR_R11B_RESOURCE_CONFIRMATION",
        "scientific_minima_preserved": True,
        "shapeit_schedule": {
            "total_calls": expected_calls["total"],
            "projected_core_hours": round(projected_phasing_core_hours, 6),
            "projected_core_hours_with_25_percent_reserve": round(
                reserved_phasing_core_hours, 6
            ),
            "allocation_core_hours": allocations["phasing"],
            "per_call_timeout_hours": resources["r2_autosomal_call_timeout_hours"],
            "simultaneous_workers_max": workers,
            "safe_rss_gib_per_worker": round(safe_rss_gib, 6),
        },
        "ibd_schedule": {
            "independent_worlds": worlds["total"],
            "tool_density_calls": expected_ibd_calls,
            "allocation_core_hours": allocations["explicit_ibd"],
            "resource_projection_complete": False,
        },
        "required_resource_evidence_before_g5_scientific": required_evidence,
        "r07_revision_complete": True,
        "scientific_g5_authorized": False,
        "next_phase": "R11B_RESOURCE_CONFIRMATION",
        "evidence_sha256": {
            CONFIGURATION_PATH.as_posix(): file_sha256(configuration_path),
            PROJECTION_PATH.as_posix(): file_sha256(projection_path),
        },
        "qualification_flags": _qualification_flags(),
        "restrictions": [
            "R07_REVISION_DOES_NOT_ITSELF_QUALIFY_SCIENTIFIC_PERFORMANCE",
            "R11B_RESOURCE_EVIDENCE_REQUIRED",
            "NEW_SCIENTIFIC_RELEASE_AND_BLIND_VALIDATION_REQUIRED",
            "NO_REAL_DATA_AUTHORIZED_BY_THIS_RECEIPT",
        ],
    }


def validate_tracked_receipt(repository_root: Path) -> dict[str, Any]:
    """Compare le reçu suivi au résultat déterministe du validateur."""
    expected = assess_r07_revision(repository_root)
    actual = load_json(repository_root.resolve(strict=True) / RECEIPT_PATH)
    _expect(actual == expected, "r07_tracked_receipt_mismatch")
    return actual


def _parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--repository-root", type=Path, default=Path.cwd())
    parser.add_argument("--validate-receipt", action="store_true")
    return parser


def main(argv: list[str] | None = None) -> int:
    arguments = _parser().parse_args(argv)
    document = (
        validate_tracked_receipt(arguments.repository_root)
        if arguments.validate_receipt
        else assess_r07_revision(arguments.repository_root)
    )
    print(json.dumps(document, ensure_ascii=False, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    sys.exit(main())
