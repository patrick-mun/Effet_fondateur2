"""Porte statique R09 ouvrant uniquement la construction R10 du micro-panel."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Any

from effet_fondateur.calibration.contracts import (
    CalibrationContractError,
    file_sha256,
    load_json,
    load_yaml,
    validate_contract_bundle,
    validate_document,
    validate_schema_set,
)
from effet_fondateur.calibration.evaluation import validate_criteria_document
from effet_fondateur.calibration.truth_router import validate_truth_policy


CONFIGURATION_FILES = (
    "contract_bundle_v1.yaml",
    "stage_criteria_v1.yaml",
    "micro_panel_v1.yaml",
    "isolation_policy_v1.yaml",
    "tool_harness_v1.yaml",
)


class ReadinessError(ValueError):
    """Échec bloquant de la porte R09."""

    def __init__(self, code: str, detail: str | None = None) -> None:
        self.code = code
        self.detail = detail
        super().__init__(code if detail is None else f"{code}: {detail}")


@dataclass(frozen=True)
class ReadinessAssessment:
    """Résultat déterministe de la validation statique R09."""

    ready_to_run: bool
    scope: str
    terms: tuple[str, ...]
    evidence_sha256: dict[str, str]
    checks: tuple[str, ...]


def _expect(condition: bool, code: str, detail: str | None = None) -> None:
    if not condition:
        raise ReadinessError(code, detail)


def _validate_micro_configuration(document: dict[str, Any]) -> None:
    _expect(document.get("release_design") == "panel_master_generic_v1", "micro_release_mismatch")
    _expect(document.get("scale") == "RESOURCE_ONLY", "micro_scale_mismatch")
    _expect(document.get("generation_phase") == "R10", "micro_generation_phase_mismatch")
    _expect(document.get("generated") is False, "micro_must_not_be_generated_in_r09")
    counts = document.get("counts", {})
    _expect(counts.get("study_identities") == 60, "micro_study_count_mismatch")
    _expect(counts.get("phasing_reference") == 90, "micro_phase_reference_mismatch")
    _expect(counts.get("ancestry_reference") == 90, "micro_ancestry_reference_mismatch")
    _expect(counts.get("null_reference") == 60, "micro_null_reference_mismatch")
    _expect(counts.get("trios") == 4, "micro_trio_count_mismatch")
    _expect(counts.get("sibships_of_four") == 4, "micro_sibship_count_mismatch")
    _expect(counts.get("multigenerational_components_of_five") == 4, "micro_multigeneration_count_mismatch")
    _expect(counts.get("unrelated") == 12, "micro_unrelated_count_mismatch")
    resources = document.get("resources", {})
    _expect(resources.get("threads_max") == 8, "micro_threads_mismatch")
    _expect(resources.get("rss_gib_max") == 16, "micro_rss_mismatch")
    _expect(resources.get("storage_gib_max") == 100, "micro_storage_mismatch")
    _expect(resources.get("wall_hours_max") == 12, "micro_wall_time_mismatch")
    _expect(resources.get("core_hours_max") == 100, "micro_core_hours_mismatch")
    _expect(document.get("scientific_performance_decision_allowed") is False, "micro_performance_forbidden")


def _validate_isolation_policy(document: dict[str, Any]) -> None:
    _expect(document.get("pipeline_forbidden_roots") == ["private_truth", "governance"], "private_roots_mismatch")
    _expect(document.get("physical_path_verification_required_from") == "R10", "path_verification_phase_mismatch")
    _expect(document.get("pipeline_process_identity_required") is True, "pipeline_identity_guard_missing")
    _expect(document.get("symlinks_and_hardlinks_forbidden") is True, "link_guard_missing")
    validate_truth_policy(document)


def _validate_tool_harness(document: dict[str, Any]) -> None:
    _expect(document.get("r09_external_execution_allowed") is False, "external_execution_enabled")
    tools = document.get("external_tools")
    expected = {"plink", "king", "bcftools", "shapeit5", "hap_ibd", "refined_ibd", "java"}
    _expect(isinstance(tools, dict) and set(tools) == expected, "tool_harness_set_mismatch")
    for tool_id, contract in tools.items():
        _expect(isinstance(contract, dict), "tool_contract_invalid", tool_id)
        _expect(contract.get("r09_mode") == "SIMULATED_CONTRACT_ONLY", "tool_r09_mode_invalid", tool_id)
        _expect(bool(contract.get("failure_cases")), "tool_failure_cases_missing", tool_id)


def assess_r09_repository(repository_root: Path) -> ReadinessAssessment:
    """Valide les contrats R09 sans générer de panel ni appeler d'outil externe."""
    calibration_root = repository_root / "calibration"
    schema_root = calibration_root / "schemas"
    configuration_root = calibration_root / "configurations"
    plan_path = calibration_root / "PLAN_DEVELOPPEMENT_PANEL_SIMULATION.md"
    schema_hashes = validate_schema_set(schema_root)
    missing_configs = [
        name for name in CONFIGURATION_FILES if not (configuration_root / name).is_file()
    ]
    if missing_configs:
        raise ReadinessError("configuration_set_incomplete", ",".join(missing_configs))

    bundle = load_yaml(configuration_root / "contract_bundle_v1.yaml")
    criteria = load_yaml(configuration_root / "stage_criteria_v1.yaml")
    micro = load_yaml(configuration_root / "micro_panel_v1.yaml")
    isolation = load_yaml(configuration_root / "isolation_policy_v1.yaml")
    tools = load_yaml(configuration_root / "tool_harness_v1.yaml")
    try:
        validate_contract_bundle(bundle, plan_path)
        validate_criteria_document(criteria)
    except CalibrationContractError as error:
        raise ReadinessError(error.code, error.detail) from error
    _validate_micro_configuration(micro)
    _validate_isolation_policy(isolation)
    _validate_tool_harness(tools)

    conventions = calibration_root / "CONVENTIONS_DEVELOPPEMENT.md"
    _expect(conventions.is_file(), "report_conventions_missing")
    evidence = {
        f"schemas/{name}": digest for name, digest in schema_hashes.items()
    }
    for name in CONFIGURATION_FILES:
        evidence[f"configurations/{name}"] = file_sha256(configuration_root / name)
    evidence["CONVENTIONS_DEVELOPPEMENT.md"] = file_sha256(conventions)
    return ReadinessAssessment(
        ready_to_run=True,
        scope="R10_MICRO_PANEL_CONSTRUCTION",
        terms=(
            "contracts",
            "data_specification",
            "sample_sizes",
            "independence",
            "tool_contracts",
            "blinding",
            "resources",
            "expected_outputs",
        ),
        evidence_sha256=evidence,
        checks=(
            "NINE_SCHEMAS_VALID",
            "TWENTY_FIVE_STAGES_COVERED",
            "FOURTEEN_TRUTH_FAMILIES_ROUTED",
            "ELEVEN_BLINDED_VIEWS_CLOSED",
            "R08_CRITERIA_HASH_MATCHES",
            "R07_MICRO_DIMENSIONS_MATCH",
            "NO_EXTERNAL_TOOL_EXECUTION_IN_R09",
            "NO_PANEL_GENERATED_IN_R09",
            "REAL_INPUTS_FORBIDDEN",
            "TRUTH_ACCESS_AFTER_CLOSURE_ONLY",
            "DYNAMIC_OUTPUTS_REQUIRE_CLOSED_MANIFEST",
            "HTML_PDF_RESERVED_FOR_R14",
        ),
    )


def validate_readiness_receipt(repository_root: Path, receipt_path: Path) -> ReadinessAssessment:
    """Compare un reçu versionné au calcul statique courant."""
    assessment = assess_r09_repository(repository_root)
    receipt = load_json(receipt_path)
    validate_document(
        receipt,
        repository_root / "calibration" / "schemas" / "readiness_receipt.schema.json",
    )
    _expect(receipt.get("ready_to_run") is True, "receipt_not_ready")
    _expect(receipt.get("scope") == assessment.scope, "receipt_scope_mismatch")
    _expect(tuple(receipt.get("terms", ())) == assessment.terms, "receipt_terms_mismatch")
    _expect(receipt.get("checks") == list(assessment.checks), "receipt_checks_mismatch")
    _expect(receipt.get("evidence_sha256") == assessment.evidence_sha256, "receipt_evidence_mismatch")
    flags = receipt.get("qualification_flags", {})
    _expect(flags.get("calibrated") is False, "calibrated_flag_must_be_false")
    _expect(flags.get("founder_effect_proven") is False, "founder_flag_must_be_false")
    _expect(flags.get("production_promotion") is False, "promotion_flag_must_be_false")
    return assessment
