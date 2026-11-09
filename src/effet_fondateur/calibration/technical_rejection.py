"""Contrôle négatif D19 : seul le refus précis d'une cible déjà présente vaut PASS."""

from __future__ import annotations

from pathlib import Path
from typing import Any

from effet_fondateur.calibration.evaluate_technical import _dataset, _stage_artifacts, compare_fam, contained_file
from effet_fondateur.calibration.materialize_ibd import _json, _sha
from effet_fondateur.calibration.technical_evaluation import GenotypeMatrix, TechnicalEvaluationError, compare_genotypes, result
from effet_fondateur.contracts import validate_json_document


def evaluate_target_rejection(*, run_root: Path, manifest: dict[str, Any], receipt: dict[str, Any], expected: GenotypeMatrix, metadata: list[dict[str, str]], target: dict[str, Any]) -> dict[str, Any]:
    """Vérifie refus V2=4, stderr exact, tentative liée et base inchangée.

    Le pipeline reste FAILED/BLOCKED. Le PASS porte uniquement sur le test du
    garde ; aucun génotype cible injecté, taux de discordance ou succès de
    préparation ne lui est attribué.
    """
    entries = [item for item in manifest["stages"] if item["stage_id"] == "04"]
    evidence = receipt.get("rejection_evidence")
    if len(entries) != 1 or not evidence:
        return result("04", {}, ["expected_target_rejection_evidence_missing"])
    entry = entries[0]
    if entry["state"] != "FAILED" or entry["last_error_code"] != 4 or entry["attempt_count"] != 1 or entry["stage_name"] != "prepare_target_variant_dataset":
        return result("04", {}, ["expected_target_rejection_state_mismatch"])
    path = contained_file(run_root, evidence["path"])
    if not path.is_file() or _sha(path) != evidence["sha256"] or path.name != "stderr.log" or not path.relative_to(run_root).as_posix().startswith("stages/attempts/04_prepare_target_variant_dataset.") or not path.parent.name.endswith(".failed"):
        raise TechnicalEvaluationError("d19_rejection_evidence_invalid")
    stage_inputs = _json(contained_file(run_root, path.parent.relative_to(run_root).as_posix() + "/stage_inputs.json"))
    validate_json_document(stage_inputs, "stage_inputs.schema.json")
    if any(stage_inputs[key] != value for key, value in (("stage_id", "04"), ("stage_name", entry["stage_name"]), ("run_id", manifest["run_id"]), ("signature", entry["signature"]), ("attempt_number", 1))):
        raise TechnicalEvaluationError("d19_rejection_attempt_binding_invalid")
    if path.read_text().strip() != "target_variant_already_present" or (run_root / "stages/04_prepare_target_variant_dataset").exists():
        return result("04", {}, ["expected_target_rejection_reason_mismatch"])
    base = _stage_artifacts(run_root, manifest, "03")
    if base is None:
        return result("04", {}, ["expected_target_rejection_base_missing"])
    comparison = compare_genotypes(expected, _dataset(base, "target_chromosome_base"), stage="04")
    family = compare_fam(base["target_chromosome_base_fam"], metadata, stage="04")
    locus_count = sum(variant == (target["chromosome"], target["position_bp"], target["project_variant_id"]) for variant in expected.variants)
    reasons = comparison["reasons"] + family["reasons"]
    if locus_count != 1:
        reasons.append("expected_target_rejection_locus_not_present_once")
    return result("04", {**comparison["metrics"], "expected_duplicate_target_rejection_count": 1,
                         "target_insertion_performed_count": 0, "target_genotype_discordance_rate": None}, sorted(set(reasons)))
