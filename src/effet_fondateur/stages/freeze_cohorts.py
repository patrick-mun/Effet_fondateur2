"""Étape 09 : résolution auditée et gel des cohortes analytiques."""

from __future__ import annotations

import argparse
import csv
import json
import math
import sys
from collections import defaultdict
from datetime import datetime
from pathlib import Path, PurePosixPath
from time import monotonic
from typing import Any, Sequence

import yaml

from effet_fondateur.audit import atomic_write_json, read_json, sha256_file
from effet_fondateur.contracts import (
    DocumentValidationError,
    TableValidationError,
    build_file_artifact,
    load_pipeline_config,
    validate_cohorts_frozen,
    validate_json_document,
    validate_tsv_table,
)
from effet_fondateur.orchestrator.state import utc_now


COHORT_IDS = (
    "controls_unrelated",
    "target_carriers_all",
    "target_carriers_independent",
    "family_noncarriers",
    "affected_noncarriers",
    "target_chromosome_all_qc",
    "genomewide_structure_reference",
)
INDEPENDENT_COHORT_IDS = {
    "controls_unrelated",
    "target_carriers_independent",
    "genomewide_structure_reference",
}
DEFAULT_REQUIRED_COHORT_IDS = (
    "controls_unrelated",
    "target_carriers_all",
    "target_carriers_independent",
    "target_chromosome_all_qc",
    "genomewide_structure_reference",
)
KINSHIP_DECISION_ID = "kinship_exclusion_approval"
POPULATION_DECISION_ID = "population_structure_exclusion_approval"
DECISION_SCOPES = (
    "controls_unrelated",
    "target_carriers_independent",
    "genomewide_structure_reference",
)
COHORT_COLUMNS = (
    "COHORT_ID",
    "SAMPLE_ID",
    "ROLE",
    "INCLUDED",
    "EXCLUSION_CODE",
    "DECISION_SOURCE",
    "INDEPENDENT_UNIT_ID",
    "FAMILY_REPRESENTATIVE",
)
SUMMARY_COLUMNS = (
    "COHORT_ID",
    "CANDIDATE_COUNT",
    "INCLUDED_COUNT",
    "EXCLUDED_COUNT",
    "INDEPENDENT_UNIT_COUNT",
    "REQUIRED_NONEMPTY",
    "COHORT_STATUS",
)
DECISION_COLUMNS = (
    "MANUAL_DECISION_ID",
    "SAMPLE_ID",
    "SOURCE_STAGE",
    "COHORT_ID",
    "DECISION",
    "REVIEW_ID",
    "REVIEWER_ROLE",
    "REVIEWED_AT",
    "REVIEWED_ARTIFACT_SHA256",
)


class CohortFreezeInputError(ValueError):
    """Signale une configuration ou une entrée de cohorte invalide."""


class CohortFreezeBlockError(RuntimeError):
    """Signale une décision ou une cohorte indispensable non résolue."""


def _build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(prog="freeze-cohorts")
    parser.add_argument("--stage-inputs", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    return parser


def _resolve_input_path(artifact: dict[str, Any], run_dir: Path) -> Path:
    artifact_path = Path(artifact["path"])
    if artifact_path.is_absolute():
        return artifact_path
    if PurePosixPath(artifact["path"]).parts[:1] == ("stages",):
        return run_dir / artifact_path
    return Path.cwd() / artifact_path


def _artifact_by_id(stage_inputs: dict[str, Any], artifact_id: str) -> dict[str, Any]:
    matches = [
        artifact
        for artifact in stage_inputs["artifacts"]
        if artifact["artifact_id"] == artifact_id
    ]
    if len(matches) != 1:
        raise CohortFreezeInputError(f"{artifact_id}_missing_or_ambiguous")
    return matches[0]


def _validated_artifact_path(artifact: dict[str, Any], run_dir: Path) -> Path:
    path = _resolve_input_path(artifact, run_dir)
    if not path.is_file() or sha256_file(path) != artifact["sha256"]:
        raise CohortFreezeInputError("declared_input_modified")
    return path


def _stable_identifier(value: Any, field: str) -> str:
    if (
        not isinstance(value, str)
        or not value
        or not value[0].isalnum()
        or any(not (character.isalnum() or character in "._-") for character in value)
    ):
        raise CohortFreezeInputError(f"invalid_parameter:{field}")
    return value


def _parameters(parameters: dict[str, Any]) -> dict[str, Any]:
    control_group_labels = parameters.get("control_group_labels", ["CONTROL"])
    if (
        not isinstance(control_group_labels, list)
        or not control_group_labels
        or len(control_group_labels) != len(set(control_group_labels))
    ):
        raise CohortFreezeInputError("invalid_parameter:control_group_labels")
    normalized_control_labels = [
        _stable_identifier(value, "control_group_labels")
        for value in control_group_labels
    ]
    required_cohorts = parameters.get(
        "required_nonempty_cohorts", list(DEFAULT_REQUIRED_COHORT_IDS)
    )
    if (
        not isinstance(required_cohorts, list)
        or len(required_cohorts) != len(set(required_cohorts))
        or any(cohort_id not in COHORT_IDS for cohort_id in required_cohorts)
    ):
        raise CohortFreezeInputError("invalid_parameter:required_nonempty_cohorts")
    decision_reviews = parameters.get("decision_reviews", {})
    if not isinstance(decision_reviews, dict) or any(
        decision_id not in {KINSHIP_DECISION_ID, POPULATION_DECISION_ID}
        for decision_id in decision_reviews
    ):
        raise CohortFreezeInputError("invalid_parameter:decision_reviews")
    return {
        "control_group_labels": normalized_control_labels,
        "required_nonempty_cohorts": required_cohorts,
        "decision_reviews": decision_reviews,
    }


def _load_target_metadata(path: Path) -> dict[str, Any]:
    try:
        metadata = yaml.safe_load(path.read_text(encoding="utf-8"))
    except (OSError, yaml.YAMLError) as error:
        raise CohortFreezeInputError("invalid_target_variant_metadata") from error
    if not isinstance(metadata, dict):
        raise CohortFreezeInputError("invalid_target_variant_metadata")
    validate_json_document(metadata, "target_variant_metadata.schema.json")
    return metadata


def _table_rows(path: Path, schema_name: str) -> list[dict[str, Any]]:
    return list(validate_tsv_table(path, schema_name).rows)


def _rows_by_sample(
    rows: list[dict[str, Any]],
    table_name: str,
) -> dict[str, dict[str, Any]]:
    by_sample = {row["SAMPLE_ID"]: row for row in rows}
    if len(by_sample) != len(rows):
        raise CohortFreezeInputError(f"duplicate_sample:{table_name}")
    return by_sample


def _validate_context(
    paths: dict[str, Path],
) -> tuple[
    list[dict[str, Any]],
    dict[str, dict[str, Any]],
    dict[str, dict[str, Any]],
    dict[str, dict[str, Any]],
    dict[str, dict[str, Any]],
    dict[str, dict[str, Any]],
]:
    sample_rows = _table_rows(paths["samples_master"], "samples_master.schema.json")
    genotype_by_sample = _rows_by_sample(
        _table_rows(paths["target_genotype_audit"], "target_genotype_audit.schema.json"),
        "target_genotype_audit",
    )
    qc_by_sample = _rows_by_sample(
        _table_rows(paths["qc_individual_metrics"], "qc_individual_metrics.schema.json"),
        "qc_individual_metrics",
    )
    proposal_by_sample = _rows_by_sample(
        _table_rows(paths["independent_set_proposal"], "independent_set_proposal.schema.json"),
        "independent_set_proposal",
    )
    scores_by_sample = _rows_by_sample(
        _table_rows(paths["population_scores"], "population_scores.schema.json"),
        "population_scores",
    )
    outliers_by_sample = _rows_by_sample(
        _table_rows(paths["population_outliers"], "population_outliers.schema.json"),
        "population_outliers",
    )
    sample_ids = {row["SAMPLE_ID"] for row in sample_rows}
    if set(genotype_by_sample) != sample_ids:
        raise CohortFreezeBlockError("target_genotype_sample_set_incomplete")
    if not set(qc_by_sample) <= sample_ids:
        raise CohortFreezeInputError("qc_sample_set_mismatch")
    panel_ids = set(proposal_by_sample)
    if not panel_ids or panel_ids != set(scores_by_sample) or panel_ids != set(outliers_by_sample):
        raise CohortFreezeInputError("population_and_kinship_sample_set_mismatch")
    if not panel_ids <= sample_ids:
        raise CohortFreezeInputError("panel_sample_absent_from_registry")
    samples_by_id = {row["SAMPLE_ID"]: row for row in sample_rows}
    for sample_id in panel_ids:
        sample = samples_by_id[sample_id]
        proposal = proposal_by_sample[sample_id]
        score = scores_by_sample[sample_id]
        if (
            proposal["FID"] != sample["FID"]
            or proposal["IID"] != sample["IID"]
            or score["FID"] != sample["FID"]
            or score["IID"] != sample["IID"]
            or score["REFERENCE_INCLUDED"] != proposal["PROPOSED_INCLUDED"]
        ):
            raise CohortFreezeInputError("cross_stage_sample_context_mismatch")
        if outliers_by_sample[sample_id]["SAMPLE_SET_ID"] != score["SAMPLE_SET_ID"]:
            raise CohortFreezeInputError("population_sample_set_id_mismatch")
    return (
        sample_rows,
        genotype_by_sample,
        qc_by_sample,
        proposal_by_sample,
        scores_by_sample,
        outliers_by_sample,
    )


def _validate_review(
    decision_id: str,
    proposed_sample_ids: set[str],
    review: Any,
    reviewed_artifact: dict[str, Any],
) -> dict[str, Any] | None:
    if not proposed_sample_ids:
        if review is not None:
            raise CohortFreezeInputError(f"unexpected_review:{decision_id}")
        return None
    if review is None:
        raise CohortFreezeBlockError(f"manual_decision_unresolved:{decision_id}")
    expected_fields = {
        "review_id",
        "reviewer_role",
        "reviewed_at",
        "reviewed_artifact_sha256",
        "approved_exclusion_sample_ids",
    }
    if not isinstance(review, dict) or set(review) != expected_fields:
        raise CohortFreezeInputError(f"invalid_review_contract:{decision_id}")
    review_id = _stable_identifier(review["review_id"], f"{decision_id}.review_id")
    reviewer_role = _stable_identifier(
        review["reviewer_role"], f"{decision_id}.reviewer_role"
    )
    reviewed_at = review["reviewed_at"]
    if not isinstance(reviewed_at, str):
        raise CohortFreezeInputError(f"invalid_reviewed_at:{decision_id}")
    try:
        parsed_date = datetime.fromisoformat(reviewed_at.replace("Z", "+00:00"))
    except ValueError as error:
        raise CohortFreezeInputError(f"invalid_reviewed_at:{decision_id}") from error
    if parsed_date.tzinfo is None:
        raise CohortFreezeInputError(f"reviewed_at_requires_timezone:{decision_id}")
    approved_ids = review["approved_exclusion_sample_ids"]
    if (
        not isinstance(approved_ids, list)
        or len(approved_ids) != len(set(approved_ids))
        or set(approved_ids) != proposed_sample_ids
    ):
        raise CohortFreezeBlockError(f"incomplete_approved_exclusions:{decision_id}")
    reviewed_sha256 = review["reviewed_artifact_sha256"]
    if reviewed_sha256 != reviewed_artifact["sha256"]:
        raise CohortFreezeBlockError(f"reviewed_artifact_mismatch:{decision_id}")
    return {
        "review_id": review_id,
        "reviewer_role": reviewer_role,
        "reviewed_at": reviewed_at,
        "reviewed_artifact_sha256": reviewed_sha256,
        "approved_exclusion_sample_ids": sorted(proposed_sample_ids),
    }


def _quality_key(
    sample_id: str,
    qc_by_sample: dict[str, dict[str, Any]],
) -> tuple[int, float, float, str]:
    qc = qc_by_sample.get(sample_id)
    if qc is None:
        return (2, math.inf, math.inf, sample_id)
    status_rank = {"PASS": 0, "ALERT": 1, "EXCLUDED": 2}[qc["QC_STATUS"]]
    heterozygosity_z = qc["HETEROZYGOSITY_Z"]
    return (
        status_rank,
        float(qc["F_MISS"]),
        abs(float(heterozygosity_z)) if heterozygosity_z is not None else math.inf,
        sample_id,
    )


def _family_representatives(
    sample_ids: set[str],
    samples_by_id: dict[str, dict[str, Any]],
    qc_by_sample: dict[str, dict[str, Any]],
) -> dict[str, str]:
    candidates_by_family: dict[str, list[str]] = defaultdict(list)
    for sample_id in sample_ids:
        candidates_by_family[samples_by_id[sample_id]["FID"]].append(sample_id)
    return {
        family_id: min(candidates, key=lambda item: _quality_key(item, qc_by_sample))
        for family_id, candidates in candidates_by_family.items()
    }


def _carrier_status(genotype: str, alt: str) -> bool:
    alleles = genotype.replace("|", "/").split("/")
    if len(alleles) != 2:
        raise CohortFreezeInputError("invalid_accepted_target_genotype")
    return alt in alleles


def _excluded_row(
    cohort_id: str,
    sample_id: str,
    role: str,
    exclusion_code: str,
    decision_source: str,
) -> dict[str, str]:
    return {
        "COHORT_ID": cohort_id,
        "SAMPLE_ID": sample_id,
        "ROLE": role,
        "INCLUDED": "false",
        "EXCLUSION_CODE": exclusion_code,
        "DECISION_SOURCE": decision_source,
        "INDEPENDENT_UNIT_ID": "",
        "FAMILY_REPRESENTATIVE": "false",
    }


def _included_row(
    cohort_id: str,
    sample_id: str,
    role: str,
    decision_source: str,
    *,
    independent_unit_id: str = "",
    family_representative: bool = False,
) -> dict[str, str]:
    return {
        "COHORT_ID": cohort_id,
        "SAMPLE_ID": sample_id,
        "ROLE": role,
        "INCLUDED": "true",
        "EXCLUSION_CODE": "",
        "DECISION_SOURCE": decision_source,
        "INDEPENDENT_UNIT_ID": independent_unit_id,
        "FAMILY_REPRESENTATIVE": str(family_representative).lower(),
    }


def _cohort_rows(
    sample_rows: list[dict[str, Any]],
    genotype_by_sample: dict[str, dict[str, Any]],
    qc_by_sample: dict[str, dict[str, Any]],
    proposal_by_sample: dict[str, dict[str, Any]],
    scores_by_sample: dict[str, dict[str, Any]],
    approved_kinship_ids: set[str],
    approved_population_ids: set[str],
    control_group_labels: set[str],
    alt: str,
) -> tuple[list[dict[str, str]], dict[str, set[str]]]:
    samples_by_id = {row["SAMPLE_ID"]: row for row in sample_rows}
    carrier_ids = {
        sample_id
        for sample_id, genotype in genotype_by_sample.items()
        if _carrier_status(genotype["GENOTYPE"], alt)
    }
    noncarrier_ids = set(samples_by_id) - carrier_ids
    carrier_family_ids = {samples_by_id[sample_id]["FID"] for sample_id in carrier_ids}
    control_candidates = {
        sample_id
        for sample_id in noncarrier_ids
        if samples_by_id[sample_id]["GROUP_LABEL"] in control_group_labels
        and samples_by_id[sample_id]["INCLUDE_GENOMEWIDE"]
    }
    independent_carrier_candidates = {
        sample_id
        for sample_id in carrier_ids
        if samples_by_id[sample_id]["INCLUDE_TARGET_CHROMOSOME"]
    }
    control_eligible = {
        sample_id
        for sample_id in control_candidates
        if sample_id in proposal_by_sample
        and proposal_by_sample[sample_id]["PROPOSED_INCLUDED"]
        and sample_id not in approved_kinship_ids
        and sample_id not in approved_population_ids
    }
    independent_carrier_eligible = {
        sample_id
        for sample_id in independent_carrier_candidates
        if sample_id in proposal_by_sample
        and proposal_by_sample[sample_id]["PROPOSED_INCLUDED"]
        and sample_id not in approved_kinship_ids
        and sample_id not in approved_population_ids
    }
    control_representatives = _family_representatives(
        control_eligible, samples_by_id, qc_by_sample
    )
    carrier_representatives = _family_representatives(
        independent_carrier_eligible, samples_by_id, qc_by_sample
    )
    controls_included = set(control_representatives.values())
    independent_carriers_included = set(carrier_representatives.values())
    all_carrier_representatives = _family_representatives(
        carrier_ids, samples_by_id, qc_by_sample
    )
    structure_candidates = {
        sample_id
        for sample_id, score in scores_by_sample.items()
        if score["REFERENCE_INCLUDED"]
    }
    structure_included = structure_candidates - approved_population_ids - approved_kinship_ids
    structure_family_representatives = _family_representatives(
        structure_included, samples_by_id, qc_by_sample
    )
    candidate_sets = {
        "controls_unrelated": control_candidates,
        "target_carriers_all": carrier_ids,
        "target_carriers_independent": independent_carrier_candidates,
        "family_noncarriers": {
            sample_id
            for sample_id in noncarrier_ids
            if samples_by_id[sample_id]["FID"] in carrier_family_ids
            and samples_by_id[sample_id]["INCLUDE_TARGET_CHROMOSOME"]
        },
        "affected_noncarriers": {
            sample_id
            for sample_id in noncarrier_ids
            if samples_by_id[sample_id]["CLINICAL_STATUS"] == "AFFECTED"
            and samples_by_id[sample_id]["INCLUDE_TARGET_CHROMOSOME"]
        },
        "target_chromosome_all_qc": {
            sample_id
            for sample_id, sample in samples_by_id.items()
            if sample["INCLUDE_TARGET_CHROMOSOME"]
        },
        "genomewide_structure_reference": structure_candidates,
    }
    included_sets = {
        "controls_unrelated": controls_included,
        "target_carriers_all": carrier_ids,
        "target_carriers_independent": independent_carriers_included,
        "family_noncarriers": candidate_sets["family_noncarriers"],
        "affected_noncarriers": candidate_sets["affected_noncarriers"],
        "target_chromosome_all_qc": candidate_sets["target_chromosome_all_qc"],
        "genomewide_structure_reference": structure_included,
    }
    rows: list[dict[str, str]] = []
    for cohort_id in COHORT_IDS:
        for sample in sample_rows:
            sample_id = sample["SAMPLE_ID"]
            carrier = sample_id in carrier_ids
            role = "CARRIER" if carrier else "NON_CARRIER"
            if cohort_id == "controls_unrelated":
                role = "CONTROL" if sample_id in control_candidates else role
                if sample_id in controls_included:
                    rows.append(
                        _included_row(
                            cohort_id,
                            sample_id,
                            "CONTROL",
                            "stage_09_family_unit_selection",
                            independent_unit_id=f"control_{sample_id}",
                            family_representative=True,
                        )
                    )
                elif sample["GROUP_LABEL"] not in control_group_labels:
                    rows.append(_excluded_row(cohort_id, sample_id, role, "not_explicit_control_group", "explicit_control_group_rule"))
                elif carrier:
                    rows.append(_excluded_row(cohort_id, sample_id, role, "target_variant_carrier", "explicit_target_genotype"))
                elif not sample["INCLUDE_GENOMEWIDE"]:
                    rows.append(_excluded_row(cohort_id, sample_id, role, "genomewide_not_included", "registry_inclusion_flags"))
                elif sample_id not in proposal_by_sample:
                    rows.append(_excluded_row(cohort_id, sample_id, role, "preliminary_qc_excluded", "stage_05_preliminary_qc"))
                elif sample_id in approved_kinship_ids or not proposal_by_sample[sample_id]["PROPOSED_INCLUDED"]:
                    rows.append(_excluded_row(cohort_id, sample_id, role, "approved_kinship_exclusion", "stage_07_independent_proposal"))
                elif sample_id in approved_population_ids:
                    rows.append(_excluded_row(cohort_id, sample_id, role, "approved_population_exclusion", "stage_08_population_outlier"))
                else:
                    rows.append(_excluded_row(cohort_id, sample_id, role, "same_family_control_represented", "stage_09_family_unit_selection"))
            elif cohort_id == "target_carriers_all":
                if carrier and sample["INCLUDE_TARGET_CHROMOSOME"]:
                    family_id = sample["FID"]
                    rows.append(
                        _included_row(
                            cohort_id,
                            sample_id,
                            "CARRIER",
                            "explicit_target_genotype",
                            independent_unit_id=f"carrier_family_{family_id}",
                            family_representative=all_carrier_representatives[family_id] == sample_id,
                        )
                    )
                elif not carrier:
                    rows.append(_excluded_row(cohort_id, sample_id, role, "target_variant_noncarrier", "explicit_target_genotype"))
                else:
                    rows.append(_excluded_row(cohort_id, sample_id, role, "target_chromosome_not_included", "registry_inclusion_flags"))
            elif cohort_id == "target_carriers_independent":
                if sample_id in independent_carriers_included:
                    rows.append(
                        _included_row(
                            cohort_id,
                            sample_id,
                            "CARRIER",
                            "stage_09_family_unit_selection",
                            independent_unit_id=f"carrier_family_{sample['FID']}",
                            family_representative=True,
                        )
                    )
                elif not carrier:
                    rows.append(_excluded_row(cohort_id, sample_id, role, "target_variant_noncarrier", "explicit_target_genotype"))
                elif not sample["INCLUDE_TARGET_CHROMOSOME"]:
                    rows.append(_excluded_row(cohort_id, sample_id, role, "target_chromosome_not_included", "registry_inclusion_flags"))
                elif sample_id not in proposal_by_sample:
                    rows.append(_excluded_row(cohort_id, sample_id, role, "preliminary_qc_excluded", "stage_05_preliminary_qc"))
                elif sample_id in approved_kinship_ids or not proposal_by_sample[sample_id]["PROPOSED_INCLUDED"]:
                    rows.append(_excluded_row(cohort_id, sample_id, "RELATED", "approved_kinship_exclusion", "stage_07_independent_proposal"))
                elif sample_id in approved_population_ids:
                    rows.append(_excluded_row(cohort_id, sample_id, role, "approved_population_exclusion", "stage_08_population_outlier"))
                else:
                    rows.append(_excluded_row(cohort_id, sample_id, role, "same_family_carrier_represented", "stage_09_family_unit_selection"))
            elif cohort_id == "family_noncarriers":
                if sample_id in included_sets[cohort_id]:
                    rows.append(_included_row(cohort_id, sample_id, "NON_CARRIER", "explicit_target_genotype", independent_unit_id=f"carrier_family_{sample['FID']}"))
                elif carrier:
                    rows.append(_excluded_row(cohort_id, sample_id, role, "target_variant_carrier", "explicit_target_genotype"))
                elif sample["FID"] not in carrier_family_ids:
                    rows.append(_excluded_row(cohort_id, sample_id, role, "family_without_target_carrier", "explicit_family_membership"))
                else:
                    rows.append(_excluded_row(cohort_id, sample_id, role, "target_chromosome_not_included", "registry_inclusion_flags"))
            elif cohort_id == "affected_noncarriers":
                if sample_id in included_sets[cohort_id]:
                    rows.append(_included_row(cohort_id, sample_id, "NON_CARRIER", "clinical_status_and_explicit_genotype"))
                elif carrier:
                    rows.append(_excluded_row(cohort_id, sample_id, role, "target_variant_carrier", "explicit_target_genotype"))
                elif sample["CLINICAL_STATUS"] != "AFFECTED":
                    rows.append(_excluded_row(cohort_id, sample_id, role, "not_affected", "explicit_clinical_status"))
                else:
                    rows.append(_excluded_row(cohort_id, sample_id, role, "target_chromosome_not_included", "registry_inclusion_flags"))
            elif cohort_id == "target_chromosome_all_qc":
                if sample_id in included_sets[cohort_id]:
                    rows.append(_included_row(cohort_id, sample_id, role, "accepted_target_genotype_and_registry"))
                else:
                    rows.append(_excluded_row(cohort_id, sample_id, role, "target_chromosome_not_included", "registry_inclusion_flags"))
            else:
                if sample_id in structure_included:
                    family_representative = structure_family_representatives[sample["FID"]] == sample_id
                    rows.append(
                        _included_row(
                            cohort_id,
                            sample_id,
                            "OTHER",
                            "stage_08_population_reference",
                            independent_unit_id=proposal_by_sample[sample_id]["COMPONENT_ID"],
                            family_representative=family_representative,
                        )
                    )
                elif sample_id not in scores_by_sample:
                    rows.append(_excluded_row(cohort_id, sample_id, role, "preliminary_qc_excluded", "stage_05_preliminary_qc"))
                elif sample_id in approved_kinship_ids or not proposal_by_sample[sample_id]["PROPOSED_INCLUDED"]:
                    rows.append(_excluded_row(cohort_id, sample_id, "RELATED", "approved_kinship_exclusion", "stage_07_independent_proposal"))
                elif sample_id in approved_population_ids:
                    rows.append(_excluded_row(cohort_id, sample_id, role, "approved_population_exclusion", "stage_08_population_outlier"))
                else:
                    rows.append(_excluded_row(cohort_id, sample_id, role, "not_population_reference", "stage_08_population_reference"))
    return rows, candidate_sets


def _summary_rows(
    cohort_rows: list[dict[str, str]],
    candidate_sets: dict[str, set[str]],
    required_cohorts: set[str],
) -> list[dict[str, str]]:
    rows = []
    for cohort_id in COHORT_IDS:
        cohort_memberships = [row for row in cohort_rows if row["COHORT_ID"] == cohort_id]
        included = [row for row in cohort_memberships if row["INCLUDED"] == "true"]
        if cohort_id in required_cohorts and not included:
            raise CohortFreezeBlockError(f"required_cohort_empty:{cohort_id}")
        units = {row["INDEPENDENT_UNIT_ID"] for row in included if row["INDEPENDENT_UNIT_ID"]}
        rows.append(
            {
                "COHORT_ID": cohort_id,
                "CANDIDATE_COUNT": str(len(candidate_sets[cohort_id])),
                "INCLUDED_COUNT": str(len(included)),
                "EXCLUDED_COUNT": str(len(candidate_sets[cohort_id]) - len(included)),
                "INDEPENDENT_UNIT_COUNT": str(len(units)),
                "REQUIRED_NONEMPTY": str(cohort_id in required_cohorts).lower(),
                "COHORT_STATUS": "POPULATED" if included else "EMPTY_ALLOWED",
            }
        )
    return rows


def _decision_rows(
    reviews: dict[str, dict[str, Any] | None],
) -> list[dict[str, str]]:
    rows = []
    source_stages = {
        KINSHIP_DECISION_ID: "07_infer_kinship",
        POPULATION_DECISION_ID: "08_analyze_population_structure",
    }
    for decision_id in (KINSHIP_DECISION_ID, POPULATION_DECISION_ID):
        review = reviews[decision_id]
        if review is None:
            continue
        for sample_id in review["approved_exclusion_sample_ids"]:
            for cohort_id in DECISION_SCOPES:
                rows.append(
                    {
                        "MANUAL_DECISION_ID": decision_id,
                        "SAMPLE_ID": sample_id,
                        "SOURCE_STAGE": source_stages[decision_id],
                        "COHORT_ID": cohort_id,
                        "DECISION": "APPROVED_EXCLUSION",
                        "REVIEW_ID": review["review_id"],
                        "REVIEWER_ROLE": review["reviewer_role"],
                        "REVIEWED_AT": review["reviewed_at"],
                        "REVIEWED_ARTIFACT_SHA256": review["reviewed_artifact_sha256"],
                    }
                )
    return rows


def _audit_parameters(parameters: dict[str, Any]) -> dict[str, Any]:
    """Retire les identifiants individuels des paramètres publiés dans l'audit."""
    review_summary = {}
    for decision_id, review in parameters["decision_reviews"].items():
        review_summary[decision_id] = {
            "review_id": review["review_id"],
            "reviewer_role": review["reviewer_role"],
            "reviewed_at": review["reviewed_at"],
            "reviewed_artifact_sha256": review["reviewed_artifact_sha256"],
            "approved_exclusion_count": len(review["approved_exclusion_sample_ids"]),
        }
    return {
        "control_group_labels": parameters["control_group_labels"],
        "required_nonempty_cohorts": parameters["required_nonempty_cohorts"],
        "decision_reviews": review_summary,
    }


def _write_tsv(path: Path, columns: tuple[str, ...], rows: list[dict[str, str]]) -> None:
    with path.open("w", encoding="utf-8", newline="") as output_file:
        writer = csv.DictWriter(output_file, fieldnames=columns, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows({column: row[column] for column in columns} for row in rows)


def _write_keep_files(
    keep_dir: Path,
    cohort_rows: list[dict[str, str]],
    samples_by_id: dict[str, dict[str, Any]],
) -> dict[str, Path]:
    keep_dir.mkdir(parents=True, exist_ok=True)
    paths = {}
    for cohort_id in COHORT_IDS:
        path = keep_dir / f"{cohort_id}.keep"
        included_sample_ids = {
            row["SAMPLE_ID"]
            for row in cohort_rows
            if row["COHORT_ID"] == cohort_id and row["INCLUDED"] == "true"
        }
        path.write_text(
            "".join(
                f"{sample['FID']}\t{sample['IID']}\n"
                for sample_id, sample in samples_by_id.items()
                if sample_id in included_sample_ids
            ),
            encoding="utf-8",
        )
        paths[cohort_id] = path
    return paths


def _output_artifacts(
    stage_inputs: dict[str, Any],
    metadata: dict[str, Any],
    sample_set_id: str | None,
    variant_set_id: str | None,
    contracts: list[tuple[str, str, Path, str | None, str]],
) -> list[dict[str, Any]]:
    producer_stage = f"{stage_inputs['stage_id']}_{stage_inputs['stage_name']}"
    published_dir = PurePosixPath(stage_inputs["published_output_dir"])
    artifacts = []
    for artifact_id, artifact_type, path, schema_name, sensitivity in contracts:
        relative_name = path.relative_to(path.parents[1]) if path.parent.name == "keep" else Path(path.name)
        artifacts.append(
            build_file_artifact(
                physical_path=path,
                published_path=str(published_dir / relative_name.as_posix()),
                artifact_id=artifact_id,
                artifact_type=artifact_type,
                media_type=(
                    "application/json"
                    if path.suffix == ".json"
                    else "text/plain"
                    if path.suffix == ".keep"
                    else "text/tab-separated-values"
                ),
                producer_stage=producer_stage,
                producer_signature=stage_inputs["signature"],
                schema_name=schema_name,
                schema_version="1.0.0" if schema_name else None,
                assembly=metadata["assembly"],
                sample_set_id=sample_set_id,
                variant_set_id=variant_set_id,
                sensitivity=sensitivity,
            )
        )
    return artifacts


def execute(stage_inputs_path: Path, output_dir: Path) -> int:
    """Valide les décisions humaines et publie les cohortes immuables."""
    started_at = utc_now()
    started_clock = monotonic()
    stage_inputs = read_json(stage_inputs_path)
    validate_json_document(stage_inputs, "stage_inputs.schema.json")
    run_dir = output_dir.parent.parent
    load_pipeline_config(run_dir / "config.resolved.yaml")
    parameters = _parameters(stage_inputs["parameters"])
    artifact_ids = (
        "config_input_target_variant_metadata",
        "samples_master",
        "target_genotype_audit",
        "qc_individual_metrics",
        "independent_set_proposal",
        "population_scores",
        "population_outliers",
    )
    input_artifacts = {artifact_id: _artifact_by_id(stage_inputs, artifact_id) for artifact_id in artifact_ids}
    paths = {artifact_id: _validated_artifact_path(artifact, run_dir) for artifact_id, artifact in input_artifacts.items()}
    metadata = _load_target_metadata(paths["config_input_target_variant_metadata"])
    (
        sample_rows,
        genotype_by_sample,
        qc_by_sample,
        proposal_by_sample,
        scores_by_sample,
        outliers_by_sample,
    ) = _validate_context(paths)
    proposed_kinship_ids = {
        sample_id
        for sample_id, row in proposal_by_sample.items()
        if not row["PROPOSED_INCLUDED"]
    }
    proposed_population_ids = {
        sample_id
        for sample_id, row in outliers_by_sample.items()
        if row["PROPOSED_EXCLUSION"]
    }
    reviews_config = parameters["decision_reviews"]
    reviews = {
        KINSHIP_DECISION_ID: _validate_review(
            KINSHIP_DECISION_ID,
            proposed_kinship_ids,
            reviews_config.get(KINSHIP_DECISION_ID),
            input_artifacts["independent_set_proposal"],
        ),
        POPULATION_DECISION_ID: _validate_review(
            POPULATION_DECISION_ID,
            proposed_population_ids,
            reviews_config.get(POPULATION_DECISION_ID),
            input_artifacts["population_outliers"],
        ),
    }
    cohort_rows, candidate_sets = _cohort_rows(
        sample_rows,
        genotype_by_sample,
        qc_by_sample,
        proposal_by_sample,
        scores_by_sample,
        proposed_kinship_ids,
        proposed_population_ids,
        set(parameters["control_group_labels"]),
        metadata["alt"],
    )
    required_cohorts = set(parameters["required_nonempty_cohorts"])
    summary_rows = _summary_rows(cohort_rows, candidate_sets, required_cohorts)
    decision_rows = _decision_rows(reviews)

    output_dir.mkdir(parents=True, exist_ok=True)
    cohorts_path = output_dir / "cohorts.frozen.tsv"
    summary_path = output_dir / "cohort_summary.tsv"
    decisions_path = output_dir / "cohort_decisions.tsv"
    _write_tsv(cohorts_path, COHORT_COLUMNS, cohort_rows)
    _write_tsv(summary_path, SUMMARY_COLUMNS, summary_rows)
    _write_tsv(decisions_path, DECISION_COLUMNS, decision_rows)
    validate_cohorts_frozen(cohorts_path, {row["SAMPLE_ID"] for row in sample_rows})
    validate_tsv_table(summary_path, "cohort_summary.schema.json")
    validate_tsv_table(decisions_path, "cohort_decisions.schema.json")
    keep_paths = _write_keep_files(
        output_dir / "keep",
        cohort_rows,
        {row["SAMPLE_ID"]: row for row in sample_rows},
    )

    report_path = output_dir / "cohort_freeze_report.json"
    report = {
        "schema_version": "1.0.0",
        "method_id": "explicit_target_genotype_cohort_freeze_v1",
        "target_variant_id": metadata["project_variant_id"],
        "sample_count": len(sample_rows),
        "cohort_count": len(COHORT_IDS),
        "control_group_labels": parameters["control_group_labels"],
        "required_nonempty_cohorts": parameters["required_nonempty_cohorts"],
        "cohort_included_counts": {
            row["COHORT_ID"]: int(row["INCLUDED_COUNT"]) for row in summary_rows
        },
        "approved_kinship_exclusion_count": len(proposed_kinship_ids),
        "approved_population_exclusion_count": len(proposed_population_ids),
        "manual_validation_required": False,
    }
    atomic_write_json(report_path, report)
    contracts = [
        ("cohorts_frozen", "cohorts_frozen", cohorts_path, "cohorts_frozen.schema.json", "sensitive_genetic"),
        ("cohort_summary", "cohort_summary", summary_path, "cohort_summary.schema.json", "internal"),
        ("cohort_decisions", "cohort_decisions", decisions_path, "cohort_decisions.schema.json", "sensitive_genetic"),
        ("cohort_freeze_report", "cohort_freeze_report", report_path, None, "internal"),
    ]
    contracts.extend(
        (
            f"cohort_keep_{cohort_id}",
            "plink_keep",
            keep_paths[cohort_id],
            None,
            "sensitive_genetic",
        )
        for cohort_id in COHORT_IDS
    )
    sample_set_id = input_artifacts["samples_master"].get("sample_set_id")
    variant_set_id = input_artifacts["target_genotype_audit"].get("variant_set_id")
    output_artifacts = _output_artifacts(
        stage_inputs, metadata, sample_set_id, variant_set_id, contracts
    )
    stage_outputs = {
        "schema_version": "1.0.0",
        "run_id": stage_inputs["run_id"],
        "stage_id": stage_inputs["stage_id"],
        "stage_name": stage_inputs["stage_name"],
        "signature": stage_inputs["signature"],
        "artifacts": output_artifacts,
    }
    validate_json_document(stage_outputs, "stage_outputs.schema.json")
    atomic_write_json(output_dir / "stage_outputs.json", stage_outputs)
    audit = {
        "schema_version": "1.0.0",
        "run_id": stage_inputs["run_id"],
        "stage_id": stage_inputs["stage_id"],
        "stage_name": stage_inputs["stage_name"],
        "method_id": "explicit_target_genotype_cohort_freeze_v1",
        "signature": stage_inputs["signature"],
        "started_at": started_at,
        "completed_at": utc_now(),
        "duration_seconds": monotonic() - started_clock,
        "inputs": list(input_artifacts.values()),
        "outputs": output_artifacts,
        "parameters": _audit_parameters(parameters),
        "tools": [],
        "counts": {
            "samples": len(sample_rows),
            "cohorts": len(COHORT_IDS),
            "memberships": len(cohort_rows),
            "manual_decision_rows": len(decision_rows),
            "approved_kinship_exclusions": len(proposed_kinship_ids),
            "approved_population_exclusions": len(proposed_population_ids),
        },
        "metrics": {
            "included_counts": report["cohort_included_counts"],
        },
        "exclusions": [],
        "warnings": [],
        "checks": [
            {"check": "target_genotypes_explicit", "status": "PASS"},
            {"check": "cross_stage_sample_sets", "status": "PASS"},
            {"check": "manual_decision_reviews", "status": "PASS"},
            {"check": "required_cohorts_nonempty", "status": "PASS"},
            {"check": "plink_keep_files", "status": "PASS"},
        ],
        "known_limits": [
            "Les groupes témoins proviennent d'une liste explicite de GROUP_LABEL configurée ; ils ne définissent jamais le génotype cible.",
            "Une approbation accepte exactement les exclusions proposées par l'artefact dont l'empreinte est revue ; changer de représentants exige un nouveau calcul amont.",
            "Les exclusions d'apparentement et de structure ne s'appliquent qu'aux cohortes analytiques indépendantes, pas à target_carriers_all ni à target_chromosome_all_qc.",
            "Le gel ne remplace pas le QC final propre à chaque cohorte de l'étape 10.",
        ],
        "expected_visualizations": [
            "cohort_membership_counts",
            "cohort_exclusion_reasons",
            "cohort_overlap",
        ],
        "manual_validation_required": False,
    }
    validate_json_document(audit, "stage_audit.schema.json")
    atomic_write_json(output_dir / "audit.json", audit)
    (output_dir / "checksums.sha256").write_text(
        "".join(
            f"{artifact['sha256']}  {Path(str(artifact['path'])).name}\n"
            for artifact in output_artifacts
        ),
        encoding="utf-8",
    )
    return 0


def main(arguments: Sequence[str] | None = None) -> int:
    """Point d'entrée autonome respectant les codes de retour V2."""
    parser = _build_parser()
    parsed_arguments = parser.parse_args(arguments)
    try:
        return execute(parsed_arguments.stage_inputs, parsed_arguments.output_dir)
    except CohortFreezeBlockError as error:
        sys.stderr.write(f"{error}\n")
        return 4
    except (
        CohortFreezeInputError,
        DocumentValidationError,
        TableValidationError,
        OSError,
        ValueError,
        json.JSONDecodeError,
    ) as error:
        parser.error(str(error))
    return 2


if __name__ == "__main__":
    raise SystemExit(main())
