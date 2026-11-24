"""D22 : oracle déclaratif des sept cohortes, sans helper de production 09."""

from __future__ import annotations

import math
from typing import Any

from effet_fondateur.calibration.structure_bench import COHORTS


def cohort_expectation(registry: list[dict[str, str]], genotypes: list[dict[str, str]], qc_rows: list[dict[str, str]],
                       proposals: list[dict[str, str]], scores: list[dict[str, str]], outliers: list[dict[str, str]],
                       parameters: dict[str, Any], alt: str) -> dict[str, Any]:
    """Reconstruit décisions, sources et unités depuis les seules entrées 09.

    Les génotypes cibles sont explicites. L'ordre des exclusions correspond
    au contrat de rôle, jamais à une appartenance de strate privée.
    """
    samples = {row["SAMPLE_ID"]: row for row in registry}
    qc = {row["SAMPLE_ID"]: row for row in qc_rows}
    proposal = {row["SAMPLE_ID"]: row for row in proposals}
    carrier = {row["SAMPLE_ID"] for row in genotypes if alt in row["GENOTYPE"].split("/")}
    families = {samples[sample]["FID"] for sample in carrier}
    controls = set(parameters["control_group_labels"])
    kin = {row["SAMPLE_ID"] for row in proposals if row["PROPOSED_INCLUDED"] == "false"}
    population = {row["SAMPLE_ID"] for row in outliers if row["PROPOSED_EXCLUSION"] == "true"}
    reference = {row["SAMPLE_ID"] for row in scores if row["REFERENCE_INCLUDED"] == "true"}
    all_scores = {row["SAMPLE_ID"] for row in scores}

    def quality(sample: str) -> tuple[int, float, float, str]:
        row = qc.get(sample)
        if row is None:
            return (2, math.inf, math.inf, sample)
        return (dict(PASS=0, ALERT=1, EXCLUDED=2)[row["QC_STATUS"]], float(row["F_MISS"]),
                abs(float(row["HETEROZYGOSITY_Z"])) if row["HETEROZYGOSITY_Z"] else math.inf, sample)

    def representatives(candidates: set[str]) -> dict[str, str]:
        result = {}
        for sample in sorted(candidates, key=quality):
            result.setdefault(samples[sample]["FID"], sample)
        return result

    control_candidates = {sample for sample, row in samples.items() if sample not in carrier and row["GROUP_LABEL"] in controls and row["INCLUDE_GENOMEWIDE"] == "true"}
    carrier_candidates = {sample for sample in carrier if samples[sample]["INCLUDE_TARGET_CHROMOSOME"] == "true"}
    independent = set(proposal) - kin - population
    control_reps = representatives(control_candidates & independent)
    carrier_reps = representatives(carrier_candidates & independent)
    all_carrier_reps = representatives(carrier)
    structure = reference - kin - population
    structure_reps = representatives(structure)
    candidates = dict(controls_unrelated=control_candidates, target_carriers_all=carrier, target_carriers_independent=carrier_candidates,
                      family_noncarriers={sample for sample, row in samples.items() if sample not in carrier and row["FID"] in families and row["INCLUDE_TARGET_CHROMOSOME"] == "true"},
                      affected_noncarriers={sample for sample, row in samples.items() if sample not in carrier and row["CLINICAL_STATUS"] == "AFFECTED" and row["INCLUDE_TARGET_CHROMOSOME"] == "true"},
                      target_chromosome_all_qc={sample for sample, row in samples.items() if row["INCLUDE_TARGET_CHROMOSOME"] == "true"}, genomewide_structure_reference=reference)
    memberships = dict(controls_unrelated=set(control_reps.values()), target_carriers_all=carrier_candidates,
                       target_carriers_independent=set(carrier_reps.values()), family_noncarriers=candidates["family_noncarriers"],
                       affected_noncarriers=candidates["affected_noncarriers"], target_chromosome_all_qc=candidates["target_chromosome_all_qc"], genomewide_structure_reference=structure)
    rows = []
    for cohort in COHORTS:
        for sample, metadata in samples.items():
            is_carrier = sample in carrier
            role = "CARRIER" if is_carrier else "NON_CARRIER"
            family = metadata["FID"]
            unit, representative = "", False
            included = sample in memberships[cohort]
            if included:
                if cohort == "controls_unrelated":
                    role, source, unit, representative = "CONTROL", "stage_09_family_unit_selection", "control_" + sample, True
                elif cohort in {"target_carriers_all", "target_carriers_independent"}:
                    source = "explicit_target_genotype" if cohort == "target_carriers_all" else "stage_09_family_unit_selection"
                    unit, representative = "carrier_family_" + family, all_carrier_reps[family] == sample if cohort == "target_carriers_all" else True
                elif cohort == "family_noncarriers":
                    source, unit = "explicit_target_genotype", "carrier_family_" + family
                elif cohort == "affected_noncarriers":
                    source = "clinical_status_and_explicit_genotype"
                elif cohort == "target_chromosome_all_qc":
                    source = "accepted_target_genotype_and_registry"
                else:
                    role, source, unit, representative = "OTHER", "stage_08_population_reference", proposal[sample]["COMPONENT_ID"], structure_reps[family] == sample
                code = ""
            else:
                if cohort == "controls_unrelated":
                    role = "CONTROL" if sample in control_candidates else role
                    conditions = [(metadata["GROUP_LABEL"] not in controls, "not_explicit_control_group", "explicit_control_group_rule"),
                                  (is_carrier, "target_variant_carrier", "explicit_target_genotype"),
                                  (metadata["INCLUDE_GENOMEWIDE"] != "true", "genomewide_not_included", "registry_inclusion_flags"),
                                  (sample not in proposal, "preliminary_qc_excluded", "stage_05_preliminary_qc"),
                                  (sample in kin, "approved_kinship_exclusion", "stage_07_independent_proposal"),
                                  (sample in population, "approved_population_exclusion", "stage_08_population_outlier"),
                                  (True, "same_family_control_represented", "stage_09_family_unit_selection")]
                elif cohort in {"target_carriers_all", "target_carriers_independent"}:
                    conditions = [(not is_carrier, "target_variant_noncarrier", "explicit_target_genotype"),
                                  (metadata["INCLUDE_TARGET_CHROMOSOME"] != "true", "target_chromosome_not_included", "registry_inclusion_flags"),
                                  (sample not in proposal, "preliminary_qc_excluded", "stage_05_preliminary_qc"),
                                  (sample in kin, "approved_kinship_exclusion", "stage_07_independent_proposal"),
                                  (sample in population, "approved_population_exclusion", "stage_08_population_outlier"),
                                  (True, "same_family_carrier_represented", "stage_09_family_unit_selection")]
                elif cohort == "family_noncarriers":
                    conditions = [(is_carrier, "target_variant_carrier", "explicit_target_genotype"),
                                  (family not in families, "family_without_target_carrier", "explicit_family_membership"),
                                  (True, "target_chromosome_not_included", "registry_inclusion_flags")]
                elif cohort == "affected_noncarriers":
                    conditions = [(is_carrier, "target_variant_carrier", "explicit_target_genotype"),
                                  (metadata["CLINICAL_STATUS"] != "AFFECTED", "not_affected", "explicit_clinical_status"),
                                  (True, "target_chromosome_not_included", "registry_inclusion_flags")]
                elif cohort == "target_chromosome_all_qc":
                    conditions = [(True, "target_chromosome_not_included", "registry_inclusion_flags")]
                else:
                    conditions = [(sample not in all_scores, "preliminary_qc_excluded", "stage_05_preliminary_qc"),
                                  (sample in kin, "approved_kinship_exclusion", "stage_07_independent_proposal"),
                                  (sample in population, "approved_population_exclusion", "stage_08_population_outlier"),
                                  (True, "not_population_reference", "stage_08_population_reference")]
                _, code, source = next(condition for condition in conditions if condition[0])
                if code == "approved_kinship_exclusion" and cohort != "controls_unrelated":
                    role = "RELATED"
            rows.append(dict(COHORT_ID=cohort, SAMPLE_ID=sample, ROLE=role, INCLUDED=str(included).lower(), EXCLUSION_CODE=code,
                             DECISION_SOURCE=source, INDEPENDENT_UNIT_ID=unit, FAMILY_REPRESENTATIVE=str(representative).lower()))
    summary, keeps = [], {}
    for cohort in COHORTS:
        included_rows = [row for row in rows if row["COHORT_ID"] == cohort and row["INCLUDED"] == "true"]
        summary.append(dict(COHORT_ID=cohort, CANDIDATE_COUNT=len(candidates[cohort]), INCLUDED_COUNT=len(included_rows),
                            EXCLUDED_COUNT=len(candidates[cohort]) - len(included_rows), INDEPENDENT_UNIT_COUNT=len({row["INDEPENDENT_UNIT_ID"] for row in included_rows if row["INDEPENDENT_UNIT_ID"]}),
                            REQUIRED_NONEMPTY=str(cohort in parameters["required_nonempty_cohorts"]).lower(), COHORT_STATUS="POPULATED" if included_rows else "EMPTY_ALLOWED"))
        keeps[cohort] = "".join(f'{row["FID"]}\t{row["IID"]}\n' for sample, row in samples.items() if sample in memberships[cohort])
    decisions = []
    for decision_id, source_stage, excluded in (("kinship_exclusion_approval", "07_infer_kinship", kin), ("population_structure_exclusion_approval", "08_analyze_population_structure", population)):
        review = parameters["decision_reviews"].get(decision_id)
        if not excluded:
            continue
        if review is None:
            raise ValueError("d22_oracle_review_missing")
        for sample in sorted(excluded):
            for cohort in ("controls_unrelated", "target_carriers_independent", "genomewide_structure_reference"):
                decisions.append(dict(MANUAL_DECISION_ID=decision_id, SAMPLE_ID=sample, SOURCE_STAGE=source_stage, COHORT_ID=cohort, DECISION="APPROVED_EXCLUSION",
                                      REVIEW_ID=review["review_id"], REVIEWER_ROLE=review["reviewer_role"], REVIEWED_AT=review["reviewed_at"], REVIEWED_ARTIFACT_SHA256=review["reviewed_artifact_sha256"]))
    return dict(rows=rows, summary=summary, decisions=decisions, keeps=keeps)
