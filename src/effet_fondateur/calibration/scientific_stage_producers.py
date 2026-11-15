"""Producteurs scientifiques propres aux étapes R13.

Chaque producteur transforme des preuves closes en observations métriques. Les
étapes absentes de ``IMPLEMENTED_STAGE_PRODUCERS`` restent explicitement
bloquantes pour G5.
"""

from __future__ import annotations

from collections.abc import Mapping
import math
from typing import Any

from effet_fondateur.calibration.contracts import STAGE_IDS
from effet_fondateur.calibration.scientific_metric_producers import (
    cluster_statistic_observation,
    exact_observation,
    inference_method_for_rule,
    proportion_observation,
)
from effet_fondateur.calibration.scientific_statistics import (
    DEFAULT_BOOTSTRAP_REPLICATES,
)


IMPLEMENTED_STAGE_PRODUCERS = (
    "00",
    "01",
    "02",
    "03",
    "04",
    "05",
    "06",
    "07",
    "08",
    "09",
    "10",
    "10A",
    "11",
    "12",
    "12A",
    "13",
    "14",
    "15",
    "16",
    "16A",
    "16B",
    "16C",
    "17",
    "18",
    "19",
)
STAGE_PRODUCER_STATE = "FULL_STAGE_PRODUCERS_FROZEN"


class ScientificStageProducerError(ValueError):
    """Preuve d'étape absente, dupliquée ou mal formée."""

    def __init__(self, code: str, detail: str | None = None) -> None:
        self.code = code
        self.detail = detail
        super().__init__(code if detail is None else f"{code}: {detail}")


def _is_sha256(value: Any) -> bool:
    return (
        isinstance(value, str)
        and len(value) == 64
        and all(character in "0123456789abcdef" for character in value)
    )


def _proportion(successes: int, trials: int) -> float:
    if trials < 1:
        raise ScientificStageProducerError("stage00_denominator_empty")
    return successes / trials


def _nonempty_text(value: Any) -> bool:
    return isinstance(value, str) and bool(value)


def _uses_forbidden_inference_source(value: Any) -> bool:
    if not _nonempty_text(value):
        return True
    normalized = value.strip().upper()
    return any(token in normalized for token in ("CLINICAL", "GROUP", "PHENOTYPE"))


def _challenge_concordance(
    records: Any,
    *,
    expected_field: str,
    observed_field: str,
    error_code: str,
) -> tuple[int, int]:
    if not isinstance(records, list) or not records:
        raise ScientificStageProducerError(error_code)
    matches = 0
    for record in records:
        if (
            not isinstance(record, Mapping)
            or not isinstance(record.get(expected_field), bool)
            or not isinstance(record.get(observed_field), bool)
        ):
            raise ScientificStageProducerError(error_code)
        matches += int(record[expected_field] == record[observed_field])
    return matches, len(records)


def _independent_records(
    records: Any,
    *,
    id_field: str,
    error_code: str,
) -> list[Mapping[str, Any]]:
    if not isinstance(records, list) or not records:
        raise ScientificStageProducerError(error_code)
    identifiers: set[str] = set()
    checked: list[Mapping[str, Any]] = []
    for record in records:
        if not isinstance(record, Mapping):
            raise ScientificStageProducerError(error_code)
        identifier = record.get(id_field)
        if not _nonempty_text(identifier) or identifier in identifiers:
            raise ScientificStageProducerError(error_code)
        identifiers.add(identifier)
        checked.append(record)
    return checked


def _flag_successes(
    records: Any,
    *,
    id_field: str,
    flag_field: str,
    error_code: str,
    success_value: bool = True,
) -> tuple[int, int]:
    checked = _independent_records(records, id_field=id_field, error_code=error_code)
    successes = 0
    for record in checked:
        flag = record.get(flag_field)
        if not isinstance(flag, bool):
            raise ScientificStageProducerError(error_code)
        successes += int(flag is success_value)
    return successes, len(checked)


def _rate(value: Any, *, error_code: str) -> float:
    if (
        isinstance(value, bool)
        or not isinstance(value, (int, float))
        or not 0.0 <= float(value) <= 1.0
    ):
        raise ScientificStageProducerError(error_code)
    return float(value)


def _finite_number(value: Any, *, error_code: str) -> float:
    if (
        isinstance(value, bool)
        or not isinstance(value, (int, float))
        or not math.isfinite(float(value))
    ):
        raise ScientificStageProducerError(error_code)
    return float(value)


def produce_stage_00_observations(evidence: Mapping[str, Any]) -> list[dict[str, Any]]:
    """Produit les quatre contrôles déterministes de l'étape 00.

    ``governance_roots`` contient exactement les huit racines attendues avec
    leur empreinte gelée et observée. ``isolation_probes`` contient les essais
    d'altération, de lecture privée et d'écriture après clôture. ``restart``
    compare les empreintes scientifiques du run initial et de sa reprise.
    """
    roots = evidence.get("governance_roots")
    if not isinstance(roots, list) or len(roots) != 8:
        raise ScientificStageProducerError("stage00_root_count_mismatch")
    root_ids: set[str] = set()
    complete_count = 0
    matching_hash_count = 0
    for root in roots:
        if not isinstance(root, Mapping):
            raise ScientificStageProducerError("stage00_root_invalid")
        root_id = root.get("root_id")
        if not isinstance(root_id, str) or not root_id or root_id in root_ids:
            raise ScientificStageProducerError("stage00_root_id_invalid")
        root_ids.add(root_id)
        present = root.get("present")
        if not isinstance(present, bool):
            raise ScientificStageProducerError("stage00_root_presence_invalid", root_id)
        expected_sha256 = root.get("expected_sha256")
        observed_sha256 = root.get("observed_sha256")
        if not _is_sha256(expected_sha256):
            raise ScientificStageProducerError("stage00_expected_hash_invalid", root_id)
        if observed_sha256 is not None and not _is_sha256(observed_sha256):
            raise ScientificStageProducerError("stage00_observed_hash_invalid", root_id)
        complete_count += int(present)
        matching_hash_count += int(present and observed_sha256 == expected_sha256)

    probes = evidence.get("isolation_probes")
    if not isinstance(probes, list) or not probes:
        raise ScientificStageProducerError("stage00_isolation_probes_missing")
    probe_ids: set[str] = set()
    denied_count = 0
    for probe in probes:
        if not isinstance(probe, Mapping):
            raise ScientificStageProducerError("stage00_isolation_probe_invalid")
        probe_id = probe.get("probe_id")
        if not isinstance(probe_id, str) or not probe_id or probe_id in probe_ids:
            raise ScientificStageProducerError("stage00_isolation_probe_id_invalid")
        probe_ids.add(probe_id)
        denied = probe.get("access_or_write_denied")
        if not isinstance(denied, bool):
            raise ScientificStageProducerError(
                "stage00_isolation_probe_result_invalid", probe_id
            )
        denied_count += int(denied)

    restart = evidence.get("restart")
    if not isinstance(restart, Mapping):
        raise ScientificStageProducerError("stage00_restart_evidence_missing")
    initial_hashes = restart.get("initial_scientific_sha256")
    resumed_hashes = restart.get("resumed_scientific_sha256")
    if (
        not isinstance(initial_hashes, Mapping)
        or not initial_hashes
        or not isinstance(resumed_hashes, Mapping)
        or set(initial_hashes) != set(resumed_hashes)
    ):
        raise ScientificStageProducerError("stage00_restart_hash_set_mismatch")
    for artifact_id, digest in initial_hashes.items():
        if (
            not isinstance(artifact_id, str)
            or not artifact_id
            or not _is_sha256(digest)
            or not _is_sha256(resumed_hashes[artifact_id])
        ):
            raise ScientificStageProducerError("stage00_restart_hash_invalid")
    restart_matches = initial_hashes == resumed_hashes

    return [
        exact_observation(
            metric_id="root_completeness",
            estimate=_proportion(complete_count, len(roots)),
            independent_n=len(roots),
            independent_unit="governance_root",
        ),
        exact_observation(
            metric_id="hash_integrity",
            estimate=_proportion(matching_hash_count, len(roots)),
            independent_n=len(roots),
            independent_unit="declared_artifact",
        ),
        exact_observation(
            metric_id="isolation",
            estimate=_proportion(denied_count, len(probes)),
            independent_n=len(probes),
            independent_unit="isolation_probe",
        ),
        exact_observation(
            metric_id="restart_identity",
            estimate=float(restart_matches),
            independent_n=1,
            independent_unit="paired_run",
        ),
    ]


def produce_stage_01_observations(evidence: Mapping[str, Any]) -> list[dict[str, Any]]:
    """Produit classification, couverture et provenance des sources."""
    sentinels = evidence.get("source_sentinels")
    if not isinstance(sentinels, list) or not sentinels:
        raise ScientificStageProducerError("stage01_source_sentinels_missing")
    positive_count = negative_count = true_positive = true_negative = 0
    for sentinel in sentinels:
        if (
            not isinstance(sentinel, Mapping)
            or not isinstance(sentinel.get("expected_valid"), bool)
            or not isinstance(sentinel.get("observed_accepted"), bool)
        ):
            raise ScientificStageProducerError("stage01_source_sentinel_invalid")
        if sentinel["expected_valid"]:
            positive_count += 1
            true_positive += int(sentinel["observed_accepted"])
        else:
            negative_count += 1
            true_negative += int(not sentinel["observed_accepted"])
    if positive_count == 0 or negative_count == 0:
        raise ScientificStageProducerError("stage01_source_sentinel_class_missing")

    autosomes = evidence.get("covered_autosomes")
    if (
        not isinstance(autosomes, list)
        or any(isinstance(value, bool) or not isinstance(value, int) for value in autosomes)
        or len(autosomes) != len(set(autosomes))
    ):
        raise ScientificStageProducerError("stage01_autosome_coverage_invalid")
    covered_autosome_count = len(set(autosomes).intersection(range(1, 23)))

    expected_target = evidence.get("expected_target")
    observed_targets = evidence.get("observed_targets")
    target_fields = ("chromosome", "position_bp", "ref", "alt", "variant_id")
    if (
        not isinstance(expected_target, Mapping)
        or not isinstance(observed_targets, list)
        or any(field not in expected_target for field in target_fields)
        or any(not isinstance(target, Mapping) for target in observed_targets)
    ):
        raise ScientificStageProducerError("stage01_target_evidence_invalid")
    target_matches = sum(
        all(target.get(field) == expected_target[field] for field in target_fields)
        for target in observed_targets
    )

    provenance = evidence.get("source_provenance")
    required_provenance = ("assembly", "license", "version", "sha256")
    if not isinstance(provenance, list) or not provenance:
        raise ScientificStageProducerError("stage01_provenance_missing")
    provenance_complete = 0
    for record in provenance:
        if not isinstance(record, Mapping):
            raise ScientificStageProducerError("stage01_provenance_invalid")
        complete = all(_nonempty_text(record.get(field)) for field in required_provenance)
        if _nonempty_text(record.get("sha256")) and not _is_sha256(record["sha256"]):
            raise ScientificStageProducerError("stage01_provenance_hash_invalid")
        provenance_complete += int(complete)

    return [
        exact_observation(
            metric_id="source_sensitivity",
            estimate=_proportion(true_positive, positive_count),
            independent_n=positive_count,
            independent_unit="deterministic_source_sentinel",
        ),
        exact_observation(
            metric_id="source_specificity",
            estimate=_proportion(true_negative, negative_count),
            independent_n=negative_count,
            independent_unit="deterministic_source_sentinel",
        ),
        exact_observation(
            metric_id="autosome_coverage",
            estimate=covered_autosome_count,
            independent_n=1,
            independent_unit="release",
        ),
        exact_observation(
            metric_id="target_coverage",
            estimate=float(target_matches == 1),
            independent_n=1,
            independent_unit="configured_target",
        ),
        exact_observation(
            metric_id="provenance",
            estimate=_proportion(provenance_complete, len(provenance)),
            independent_n=len(provenance),
            independent_unit="accepted_source",
        ),
    ]


def produce_stage_02_observations(evidence: Mapping[str, Any]) -> list[dict[str, Any]]:
    """Produit les contrôles de bijection, pedigree et non-inférence clinique."""
    rows = evidence.get("registry_rows")
    id_fields = ("source_sample_id", "sample_id", "identity_id")
    if not isinstance(rows, list) or not rows:
        raise ScientificStageProducerError("stage02_registry_missing")
    valid_rows = 0
    columns: dict[str, set[str]] = {field: set() for field in id_fields}
    for row in rows:
        if not isinstance(row, Mapping):
            raise ScientificStageProducerError("stage02_registry_row_invalid")
        values = [row.get(field) for field in id_fields]
        if all(_nonempty_text(value) for value in values):
            valid_rows += 1
            for field, value in zip(id_fields, values, strict=True):
                columns[field].add(value)
    bijective = valid_rows == len(rows) and all(
        len(values) == len(rows) for values in columns.values()
    )

    collision_matches, collision_count = _challenge_concordance(
        evidence.get("identity_challenges"),
        expected_field="expected_detected",
        observed_field="observed_detected",
        error_code="stage02_identity_challenges_invalid",
    )
    pedigree_matches, pedigree_count = _challenge_concordance(
        evidence.get("pedigree_challenges"),
        expected_field="expected_flagged",
        observed_field="observed_flagged",
        error_code="stage02_pedigree_challenges_invalid",
    )

    assignments = evidence.get("genotype_assignments")
    if not isinstance(assignments, list) or not assignments:
        raise ScientificStageProducerError("stage02_genotype_assignments_missing")
    compliant_assignments = 0
    for assignment in assignments:
        if not isinstance(assignment, Mapping) or not _nonempty_text(
            assignment.get("genotype_source")
        ):
            raise ScientificStageProducerError("stage02_genotype_assignment_invalid")
        compliant_assignments += int(
            not _uses_forbidden_inference_source(assignment["genotype_source"])
        )

    return [
        exact_observation(
            metric_id="registry_bijection",
            estimate=float(bijective),
            independent_n=len(rows),
            independent_unit="study_identity",
        ),
        exact_observation(
            metric_id="collision_detection",
            estimate=_proportion(collision_matches, collision_count),
            independent_n=collision_count,
            independent_unit="deterministic_identity_challenge",
        ),
        exact_observation(
            metric_id="pedigree_detection",
            estimate=_proportion(pedigree_matches, pedigree_count),
            independent_n=pedigree_count,
            independent_unit="pedigree_component",
        ),
        exact_observation(
            metric_id="no_clinical_inference",
            estimate=_proportion(compliant_assignments, len(assignments)),
            independent_n=len(assignments),
            independent_unit="release",
        ),
    ]


def produce_stage_03_observations(evidence: Mapping[str, Any]) -> list[dict[str, Any]]:
    """Produit les fidélités de conversion à partir d'indices et empreintes clos."""
    expected_samples = evidence.get("expected_sample_ids")
    observed_samples = evidence.get("observed_sample_ids")
    if (
        not isinstance(expected_samples, list)
        or not expected_samples
        or any(not _nonempty_text(value) for value in expected_samples)
        or not isinstance(observed_samples, list)
        or any(not _nonempty_text(value) for value in observed_samples)
    ):
        raise ScientificStageProducerError("stage03_sample_index_invalid")
    sample_matches = sum(
        expected == observed
        for expected, observed in zip(expected_samples, observed_samples, strict=False)
    )
    sample_denominator = max(len(expected_samples), len(observed_samples))

    variant_index = evidence.get("variant_index")
    if not isinstance(variant_index, Mapping):
        raise ScientificStageProducerError("stage03_variant_index_missing")
    expected_variant_count = variant_index.get("expected_count")
    observed_variant_count = variant_index.get("observed_count")
    if (
        isinstance(expected_variant_count, bool)
        or not isinstance(expected_variant_count, int)
        or expected_variant_count < 1
        or isinstance(observed_variant_count, bool)
        or not isinstance(observed_variant_count, int)
        or not _is_sha256(variant_index.get("expected_sha256"))
        or not _is_sha256(variant_index.get("observed_sha256"))
    ):
        raise ScientificStageProducerError("stage03_variant_index_invalid")
    variants_match = (
        expected_variant_count == observed_variant_count
        and variant_index["expected_sha256"] == variant_index["observed_sha256"]
    )

    genotype_worlds = evidence.get("genotype_worlds")
    if not isinstance(genotype_worlds, list) or not genotype_worlds:
        raise ScientificStageProducerError("stage03_genotype_worlds_missing")
    matching_worlds = 0
    for world in genotype_worlds:
        if (
            not isinstance(world, Mapping)
            or not _is_sha256(world.get("expected_sha256"))
            or not _is_sha256(world.get("observed_sha256"))
        ):
            raise ScientificStageProducerError("stage03_genotype_world_invalid")
        matching_worlds += int(world["expected_sha256"] == world["observed_sha256"])

    missingness = evidence.get("missingness")
    if (
        not isinstance(missingness, Mapping)
        or isinstance(missingness.get("expected_count"), bool)
        or not isinstance(missingness.get("expected_count"), int)
        or isinstance(missingness.get("observed_count"), bool)
        or not isinstance(missingness.get("observed_count"), int)
    ):
        raise ScientificStageProducerError("stage03_missingness_invalid")
    missingness_difference = abs(
        missingness["observed_count"] - missingness["expected_count"]
    )

    allele_sentinels = evidence.get("allele_sentinels")
    if not isinstance(allele_sentinels, list) or not allele_sentinels:
        raise ScientificStageProducerError("stage03_allele_sentinels_missing")
    allele_matches = 0
    for sentinel in allele_sentinels:
        if not isinstance(sentinel, Mapping):
            raise ScientificStageProducerError("stage03_allele_sentinel_invalid")
        fields = ("expected_ref", "expected_alt", "observed_ref", "observed_alt")
        if any(not _nonempty_text(sentinel.get(field)) for field in fields):
            raise ScientificStageProducerError("stage03_allele_sentinel_invalid")
        allele_matches += int(
            sentinel["expected_ref"] == sentinel["observed_ref"]
            and sentinel["expected_alt"] == sentinel["observed_alt"]
        )

    return [
        exact_observation(
            metric_id="sample_concordance",
            estimate=_proportion(sample_matches, sample_denominator),
            independent_n=sample_denominator,
            independent_unit="study_identity",
        ),
        exact_observation(
            metric_id="variant_concordance",
            estimate=float(variants_match),
            independent_n=expected_variant_count,
            independent_unit="observed_variant",
        ),
        exact_observation(
            metric_id="genotype_concordance",
            estimate=_proportion(matching_worlds, len(genotype_worlds)),
            independent_n=len(genotype_worlds),
            independent_unit="scenario_world",
        ),
        exact_observation(
            metric_id="missingness_count",
            estimate=missingness_difference,
            independent_n=1,
            independent_unit="scenario_world",
        ),
        exact_observation(
            metric_id="allele_integrity",
            estimate=_proportion(allele_matches, len(allele_sentinels)),
            independent_n=len(allele_sentinels),
            independent_unit="allele_sentinel",
        ),
    ]


def produce_stage_04_observations(evidence: Mapping[str, Any]) -> list[dict[str, Any]]:
    """Produit concordance, manque, provenance moléculaire et Mendel."""
    calls = evidence.get("target_calls")
    if not isinstance(calls, list) or not calls:
        raise ScientificStageProducerError("stage04_target_calls_missing")
    nonmissing_count = concordant_count = missingness_preserved = provenance_count = 0
    for call in calls:
        if (
            not isinstance(call, Mapping)
            or not isinstance(call.get("expected_missing"), bool)
            or not isinstance(call.get("observed_missing"), bool)
        ):
            raise ScientificStageProducerError("stage04_target_call_invalid")
        missingness_preserved += int(
            call["expected_missing"] == call["observed_missing"]
        )
        if not call["expected_missing"]:
            nonmissing_count += 1
            concordant_count += int(call.get("expected_call") == call.get("observed_call"))
        provenance_count += int(
            _nonempty_text(call.get("molecular_source"))
            and not _uses_forbidden_inference_source(call.get("inference_source"))
        )
    if nonmissing_count == 0:
        raise ScientificStageProducerError("stage04_nonmissing_calls_absent")

    mendelian_matches, trio_count = _challenge_concordance(
        evidence.get("mendelian_challenges"),
        expected_field="expected_flagged",
        observed_field="observed_flagged",
        error_code="stage04_mendelian_challenges_invalid",
    )
    return [
        exact_observation(
            metric_id="target_call_concordance",
            estimate=_proportion(concordant_count, nonmissing_count),
            independent_n=len(calls),
            independent_unit="study_identity",
        ),
        exact_observation(
            metric_id="missingness_preservation",
            estimate=_proportion(missingness_preserved, len(calls)),
            independent_n=len(calls),
            independent_unit="study_identity",
        ),
        exact_observation(
            metric_id="molecular_provenance",
            estimate=_proportion(provenance_count, len(calls)),
            independent_n=len(calls),
            independent_unit="target_call",
        ),
        exact_observation(
            metric_id="mendelian_detection",
            estimate=_proportion(mendelian_matches, trio_count),
            independent_n=trio_count,
            independent_unit="trio",
        ),
    ]


def produce_stage_05_observations(
    evidence: Mapping[str, Any],
    *,
    evaluation_secret: bytes,
    release_id: str,
    bootstrap_replicates: int = DEFAULT_BOOTSTRAP_REPLICATES,
) -> list[dict[str, Any]]:
    """Produit les performances du QC préliminaire et leurs incertitudes."""
    high_successes, high_count = _flag_successes(
        evidence.get("high_defect_units"),
        id_field="unit_id",
        flag_field="observed_flagged",
        error_code="stage05_high_defect_units_invalid",
    )
    low_successes, low_count = _flag_successes(
        evidence.get("low_defect_units"),
        id_field="unit_id",
        flag_field="observed_flagged",
        error_code="stage05_low_defect_units_invalid",
    )
    negative_successes, negative_count = _flag_successes(
        evidence.get("negative_units"),
        id_field="unit_id",
        flag_field="observed_flagged",
        success_value=False,
        error_code="stage05_negative_units_invalid",
    )

    scenarios = _independent_records(
        evidence.get("scenario_units"),
        id_field="scenario_unit_id",
        error_code="stage05_scenario_units_invalid",
    )
    missingness_errors: list[float] = []
    call_rate_errors: list[float] = []
    cluster_ids: list[str] = []
    evaluable_count = 0
    for scenario in scenarios:
        expected_missingness = _rate(
            scenario.get("expected_missingness"),
            error_code="stage05_scenario_rate_invalid",
        )
        observed_missingness = _rate(
            scenario.get("observed_missingness"),
            error_code="stage05_scenario_rate_invalid",
        )
        expected_call_rate = _rate(
            scenario.get("expected_call_rate"),
            error_code="stage05_scenario_rate_invalid",
        )
        observed_call_rate = _rate(
            scenario.get("observed_call_rate"),
            error_code="stage05_scenario_rate_invalid",
        )
        if not isinstance(scenario.get("evaluable"), bool):
            raise ScientificStageProducerError("stage05_scenario_evaluability_invalid")
        missingness_errors.append(observed_missingness - expected_missingness)
        call_rate_errors.append(observed_call_rate - expected_call_rate)
        cluster_ids.append(scenario["scenario_unit_id"])
        evaluable_count += int(scenario["evaluable"])

    bootstrap_arguments = {
        "cluster_ids": cluster_ids,
        "independent_unit": "scenario_unit",
        "statistic_id": "MEAN_ABSOLUTE",
        "evaluation_secret": evaluation_secret,
        "release_id": release_id,
        "stage_id": "05",
        "replicates": bootstrap_replicates,
    }
    return [
        proportion_observation(
            metric_id="high_defect_sensitivity",
            successes=high_successes,
            trials=high_count,
            independent_unit="injected_high_defect_unit",
            method="WILSON_ONE_SIDED",
        ),
        proportion_observation(
            metric_id="low_defect_sensitivity",
            successes=low_successes,
            trials=low_count,
            independent_unit="injected_low_defect_unit",
            method="WILSON_ONE_SIDED",
        ),
        proportion_observation(
            metric_id="specificity",
            successes=negative_successes,
            trials=negative_count,
            independent_unit="negative_unit",
            method="WILSON_ONE_SIDED",
        ),
        cluster_statistic_observation(
            metric_id="missingness_error",
            values=missingness_errors,
            **bootstrap_arguments,
        ),
        cluster_statistic_observation(
            metric_id="call_rate_error",
            values=call_rate_errors,
            **bootstrap_arguments,
        ),
        proportion_observation(
            metric_id="evaluability",
            successes=evaluable_count,
            trials=len(scenarios),
            independent_unit="eligible_scenario_unit",
            method="WILSON_ONE_SIDED",
        ),
    ]


def produce_stage_06_observations(evidence: Mapping[str, Any]) -> list[dict[str, Any]]:
    """Produit fidélité du panel de parenté, densités et LD résiduel."""
    variants = _independent_records(
        evidence.get("variant_units"),
        id_field="variant_id",
        error_code="stage06_variant_units_invalid",
    )
    eligible_count = selected_count = selected_eligible_count = 0
    for variant in variants:
        eligible = variant.get("expected_eligible")
        selected = variant.get("observed_selected")
        if not isinstance(eligible, bool) or not isinstance(selected, bool):
            raise ScientificStageProducerError("stage06_variant_classification_invalid")
        eligible_count += int(eligible)
        selected_count += int(selected)
        selected_eligible_count += int(eligible and selected)
    if eligible_count == 0 or selected_count == 0:
        raise ScientificStageProducerError("stage06_variant_class_empty")

    chromosome_counts = evidence.get("chromosome_variant_counts")
    density_contexts = {
        "D_MAX": "density_dmax",
        "D_MID": "density_dmid",
        "D_LOW": "density_dlow",
    }
    if not isinstance(chromosome_counts, Mapping):
        raise ScientificStageProducerError("stage06_chromosome_counts_missing")
    density_observations: list[dict[str, Any]] = []
    for density, context_id in density_contexts.items():
        counts = chromosome_counts.get(density)
        if (
            not isinstance(counts, Mapping)
            or set(counts) != set(range(1, 23))
            or any(
                isinstance(count, bool) or not isinstance(count, int) or count < 0
                for count in counts.values()
            )
        ):
            raise ScientificStageProducerError("stage06_chromosome_counts_invalid", density)
        density_observations.append(
            exact_observation(
                metric_id="density_minimum",
                estimate=min(counts.values()),
                independent_n=22,
                independent_unit="autosome",
                context_id=context_id,
                context={"density": density},
            )
        )

    ld_blocks = _independent_records(
        evidence.get("ld_pair_blocks"),
        id_field="block_id",
        error_code="stage06_ld_pair_blocks_invalid",
    )
    exceeding_count = 0
    for block in ld_blocks:
        exceeds = block.get("exceeds_pruning_threshold")
        if not isinstance(exceeds, bool):
            raise ScientificStageProducerError("stage06_ld_pair_block_invalid")
        exceeding_count += int(exceeds)

    violations = evidence.get("forbidden_variant_uses")
    if not isinstance(violations, list) or any(
        not isinstance(violation, Mapping) for violation in violations
    ):
        raise ScientificStageProducerError("stage06_forbidden_variant_uses_invalid")

    return [
        proportion_observation(
            metric_id="eligibility_recall",
            successes=selected_eligible_count,
            trials=eligible_count,
            independent_unit="eligible_variant",
            method="WILSON_ONE_SIDED",
        ),
        proportion_observation(
            metric_id="eligibility_precision",
            successes=selected_eligible_count,
            trials=selected_count,
            independent_unit="selected_variant",
            method="WILSON_ONE_SIDED",
        ),
        exact_observation(
            metric_id="chromosome_coverage",
            estimate=len(set(chromosome_counts["D_MAX"])),
            independent_n=1,
            independent_unit="release",
        ),
        *density_observations,
        proportion_observation(
            metric_id="residual_ld",
            successes=exceeding_count,
            trials=len(ld_blocks),
            independent_unit="variant_pair_block",
            method="WILSON_ONE_SIDED",
        ),
        exact_observation(
            metric_id="forbidden_variant_use",
            estimate=len(violations),
            independent_n=1,
            independent_unit="release",
        ),
    ]


def produce_stage_07_observations(
    evidence: Mapping[str, Any],
    *,
    evaluation_secret: bytes,
    release_id: str,
    bootstrap_replicates: int = DEFAULT_BOOTSTRAP_REPLICATES,
) -> list[dict[str, Any]]:
    """Produit classification et erreur de parenté, grappées par composante."""
    pairs = _independent_records(
        evidence.get("designated_pairs"),
        id_field="pair_id",
        error_code="stage07_designated_pairs_invalid",
    )
    related_count = related_correct = observed_related_count = observed_related_correct = 0
    unrelated_count = unrelated_correct = cryptic_count = cryptic_detected = 0
    kinship_errors: list[float] = []
    component_ids: list[str] = []
    for pair in pairs:
        expected_class = pair.get("expected_class")
        observed_class = pair.get("observed_class")
        component_id = pair.get("component_id")
        cryptic = pair.get("cryptic")
        if (
            not _nonempty_text(expected_class)
            or not _nonempty_text(observed_class)
            or not _nonempty_text(component_id)
            or not isinstance(cryptic, bool)
        ):
            raise ScientificStageProducerError("stage07_pair_classification_invalid")
        expected_related = expected_class != "UNRELATED"
        observed_related = observed_class != "UNRELATED"
        exact_class = expected_class == observed_class
        related_count += int(expected_related)
        related_correct += int(expected_related and exact_class)
        observed_related_count += int(observed_related)
        observed_related_correct += int(observed_related and exact_class)
        unrelated_count += int(not expected_related)
        unrelated_correct += int(not expected_related and not observed_related)
        cryptic_count += int(cryptic)
        cryptic_detected += int(cryptic and observed_related)
        expected_kinship = _finite_number(
            pair.get("expected_kinship"), error_code="stage07_kinship_invalid"
        )
        observed_kinship = _finite_number(
            pair.get("observed_kinship"), error_code="stage07_kinship_invalid"
        )
        kinship_errors.append(observed_kinship - expected_kinship)
        component_ids.append(component_id)
    if min(related_count, observed_related_count, unrelated_count, cryptic_count) == 0:
        raise ScientificStageProducerError("stage07_pair_class_missing")

    return [
        proportion_observation(
            metric_id="related_class_sensitivity",
            successes=related_correct,
            trials=related_count,
            independent_unit="designated_related_pair",
            method="WILSON_ONE_SIDED",
        ),
        proportion_observation(
            metric_id="related_class_ppv",
            successes=observed_related_correct,
            trials=observed_related_count,
            independent_unit="designated_related_pair",
            method="WILSON_ONE_SIDED",
        ),
        proportion_observation(
            metric_id="unrelated_specificity",
            successes=unrelated_correct,
            trials=unrelated_count,
            independent_unit="designated_unrelated_pair",
            method="WILSON_ONE_SIDED",
        ),
        cluster_statistic_observation(
            metric_id="kinship_mae",
            values=kinship_errors,
            cluster_ids=component_ids,
            independent_unit="designated_pair_clustered_by_component",
            statistic_id="MEAN_ABSOLUTE",
            evaluation_secret=evaluation_secret,
            release_id=release_id,
            stage_id="07",
            replicates=bootstrap_replicates,
        ),
        proportion_observation(
            metric_id="cryptic_detection",
            successes=cryptic_detected,
            trials=cryptic_count,
            independent_unit="cryptic_related_pair",
            method="WILSON_ONE_SIDED",
        ),
    ]


def produce_stage_08_observations(
    evidence: Mapping[str, Any],
    *,
    evaluation_secret: bytes,
    release_id: str,
    bootstrap_replicates: int = DEFAULT_BOOTSTRAP_REPLICATES,
) -> list[dict[str, Any]]:
    """Produit séparation PCA, atypiques et stabilité Procrustes."""
    fit_subjects = _independent_records(
        evidence.get("fit_subjects"),
        id_field="subject_id",
        error_code="stage08_fit_subjects_invalid",
    )
    fit_relative_count = 0
    for subject in fit_subjects:
        if not isinstance(subject.get("is_relative"), bool):
            raise ScientificStageProducerError("stage08_fit_subject_invalid")
        fit_relative_count += int(subject["is_relative"])

    representatives = _independent_records(
        evidence.get("representatives"),
        id_field="representative_id",
        error_code="stage08_representatives_invalid",
    )
    pure_count = pure_correct = outlier_count = outlier_detected = 0
    nonoutlier_count = nonoutlier_correct = 0
    for representative in representatives:
        boolean_fields = (
            "pure_state",
            "centroid_correct",
            "expected_outlier",
            "observed_outlier",
        )
        if any(not isinstance(representative.get(field), bool) for field in boolean_fields):
            raise ScientificStageProducerError("stage08_representative_invalid")
        pure_count += int(representative["pure_state"])
        pure_correct += int(
            representative["pure_state"] and representative["centroid_correct"]
        )
        outlier_count += int(representative["expected_outlier"])
        outlier_detected += int(
            representative["expected_outlier"] and representative["observed_outlier"]
        )
        nonoutlier_count += int(not representative["expected_outlier"])
        nonoutlier_correct += int(
            not representative["expected_outlier"]
            and not representative["observed_outlier"]
        )
    if min(pure_count, outlier_count, nonoutlier_count) == 0:
        raise ScientificStageProducerError("stage08_representative_class_missing")

    shifts = _independent_records(
        evidence.get("procrustes_shifts"),
        id_field="representative_id",
        error_code="stage08_procrustes_shifts_invalid",
    )
    shift_values = [
        _finite_number(
            shift.get("shift_standard_deviation"),
            error_code="stage08_procrustes_shift_invalid",
        )
        for shift in shifts
    ]
    if any(value < 0.0 for value in shift_values):
        raise ScientificStageProducerError("stage08_procrustes_shift_invalid")
    shift_ids = [shift["representative_id"] for shift in shifts]
    representative_ids = {
        representative["representative_id"] for representative in representatives
    }
    if not set(shift_ids).issubset(representative_ids):
        raise ScientificStageProducerError("stage08_procrustes_identity_unknown")
    bootstrap_arguments = {
        "values": shift_values,
        "cluster_ids": shift_ids,
        "independent_unit": "independent_representative",
        "evaluation_secret": evaluation_secret,
        "release_id": release_id,
        "stage_id": "08",
        "replicates": bootstrap_replicates,
    }
    return [
        exact_observation(
            metric_id="fit_relative_count",
            estimate=fit_relative_count,
            independent_n=len(fit_subjects),
            independent_unit="fitted_reference",
        ),
        proportion_observation(
            metric_id="pure_centroid_accuracy",
            successes=pure_correct,
            trials=pure_count,
            independent_unit="independent_representative",
            method="WILSON_ONE_SIDED",
        ),
        proportion_observation(
            metric_id="outlier_sensitivity",
            successes=outlier_detected,
            trials=outlier_count,
            independent_unit="independent_representative",
            method="WILSON_ONE_SIDED",
        ),
        proportion_observation(
            metric_id="outlier_specificity",
            successes=nonoutlier_correct,
            trials=nonoutlier_count,
            independent_unit="independent_representative",
            method="WILSON_ONE_SIDED",
        ),
        cluster_statistic_observation(
            metric_id="procrustes_median",
            statistic_id="MEDIAN",
            **bootstrap_arguments,
        ),
        cluster_statistic_observation(
            metric_id="procrustes_p95",
            statistic_id="QUANTILE_95",
            **bootstrap_arguments,
        ),
    ]


def produce_stage_09_observations(evidence: Mapping[str, Any]) -> list[dict[str, Any]]:
    """Produit les contrôles déterministes du gel signé des cohortes."""
    assignments = _independent_records(
        evidence.get("cohort_assignments"),
        id_field="identity_id",
        error_code="stage09_cohort_assignments_invalid",
    )
    assignment_matches = 0
    for assignment in assignments:
        expected_cohort = assignment.get("expected_cohort")
        observed_cohort = assignment.get("observed_cohort")
        if not _nonempty_text(expected_cohort) or not _nonempty_text(observed_cohort):
            raise ScientificStageProducerError("stage09_cohort_assignment_invalid")
        assignment_matches += int(expected_cohort == observed_cohort)

    components = _independent_records(
        evidence.get("family_components"),
        id_field="component_id",
        error_code="stage09_family_components_invalid",
    )
    independent_components = 0
    for component in components:
        preserved = component.get("independent_partition_preserved")
        if not isinstance(preserved, bool):
            raise ScientificStageProducerError("stage09_family_component_invalid")
        independent_components += int(preserved)

    mandatory_cohorts = _independent_records(
        evidence.get("mandatory_cohorts"),
        id_field="cohort_id",
        error_code="stage09_mandatory_cohorts_invalid",
    )
    nonempty_cohorts = 0
    for cohort in mandatory_cohorts:
        member_count = cohort.get("member_count")
        if (
            isinstance(member_count, bool)
            or not isinstance(member_count, int)
            or member_count < 0
        ):
            raise ScientificStageProducerError("stage09_mandatory_cohort_invalid")
        nonempty_cohorts += int(member_count > 0)

    review = evidence.get("review")
    if (
        not isinstance(review, Mapping)
        or not isinstance(review.get("required"), bool)
        or not isinstance(review.get("signed"), bool)
        or not _is_sha256(review.get("expected_sha256"))
    ):
        raise ScientificStageProducerError("stage09_review_invalid")
    observed_review_sha256 = review.get("observed_sha256")
    if observed_review_sha256 is not None and not _is_sha256(observed_review_sha256):
        raise ScientificStageProducerError("stage09_review_hash_invalid")
    signature_valid = (
        not review["required"]
        or review["signed"]
        and observed_review_sha256 == review["expected_sha256"]
    )
    review_observation = exact_observation(
        metric_id="signed_review",
        estimate=float(signature_valid),
        independent_n=1,
        independent_unit="release",
    )
    if review["required"] and not review["signed"]:
        review_observation["structural_not_evaluable_reason"] = "HUMAN_REVIEW_UNSIGNED"

    return [
        exact_observation(
            metric_id="cohort_concordance",
            estimate=_proportion(assignment_matches, len(assignments)),
            independent_n=len(assignments),
            independent_unit="study_identity",
        ),
        exact_observation(
            metric_id="family_independence",
            estimate=_proportion(independent_components, len(components)),
            independent_n=len(components),
            independent_unit="pedigree_component",
        ),
        exact_observation(
            metric_id="mandatory_nonempty",
            estimate=_proportion(nonempty_cohorts, len(mandatory_cohorts)),
            independent_n=len(mandatory_cohorts),
            independent_unit="mandatory_cohort",
        ),
        review_observation,
    ]


def produce_stage_10_observations(evidence: Mapping[str, Any]) -> list[dict[str, Any]]:
    """Produit les décisions QC finales et les audits de cible et cohorte."""
    high_successes, high_count = _flag_successes(
        evidence.get("high_defect_units"),
        id_field="unit_id",
        flag_field="observed_flagged",
        error_code="stage10_high_defect_units_invalid",
    )
    low_successes, low_count = _flag_successes(
        evidence.get("low_defect_units"),
        id_field="unit_id",
        flag_field="observed_flagged",
        error_code="stage10_low_defect_units_invalid",
    )
    negative_successes, negative_count = _flag_successes(
        evidence.get("negative_units"),
        id_field="unit_id",
        flag_field="observed_flagged",
        success_value=False,
        error_code="stage10_negative_units_invalid",
    )
    target_matches, target_count = _challenge_concordance(
        evidence.get("target_exception_checks"),
        expected_field="expected_preserved",
        observed_field="observed_preserved",
        error_code="stage10_target_exception_checks_invalid",
    )
    changes = _independent_records(
        evidence.get("cohort_changes"),
        id_field="change_id",
        error_code="stage10_cohort_changes_invalid",
    )
    audited_count = 0
    for change in changes:
        if not isinstance(change.get("audited"), bool):
            raise ScientificStageProducerError("stage10_cohort_change_invalid")
        audited_count += int(change["audited"])
    return [
        proportion_observation(
            metric_id="high_defect_sensitivity",
            successes=high_successes,
            trials=high_count,
            independent_unit="injected_high_defect_unit",
            method="WILSON_ONE_SIDED",
        ),
        proportion_observation(
            metric_id="low_defect_sensitivity",
            successes=low_successes,
            trials=low_count,
            independent_unit="injected_low_defect_unit",
            method="WILSON_ONE_SIDED",
        ),
        proportion_observation(
            metric_id="specificity",
            successes=negative_successes,
            trials=negative_count,
            independent_unit="negative_unit",
            method="WILSON_ONE_SIDED",
        ),
        exact_observation(
            metric_id="target_exception_integrity",
            estimate=_proportion(target_matches, target_count),
            independent_n=target_count,
            independent_unit="configured_target_exception",
        ),
        exact_observation(
            metric_id="cohort_audit",
            estimate=_proportion(audited_count, len(changes)),
            independent_n=len(changes),
            independent_unit="cohort_change",
        ),
    ]


def produce_stage_10a_observations(evidence: Mapping[str, Any]) -> list[dict[str, Any]]:
    """Produit les contrôles structurels du panel autosomal final."""
    expected_subjects = evidence.get("expected_subject_ids")
    observed_subjects = evidence.get("observed_subject_ids")
    if (
        not isinstance(expected_subjects, list)
        or not expected_subjects
        or any(not _nonempty_text(subject_id) for subject_id in expected_subjects)
        or len(set(expected_subjects)) != len(expected_subjects)
        or not isinstance(observed_subjects, list)
        or any(not _nonempty_text(subject_id) for subject_id in observed_subjects)
        or len(set(observed_subjects)) != len(observed_subjects)
    ):
        raise ScientificStageProducerError("stage10a_subject_index_invalid")
    subject_denominator = max(len(expected_subjects), len(observed_subjects))
    subject_matches = len(set(expected_subjects).intersection(observed_subjects))
    order_matches = sum(
        expected == observed
        for expected, observed in zip(expected_subjects, observed_subjects, strict=False)
    )

    chromosome_counts = evidence.get("autosome_variant_counts")
    if (
        not isinstance(chromosome_counts, Mapping)
        or any(
            isinstance(chromosome, bool)
            or not isinstance(chromosome, int)
            or isinstance(count, bool)
            or not isinstance(count, int)
            or count < 0
            for chromosome, count in chromosome_counts.items()
        )
    ):
        raise ScientificStageProducerError("stage10a_chromosome_counts_invalid")
    structural_missingness = evidence.get("structural_missing_genotype_count")
    duplicate_count = evidence.get("duplicate_variant_count")
    for value in (structural_missingness, duplicate_count):
        if isinstance(value, bool) or not isinstance(value, int) or value < 0:
            raise ScientificStageProducerError("stage10a_count_invalid")

    return [
        exact_observation(
            metric_id="autosome_count",
            estimate=len(chromosome_counts),
            independent_n=1,
            independent_unit="release",
        ),
        exact_observation(
            metric_id="subject_completeness",
            estimate=_proportion(subject_matches, subject_denominator),
            independent_n=subject_denominator,
            independent_unit="study_identity",
        ),
        exact_observation(
            metric_id="order_identity",
            estimate=_proportion(order_matches, subject_denominator),
            independent_n=subject_denominator,
            independent_unit="study_identity",
        ),
        exact_observation(
            metric_id="structural_missingness",
            estimate=structural_missingness,
            independent_n=1,
            independent_unit="release",
        ),
        exact_observation(
            metric_id="duplicates",
            estimate=duplicate_count,
            independent_n=1,
            independent_unit="release",
        ),
        exact_observation(
            metric_id="empty_chromosomes",
            estimate=sum(count == 0 for count in chromosome_counts.values()),
            independent_n=1,
            independent_unit="release",
        ),
    ]


def produce_stage_11_observations(evidence: Mapping[str, Any]) -> list[dict[str, Any]]:
    """Produit exactitude de cible, interpolation et faisabilité des fenêtres."""
    expected_target = evidence.get("expected_target")
    observed_target = evidence.get("observed_target")
    target_fields = ("chromosome", "position_bp", "ref", "alt", "variant_id")
    if not isinstance(expected_target, Mapping) or not isinstance(observed_target, Mapping):
        raise ScientificStageProducerError("stage11_target_invalid")
    if any(field not in expected_target or field not in observed_target for field in target_fields):
        raise ScientificStageProducerError("stage11_target_invalid")
    target_matches = all(
        expected_target[field] == observed_target[field] for field in target_fields
    )

    pairs = _independent_records(
        evidence.get("target_age_pairs"),
        id_field="pair_id",
        error_code="stage11_target_age_pairs_invalid",
    )
    maximum_cm_error = 0.0
    maximum_bp_error = 0.0
    extrapolation_count = 0
    for pair in pairs:
        expected_cm = _finite_number(
            pair.get("expected_cm"), error_code="stage11_cm_value_invalid"
        )
        observed_cm = _finite_number(
            pair.get("observed_cm"), error_code="stage11_cm_value_invalid"
        )
        expected_start = _finite_number(
            pair.get("expected_start_bp"), error_code="stage11_bp_value_invalid"
        )
        observed_start = _finite_number(
            pair.get("observed_start_bp"), error_code="stage11_bp_value_invalid"
        )
        expected_end = _finite_number(
            pair.get("expected_end_bp"), error_code="stage11_bp_value_invalid"
        )
        observed_end = _finite_number(
            pair.get("observed_end_bp"), error_code="stage11_bp_value_invalid"
        )
        if not isinstance(pair.get("extrapolated"), bool):
            raise ScientificStageProducerError("stage11_extrapolation_invalid")
        maximum_cm_error = max(maximum_cm_error, abs(observed_cm - expected_cm))
        maximum_bp_error = max(
            maximum_bp_error,
            abs(observed_start - expected_start),
            abs(observed_end - expected_end),
        )
        extrapolation_count += int(pair["extrapolated"])

    windows = _independent_records(
        evidence.get("evaluable_target_windows"),
        id_field="window_id",
        error_code="stage11_target_windows_invalid",
    )
    compliant_windows = 0
    for window in windows:
        if not isinstance(window.get("density_compliant"), bool):
            raise ScientificStageProducerError("stage11_target_window_invalid")
        compliant_windows += int(window["density_compliant"])

    infeasible_matches, infeasible_count = _challenge_concordance(
        evidence.get("infeasible_target_age_pairs"),
        expected_field="expected_infeasible",
        observed_field="observed_infeasible",
        error_code="stage11_infeasible_pairs_invalid",
    )
    return [
        exact_observation(
            metric_id="target_identity",
            estimate=float(target_matches),
            independent_n=1,
            independent_unit="configured_target",
        ),
        exact_observation(
            metric_id="cm_error",
            estimate=maximum_cm_error,
            independent_n=len(pairs),
            independent_unit="target_age_pair",
        ),
        exact_observation(
            metric_id="bp_error",
            estimate=maximum_bp_error,
            independent_n=len(pairs),
            independent_unit="target_age_pair",
        ),
        exact_observation(
            metric_id="extrapolation",
            estimate=extrapolation_count,
            independent_n=len(pairs),
            independent_unit="target_age_pair",
        ),
        exact_observation(
            metric_id="density",
            estimate=_proportion(compliant_windows, len(windows)),
            independent_n=len(windows),
            independent_unit="evaluable_target_window",
        ),
        exact_observation(
            metric_id="infeasibility_classification",
            estimate=_proportion(infeasible_matches, infeasible_count),
            independent_n=infeasible_count,
            independent_unit="infeasible_target_age_pair",
        ),
    ]


REMAINING_STAGE_IDS = (
    "12",
    "12A",
    "13",
    "14",
    "15",
    "16",
    "16A",
    "16B",
    "16C",
    "17",
    "18",
    "19",
)
REMAINING_CLUSTER_STATISTICS = {
    ("13", "median_jaccard"): "MEDIAN",
    ("13", "normalized_boundary_error"): "MEDIAN_ABSOLUTE",
    ("14", "relative_median_bias"): "MEDIAN_ABSOLUTE",
    ("14", "mdape"): "MEDIAN_ABSOLUTE",
    ("14", "interval_width_ratio"): "MEDIAN",
    ("14", "leave_one_family_out"): "QUANTILE_95",
    ("15", "r2_mae"): "MEAN_ABSOLUTE",
    ("15", "dprime_mae"): "MEAN_ABSOLUTE",
    ("15", "r2_bias"): "ABSOLUTE_MEAN",
    ("15", "dprime_bias"): "ABSOLUTE_MEAN",
    ("15", "density_delta"): "MEDIAN_ABSOLUTE",
    ("16", "boundary_mae"): "MEAN_ABSOLUTE",
    ("16", "length_error"): "MEAN_ABSOLUTE",
    ("16", "burden_error"): "MEAN_ABSOLUTE",
    ("16A", "mixture_rank_correlation"): "MEAN",
    ("16A", "procrustes_median"): "MEDIAN",
    ("16A", "procrustes_p95"): "QUANTILE_95",
    ("16C", "boundary_error"): "MEAN_ABSOLUTE",
    ("16C", "tool_jaccard"): "MEDIAN",
    ("17", "signed_delta"): "MEAN",
    ("17", "monotonicity"): "MEAN",
    ("17", "equivalence"): "MEAN",
    ("17", "tool_seed_convergence"): "STANDARD_DEVIATION",
}


def _cell_key(metric_id: str, context_id: str | None) -> str:
    return metric_id if context_id is None else f"{metric_id}::{context_id}"


def _expanded_rules(stage_rules: Mapping[str, Any]) -> list[tuple[Mapping[str, Any], Mapping[str, Any], str | None]]:
    rules = stage_rules.get("rules")
    if not isinstance(rules, list) or not rules:
        raise ScientificStageProducerError("remaining_stage_rules_invalid")
    expanded: list[tuple[Mapping[str, Any], Mapping[str, Any], str | None]] = []
    for rule in rules:
        if not isinstance(rule, Mapping) or not _nonempty_text(rule.get("metric_id")):
            raise ScientificStageProducerError("remaining_stage_rule_invalid")
        if rule.get("operator") != "contextual":
            expanded.append((rule, rule, None))
            continue
        contexts = rule.get("contexts")
        if not isinstance(contexts, list) or not contexts:
            raise ScientificStageProducerError("remaining_stage_contexts_invalid")
        for context_rule in contexts:
            if not isinstance(context_rule, Mapping) or not _nonempty_text(
                context_rule.get("context_id")
            ):
                raise ScientificStageProducerError("remaining_stage_context_invalid")
            expanded.append((rule, context_rule, context_rule["context_id"]))
    return expanded


def produce_remaining_stage_observations(
    stage_id: str,
    stage_rules: Mapping[str, Any],
    evidence: Mapping[str, Any],
    *,
    evaluation_secret: bytes,
    release_id: str,
    bootstrap_replicates: int = DEFAULT_BOOTSTRAP_REPLICATES,
) -> list[dict[str, Any]]:
    """Produit les métriques `12`–`19` depuis des cellules closes non arrondies.

    Une cellule contient soit les résultats binaires de chaque unité
    indépendante, soit les valeurs par grappe, soit les contrôles déterministes
    élémentaires. La méthode et la statistique sont figées par le couple
    étape–métrique et ne peuvent pas être choisies après lecture du résultat.
    """
    if stage_id not in REMAINING_STAGE_IDS:
        raise ScientificStageProducerError("remaining_stage_id_invalid", stage_id)
    cells = evidence.get("metric_cells")
    if not isinstance(cells, Mapping):
        raise ScientificStageProducerError("remaining_stage_cells_missing", stage_id)
    expanded = _expanded_rules(stage_rules)
    expected_keys = {
        _cell_key(primary["metric_id"], context_id)
        for primary, _, context_id in expanded
    }
    if set(cells) != expected_keys:
        raise ScientificStageProducerError("remaining_stage_cell_coverage_mismatch", stage_id)

    observations: list[dict[str, Any]] = []
    for primary, selected, context_id in expanded:
        metric_id = primary["metric_id"]
        cell = cells[_cell_key(metric_id, context_id)]
        if not isinstance(cell, Mapping):
            raise ScientificStageProducerError("remaining_stage_cell_invalid", metric_id)
        context_fields: dict[str, Any] = {}
        if context_id is not None:
            context = cell.get("context")
            if not isinstance(context, Mapping) or not context:
                raise ScientificStageProducerError("remaining_stage_context_evidence_invalid", metric_id)
            context_fields = {"context_id": context_id, "context": context}

        method = inference_method_for_rule(stage_id, primary)
        if method == "DETERMINISTIC_EXACT":
            threshold = selected.get("value")
            audit_units = cell.get("audit_units")
            if not isinstance(audit_units, list) or not audit_units:
                raise ScientificStageProducerError("remaining_exact_audit_units_invalid", metric_id)
            if threshold == 0:
                violations = cell.get("violations")
                if not isinstance(violations, list):
                    raise ScientificStageProducerError("remaining_exact_violations_invalid", metric_id)
                estimate = float(len(violations))
            elif threshold == 1 or threshold == 1.0:
                checks = cell.get("checks")
                if (
                    not isinstance(checks, list)
                    or not checks
                    or any(not isinstance(check, bool) for check in checks)
                ):
                    raise ScientificStageProducerError("remaining_exact_checks_invalid", metric_id)
                estimate = _proportion(sum(checks), len(checks))
                audit_units = checks
            else:
                items = cell.get("items")
                if not isinstance(items, list):
                    raise ScientificStageProducerError("remaining_exact_items_invalid", metric_id)
                estimate = float(len(items))
            observations.append(
                exact_observation(
                    metric_id=metric_id,
                    estimate=estimate,
                    independent_n=len(audit_units),
                    independent_unit=primary["independent_unit"],
                    **context_fields,
                )
            )
            continue

        if method in {"WILSON_ONE_SIDED", "CLOPPER_PEARSON_ONE_SIDED"}:
            outcomes = cell.get("outcomes")
            if (
                not isinstance(outcomes, list)
                or not outcomes
                or any(not isinstance(outcome, bool) for outcome in outcomes)
            ):
                raise ScientificStageProducerError("remaining_proportion_outcomes_invalid", metric_id)
            observations.append(
                proportion_observation(
                    metric_id=metric_id,
                    successes=sum(outcomes),
                    trials=len(outcomes),
                    independent_unit=primary["independent_unit"],
                    method=method,
                    **context_fields,
                )
            )
            continue

        statistic_id = REMAINING_CLUSTER_STATISTICS.get((stage_id, metric_id))
        values = cell.get("values")
        cluster_ids = cell.get("cluster_ids")
        if statistic_id is None:
            raise ScientificStageProducerError("remaining_cluster_statistic_unbound", f"{stage_id}:{metric_id}")
        if not isinstance(values, list) or not isinstance(cluster_ids, list):
            raise ScientificStageProducerError("remaining_cluster_values_invalid", metric_id)
        observations.append(
            cluster_statistic_observation(
                metric_id=metric_id,
                values=values,
                cluster_ids=cluster_ids,
                independent_unit=primary["independent_unit"],
                statistic_id=statistic_id,
                evaluation_secret=evaluation_secret,
                release_id=release_id,
                stage_id=stage_id,
                replicates=bootstrap_replicates,
                **context_fields,
            )
        )
    return observations


def stage_producer_audit() -> dict[str, Any]:
    """Expose l'avancement sans assimiler un producteur partiel à G5."""
    missing = [stage_id for stage_id in STAGE_IDS if stage_id not in IMPLEMENTED_STAGE_PRODUCERS]
    return {
        "state": STAGE_PRODUCER_STATE,
        "implemented_stage_count": len(IMPLEMENTED_STAGE_PRODUCERS),
        "implemented_stages": list(IMPLEMENTED_STAGE_PRODUCERS),
        "missing_stage_count": len(missing),
        "missing_stages": missing,
        "stage_producers_complete": not missing,
        "scientific_evaluation_authorized": False,
    }
