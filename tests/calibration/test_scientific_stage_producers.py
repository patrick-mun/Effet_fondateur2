from __future__ import annotations

import copy
from pathlib import Path

import pytest

from effet_fondateur.calibration.r13_scientific_evaluator import (
    evaluate_stage_observations,
    load_metric_rule_bundle,
)
from effet_fondateur.calibration.scientific_stage_producers import (
    ScientificStageProducerError,
    produce_stage_00_observations,
    produce_stage_01_observations,
    produce_stage_02_observations,
    produce_stage_03_observations,
    produce_stage_04_observations,
    produce_stage_05_observations,
    produce_stage_06_observations,
    produce_stage_07_observations,
    produce_stage_08_observations,
    produce_stage_09_observations,
    produce_stage_10_observations,
    produce_stage_10a_observations,
    produce_stage_11_observations,
    produce_remaining_stage_observations,
    REMAINING_STAGE_IDS,
    stage_producer_audit,
)


REPOSITORY_ROOT = Path(__file__).resolve().parents[2]


def _context_witness(condition: object) -> dict[str, object]:
    assert isinstance(condition, dict)
    if "any_of" in condition:
        return _context_witness(condition["any_of"][0])
    if "all_of" in condition:
        merged: dict[str, object] = {}
        for branch in condition["all_of"]:
            merged.update(_context_witness(branch))
        return merged
    witness: dict[str, object] = {}
    for key, value in condition.items():
        target_key = key.removesuffix("_min").removesuffix("_max")
        witness[target_key] = value[0] if isinstance(value, list) else value
    return witness


def _passing_remaining_evidence(stage_rules: dict[str, object]) -> dict[str, object]:
    cells: dict[str, object] = {}
    for primary in stage_rules["rules"]:  # type: ignore[index]
        contexts = primary.get("contexts") if primary["operator"] == "contextual" else [None]
        for context_rule in contexts:
            selected = primary if context_rule is None else context_rule
            context_id = None if context_rule is None else context_rule["context_id"]
            key = primary["metric_id"] if context_id is None else f"{primary['metric_id']}::{context_id}"
            minimum_n = primary["minimum_independent_n"]
            cell: dict[str, object] = {}
            if context_rule is not None:
                cell["context"] = _context_witness(context_rule["when"])
            operator = selected["operator"]
            unit = primary["unit"]
            deterministic = primary.get("deterministic") is True or operator == "exact"
            if deterministic:
                cell["audit_units"] = [f"AUDIT-{index}" for index in range(minimum_n)]
                threshold = selected["value"]
                if threshold == 0:
                    cell["violations"] = []
                elif threshold == 1 or threshold == 1.0:
                    cell["checks"] = [True] * minimum_n
                else:
                    cell["items"] = [f"ITEM-{index}" for index in range(int(threshold))]
            elif unit == "proportion":
                precision_n = max(minimum_n, 500)
                cell["outcomes"] = [operator == "min"] * precision_n
            else:
                value = 1.0 if operator == "min" else 0.0
                cell["values"] = [value] * minimum_n
                cell["cluster_ids"] = [f"CLUSTER-{index}" for index in range(minimum_n)]
            cells[key] = cell
    return {"metric_cells": cells}


def _evidence() -> dict[str, object]:
    roots = [
        {
            "root_id": f"ROOT-{index}",
            "present": True,
            "expected_sha256": f"{index:x}" * 64,
            "observed_sha256": f"{index:x}" * 64,
        }
        for index in range(1, 9)
    ]
    hashes = {f"ARTIFACT-{index}": f"{index:x}" * 64 for index in range(1, 4)}
    return {
        "governance_roots": roots,
        "isolation_probes": [
            {"probe_id": "PRIVATE-READ", "access_or_write_denied": True},
            {"probe_id": "CLOSED-WRITE", "access_or_write_denied": True},
            {"probe_id": "TAMPER", "access_or_write_denied": True},
        ],
        "restart": {
            "initial_scientific_sha256": hashes,
            "resumed_scientific_sha256": dict(hashes),
        },
    }


def _stage_01_evidence() -> dict[str, object]:
    target = {
        "chromosome": 1,
        "position_bp": 123_456,
        "ref": "A",
        "alt": "G",
        "variant_id": "TARGET-1",
    }
    return {
        "source_sentinels": [
            {"expected_valid": True, "observed_accepted": True},
            {"expected_valid": False, "observed_accepted": False},
        ],
        "covered_autosomes": list(range(1, 23)),
        "expected_target": target,
        "observed_targets": [dict(target)],
        "source_provenance": [
            {
                "assembly": "GRCh38",
                "license": "synthetic-test-data",
                "version": "1.0",
                "sha256": "a" * 64,
            }
        ],
    }


def _stage_02_evidence() -> dict[str, object]:
    return {
        "registry_rows": [
            {
                "source_sample_id": f"SOURCE-{index}",
                "sample_id": f"SAMPLE-{index}",
                "identity_id": f"IDENTITY-{index}",
            }
            for index in range(4_000)
        ],
        "identity_challenges": [
            {"expected_detected": True, "observed_detected": True}
        ],
        "pedigree_challenges": [
            {"expected_flagged": True, "observed_flagged": True}
            for _ in range(200)
        ],
        "genotype_assignments": [{"genotype_source": "MOLECULAR_ASSAY"}],
    }


def _stage_03_evidence() -> dict[str, object]:
    sample_ids = [f"SAMPLE-{index}" for index in range(4_000)]
    return {
        "expected_sample_ids": sample_ids,
        "observed_sample_ids": list(sample_ids),
        "variant_index": {
            "expected_count": 100_000,
            "observed_count": 100_000,
            "expected_sha256": "b" * 64,
            "observed_sha256": "b" * 64,
        },
        "genotype_worlds": [
            {"expected_sha256": "c" * 64, "observed_sha256": "c" * 64}
        ],
        "missingness": {"expected_count": 100, "observed_count": 100},
        "allele_sentinels": [
            {
                "expected_ref": "A",
                "expected_alt": "G",
                "observed_ref": "A",
                "observed_alt": "G",
            }
        ],
    }


def _stage_04_evidence() -> dict[str, object]:
    return {
        "target_calls": [
            {
                "expected_missing": False,
                "observed_missing": False,
                "expected_call": "0/1",
                "observed_call": "0/1",
                "molecular_source": "SYNTHETIC_TRUTH",
                "inference_source": "MOLECULAR",
            }
            for _ in range(4_000)
        ],
        "mendelian_challenges": [
            {"expected_flagged": True, "observed_flagged": True}
            for _ in range(200)
        ],
    }


def _stage_05_evidence() -> dict[str, object]:
    return {
        "high_defect_units": [
            {"unit_id": f"HIGH-{index}", "observed_flagged": True}
            for index in range(500)
        ],
        "low_defect_units": [
            {"unit_id": f"LOW-{index}", "observed_flagged": True}
            for index in range(500)
        ],
        "negative_units": [
            {"unit_id": f"NEGATIVE-{index}", "observed_flagged": False}
            for index in range(3_000)
        ],
        "scenario_units": [
            {
                "scenario_unit_id": f"SCENARIO-{index}",
                "expected_missingness": 0.01,
                "observed_missingness": 0.01,
                "expected_call_rate": 0.99,
                "observed_call_rate": 0.99,
                "evaluable": True,
            }
            for index in range(500)
        ],
    }


def _stage_06_evidence() -> dict[str, object]:
    return {
        "variant_units": [
            {
                "variant_id": f"VARIANT-{index}",
                "expected_eligible": True,
                "observed_selected": True,
            }
            for index in range(500)
        ],
        "chromosome_variant_counts": {
            density: {
                chromosome: minimum for chromosome in range(1, 23)
            }
            for density, minimum in (
                ("D_MAX", 5_000),
                ("D_MID", 3_000),
                ("D_LOW", 1_000),
            )
        },
        "ld_pair_blocks": [
            {"block_id": f"BLOCK-{index}", "exceeds_pruning_threshold": False}
            for index in range(500)
        ],
        "forbidden_variant_uses": [],
    }


def _stage_07_evidence() -> dict[str, object]:
    related = [
        {
            "pair_id": f"RELATED-{index}",
            "expected_class": "DEGREE_1",
            "observed_class": "DEGREE_1",
            "component_id": f"RELATED-COMPONENT-{index}",
            "cryptic": True,
            "expected_kinship": 0.25,
            "observed_kinship": 0.25,
        }
        for index in range(500)
    ]
    unrelated = [
        {
            "pair_id": f"UNRELATED-{index}",
            "expected_class": "UNRELATED",
            "observed_class": "UNRELATED",
            "component_id": f"UNRELATED-COMPONENT-{index}",
            "cryptic": False,
            "expected_kinship": 0.0,
            "observed_kinship": 0.0,
        }
        for index in range(3_000)
    ]
    return {"designated_pairs": related + unrelated}


def _stage_08_evidence() -> dict[str, object]:
    pure_representatives = [
        {
            "representative_id": f"PURE-{index}",
            "pure_state": True,
            "centroid_correct": True,
            "expected_outlier": False,
            "observed_outlier": False,
        }
        for index in range(500)
    ]
    outlier_representatives = [
        {
            "representative_id": f"OUTLIER-{index}",
            "pure_state": False,
            "centroid_correct": False,
            "expected_outlier": True,
            "observed_outlier": True,
        }
        for index in range(500)
    ]
    return {
        "fit_subjects": [{"subject_id": "FIT-1", "is_relative": False}],
        "representatives": pure_representatives + outlier_representatives,
        "procrustes_shifts": [
            {
                "representative_id": f"PURE-{index}",
                "shift_standard_deviation": 0.0,
            }
            for index in range(500)
        ],
    }


def _stage_09_evidence() -> dict[str, object]:
    return {
        "cohort_assignments": [
            {
                "identity_id": f"IDENTITY-{index}",
                "expected_cohort": "CASE",
                "observed_cohort": "CASE",
            }
            for index in range(4_000)
        ],
        "family_components": [
            {
                "component_id": f"COMPONENT-{index}",
                "independent_partition_preserved": True,
            }
            for index in range(200)
        ],
        "mandatory_cohorts": [{"cohort_id": "CASE", "member_count": 4_000}],
        "review": {
            "required": True,
            "signed": True,
            "expected_sha256": "d" * 64,
            "observed_sha256": "d" * 64,
        },
    }


def _stage_10_evidence() -> dict[str, object]:
    return {
        "high_defect_units": [
            {"unit_id": f"HIGH-{index}", "observed_flagged": True}
            for index in range(500)
        ],
        "low_defect_units": [
            {"unit_id": f"LOW-{index}", "observed_flagged": True}
            for index in range(500)
        ],
        "negative_units": [
            {"unit_id": f"NEGATIVE-{index}", "observed_flagged": False}
            for index in range(3_000)
        ],
        "target_exception_checks": [
            {"expected_preserved": True, "observed_preserved": True}
        ],
        "cohort_changes": [{"change_id": "CHANGE-1", "audited": True}],
    }


def _stage_10a_evidence() -> dict[str, object]:
    subject_ids = [f"SUBJECT-{index}" for index in range(4_000)]
    return {
        "expected_subject_ids": subject_ids,
        "observed_subject_ids": list(subject_ids),
        "autosome_variant_counts": {
            chromosome: 1_000 for chromosome in range(1, 23)
        },
        "structural_missing_genotype_count": 0,
        "duplicate_variant_count": 0,
    }


def _stage_11_evidence() -> dict[str, object]:
    target = {
        "chromosome": 1,
        "position_bp": 123_456,
        "ref": "A",
        "alt": "G",
        "variant_id": "TARGET-1",
    }
    return {
        "expected_target": target,
        "observed_target": dict(target),
        "target_age_pairs": [
            {
                "pair_id": f"PAIR-{index}",
                "expected_cm": 1.25,
                "observed_cm": 1.25,
                "expected_start_bp": 100_000,
                "observed_start_bp": 100_000,
                "expected_end_bp": 200_000,
                "observed_end_bp": 200_000,
                "extrapolated": False,
            }
            for index in range(10)
        ],
        "evaluable_target_windows": [
            {"window_id": "WINDOW-1", "density_compliant": True}
        ],
        "infeasible_target_age_pairs": [
            {"expected_infeasible": True, "observed_infeasible": True}
        ],
    }


def test_stage_00_producer_passes_only_from_matching_audit_evidence() -> None:
    bundle = load_metric_rule_bundle(REPOSITORY_ROOT)
    observations = produce_stage_00_observations(_evidence())

    decision = evaluate_stage_observations(
        "00",
        bundle["stages"]["00"],
        observations,
    )

    assert decision["status"] == "PASS"
    assert decision["metric_decision_count"] == 4


def test_stage_00_producer_exposes_hash_tampering_as_point_failure() -> None:
    bundle = load_metric_rule_bundle(REPOSITORY_ROOT)
    evidence = copy.deepcopy(_evidence())
    evidence["governance_roots"][0]["observed_sha256"] = "f" * 64  # type: ignore[index]

    decision = evaluate_stage_observations(
        "00",
        bundle["stages"]["00"],
        produce_stage_00_observations(evidence),
    )

    assert decision["status"] == "FAIL_POINT_ESTIMATE"
    hash_decision = next(
        item for item in decision["metric_decisions"] if item["metric_id"] == "hash_integrity"
    )
    assert hash_decision["estimate"] == pytest.approx(7 / 8)


def test_stage_00_producer_rejects_missing_governance_root() -> None:
    evidence = _evidence()
    evidence["governance_roots"].pop()  # type: ignore[union-attr]

    with pytest.raises(ScientificStageProducerError, match="root_count_mismatch"):
        produce_stage_00_observations(evidence)


@pytest.mark.parametrize(
    ("stage_id", "producer", "evidence_factory", "metric_count"),
    [
        ("01", produce_stage_01_observations, _stage_01_evidence, 5),
        ("02", produce_stage_02_observations, _stage_02_evidence, 4),
        ("03", produce_stage_03_observations, _stage_03_evidence, 5),
        ("04", produce_stage_04_observations, _stage_04_evidence, 4),
    ],
)
def test_stage_01_to_04_producers_pass_complete_closed_evidence(
    stage_id: str,
    producer: object,
    evidence_factory: object,
    metric_count: int,
) -> None:
    bundle = load_metric_rule_bundle(REPOSITORY_ROOT)

    observations = producer(evidence_factory())  # type: ignore[operator]
    decision = evaluate_stage_observations(
        stage_id,
        bundle["stages"][stage_id],
        observations,
    )

    assert decision["status"] == "PASS"
    assert decision["metric_decision_count"] == metric_count


def test_stage_01_producer_exposes_incomplete_autosome_coverage() -> None:
    evidence = _stage_01_evidence()
    evidence["covered_autosomes"].pop()  # type: ignore[union-attr]
    bundle = load_metric_rule_bundle(REPOSITORY_ROOT)

    decision = evaluate_stage_observations(
        "01",
        bundle["stages"]["01"],
        produce_stage_01_observations(evidence),
    )

    assert decision["status"] == "FAIL_POINT_ESTIMATE"


def test_stage_02_producer_exposes_non_bijective_registry() -> None:
    evidence = _stage_02_evidence()
    rows = evidence["registry_rows"]
    rows[1]["identity_id"] = rows[0]["identity_id"]  # type: ignore[index]
    bundle = load_metric_rule_bundle(REPOSITORY_ROOT)

    decision = evaluate_stage_observations(
        "02",
        bundle["stages"]["02"],
        produce_stage_02_observations(evidence),
    )

    assert decision["status"] == "FAIL_POINT_ESTIMATE"


def test_stage_03_producer_exposes_allele_swap() -> None:
    evidence = _stage_03_evidence()
    evidence["allele_sentinels"][0]["observed_ref"] = "G"  # type: ignore[index]
    evidence["allele_sentinels"][0]["observed_alt"] = "A"  # type: ignore[index]
    bundle = load_metric_rule_bundle(REPOSITORY_ROOT)

    decision = evaluate_stage_observations(
        "03",
        bundle["stages"]["03"],
        produce_stage_03_observations(evidence),
    )

    assert decision["status"] == "FAIL_POINT_ESTIMATE"


def test_stage_04_producer_exposes_clinical_inference() -> None:
    evidence = _stage_04_evidence()
    evidence["target_calls"][0]["inference_source"] = "clinical_status"  # type: ignore[index]
    bundle = load_metric_rule_bundle(REPOSITORY_ROOT)

    decision = evaluate_stage_observations(
        "04",
        bundle["stages"]["04"],
        produce_stage_04_observations(evidence),
    )

    assert decision["status"] == "FAIL_POINT_ESTIMATE"


def test_stage_05_producer_passes_with_prespecified_uncertainty() -> None:
    bundle = load_metric_rule_bundle(REPOSITORY_ROOT)

    observations = produce_stage_05_observations(
        _stage_05_evidence(),
        evaluation_secret=b"stage-05-test-evaluation-secret-32-bytes",
        release_id="stage-producer-test",
        bootstrap_replicates=200,
    )
    decision = evaluate_stage_observations(
        "05",
        bundle["stages"]["05"],
        observations,
    )

    assert decision["status"] == "PASS"
    assert decision["metric_decision_count"] == 6
    assert next(
        observation
        for observation in observations
        if observation["metric_id"] == "missingness_error"
    )["producer"]["bootstrap_replicates"] == 200


def test_stage_05_producer_exposes_sensitivity_below_threshold() -> None:
    evidence = _stage_05_evidence()
    for record in evidence["high_defect_units"][:51]:  # type: ignore[index]
        record["observed_flagged"] = False
    bundle = load_metric_rule_bundle(REPOSITORY_ROOT)

    decision = evaluate_stage_observations(
        "05",
        bundle["stages"]["05"],
        produce_stage_05_observations(
            evidence,
            evaluation_secret=b"stage-05-test-evaluation-secret-32-bytes",
            release_id="stage-producer-test",
            bootstrap_replicates=100,
        ),
    )

    assert decision["status"] == "FAIL_POINT_ESTIMATE"


def test_stage_06_producer_passes_all_frozen_density_contexts() -> None:
    bundle = load_metric_rule_bundle(REPOSITORY_ROOT)

    observations = produce_stage_06_observations(_stage_06_evidence())
    decision = evaluate_stage_observations(
        "06",
        bundle["stages"]["06"],
        observations,
    )

    assert decision["status"] == "PASS"
    assert decision["metric_decision_count"] == 8
    assert {
        observation.get("context_id")
        for observation in observations
        if observation["metric_id"] == "density_minimum"
    } == {"density_dmax", "density_dmid", "density_dlow"}


def test_stage_06_producer_exposes_density_below_minimum() -> None:
    evidence = _stage_06_evidence()
    evidence["chromosome_variant_counts"]["D_LOW"][22] = 999  # type: ignore[index]
    bundle = load_metric_rule_bundle(REPOSITORY_ROOT)

    decision = evaluate_stage_observations(
        "06",
        bundle["stages"]["06"],
        produce_stage_06_observations(evidence),
    )

    assert decision["status"] == "FAIL_POINT_ESTIMATE"


def test_stage_07_producer_passes_with_component_cluster_bootstrap() -> None:
    bundle = load_metric_rule_bundle(REPOSITORY_ROOT)

    observations = produce_stage_07_observations(
        _stage_07_evidence(),
        evaluation_secret=b"stage-07-test-evaluation-secret-32-bytes",
        release_id="stage-producer-test",
        bootstrap_replicates=100,
    )
    decision = evaluate_stage_observations(
        "07",
        bundle["stages"]["07"],
        observations,
    )

    assert decision["status"] == "PASS"
    kinship = next(
        observation
        for observation in observations
        if observation["metric_id"] == "kinship_mae"
    )
    assert kinship["independent_n"] == 3_500


def test_stage_08_producer_passes_with_independent_representatives() -> None:
    bundle = load_metric_rule_bundle(REPOSITORY_ROOT)

    observations = produce_stage_08_observations(
        _stage_08_evidence(),
        evaluation_secret=b"stage-08-test-evaluation-secret-32-bytes",
        release_id="stage-producer-test",
        bootstrap_replicates=100,
    )
    decision = evaluate_stage_observations(
        "08",
        bundle["stages"]["08"],
        observations,
    )

    assert decision["status"] == "PASS"
    assert decision["metric_decision_count"] == 6


def test_stage_08_producer_rejects_unknown_procrustes_identity() -> None:
    evidence = _stage_08_evidence()
    evidence["procrustes_shifts"][0]["representative_id"] = "UNKNOWN"  # type: ignore[index]

    with pytest.raises(ScientificStageProducerError, match="identity_unknown"):
        produce_stage_08_observations(
            evidence,
            evaluation_secret=b"stage-08-test-evaluation-secret-32-bytes",
            release_id="stage-producer-test",
            bootstrap_replicates=10,
        )


def test_stage_09_producer_passes_signed_cohort_freeze() -> None:
    bundle = load_metric_rule_bundle(REPOSITORY_ROOT)

    decision = evaluate_stage_observations(
        "09",
        bundle["stages"]["09"],
        produce_stage_09_observations(_stage_09_evidence()),
    )

    assert decision["status"] == "PASS"


def test_stage_09_unsigned_review_requires_review_instead_of_passing() -> None:
    evidence = _stage_09_evidence()
    evidence["review"]["signed"] = False  # type: ignore[index]
    evidence["review"]["observed_sha256"] = None  # type: ignore[index]
    bundle = load_metric_rule_bundle(REPOSITORY_ROOT)

    decision = evaluate_stage_observations(
        "09",
        bundle["stages"]["09"],
        produce_stage_09_observations(evidence),
    )

    assert decision["status"] == "REVIEW_REQUIRED"


def test_stage_10_producer_passes_final_qc_contract() -> None:
    bundle = load_metric_rule_bundle(REPOSITORY_ROOT)

    decision = evaluate_stage_observations(
        "10",
        bundle["stages"]["10"],
        produce_stage_10_observations(_stage_10_evidence()),
    )

    assert decision["status"] == "PASS"


def test_stage_10a_producer_passes_structural_panel_contract() -> None:
    bundle = load_metric_rule_bundle(REPOSITORY_ROOT)

    decision = evaluate_stage_observations(
        "10A",
        bundle["stages"]["10A"],
        produce_stage_10a_observations(_stage_10a_evidence()),
    )

    assert decision["status"] == "PASS"
    assert decision["metric_decision_count"] == 6


def test_stage_10a_producer_exposes_subject_order_change() -> None:
    evidence = _stage_10a_evidence()
    observed = evidence["observed_subject_ids"]
    observed[0], observed[1] = observed[1], observed[0]  # type: ignore[index]
    bundle = load_metric_rule_bundle(REPOSITORY_ROOT)

    decision = evaluate_stage_observations(
        "10A",
        bundle["stages"]["10A"],
        produce_stage_10a_observations(evidence),
    )

    assert decision["status"] == "FAIL_POINT_ESTIMATE"


def test_stage_11_producer_passes_target_map_and_window_contract() -> None:
    bundle = load_metric_rule_bundle(REPOSITORY_ROOT)

    decision = evaluate_stage_observations(
        "11",
        bundle["stages"]["11"],
        produce_stage_11_observations(_stage_11_evidence()),
    )

    assert decision["status"] == "PASS"
    assert decision["metric_decision_count"] == 6


def test_stage_11_producer_exposes_map_error_above_one_bp() -> None:
    evidence = _stage_11_evidence()
    evidence["target_age_pairs"][0]["observed_start_bp"] = 100_002  # type: ignore[index]
    bundle = load_metric_rule_bundle(REPOSITORY_ROOT)

    decision = evaluate_stage_observations(
        "11",
        bundle["stages"]["11"],
        produce_stage_11_observations(evidence),
    )

    assert decision["status"] == "FAIL_POINT_ESTIMATE"


@pytest.mark.parametrize("stage_id", REMAINING_STAGE_IDS)
def test_remaining_stage_producers_cover_and_pass_frozen_rules(stage_id: str) -> None:
    bundle = load_metric_rule_bundle(REPOSITORY_ROOT)
    stage_rules = bundle["stages"][stage_id]

    observations = produce_remaining_stage_observations(
        stage_id,
        stage_rules,
        _passing_remaining_evidence(stage_rules),
        evaluation_secret=b"remaining-stage-evaluation-secret-32-bytes",
        release_id="stage-producer-test",
        bootstrap_replicates=10,
    )
    decision = evaluate_stage_observations(stage_id, stage_rules, observations)

    assert decision["status"] == "PASS"


def test_stage_producer_audit_remains_partial() -> None:
    audit = stage_producer_audit()

    assert audit["implemented_stages"] == [
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
    ]
    assert audit["implemented_stage_count"] == 25
    assert audit["missing_stage_count"] == 0
    assert audit["stage_producers_complete"] is True
    assert audit["scientific_evaluation_authorized"] is False
