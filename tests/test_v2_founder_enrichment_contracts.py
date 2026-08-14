import csv
import gzip
from pathlib import Path

import pytest

from effet_fondateur.contracts import TableValidationError, validate_json_document, validate_tsv_table
from effet_fondateur.contracts.documents import DocumentValidationError
from effet_fondateur.stages.evaluate_founder_haplotype_enrichment import (
    _load_target_metadata,
)
from effet_fondateur.founder_enrichment.publication import write_tsv


def _write(path: Path, columns: list[str], row: list[str]) -> None:
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(columns)
        writer.writerow(row)


def test_founder_enrichment_unit_contract_is_strict(tmp_path: Path) -> None:
    path = tmp_path / "units.tsv"
    columns = ["UNIT_ID", "FAMILY_ID", "SAMPLE_ID", "HAPLOTYPE_ID", "ROLE", "SELECTION_SOURCE", "STATUS", "EXCLUSION_CODE"]
    _write(path, columns, ["unit_1", "family_1", "sample_1", "H1", "INDEPENDENT_CARRIER_FAMILY", "STEP13_PRESELECTED_REPRESENTATIVE", "INCLUDED", ""])
    assert validate_tsv_table(path, "founder_haplotype_units.schema.json").row_count == 1
    _write(path, columns, ["unit_1", "family_1", "sample_1", "H1", "INDEPENDENT_CARRIER_FAMILY", "BEST_A_POSTERIORI", "INCLUDED", ""])
    with pytest.raises(TableValidationError):
        validate_tsv_table(path, "founder_haplotype_units.schema.json")


def test_founder_enrichment_summary_forbids_scientific_overclaim() -> None:
    sha = "a" * 64
    summary = {
        "schema_version": "1.0.0", "method_id": "target_centered_empirical_haplotype_sharing_v1",
        "primary_statistic": "total_shared_cm", "status": "NOT_CLASSIFIED", "independent_family_count": 3,
        "observed": {"evaluation_status": "EVALUATED", "left_shared_cm": 0.4, "right_shared_cm": 0.9, "total_shared_cm": 1.3, "left_marker_count": 5, "right_marker_count": 9},
        "null_results": [{"source": "EXTERNAL", "stratum": "ALL", "requested_draws": 100000, "attempted_draws": 100010, "evaluable_draws": 100000, "non_evaluable_draws": 10, "exceedance_count": 4, "empirical_probability": 5 / 100001, "exact_probability": None, "interval_low": 0.0, "interval_high": 0.0001}],
        "classification_threshold": None, "random_seed": 161602026,
        "provenance": {"assembly": "GRCh38", "target_variant_id": "target", "target_ref": "C", "target_alt": "A", "map_sha256": sha, "study_bcf_sha256": sha, "reference_vcf_sha256": sha, "step13_summary_sha256": sha, "step16a_summary_sha256": sha},
        "interpretation": {"ibs_only": True, "ibd_proven": False, "founder_effect_proven": False, "geographic_origin_inferred": False, "composite_score_calculated": False, "statement": "Partage IBS centré cible, pas preuve IBD."},
    }
    validate_json_document(summary, "founder_haplotype_enrichment_summary.schema.json")
    summary["interpretation"]["ibd_proven"] = True
    with pytest.raises(DocumentValidationError):
        validate_json_document(summary, "founder_haplotype_enrichment_summary.schema.json")


def test_null_draw_contract_can_be_validated_while_gzipped(tmp_path: Path) -> None:
    path = tmp_path / "draws.tsv.gz"
    columns = ["NULL_SOURCE", "STRATUM", "DRAW_INDEX", "ATTEMPT_INDEX", "UNIT_COUNT", "EVALUATION_STATUS", "LEFT_SHARED_CM", "RIGHT_SHARED_CM", "TOTAL_SHARED_CM", "LEFT_MARKER_COUNT", "RIGHT_MARKER_COUNT", "NON_EVALUABLE_REASON"]
    with gzip.open(path, "wt", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(columns)
        writer.writerow(["EXTERNAL", "ALL", "1", "1", "3", "EVALUATED", "0.4", "0.9", "1.3", "5", "9", ""])
    assert validate_tsv_table(path, "founder_haplotype_null_draws.schema.json").row_count == 1


def test_founder_enrichment_loads_target_metadata_from_yaml(tmp_path: Path) -> None:
    path = tmp_path / "target_variant.yaml"
    path.write_text(
        "\n".join(
            (
                "schema_version: 1.0.0",
                "assembly: GRCh38",
                "gene: DOCK6",
                "chromosome: 19",
                "position_bp: 11222390",
                "ref: C",
                "alt: A",
                "transcript: NM_020812.4",
                "hgvs_c: c.1833-1G>T",
                "hgvs_p: null",
                "project_variant_id: DOCK6_GRCh38_19_11222390_C_A",
                "rsid: null",
                "reference_validation:",
                "  status: CONFIRMED",
                "  source_id: test_reference",
                "  confirmed_at: '2026-08-13T12:00:00+04:00'",
                "annotation_validation:",
                "  status: CONFIRMED",
                "  source_id: test_annotation",
                "  confirmed_at: '2026-08-13T12:00:00+04:00'",
                "",
            )
        ),
        encoding="utf-8",
    )

    target = _load_target_metadata(path)

    assert target["project_variant_id"] == "DOCK6_GRCh38_19_11222390_C_A"
    assert target["position_bp"] == 11222390


def test_founder_enrichment_writes_contract_boolean_values(tmp_path: Path) -> None:
    path = tmp_path / "variant_audit.tsv"
    columns = [
        "VARIANT_ORDER", "VARIANT_ID", "CHROMOSOME", "POSITION_BP",
        "POSITION_CM", "IS_TARGET", "OBSERVED_USE", "INTERNAL_NULL_USE",
        "EXTERNAL_NULL_USE", "DECISION_REASON",
    ]
    write_tsv(
        path,
        columns,
        [
            {
                "VARIANT_ORDER": 1,
                "VARIANT_ID": "target",
                "CHROMOSOME": 19,
                "POSITION_BP": 11222390,
                "POSITION_CM": "1.25",
                "IS_TARGET": True,
                "OBSERVED_USE": "ANCHOR_ONLY",
                "INTERNAL_NULL_USE": "ANCHOR_ONLY",
                "EXTERNAL_NULL_USE": "ANCHOR_ONLY",
                "DECISION_REASON": "TARGET_EXCLUDED_FROM_SIGNATURE",
            }
        ],
    )

    assert "\ttrue\t" in path.read_text(encoding="utf-8")
    assert validate_tsv_table(
        path, "founder_haplotype_variant_audit.schema.json"
    ).row_count == 1
