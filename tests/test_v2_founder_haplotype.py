import csv
import json
import sys
from pathlib import Path

from effet_fondateur.contracts import validate_tsv_table
from effet_fondateur.founder import infer_target_centered_ibs


def _write_tsv(path: Path, columns: list[str], rows: list[list[str]]) -> None:
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(columns)
        writer.writerows(rows)


def _write_fake_bcftools(path: Path) -> None:
    script = f'''#!{sys.executable}
import sys

arguments = sys.argv[1:]
if arguments == ["--version"]:
    print("bcftools 1.21")
elif arguments[:2] == ["query", "--list-samples"]:
    print("carrier_1\\ncarrier_2\\ncarrier_3\\ncontrol_1")
elif arguments[:2] == ["query", "--format"]:
    print("v80\\t80\\t0|1\\t0|0\\t1|0\\t0|1")
    print("v90\\t90\\t1|0\\t1|1\\t0|1\\t1|0")
    print("target\\t100\\t1|0\\t1|0\\t0|1\\t0|0")
    print("v110\\t110\\t1|0\\t1|1\\t0|1\\t1|0")
    print("v120\\t120\\t1|0\\t1|0\\t0|0\\t0|0")
else:
    raise SystemExit(9)
'''
    path.write_text(script, encoding="utf-8")
    path.chmod(0o755)


def _prepare_inputs(tmp_path: Path) -> dict[str, Path]:
    bcf_path = tmp_path / "phased.bcf"
    bcf_path.write_bytes(b"synthetic")
    map_path = tmp_path / "target_genetic_map.tsv"
    map_columns = [
        "DATASET_ID", "VARIANT_ORDER", "VARIANT_ID", "CHROMOSOME",
        "POSITION_BP", "POSITION_CM", "MAP_STATUS", "LEFT_MAP_BP",
        "LEFT_MAP_CM", "RIGHT_MAP_BP", "RIGHT_MAP_CM", "LOCAL_RATE_CM_PER_MB",
        "IS_TARGET_VARIANT",
    ]
    map_rows = []
    for order, (variant_id, position, position_cm) in enumerate(
        [("v80", 80, "0.8"), ("v90", 90, "0.9"), ("target", 100, "1.0"), ("v110", 110, "1.1"), ("v120", 120, "1.2")],
        start=1,
    ):
        map_rows.append([
            "dataset", str(order), variant_id, "19", str(position), position_cm,
            "EXACT", str(position), position_cm, str(position), position_cm, "1",
            "true" if variant_id == "target" else "false",
        ])
    _write_tsv(map_path, map_columns, map_rows)
    carrier_path = tmp_path / "carrier_haplotypes.tsv"
    carrier_columns = [
        "SAMPLE_ORDER", "SAMPLE_ID", "TARGET_VARIANT_ID", "EXPLICIT_GENOTYPE",
        "PHASED_GT", "ALT_COPY_COUNT", "CARRIER_HAPLOTYPE", "PHASE_CONFIDENCE",
        "CONFIDENCE_STATUS", "RELIABILITY_STATUS", "UNRELIABLE_REASON",
    ]
    _write_tsv(carrier_path, carrier_columns, [
        ["1", "carrier_1", "target", "A/G", "1|0", "1", "H1", "0.99", "SCORED_PASS", "PASS", ""],
        ["2", "carrier_2", "target", "A/G", "1|0", "1", "H1", "0.98", "SCORED_PASS", "PASS", ""],
        ["3", "carrier_3", "target", "A/G", "0|1", "1", "H2", "0.97", "SCORED_PASS", "PASS", ""],
        ["4", "control_1", "target", "A/A", "0|0", "0", "NONE", "", "NOT_APPLICABLE_HOMOZYGOUS", "PASS", ""],
    ])
    cohorts_path = tmp_path / "cohorts_frozen.tsv"
    cohort_columns = [
        "COHORT_ID", "SAMPLE_ID", "ROLE", "INCLUDED", "EXCLUSION_CODE",
        "DECISION_SOURCE", "INDEPENDENT_UNIT_ID", "FAMILY_REPRESENTATIVE",
    ]
    cohort_rows = [
        ["target_carriers_independent", f"carrier_{index}", "CARRIER", "true", "", "test", f"unit_{index}", "true"]
        for index in range(1, 4)
    ] + [
        ["controls_unrelated", "control_1", "CONTROL", "true", "", "test", "control_unit_1", "true"],
        ["target_chromosome_all_qc", "control_1", "NON_CARRIER", "true", "", "test", "", "false"],
    ]
    _write_tsv(cohorts_path, cohort_columns, cohort_rows)
    samples_path = tmp_path / "samples.master.tsv"
    sample_columns = [
        "SAMPLE_ID", "SOURCE_FILE", "FID", "IID", "PID", "MID", "SEX",
        "CLINICAL_STATUS", "GROUP_LABEL", "TARGET_GENOTYPE",
        "TARGET_GENOTYPE_SOURCE", "ARRAY_BATCH", "INCLUDE_GENOMEWIDE",
        "INCLUDE_TARGET_CHROMOSOME", "NOTES_CODE",
    ]
    sample_rows = [
        [sample_id, f"{sample_id}.txt", family_id, sample_id, "0", "0", "UNKNOWN", "UNKNOWN", group, genotype, "explicit", "batch", "true", "true", ""]
        for sample_id, family_id, group, genotype in [
            ("carrier_1", "family_1", "CASE", "A/G"),
            ("carrier_2", "family_2", "CASE", "A/G"),
            ("carrier_3", "family_3", "CASE", "A/G"),
            ("control_1", "family_4", "CONTROL", "A/A"),
        ]
    ]
    _write_tsv(samples_path, sample_columns, sample_rows)
    bcftools_path = tmp_path / "bcftools"
    _write_fake_bcftools(bcftools_path)
    return {
        "bcf": bcf_path, "map": map_path, "carriers": carrier_path,
        "cohorts": cohorts_path, "samples": samples_path, "bcftools": bcftools_path,
    }


def test_exact_ibs_publishes_conservative_founder_candidate(tmp_path: Path) -> None:
    paths = _prepare_inputs(tmp_path)

    result = infer_target_centered_ibs(
        phased_bcf_path=paths["bcf"], carrier_haplotypes_path=paths["carriers"],
        cohorts_path=paths["cohorts"], samples_master_path=paths["samples"],
        genetic_map_path=paths["map"], output_dir=tmp_path / "output",
        bcftools_command=str(paths["bcftools"]), timeout_seconds=30,
        minimum_independent_carriers=3, minimum_flank_markers=1,
    )

    assert result.status == "SUPPORTED_IBS_CANDIDATE"
    segment_rows = validate_tsv_table(
        result.segments_path, "founder_segments.schema.json"
    ).rows
    assert len(segment_rows) == 3
    assert {row["LEFT_BOUND_BP"] for row in segment_rows} == {"80"}
    assert {row["RIGHT_BOUND_BP"] for row in segment_rows} == {"110"}
    consensus = validate_tsv_table(
        result.consensus_path, "founder_consensus.schema.json"
    ).rows[0]
    assert consensus["INTERPRETATION"] == "IBS_SHARED_CANDIDATE"
    assert consensus["BACKGROUND_HAPLOTYPE_COUNT"] == "2"
    assert consensus["MATCHING_BACKGROUND_HAPLOTYPE_COUNT"] == "1"
    assert consensus["BACKGROUND_MATCH_FREQUENCY"] == "0.5"
    sharing = validate_tsv_table(
        result.sharing_path, "founder_sharing_matrix.schema.json"
    )
    assert sharing.row_count == 6
    summary = json.loads(result.summary_path.read_text(encoding="utf-8"))
    assert summary["ibd_claimed"] is False


def test_exact_ibs_refuses_founder_conclusion_when_carriers_are_insufficient(
    tmp_path: Path,
) -> None:
    paths = _prepare_inputs(tmp_path)

    result = infer_target_centered_ibs(
        phased_bcf_path=paths["bcf"], carrier_haplotypes_path=paths["carriers"],
        cohorts_path=paths["cohorts"], samples_master_path=paths["samples"],
        genetic_map_path=paths["map"], output_dir=tmp_path / "output",
        bcftools_command=str(paths["bcftools"]), timeout_seconds=30,
        minimum_independent_carriers=4, minimum_flank_markers=1,
    )

    assert result.status == "INSUFFICIENT_CARRIERS"
    consensus = validate_tsv_table(
        result.consensus_path, "founder_consensus.schema.json"
    ).rows[0]
    assert consensus["INTERPRETATION"] == "NO_FOUNDER_CONCLUSION"
