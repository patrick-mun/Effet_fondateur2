"""Rapports D23 : agrégats fermés, illustration et mathématiques, sans outil réel."""
from copy import deepcopy
import json
from pathlib import Path
import re

import numpy as np
import pytest

from effet_fondateur.calibration.evaluate_phasing import STAGES, acceptance_results, validate_evaluation
from effet_fondateur.calibration.phasing_bench import CASES, METHOD_ID, load_configuration
from effet_fondateur.calibration.phasing_oracle import PhasePanel, score_phase
from effet_fondateur.calibration.qc_bench import file_inventory, write_json
from effet_fondateur.calibration.report_phasing import build_reports
from effet_fondateur.calibration.technical_evaluation import GenotypeMatrix, compare_genotypes


def aggregate_documents() -> dict[str, dict]:
    """Petite vérité artificielle, avec tous les champs autorisés de chaque rapport."""
    sample_count = 40
    calls = np.tile(np.array([[[0, 1]]], dtype=np.int8), (7, sample_count, 1))
    panel = PhasePanel(tuple(f"SYNTHETIC_{index}" for index in range(sample_count)), tuple((5, index+1, "A", "G") for index in range(7)), calls,
                       np.ones((7, sample_count), dtype=bool), np.full((7, sample_count), "", dtype=object))
    performance = score_phase(panel, panel, founder_indices=list(range(32)), unit_indices=[index//3 for index in range(24)]+list(range(8,24)), target=(5,4),
                              rare_maf_upper_bound=.05, density_gap_bp=100_000, bootstrap_seed=2305, bootstrap_replicates=20)
    performance.update(status="PASS", sample_order_correct=True, variant_order_correct=True, reference_membership_correct=True, target_retained=True, observed_genotype_preservation_errors=0,
                       map_coordinate_error_count=0, published_manual_review_required=True, reliable_carrier_count=20, unreliable_carrier_count=6,
                       chromosome_counts=[dict(chromosome=5, observed=7)])
    matrix = GenotypeMatrix(tuple(("F"+str(index), sample) for index, sample in enumerate(panel.samples)), tuple((5,index+1,"m"+str(index)) for index in range(7)), calls)
    input_row = dict(status="PASS", concordance=compare_genotypes(matrix,matrix,stage="10A"), source_markers=7, expected_markers=7, observed_markers=7,
                     excluded_incomplete=0, excluded_monomorphic=0, chromosome_counts=[dict(chromosome=chromosome, expected=7 if chromosome==5 else 0, observed=7 if chromosome==5 else 0) for chromosome in range(1,23)])
    region_row = dict(status="PASS", concordance=compare_genotypes(matrix,matrix,stage="11"), expected_markers=7, observed_markers=7,
                      map_coordinate_error_count=0, map_order_correct=True, map_max_absolute_error_cm=0., target_retained=True,
                      map_plot=dict(position_bp=[1,2,3], expected_cm=[.1,.2,.3], observed_cm=[.1,.2,.3]))
    result = {}
    for stage in STAGES:
        row = input_row if stage == "10A" else region_row if stage == "11" else performance
        result[stage] = dict(schema_version="1.0.0", method_id=METHOD_ID, stage_id=stage, evaluations=[dict(profile=case, **deepcopy(row)) for case in CASES],
                             status="PASS", performance_qualification="REVIEW_REQUIRED" if stage.startswith("12") else "NOT_APPLICABLE", calibrated=False, founder_effect_proven=False)
        for evaluation in result[stage]["evaluations"]:
            evaluation["criteria"] = acceptance_results(stage, evaluation)
    return result


def closed_fixture(workspace: Path) -> Path:
    root = workspace / "evaluation"
    root.mkdir()
    for stage, document in aggregate_documents().items():
        validate_evaluation(document)
        write_json(root / (stage + ".evaluation.json"), document)
    config = load_configuration(Path(__file__).resolve().parents[2] / "calibration/configurations/d23_phasing_bench_v1.json")
    manifest = dict(schema_version="1.0.0", method_id=METHOD_ID, closure_state="CLOSED", created_at="2026-10-30T12:00:00Z",
                    preparation_manifest_sha256=config["preparation_manifest_sha256"], source_run_sha256="a"*64, regional_source_sha256=["b"*64]*3,
                    configuration_sha256="c"*64, implementation_sha256="d"*64, oracle_implementation_sha256="e"*64, config=config,
                    acceptance_criteria_sha256="a"*64, prespecification_sha256="b"*64,
                    stages=[dict(stage_id=stage, status="PASS", performance_qualification="REVIEW_REQUIRED" if stage.startswith("12") else "NOT_APPLICABLE") for stage in STAGES],
                    independent_units=24, limitations=["PERFORMANCE_THRESHOLDS_NOT_SEALED"], artifacts=file_inventory(root), calibrated=False, founder_effect_proven=False)
    validate_evaluation(manifest, manifest=True)
    write_json(root / "run_index.json", manifest)
    return root


def test_four_reports_aggregate_only_illustrated_and_no_requalification(tmp_path: Path) -> None:
    workspace = tmp_path / "calibration_workspace"
    workspace.mkdir()
    evaluation = closed_fixture(workspace)
    original = (evaluation / "run_index.json").read_bytes()
    template = tmp_path / "template.html"
    template.write_text('<html><head><title>Template</title></head><body></body></html>')
    paths = build_reports(evaluation_root=evaluation, template_path=template, output_root=workspace / "reports")
    assert len(paths) == 4
    for path in paths:
        html = path.read_text()
        assert len(re.findall('<section id="s[0-9]+">', html)) == 14
        assert html.count('role="math"') == 2
        assert 2 <= html.count('<figure class="scientific-figure">') <= 4
        data = (path.parent / "report_data.json").read_text()
        assert "SYNTHETIC_" not in html + data
        assert "calibrated=false" in html and "REVIEW_REQUIRED" in html
    assert (evaluation / "run_index.json").read_bytes() == original
    with pytest.raises(ValueError, match="run_exists"):
        build_reports(evaluation_root=evaluation, template_path=template, output_root=workspace / "reports")


def test_evaluation_schema_refuses_individual_fields_and_false_qualification() -> None:
    document = aggregate_documents()["12"]
    document["evaluations"][0]["sample_ids"] = ["SYNTHETIC_NEVER_EXPORT"]
    with pytest.raises(ValueError, match="aggregate_evaluation_schema_invalid"):
        validate_evaluation(document)
    document = aggregate_documents()["12"]
    document["performance_qualification"] = "PASS"
    with pytest.raises(ValueError, match="aggregate_evaluation_schema_invalid"):
        validate_evaluation(document)
