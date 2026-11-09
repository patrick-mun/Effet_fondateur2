"""Préservation du manque, des axes et de la référence : aucun exécutable réel."""
from dataclasses import replace
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pytest

from effet_fondateur.calibration.evaluate_phasing import evaluate_phase_outputs
from effet_fondateur.calibration.phasing_oracle import PhasePanel
from effet_fondateur.calibration.qc_bench import write_json
from effet_fondateur.calibration.simulation.measurement_errors import MeasurementPanel


@pytest.mark.parametrize("case,expected_errors", [("CORRECT_REMASK",0),("IMPUTED_AS_OBSERVED",1),("PARTIAL_REMASK",1),("UNPHASED",0),("REFERENCE_CHANGED",0),("LOCUS_ORDER_CHANGED",0),("LOCUS_ID_CHANGED",0)])
def test_publication_oracle_detects_non_scientific_conformance_defects(tmp_path: Path, monkeypatch: pytest.MonkeyPatch, case: str, expected_errors: int) -> None:
    calls = np.tile(np.array([[[0,1],[0,1]]], dtype=np.int8), (7,1,1))
    samples = ("SYNTHETIC1", "SYNTHETIC2")
    truth = PhasePanel(samples, tuple((5,index+1,"A","G") for index in range(7)), calls, np.ones((7,2), dtype=bool), np.full((7,2),"",dtype=object),tuple("m"+str(index) for index in range(7)))
    observation = calls.copy()
    observation[1,0] = -1
    variants = tuple(("chr5",str(index+1),"m"+str(index),"A","G",".","PASS",".","GT") for index in range(7))
    prediction = observation.copy()
    if case == "IMPUTED_AS_OBSERVED":
        prediction[1,0] = (0,1)
    if case == "PARTIAL_REMASK":
        prediction[1,0] = (0,-1)
    predicted = replace(truth, alleles=prediction, phased=truth.phased.copy())
    if case == "UNPHASED":
        predicted.phased[2,0] = False
    if case == "LOCUS_ORDER_CHANGED":
        predicted = replace(predicted, loci=predicted.loci[::-1], alleles=predicted.alleles[::-1], phased=predicted.phased[::-1], phase_sets=predicted.phase_sets[::-1])
    if case == "LOCUS_ID_CHANGED":
        predicted = replace(predicted, variant_ids=("OTHER_ID",)+predicted.variant_ids[1:])
    monkeypatch.setattr("effet_fondateur.calibration.evaluate_phasing._read_phased", lambda path, timeout: predicted)
    header = "REF_OTHER\n" if case == "REFERENCE_CHANGED" else "REF1\nREF2\n"
    monkeypatch.setattr("effet_fondateur.calibration.evaluate_phasing.subprocess.run", lambda *args,**kwargs: SimpleNamespace(returncode=0, stdout=header))
    manifest = tmp_path / "reference.json"
    write_json(manifest, dict(sample_selection=dict(membership_sha256="a"*64)))
    summary = tmp_path / "summary.json"
    write_json(summary, dict(manual_validation_required=True, reliable_carrier_count=1, unreliable_carrier_count=1))
    paths = dict(shapeit5_final_bcf=tmp_path/"unused.bcf", reference_window_manifest=manifest, reference_window_vcf=tmp_path/"reference.vcf.gz",phase_target_region_summary=summary)
    registry = [dict(SAMPLE_ID=sample, FID="F"+str(index), PID="0",MID="0") for index,sample in enumerate(samples)]
    config = dict(command_timeout_seconds=1, rare_maf_upper_bound=.05, density_gap_bp=100_000, bootstrap_seed=2305, bootstrap_replicates=20)
    result = evaluate_phase_outputs(paths,"12",truth,registry,config,dict(chromosome=5,position_bp=4),list(samples),dict(membership_sha256="a"*64,sample_ids=["REF1","REF2"]),MeasurementPanel(samples,variants,observation),tmp_path)
    assert result["observed_genotype_preservation_errors"] == expected_errors
    assert result["status"] == ("PASS" if case == "CORRECT_REMASK" else "FAIL")
