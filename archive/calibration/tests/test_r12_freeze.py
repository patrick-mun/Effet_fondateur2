import gzip
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pytest
import yaml

from effet_fondateur.audit import sha256_file
from effet_fondateur.calibration.r12_freeze import (
    COMPATIBILITY,
    R12FreezeError,
    _ibd_pairs,
    _switch_error,
    _validate_compatibility,
)


REPOSITORY_ROOT = Path(__file__).resolve().parents[2]


def test_r12_compatibility_document_is_generic_and_unqualified() -> None:
    document = yaml.safe_load((REPOSITORY_ROOT / COMPATIBILITY).read_text())

    assert document["freeze_phase"] == "R12"
    assert document["target_contract"] == {
        "selection": "PROFILE_CONFIGURED_AUTOSOMAL",
        "hard_coded": False,
        "autosomal_only": True,
    }
    assert document["qualification_flags"] == {
        "calibrated": False,
        "founder_effect_proven": False,
        "production_promotion": False,
    }


def test_switch_error_counts_orientation_changes() -> None:
    truth = np.array([[[0, 0, 0, 0], [1, 1, 1, 1]]], dtype=np.int8)
    observed = np.array(
        [[[0, 1], [0, 1], [1, 0], [1, 0]]], dtype=np.int8
    )

    switches, transitions, rate = _switch_error(
        truth, observed, slice(0, truth.shape[2])
    )

    assert switches == 1
    assert transitions == 3
    assert rate == pytest.approx(1 / 3)


def test_ibd_pairs_counts_segments_and_unique_pairs(tmp_path: Path) -> None:
    output = tmp_path / "segments.ibd.gz"
    with gzip.open(output, "wt", encoding="utf-8") as handle:
        handle.write("S1\t1\tS2\t2\t7\t100\t200\t2.0\n")
        handle.write("S2\t1\tS1\t2\t7\t300\t500\t3.0\n")
        handle.write("S1\t1\tS3\t1\t7\t700\t900\t2.0\n")

    segment_count, pairs = _ibd_pairs(output)

    assert segment_count == 3
    assert pairs == {("S1", "S2"), ("S1", "S3")}


def test_tool_compatibility_rejects_a_version_drift(tmp_path: Path) -> None:
    hap_jar = tmp_path / "hap.jar"
    refined_jar = tmp_path / "refined.jar"
    hap_jar.write_bytes(b"hap")
    refined_jar.write_bytes(b"refined")
    compatibility = {
        "freeze_phase": "R12",
        "tools": {
            "hap_ibd": {"sha256": sha256_file(hap_jar)},
            "refined_ibd": {"sha256": sha256_file(refined_jar)},
        },
        "qualification_flags": {
            "calibrated": False,
            "founder_effect_proven": False,
            "production_promotion": False,
        },
    }
    tools = {
        "probe_text": {
            "plink": "PLINK v1.9.0-b.7.11 64-bit",
            "king": "KING 2.3.1",
            "bcftools": "bcftools 1.21",
            "java": 'openjdk version "17.0.21"',
            "hap_ibd": "hap-ibd.jar [ version 1.0, 15Jun23.92f ]",
            "refined_ibd": "refined-ibd.17Jan20.102.jar",
        },
        "shapeit_probe": SimpleNamespace(
            phase_common_version="5.1.1", phase_rare_version="5.1.1"
        ),
        "hap_ibd_jar": str(hap_jar),
        "refined_ibd_jar": str(refined_jar),
    }

    with pytest.raises(
        R12FreezeError, match="r12_tool_version_incompatible:king"
    ):
        _validate_compatibility(compatibility, tools)
