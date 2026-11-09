"""Contrastes D23b analytiques, sans données ni outils externes."""
import numpy as np
import pytest

from effet_fondateur.calibration.evaluate_phasing_sensitivity import _contrasts, paired_rate_contrast


def test_paired_contrast_reconciles_counts_and_is_reproducible() -> None:
    first = [{"numerator": 2, "denominator": 10}, {"numerator": 1, "denominator": 5}]
    second = [{"numerator": 1, "denominator": 10}, {"numerator": 0, "denominator": 5}]
    result = paired_rate_contrast(first, second, seed=2311, replicates=100)
    assert result == paired_rate_contrast(first, second, seed=2311, replicates=100)
    assert result["estimate"] == pytest.approx(2 / 15)
    assert result["bootstrap_valid_replicates"] == 100
    assert result["bootstrap_excluded_zero_denominator"] == 0
    assert np.allclose(result["interval_95"], [0.1, 0.2])


def test_paired_contrast_counts_zero_denominator_draws() -> None:
    first = [{"numerator": 0, "denominator": 0}, {"numerator": 1, "denominator": 2}]
    second = [{"numerator": 0, "denominator": 1}, {"numerator": 0, "denominator": 0}]
    result = paired_rate_contrast(first, second, seed=1, replicates=100)
    assert result["status"] == "EVALUABLE"
    assert result["bootstrap_valid_replicates"] + result["bootstrap_excluded_zero_denominator"] == 100
    with pytest.raises(ValueError, match="paired_counts_invalid"):
        paired_rate_contrast(first, second[:-1], seed=1, replicates=10)


def test_failed_source_cell_makes_contrasts_not_evaluable() -> None:
    modes = ("WITH_PEDIGREE", "WITHOUT_PEDIGREE")
    performance = {
        **{("AUTOSOMAL", profile, "WITHOUT_PEDIGREE_BY_CONTRACT"): {"status": "FAIL"}
           for profile in ("DENSITY_25", "DENSITY_50", "DENSITY_75", "DENSITY_100")},
        **{("REGIONAL", profile, mode): {"status": "FAIL"}
           for profile in ("DENSITY_25", "DENSITY_50", "DENSITY_75", "DENSITY_100") for mode in modes},
    }
    rows = _contrasts(performance, {"bootstrap_seed": 2311, "bootstrap_replicates": 20})
    assert len(rows) == 65
    assert {row["status"] for row in rows} == {"NOT_EVALUABLE"}
    assert {row["reason"] for row in rows} == {"SOURCE_CELL_NOT_EVALUABLE"}
