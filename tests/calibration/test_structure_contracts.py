"""Régression D22 : effectif de référence canonique >= 2, sans outil externe."""

from pathlib import Path

import pytest

from effet_fondateur.calibration.qc_bench import write_table
from effet_fondateur.contracts import TableValidationError, validate_tsv_table


def eigenvalue_fixture(path: Path, count: str) -> None:
    """Une ligne minimale isole la validation du compte, pas le calcul PCA."""
    row = dict(COMPONENT="PC1", EIGENVALUE="1", EXPLAINED_VARIANCE_RATIO="1", CUMULATIVE_EXPLAINED_VARIANCE_RATIO="1",
               SINGULAR_VALUE="1", REFERENCE_SAMPLE_COUNT=count, INFORMATIVE_VARIANT_COUNT="10")
    write_table(path, tuple(row), [row])


@pytest.mark.parametrize("count", ["2", "9", "10", "19", "20", "99", "100", "108", "199", "200", "1000"])
def test_population_eigenvalues_accept_all_reference_counts_at_least_two(tmp_path, count):
    path = tmp_path / "eigenvalues.tsv"
    eigenvalue_fixture(path, count)
    validate_tsv_table(path, "population_eigenvalues.schema.json")


@pytest.mark.parametrize("count", ["0", "1", "-2", "2.0", "02", "", "abc"])
def test_population_eigenvalues_refuse_invalid_reference_counts(tmp_path, count):
    path = tmp_path / "eigenvalues.tsv"
    eigenvalue_fixture(path, count)
    with pytest.raises(TableValidationError):
        validate_tsv_table(path, "population_eigenvalues.schema.json")
