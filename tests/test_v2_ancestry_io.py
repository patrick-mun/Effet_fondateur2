from pathlib import Path

import numpy as np
import pytest

from effet_fondateur.ancestry import (
    AncestryAnalysisError,
    parse_vcf_query_panel,
    read_bim_variants,
    read_plink_raw_panel,
)
from effet_fondateur.ancestry.io import _dosage_column_indexes


class _CountingHeader(list[str]):
    """Compte les parcours complets pour protéger la complexité du parseur."""

    def __init__(self, values: list[str]) -> None:
        super().__init__(values)
        self.iteration_count = 0

    def __iter__(self):  # type: ignore[no-untyped-def]
        self.iteration_count += 1
        return super().__iter__()


def test_plink_raw_is_mapped_to_master_ids_and_bim_a1(tmp_path: Path) -> None:
    bim = tmp_path / "panel.bim"
    raw = tmp_path / "panel.raw"
    bim.write_text("2 v1 0 100 G A\n2 v2 0 200 T C\n", encoding="utf-8")
    raw.write_text(
        "FID IID PAT MAT SEX PHENOTYPE v1_G v2_T\nF I 0 0 1 -9 1 NA\n",
        encoding="utf-8",
    )
    variants = read_bim_variants(bim)
    panel = read_plink_raw_panel(raw, variants, {("F", "I"): "STUDY"})
    assert (variants[0].ref, variants[0].alt) == ("A", "G")
    assert panel.sample_ids == ("STUDY",)
    assert panel.alt_dosages[0, 0] == 1
    assert np.isnan(panel.alt_dosages[0, 1])


def test_plink_raw_header_is_indexed_once_for_large_panels(tmp_path: Path) -> None:
    bim = tmp_path / "panel.bim"
    bim.write_text(
        "".join(f"1 v{index} 0 {index + 1} G A\n" for index in range(10_000)),
        encoding="utf-8",
    )
    variants = read_bim_variants(bim)
    header = _CountingHeader(
        ["FID", "IID", "PAT", "MAT", "SEX", "PHENOTYPE"]
        + [f"v{index}_G" for index in range(10_000)]
    )

    indexes = _dosage_column_indexes(header, variants)

    assert indexes == list(range(6, 10_006))
    assert header.iteration_count == 1


def test_plink_raw_rejects_duplicate_expected_dosage_column(tmp_path: Path) -> None:
    bim = tmp_path / "panel.bim"
    bim.write_text("1 v1 0 1 G A\n", encoding="utf-8")
    variants = read_bim_variants(bim)

    with pytest.raises(AncestryAnalysisError, match="missing_or_ambiguous"):
        _dosage_column_indexes(
            ["FID", "IID", "PAT", "MAT", "SEX", "PHENOTYPE", "v1_G", "v1_G"],
            variants,
        )


def test_vcf_query_decomposes_any_configured_local_haplotypes() -> None:
    panel = parse_vcf_query_panel(
        ["chr11\t900\tlocal1\tC\tA\t1|0\t0|1\n"],
        ("S1", "S2"),
        haplotypes=True,
        phased_required=True,
    )
    assert panel.sample_ids == ("S1:H1", "S1:H2", "S2:H1", "S2:H2")
    np.testing.assert_array_equal(panel.alt_dosages[:, 0], [1, 0, 0, 1])


def test_unphased_local_query_blocks() -> None:
    with pytest.raises(AncestryAnalysisError, match="not_phased"):
        parse_vcf_query_panel(
            ["chr4\t10\tv\tA\tG\t0/1\n"],
            ("S",),
            haplotypes=True,
            phased_required=True,
        )
