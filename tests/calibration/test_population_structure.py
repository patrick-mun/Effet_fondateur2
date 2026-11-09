"""Tests unitaires des strates et mélanges D12."""

from __future__ import annotations

from collections import Counter
import csv
import io
import json

import pytest

from effet_fondateur.calibration.simulation.population_structure import (
    MixtureSpecification,
    PopulationSourceDiploid,
    PopulationStructureError,
    build_population_panel,
    canonical_ancestry_truth_bytes,
    parse_d11_population_provenance,
    parse_mixture_specifications,
    parse_observed_variant_keys,
    transform_population_vcf,
)


SOURCES = (
    PopulationSourceDiploid("SYN_CAL_000001", "POP_A"),
    PopulationSourceDiploid("SYN_CAL_000002", "POP_A"),
    PopulationSourceDiploid("SYN_CAL_000003", "POP_B"),
    PopulationSourceDiploid("SYN_CAL_000004", "POP_B"),
    PopulationSourceDiploid("SYN_CAL_000005", "POP_C"),
)


def _mixtures() -> tuple[MixtureSpecification, ...]:
    return (MixtureSpecification("ADMIX_A_B", "POP_A", "POP_B", 2),)


def _panel(seed: int = 1002):
    return build_population_panel(
        SOURCES,
        _mixtures(),
        random_seed=seed,
        synthetic_identifier_prefix="SYN_CAL_D12_",
    )


def test_d11_provenance_is_reduced_to_synthetic_population_assignments() -> None:
    provenance = io.StringIO(
        "synthetic_identifier\tsynthetic_haplotype_index\tsource_identifier\t"
        "source_haplotype_index\tsource_population\n"
        "SYN_CAL_000001\t1\tDONOR_A\t1\tPOP_A\n"
        "SYN_CAL_000001\t2\tDONOR_B\t2\tPOP_A\n"
    )
    assert parse_d11_population_provenance(provenance) == (
        PopulationSourceDiploid("SYN_CAL_000001", "POP_A"),
    )


def test_population_panel_is_deterministic_and_uses_each_haplotype_once() -> None:
    first = _panel(1002)
    assert first == _panel(1002)
    observed = Counter(
        (origin.sample_id, origin.haplotype_index)
        for diploid in first
        for origin in (diploid.haplotype_1, diploid.haplotype_2)
    )
    expected = Counter(
        (source.sample_id, index) for source in SOURCES for index in (1, 2)
    )
    assert observed == expected
    assert len(first) == len(SOURCES)


def test_reciprocal_pair_produces_two_exact_fifty_fifty_mixtures() -> None:
    admixed = [diploid for diploid in _panel() if diploid.stratum_id == "ADMIX_A_B"]
    assert len(admixed) == 2
    assert all(
        diploid.mixture_proportions == (("POP_A", 0.5), ("POP_B", 0.5))
        for diploid in admixed
    )
    assert sum(len(diploid.mixture_proportions) == 1 for diploid in _panel()) == 3


def test_configuration_rejects_unequal_mixture_before_d14() -> None:
    config = {
        "schema_version": "1.0.0",
        "population_structure_config_id": "d12.test",
        "config_state": "FINALIZED",
        "method": {
            "method_id": "reciprocal_whole_haplotype_mixture_v1",
            "method_version": "1.0.0",
            "random_seed": 1002,
            "pure_policy": "RETAIN_UNUSED_D11_DIPLOIDS",
            "mixture_unit": "RECIPROCAL_WHOLE_HAPLOTYPE_PAIR",
            "supported_proportions": "PURE_OR_TWO_WAY_EQUAL",
        },
        "mixtures": [
            {
                "mixture_id": "ADMIX_A_B",
                "individual_count": 2,
                "components": [
                    {"stratum_id": "POP_A", "proportion": 0.75},
                    {"stratum_id": "POP_B", "proportion": 0.25},
                ],
            }
        ],
    }
    with pytest.raises(
        PopulationStructureError,
        match="population_structure_mixture_proportion_unsupported",
    ):
        parse_mixture_specifications(config)


def test_ancestry_truth_has_canonical_json_proportions() -> None:
    rows = canonical_ancestry_truth_bytes(_panel()).decode().splitlines()
    assert rows[0] == "sample_id\tstratum_id\tmixture_proportions"
    records = csv.DictReader(rows, delimiter="\t")
    decoded = [json.loads(row["mixture_proportions"]) for row in records]
    assert all(abs(sum(proportions.values()) - 1.0) < 1e-12 for proportions in decoded)


def test_vcf_transformation_keeps_dense_haplotypes_and_d11b_marker_universe() -> None:
    source = io.StringIO(
        "##fileformat=VCFv4.2\n"
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t"
        + "\t".join(source.sample_id for source in SOURCES)
        + "\n"
        "chr19\t100\tv1\tA\tG\t.\tPASS\t.\tGT\t0|0\t0|1\t1|0\t1|1\t0|1\n"
        "chr19\t200\tv2\tC\tT\t.\tPASS\t.\tGT\t0|1\t1|0\t0|0\t1|1\t0|1\n"
    )
    observed_keys = parse_observed_variant_keys(
        io.StringIO(
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tX\n"
            "chr19\t200\tCHIP-2\tC\tT\t.\tPASS\t.\tGT\t0/1\n"
        )
    )
    dense = io.StringIO()
    observed = io.StringIO()
    stats = transform_population_vcf(
        source,
        dense,
        observed,
        panel=_panel(),
        observed_variant_keys=observed_keys,
        chromosome=19,
    )
    assert stats.dense_variant_count == 2
    assert stats.observed_variant_count == 1
    assert stats.allele_count_difference == 0
    assert "\t100\t" not in observed.getvalue()
    assert "\t200\t" in observed.getvalue()
    assert "\tCHIP-2\t" in observed.getvalue()
    assert "POP_A" not in observed.getvalue()
