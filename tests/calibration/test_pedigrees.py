"""Tests unitaires des pedigrees et transmissions D13."""

from __future__ import annotations

from collections import Counter
import io

import pytest

from effet_fondateur.calibration.simulation.pedigrees import (
    PedigreeFamilySpecification,
    PedigreeSimulationError,
    SourceIndividual,
    build_pedigree_panel,
    canonical_transmission_bytes,
    parse_d12_source_truth,
    parse_pedigree_configuration,
    transform_pedigree_vcf,
)


SOURCES = tuple(
    SourceIndividual(
        f"SYN_CAL_D12_{index:06d}",
        "POP_A",
        (("POP_A", 1.0),),
        "POP_A",
        "POP_A",
    )
    for index in range(1, 9)
)
SPECIFICATIONS = (PedigreeFamilySpecification("PED_A", "POP_A", 1),)


def _panel(seed: int = 1003):
    return build_pedigree_panel(
        SOURCES,
        SPECIFICATIONS,
        unrelated_pair_count=1,
        random_seed=seed,
        synthetic_identifier_prefix="SYN_CAL_D13_",
    )


def test_d12_ancestry_and_haplotype_provenance_are_bound() -> None:
    records = parse_d12_source_truth(
        io.StringIO(
            "sample_id\tstratum_id\tmixture_proportions\n"
            'SYN_CAL_D12_000001\tADMIX_A_B\t"{""A"":0.5,""B"":0.5}"\n'
        ),
        io.StringIO(
            "sample_id\toutput_haplotype_index\td11_sample_id\t"
            "d11_haplotype_index\tsource_stratum_id\n"
            "SYN_CAL_D12_000001\t1\tSYN_CAL_D11_000001\t1\tA\n"
            "SYN_CAL_D12_000001\t2\tSYN_CAL_D11_000002\t2\tB\n"
        ),
    )
    assert records == (
        SourceIndividual(
            "SYN_CAL_D12_000001",
            "ADMIX_A_B",
            (("A", 0.5), ("B", 0.5)),
            "A",
            "B",
        ),
    )


def test_descendant_ancestry_follows_transmitted_haplotypes() -> None:
    sources = tuple(
        SourceIndividual(
            f"SYN_CAL_D12_{index:06d}",
            "ADMIX_A_B",
            (("A", 0.5), ("B", 0.5)),
            "A",
            "B",
        )
        for index in range(1, 9)
    )
    panel = build_pedigree_panel(
        sources,
        (PedigreeFamilySpecification("PED_ADMIX", "ADMIX_A_B", 1),),
        unrelated_pair_count=1,
        random_seed=1003,
        synthetic_identifier_prefix="SYN_CAL_D13_",
    )
    descendants = [
        individual for individual in panel.individuals if not individual.is_founder
    ]
    source_by_id = {source.sample_id: source for source in sources}
    for descendant in descendants:
        haplotype_strata = []
        for origin in (descendant.haplotype_1, descendant.haplotype_2):
            source = source_by_id[origin.d12_sample_id]
            haplotype_strata.append(
                (
                    source.haplotype_1_stratum_id,
                    source.haplotype_2_stratum_id,
                )[origin.d12_haplotype_index - 1]
            )
        expected = Counter(haplotype_strata)
        assert descendant.mixture_proportions == tuple(
            sorted((stratum, count / 2) for stratum, count in expected.items())
        )
        assert descendant.stratum_id == (
            haplotype_strata[0]
            if len(expected) == 1
            else "ADMIX_" + "_".join(sorted(expected))
        )


def test_first_cousin_template_is_deterministic_and_complete() -> None:
    first = _panel()
    assert first == _panel()
    assert first.source_individual_count == 8
    assert len(first.individuals) == 12
    assert len(first.transmissions) == 8
    assert first.family_count == 1
    assert first.trio_count == 4


def test_relationship_classes_and_coefficients_are_derived_from_pedigree() -> None:
    counts = Counter(record.relationship_class for record in _panel().kinship_records)
    assert counts == {
        "UNRELATED": 13,
        "FIRST_DEGREE": 9,
        "SECOND_DEGREE": 6,
        "THIRD_DEGREE": 1,
    }
    coefficients = {
        record.relationship_class: record.kinship_coefficient
        for record in _panel().kinship_records
    }
    assert coefficients == {
        "UNRELATED": 0.0,
        "FIRST_DEGREE": 0.25,
        "SECOND_DEGREE": 0.125,
        "THIRD_DEGREE": 0.0625,
    }


def test_every_descendant_copy_matches_the_declared_parental_copy() -> None:
    panel = _panel()
    by_id = {individual.sample_id: individual for individual in panel.individuals}
    for record in panel.transmissions:
        child = by_id[record.child_id]
        parent = by_id[record.parent_id]
        child_origin = (child.haplotype_1, child.haplotype_2)[
            record.child_haplotype_index - 1
        ]
        parent_origin = (parent.haplotype_1, parent.haplotype_2)[
            record.parental_haplotype_index - 1
        ]
        assert child_origin == parent_origin


def test_transmission_truth_explicitly_records_no_recombination() -> None:
    rows = canonical_transmission_bytes(_panel()).decode().splitlines()
    assert rows[0].endswith("recombination_count")
    assert all(row.endswith("\t0") for row in rows[1:])


def test_configuration_rejects_recombination_before_d14() -> None:
    config = {
        "schema_version": "1.0.0",
        "pedigree_config_id": "d13.test",
        "config_state": "FINALIZED",
        "method": {
            "method_id": "first_cousin_whole_haplotype_transmission_v1",
            "method_version": "1.0.0",
            "random_seed": 1003,
            "family_template": "FIRST_COUSIN_PEDIGREE",
            "source_selection": "SAME_STRATUM_WITHOUT_REPLACEMENT",
            "recombination_policy": "ENABLED",
        },
        "family_groups": [
            {
                "family_group_id": "PED_A",
                "stratum_id": "POP_A",
                "family_count": 1,
            }
        ],
        "unrelated_pair_count": 1,
    }
    with pytest.raises(PedigreeSimulationError, match="pedigree_method_invalid"):
        parse_pedigree_configuration(config)


def test_vcf_descendants_are_mendelian_and_observation_keeps_chip_id() -> None:
    source = io.StringIO(
        "##fileformat=VCFv4.2\n"
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t"
        + "\t".join(source.sample_id for source in SOURCES)
        + "\n"
        "chr19\t100\tv1\tA\tG\t.\tPASS\t.\tGT\t"
        "0|0\t0|1\t1|0\t1|1\t0|1\t1|0\t0|0\t1|1\n"
        "chr19\t200\tv2\tC\tT\t.\tPASS\t.\tGT\t"
        "0|1\t1|0\t0|0\t1|1\t0|1\t1|0\t0|0\t1|1\n"
    )
    dense = io.StringIO()
    observed = io.StringIO()
    stats = transform_pedigree_vcf(
        source,
        dense,
        observed,
        panel=_panel(),
        observed_variant_keys={("chr19", 200, "C", "T"): "CHIP-2"},
        chromosome=19,
    )
    assert stats.dense_variant_count == 2
    assert stats.observed_variant_count == 1
    assert stats.output_individual_count == 12
    assert stats.founder_haplotype_copy_error_count == 0
    assert stats.founder_allele_count_difference == 0
    assert stats.mendelian_error_count == 0
    assert "\tCHIP-2\t" in observed.getvalue()
    assert "POP_A" not in observed.getvalue()
