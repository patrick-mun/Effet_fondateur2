"""Tests unitaires du réappariement et de la transformation D11."""

from __future__ import annotations

from collections import Counter
from copy import deepcopy
import io

import pytest

from effet_fondateur.calibration.simulation.base_haplotypes import (
    BaseHaplotypeError,
    HaplotypeOrigin,
    SyntheticDiploid,
    canonical_provenance_bytes,
    pair_donor_haplotypes,
    parse_genetic_map,
    parse_population_assignments,
    transform_phased_vcf,
    validate_base_haplotype_manifest,
    validate_pairings,
)


DONORS = ("DONOR_A", "DONOR_B", "DONOR_C", "DONOR_D")
POPULATIONS = {
    "DONOR_A": "POP_1",
    "DONOR_B": "POP_1",
    "DONOR_C": "POP_2",
    "DONOR_D": "POP_2",
}
SHA = "a" * 64


def _pairings(seed: int = 1002) -> tuple[SyntheticDiploid, ...]:
    return pair_donor_haplotypes(
        DONORS,
        POPULATIONS,
        random_seed=seed,
        synthetic_identifier_prefix="SYN_CAL_",
    )


def _map_points():
    return parse_genetic_map(
        [
            "CHROMOSOME\tPOSITION_BP\tPOSITION_CM\n",
            "19\t50\t0.0\n",
            "19\t200\t1.0\n",
            "19\t400\t3.0\n",
        ],
        expected_chromosome=19,
    )


def _vcf_lines() -> list[str]:
    return [
        "##fileformat=VCFv4.2\n",
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tDONOR_A\tDONOR_B\tDONOR_C\tDONOR_D\n",
        "chr19\t100\tv1\tA\tG\t.\tPASS\t.\tGT\t0|0\t0|1\t1|0\t1|1\n",
        "chr19\t200\tv2\tC\tT\t.\tPASS\t.\tGT\t0|1\t0|1\t1|0\t1|0\n",
        "chr19\t300\tv3\tG\tA\t.\tPASS\t.\tGT\t1|1\t0|0\t1|0\t0|1\n",
    ]


def _transform(pairings: tuple[SyntheticDiploid, ...] | None = None):
    phased = io.StringIO()
    unphased = io.StringIO()
    genetic_map = io.StringIO()
    stats = transform_phased_vcf(
        _vcf_lines(),
        phased,
        unphased,
        genetic_map,
        pairings=pairings or _pairings(),
        genetic_map_points=_map_points(),
        chromosome=19,
        start_bp=100,
        end_bp=300,
    )
    return stats, phased.getvalue(), unphased.getvalue(), genetic_map.getvalue()


def test_population_assignments_are_parsed() -> None:
    assignments = parse_population_assignments(
        [
            "FamilyID SampleID FatherID MotherID Population\n",
            "F1 DONOR_A 0 0 POP_1\n",
            "F2 DONOR_B 0 0 POP_1\n",
        ]
    )
    assert assignments == {"DONOR_A": "POP_1", "DONOR_B": "POP_1"}


def test_pairing_is_deterministic_for_same_seed() -> None:
    assert _pairings(1234) == _pairings(1234)


def test_pairing_changes_with_seed() -> None:
    assert _pairings(1234) != _pairings(4321)


def test_every_source_haplotype_is_used_exactly_once() -> None:
    observed = Counter(
        (origin.source_identifier, origin.source_haplotype_index)
        for pairing in _pairings()
        for origin in (pairing.haplotype_1, pairing.haplotype_2)
    )
    expected = Counter(
        (donor, haplotype_index)
        for donor in DONORS
        for haplotype_index in (1, 2)
    )
    assert observed == expected


def test_no_synthetic_pair_uses_the_same_source_donor() -> None:
    for pairing in _pairings():
        assert pairing.haplotype_1.source_identifier != pairing.haplotype_2.source_identifier


def test_pairs_stay_within_public_population_stratum() -> None:
    for pairing in _pairings():
        assert pairing.haplotype_1.source_population == pairing.source_population
        assert pairing.haplotype_2.source_population == pairing.source_population


def test_singleton_population_is_rejected() -> None:
    populations = dict(POPULATIONS)
    populations["DONOR_D"] = "SINGLETON"
    with pytest.raises(BaseHaplotypeError, match="population_donor_count_insufficient"):
        pair_donor_haplotypes(
            DONORS,
            populations,
            random_seed=1,
            synthetic_identifier_prefix="SYN_CAL_",
        )


def test_missing_population_is_rejected() -> None:
    populations = dict(POPULATIONS)
    del populations["DONOR_D"]
    with pytest.raises(BaseHaplotypeError, match="population_assignment_missing_for_donor"):
        pair_donor_haplotypes(
            DONORS,
            populations,
            random_seed=1,
            synthetic_identifier_prefix="SYN_CAL_",
        )


def test_same_source_pair_is_rejected_by_runtime_audit() -> None:
    pairings = list(_pairings())
    first = pairings[0]
    pairings[0] = SyntheticDiploid(
        first.synthetic_identifier,
        first.haplotype_1,
        HaplotypeOrigin(
            first.haplotype_1.source_identifier,
            2,
            first.haplotype_1.source_population,
        ),
        first.source_population,
    )
    with pytest.raises(BaseHaplotypeError, match="same_source_pair_forbidden"):
        validate_pairings(pairings, DONORS, POPULATIONS)


def test_private_provenance_is_canonical_and_complete() -> None:
    payload = canonical_provenance_bytes(_pairings())
    rows = payload.decode("utf-8").splitlines()
    assert rows[0].startswith("synthetic_identifier\t")
    assert len(rows) == 1 + 2 * len(DONORS)
    assert payload == canonical_provenance_bytes(tuple(reversed(_pairings())))


def test_genetic_map_requires_monotonic_bp_and_cm() -> None:
    with pytest.raises(BaseHaplotypeError, match="genetic_map_not_monotonic"):
        parse_genetic_map(
            [
                "CHROMOSOME\tPOSITION_BP\tPOSITION_CM\n",
                "19\t100\t1.0\n",
                "19\t200\t0.5\n",
            ],
            expected_chromosome=19,
        )


def test_vcf_transformation_has_exact_dimensions_and_concordance() -> None:
    stats, phased, unphased, genetic_map = _transform()
    assert stats.variant_count == 3
    assert stats.synthetic_individual_count == 4
    assert stats.source_haplotype_count == 8
    assert stats.synthetic_haplotype_count == 8
    assert stats.allele_count_difference == 0
    assert stats.genotype_haplotype_discordance_count == 0
    assert sum(line.startswith("chr19\t") for line in phased.splitlines()) == 3
    assert sum(line.startswith("chr19\t") for line in unphased.splitlines()) == 3
    assert genetic_map.splitlines() == [
        "CHROMOSOME\tPOSITION_BP\tPOSITION_CM",
        "19\t100\t0.333333333333",
        "19\t200\t1",
        "19\t300\t2",
    ]


def test_distinct_biallelic_variants_may_share_a_position() -> None:
    lines = _vcf_lines()
    lines.insert(
        3,
        "chr19\t100\tv1b\tA\tC\t.\tPASS\t.\tGT\t0|1\t0|0\t1|0\t1|1\n",
    )
    phased = io.StringIO()
    stats = transform_phased_vcf(
        lines,
        phased,
        io.StringIO(),
        io.StringIO(),
        pairings=_pairings(),
        genetic_map_points=_map_points(),
        chromosome=19,
        start_bp=100,
        end_bp=300,
    )
    assert stats.variant_count == 4
    assert phased.getvalue().count("chr19\t100\t") == 2


def test_vcf_output_contains_only_synthetic_sample_names() -> None:
    _, phased, _, _ = _transform()
    header = next(line for line in phased.splitlines() if line.startswith("#CHROM"))
    assert all(donor not in header for donor in DONORS)
    assert all(pairing.synthetic_identifier in header for pairing in _pairings())


def test_haplotype_sequences_are_preserved_as_a_multiset() -> None:
    _, phased, _, _ = _transform()
    source_records = [line.rstrip().split("\t") for line in _vcf_lines() if line.startswith("chr19")]
    source_sequences: list[tuple[str, ...]] = []
    for sample_index in range(4):
        genotypes = [record[9 + sample_index] for record in source_records]
        for haplotype_index in (0, 1):
            source_sequences.append(
                tuple(genotype.split("|")[haplotype_index] for genotype in genotypes)
            )
    output_records = [line.split("\t") for line in phased.splitlines() if line.startswith("chr19")]
    output_sequences: list[tuple[str, ...]] = []
    for sample_index in range(4):
        genotypes = [record[9 + sample_index] for record in output_records]
        for haplotype_index in (0, 1):
            output_sequences.append(
                tuple(genotype.split("|")[haplotype_index] for genotype in genotypes)
            )
    assert Counter(output_sequences) == Counter(source_sequences)


@pytest.mark.parametrize(
    ("replacement", "error"),
    [
        ("0/1", "vcf_genotype_not_phased_diploid"),
        ("0|.", "vcf_genotype_not_complete_biallelic"),
        ("0|2", "vcf_genotype_not_complete_biallelic"),
    ],
)
def test_invalid_source_genotypes_are_rejected(replacement: str, error: str) -> None:
    lines = _vcf_lines()
    lines[2] = lines[2].replace("0|0", replacement, 1)
    with pytest.raises(BaseHaplotypeError, match=error):
        transform_phased_vcf(
            lines,
            io.StringIO(),
            io.StringIO(),
            io.StringIO(),
            pairings=_pairings(),
            genetic_map_points=_map_points(),
            chromosome=19,
            start_bp=100,
            end_bp=300,
        )


def test_map_extrapolation_is_rejected() -> None:
    with pytest.raises(BaseHaplotypeError, match="genetic_map_domain_not_covered"):
        transform_phased_vcf(
            _vcf_lines(),
            io.StringIO(),
            io.StringIO(),
            io.StringIO(),
            pairings=_pairings(),
            genetic_map_points=_map_points(),
            chromosome=19,
            start_bp=1,
            end_bp=300,
        )


def _manifest(stats) -> dict[str, object]:
    return {
        "schema_version": "1.0.0",
        "base_haplotype_manifest_id": "d11.manifest",
        "manifest_state": "FINALIZED",
        "panel": {
            "role": "CALIBRATION",
            "individual_count": 4,
            "haplotype_count": 8,
        },
        "source": {
            "vcf_sha256": SHA,
            "index_sha256": SHA,
            "genetic_map_sha256": SHA,
        },
        "partition": {
            "partition_manifest_sha256": SHA,
            "membership_sha256": SHA,
        },
        "method": {
            "method_id": "population_stratified_whole_haplotype_repairing_v1",
            "method_version": "1.0.0",
        },
        "confidentiality": {
            "raw_donor_identifiers_in_manifest": False,
            "local_source_paths_in_manifest": False,
            "donor_provenance_private": True,
            "artifacts_version_control_allowed": False,
            "aggregate_counts_publishable": True,
        },
        "artifacts": [
            {"role": "PHASED_HAPLOTYPES"},
            {"role": "UNPHASED_GENOTYPES"},
            {"role": "DONOR_PROVENANCE"},
            {"role": "GENETIC_MAP"},
        ],
        "integrity": {
            "source_donor_count": 4,
            "synthetic_individual_count": 4,
            "source_haplotype_count": 8,
            "synthetic_haplotype_count": 8,
            "variant_count": stats.variant_count,
            "population_stratum_count": 2,
            "haplotype_copy_error_count": 0,
            "same_source_pair_count": 0,
            "allele_count_difference": 0,
            "genotype_haplotype_discordance_count": 0,
            "genetic_map_coordinate_error_count": 0,
            "genetic_map_domain_covered": True,
            "mendelian_evaluation": "DEFERRED_TO_D13_NO_PEDIGREE",
            "checks": [
                "SOURCE_HASHES_MATCHED",
                "PARTITION_BINDING_MATCHED",
                "HAPLOTYPE_COPY_BIJECTION",
                "SOURCE_PAIRS_DISJOINT",
                "ALLELE_FREQUENCIES_EXACT",
                "HAPLOTYPE_LD_EXACT_BY_PERMUTATION",
                "GENOTYPE_HAPLOTYPE_CONCORDANCE_100_PERCENT",
                "GENETIC_MAP_COORDINATES_VALID",
                "RAW_DONOR_IDENTIFIERS_NOT_IN_MANIFEST",
            ],
        },
    }


def test_runtime_manifest_audit_accepts_aggregate_contract() -> None:
    stats, _, _, _ = _transform()
    audit = validate_base_haplotype_manifest(
        _manifest(stats), pairings=_pairings(), transformation_stats=stats
    )
    assert audit.individual_count == 4
    assert audit.variant_count == 3


def test_runtime_manifest_audit_rejects_raw_donor_identifier() -> None:
    stats, _, _, _ = _transform()
    manifest = _manifest(stats)
    manifest["unexpected_note"] = f"derived from {DONORS[0]}"
    with pytest.raises(BaseHaplotypeError, match="raw_donor_identifier_in_manifest"):
        validate_base_haplotype_manifest(
            manifest, pairings=_pairings(), transformation_stats=stats
        )


def test_runtime_manifest_audit_rejects_nonzero_discordance() -> None:
    stats, _, _, _ = _transform()
    manifest = deepcopy(_manifest(stats))
    manifest["integrity"]["genotype_haplotype_discordance_count"] = 1
    with pytest.raises(BaseHaplotypeError, match="base_haplotype_manifest_integrity_mismatch"):
        validate_base_haplotype_manifest(
            manifest, pairings=_pairings(), transformation_stats=stats
        )
