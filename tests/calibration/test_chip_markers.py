"""Tests unitaires de l'harmonisation puce–1000G D11b."""

from __future__ import annotations

import io

import pytest

from effet_fondateur.calibration.simulation.chip_markers import (
    ChipMarker,
    ChipMarkerError,
    ReferenceVariant,
    filter_vcf_to_chip_markers,
    harmonize_chip_markers,
    parse_chip_marker_manifest,
    parse_reference_vcf_variants,
    write_marker_map,
)


HEADER = (
    "MARKER_ID\tARRAY_DESIGN\tDESIGN_VERSION\tASSEMBLY\tCHROMOSOME\t"
    "POSITION_BP\tALLELE_A\tALLELE_B\tSTRAND\n"
)


def _marker_line(
    marker_id: str,
    position_bp: int,
    allele_a: str,
    allele_b: str,
    strand: str = "+",
) -> str:
    return (
        f"{marker_id}\tCytoScan750K_Accel\tNA36.r1\tGRCh38\t19\t"
        f"{position_bp}\t{allele_a}\t{allele_b}\t{strand}\n"
    )


def _variant(position_bp: int, reference: str, alternate: str) -> ReferenceVariant:
    return ReferenceVariant(19, position_bp, f"v{position_bp}", reference, alternate)


def test_official_normalized_marker_manifest_is_parsed() -> None:
    markers = parse_chip_marker_manifest(
        [HEADER, _marker_line("AX-1", 100, "A", "G")]
    )
    assert markers == (ChipMarker("AX-1", 19, 100, "A", "G", "+"),)


@pytest.mark.parametrize(
    ("field", "replacement", "error_code"),
    [
        ("CytoScan750K_Accel", "CytoScan750K_Array", "chip_marker_design_mismatch"),
        ("NA36.r1", "NA36.r7", "chip_marker_design_version_mismatch"),
        ("GRCh38", "GRCh37", "chip_marker_assembly_mismatch"),
    ],
)
def test_marker_manifest_refuses_wrong_design_version_or_assembly(
    field: str,
    replacement: str,
    error_code: str,
) -> None:
    with pytest.raises(ChipMarkerError) as captured:
        parse_chip_marker_manifest(
            [HEADER, _marker_line("AX-1", 100, "A", "G").replace(field, replacement)]
        )
    assert captured.value.code == error_code


def test_reference_vcf_parser_keeps_only_biallelic_snvs_in_domain() -> None:
    variants = parse_reference_vcf_variants(
        [
            "##fileformat=VCFv4.2\n",
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n",
            "chr19\t100\tv1\tA\tG\t.\tPASS\t.\n",
            "chr19\t150\tv2\tA\tAT\t.\tPASS\t.\n",
            "chr19\t300\tv3\tC\tT\t.\tPASS\t.\n",
        ],
        chromosome=19,
        start_bp=90,
        end_bp=200,
    )
    assert variants == (ReferenceVariant(19, 100, "v1", "A", "G"),)


def test_harmonization_accepts_direct_swapped_and_complement_relations() -> None:
    markers = (
        ChipMarker("AX-1", 19, 100, "A", "G", "+"),
        ChipMarker("AX-2", 19, 200, "G", "A", "+"),
        ChipMarker("AX-3", 19, 300, "T", "C", "-"),
        ChipMarker("AX-4", 19, 400, "C", "T", "-"),
    )
    variants = tuple(_variant(position, "A", "G") for position in (100, 200, 300, 400))
    matches, audit = harmonize_chip_markers(
        markers,
        variants,
        chromosome=19,
        start_bp=100,
        end_bp=400,
    )
    assert [match.allele_relation for match in matches] == [
        "DIRECT",
        "SWAPPED",
        "COMPLEMENT",
        "COMPLEMENT_SWAPPED",
    ]
    assert audit.matched_marker_count == 4
    assert audit.excluded_marker_count == 0


def test_harmonization_excludes_every_ambiguous_case_with_reason() -> None:
    markers = (
        ChipMarker("DUP-1", 19, 100, "A", "G", "+"),
        ChipMarker("DUP-2", 19, 100, "A", "G", "+"),
        ChipMarker("PAL", 19, 200, "A", "T", "UNKNOWN"),
        ChipMarker("MISSING", 19, 300, "A", "G", "+"),
        ChipMarker("MULTI", 19, 400, "A", "G", "+"),
        ChipMarker("ALLELES", 19, 500, "A", "C", "+"),
    )
    variants = (
        _variant(100, "A", "G"),
        _variant(200, "A", "T"),
        _variant(400, "A", "G"),
        ReferenceVariant(19, 400, "second", "A", "C"),
        _variant(500, "A", "G"),
    )
    matches, audit = harmonize_chip_markers(
        markers,
        variants,
        chromosome=19,
        start_bp=100,
        end_bp=500,
    )
    assert matches == ()
    assert dict(audit.exclusion_counts) == {
        "ALLELES_INCOMPATIBLE_OR_AMBIGUOUS": 1,
        "DUPLICATE_CHIP_POSITION": 2,
        "NO_REFERENCE_VARIANT": 1,
        "PALINDROMIC_STRAND_AMBIGUOUS": 1,
        "REFERENCE_POSITION_AMBIGUOUS": 1,
    }


def test_marker_map_and_observed_vcf_preserve_reference_orientation() -> None:
    markers = (ChipMarker("AX-1", 19, 100, "G", "A", "+"),)
    variants = (_variant(100, "A", "G"),)
    matches, _audit = harmonize_chip_markers(
        markers,
        variants,
        chromosome=19,
        start_bp=100,
        end_bp=100,
    )
    marker_map = io.StringIO()
    assert write_marker_map(matches, marker_map) == 1
    assert "AX-1\t19\t100\tA\tG\tSWAPPED" in marker_map.getvalue()

    observed = io.StringIO()
    written = filter_vcf_to_chip_markers(
        [
            "##fileformat=VCFv4.2\n",
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSYN\n",
            "chr19\t100\tv1\tA\tG\t.\tPASS\t.\tGT\t0/1\n",
            "chr19\t200\tv2\tC\tT\t.\tPASS\t.\tGT\t0/0\n",
        ],
        observed,
        matches,
    )
    assert written == 1
    assert "chr19\t100\tAX-1\tA\tG" in observed.getvalue()
    assert "\t200\t" not in observed.getvalue()
