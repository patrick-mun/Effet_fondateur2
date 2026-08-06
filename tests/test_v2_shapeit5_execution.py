from effet_fondateur.phasing.execution import (
    PhasedVariant,
    _carrier_rows,
    _mendel_errors,
    _transmissions,
)


def _target(*genotypes: str, confidences: tuple[float | None, ...]) -> PhasedVariant:
    return PhasedVariant(
        "chr19",
        100_000,
        "target_GRCh38_1_100000_A_G",
        "A",
        "G",
        genotypes,
        confidences,
    )


def test_carrier_assignment_marks_low_singleton_confidence_unreliable() -> None:
    target = _target("0|1", "0|0", confidences=(0.7, None))
    explicit = (
        {"SAMPLE_ID": "carrier", "GENOTYPE": "A/G"},
        {"SAMPLE_ID": "noncarrier", "GENOTYPE": "A/A"},
    )

    rows = _carrier_rows(target, ["carrier", "noncarrier"], explicit, 0.9)

    assert rows[0]["CARRIER_HAPLOTYPE"] == "H2"
    assert rows[0]["CONFIDENCE_STATUS"] == "SCORED_LOW"
    assert rows[0]["RELIABILITY_STATUS"] == "UNRELIABLE"
    assert rows[1]["CONFIDENCE_STATUS"] == "NOT_APPLICABLE_HOMOZYGOUS"


def test_trio_transmission_identifies_direct_child_haplotype_orientation() -> None:
    target = _target("0|0", "1|1", "0|1", confidences=(None, None, 0.95))
    samples = ["father", "mother", "child"]
    pedigree = [("child", "father", "mother")]

    rows = _transmissions(target, samples, pedigree)

    assert rows[0]["TRANSMISSION_STATUS"] == "DIRECT"
    assert rows[0]["PATERNAL_CHILD_HAPLOTYPE"] == "H1"
    assert rows[0]["MATERNAL_CHILD_HAPLOTYPE"] == "H2"
    assert _mendel_errors([target], samples, pedigree) == 0


def test_mendel_error_is_detected_before_phasing() -> None:
    target = _target("0/0", "0/0", "1/1", confidences=(None, None, None))

    assert _mendel_errors(
        [target], ["father", "mother", "child"], [("child", "father", "mother")]
    ) == 1
