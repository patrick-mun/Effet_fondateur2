import pytest

from effet_fondateur.phasing.execution import (
    PhasedVariant,
    Shapeit5ExecutionBlockError,
    _build_missing_genotype_mask,
    _carrier_rows,
    _mendel_check,
    _mendel_errors,
    _remask_vcf_gt_text,
    _same_unphased_genotype,
    _transmissions,
    _verify_missing_mask_restored,
    _verify_observed_genotypes_preserved,
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


def test_missing_pedigree_genotype_is_not_evaluated() -> None:
    target = _target("0/0", "1/1", "./.", confidences=(None, None, None))

    summary = _mendel_check(
        [target], ["father", "mother", "child"], [("child", "father", "mother")]
    )

    assert summary.error_count == 0
    assert summary.evaluable_record_count == 0
    assert summary.not_evaluated_record_count == 1
    assert _same_unphased_genotype("./.", ".|.")


def test_invalid_phasing_allele_remains_blocking() -> None:
    target = _target("0/0", "1/1", "0/2", confidences=(None, None, None))

    with pytest.raises(Shapeit5ExecutionBlockError, match="unsupported_phasing_genotype"):
        _mendel_check(
            [target],
            ["father", "mother", "child"],
            [("child", "father", "mother")],
        )


def test_missing_gt_completed_then_remasked_without_changing_other_format() -> None:
    input_variant = _target("./.", "0/1", confidences=(None, None))
    phased_variant = _target("1|0", "1|0", confidences=(0.8, 0.9))
    samples = ["missing", "observed"]
    mask = _build_missing_genotype_mask([input_variant], samples)

    assert _verify_observed_genotypes_preserved(
        [input_variant], [phased_variant], samples
    ) == 1
    vcf = (
        "##fileformat=VCFv4.2\n"
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tmissing\tobserved\n"
        "chr19\t100000\ttarget_GRCh38_1_100000_A_G\tA\tG\t.\tPASS\t.\tGT:DP:AD\t1|0:12:5,7\t1|0:9:4,5\n"
    )

    remasked, count = _remask_vcf_gt_text(vcf, samples, mask)

    assert count == 1
    assert ".|.:12:5,7" in remasked
    assert "1|0:9:4,5" in remasked
    restored_variant = _target(".|.", "1|0", confidences=(None, None))
    _verify_missing_mask_restored([input_variant], [restored_variant], samples)


def test_observed_gt_phase_separator_change_is_preserved() -> None:
    assert _verify_observed_genotypes_preserved(
        [_target("0/1", confidences=(None,))],
        [_target("1|0", confidences=(0.95,))],
        ["sample"],
    ) == 0


def test_observed_gt_allele_change_blocks() -> None:
    with pytest.raises(
        Shapeit5ExecutionBlockError, match="shapeit5_observed_genotype_modified"
    ):
        _verify_observed_genotypes_preserved(
            [_target("0/1", confidences=(None,))],
            [_target("1|1", confidences=(0.95,))],
            ["sample"],
        )


@pytest.mark.parametrize("missing_coordinate", ["absent_sample", "absent_variant"])
def test_partial_mask_restoration_blocks(missing_coordinate: str) -> None:
    key = ("chr19", 100_000, "target_GRCh38_1_100000_A_G", "A", "G")
    sample = "absent" if missing_coordinate == "absent_sample" else "sample"
    if missing_coordinate == "absent_variant":
        key = ("chr19", 999, "absent", "A", "G")
    mask = frozenset({(key, sample)})
    vcf = (
        "##fileformat=VCFv4.2\n"
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tsample\n"
        "chr19\t100000\ttarget_GRCh38_1_100000_A_G\tA\tG\t.\tPASS\t.\tGT\t0|1\n"
    )

    with pytest.raises(Shapeit5ExecutionBlockError):
        _remask_vcf_gt_text(vcf, ["sample"], mask)


def test_missing_target_gt_remains_forbidden() -> None:
    with pytest.raises(
        Shapeit5ExecutionBlockError, match="target_phasing_genotype_missing"
    ):
        _carrier_rows(
            _target(".|.", confidences=(None,)),
            ["sample"],
            [{"SAMPLE_ID": "sample", "GENOTYPE": "A/G"}],
            0.9,
        )
