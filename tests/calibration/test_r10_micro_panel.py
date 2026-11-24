"""Tests ciblés du constructeur transversal R10 sans outil scientifique."""

from pathlib import Path

import numpy as np

from effet_fondateur.calibration.contracts import STAGE_IDS, load_yaml
from effet_fondateur.calibration.micro_contract import (
    CORE_STAGES,
    STAGE_DEPENDENCIES,
    STAGE_PRIMARY_OUTPUT,
    validate_micro_topology,
)
from effet_fondateur.calibration.micro_panel import (
    MicroPanelError,
    _build_haplotypes,
    _build_subjects,
    _build_variants,
    _observations,
    _reference_pools,
    _sandbox_profile,
    _seed128,
    _write_plink_dataset,
    _write_release_checksum_index,
    _write_vcf,
    _verify_release_checksum_index,
)
from effet_fondateur.calibration.micro_stage_worker import execute_stage


ROOT = Path(__file__).resolve().parents[2]
CONFIGURATION = ROOT / "calibration/configurations/micro_panel_v1.yaml"


def test_micro_topology_covers_all_stages_and_is_acyclic() -> None:
    validate_micro_topology()
    assert tuple(STAGE_DEPENDENCIES) == STAGE_IDS
    assert set(STAGE_PRIMARY_OUTPUT) == set(STAGE_IDS)
    assert "17" in CORE_STAGES
    assert STAGE_DEPENDENCIES["16C"] == ("12", "12A", "13", "16A")
    assert len(STAGE_DEPENDENCIES["17"]) == 16


def test_micro_subject_composition_is_exact() -> None:
    subjects = _build_subjects()
    counts: dict[str, int] = {}
    for subject in subjects:
        counts[subject.component_type] = counts.get(subject.component_type, 0) + 1
    assert len(subjects) == 60
    assert counts == {
        "TRIO": 12,
        "SIBSHIP_4": 16,
        "MULTIGENERATIONAL_5": 20,
        "UNRELATED": 12,
    }
    assert len({subject.master_id for subject in subjects}) == 60


def test_configured_target_catalog_and_22_autosomes_are_materialized() -> None:
    variants, target = _build_variants(load_yaml(CONFIGURATION))
    assert target["role"] == "STANDARD"
    assert target["variant_id"] == "SYN-TARGET-STANDARD"
    assert {variant.chromosome for variant in variants} == set(range(1, 23))
    assert any(variant.variant_id == target["variant_id"] for variant in variants)
    assert len(variants) > 100


def test_haplotypes_and_measurement_copy_are_deterministic() -> None:
    configuration = load_yaml(CONFIGURATION)
    subjects = _build_subjects()
    variants, target = _build_variants(configuration)
    first = _build_haplotypes(subjects, variants, target)
    second = _build_haplotypes(subjects, variants, target)
    assert np.array_equal(first, second)
    observed, events = _observations(first)
    assert first.shape == (60, 2, len(variants))
    assert observed.shape == (60, len(variants))
    assert len(events) == 3
    assert int(np.sum(observed == -1)) == 2


def test_seed_derivation_is_domain_separated() -> None:
    assert _seed128("RESOURCE_ONLY", "TARGET", "MASTER") == _seed128(
        "RESOURCE_ONLY", "TARGET", "MASTER"
    )
    assert _seed128("RESOURCE_ONLY", "TARGET", "MASTER") != _seed128(
        "RESOURCE_ONLY", "IBD", "MASTER"
    )


def test_six_reference_pools_are_disjoint_and_have_r07_counts() -> None:
    pools = _reference_pools()
    members = [member for pool in pools for member in pool["member_ids"]]
    by_type = {
        reference_type: sum(
            len(pool["member_ids"])
            for pool in pools
            if pool["reference_type"] == reference_type
        )
        for reference_type in {"PHASING", "ANCESTRY", "NULL"}
    }
    assert len(pools) == 6
    assert len(members) == len(set(members)) == 240
    assert by_type == {"PHASING": 90, "ANCESTRY": 90, "NULL": 60}


def test_plink_bed_is_variant_major_and_dimensioned(tmp_path: Path) -> None:
    configuration = load_yaml(CONFIGURATION)
    subjects = _build_subjects()
    variants, target = _build_variants(configuration)
    haplotypes = _build_haplotypes(subjects, variants, target)
    observed, _ = _observations(haplotypes)
    pseudonyms = [f"SYN-{index:03d}" for index in range(len(subjects))]
    bed, bim, fam = _write_plink_dataset(tmp_path, pseudonyms, variants, observed)
    assert bed.read_bytes()[:3] == bytes((0x6C, 0x1B, 0x01))
    assert bed.stat().st_size == 3 + len(variants) * 15
    assert len(bim.read_text().splitlines()) == len(variants)
    assert len(fam.read_text().splitlines()) == 60


def test_gzip_vcf_is_byte_reproducible(tmp_path: Path) -> None:
    configuration = load_yaml(CONFIGURATION)
    subjects = _build_subjects()
    variants, target = _build_variants(configuration)
    haplotypes = _build_haplotypes(subjects, variants, target)
    observed, _ = _observations(haplotypes)
    pseudonyms = [f"SYN-{index:03d}" for index in range(len(subjects))]
    first = tmp_path / "first.vcf.gz"
    second = tmp_path / "second.vcf.gz"
    _write_vcf(first, pseudonyms, variants, observed)
    _write_vcf(second, pseudonyms, variants, observed)
    assert first.read_bytes() == second.read_bytes()


def test_release_checksum_index_is_exhaustive(tmp_path: Path) -> None:
    (tmp_path / "governance").mkdir()
    (tmp_path / "a.txt").write_text("a\n")
    (tmp_path / "governance/b.json").write_text("{}\n")
    _write_release_checksum_index(tmp_path)
    _verify_release_checksum_index(tmp_path)
    (tmp_path / "a.txt").write_text("changed\n")
    try:
        _verify_release_checksum_index(tmp_path)
    except MicroPanelError as error:
        assert error.code == "release_checksum_mismatch"
    else:
        raise AssertionError("Une altération doit invalider l'index de release")


def test_worker_traverses_declared_dependency_without_truth(tmp_path: Path) -> None:
    release_root = tmp_path / "release"
    input_root = release_root / "blinded_inputs"
    run_root = release_root / "runs/RUN"
    input_root.mkdir(parents=True)
    run_root.mkdir(parents=True)
    (input_root / "inventory.json").write_text('{"synthetic": true}\n')
    execute_stage(
        release_root=release_root,
        input_root=Path("blinded_inputs"),
        run_root=Path("runs/RUN"),
        stage_id="00",
        release_id="release-test",
    )
    execute_stage(
        release_root=release_root,
        input_root=Path("blinded_inputs"),
        run_root=Path("runs/RUN"),
        stage_id="01",
        release_id="release-test",
    )
    assert (run_root / "00/stage_outputs.json").is_file()
    assert (run_root / "01/source_qc.tsv").is_file()
    assert "private_truth" not in (run_root / "01/stage_inputs.json").read_text()


def test_sandbox_profile_denies_private_governance_and_real_inputs(tmp_path: Path) -> None:
    repository = tmp_path / "repository"
    release = repository / "calibration_workspace/releases/release"
    profile = _sandbox_profile(repository, release)
    assert str(release / "private_truth") in profile
    assert str(release / "governance") in profile
    assert str(repository / "data/input") in profile
    assert "deny file-read*" in profile
    assert "deny file-write*" in profile
