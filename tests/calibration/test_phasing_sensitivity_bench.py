"""Contrats D23b sans outil génétique externe."""
import json
from pathlib import Path

import numpy as np
import pytest

from effet_fondateur.calibration.phasing_sensitivity_bench import _schemas_dir, density_profile, load_configuration
from effet_fondateur.calibration.simulation.measurement_errors import MeasurementPanel


def _panel() -> MeasurementPanel:
    variants = tuple(
        (f"chr{chromosome}", str(position), f"m{chromosome}_{position}", "A", "G", ".", "PASS", ".", "GT")
        for chromosome in range(1, 23) for position in range(1, 9)
    )
    alleles = np.zeros((len(variants), 2, 2), dtype=np.int8)
    return MeasurementPanel(("S1", "S2"), variants, alleles)


def test_density_profiles_are_nested_by_chromosome_rank_and_keep_target() -> None:
    panel = _panel()
    target = {"chromosome": 5, "position_bp": 2}
    definitions = [
        {"profile_id": "DENSITY_25", "rank_residues": [0], "rank_modulus": 4},
        {"profile_id": "DENSITY_50", "rank_residues": [0, 2], "rank_modulus": 4},
        {"profile_id": "DENSITY_75", "rank_residues": [0, 1, 2], "rank_modulus": 4},
        {"profile_id": "DENSITY_100", "rank_residues": [0, 1, 2, 3], "rank_modulus": 4},
    ]
    profiles = [density_profile(panel, definition, target) for definition in definitions]
    loci = [{row[:2] for row in profile.variants} for profile in profiles]
    assert loci[0] < loci[1] < loci[2] < loci[3]
    assert [len(profile.variants) for profile in profiles] == [45, 89, 133, 176]
    assert all(("chr5", "2") in values for values in loci)
    assert {row[1] for row in profiles[0].variants if row[0] == "chr5"} == {"1", "2", "6"}


def test_sensitivity_configuration_is_closed(tmp_path: Path) -> None:
    source = Path(__file__).resolve().parents[2] / "calibration/schemas/phasing_sensitivity_config.schema.json"
    config = {
        "schema_version": "1.0.0",
        "method_id": "paired_density_pedigree_phasing_bench_v1",
        "preparation_manifest_sha256": "0" * 64,
        "density_profiles": [
            {"profile_id": "DENSITY_25", "rank_residues": [0], "rank_modulus": 4},
            {"profile_id": "DENSITY_50", "rank_residues": [0, 2], "rank_modulus": 4},
            {"profile_id": "DENSITY_75", "rank_residues": [0, 1, 2], "rank_modulus": 4},
            {"profile_id": "DENSITY_100", "rank_residues": [0, 1, 2, 3], "rank_modulus": 4},
        ],
        "regional_pedigree_modes": ["WITH_PEDIGREE", "WITHOUT_PEDIGREE"],
        "rare_maf_upper_bound": 0.05,
        "map_absolute_tolerance_cm": 1e-8,
        "phasing_seed": 2303,
        "bootstrap_seed": 2311,
        "bootstrap_replicates": 2000,
        "density_gap_bp": 100000,
        "target_flank_bp": 5000000,
        "command_timeout_seconds": 1800,
        "minimum_total_variants": 10000,
        "minimum_variants_per_autosome": 100,
    }
    assert source.is_file()
    path = tmp_path / "config.json"
    path.write_text(json.dumps(config))
    assert load_configuration(path) == config
    path.write_text(json.dumps(config | {"unexpected": True}))
    with pytest.raises(ValueError):
        load_configuration(path)


def test_sensitivity_schemas_resolve_from_calibration_directory() -> None:
    assert (_schemas_dir() / "phasing_sensitivity_config.schema.json").is_file()
    assert (_schemas_dir() / "phasing_sensitivity_run.schema.json").is_file()
