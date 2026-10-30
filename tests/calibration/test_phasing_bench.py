"""Contrôles D23 sans outil externe ni donnée réelle."""
from pathlib import Path

import numpy as np
import pytest

from effet_fondateur.calibration.phasing_bench import fixture_inputs, load_configuration, observed_profile, write_observation
from effet_fondateur.calibration.simulation.measurement_errors import MeasurementPanel, read_panel


def panel() -> MeasurementPanel:
    variants = tuple((f"chr{chromosome}", str(position), f"m{chromosome}_{position}", "A", "G", ".", "PASS", ".", "GT")
                     for chromosome in range(1, 23) for position in range(1, 7))
    return MeasurementPanel(("S1", "S2"), variants, np.tile(np.array([[[0, 1], [1, 1]]], dtype=np.int8), (len(variants), 1, 1)))


def test_profiles_are_paired_reproducible_and_target_protected(tmp_path: Path) -> None:
    original = panel()
    config = dict(missingness_seed=2304, missingness_fraction=.5)
    target = dict(chromosome=5, position_bp=2)
    sparse = observed_profile(original, "DENSITY_50", config, target)
    assert len(sparse.variants) == 67
    assert any(row[:2] == ("chr5", "2") for row in sparse.variants)
    missing = observed_profile(original, "MISSING_1PCT", config, target)
    np.testing.assert_array_equal(missing.alleles, observed_profile(original, "MISSING_1PCT", config, target).alleles)
    assert np.any(missing.alleles < 0)
    assert np.all(original.alleles >= 0)
    index = next(index for index, row in enumerate(original.variants) if row[:2] == ("chr5", "2"))
    np.testing.assert_array_equal(missing.alleles[index], original.alleles[index])
    path = tmp_path / "observed.vcf"
    write_observation(path, missing)
    with path.open() as handle:
        recovered = read_panel(handle, phased=False)
    np.testing.assert_array_equal(recovered.alleles, missing.alleles)
    assert "|" not in path.read_text()


def test_configuration_closed_and_seeds_distinct(tmp_path: Path) -> None:
    import json
    source = Path(__file__).resolve().parents[2] / "calibration/configurations/d23_phasing_bench_v1.json"
    config = load_configuration(source)
    for key, value in (("phasing_seed", config["missingness_seed"]), ("unexpected", True), ("bootstrap_replicates", True)):
        altered = dict(config, **{key: value})
        path = tmp_path / "config.json"
        path.write_text(json.dumps(altered))
        with pytest.raises(ValueError):
            load_configuration(path)


def test_unknown_profile_rejected() -> None:
    with pytest.raises(ValueError, match="profile_unknown"):
        observed_profile(panel(), "OTHER", {}, {})


def test_new_fixture_creation_has_canonical_tables_and_bound_dataset_sets(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    """Régression des deux défauts d'intégration, PLINK intégralement simulé."""
    import json
    import yaml
    from effet_fondateur.calibration.export_blinded_inputs import SAMPLE_COLUMNS, pipeline_configuration
    from effet_fondateur.calibration.qc_bench import write_table
    from effet_fondateur.contracts import validate_tsv_table
    source = panel()
    preparation, root = tmp_path / "preparation", tmp_path / "run"
    metadata = preparation / "blinded_inputs/metadata"
    metadata.mkdir(parents=True)
    root.mkdir()
    registry = [dict(zip(SAMPLE_COLUMNS, (sample,sample+".txt","F1",sample,"0","0","UNKNOWN","UNKNOWN","UNSPECIFIED","","","","true","true",""),strict=True)) for sample in source.samples]
    write_table(metadata / "samples.tsv", SAMPLE_COLUMNS, registry)
    genotype_rows = [dict(SAMPLE_ID=sample, GENOTYPE="A/G", GENOTYPE_SOURCE="SIMULATION", MEASUREMENT_METHOD="NOT_MEASURED") for sample in source.samples]
    write_table(preparation / "blinded_inputs/target_genotypes.tsv", tuple(genotype_rows[0]), genotype_rows)
    target = dict(chromosome=5, position_bp=2, project_variant_id="target", ref="A", alt="G")
    config = pipeline_configuration(target)
    for name in ("target_variant_metadata", "reference_panel_catalog", "ancestry_reference_catalog", "genetic_map_catalog", "phasing_reference_selection"):
        path = tmp_path / (name + ".json")
        path.write_text('{}')
        config["inputs"][name] = str(path)
    (root / "config.resolved.yaml").write_text(yaml.safe_dump(config,sort_keys=False))
    def fake_command(command: list[str], run_root: Path, timeout: int) -> int:
        prefix = Path(command[command.index("--out")+1])
        prefix.with_suffix(".bed").write_bytes(b"\x6c\x1b\x01")
        prefix.with_suffix(".fam").write_text("S1 S1 0 0 0 -9\nS2 S2 0 0 0 -9\n")
        variants = [row for row in source.variants if "--chr" not in command or row[0] == "chr5"]
        prefix.with_suffix(".bim").write_text(''.join(f'{row[0].removeprefix("chr")} {row[2]} 0 {row[1]} {row[4]} {row[3]}\n' for row in variants))
        return 0
    monkeypatch.setattr("effet_fondateur.calibration.phasing_bench._command",fake_command)
    artifacts = fixture_inputs(root,preparation,source,dict(command_timeout_seconds=1),target)
    audit = next(item for item in artifacts if item["artifact_id"] == "target_genotype_audit")
    validate_tsv_table(Path(audit["path"]),"target_genotype_audit.schema.json")
    for prefix in ("genomewide_pre_qc","target_chromosome_all_qc"):
        descriptor = json.loads(Path(next(item for item in artifacts if item["artifact_id"] == prefix+"_dataset")["path"]).read_text())
        for item in artifacts:
            if item["artifact_id"].startswith(prefix+"_"):
                assert item["sample_set_id"] == descriptor["sample_set_id"]
                assert item["variant_set_id"] == descriptor["variant_set_id"]
