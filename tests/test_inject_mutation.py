import json
from pathlib import Path

import pytest

from simulation_genotype_famille.inject_mutation import (
    MutationInjectionError,
    inject_mutation,
    load_mutation_metadata,
)


def write_source_dataset(directory: Path) -> tuple[Path, Path, Path]:
    source_prefix = directory / "source" / "genotype_data"
    source_prefix.parent.mkdir()
    source_prefix.with_suffix(".map").write_text(
        "19\trs1\t0\t100\n19\trs2\t0\t300\n",
        encoding="utf-8",
    )
    source_prefix.with_suffix(".ped").write_text(
        "F1\tCASE\t0\t0\t1\t2\tA\tA\tC\tT\n"
        "F2\tCONTROL\t0\t0\t2\t1\tA\tG\tC\tC\n",
        encoding="utf-8",
    )
    groups_path = source_prefix.parent / "groupes.txt"
    groups_path.write_text("CASE\tATTEINT\nCONTROL\tSAINS\n", encoding="utf-8")
    (source_prefix.parent / "cas.txt").write_text("F1\tCASE\n", encoding="utf-8")
    (source_prefix.parent / "temoins.txt").write_text(
        "F2\tCONTROL\n", encoding="utf-8"
    )
    config_path = directory / "mutation.json"
    config_path.write_text(
        json.dumps(
            {
                "gene": "DOCK6",
                "mutation_name": "Variant test",
                "chromosome": 19,
                "position_bp": 200,
                "reference_allele": "G",
                "alternate_allele": "T",
                "hgvs_c": "c.1G>T",
                "hgvs_p": "p.Gly1Val",
                "transcript": "NM_TEST.1",
                "disease": "Maladie test",
            }
        ),
        encoding="utf-8",
    )
    return source_prefix, groups_path, config_path


def test_injection_adds_sorted_marker_genotypes_metadata_and_audit(tmp_path):
    source_prefix, groups_path, config_path = write_source_dataset(tmp_path)
    output_prefix = tmp_path / "output" / "genotype_data"
    source_ped_before = source_prefix.with_suffix(".ped").read_text()

    report = inject_mutation(
        source_prefix=source_prefix,
        output_prefix=output_prefix,
        mutation_config_path=config_path,
        groups_path=groups_path,
        group_genotype_values=["ATTEINT=ALT/ALT", "SAINS=REF/REF"],
    )

    assert output_prefix.with_suffix(".map").read_text().splitlines() == [
        "19\trs1\t0\t100",
        "19\trsMUT_200\t0\t200",
        "19\trs2\t0\t300",
    ]
    output_rows = [
        line.split()
        for line in output_prefix.with_suffix(".ped").read_text().splitlines()
    ]
    assert output_rows[0][6:] == ["A", "A", "T", "T", "C", "T"]
    assert output_rows[1][6:] == ["A", "G", "G", "G", "C", "C"]
    assert source_prefix.with_suffix(".ped").read_text() == source_ped_before
    assert (output_prefix.parent / "groupes.txt").read_text() == groups_path.read_text()
    assert (output_prefix.parent / "cas.txt").is_file()
    assert (output_prefix.parent / "temoins.txt").is_file()
    assert report["genotype_counts"] == {"ALT/ALT": 1, "REF/REF": 1}

    mutation_info = json.loads((output_prefix.parent / "mutation_info.json").read_text())
    project_info = json.loads((output_prefix.parent / "project_info.json").read_text())
    audit_lines = (output_prefix.parent / "mutation_genotype_audit.tsv").read_text().splitlines()
    assert mutation_info["hgvs_c"] == "c.1G>T"
    assert project_info["HGVS p."] == "p.Gly1Val"
    assert audit_lines[1].endswith("ALT/ALT\tT\tT\tgroup:ATTEINT")
    assert report["source_sha256"]["ped"] != report["output_sha256"]["ped"]


def test_individual_assignment_overrides_group_rule(tmp_path):
    source_prefix, groups_path, config_path = write_source_dataset(tmp_path)
    genotypes_path = tmp_path / "genotypes.tsv"
    genotypes_path.write_text(
        "FID\tIID\tGENOTYPE\nF1\tCASE\tREF/ALT\n", encoding="utf-8"
    )

    report = inject_mutation(
        source_prefix=source_prefix,
        output_prefix=tmp_path / "output" / "genotype_data",
        mutation_config_path=config_path,
        groups_path=groups_path,
        group_genotype_values=["ATTEINT=ALT/ALT", "SAINS=REF/REF"],
        individual_genotypes_path=genotypes_path,
    )

    assert report["genotype_counts"] == {"REF/ALT": 1, "REF/REF": 1}
    audit = (tmp_path / "output" / "mutation_genotype_audit.tsv").read_text()
    assert "F1\tCASE\tATTEINT\tREF/ALT\tG\tT\tindividual_override" in audit


def test_injection_requires_an_explicit_genotype_for_every_individual(tmp_path):
    source_prefix, groups_path, config_path = write_source_dataset(tmp_path)

    with pytest.raises(MutationInjectionError, match="F2/CONTROL"):
        inject_mutation(
            source_prefix=source_prefix,
            output_prefix=tmp_path / "output" / "genotype_data",
            mutation_config_path=config_path,
            groups_path=groups_path,
            group_genotype_values=["ATTEINT=ALT/ALT"],
        )


def test_injection_refuses_overwrite_and_source_directory(tmp_path):
    source_prefix, groups_path, config_path = write_source_dataset(tmp_path)
    output_prefix = tmp_path / "output" / "genotype_data"
    arguments = {
        "source_prefix": source_prefix,
        "mutation_config_path": config_path,
        "groups_path": groups_path,
        "group_genotype_values": ["ATTEINT=ALT/ALT", "SAINS=REF/REF"],
    }
    inject_mutation(output_prefix=output_prefix, **arguments)

    with pytest.raises(MutationInjectionError, match="Sorties déjà présentes"):
        inject_mutation(output_prefix=output_prefix, **arguments)
    first_ped = output_prefix.with_suffix(".ped").read_bytes()
    first_map = output_prefix.with_suffix(".map").read_bytes()
    inject_mutation(output_prefix=output_prefix, force=True, **arguments)
    assert output_prefix.with_suffix(".ped").read_bytes() == first_ped
    assert output_prefix.with_suffix(".map").read_bytes() == first_map
    with pytest.raises(MutationInjectionError, match="distinct du dossier source"):
        inject_mutation(output_prefix=source_prefix.parent / "injected", **arguments)


def test_mutation_metadata_rejects_non_plink_marker_and_invalid_allele(tmp_path):
    config_path = tmp_path / "mutation.json"
    base_config = {
        "gene": "DOCK6",
        "mutation_name": "Variant test",
        "chromosome": 19,
        "position_bp": 200,
        "reference_allele": "G",
        "alternate_allele": "T",
    }
    config_path.write_text(
        json.dumps({**base_config, "marker_id": "variant_200"}), encoding="utf-8"
    )
    with pytest.raises(MutationInjectionError, match="commencer par 'rsMUT'"):
        load_mutation_metadata(config_path)

    config_path.write_text(
        json.dumps({**base_config, "alternate_allele": "DEL"}), encoding="utf-8"
    )
    with pytest.raises(MutationInjectionError, match="bases génomiques simples"):
        load_mutation_metadata(config_path)
