"""Préparation D20 : stress QC appariés IID, distincts du panel biologique.

Ce module ne lance aucun outil. Les sorties contiennent uniquement des
sujets synthétiques sans correspondance réelle, confinés au workspace ignoré.
La géométrie artificielle est autosomique ; elle ne simule pas la puce 1000G.
"""

from __future__ import annotations

import argparse
from copy import deepcopy
import csv
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
import shutil
import tempfile
from typing import Any

import numpy as np
import yaml

from effet_fondateur.calibration.evaluate_technical import contained_file, raw_table, validate_schema
from effet_fondateur.calibration.export_blinded_inputs import _write_bundle
from effet_fondateur.calibration.materialize_ibd import _json, _sha
from effet_fondateur.calibration.qc_oracle import COHORTS, preliminary_expectation
from effet_fondateur.calibration.simulation.measurement_errors import MeasurementPanel
from effet_fondateur.calibration.technical_evaluation import BASES, GenotypeMatrix, TechnicalEvaluationError, observation_matrix
from effet_fondateur.calibration.technical_scenarios import new_workspace_run, technical_configuration, validate_scenario_config
from effet_fondateur.contracts import load_pipeline_config, validate_cohorts_frozen


METHOD_ID = "iid_autosomal_qc_component_bench_v1"
CASES = ("CLEAN", "SAMPLE_MISSING", "VARIANT_MISSING", "HET_OUTLIER", "BATCH_MISSING",
         "DUPLICATE", "RARE_VARIANT", "CONTROL_HWE", "TARGET_EXCEPTION", "LOCAL_SAMPLE_MISSING", "CONTROL_VARIANT_MISSING",
         "CARRIER_VARIANT_MISSING", "LOCAL_VARIANT_MISSING", "POST05_SAMPLE_MISSING")
COHORT_COLUMNS = ("COHORT_ID", "SAMPLE_ID", "ROLE", "INCLUDED", "EXCLUSION_CODE", "DECISION_SOURCE", "INDEPENDENT_UNIT_ID", "FAMILY_REPRESENTATIVE")
PRODUCER = "D20_COHORT_CONTRACT_FIXTURE"


def write_json(path: Path, value: Any) -> None:
    """Sérialise sans NaN ; le répertoire parent est créé par l'appelant."""
    path.write_text(json.dumps(value, indent=2, allow_nan=False) + "\n", encoding="utf-8")


def write_table(path: Path, columns: tuple[str, ...], rows: list[dict[str, Any]]) -> None:
    """Écrit une table explicite sans colonnes supplémentaires silencieuses."""
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def file_inventory(root: Path, excluded: tuple[str, ...] = ()) -> list[dict[str, str]]:
    """Inventorie tous les fichiers d'un run, hors document auto-référent."""
    return [dict(path=path.relative_to(root).as_posix(), sha256=_sha(path))
            for path in sorted(root.rglob("*")) if path.is_file() and path.relative_to(root).as_posix() not in excluded]


def checked_inventory(root: Path, document_name: str, state_key: str, state: str) -> dict[str, Any]:
    """Contrôle clôture, confinement, exhaustivité et empreintes avant lecture."""
    document_path = contained_file(root, document_name)
    if not document_path.is_file():
        raise TechnicalEvaluationError("d20_closure_document_missing")
    document = _json(document_path)
    validate_schema(document, "qc_bench_document")
    if document.get(state_key) != state or document.get("method_id") != METHOD_ID:
        raise TechnicalEvaluationError("d20_document_state_invalid")
    artifacts = document.get("artifacts", [])
    if not artifacts or len({item["path"] for item in artifacts}) != len(artifacts):
        raise TechnicalEvaluationError("d20_inventory_invalid")
    for item in artifacts:
        path = contained_file(root, item["path"])
        if not path.is_file() or _sha(path) != item["sha256"]:
            raise TechnicalEvaluationError("d20_artifact_checksum_mismatch")
    if file_inventory(root, (document_name,)) != artifacts:
        raise TechnicalEvaluationError("d20_uninventoried_file")
    return document


def generate_cases(configuration: dict[str, Any]) -> list[tuple[str, MeasurementPanel, list[tuple[str, str]]]]:
    """Fond reproductible de 1 408 loci ; quatorze cas sans retoucher seuils.

    Le troisième élément est la vérité de paire clonée, pas une parenté IBD.
    Les rotations de manques de lot sont déterministes après le tirage initial.
    """
    validate_schema(configuration, "qc_bench_config")
    sample_count = configuration["sample_count"]
    per_chromosome = configuration["markers_per_autosome"]
    variants = tuple((str(chromosome), str(1_000_000 + marker * 10_000), f"qc_chr{chromosome}_{marker:03d}", "A", "G", ".", "PASS", ".", "GT")
                     for chromosome in range(1, 23) for marker in range(per_chromosome))
    baseline = np.sort(np.random.default_rng(configuration["seed"]).binomial(
        1, configuration["allele_frequency"], (len(variants), sample_count, 2)).astype(np.int8), axis=2)
    output = []
    for case_index, case in enumerate(CASES):
        samples = tuple("S" + hashlib.sha256(f'{METHOD_ID}:{configuration["seed"]}:{case_index}:{index}'.encode()).hexdigest()[:16] for index in range(sample_count))
        calls = baseline.copy()
        pairs = []
        if case == "SAMPLE_MISSING":
            calls[:int(len(variants) * .3), 0] = -1
        elif case == "VARIANT_MISSING":
            calls[0, :sample_count // 4] = -1
        elif case == "HET_OUTLIER":
            calls[:, 0] = (0, 1)
        elif case == "BATCH_MISSING":
            batch_b = np.arange(1, sample_count, 2)
            for marker in range(len(variants)):
                if marker % 5:
                    calls[marker, batch_b[(marker + np.arange(5)) % len(batch_b)]] = -1
        elif case == "DUPLICATE":
            calls[:, -1] = calls[:, -2]
            pairs = [(samples[-2], samples[-1])]
        elif case == "RARE_VARIANT":
            calls[[12, 4 * per_chromosome + 12]] = 0
            calls[[12, 4 * per_chromosome + 12], 0] = (0, 1)
        elif case == "CONTROL_HWE":
            calls[13, :sample_count // 2] = (0, 1)
        elif case == "LOCAL_SAMPLE_MISSING":
            first = 4 * per_chromosome
            calls[first:first + int(per_chromosome * .3), 0] = -1
        elif case == "CONTROL_VARIANT_MISSING":
            calls[14, :8] = -1
        elif case == "CARRIER_VARIANT_MISSING":
            calls[15, sample_count // 2:sample_count // 2 + 8] = -1
        elif case == "LOCAL_VARIANT_MISSING":
            calls[4 * per_chromosome, :sample_count // 4] = -1
        elif case == "POST05_SAMPLE_MISSING":
            calls[:200, :10] = -1
            calls[-130:, [20, 100]] = -1
        output.append((case, MeasurementPanel(samples, variants, calls), pairs))
    return output


def with_target(matrix: GenotypeMatrix, target: dict[str, Any], genotypes: list[str]) -> GenotypeMatrix:
    """Ajoute la cible aux seuls loci de son chromosome ; tri déterministe bp/ID."""
    selected = [index for index, variant in enumerate(matrix.variants) if variant[0] == target["chromosome"]]
    variants = [matrix.variants[index] for index in selected] + [(target["chromosome"], target["position_bp"], target["project_variant_id"])]
    calls = np.array([[sorted(BASES.index(base) for base in genotype.split("/")) for genotype in genotypes]], dtype=np.int8)
    alleles = np.concatenate((matrix.alleles[selected], calls))
    order = sorted(range(len(variants)), key=lambda index: variants[index])
    return GenotypeMatrix(matrix.samples, tuple(variants[index] for index in order), alleles[order])


def _cohort_fixture(path: Path, matrix: GenotypeMatrix, preliminary: GenotypeMatrix, genotypes: list[str], target: dict[str, Any], case: str) -> dict[str, list[str]]:
    """Éligibilité connue avant run, génotypes explicites, clone non indépendant."""
    eligible = {sample for _, sample in preliminary.samples}
    memberships = {COHORTS[0]: [sample for (_, sample), genotype in zip(matrix.samples, genotypes, strict=True)
                               if genotype == f'{target["ref"]}/{target["ref"]}' and sample in eligible],
                   COHORTS[1]: [sample for (_, sample), genotype in zip(matrix.samples, genotypes, strict=True)
                               if genotype != f'{target["ref"]}/{target["ref"]}' and sample in eligible],
                   COHORTS[2]: [sample for _, sample in matrix.samples]}
    if case == "DUPLICATE":
        memberships[COHORTS[1]].remove(matrix.samples[-1][1])
    rows = []
    for cohort, samples in memberships.items():
        if not samples:
            raise TechnicalEvaluationError("d20_empty_prespecified_cohort")
        for sample in samples:
            rows.append(dict(zip(COHORT_COLUMNS, (cohort, sample, "CONTROL" if cohort == COHORTS[0] else "CARRIER" if cohort == COHORTS[1] else "OTHER",
                                                  "true", "", PRODUCER, sample, "false"), strict=True)))
        (path / f"{cohort}.keep").write_text("".join(f"{sample}\t{sample}\n" for sample in samples), encoding="utf-8")
    write_table(path / "cohorts.tsv", COHORT_COLUMNS, rows)
    return memberships


def prepare_qc_bench(*, configuration_path: Path, workspace_root: Path, run_id: str) -> dict[str, Any]:
    """Publie entrées aveugles, fixtures publiques et vérité privée immuables."""
    destination = new_workspace_run(workspace_root, run_id)
    config = _json(configuration_path)
    validate_schema(config, "qc_bench_config")
    reference = contained_file(configuration_path.absolute().parent, config["target_configuration"])
    if _sha(reference) != config["target_configuration_sha256"]:
        raise TechnicalEvaluationError("d20_target_reference_checksum_mismatch")
    target_config = _json(reference)
    validate_scenario_config(target_config)
    temporary = Path(tempfile.mkdtemp(prefix=f".{run_id}.", dir=destination.parent))
    try:
        private = temporary / "private_truth"
        private.mkdir(mode=0o700)
        (temporary / "blinded_inputs").mkdir(mode=0o700)
        shutil.copyfile(configuration_path, private / "configuration.json")
        shutil.copyfile(reference, private / "target_configuration.json")
        write_json(temporary / "blinded_inputs/qc_execution.json", dict(final_parameters=config["final_parameters"],
                                                                      cohort_producer=PRODUCER, scope="ISOLATED_FINAL_QC_NOT_STAGE_09"))
        bundles = []
        for index, (case, panel, pairs) in enumerate(generate_cases(config)):
            bundle_id = f"B{index:02d}"
            bundle = temporary / "blinded_inputs/bundles" / bundle_id
            _write_bundle(bundle, panel, dict(target=target_config["target_metadata"]))
            metadata = raw_table(bundle / "metadata/samples.tsv")
            batches = ["LOT_A" if rank % 2 == 0 else "LOT_B" for rank in range(len(panel.samples))]
            for row, batch in zip(metadata, batches, strict=True):
                row.update(ARRAY_BATCH=batch, INCLUDE_GENOMEWIDE="true")
            write_table(bundle / "metadata/samples.tsv", tuple(metadata[0]), metadata)
            matrix = observation_matrix(panel, metadata)
            expected = preliminary_expectation(matrix, batches, config["preliminary_parameters"])
            target = target_config["target_metadata"]
            genotypes = [f'{target["ref"]}/{target["ref"]}' if rank < len(panel.samples) // 2 else f'{target["ref"]}/{target["alt"]}' for rank in range(len(panel.samples))]
            if case == "TARGET_EXCEPTION":
                genotypes = [f'{target["ref"]}/{target["ref"]}'] * (len(panel.samples) - 1) + [f'{target["ref"]}/{target["alt"]}']
            write_table(bundle / "target_genotypes.tsv", ("SAMPLE_ID", "GENOTYPE", "GENOTYPE_SOURCE", "MEASUREMENT_METHOD"),
                        [dict(SAMPLE_ID=sample, GENOTYPE=genotype, GENOTYPE_SOURCE="EXPLICIT_SYNTHETIC_D20", MEASUREMENT_METHOD="DETERMINISTIC_FIXTURE_NOT_A_MEASUREMENT")
                         for sample, genotype in zip(panel.samples, genotypes, strict=True)])
            (bundle / "target_metadata.yaml").write_text(yaml.safe_dump(target, sort_keys=False), encoding="utf-8")
            pipeline = technical_configuration(target_config)
            pipeline["project"]["project_id"] = "d20_qc_component_bench"
            pipeline["tools"]["plink"] = config["plink_executable"]
            pipeline["stages"]["validate_sources"]["parameters"]["required_autosomes"] = list(range(1, 23))
            pipeline["stages"]["build_sample_registry"]["parameters"]["manual_approval"] = deepcopy(config["manual_approval"])
            pipeline["stages"]["qc_preliminary"] = dict(enabled=True, parameters=config["preliminary_parameters"])
            (bundle / "pipeline.yaml").write_text(yaml.safe_dump(pipeline, sort_keys=False), encoding="utf-8")
            load_pipeline_config(bundle / "pipeline.yaml")
            fixture = temporary / "blinded_inputs/cohort_fixtures" / bundle_id
            fixture.mkdir(parents=True)
            memberships = _cohort_fixture(fixture, matrix, expected.retained, genotypes, target, case)
            validate_cohorts_frozen(fixture / "cohorts.tsv", set(panel.samples))
            write_json(private / f"{bundle_id}.json", dict(case=case, duplicate_pairs=pairs, cohorts=memberships,
                                                           target_genotypes=genotypes, sample_count=len(panel.samples), variant_count=len(panel.variants)))
            bundles.append(dict(bundle_id=bundle_id, sample_count=len(panel.samples), marker_count=len(panel.variants)))
        manifest = dict(schema_version="1.0.0", method_id=METHOD_ID, manifest_state="FINALIZED", created_at=datetime.now(timezone.utc).isoformat(),
                        configuration_sha256=_sha(configuration_path), bundles=bundles, artifacts=file_inventory(temporary),
                        calibrated=False, founder_effect_proven=False)
        write_json(temporary / "panel_manifest.json", manifest)
        checked_inventory(temporary, "panel_manifest.json", "manifest_state", "FINALIZED")
        for path in temporary.rglob("*"):
            path.chmod(0o700 if path.is_dir() else 0o600)
        temporary.rename(destination)
    except Exception:
        shutil.rmtree(temporary)
        raise
    return manifest


def main() -> int:
    """CLI de préparation seulement ; stdout réduit aux effectifs et au statut."""
    parser = argparse.ArgumentParser(description="Préparation du banc D20, sans outil externe")
    parser.add_argument("--configuration-path", type=Path, required=True)
    parser.add_argument("--workspace-root", type=Path, required=True)
    parser.add_argument("--run-id", required=True)
    document = prepare_qc_bench(**vars(parser.parse_args()))
    print(json.dumps(dict(bundles=len(document["bundles"]), state=document["manifest_state"])))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
