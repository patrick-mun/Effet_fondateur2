"""D22 : préparation aveugle de structure et fixtures amont explicites.

Le banc IID est distinct des panels 1000G et du FAIL D21. Il teste 08/09
comme composants ; aucune exécution de 04–07 n'est fabriquée.
"""

from __future__ import annotations

import argparse
from copy import deepcopy
from datetime import datetime, timezone
import hashlib
from pathlib import Path
import shutil
import tempfile
from typing import Any

import numpy as np
import yaml

from effet_fondateur.calibration.evaluate_technical import contained_file, raw_table, validate_schema
from effet_fondateur.calibration.export_blinded_inputs import _write_bundle
from effet_fondateur.calibration.materialize_ibd import _json, _sha
from effet_fondateur.calibration.qc_bench import file_inventory, write_json, write_table
from effet_fondateur.calibration.simulation.measurement_errors import MeasurementPanel
from effet_fondateur.calibration.technical_evaluation import TechnicalEvaluationError
from effet_fondateur.calibration.technical_scenarios import new_workspace_run, technical_configuration, validate_scenario_config
from effet_fondateur.contracts import load_pipeline_config, validate_tsv_table


METHOD_ID = "iid_structure_component_bench_v1"
PRODUCER = "D22_UPSTREAM_CONTRACT_FIXTURE"
CASES = ("BASELINE", "PROJECTED_STRESS", "LABEL_PERMUTATION", "MISSING_5PCT")
UPSTREAM_STAGES = ("initialize_run", "validate_sources", "build_sample_registry", "convert_acpa")
COHORTS = ("controls_unrelated", "target_carriers_all", "target_carriers_independent", "family_noncarriers",
           "affected_noncarriers", "target_chromosome_all_qc", "genomewide_structure_reference")


def validate_configuration(config: dict[str, Any]) -> None:
    """Valide le domaine versionné et refuse les graines confondues."""
    validate_schema(config, "structure_bench_config")
    if config["seed"] == config["missingness_seed"]:
        raise TechnicalEvaluationError("d22_seed_collision")


def checked_document(root: Path, name: str, state_key: str, state: str) -> dict[str, Any]:
    """Contrôle clôture, méthode, inventaire exhaustif et confinement physique."""
    document = _json(contained_file(root, name))
    validate_schema(document, "structure_bench_document")
    if document.get("schema_version") != "1.0.0" or document.get("method_id") != METHOD_ID or document.get(state_key) != state or document.get("calibrated") is not False or document.get("founder_effect_proven") is not False:
        raise TechnicalEvaluationError("d22_document_state_invalid")
    artifacts = document.get("artifacts", [])
    if not artifacts or len({row["path"] for row in artifacts}) != len(artifacts):
        raise TechnicalEvaluationError("d22_inventory_invalid")
    for row in artifacts:
        if _sha(contained_file(root, row["path"])) != row["sha256"]:
            raise TechnicalEvaluationError("d22_artifact_changed")
    if file_inventory(root, (name,)) != artifacts:
        raise TechnicalEvaluationError("d22_uninventoried_file")
    return document


def generate_cases(config: dict[str, Any]) -> tuple[list[MeasurementPanel], dict[str, Any]]:
    """Tire des copies IID ; retourne quatre observations et une vérité privée.

    Les strates sont mathématiques. Le mélange 50/50 tire l'origine de chaque
    copie à chaque locus, sans segments ni population géographique inférée.
    """
    validate_configuration(config)
    rng = np.random.default_rng(config["seed"])
    count, per_group = 22 * config["markers_per_autosome"], config["pure_per_group"]
    groups = ["A"] * per_group + ["B"] * per_group + ["MIXED"] * config["admixed_count"] + ["RARE"] * config["rare_count"]
    reference_count = len(groups)
    frequencies = config["frequencies"]
    calls = []
    for group in groups:
        if group == "MIXED":
            origins = rng.integers(0, 2, (count, 2))
            probability = np.array(frequencies)[origins]
        elif group == "RARE":
            probability = np.broadcast_to(np.where(np.arange(count) < count // 2, frequencies[1], frequencies[0])[:, None], (count, 2))
        else:
            probability = frequencies[int(group == "B")]
        calls.append(rng.binomial(1, probability, size=(count, 2)).astype(np.int8))
    donors = projected_donors(config)
    calls.extend(calls[index].copy() for index in donors)
    groups.extend(groups[index] for index in donors)
    baseline = np.sort(np.stack(calls, axis=1), axis=2)
    baseline[0, :reference_count] = 0
    baseline[0, reference_count:] = (0, 1)
    baseline[1, :reference_count] = -1
    variants = tuple((str(chromosome), str(1_000_000 + offset * 10_000), f"struct_chr{chromosome}_{offset:04d}", "A", "G", ".", "PASS", ".", "GT")
                     for chromosome in range(1, 23) for offset in range(config["markers_per_autosome"]))
    samples = tuple("S" + hashlib.sha256(f'{METHOD_ID}:{config["seed"]}:{index}'.encode()).hexdigest()[:16] for index in range(len(groups)))
    panels = []
    for case in CASES:
        observed = baseline.copy()
        if case == "PROJECTED_STRESS":
            projected = observed[:, reference_count:]
            projected[projected >= 0] = 1 - projected[projected >= 0]
        elif case == "MISSING_5PCT":
            missing = np.random.default_rng(config["missingness_seed"]).random(observed.shape[:2]) < config["missingness_fraction"]
            observed[missing] = -1
        panels.append(MeasurementPanel(samples, variants, np.sort(observed, axis=2)))
    return panels, dict(groups=groups, reference_count=reference_count, projected_donors=donors,
                        injected_outlier_indices=[index for index, group in enumerate(groups) if group == "RARE"])


def projected_donors(config: dict[str, Any]) -> list[int]:
    """Décrit les copies projetées, sans déduire leur identité des génotypes."""
    per_group = config["pure_per_group"]
    reference_count = 2 * per_group + config["admixed_count"] + config["rare_count"]
    return list(range(config["projected_per_group"])) + list(range(per_group, per_group + config["projected_per_group"])) + list(range(reference_count - config["projected_rare_count"], reference_count))


def registry_rows(panel: MeasurementPanel, case: str) -> list[dict[str, str]]:
    """Registre neutre : rôles techniques indépendants des strates et génotypes."""
    from effet_fondateur.calibration.export_blinded_inputs import SAMPLE_COLUMNS
    rows = []
    for index, sample in enumerate(panel.samples):
        alternate = (index + int(case == "LABEL_PERMUTATION")) % 2
        family = "F" + hashlib.sha256(f'{METHOD_ID}:family:{index // 6}'.encode()).hexdigest()[:16]
        values = (sample, sample + ".txt", family, sample, "0", "0", "UNKNOWN", "AFFECTED" if (index + int(case == "LABEL_PERMUTATION")) % 5 == 0 else "UNKNOWN",
                  "CONTROL" if alternate == 0 else "OTHER", "", "", "LOT_A" if alternate == 0 else "LOT_B", "true", "true", "")
        rows.append(dict(zip(SAMPLE_COLUMNS, values, strict=True)))
    return rows


def upstream_configuration(config: dict[str, Any], target_config: dict[str, Any]) -> dict[str, Any]:
    """Profil 00–03 seulement ; les autres entrées portent leur vraie provenance."""
    result = technical_configuration(target_config)
    result["project"]["project_id"] = "d22_structure_component_bench"
    result["stages"].pop("prepare_target_variant_dataset")
    result["inputs"].update(target_variant_metadata=None, target_variant_genotypes=None)
    result["stages"]["validate_sources"]["parameters"]["required_autosomes"] = list(range(1, 23))
    result["stages"]["build_sample_registry"]["parameters"]["manual_approval"] = deepcopy(config["manual_approval"])
    result["tools"]["plink"] = config["plink_executable"]
    return result


def fixture_tables(panel: MeasurementPanel, registry: list[dict[str, str]], config: dict[str, Any], target: dict[str, Any]) -> dict[str, list[dict[str, Any]]]:
    """Contrats 04/05/07 explicitement simulés, sans fonction amont prétendue."""
    reference_count = 2 * config["pure_per_group"] + config["admixed_count"] + config["rare_count"]
    qc, proposal, genotypes = [], [], []
    cycle = (f'{target["ref"]}/{target["ref"]}', f'{target["ref"]}/{target["alt"]}', f'{target["alt"]}/{target["alt"]}')
    for index, row in enumerate(registry):
        base = {key: row[key] for key in ("SAMPLE_ID", "FID", "IID")}
        missing = int((panel.alleles[:, index, 0] < 0).sum())
        rate = f'{missing / len(panel.variants):.12f}'
        qc.append(dict(**base, ARRAY_BATCH=row["ARRAY_BATCH"], N_MISS=str(missing), N_GENO=str(len(panel.variants)), F_MISS=rate,
                       OBSERVED_HETEROZYGOSITY="", EXPECTED_HETEROZYGOSITY="", INBREEDING_F="", HETEROZYGOSITY_Z=str(index % 7), QC_STATUS="PASS", EXCLUSION_CODE=""))
        included = index < reference_count
        representative = index if included else projected_donors(config)[index - reference_count]
        proposal.append(dict(**base, COMPONENT_ID=f"kinship_component_{representative+1:06d}", PROPOSED_INCLUDED=str(included).lower(),
                             PROPOSAL_STATUS="RETAIN" if included else "EXCLUDE_RELATED", REPRESENTATIVE_SAMPLE_ID="" if included else panel.samples[representative],
                             RELATION_BASIS="" if included else "KING_RELATEDNESS", F_MISS=rate, ABS_HETEROZYGOSITY_Z=str(index % 7),
                             DECISION_SOURCE="stage_07_quality_aware_greedy_maximal", MANUAL_VALIDATION_REQUIRED=str(not included).lower()))
        genotypes.append(dict(**base, GENOTYPE=cycle[index % 3], GENOTYPE_SOURCE="EXPLICIT_SYNTHETIC_D22", MEASUREMENT_METHOD="DETERMINISTIC_FIXTURE_NOT_A_MEASUREMENT", STATUS="ACCEPTED"))
    return dict(qc_individual_metrics=qc, independent_set_proposal=proposal, target_genotype_audit=genotypes)


def prepare_structure_bench(*, configuration_path: Path, workspace_root: Path, run_id: str) -> dict[str, Any]:
    """Publie les quatre cas sans lancer PLINK ni ouvrir d'autres panels."""
    destination = new_workspace_run(workspace_root, run_id)
    config = _json(configuration_path)
    panels, truth = generate_cases(config)
    reference = contained_file(configuration_path.absolute().parent, config["target_configuration"])
    if _sha(reference) != config["target_configuration_sha256"]:
        raise TechnicalEvaluationError("d22_target_configuration_changed")
    target_config = _json(reference)
    validate_scenario_config(target_config)
    temporary = Path(tempfile.mkdtemp(prefix=f".{run_id}.", dir=destination.parent))
    try:
        private = temporary / "private_truth"
        private.mkdir(mode=0o700)
        shutil.copyfile(configuration_path, private / "configuration.json")
        write_json(private / "structure_truth.json", truth)
        execution = dict(producer=PRODUCER, scope="ISOLATED_COMPONENTS_08_09_NOT_04_07", population_parameters=config["population_parameters"],
                         fixture_review_policy=config["fixture_review_policy"], freeze_parameters=dict(control_group_labels=["CONTROL"], required_nonempty_cohorts=list(COHORTS)))
        bundles = []
        for index, (case, panel) in enumerate(zip(CASES, panels, strict=True)):
            bundle_id = f"B{index:02d}"
            bundle = temporary / "blinded_inputs/bundles" / bundle_id
            _write_bundle(bundle, panel, dict(target=target_config["target_metadata"]))
            registry = registry_rows(panel, case)
            write_table(bundle / "metadata/samples.tsv", tuple(registry[0]), registry)
            (bundle / "pipeline.yaml").write_text(yaml.safe_dump(upstream_configuration(config, target_config), sort_keys=False), encoding="utf-8")
            load_pipeline_config(bundle / "pipeline.yaml")
            fixture = bundle / "contract_fixtures"
            fixture.mkdir()
            for name, rows in fixture_tables(panel, registry, config, target_config["target_metadata"]).items():
                path = fixture / (name + ".tsv")
                write_table(path, tuple(rows[0]), rows)
                validate_tsv_table(path, name + ".schema.json")
            (fixture / "target_metadata.yaml").write_text(yaml.safe_dump(target_config["target_metadata"], sort_keys=False), encoding="utf-8")
            write_json(private / f"{bundle_id}.json", dict(case=case))
            bundles.append(dict(bundle_id=bundle_id, sample_count=len(panel.samples), marker_count=len(panel.variants)))
        write_json(temporary / "blinded_inputs/execution.json", execution)
        document = dict(schema_version="1.0.0", method_id=METHOD_ID, manifest_state="FINALIZED", created_at=datetime.now(timezone.utc).isoformat(),
                        configuration_sha256=_sha(configuration_path), bundles=bundles, artifacts=file_inventory(temporary), calibrated=False, founder_effect_proven=False)
        write_json(temporary / "panel_manifest.json", document)
        checked_document(temporary, "panel_manifest.json", "manifest_state", "FINALIZED")
        for path in temporary.rglob("*"):
            path.chmod(0o700 if path.is_dir() else 0o600)
        temporary.rename(destination)
    except Exception:
        shutil.rmtree(temporary)
        raise
    return document


def main() -> int:
    """CLI de préparation uniquement ; sortie utilisateur sans identifiants."""
    parser = argparse.ArgumentParser(description="Préparation D22 de structure et gel")
    for name in ("configuration-path", "workspace-root"):
        parser.add_argument("--" + name, type=Path, required=True)
    parser.add_argument("--run-id", required=True)
    document = prepare_structure_bench(**vars(parser.parse_args()))
    print(dict(state=document["manifest_state"], cases=len(document["bundles"])))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
