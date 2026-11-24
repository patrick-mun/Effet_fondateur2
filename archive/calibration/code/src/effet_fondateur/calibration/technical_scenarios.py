"""D19 : dérivés techniques 00–04, sans modification des observations D18.

La cible technique est indépendante de la cible biologique D14–D18. Son
annotation est publique ; ses génotypes sont explicitement synthétiques,
répartis par rang canonique, sans tirage ni inférence clinique. Le négatif
ajoute cette cible à une copie des sondes pour éprouver le refus de doublon.
"""

from __future__ import annotations

import argparse
from copy import deepcopy
import csv
from datetime import datetime, timezone
import json
from pathlib import Path
import re
import shutil
import tempfile
from typing import Any

import numpy as np
import yaml

from effet_fondateur.calibration.evaluate_technical import contained_file, raw_table, validate_schema
from effet_fondateur.calibration.export_blinded_inputs import load_checked_export, pipeline_configuration
from effet_fondateur.calibration.materialize_ibd import _json, _sha, _text
from effet_fondateur.calibration.materialize_measurement_errors import deterministic_gzip
from effet_fondateur.calibration.simulation.measurement_errors import MeasurementPanel, read_panel, write_panel
from effet_fondateur.calibration.technical_evaluation import TechnicalEvaluationError
from effet_fondateur.contracts import load_pipeline_config, validate_json_document
from effet_fondateur.contracts.samples import validate_samples_master


METHOD_ID = "d19_technical_target_fixture_v1"
MANIFEST_PATH = "private_truth/d19_inputs/technical_panel_manifest.json"
GENOTYPE_SOURCE = "SYNTHETIC_EXPLICIT_TARGET_D19"
MEASUREMENT_METHOD = "DETERMINISTIC_FIXTURE_NOT_A_MEASUREMENT"
SCENARIOS = ("TARGET_ABSENT_INSERTION", "TARGET_PRESENT_REJECTION")


def new_workspace_run(workspace_root: Path, run_id: str) -> Path:
    """Confine une nouvelle publication locale avant toute lecture/écriture."""
    root = workspace_root.absolute()
    parts = tuple(part.casefold() for part in root.parts)
    if root.name != "calibration_workspace" or ".." in parts or "dock6_real" in parts or any(
        left == "data" and right == "input" for left, right in zip(parts, parts[1:])
    ) or any(path.is_symlink() for path in (root, *root.parents)) or not re.fullmatch(r"[A-Za-z0-9][A-Za-z0-9._-]*", run_id):
        raise TechnicalEvaluationError("d19_workspace_or_run_invalid")
    target = root / run_id
    if target.exists() or target.is_symlink():
        raise TechnicalEvaluationError("d19_run_exists")
    return target


def validate_scenario_config(config: dict[str, Any]) -> None:
    """Valide annotation publique et règle synthétique, sans revendication clinique."""
    validate_schema(config, "technical_scenario_config")
    metadata = config["target_metadata"]
    validate_json_document(metadata, "target_variant_metadata.schema.json")
    if metadata["ref"] == metadata["alt"] or any(len(metadata[key]) != 1 for key in ("ref", "alt")) or any(
        metadata[key]["source_id"] != config["source"]["source_id"] for key in ("reference_validation", "annotation_validation")
    ) or not config["source"]["url"].startswith("https://www.ncbi.nlm.nih.gov/clinvar/"):
        raise TechnicalEvaluationError("d19_public_target_annotation_invalid")


def technical_configuration(config: dict[str, Any]) -> dict[str, Any]:
    """Profil borné 00–04 ; l'approbation fournie n'est jamais créée par le runner."""
    validate_scenario_config(config)
    metadata = config["target_metadata"]
    result = pipeline_configuration(metadata)
    result["project"]["project_id"] = "technical_03_04_only"
    result["target"] = {key: metadata[key] for key in ("gene", "chromosome", "position_bp", "ref", "alt", "transcript", "project_variant_id", "rsid")}
    result["inputs"].update(target_variant_metadata="target_metadata.yaml", target_variant_genotypes="target_genotypes.tsv")
    result["tools"]["plink"] = config["plink_executable"]
    result["stages"]["build_sample_registry"]["parameters"] = dict(manual_approval=deepcopy(config["manual_approval"]))
    result["stages"].update(convert_acpa=dict(enabled=True, parameters=dict(marker_mode="intersection", plink_timeout_seconds=300)),
                            prepare_target_variant_dataset=dict(enabled=True, parameters=dict(plink_timeout_seconds=300)))
    return result


def target_expectation(metadata: list[dict[str, str]], config: dict[str, Any], bundle_id: str) -> dict[str, Any]:
    """Attente privée figée avant exécution : RR, RA, AA selon rang modulo 3.

    Le rang est celui du registre aveugle, sans lien avec sexe, famille,
    ascendance, clinique ou observation d'un vrai porteur. Ce stress ne
    représente pas une fréquence populationnelle et n'ajoute aucune graine.
    """
    target = config["target_metadata"]
    ref, alt = target["ref"], target["alt"]
    cycle = (f"{ref}/{ref}", f"{ref}/{alt}", f"{alt}/{alt}")
    result = dict(schema_version="1.0.0", assembly="GRCh38", bundle_id=bundle_id,
                  genotype_source="EXPLICIT_SYNTHETIC_TARGET_EXPECTATION",
                  target=dict(chromosome=target["chromosome"], position_bp=target["position_bp"], variant_id=target["project_variant_id"], ref=ref, alt=alt),
                  genotypes=[dict(sample_id=row["SAMPLE_ID"], genotype=cycle[index % 3]) for index, row in enumerate(metadata)])
    validate_schema(result, "technical_target_expectation")
    return result


def _write_target_files(bundle: Path, private: Path, config: dict[str, Any], expectation: dict[str, Any]) -> None:
    """Sérialise la vérité et l'entrée avec des provenances non expérimentales."""
    (private / f'{expectation["bundle_id"]}.target.json').write_text(json.dumps(expectation, indent=2) + "\n", encoding="utf-8")
    (bundle / "target_metadata.yaml").write_text(yaml.safe_dump(config["target_metadata"], sort_keys=False), encoding="utf-8")
    with (bundle / "target_genotypes.tsv").open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(("SAMPLE_ID", "GENOTYPE", "GENOTYPE_SOURCE", "MEASUREMENT_METHOD"))
        for row in expectation["genotypes"]:
            writer.writerow((row["sample_id"], row["genotype"], GENOTYPE_SOURCE, MEASUREMENT_METHOD))
    (bundle / "pipeline.yaml").write_text(yaml.safe_dump(technical_configuration(config), sort_keys=False), encoding="utf-8")


def _add_observed_target(bundle: Path, panel: MeasurementPanel, expectation: dict[str, Any]) -> None:
    """Négatif explicite : copie enrichie, jamais suppression d'une sonde D18."""
    target = expectation["target"]
    by_sample = {row["sample_id"]: row["genotype"] for row in expectation["genotypes"]}
    calls = np.array([[tuple((target["ref"], target["alt"]).index(base) for base in by_sample[sample].split("/")) for sample in panel.samples]], dtype=np.int8)
    variant = (str(target["chromosome"]), str(target["position_bp"]), target["variant_id"], target["ref"], target["alt"], ".", "PASS", ".", "GT")
    variants = panel.variants + (variant,)
    alleles = np.concatenate((panel.alleles, calls), axis=0)
    order = sorted(range(len(variants)), key=lambda index: (int(variants[index][0].removeprefix("chr")), int(variants[index][1]), variants[index][2]))
    observed = MeasurementPanel(panel.samples, tuple(variants[index] for index in order), alleles[order])
    with deterministic_gzip(bundle / "observed.vcf.gz") as handle:
        write_panel(handle, observed, observed.alleles, np.ones(len(variants), dtype=bool), phased=False)
    for sample in panel.samples:
        with (bundle / "acpa_samples" / f"{sample}.txt").open("a", encoding="utf-8", newline="") as handle:
            csv.writer(handle, delimiter="\t", lineterminator="\n").writerow((target["variant_id"], by_sample[sample].replace("/", ""), "", target["chromosome"], target["position_bp"]))


def prepare_technical_panel(*, d18_root: Path, configuration_path: Path, workspace_root: Path, run_id: str, created_at: str, scenario: str = SCENARIOS[0]) -> dict[str, Any]:
    """Publie un dérivé neuf et fermé ; les sources D18 restent inchangées."""
    target_root = new_workspace_run(workspace_root, run_id)
    if d18_root.absolute().parent != target_root.parent or scenario not in SCENARIOS:
        raise TechnicalEvaluationError("d19_panel_arguments_invalid")
    config = _json(configuration_path)
    validate_scenario_config(config)
    load_checked_export(d18_root)
    source_manifest = d18_root / "private_truth/d18/blinded_export_manifest.json"
    source_sha, config_sha = _sha(source_manifest), _sha(configuration_path)
    inventory = _json(d18_root / "blinded_inputs/inventory.json")
    temporary = Path(tempfile.mkdtemp(prefix=f".{run_id}.", dir=target_root.parent))
    try:
        private = temporary / "private_truth/d19_inputs"
        private.mkdir(parents=True, mode=0o700)
        shutil.copyfile(configuration_path, private / "configuration.json")
        result_inventory: dict[str, Any] = dict(schema_version="1.0.0", assembly="GRCh38", bundles=[], artifacts=[])
        # Un négatif représentatif suffit au garde déterministe de doublon ;
        # les huit profils sont tous conservés dans le scénario principal.
        bundles = inventory["bundles"] if scenario == SCENARIOS[0] else inventory["bundles"][:1]
        for entry in bundles:
            bundle_id = entry["bundle_id"]
            source = d18_root / "blinded_inputs/bundles" / bundle_id
            bundle = temporary / "blinded_inputs/bundles" / bundle_id
            shutil.copytree(source, bundle)
            metadata = raw_table(bundle / "metadata/samples.tsv")
            validate_samples_master(bundle / "metadata/samples.tsv", bundle / "acpa_samples")
            with _text(bundle / "observed.vcf.gz") as handle:
                panel = read_panel(handle, phased=False)
            target = config["target_metadata"]
            if tuple(row["SAMPLE_ID"] for row in metadata) != panel.samples or any(
                row[2] == target["project_variant_id"] or (int(row[0].removeprefix("chr")), int(row[1])) == (target["chromosome"], target["position_bp"]) for row in panel.variants
            ) or any(int(row[0].removeprefix("chr")) != target["chromosome"] for row in panel.variants):
                raise TechnicalEvaluationError("d19_target_not_independent_of_chip")
            expectation = target_expectation(metadata, config, bundle_id)
            _write_target_files(bundle, private, config, expectation)
            if scenario == SCENARIOS[1]:
                _add_observed_target(bundle, panel, expectation)
            load_pipeline_config(bundle / "pipeline.yaml")
            result_inventory["bundles"].append(dict(bundle_id=bundle_id, sample_count=entry["sample_count"], marker_count=entry["marker_count"] + int(scenario == SCENARIOS[1])))
            for source_artifact in inventory["artifacts"]:
                if source_artifact["bundle_id"] != bundle_id or source_artifact["role"] != "ACPA_SOURCE":
                    continue
                path = temporary / "blinded_inputs" / source_artifact["path"]
                result_inventory["artifacts"].append(dict(bundle_id=bundle_id, role="ACPA_SOURCE", path=source_artifact["path"], sha256=_sha(path), size_bytes=path.stat().st_size))
        manifest = dict(schema_version="1.0.0", method_id=METHOD_ID, manifest_state="FINALIZED", created_at=created_at,
                        source_d18_manifest_sha256=source_sha, configuration_sha256=config_sha, scenario=scenario,
                        inventory=result_inventory, artifacts=[dict(path=path.relative_to(temporary).as_posix(), sha256=_sha(path)) for path in sorted(temporary.rglob("*")) if path.is_file()],
                        calibrated=False, founder_effect_proven=False)
        validate_schema(manifest, "technical_panel_manifest")
        (temporary / MANIFEST_PATH).write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")
        load_checked_panel(temporary)
        validate_derivation(temporary, d18_root)
        load_checked_export(d18_root)
        if _sha(source_manifest) != source_sha or _sha(configuration_path) != config_sha:
            raise TechnicalEvaluationError("d19_panel_source_changed")
        for path in temporary.rglob("*"):
            path.chmod(0o700 if path.is_dir() else 0o600)
        temporary.rename(target_root)
    except Exception:
        shutil.rmtree(temporary)
        raise
    return manifest


def validate_derivation(run_root: Path, d18_root: Path) -> None:
    """Vérifie la copie depuis D18, puis le seul ajout préspecifié du négatif.

    Appelée à la préparation et après clôture pour l'évaluation, jamais par
    le processus pipeline. L'attente cible et sa sérialisation doivent toutes
    deux correspondre à la règle indépendante du registre neutre.
    """
    manifest = load_checked_panel(run_root)
    original = _json(d18_root / "blinded_inputs/inventory.json")
    if _sha(d18_root / "private_truth/d18/blinded_export_manifest.json") != manifest["source_d18_manifest_sha256"]:
        raise TechnicalEvaluationError("d19_derivation_source_binding_invalid")
    config = _json(run_root / "private_truth/d19_inputs/configuration.json")
    target = config["target_metadata"]
    negative = manifest["scenario"] == SCENARIOS[1]
    entries = original["bundles"][:1] if negative else original["bundles"]
    if [entry["bundle_id"] for entry in entries] != [entry["bundle_id"] for entry in manifest["inventory"]["bundles"]]:
        raise TechnicalEvaluationError("d19_derivation_population_invalid")
    for entry, derived in zip(entries, manifest["inventory"]["bundles"], strict=True):
        bundle_id = entry["bundle_id"]
        source = d18_root / "blinded_inputs/bundles" / bundle_id
        bundle = run_root / "blinded_inputs/bundles" / bundle_id
        if derived["sample_count"] != entry["sample_count"] or derived["marker_count"] != entry["marker_count"] + int(negative) or _sha(source / "metadata/samples.tsv") != _sha(bundle / "metadata/samples.tsv"):
            raise TechnicalEvaluationError("d19_derivation_registry_changed")
        metadata = raw_table(bundle / "metadata/samples.tsv")
        expectation = target_expectation(metadata, config, bundle_id)
        if _json(run_root / "private_truth/d19_inputs" / f"{bundle_id}.target.json") != expectation:
            raise TechnicalEvaluationError("d19_derivation_target_truth_invalid")
        input_rows = raw_table(bundle / "target_genotypes.tsv")
        wanted_rows = [dict(SAMPLE_ID=row["sample_id"], GENOTYPE=row["genotype"], GENOTYPE_SOURCE=GENOTYPE_SOURCE, MEASUREMENT_METHOD=MEASUREMENT_METHOD) for row in expectation["genotypes"]]
        if input_rows != wanted_rows or yaml.safe_load((bundle / "target_metadata.yaml").read_text()) != target:
            raise TechnicalEvaluationError("d19_derivation_target_input_invalid")
        for row, wanted in zip(metadata, wanted_rows, strict=True):
            source_bytes = (source / "acpa_samples" / row["SOURCE_FILE"]).read_bytes()
            appended = ("\t".join((target["project_variant_id"], wanted["GENOTYPE"].replace("/", ""), "", str(target["chromosome"]), str(target["position_bp"]))) + "\n").encode() if negative else b""
            if (bundle / "acpa_samples" / row["SOURCE_FILE"]).read_bytes() != source_bytes + appended:
                raise TechnicalEvaluationError("d19_derivation_acpa_changed")
        if not negative:
            if _sha(source / "observed.vcf.gz") != _sha(bundle / "observed.vcf.gz"):
                raise TechnicalEvaluationError("d19_derivation_observation_changed")
            continue
        with _text(source / "observed.vcf.gz") as handle:
            before = read_panel(handle, phased=False)
        with _text(bundle / "observed.vcf.gz") as handle:
            after = read_panel(handle, phased=False)
        retained = [index for index, row in enumerate(after.variants) if row[2] != target["project_variant_id"]]
        added = [index for index, row in enumerate(after.variants) if row[2] == target["project_variant_id"]]
        if after.samples != before.samples or tuple(after.variants[index] for index in retained) != before.variants or not np.array_equal(after.alleles[retained], before.alleles) or len(added) != 1:
            raise TechnicalEvaluationError("d19_derivation_negative_background_changed")
        variant = after.variants[added[0]]
        expected_calls = np.array([[(target["ref"], target["alt"]).index(base) for base in row["GENOTYPE"].split("/")] for row in wanted_rows], dtype=np.int8)
        if variant[:5] != (str(target["chromosome"]), str(target["position_bp"]), target["project_variant_id"], target["ref"], target["alt"]) or not np.array_equal(after.alleles[added[0]], expected_calls):
            raise TechnicalEvaluationError("d19_derivation_negative_target_invalid")


def load_checked_panel(run_root: Path) -> dict[str, Any]:
    """Contrôle inventaire fermé, configurations bornées et toutes les empreintes."""
    manifest = _json(contained_file(run_root, MANIFEST_PATH))
    validate_schema(manifest, "technical_panel_manifest")
    files = {MANIFEST_PATH}
    for artifact in manifest["artifacts"]:
        if artifact["path"] in files:
            raise TechnicalEvaluationError("d19_panel_duplicate_artifact")
        path = contained_file(run_root, artifact["path"])
        if not path.is_file() or _sha(path) != artifact["sha256"]:
            raise TechnicalEvaluationError("d19_panel_checksum_mismatch")
        files.add(artifact["path"])
    directories = {parent.as_posix() for name in files for parent in Path(name).parents if parent != Path(".")}
    if {path.relative_to(run_root).as_posix() for path in run_root.rglob("*")} != files | directories:
        raise TechnicalEvaluationError("d19_panel_tree_mismatch")
    config_path = contained_file(run_root, "private_truth/d19_inputs/configuration.json")
    if _sha(config_path) != manifest["configuration_sha256"]:
        raise TechnicalEvaluationError("d19_panel_configuration_mismatch")
    config = _json(config_path)
    expected_config = technical_configuration(config)
    inventory = manifest["inventory"]
    bundle_ids = [item["bundle_id"] for item in inventory["bundles"]]
    if len(bundle_ids) != len(set(bundle_ids)):
        raise TechnicalEvaluationError("d19_panel_duplicate_bundle")
    for bundle in inventory["bundles"]:
        root = run_root / "blinded_inputs/bundles" / bundle["bundle_id"]
        if load_pipeline_config(root / "pipeline.yaml") != expected_config:
            raise TechnicalEvaluationError("d19_panel_scope_invalid")
        selected = [item for item in inventory["artifacts"] if item["bundle_id"] == bundle["bundle_id"]]
        if len(selected) != bundle["sample_count"] or len({item["path"] for item in selected}) != len(selected):
            raise TechnicalEvaluationError("d19_panel_source_count_invalid")
        for item in selected:
            path = contained_file(run_root, "blinded_inputs/" + item["path"])
            if not item["path"].startswith(f'bundles/{bundle["bundle_id"]}/acpa_samples/') or _sha(path) != item["sha256"] or path.stat().st_size != item["size_bytes"]:
                raise TechnicalEvaluationError("d19_panel_inventory_mismatch")
    if any(item["bundle_id"] not in bundle_ids for item in inventory["artifacts"]):
        raise TechnicalEvaluationError("d19_panel_unknown_bundle")
    return manifest


def main() -> int:
    parser = argparse.ArgumentParser(description="Préparation des scénarios D19 03–04")
    for name in ("d18-root", "configuration-path", "workspace-root"):
        parser.add_argument("--" + name, type=Path, required=True)
    parser.add_argument("--run-id", required=True)
    parser.add_argument("--scenario", choices=SCENARIOS, default=SCENARIOS[0])
    manifest = prepare_technical_panel(**vars(parser.parse_args()), created_at=datetime.now(timezone.utc).isoformat())
    print(json.dumps(dict(bundles=len(manifest["inventory"]["bundles"]), scenario=manifest["scenario"])))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
