"""D18 : exports ACPA/VCF observés, registre neutre et inventaire sans vérité.

Le seul contenu génotypique consommé est OBSERVED_MEASUREMENT_PANEL D17.
Les empreintes des autres artefacts D17 sont contrôlées sans ouvrir les
registres ni les candidats phasés. La correspondance des identifiants et des
profils reste privée. Le brouillon de configuration n'active que 00–02 : ni
approbation humaine fictive, ni annotation moléculaire CONFIRMED inventée.
"""

from __future__ import annotations

import argparse
import csv
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
import random
import re
import shutil
import tempfile
from typing import Any

import jsonschema
import numpy as np
from numpy.typing import NDArray
import yaml

from effet_fondateur.calibration.materialize_ibd import _json, _safe_file, _sha, _text
from effet_fondateur.calibration.materialize_measurement_errors import deterministic_gzip, load_checked_run
from effet_fondateur.calibration.simulation.ibd import IbdSimulationError
from effet_fondateur.calibration.simulation.measurement_errors import MeasurementError, MeasurementPanel, read_panel, write_panel
from effet_fondateur.contracts import load_pipeline_config
from effet_fondateur.contracts.samples import validate_samples_master
from effet_fondateur.stages.validate_sources import inspect_source


METHOD_ID = "observed_only_acpa_blinding_v1"
ACPA_COLUMNS = ("Probe Set ID", "Forward Strand Base Calls", "dbSNP RS ID", "Chromosome", "Chromosomal Position")
SAMPLE_COLUMNS = ("SAMPLE_ID", "SOURCE_FILE", "FID", "IID", "PID", "MID", "SEX", "CLINICAL_STATUS", "GROUP_LABEL", "TARGET_GENOTYPE", "TARGET_GENOTYPE_SOURCE", "ARRAY_BATCH", "INCLUDE_GENOMEWIDE", "INCLUDE_TARGET_CHROMOSOME", "NOTES_CODE")


class BlindedExportError(ValueError):
    """Violation D18 sans donnée individuelle dans le message."""


def validate_schema(document: object, name: str) -> None:
    """Valide un contrat D18 sans accepter de champ supplémentaire."""
    schema = _json(Path(__file__).resolve().parents[3] / "calibration/schemas" / f"{name}.schema.json")
    try:
        jsonschema.Draft202012Validator(schema, format_checker=jsonschema.FormatChecker()).validate(document)
    except jsonschema.ValidationError as error:
        raise BlindedExportError(f"d18_schema_invalid:{name}") from error


def _opaque_id(config: dict[str, Any], prefix: str, *parts: str) -> str:
    """Pseudonyme déterministe à domaine distinct ; pas une garantie cryptographique."""
    payload = f'{METHOD_ID}:{config["panel_role"]}:{config["blinding_seed"]}:' + ":".join(parts)
    return prefix + hashlib.sha256(payload.encode()).hexdigest()[:16]


def blind_panel(panel: MeasurementPanel, config: dict[str, Any], bundle_id: str) -> tuple[MeasurementPanel, list[dict[str, str]]]:
    """Permute et renomme les colonnes sans changer une seule cellule observée."""
    labels = [_opaque_id(config, "S", bundle_id, sample) for sample in panel.samples]
    if len(set(labels)) != len(labels):
        raise BlindedExportError("d18_pseudonym_collision")
    order = list(range(len(labels)))
    seed = int(_opaque_id(config, "B", bundle_id, "order")[1:], 16)
    random.Random(seed).shuffle(order)
    mapping = [dict(source_sample_id=panel.samples[index], blinded_sample_id=labels[index]) for index in order]
    result = MeasurementPanel(tuple(labels[index] for index in order), panel.variants, panel.alleles[:, order, :].copy())
    return result, mapping


def pipeline_configuration(target: dict[str, Any]) -> dict[str, Any]:
    """Profil 00–02 strictement régional, sans filiation ni approbation prétendue."""
    return dict(
        schema_version="1.0.0", project=dict(project_id="input_validation", assembly="GRCh38"),
        target=dict(gene="SYNTHETIC_LOCUS", chromosome=target["chromosome"], position_bp=target["position_bp"],
                    ref=target["ref"], alt=target["alt"], transcript=None, project_variant_id="analysis_target", rsid=None),
        inputs=dict(acpa_samples_dir="acpa_samples", sample_metadata="metadata/samples.tsv",
                    target_variant_metadata=None, target_variant_genotypes=None, genetic_map=None,
                    genetic_map_catalog=None, reference_panel_catalog=None),
        tools=dict(plink=None, king=None, bcftools=None, rscript=None, phasing_adapter=None, local_ibd_adapter=None, explicit_ibd_adapters=None),
        policies=dict(overwrite_inputs=False, allow_clinical_genotype_inference=False,
                      pseudonymize_figures=True, stop_on_critical_warning=True),
        stages=dict(initialize_run=dict(enabled=True, parameters={}),
                    validate_sources=dict(enabled=True, parameters=dict(allowed_extensions=[".txt"], required_autosomes=[target["chromosome"]])),
                    build_sample_registry=dict(enabled=True, parameters={})))


def _base_call(pair: NDArray[np.int8], ref: str, alt: str) -> str:
    return "NoCall" if pair[0] < 0 else "".join((ref, alt)[int(allele)] for allele in pair)


def _write_bundle(root: Path, panel: MeasurementPanel, config: dict[str, Any]) -> None:
    """Sérialise uniquement les colonnes autorisées, sans recopier d'en-tête source."""
    sources = root / "acpa_samples"
    sources.mkdir(parents=True, mode=0o700)
    (root / "metadata").mkdir(mode=0o700)
    with deterministic_gzip(root / "observed.vcf.gz") as handle:
        write_panel(handle, panel, panel.alleles, np.ones(len(panel.variants), dtype=bool), phased=False)
    for sample_index, sample in enumerate(panel.samples):
        with (sources / f"{sample}.txt").open("w", encoding="utf-8", newline="") as handle:
            handle.write("# Array Type Name: CytoScan 750K Accel Array\n# UCSC Genomic Version: hg38\n")
            writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
            writer.writerow(ACPA_COLUMNS)
            for index, row in enumerate(panel.variants):
                writer.writerow((row[2], _base_call(panel.alleles[index, sample_index], row[3], row[4]),
                                 "", row[0].removeprefix("chr"), row[1]))
    with (root / "metadata/samples.tsv").open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(SAMPLE_COLUMNS)
        for sample in panel.samples:
            # FID individuel = unité technique, pas une affirmation d'absence
            # de parenté. PID/MID=0 signifient filiation non renseignée.
            writer.writerow((sample, f"{sample}.txt", sample, sample, "0", "0", "UNKNOWN", "UNKNOWN", "UNSPECIFIED",
                             "", "", "", "false", "true", ""))
    (root / "pipeline.yaml").write_text(yaml.safe_dump(pipeline_configuration(config["target"]), sort_keys=False), encoding="utf-8")


def audit_bundle(root: Path, expected: MeasurementPanel, config: dict[str, Any]) -> dict[str, int]:
    """Relit VCF, registre, configuration et chaque appel ACPA avant publication.

    Les parseurs de structure 01 et de registre 02 sont ceux du pipeline. Une
    seconde voie décode chaque ACPA en allèles et contrôle le manque explicite.
    Aucun PLINK ni pipeline complet n'est exécuté.
    """
    with _text(root / "observed.vcf.gz") as handle:
        exported = read_panel(handle, phased=False)
    expected_variants = tuple((*row[:5], ".", "PASS", ".", "GT") for row in expected.variants)
    if exported.samples != expected.samples or exported.variants != expected_variants or not np.array_equal(exported.alleles, expected.alleles):
        raise BlindedExportError("d18_vcf_roundtrip_mismatch")
    metadata = validate_samples_master(root / "metadata/samples.tsv", root / "acpa_samples")
    if tuple(row["SAMPLE_ID"] for row in metadata.rows) != expected.samples:
        raise BlindedExportError("d18_sample_order_mismatch")
    with (root / "metadata/samples.tsv").open(encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    for row, sample in zip(rows, expected.samples, strict=True):
        expected_values = (sample, f"{sample}.txt", sample, sample, "0", "0", "UNKNOWN", "UNKNOWN", "UNSPECIFIED", "", "", "", "false", "true", "")
        if tuple(row[name] for name in SAMPLE_COLUMNS) != expected_values:
            raise BlindedExportError("d18_non_neutral_metadata")
    if load_pipeline_config(root / "pipeline.yaml") != pipeline_configuration(config["target"]):
        raise BlindedExportError("d18_pipeline_config_mismatch")
    missing = 0
    for sample_index, sample in enumerate(expected.samples):
        path = root / "acpa_samples" / f"{sample}.txt"
        inspection = inspect_source(path, path.name, _sha(path), "GRCh38", str(config["target"]["chromosome"]),
                                    (str(config["target"]["chromosome"]),), (".txt",))
        if not inspection.passed or inspection.data_row_count != len(expected.variants):
            raise BlindedExportError("d18_pipeline_source_validation_failed")
        with path.open(encoding="utf-8") as handle:
            if next(handle) != "# Array Type Name: CytoScan 750K Accel Array\n" or next(handle) != "# UCSC Genomic Version: hg38\n":
                raise BlindedExportError("d18_acpa_header_mismatch")
            reader = csv.DictReader(handle, delimiter="\t")
            if tuple(reader.fieldnames or ()) != ACPA_COLUMNS:
                raise BlindedExportError("d18_acpa_columns_invalid")
            actual = list(reader)
        if len(actual) != len(expected.variants):
            raise BlindedExportError("d18_acpa_row_count_mismatch")
        for index, (row, variant) in enumerate(zip(actual, expected.variants, strict=True)):
            if set(row) != set(ACPA_COLUMNS) or any(value is None for value in row.values()):
                raise BlindedExportError("d18_acpa_extra_field")
            calls = row["Forward Strand Base Calls"]
            if row["Probe Set ID"] != variant[2] or row["Chromosome"] != variant[0].removeprefix("chr") or row["Chromosomal Position"] != variant[1] or row["dbSNP RS ID"]:
                raise BlindedExportError("d18_acpa_locus_mismatch")
            if calls == "NoCall":
                decoded = (-1, -1)
                missing += 1
            else:
                try:
                    decoded = tuple(sorted((variant[3], variant[4]).index(base) for base in calls))
                except ValueError as error:
                    raise BlindedExportError("d18_acpa_allele_invalid") from error
            if tuple(expected.alleles[index, sample_index]) != decoded:
                raise BlindedExportError("d18_acpa_genotype_mismatch")
    return dict(sample_count=len(expected.samples), marker_count=len(expected.variants),
                missing_call_count=missing, genotype_mismatch_count=0, truth_leak_count=0)


def _artifact(path: Path, root: Path, bundle_id: str, role: str, count: int) -> dict[str, Any]:
    return dict(bundle_id=bundle_id, role=role, path=path.relative_to(root).as_posix(),
                sha256=_sha(path), size_bytes=path.stat().st_size, record_count=count)


def validate_inventory(document: dict[str, Any]) -> None:
    """Vérifie la bijection inventaire/bundles et les chemins exacts autorisés."""
    validate_schema(document, "blinded_input_inventory")
    bundles = {item["bundle_id"]: item for item in document["bundles"]}
    artifacts = document["artifacts"]
    if len(bundles) != len(document["bundles"]) or len({item["path"] for item in artifacts}) != len(artifacts):
        raise BlindedExportError("d18_inventory_duplicate")
    for bundle_id, bundle in bundles.items():
        selected = [item for item in artifacts if item["bundle_id"] == bundle_id]
        roles = [item["role"] for item in selected]
        if len(selected) != bundle["sample_count"] + 3 or any(roles.count(role) != 1 for role in ("OBSERVED_VCF", "SAMPLE_METADATA", "PIPELINE_CONFIG")) or roles.count("ACPA_SOURCE") != bundle["sample_count"]:
            raise BlindedExportError("d18_inventory_role_mismatch")
        for item in selected:
            if not item["path"].startswith(f"bundles/{bundle_id}/"):
                raise BlindedExportError("d18_inventory_bundle_mismatch")
            suffix = item["path"].removeprefix(f"bundles/{bundle_id}/")
            required = {"OBSERVED_VCF": "observed.vcf.gz", "PIPELINE_CONFIG": "pipeline.yaml", "SAMPLE_METADATA": "metadata/samples.tsv"}
            if item["role"] in required and suffix != required[item["role"]] or item["role"] == "ACPA_SOURCE" and not suffix.startswith("acpa_samples/"):
                raise BlindedExportError("d18_inventory_role_path_mismatch")
            count = bundle["sample_count"] if item["role"] == "SAMPLE_METADATA" else 1 if item["role"] == "PIPELINE_CONFIG" else bundle["marker_count"]
            if item["record_count"] != count:
                raise BlindedExportError("d18_inventory_record_count_mismatch")
    if any(item["bundle_id"] not in bundles for item in artifacts):
        raise BlindedExportError("d18_inventory_unknown_bundle")


def _validate_file_tree(root: Path, inventory: dict[str, Any]) -> None:
    """Refuse fichiers, répertoires et liens non inventoriés, même vides."""
    files = {root / item["path"] for item in inventory["artifacts"]} | {root / "inventory.json"}
    directories = {root / "bundles"}
    for bundle in inventory["bundles"]:
        parent = root / "bundles" / bundle["bundle_id"]
        directories.update((parent, parent / "metadata", parent / "acpa_samples"))
    if set(root.rglob("*")) != files | directories or any(path.is_symlink() for path in (root, *root.rglob("*"))):
        raise BlindedExportError("d18_uninventoried_file")


def _validate_mapping(mapping: list[dict[str, Any]], inventory: dict[str, Any]) -> None:
    """Exige une bijection privée par profil et aucun pseudonyme inter-lots partagé."""
    validate_schema(mapping, "blinding_map")
    bundles = {item["bundle_id"]: item["sample_count"] for item in inventory["bundles"]}
    if len(mapping) != len(bundles) or {item["bundle_id"] for item in mapping} != set(bundles) or len({item["profile_id"] for item in mapping}) != len(mapping):
        raise BlindedExportError("d18_mapping_bundle_mismatch")
    all_labels = []
    for item in mapping:
        source_ids = [sample["source_sample_id"] for sample in item["samples"]]
        labels = [sample["blinded_sample_id"] for sample in item["samples"]]
        if len(source_ids) != bundles[item["bundle_id"]] or len(set(source_ids)) != len(source_ids) or set(source_ids) & set(labels):
            raise BlindedExportError("d18_mapping_sample_mismatch")
        all_labels.extend(labels)
    if len(set(all_labels)) != len(all_labels):
        raise BlindedExportError("d18_mapping_pseudonym_duplicate")


def _export(
    *, d17_root: Path, config_path: Path, workspace_root: Path,
    run_id: str, created_at: str,
) -> dict[str, Any]:
    workspace_root = workspace_root.absolute()
    parts = tuple(part.casefold() for part in workspace_root.parts)
    if workspace_root.name != "calibration_workspace" or ".." in parts or "dock6_real" in parts or any(
        left == "data" and right == "input" for left, right in zip(parts, parts[1:])
    ) or any(parent.is_symlink() for parent in (workspace_root, *workspace_root.parents)):
        raise BlindedExportError("d18_workspace_invalid")
    if not re.fullmatch(r"[A-Za-z0-9][A-Za-z0-9._-]*", run_id) or d17_root.absolute().parent != workspace_root:
        raise BlindedExportError("d18_run_or_input_root_invalid")
    target = workspace_root / run_id
    if target.exists() or target.is_symlink():
        raise BlindedExportError("d18_run_exists")
    config = _json(config_path)
    validate_schema(config, "blinded_export_config")
    source = load_checked_run(d17_root)
    if config["panel_role"] != source["panel_role"] or config["blinding_seed"] == source["configuration"]["measurement_error_seed"]:
        raise BlindedExportError("d18_role_or_seed_invalid")
    if config["target"]["ref"] == config["target"]["alt"] or any(config["target"][key] != source["domain"][field] for key, field in (("chromosome", "chromosome"), ("position_bp", "target_position_bp"))):
        raise BlindedExportError("d18_target_mismatch")
    manifest_path = d17_root / "private_truth/d17/measurement_error_manifest.json"
    source_sha, config_sha = _sha(manifest_path), _sha(config_path)
    temporary = Path(tempfile.mkdtemp(prefix=f".{run_id}.", dir=workspace_root))
    try:
        blind_root, private = temporary / "blinded_inputs", temporary / "private_truth/d18"
        blind_root.mkdir(mode=0o700)
        private.mkdir(parents=True, mode=0o700)
        shutil.copyfile(config_path, private / "configuration.json")
        inventory: dict[str, Any] = dict(schema_version="1.0.0", scope="REGIONAL_SOURCE_VALIDATION_ONLY", assembly=config["assembly"], bundles=[], artifacts=[])
        mapping = []
        audits = []
        for profile in sorted(source["profiles"], key=lambda item: _opaque_id(config, "B", item["profile_id"])):
            profile_id = profile["profile_id"]
            bundle_id = _opaque_id(config, "B", profile_id)
            if any(item["bundle_id"] == bundle_id for item in inventory["bundles"]):
                raise BlindedExportError("d18_bundle_collision")
            artifact = next(item for item in source["artifacts"] if item["profile_id"] == profile_id and item["role"] == "OBSERVED_MEASUREMENT_PANEL")
            with _text(d17_root / artifact["path"]) as handle:
                panel = read_panel(handle, phased=False)
            if len(panel.samples) != profile["sample_count"] or len(panel.variants) != profile["retained_marker_count"]:
                raise BlindedExportError("d18_source_count_mismatch")
            if any(int(row[0].removeprefix("chr")) != config["target"]["chromosome"] for row in panel.variants):
                raise BlindedExportError("d18_source_scope_mismatch")
            for row in panel.variants:
                if not re.fullmatch(r"[A-Za-z0-9][A-Za-z0-9._-]*", row[2]):
                    raise BlindedExportError("d18_probe_identifier_invalid")
                if int(row[1]) == config["target"]["position_bp"] and (row[3], row[4]) != (config["target"]["ref"], config["target"]["alt"]):
                    raise BlindedExportError("d18_target_alleles_mismatch")
            blinded, sample_mapping = blind_panel(panel, config, bundle_id)
            bundle_root = blind_root / "bundles" / bundle_id
            _write_bundle(bundle_root, blinded, config)
            audits.append(dict(bundle_id=bundle_id, **audit_bundle(bundle_root, blinded, config)))
            mapping.append(dict(profile_id=profile_id, bundle_id=bundle_id, source_observed_sha256=artifact["sha256"], samples=sample_mapping))
            inventory["bundles"].append(dict(bundle_id=bundle_id, sample_count=len(panel.samples), marker_count=len(panel.variants)))
            for path, role, count in ((bundle_root / "observed.vcf.gz", "OBSERVED_VCF", len(panel.variants)),
                                      (bundle_root / "metadata/samples.tsv", "SAMPLE_METADATA", len(panel.samples)),
                                      (bundle_root / "pipeline.yaml", "PIPELINE_CONFIG", 1)):
                inventory["artifacts"].append(_artifact(path, blind_root, bundle_id, role, count))
            for sample in blinded.samples:
                inventory["artifacts"].append(_artifact(bundle_root / "acpa_samples" / f"{sample}.txt", blind_root, bundle_id, "ACPA_SOURCE", len(panel.variants)))
        validate_inventory(inventory)
        _validate_mapping(mapping, inventory)
        (blind_root / "inventory.json").write_text(json.dumps(inventory, indent=2) + "\n", encoding="utf-8")
        _validate_file_tree(blind_root, inventory)
        (private / "blinding_map.json").write_text(json.dumps(mapping, indent=2) + "\n", encoding="utf-8")
        if _sha(manifest_path) != source_sha or _sha(config_path) != config_sha:
            raise BlindedExportError("d18_source_modified")
        load_checked_run(d17_root)
        manifest = dict(schema_version="1.0.0", manifest_state="FINALIZED", method_id=METHOD_ID,
                        blinded_export_manifest_id=f"{run_id}.d18", created_at=created_at,
                        panel_role=config["panel_role"], scope="REGIONAL_SOURCE_VALIDATION_ONLY",
                        source_manifest_sha256=source_sha, configuration_sha256=config_sha,
                        blinding_seed=config["blinding_seed"],
                        inventory_sha256=_sha(blind_root / "inventory.json"),
                        blinding_map_sha256=_sha(private / "blinding_map.json"),
                        bundle_count=len(inventory["bundles"]), file_count=len(inventory["artifacts"]),
                        audits=audits, biological_truth_modified=False,
                        calibrated=False, founder_effect_proven=False)
        validate_schema(manifest, "blinded_export_manifest")
        (private / "blinded_export_manifest.json").write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")
        for path in temporary.rglob("*"):
            path.chmod(0o700 if path.is_dir() else 0o600)
        temporary.rename(target)
    except Exception:
        shutil.rmtree(temporary)
        raise
    return manifest


def export_blinded_inputs(
    *, d17_root: Path, config_path: Path, workspace_root: Path, run_id: str, created_at: str,
) -> dict[str, Any]:
    """Exporte D18 dans un run neuf ; aucun candidat phasé n'est décodé."""
    try:
        return _export(d17_root=d17_root, config_path=config_path, workspace_root=workspace_root,
                       run_id=run_id, created_at=created_at)
    except (IbdSimulationError, MeasurementError) as error:
        raise BlindedExportError("d18_source_integrity_error") from error


def load_checked_export(run_root: Path) -> dict[str, Any]:
    """Vérifie manifest, inventaire, correspondance privée et absence d'extra."""
    private = run_root / "private_truth/d18"
    manifest = _json(private / "blinded_export_manifest.json")
    validate_schema(manifest, "blinded_export_manifest")
    for path, sha in ((run_root / "blinded_inputs/inventory.json", manifest["inventory_sha256"]),
                      (private / "blinding_map.json", manifest["blinding_map_sha256"]),
                      (private / "configuration.json", manifest["configuration_sha256"])):
        _safe_file(path)
        if _sha(path) != sha:
            raise BlindedExportError("d18_manifest_binding_mismatch")
    inventory = _json(run_root / "blinded_inputs/inventory.json")
    validate_inventory(inventory)
    config = _json(private / "configuration.json")
    validate_schema(config, "blinded_export_config")
    if any(config[field] != manifest[field] for field in ("panel_role", "method_id", "blinding_seed")):
        raise BlindedExportError("d18_configuration_binding_mismatch")
    mapping = json.loads((private / "blinding_map.json").read_text(encoding="utf-8"))
    _validate_mapping(mapping, inventory)
    root = run_root / "blinded_inputs"
    _validate_file_tree(root, inventory)
    for item in mapping:
        path = root / "bundles" / item["bundle_id"] / "metadata/samples.tsv"
        with path.open(encoding="utf-8") as handle:
            actual_labels = [row["SAMPLE_ID"] for row in csv.DictReader(handle, delimiter="\t")]
        if actual_labels != [sample["blinded_sample_id"] for sample in item["samples"]]:
            raise BlindedExportError("d18_mapping_registry_mismatch")
    for item in inventory["artifacts"]:
        path = root / item["path"]
        _safe_file(path)
        if _sha(path) != item["sha256"] or path.stat().st_size != item["size_bytes"]:
            raise BlindedExportError("d18_export_artifact_modified")
    if manifest["bundle_count"] != len(inventory["bundles"]) or manifest["file_count"] != len(inventory["artifacts"]):
        raise BlindedExportError("d18_inventory_count_mismatch")
    if sorted((item["bundle_id"], item["sample_count"], item["marker_count"]) for item in manifest["audits"]) != sorted(
        (item["bundle_id"], item["sample_count"], item["marker_count"]) for item in inventory["bundles"]
    ):
        raise BlindedExportError("d18_audit_inventory_mismatch")
    return manifest


def main() -> int:
    """CLI D18, effectifs agrégés uniquement."""
    parser = argparse.ArgumentParser(description="Export aveugle D18")
    for name in ("d17-root", "config-path", "workspace-root"):
        parser.add_argument("--" + name, type=Path, required=True)
    parser.add_argument("--run-id", required=True)
    manifest = export_blinded_inputs(**vars(parser.parse_args()), created_at=datetime.now(timezone.utc).isoformat())
    print(json.dumps(dict(bundles=manifest["bundle_count"], files=manifest["file_count"], status="PASS")))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
