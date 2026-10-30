"""D19 : évaluation des artefacts 00–04 après reçu de clôture empreinté.

Aucune exécution n'est déclenchée par cet évaluateur. Un reçu manquant rend
l'étape non évaluable ; une sortie annoncée réussie mais absente ou corrompue
est un défaut FAIL. Les attentes de conversion viennent du VCF observé D18.
Les attentes de cible 04 sont optionnelles, privées et fournies séparément.
"""

from __future__ import annotations

import argparse
import csv
from datetime import datetime, timezone
import json
from pathlib import Path, PurePosixPath
import re
import shutil
import tempfile
from typing import Any

import jsonschema
import numpy as np

from effet_fondateur.calibration.export_blinded_inputs import load_checked_export
from effet_fondateur.calibration.materialize_ibd import _json, _safe_file, _sha, _text
from effet_fondateur.calibration.simulation.measurement_errors import read_panel
from effet_fondateur.calibration.technical_evaluation import (
    BASES, METHOD_ID, STAGES, GenotypeMatrix, TechnicalEvaluationError,
    compare_genotypes, evaluate_target, observation_matrix, read_plink,
    result, row_result, unavailable,
)
from effet_fondateur.contracts import load_pipeline_config, validate_json_document, validate_tsv_table


STAGE_NAMES = dict(zip(STAGES, ("initialize_run", "validate_sources", "build_sample_registry", "convert_acpa", "prepare_target_variant_dataset"), strict=True))


def validate_schema(document: object, name: str) -> None:
    """Contrats D19 stricts ; les sorties ne contiennent que des agrégats."""
    schema = _json(Path(__file__).resolve().parents[3] / "calibration/schemas" / f"{name}.schema.json")
    try:
        jsonschema.Draft202012Validator(schema, format_checker=jsonschema.FormatChecker()).validate(document)
    except jsonschema.ValidationError as error:
        raise TechnicalEvaluationError(f"d19_schema_invalid_{name}") from error


def contained_file(root: Path, relative: str) -> Path:
    """Confine avant lecture, y compris quand le fichier déclaré est absent."""
    relative_path = PurePosixPath(relative)
    if relative_path.is_absolute() or ".." in relative_path.parts or "\\" in relative:
        raise TechnicalEvaluationError("d19_unsafe_path")
    path = root / relative
    lowered = [part.casefold() for part in path.absolute().parts]
    if "dock6_real" in lowered or any(left == "data" and right == "input" for left, right in zip(lowered, lowered[1:])) or any(
        parent.is_symlink() for parent in (path, *path.parents)
    ) or not path.resolve().is_relative_to(root.resolve()) or path.exists() and (not path.is_file() or path.stat().st_nlink != 1):
        raise TechnicalEvaluationError("d19_unsafe_path")
    return path


def raw_table(path: Path) -> list[dict[str, str]]:
    """Lit des chaînes canoniques après contrôle structurel des en-têtes/lignes."""
    with path.open(encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if not reader.fieldnames or len(set(reader.fieldnames)) != len(reader.fieldnames):
            raise TechnicalEvaluationError("d19_table_header_invalid")
        rows = list(reader)
    if any(None in row or any(value is None for value in row.values()) for row in rows):
        raise TechnicalEvaluationError("d19_table_columns_invalid")
    return rows


def _stage_artifacts(run_root: Path, manifest: dict[str, Any], stage: str) -> dict[str, Path] | None:
    """Contrôle manifest → audit/sorties → artefacts, identités et schémas."""
    matches = [item for item in manifest["stages"] if item["stage_id"] == stage]
    if not matches:
        return None
    if len(matches) != 1:
        raise TechnicalEvaluationError("d19_duplicate_stage")
    entry = matches[0]
    if entry["stage_name"] != STAGE_NAMES[stage]:
        raise TechnicalEvaluationError("d19_stage_identity_mismatch")
    if entry["state"] in {"PENDING", "SKIPPED", "BLOCKED"}:
        return None
    if entry["state"] not in {"SUCCEEDED", "CACHED"}:
        raise TechnicalEvaluationError("d19_stage_not_successful")
    directory = f"stages/{stage}_{STAGE_NAMES[stage]}"
    paths = {"audit": contained_file(run_root, entry["audit_path"] or "invalid"),
             "outputs": contained_file(run_root, directory + "/stage_outputs.json")}
    for kind, schema in (("audit", "stage_audit.schema.json"), ("outputs", "stage_outputs.schema.json")):
        path = paths[kind]
        if not path.is_file() or _sha(path) != entry["audit_sha256" if kind == "audit" else "stage_outputs_sha256"]:
            raise TechnicalEvaluationError("d19_stage_checksum_mismatch")
        document = _json(path)
        validate_json_document(document, schema)
        if any(document[key] != expected for key, expected in (
            ("run_id", manifest["run_id"]), ("stage_id", stage),
            ("stage_name", STAGE_NAMES[stage]), ("signature", entry["signature"]))):
            raise TechnicalEvaluationError("d19_stage_identity_mismatch")
    outputs, audit = _json(paths["outputs"]), _json(paths["audit"])
    if audit["outputs"] != outputs["artifacts"]:
        raise TechnicalEvaluationError("d19_audit_output_binding_mismatch")
    artifacts: dict[str, Path] = {}
    for artifact in outputs["artifacts"]:
        if artifact["artifact_id"] in artifacts or artifact["producer_signature"] != entry["signature"] or artifact["producer_stage"] != f"{stage}_{STAGE_NAMES[stage]}" or not artifact["path"].startswith(directory + "/"):
            raise TechnicalEvaluationError("d19_artifact_identity_invalid")
        path = contained_file(run_root, artifact["path"])
        if not path.is_file() or _sha(path) != artifact["sha256"]:
            raise TechnicalEvaluationError("d19_artifact_checksum_mismatch")
        if artifact["schema_name"]:
            schema_path = Path(__file__).resolve().parents[3] / "schemas" / artifact["schema_name"]
            if Path(artifact["schema_name"]).name != artifact["schema_name"]:
                raise TechnicalEvaluationError("d19_artifact_schema_path_invalid")
            schema = _json(schema_path)
            if "x-tsv" in schema:
                validate_tsv_table(path, artifact["schema_name"])
            else:
                validate_json_document(_json(path), artifact["schema_name"])
        artifacts[artifact["artifact_id"]] = path
    return artifacts


def _required(artifacts: dict[str, Path], name: str) -> Path:
    if name not in artifacts:
        raise TechnicalEvaluationError("d19_required_artifact_missing")
    return artifacts[name]


def _dataset(artifacts: dict[str, Path], prefix: str) -> GenotypeMatrix:
    """Vérifie aussi le descripteur PLINK ; les hachages ne remplacent pas les GT."""
    descriptor_path = _required(artifacts, prefix + "_dataset")
    descriptor = _json(descriptor_path)
    validate_json_document(descriptor, "plink_dataset.schema.json")
    for suffix in ("bed", "bim", "fam"):
        path = _required(artifacts, prefix + "_" + suffix)
        if descriptor["files"][suffix] != dict(name=path.name, sha256=_sha(path)):
            raise TechnicalEvaluationError("d19_dataset_descriptor_checksum_mismatch")
    bed = _required(artifacts, prefix + "_bed")
    if any(_required(artifacts, prefix + "_" + suffix) != bed.with_suffix("." + suffix) for suffix in ("bim", "fam")):
        raise TechnicalEvaluationError("d19_dataset_file_prefix_mismatch")
    matrix = read_plink(bed.with_suffix(""))
    if descriptor["sample_count"] != len(matrix.samples) or descriptor["variant_count"] != len(matrix.variants) or descriptor["assembly"] != "GRCh38":
        raise TechnicalEvaluationError("d19_dataset_descriptor_count_mismatch")
    return matrix


def _target_truth(path: Path, bundle_id: str, metadata: list[dict[str, str]]) -> GenotypeMatrix:
    document = _json(path)
    validate_schema(document, "technical_target_expectation")
    if document["bundle_id"] != bundle_id:
        raise TechnicalEvaluationError("d19_target_truth_bundle_mismatch")
    rows = document["genotypes"]
    if len({row["sample_id"] for row in rows}) != len(rows) or {row["sample_id"] for row in rows} != {row["SAMPLE_ID"] for row in metadata}:
        raise TechnicalEvaluationError("d19_target_truth_sample_mismatch")
    by_sample = {row["sample_id"]: row["genotype"] for row in rows}
    target = document["target"]
    if target["ref"] == target["alt"]:
        raise TechnicalEvaluationError("d19_target_truth_alleles_invalid")
    pairs = []
    for row in metadata:
        genotype = by_sample[row["SAMPLE_ID"]]
        if any(base not in (target["ref"], target["alt"]) for base in genotype.split("/")):
            raise TechnicalEvaluationError("d19_target_truth_alleles_invalid")
        pairs.append(sorted(BASES.index(base) for base in genotype.split("/")))
    return GenotypeMatrix(tuple((row["FID"], row["IID"]) for row in metadata),
                          ((target["chromosome"], target["position_bp"], target["variant_id"]),),
                          np.array([pairs], dtype=np.int8))


def compare_fam(path: Path, metadata: list[dict[str, str]], *, stage: str) -> dict[str, Any]:
    """Contrôle les six champs FAM, sans confondre parent inconnu et indépendance.

    Codes PLINK : sexe 1/2/0 et phénotype 1/2/-9 pour sain/atteint/inconnu.
    Aucune valeur n'est déduite d'une appartenance au panel.
    """
    columns = ("FID", "IID", "PID", "MID", "sex", "phenotype")
    wanted = [dict(zip(columns, (row["FID"], row["IID"], row["PID"], row["MID"],
                               {"MALE": "1", "FEMALE": "2"}.get(row["SEX"], "0"),
                               {"UNAFFECTED": "1", "AFFECTED": "2"}.get(row["CLINICAL_STATUS"], "-9")), strict=True)) for row in metadata]
    rows = [line.split() for line in path.read_text().splitlines() if line.strip()]
    if any(len(row) != 6 for row in rows):
        raise TechnicalEvaluationError("d19_fam_columns_invalid")
    actual = [dict(zip(columns, row, strict=True)) for row in rows]
    return row_result(stage, wanted, actual, keys=("FID", "IID"), fields=columns)


def evaluate_pipeline_run(
    *, bundle_root: Path, inventory: dict[str, Any], bundle_id: str,
    run_root: Path | None, receipt: dict[str, Any] | None,
    target_truth_path: Path | None = None,
    expected_target_rejection: bool = False,
) -> list[dict[str, Any]]:
    """Produit cinq verdicts indépendants, sans exécuter ni modifier le pipeline."""
    if run_root is None or receipt is None:
        return [unavailable(stage, "pipeline_run_not_provided") for stage in STAGES]
    manifest_path = contained_file(run_root, "manifest.json")
    if not manifest_path.is_file() or _sha(manifest_path) != receipt["manifest_sha256"]:
        return [result(stage, {}, ["closed_run_manifest_checksum_mismatch"]) for stage in STAGES]
    try:
        manifest = _json(manifest_path)
        validate_json_document(manifest, "run_manifest.schema.json")
    except ValueError:
        return [result(stage, {}, ["run_manifest_schema_invalid"]) for stage in STAGES]
    if any(item["state"] == "RUNNING" for item in manifest["stages"]):
        raise TechnicalEvaluationError("d19_run_not_closed")
    expected_config = load_pipeline_config(bundle_root / "pipeline.yaml")
    metadata = raw_table(bundle_root / "metadata/samples.tsv")
    with _text(bundle_root / "observed.vcf.gz") as handle:
        panel = read_panel(handle, phased=False)
    expected = observation_matrix(panel, metadata)
    verdicts = []
    for stage in STAGES:
        try:
            if stage == "04" and expected_target_rejection:
                from effet_fondateur.calibration.technical_rejection import evaluate_target_rejection
                if not expected_config["stages"].get(STAGE_NAMES[stage], {}).get("enabled", False):
                    verdicts.append(result(stage, {}, ["executed_stage_not_authorized_by_configuration"]))
                    continue
                verdicts.append(evaluate_target_rejection(run_root=run_root, manifest=manifest, receipt=receipt,
                                                          expected=expected, metadata=metadata, target=expected_config["target"]))
                continue
            artifacts = _stage_artifacts(run_root, manifest, stage)
            if artifacts is None:
                verdicts.append(unavailable(stage, "stage_not_executed"))
                continue
            if not expected_config["stages"].get(STAGE_NAMES[stage], {}).get("enabled", False):
                verdicts.append(result(stage, {}, ["executed_stage_not_authorized_by_configuration"]))
                continue
            if stage == "00":
                config_path = contained_file(run_root, "config.resolved.yaml")
                environment_path = contained_file(run_root, "environment.json")
                checksums = int(not config_path.is_file() or _sha(config_path) != manifest["config_sha256"])
                checksums += int(not environment_path.is_file() or _sha(environment_path) != receipt["environment_sha256"])
                if checksums:
                    verdicts.append(result(stage, {"checksum_mismatch_count": checksums}, ["initialization_checksum_mismatch"]))
                    continue
                config = load_pipeline_config(config_path)
                environment = _json(environment_path)
                mismatch = int(config != expected_config)
                summary = _json(_required(artifacts, "initialization_summary"))
                expected_summary = dict(schema_version="1.0.0", project_id=config["project"]["project_id"], assembly=config["project"]["assembly"],
                                        enabled_stage_count=sum(bool(item["enabled"]) for item in config["stages"].values()),
                                        target_definition_complete=all(config["target"][key] is not None for key in ("chromosome", "position_bp", "ref", "alt", "project_variant_id")))
                summary_mismatch = int(summary != expected_summary)
                environment_invalid = int(not environment.get("python", {}).get("version", "").startswith("3.12.") or not environment.get("pipeline_version"))
                verdicts.append(result(stage, dict(checksum_mismatch_count=0, configuration_mismatch_count=mismatch,
                                                  initialization_summary_mismatch_count=summary_mismatch,
                                                  environment_validation_error_count=environment_invalid),
                                       [name for name, value in (("configuration_mismatch", mismatch), ("environment_invalid", environment_invalid), ("initialization_summary_mismatch", summary_mismatch)) if value]))
            elif stage == "01":
                expected_rows = []
                for artifact in inventory["artifacts"]:
                    if artifact["bundle_id"] == bundle_id and artifact["role"] == "ACPA_SOURCE":
                        expected_rows.append(dict(SOURCE_FILE=Path(artifact["path"]).name, SIZE_BYTES=str(artifact["size_bytes"]),
                                                  SHA256=artifact["sha256"], ENCODING="utf-8-sig", ARRAY_TYPE="CytoScan 750K Accel Array",
                                                  GENOMIC_BUILD="hg38", DATA_ROW_COUNT=str(len(panel.variants)),
                                                  AUTOSOMES_PRESENT=",".join(str(value) for value in sorted({int(row[0].removeprefix("chr")) for row in panel.variants})),
                                                  TARGET_CHROMOSOME_PRESENT=str(any(int(row[0].removeprefix("chr")) == expected_config["target"]["chromosome"] for row in panel.variants)).lower(), STATUS="PASS"))
                expected_rows.sort(key=lambda row: row["SOURCE_FILE"])
                actual = raw_table(_required(artifacts, "source_inventory"))
                verdicts.append(row_result(stage, expected_rows, actual, keys=("SOURCE_FILE",), fields=tuple(expected_rows[0])))
            elif stage == "02":
                actual = raw_table(_required(artifacts, "samples_master"))
                verdicts.append(row_result(stage, metadata, actual, keys=("SAMPLE_ID",), fields=tuple(metadata[0])))
            elif stage == "03":
                prefixes = ("genomewide_base", "target_chromosome_base")
                comparisons = [compare_genotypes(expected, _dataset(artifacts, prefix)) for prefix in prefixes]
                families = [compare_fam(_required(artifacts, prefix + "_fam"), metadata, stage=stage) for prefix in prefixes]
                metrics = {f"{scope}_{name}": value for scope, comparison in zip(("regional_base", "target_chromosome"), comparisons, strict=True) for name, value in comparison["metrics"].items()}
                metrics.update({f"{scope}_fam_{name}": value for scope, comparison in zip(("regional_base", "target_chromosome"), families, strict=True) for name, value in comparison["metrics"].items()})
                alignment = raw_table(_required(artifacts, "sample_alignment"))
                wanted = [dict(ROW_INDEX=str(index), **{key: row[key] for key in ("SAMPLE_ID", "FID", "IID", "SOURCE_FILE")}) for index, row in enumerate(metadata, 1)]
                alignment_result = row_result(stage, wanted, alignment, keys=("SAMPLE_ID",), fields=tuple(wanted[0]))
                reasons = sorted({reason for item in (*comparisons, *families, alignment_result) for reason in item["reasons"]})
                metrics.update({f"alignment_{name}": value for name, value in alignment_result["metrics"].items()})
                verdicts.append(result(stage, metrics, reasons))
            else:
                if target_truth_path is None:
                    verdicts.append(unavailable(stage, "independent_target_truth_not_provided"))
                    continue
                target_document = _json(target_truth_path)
                target = _target_truth(target_truth_path, bundle_id, metadata)
                target_config = expected_config["target"]
                if target.variants[0] != (target_config["chromosome"], target_config["position_bp"], target_config["project_variant_id"]) or any(
                    target_document["target"][key] != target_config[key] for key in ("ref", "alt")
                ):
                    raise TechnicalEvaluationError("d19_target_truth_locus_mismatch")
                comparison = evaluate_target(expected, target, _dataset(artifacts, "target_variant"))
                family = compare_fam(_required(artifacts, "target_variant_fam"), metadata, stage=stage)
                metrics = {**comparison["metrics"], **{f"fam_{name}": value for name, value in family["metrics"].items()}}
                verdicts.append(result(stage, metrics, sorted(set(comparison["reasons"] + family["reasons"]))))
        except TechnicalEvaluationError as error:
            if str(error) == "d19_unsafe_path":
                raise
            verdicts.append(result(stage, {}, [str(error)]))
        except (ValueError, KeyError, OSError):
            verdicts.append(result(stage, {}, ["stage_contract_or_artifact_invalid"]))
    return verdicts


def evaluate_technical(
    *, d18_root: Path, run_index_path: Path, workspace_root: Path,
    run_id: str, created_at: str,
    technical_panel_root: Path | None = None,
) -> dict[str, Any]:
    """Écrit cinq évaluations distinctes et un manifest agrégé dans un run neuf."""
    workspace_root = workspace_root.absolute()
    parts = tuple(part.casefold() for part in workspace_root.parts)
    if workspace_root.name != "calibration_workspace" or ".." in parts or "dock6_real" in parts or any(
        left == "data" and right == "input" for left, right in zip(parts, parts[1:])
    ) or any(path.is_symlink() for path in (workspace_root, *workspace_root.parents)):
        raise TechnicalEvaluationError("d19_workspace_invalid")
    if d18_root.absolute().parent != workspace_root or not re.fullmatch(r"[A-Za-z0-9][A-Za-z0-9._-]*", run_id):
        raise TechnicalEvaluationError("d19_run_or_source_root_invalid")
    target_root = workspace_root / run_id
    if target_root.exists() or target_root.is_symlink():
        raise TechnicalEvaluationError("d19_run_exists")
    _safe_file(run_index_path)
    if not run_index_path.resolve().is_relative_to(workspace_root.resolve()):
        raise TechnicalEvaluationError("d19_run_index_outside_workspace")
    index = _json(run_index_path)
    index_sha = _sha(run_index_path)
    validate_schema(index, "technical_run_index")
    if len({entry["bundle_id"] for entry in index["runs"]}) != len(index["runs"]) or len({entry["run_path"] for entry in index["runs"]}) != len(index["runs"]):
        raise TechnicalEvaluationError("d19_run_index_duplicate")
    # Refuser un processus encore actif avant tout accès aux correspondances D18.
    for entry in index["runs"]:
        path = contained_file(workspace_root, entry["run_path"] + "/manifest.json")
        if path.is_file() and _sha(path) == entry["manifest_sha256"]:
            try:
                closed_manifest = _json(path)
                validate_json_document(closed_manifest, "run_manifest.schema.json")
            except ValueError:
                continue  # L'évaluateur publiera FAIL pour le contrat invalide.
            if any(item["state"] == "RUNNING" for item in closed_manifest["stages"]):
                raise TechnicalEvaluationError("d19_run_not_closed")
    load_checked_export(d18_root)
    d18_manifest_path = d18_root / "private_truth/d18/blinded_export_manifest.json"
    d18_sha = _sha(d18_manifest_path)
    if index["d18_manifest_sha256"] != d18_sha:
        raise TechnicalEvaluationError("d19_export_binding_mismatch")
    inventory = _json(d18_root / "blinded_inputs/inventory.json")
    input_root = d18_root
    panel_manifest = None
    if technical_panel_root is not None:
        from effet_fondateur.calibration.technical_scenarios import MANIFEST_PATH, load_checked_panel, validate_derivation
        if technical_panel_root.absolute().parent != workspace_root:
            raise TechnicalEvaluationError("d19_technical_panel_root_invalid")
        panel_manifest = load_checked_panel(technical_panel_root)
        validate_derivation(technical_panel_root, d18_root)
        if panel_manifest["source_d18_manifest_sha256"] != d18_sha or index.get("technical_panel_manifest_sha256") != _sha(technical_panel_root / MANIFEST_PATH):
            raise TechnicalEvaluationError("d19_technical_panel_binding_mismatch")
        input_root = technical_panel_root
        inventory = panel_manifest["inventory"]
    elif "technical_panel_manifest_sha256" in index:
        raise TechnicalEvaluationError("d19_technical_panel_required")
    if {entry["bundle_id"] for entry in index["runs"]} - {entry["bundle_id"] for entry in inventory["bundles"]}:
        raise TechnicalEvaluationError("d19_unknown_bundle")
    by_bundle = {entry["bundle_id"]: entry for entry in index["runs"]}
    outputs = {stage: [] for stage in STAGES}
    for bundle in inventory["bundles"]:
        bundle_id = bundle["bundle_id"]
        receipt = by_bundle.get(bundle_id)
        run_root = contained_file(workspace_root, receipt["run_path"] + "/manifest.json").parent if receipt else None
        # Les fichiers réellement lus doivent rester immuables pendant la comparaison.
        snapshot = _run_snapshot(run_root) if run_root else {}
        target_truth = None
        if receipt and receipt["target_truth_path"] is not None:
            if "private_truth" not in PurePosixPath(receipt["target_truth_path"]).parts:
                raise TechnicalEvaluationError("d19_target_truth_must_be_private")
            target_truth = contained_file(workspace_root, receipt["target_truth_path"])
            if panel_manifest is not None and target_truth.absolute() != (input_root / "private_truth/d19_inputs" / f"{bundle_id}.target.json").absolute():
                raise TechnicalEvaluationError("d19_technical_target_truth_binding_mismatch")
            if not target_truth.is_file() or _sha(target_truth) != receipt["target_truth_sha256"]:
                raise TechnicalEvaluationError("d19_target_truth_checksum_mismatch")
        verdicts = evaluate_pipeline_run(bundle_root=input_root / "blinded_inputs/bundles" / bundle_id,
                                         inventory=inventory, bundle_id=bundle_id,
                                         run_root=run_root, receipt=receipt, target_truth_path=target_truth,
                                         expected_target_rejection=bool(panel_manifest and panel_manifest["scenario"] == "TARGET_PRESENT_REJECTION"))
        for verdict in verdicts:
            outputs[verdict["stage_id"]].append(dict(bundle_id=bundle_id, **verdict))
        if run_root and snapshot != _run_snapshot(run_root):
            raise TechnicalEvaluationError("d19_pipeline_changed_during_evaluation")
        if target_truth and _sha(target_truth) != receipt["target_truth_sha256"]:
            raise TechnicalEvaluationError("d19_target_truth_changed_during_evaluation")
    load_checked_export(d18_root)
    if panel_manifest is not None:
        load_checked_panel(input_root)
        if _sha(input_root / MANIFEST_PATH) != index["technical_panel_manifest_sha256"]:
            raise TechnicalEvaluationError("d19_technical_panel_changed")
    if _sha(run_index_path) != index_sha or _sha(d18_manifest_path) != d18_sha:
        raise TechnicalEvaluationError("d19_sources_changed_during_evaluation")
    temporary = Path(tempfile.mkdtemp(prefix=f".{run_id}.", dir=workspace_root))
    try:
        private = temporary / "private_truth/d19"
        private.mkdir(parents=True, mode=0o700)
        artifacts = []
        stage_summaries = []
        for stage, evaluations in outputs.items():
            document = dict(schema_version="1.0.0", method_id=METHOD_ID, report_id=f"CAL-REP-{stage}", evaluations=evaluations)
            validate_schema(document, "technical_evaluation")
            path = private / f"evaluation_{stage}.json"
            path.write_text(json.dumps(document, indent=2) + "\n", encoding="utf-8")
            artifacts.append(dict(stage_id=stage, path=path.relative_to(temporary).as_posix(), sha256=_sha(path)))
            statuses = [item["status"] for item in evaluations]
            aggregate = "FAIL" if "FAIL" in statuses else "NOT_EVALUABLE" if "NOT_EVALUABLE" in statuses else "PASS"
            stage_summaries.append(dict(stage_id=stage, status=aggregate, pass_count=statuses.count("PASS"),
                                        fail_count=statuses.count("FAIL"), not_evaluable_count=statuses.count("NOT_EVALUABLE")))
        manifest = dict(schema_version="1.0.0", manifest_state="FINALIZED", method_id=METHOD_ID,
                        technical_evaluation_manifest_id=f"{run_id}.d19", created_at=created_at,
                        d18_manifest_sha256=d18_sha, run_index_sha256=index_sha,
                        bundle_count=len(inventory["bundles"]), stages=stage_summaries, artifacts=artifacts,
                        scope="REGIONAL_TECHNICAL_PARTIAL_PIPELINE", calibrated=False, founder_effect_proven=False)
        if panel_manifest is not None:
            panel_config = _json(input_root / "private_truth/d19_inputs/configuration.json")
            public_target = {key: panel_config["target_metadata"][key] for key in ("gene", "chromosome", "position_bp", "ref", "alt", "transcript", "hgvs_c", "hgvs_p", "rsid")}
            public_target.update(source_url=panel_config["source"]["url"], source_id=panel_config["source"]["source_id"], genotype_rule=panel_config["genotype_rule"], decision_id=panel_config["manual_approval"]["decision_id"])
            manifest.update(technical_panel_manifest_sha256=index["technical_panel_manifest_sha256"], scenario=panel_manifest["scenario"], technical_target=public_target)
        validate_schema(manifest, "technical_evaluation_manifest")
        (private / "technical_evaluation_manifest.json").write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")
        for path in temporary.rglob("*"):
            path.chmod(0o700 if path.is_dir() else 0o600)
        temporary.rename(target_root)
    except Exception:
        shutil.rmtree(temporary)
        raise
    return manifest


def _run_snapshot(run_root: Path) -> dict[str, str]:
    """Empreinte fermée d'un run diagnostique ; aucun contenu individuel publié."""
    snapshot = {}
    for path in sorted(run_root.rglob("*")):
        if path.is_symlink():
            raise TechnicalEvaluationError("d19_unsafe_path")
        if path.is_file():
            relative = path.relative_to(run_root).as_posix()
            snapshot[relative] = _sha(contained_file(run_root, relative))
    return snapshot


def load_checked_evaluation(run_root: Path) -> dict[str, Any]:
    """Contrôle les cinq résultats, leurs agrégats et leur population de lots."""
    manifest = _json(contained_file(run_root, "private_truth/d19/technical_evaluation_manifest.json"))
    validate_schema(manifest, "technical_evaluation_manifest")
    if [item["stage_id"] for item in manifest["stages"]] != list(STAGES) or [item["stage_id"] for item in manifest["artifacts"]] != list(STAGES):
        raise TechnicalEvaluationError("d19_evaluation_stage_set_invalid")
    bundles = None
    for artifact, summary in zip(manifest["artifacts"], manifest["stages"], strict=True):
        stage = artifact["stage_id"]
        if artifact["path"] != f"private_truth/d19/evaluation_{stage}.json":
            raise TechnicalEvaluationError("d19_evaluation_artifact_path_invalid")
        path = contained_file(run_root, artifact["path"])
        if not path.is_file() or _sha(path) != artifact["sha256"]:
            raise TechnicalEvaluationError("d19_evaluation_checksum_mismatch")
        document = _json(path)
        validate_schema(document, "technical_evaluation")
        if document["report_id"] != f"CAL-REP-{stage}":
            raise TechnicalEvaluationError("d19_evaluation_report_binding_invalid")
        evaluations = document["evaluations"]
        population = [item["bundle_id"] for item in evaluations]
        if len(set(population)) != manifest["bundle_count"] or len(population) != manifest["bundle_count"] or bundles is not None and population != bundles:
            raise TechnicalEvaluationError("d19_evaluation_population_invalid")
        bundles = population
        statuses = [item["status"] for item in evaluations]
        aggregate = "FAIL" if "FAIL" in statuses else "NOT_EVALUABLE" if "NOT_EVALUABLE" in statuses else "PASS"
        if summary != dict(stage_id=stage, status=aggregate, pass_count=statuses.count("PASS"), fail_count=statuses.count("FAIL"), not_evaluable_count=statuses.count("NOT_EVALUABLE")):
            raise TechnicalEvaluationError("d19_evaluation_summary_mismatch")
        if any((item["status"] == "PASS") != (not item["reasons"]) for item in evaluations):
            raise TechnicalEvaluationError("d19_evaluation_reason_status_invalid")
    return manifest


def main() -> int:
    parser = argparse.ArgumentParser(description="Évaluateurs techniques D19, lecture seule du pipeline")
    for name in ("d18-root", "run-index-path", "workspace-root"):
        parser.add_argument("--" + name, type=Path, required=True)
    parser.add_argument("--run-id", required=True)
    parser.add_argument("--technical-panel-root", type=Path)
    manifest = evaluate_technical(**vars(parser.parse_args()), created_at=datetime.now(timezone.utc).isoformat())
    print(json.dumps(manifest["stages"]))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
