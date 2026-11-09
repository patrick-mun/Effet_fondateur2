"""D17 : profils de mesure privés, appariés, audités et publiés atomiquement."""

from __future__ import annotations

import argparse
from collections.abc import Iterator
from contextlib import contextmanager
from datetime import datetime, timezone
import gzip
import io
import json
from pathlib import Path
import platform
import re
import shutil
import tempfile
from typing import Any, TextIO

import jsonschema
import numpy as np

from effet_fondateur.calibration.materialize_ibd import _artifact, _json, _safe_file, _sha, _text
from effet_fondateur.calibration.materialize_roh import _validate_schema as validate_roh_schema
from effet_fondateur.calibration.simulation.ibd import IbdSimulationError
from effet_fondateur.calibration.simulation.roh import RohSimulationError
from effet_fondateur.calibration.simulation.measurement_errors import (
    METHOD_ID, MeasurementError, MeasurementPanel, PerturbedPanel, inject_profile,
    derive_seed, project_dense_truth, quota, read_panel, replay_events, write_panel,
)


def validate_schema(document: object, name: str) -> None:
    """Valide un contrat D17 strict, avec date ISO et sans clé JSON dupliquée."""
    schema = _json(Path(__file__).resolve().parents[3] / "calibration/schemas" / f"{name}.schema.json")
    try:
        jsonschema.Draft202012Validator(schema, format_checker=jsonschema.FormatChecker()).validate(document)
    except jsonschema.ValidationError as error:
        raise MeasurementError(f"d17_schema_invalid:{name}") from error


def validate_configuration(document: dict[str, Any]) -> None:
    """Exige des identifiants uniques et exactement un contrôle sans défaut."""
    validate_schema(document, "measurement_error_config")
    profiles = document["profiles"]
    if len({item["profile_id"] for item in profiles}) != len(profiles):
        raise MeasurementError("d17_duplicate_profile")
    controls = [item for item in profiles if item["marker_retention_fraction"] == 1 and
                all(item[name] == 0 for name in ("missing_rate", "genotyping_error_rate", "phase_point_rate", "phase_switch_rate"))]
    if len(controls) != 1 or controls[0]["profile_id"] != "control":
        raise MeasurementError("d17_one_intact_control_required")


@contextmanager
def deterministic_gzip(path: Path) -> Iterator[TextIO]:
    """Texte UTF-8, gzip sans horodatage ni nom de fichier local."""
    with path.open("wb") as raw:
        with gzip.GzipFile(filename="", fileobj=raw, mode="wb", mtime=0) as zipped:
            with io.TextIOWrapper(zipped, encoding="utf-8", newline="") as handle:
                yield handle


def _audit_exports(
    panel: MeasurementPanel, result: PerturbedPanel, root: Path,
) -> dict[str, Any]:
    """Relit le registre et les deux VCF, puis compare toutes les cellules."""
    with _text(root / "measurement_errors.jsonl.gz") as handle:
        events = [json.loads(line) for line in handle]
    event_schema = _json(Path(__file__).resolve().parents[3] / "calibration/schemas/measurement_error_event.schema.json")
    validator = jsonschema.Draft202012Validator(event_schema)
    for event in events:
        if not validator.is_valid(event):
            raise MeasurementError("d17_event_schema_invalid")
    if events != list(result.events):
        raise MeasurementError("d17_registry_export_mismatch")
    replayed, retained = replay_events(panel, events, profile_id=result.summary["profile_id"])
    if not np.array_equal(replayed, result.alleles) or not np.array_equal(retained, result.retained):
        raise MeasurementError("d17_unregistered_change")
    for filename, phased in (("candidate.phased.vcf.gz", True), ("observed.vcf.gz", False)):
        with _text(root / filename) as handle:
            exported = read_panel(handle, phased=phased)
        expected_variants = tuple((*row[:5], ".", "PASS", ".", "GT") for index, row in enumerate(panel.variants) if retained[index])
        expected = replayed[retained] if phased else np.sort(replayed[retained], axis=2)
        if exported.samples != panel.samples or exported.variants != expected_variants or not np.array_equal(exported.alleles, expected):
            raise MeasurementError("d17_vcf_export_mismatch")
    return dict(status="PASS", unregistered_change_count=0,
                phased_export_mismatch_count=0, unphased_export_mismatch_count=0)


def _write_profile(
    root: Path, panel: MeasurementPanel, result: PerturbedPanel,
) -> dict[str, Any]:
    """Écrit puis audite un profil ; les exports de phase restent privés."""
    root.mkdir(parents=True, mode=0o700)
    for filename, phased in (("candidate.phased.vcf.gz", True), ("observed.vcf.gz", False)):
        with deterministic_gzip(root / filename) as handle:
            write_panel(handle, panel, result.alleles, result.retained, phased=phased)
    with deterministic_gzip(root / "measurement_errors.jsonl.gz") as handle:
        for event in result.events:
            handle.write(json.dumps(event, separators=(",", ":")) + "\n")
    return {**result.summary, "audit": _audit_exports(panel, result, root)}


def _materialize(
    *, d16_root: Path, config_path: Path, workspace_root: Path,
    run_id: str, created_at: str,
) -> dict[str, Any]:
    if not re.fullmatch(r"[A-Za-z0-9][A-Za-z0-9._-]*", run_id):
        raise MeasurementError("d17_run_id_invalid")
    workspace_root = workspace_root.absolute()
    parts = tuple(part.casefold() for part in workspace_root.parts)
    if workspace_root.name != "calibration_workspace" or ".." in parts or "dock6_real" in parts or any(
        left == "data" and right == "input" for left, right in zip(parts, parts[1:])
    ) or any(
        parent.is_symlink() for parent in (workspace_root, *workspace_root.parents)
    ):
        raise MeasurementError("d17_workspace_invalid")
    if d16_root.absolute().parent != workspace_root:
        raise MeasurementError("d17_input_root_invalid")
    workspace_root.mkdir(parents=True, exist_ok=True)
    target = workspace_root / run_id
    if target.exists() or target.is_symlink():
        raise MeasurementError("d17_run_exists")
    config = _json(config_path)
    validate_configuration(config)
    source_path = d16_root / "private_truth/d16/roh_simulation_manifest.json"
    source = _json(source_path)
    validate_roh_schema(source, "roh_simulation_manifest")
    if source["panel_role"] != config["panel_role"]:
        raise MeasurementError("d17_panel_role_mismatch")
    source_sha, config_sha = _sha(source_path), _sha(config_path)
    source_paths = {item["role"]: _artifact(d16_root, source, item["role"]) for item in source["artifacts"]}
    if len(source_paths) != 8:
        raise MeasurementError("d17_source_artifacts_invalid")
    with _text(source_paths["OBSERVED_ROH_PANEL"]) as handle:
        observed = read_panel(handle, phased=False)
    with _text(source_paths["DENSE_PHASED_ROH_PANEL"]) as handle:
        baseline = project_dense_truth(handle, observed)
    if len(baseline.variants) != source["integrity"]["observed_variant_count"]:
        raise MeasurementError("d17_source_marker_count_mismatch")
    domain = source["domain"]
    if any(int(row[0].removeprefix("chr")) != domain["chromosome"] or
           not domain["start_bp"] <= int(row[1]) <= domain["end_bp"] for row in baseline.variants):
        raise MeasurementError("d17_observed_domain_mismatch")
    snapshot = baseline.alleles.copy()
    temporary = Path(tempfile.mkdtemp(prefix=f".{run_id}.", dir=workspace_root))
    try:
        private = temporary / "private_truth/d17"
        private.mkdir(parents=True, mode=0o700)
        artifacts: list[dict[str, Any]] = []

        def artifact(path: Path, role: str, count: int, profile_id: str | None = None) -> None:
            path.chmod(0o600)
            artifacts.append(dict(role=role, profile_id=profile_id, path=path.relative_to(temporary).as_posix(),
                                  sha256=_sha(path), size_bytes=path.stat().st_size, record_count=count))

        shutil.copyfile(config_path, private / "configuration.json")
        artifact(private / "configuration.json", "CONFIGURATION", 1)
        with deterministic_gzip(private / "baseline.phased.vcf.gz") as handle:
            write_panel(handle, baseline, baseline.alleles, np.ones(len(baseline.variants), dtype=bool), phased=True)
        artifact(private / "baseline.phased.vcf.gz", "BASELINE_PHASED_TRUTH", len(baseline.variants))
        profiles = []
        for profile in config["profiles"]:
            result = inject_profile(baseline, profile, parent_seed=config["measurement_error_seed"],
                                    panel_role=config["panel_role"], target_chromosome=domain["chromosome"],
                                    target_bp=domain["target_position_bp"])
            root = private / "profiles" / profile["profile_id"]
            profiles.append(_write_profile(root, baseline, result))
            for filename, role, count in (
                ("observed.vcf.gz", "OBSERVED_MEASUREMENT_PANEL", result.summary["retained_marker_count"]),
                ("candidate.phased.vcf.gz", "CANDIDATE_PHASED_MEASUREMENT_PANEL", result.summary["retained_marker_count"]),
                ("measurement_errors.jsonl.gz", "MEASUREMENT_ERROR_REGISTRY", len(result.events)),
            ):
                artifact(root / filename, role, count, profile["profile_id"])
        if not np.array_equal(snapshot, baseline.alleles):
            raise MeasurementError("d17_baseline_modified")
        if _sha(source_path) != source_sha or _sha(config_path) != config_sha:
            raise MeasurementError("d17_input_modified_during_run")
        for role, path in source_paths.items():
            if _sha(path) != next(item["sha256"] for item in source["artifacts"] if item["role"] == role):
                raise MeasurementError("d17_biological_source_modified")
        manifest = dict(schema_version="1.0.0", manifest_state="FINALIZED",
                        measurement_error_manifest_id=f"{run_id}.d17", created_at=created_at,
                        method_id=METHOD_ID, panel_role=config["panel_role"],
                        scope="REGIONAL_PAIRED_MEASUREMENT_STRESS_TEST", domain=domain,
                        configuration=dict(configuration_id=config["configuration_id"], sha256=config_sha,
                                           measurement_error_seed=config["measurement_error_seed"]),
                        source=dict(stage="D16", manifest_sha256=source_sha,
                                    roh_simulation_manifest_id=source["roh_simulation_manifest_id"],
                                    artifacts=[dict(role=item["role"], sha256=item["sha256"]) for item in source["artifacts"]]),
                        software=dict(python=platform.python_version(), numpy=np.__version__),
                        sample_count=len(baseline.samples), source_marker_count=len(baseline.variants),
                        profiles=profiles, artifacts=artifacts, biological_truth_modified=False,
                        calibrated=False, gamma_inputs_modified=False, founder_effect_proven=False)
        validate_manifest(manifest)
        manifest_path = private / "measurement_error_manifest.json"
        manifest_path.write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")
        manifest_path.chmod(0o600)
        for directory in temporary.rglob("*"):
            if directory.is_dir():
                directory.chmod(0o700)
        temporary.rename(target)
    except Exception:
        shutil.rmtree(temporary)
        raise
    return manifest


def validate_manifest(manifest: dict[str, Any]) -> None:
    """Contrôle schéma, quotas, dénominateurs et bijection profils/artefacts."""
    validate_schema(manifest, "measurement_error_manifest")
    profile_ids = [item["profile_id"] for item in manifest["profiles"]]
    if len(set(profile_ids)) != len(profile_ids):
        raise MeasurementError("d17_duplicate_profile")
    validate_configuration(dict(schema_version="1.0.0", configuration_state="FINALIZED",
                                configuration_id=manifest["configuration"]["configuration_id"],
                                method_id=manifest["method_id"], panel_role=manifest["panel_role"],
                                measurement_error_seed=manifest["configuration"]["measurement_error_seed"],
                                target_policy="KEEP_IF_OBSERVED_NO_GENOTYPE_PROTECTION",
                                profiles=[item["requested"] for item in manifest["profiles"]]))
    expected = {(None, "BASELINE_PHASED_TRUTH"), (None, "CONFIGURATION")} | {
        (profile_id, role) for profile_id in profile_ids for role in (
            "OBSERVED_MEASUREMENT_PANEL", "CANDIDATE_PHASED_MEASUREMENT_PANEL", "MEASUREMENT_ERROR_REGISTRY")}
    actual = {(item["profile_id"], item["role"]) for item in manifest["artifacts"]}
    if expected != actual or len(actual) != len(manifest["artifacts"]) or len({item["path"] for item in manifest["artifacts"]}) != len(actual):
        raise MeasurementError("d17_artifact_inventory_invalid")
    for item in manifest["profiles"]:
        rates, counts = item["requested"], item["event_counts"]
        if rates["profile_id"] != item["profile_id"] or item["source_marker_count"] != manifest["source_marker_count"] or item["sample_count"] != manifest["sample_count"]:
            raise MeasurementError("d17_profile_identity_or_count_mismatch")
        if item["retained_cell_count"] != item["retained_marker_count"] * item["sample_count"]:
            raise MeasurementError("d17_cell_denominator_invalid")
        if not 0 <= item["eligible_switch_boundary_count"] <= item["heterozygote_count"] <= item["retained_cell_count"]:
            raise MeasurementError("d17_phase_denominator_invalid")
        if item["retained_marker_count"] != max(int(item["target_observed"]), quota(rates["marker_retention_fraction"], item["source_marker_count"])):
            raise MeasurementError("d17_density_quota_mismatch")
        for operation, seed in item["seeds"].items():
            if seed != derive_seed(manifest["configuration"]["measurement_error_seed"], manifest["panel_role"], item["profile_id"], operation):
                raise MeasurementError("d17_derived_seed_mismatch")
        if counts["MARKER_DROP"] != item["source_marker_count"] - item["retained_marker_count"] or item["event_count"] != sum(counts.values()):
            raise MeasurementError("d17_event_count_invalid")
        for kind, rate, denominator in (
            ("GENOTYPE", "genotyping_error_rate", "retained_cell_count"),
            ("MISSING_GT", "missing_rate", "retained_cell_count"),
            ("PHASE_POINT", "phase_point_rate", "heterozygote_count"),
            ("PHASE_SWITCH_BOUNDARY", "phase_switch_rate", "eligible_switch_boundary_count"),
        ):
            if counts[kind] != quota(rates[rate], item[denominator]):
                raise MeasurementError("d17_quota_mismatch")
    for artifact in manifest["artifacts"]:
        if artifact["profile_id"] is None:
            expected_count = 1 if artifact["role"] == "CONFIGURATION" else manifest["source_marker_count"]
        else:
            profile = next(item for item in manifest["profiles"] if item["profile_id"] == artifact["profile_id"])
            expected_count = profile["event_count"] if artifact["role"] == "MEASUREMENT_ERROR_REGISTRY" else profile["retained_marker_count"]
        if artifact["record_count"] != expected_count:
            raise MeasurementError("d17_artifact_record_count_mismatch")


def materialize_measurement_errors(
    *, d16_root: Path, config_path: Path, workspace_root: Path,
    run_id: str, created_at: str,
) -> dict[str, Any]:
    """Produit D17 depuis D16 sans appeler le pipeline ni lire de données d'étude."""
    try:
        return _materialize(d16_root=d16_root, config_path=config_path, workspace_root=workspace_root,
                            run_id=run_id, created_at=created_at)
    except (IbdSimulationError, RohSimulationError, jsonschema.ValidationError) as error:
        raise MeasurementError("d17_input_integrity_error") from error


def load_checked_run(run_root: Path) -> dict[str, Any]:
    """Relit un manifest et vérifie tous ses artefacts sans recalcul de métrique."""
    manifest = _json(run_root / "private_truth/d17/measurement_error_manifest.json")
    validate_manifest(manifest)
    for item in manifest["artifacts"]:
        path = run_root / item["path"]
        _safe_file(path)
        if run_root.resolve() not in path.resolve().parents or _sha(path) != item["sha256"] or path.stat().st_size != item["size_bytes"]:
            raise MeasurementError("d17_artifact_hash_or_path_invalid")
    configuration = next(item for item in manifest["artifacts"] if item["role"] == "CONFIGURATION")
    config = _json(run_root / configuration["path"])
    validate_configuration(config)
    if (configuration["sha256"] != manifest["configuration"]["sha256"] or
        config["profiles"] != [item["requested"] for item in manifest["profiles"]] or
        config["panel_role"] != manifest["panel_role"] or
        config["measurement_error_seed"] != manifest["configuration"]["measurement_error_seed"] or
        config["configuration_id"] != manifest["configuration"]["configuration_id"]):
        raise MeasurementError("d17_configuration_binding_mismatch")
    return manifest


def main() -> int:
    """CLI : affiche uniquement effectifs et statut technique agrégés."""
    parser = argparse.ArgumentParser(description="Profils de mesure D17")
    for name in ("d16-root", "config-path", "workspace-root"):
        parser.add_argument("--" + name, type=Path, required=True)
    parser.add_argument("--run-id", required=True)
    result = materialize_measurement_errors(**vars(parser.parse_args()), created_at=datetime.now(timezone.utc).isoformat())
    print(json.dumps(dict(profiles=len(result["profiles"]), samples=result["sample_count"],
                         markers=result["source_marker_count"], status="PASS")))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
