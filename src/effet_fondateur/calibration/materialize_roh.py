"""D16 : matérialisation atomique des ROH et conservation des vérités D15."""

from __future__ import annotations

import argparse
import csv
from dataclasses import asdict
from datetime import datetime, timezone
import gzip
import hashlib
import io
import json
from pathlib import Path
import re
import shutil
import tempfile
from typing import Any

import jsonschema

# Les lectures empreintées sont identiques à celles de D15. Aucune logique
# scientifique d'injection IBD n'est exécutée lors de leur réutilisation.
from effet_fondateur.calibration.materialize_ibd import (
    _artifact, _json, _safe_file, _sha, _text,
)
from effet_fondateur.calibration.simulation.founder_recombination import parse_genetic_map
from effet_fondateur.calibration.simulation.population_structure import parse_observed_variant_keys
from effet_fondateur.calibration.simulation.ibd import IbdSimulationError
from effet_fondateur.calibration.simulation.roh import (
    METHOD_ID, RohInjection, RohSimulationError, build_roh_injections,
    derive_roh_ibd_edges, transform_roh_vcf,
    validate_roh_density,
)


def _validate_schema(document: object, name: str) -> None:
    schema_path = Path(__file__).resolve().parents[3] / "calibration/schemas" / f"{name}.schema.json"
    try:
        jsonschema.Draft202012Validator(
            _json(schema_path), format_checker=jsonschema.FormatChecker()
        ).validate(document)
    except jsonschema.ValidationError as error:
        raise RohSimulationError(f"d16_schema_invalid:{name}") from error


def _load_inputs(
    *, d15_root: Path, d14_root: Path, d11a_root: Path,
    config_path: Path, workspace_root: Path,
) -> tuple[dict[str, Any], dict[str, dict[str, Any]], dict[str, Path]]:
    """Valide configuration, manifests, rôles et liaisons avant les génotypes."""
    for root in (d15_root, d14_root, d11a_root):
        if root.absolute().parent != workspace_root:
            raise RohSimulationError("d16_input_root_invalid")
    paths = {
        "D15": d15_root / "private_truth/d15/ibd_simulation_manifest.json",
        "D14": d14_root / "private_truth/d14/founder_simulation_manifest.json",
        "D11A": d11a_root / "private_truth/d11/base_haplotype_manifest.json",
    }
    config = _json(config_path)
    _validate_schema(config, "roh_simulation_config")
    documents = {key: _json(path) for key, path in paths.items()}
    _validate_schema(documents["D15"], "ibd_simulation_manifest")
    if any(document.get("manifest_state") != "FINALIZED" for document in documents.values()):
        raise RohSimulationError("d16_inputs_not_finalized")
    for key in ("D14", "D11A"):
        if documents["D15"]["inputs"][key] != _sha(paths[key]):
            raise RohSimulationError("d16_lineage_binding_mismatch")
    role = config["panel_role"]
    if documents["D15"]["panel_role"] != role or any(
        documents[key]["panel"]["role"] != role for key in ("D14", "D11A")
    ):
        raise RohSimulationError("d16_panel_role_mismatch")
    for key in ("D14", "D11A"):
        if any(documents[key]["domain"][field] != documents["D15"]["domain"][field]
               for field in ("chromosome", "start_bp", "end_bp")):
            raise RohSimulationError("d16_domain_mismatch")
    return config, documents, paths


def _read_plan_inputs(
    *, d15_root: Path, d14_root: Path, d11a_root: Path,
    documents: dict[str, dict[str, Any]],
) -> dict[str, Any]:
    """Charge la géométrie et les vérités ; les données restent synthétiques."""
    d15, d14, d11 = (documents[key] for key in ("D15", "D14", "D11A"))
    dense = _artifact(d15_root, d15, "DENSE_PHASED_IBD_PANEL")
    observed = _artifact(d15_root, d15, "OBSERVED_IBD_PANEL")
    truth_path = _artifact(d15_root, d15, "IBD_INJECTION_TRUTH")
    ibd_truth = json.loads(truth_path.read_text(encoding="utf-8"))
    _validate_schema(ibd_truth, "ibd_injection_truth")
    pedigree_path = _artifact(d14_root, d14, "PEDIGREE_TRUTH")
    carriers_path = _artifact(d14_root, d14, "FOUNDER_CARRIER_TRUTH")
    transmissions_path = _artifact(d14_root, d14, "MENDELIAN_TRANSMISSION_TRUTH")
    map_path = _artifact(d11a_root, d11, "GENETIC_MAP")
    with _text(pedigree_path) as handle:
        pedigree = list(csv.DictReader(handle, delimiter="\t"))
    with _text(carriers_path) as handle:
        carriers = {row["sample_id"] for row in csv.DictReader(handle, delimiter="\t")}
    participants = {sample for item in ibd_truth for sample in
                    (item["sample_1_id"], item["sample_2_id"], item["template_id"]) if sample}
    eligible = [row["sample_id"] for row in pedigree
                if row["family_unit_id"].startswith("UNRELATED_")
                and row["sample_id"] not in carriers | participants]
    with _text(map_path) as handle:
        points = parse_genetic_map(handle, d15["domain"]["chromosome"])
    with _text(observed) as handle:
        keys = parse_observed_variant_keys(handle)
    with _text(transmissions_path) as handle:
        transmissions = tuple((row["child_id"], row["parent_id"],
                               int(row["parental_haplotype_index"]),
                               int(row["child_haplotype_index"]))
                              for row in csv.DictReader(handle, delimiter="\t"))
    return {"dense_path": dense, "ibd_truth_path": truth_path,
            "ibd_truth": ibd_truth, "eligible": eligible, "map_points": points,
            "observed_keys": keys, "transmissions": transmissions,
            "pedigree_path": pedigree_path, "transmissions_path": transmissions_path}


def _write_panels(
    private: Path, plan_inputs: dict[str, Any], injections: tuple[RohInjection, ...],
) -> dict[str, int]:
    """Écrit deux flux gzip reproductibles avec contrôle intégré des copies."""
    with _text(plan_inputs["dense_path"]) as source, (
        private / "roh_panel.phased.vcf.gz"
    ).open("wb") as raw_dense, (
        private / "observed_roh_panel.vcf.gz"
    ).open("wb") as raw_observed:
        with gzip.GzipFile(filename="", fileobj=raw_dense, mode="wb", mtime=0) as zipped_dense, gzip.GzipFile(filename="", fileobj=raw_observed, mode="wb", mtime=0) as zipped_observed:
            with io.TextIOWrapper(zipped_dense, encoding="utf-8", newline="") as dense, io.TextIOWrapper(zipped_observed, encoding="utf-8", newline="") as observed:
                return transform_roh_vcf(
                    source, dense, observed, injections=injections,
                    observed_keys=plan_inputs["observed_keys"],
                    transmissions=plan_inputs["transmissions"],
                    ibd_truth=plan_inputs["ibd_truth"],
                )


def _materialize_roh(
    *, d15_root: Path, d14_root: Path, d11a_root: Path, config_path: Path,
    workspace_root: Path, run_id: str, created_at: str,
) -> dict[str, Any]:
    """Produit un run D16 neuf avec vérités ROH/IBD distinctes et contrôlées.

    Une copie locale H1→H2 réalise l'autozygotie synthétique sans redéfinir le
    pedigree global. Les sujets concernés ne participent à aucune transmission.
    Les bornes et seuils observables sont fixés dans la configuration versionnée.
    """
    if not re.fullmatch(r"[A-Za-z0-9][A-Za-z0-9._-]*", run_id):
        raise RohSimulationError("d16_run_id_invalid")
    workspace_root = workspace_root.absolute()
    if workspace_root.name != "calibration_workspace" or any(
        path.is_symlink() for path in (workspace_root, *workspace_root.parents)
    ):
        raise RohSimulationError("d16_workspace_invalid")
    # Le garde commun refuse aussi les liens, traversées et chemins d'étude.
    _safe_file(config_path)
    workspace_root.mkdir(parents=True, exist_ok=True)
    target_root = workspace_root / run_id
    if target_root.exists() or target_root.is_symlink():
        raise RohSimulationError("d16_run_exists")
    config, documents, input_paths = _load_inputs(
        d15_root=d15_root, d14_root=d14_root, d11a_root=d11a_root,
        config_path=config_path, workspace_root=workspace_root,
    )
    plan_inputs = _read_plan_inputs(d15_root=d15_root, d14_root=d14_root,
                                  d11a_root=d11a_root, documents=documents)
    parent_seed = documents["D15"]["seed"]["parent_recombination_seed"]
    role = config["panel_role"]
    derived_seed = int.from_bytes(hashlib.sha256(
        f"{METHOD_ID}:D16:{parent_seed}:{role}".encode()
    ).digest()[:4], "big")
    injections = build_roh_injections(
        scenarios=config["scenarios"], eligible_samples=plan_inputs["eligible"],
        ibd_truth=plan_inputs["ibd_truth"], map_points=plan_inputs["map_points"],
        chromosome=documents["D15"]["domain"]["chromosome"],
        target_bp=documents["D15"]["domain"]["target_position_bp"],
        random_seed=derived_seed,
    )
    if any(item.is_positive and item.length_bp < config["minimum_length_kb"] * 1000
           for item in injections):
        raise RohSimulationError("d16_length_below_configured_minimum")
    density = validate_roh_density(
        injections, [key[1] for key in plan_inputs["observed_keys"]],
        minimum_markers=config["minimum_observed_markers"],
        maximum_density_kb=config["maximum_density_kb_per_marker"],
        maximum_gap_kb=config["maximum_internal_gap_kb"],
        minimum_flanks=config["minimum_flank_markers_each_side"],
    )
    induced_ibd = derive_roh_ibd_edges(
        injections, plan_inputs["ibd_truth"], plan_inputs["map_points"],
        replication_id=run_id,
    )
    _validate_schema(induced_ibd, "roh_induced_ibd_truth")
    temporary = Path(tempfile.mkdtemp(prefix=f".{run_id}.", dir=workspace_root))
    try:
        private = temporary / "private_truth/d16"
        private.mkdir(parents=True, mode=0o700)
        integrity = _write_panels(private, plan_inputs, injections)
        for field in ("dense_variant_count", "observed_variant_count"):
            if integrity[field] != documents["D15"]["integrity"][field]:
                raise RohSimulationError("d16_variant_count_mismatch")
        truth = [{"replication_id": run_id, **asdict(item)} for item in injections]
        _validate_schema(truth, "roh_injection_truth")
        (private / "roh_truth.json").write_text(json.dumps(truth, indent=2) + "\n", encoding="utf-8")
        (private / "observability_audit.json").write_text(json.dumps(density, indent=2) + "\n", encoding="utf-8")
        (private / "roh_induced_ibd_truth.json").write_text(json.dumps(induced_ibd, indent=2) + "\n", encoding="utf-8")
        for key, filename in (("ibd_truth_path", "ibd_truth.json"),
                              ("pedigree_path", "pedigree_truth.tsv.gz"),
                              ("transmissions_path", "mendelian_transmission_truth.tsv.gz")):
            shutil.copyfile(plan_inputs[key], private / filename)
        artifact_specs = (
            ("roh_panel.phased.vcf.gz", "DENSE_PHASED_ROH_PANEL", integrity["dense_variant_count"]),
            ("observed_roh_panel.vcf.gz", "OBSERVED_ROH_PANEL", integrity["observed_variant_count"]),
            ("roh_truth.json", "ROH_INJECTION_TRUTH", len(truth)),
            ("observability_audit.json", "ROH_OBSERVABILITY_AUDIT", len(density)),
            ("roh_induced_ibd_truth.json", "ROH_INDUCED_IBD_TRUTH", len(induced_ibd)),
            ("ibd_truth.json", "PRESERVED_IBD_TRUTH", len(plan_inputs["ibd_truth"])),
            ("pedigree_truth.tsv.gz", "PRESERVED_PEDIGREE_TRUTH", documents["D14"]["panel"]["individual_count"]),
            ("mendelian_transmission_truth.tsv.gz", "PRESERVED_TRANSMISSION_TRUTH", len(plan_inputs["transmissions"])),
        )
        artifacts = []
        for filename, artifact_role, count in artifact_specs:
            path = private / filename
            path.chmod(0o600)
            artifacts.append({"role": artifact_role, "path": path.relative_to(temporary).as_posix(),
                              "sha256": _sha(path), "size_bytes": path.stat().st_size,
                              "record_count": count})
        manifest = {
            "schema_version": "1.0.0", "manifest_state": "FINALIZED",
            "roh_simulation_manifest_id": f"{run_id}.d16", "created_at": created_at,
            "method_id": METHOD_ID, "panel_role": role,
            "scope": "REGIONAL_TECHNICAL_INJECTION", "domain": documents["D15"]["domain"],
            "configuration": {"configuration_id": config["configuration_id"],
                              "configuration_sha256": _sha(config_path)},
            "inputs": {key: _sha(path) for key, path in input_paths.items()},
            "seed": {"parent_recombination_seed": parent_seed, "derived_seed": derived_seed,
                     "derivation": "SHA256(method_id:D16:parent_seed:panel_role)[0:4]_BIG_ENDIAN"},
            "truth_policy": "INJECTED_AUTOZYGOUS_ORIGIN_BACKGROUND_ROH_UNKNOWN",
            "positive_roh_count": sum(item.is_positive for item in injections),
            "negative_roh_count": sum(not item.is_positive for item in injections),
            "ibd_coupled_roh_count": sum(item.source_ibd_scenario_id is not None for item in injections),
            "induced_ibd_edge_count": len(induced_ibd),
            "observability": density,
            "artifacts": artifacts, "integrity": integrity,
            "calibrated": False, "gamma_inputs_modified": False,
            "founder_effect_proven": False,
        }
        _validate_schema(manifest, "roh_simulation_manifest")
        manifest_path = private / "roh_simulation_manifest.json"
        manifest_path.write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")
        manifest_path.chmod(0o600)
        temporary.rename(target_root)
    except Exception:
        shutil.rmtree(temporary)
        raise
    return manifest


def materialize_roh(
    *, d15_root: Path, d14_root: Path, d11a_root: Path, config_path: Path,
    workspace_root: Path, run_id: str, created_at: str,
) -> dict[str, Any]:
    """Matérialise D16 et attribue les erreurs de lecture communes à cette étape."""
    try:
        return _materialize_roh(
            d15_root=d15_root, d14_root=d14_root, d11a_root=d11a_root,
            config_path=config_path, workspace_root=workspace_root,
            run_id=run_id, created_at=created_at,
        )
    except IbdSimulationError as error:
        raise RohSimulationError("d16_input_integrity_error") from error


def main() -> int:
    """CLI D16 : arguments explicites et compteurs agrégés uniquement."""
    parser = argparse.ArgumentParser(description="Injection ROH D16")
    for argument in ("d15-root", "d14-root", "d11a-root", "config-path", "workspace-root"):
        parser.add_argument("--" + argument, type=Path, required=True)
    parser.add_argument("--run-id", required=True)
    result = materialize_roh(**vars(parser.parse_args()), created_at=datetime.now(timezone.utc).isoformat())
    print(json.dumps(result["integrity"]))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
