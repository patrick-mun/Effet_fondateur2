"""Matérialise D15 sur un run D14 lié au reçu D14b, sans outil externe."""

from __future__ import annotations

import argparse
from dataclasses import asdict
from datetime import datetime, timezone
import csv
import gzip
import hashlib
import io
import json
from pathlib import Path
import re
import shutil
import tempfile
from typing import Any, TextIO

import jsonschema

from effet_fondateur.calibration.simulation.founder_recombination import parse_genetic_map
from effet_fondateur.calibration.simulation.ibd import (
    METHOD_ID, IbdSimulationError, build_ibd_injections, transform_ibd_vcf,
)
from effet_fondateur.calibration.simulation.population_structure import parse_observed_variant_keys


def _sha(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1048576), b""):
            digest.update(block)
    return digest.hexdigest()


def _safe_file(path: Path) -> None:
    """Refuse données d'étude, liens, traversées et fichiers absents avant lecture."""
    parts = tuple(part.casefold() for part in path.absolute().parts)
    if ".." in parts or "dock6_real" in parts or any(
        left == "data" and right == "input" for left, right in zip(parts, parts[1:])
    ):
        raise IbdSimulationError("d15_real_or_traversing_path_forbidden")
    if not path.is_file() or path.stat().st_nlink != 1:
        raise IbdSimulationError("d15_file_invalid")
    if any(parent.is_symlink() for parent in (path, *path.parents)):
        raise IbdSimulationError("d15_symlink_forbidden")


def _json(path: Path) -> dict[str, Any]:
    _safe_file(path)
    def unique(pairs: list[tuple[str, Any]]) -> dict[str, Any]:
        result: dict[str, Any] = {}
        for key, value in pairs:
            if key in result:
                raise IbdSimulationError("d15_duplicate_json_key")
            result[key] = value
        return result
    document = json.loads(path.read_text(encoding="utf-8"), object_pairs_hook=unique)
    if not isinstance(document, dict):
        raise IbdSimulationError("d15_json_object_required")
    return document


def _artifact(root: Path, document: dict[str, Any], role: str) -> Path:
    matches = [item for item in document["artifacts"] if item["role"] == role]
    if len(matches) != 1:
        raise IbdSimulationError("d15_artifact_role_invalid")
    item = matches[0]
    path = root / item["path"]
    _safe_file(path)
    if root.resolve() not in path.resolve().parents or _sha(path) != item["sha256"]:
        raise IbdSimulationError("d15_artifact_hash_or_path_invalid")
    return path


def _text(path: Path) -> TextIO:
    return gzip.open(path, "rt", encoding="utf-8") if path.suffix == ".gz" else path.open(encoding="utf-8")


def materialize_ibd(
    *, d14_root: Path, d11a_root: Path, design_root: Path, replay_root: Path,
    workspace_root: Path, run_id: str, created_at: str,
) -> dict[str, Any]:
    """Injecte une paire positive et une paire négative par longueur D14b.

    La graine D15 est dérivée par SHA-256 depuis la graine D14, le rôle de
    panel et l'identifiant de méthode. Matrices et vingt sujets sont disjoints.
    Tous les sujets modifiés sont hors famille et non porteurs D14 : leurs
    transmissions restent inchangées. Les résultats régionaux préparent 16C
    mais ne constituent pas sa qualification autosomique.
    """
    if not re.fullmatch(r"[A-Za-z0-9][A-Za-z0-9._-]*", run_id):
        raise IbdSimulationError("d15_run_id_invalid")
    workspace_root = workspace_root.absolute()
    if workspace_root.name != "calibration_workspace" or any(
        parent.is_symlink() for parent in (workspace_root, *workspace_root.parents)
    ):
        raise IbdSimulationError("d15_workspace_invalid")
    workspace_root.mkdir(parents=True, exist_ok=True)
    target_root = workspace_root / run_id
    if target_root.exists():
        raise IbdSimulationError("d15_run_exists")
    roots = {"D14": d14_root, "D11A": d11a_root, "DESIGN": design_root, "REPLAY": replay_root}
    for root in roots.values():
        if root.absolute().parent != workspace_root:
            raise IbdSimulationError("d15_input_root_invalid")
    paths = {
        "D14": d14_root / "private_truth/d14/founder_simulation_manifest.json",
        "D11A": d11a_root / "private_truth/d11/base_haplotype_manifest.json",
        "DESIGN": design_root / "private_truth/d14b/simulation_window_manifest.json",
        "REPLAY": replay_root / "private_truth/d14b_replay/simulation_window_replay.json",
    }
    docs = {key: _json(path) for key, path in paths.items()}
    d14, d11, design, replay = (docs[key] for key in ("D14", "D11A", "DESIGN", "REPLAY"))
    if any(docs[key].get("manifest_state") != "FINALIZED" for key in ("D14", "D11A", "DESIGN")) or replay.get("receipt_state") != "FINALIZED":
        raise IbdSimulationError("d15_inputs_not_finalized")
    stages = {item["stage"]: item["manifest_sha256"] for item in replay["stages"]}
    if stages.get("D14") != _sha(paths["D14"]) or stages.get("D11A") != _sha(paths["D11A"]) or replay["design"]["manifest_sha256"] != _sha(paths["DESIGN"]):
        raise IbdSimulationError("d15_replay_binding_mismatch")
    role = d14["panel"]["role"]
    if role not in {"CALIBRATION", "VALIDATION"} or role != design["panel"]["role"] or role != d11["panel"]["role"]:
        raise IbdSimulationError("d15_panel_role_mismatch")
    primary = [item for item in design["windows"] if item["purpose"] == "PRIMARY"]
    if len(primary) != 1 or any(item["status"] != "EVALUABLE" for item in primary[0]["ibd_scenarios"]):
        raise IbdSimulationError("d15_primary_window_invalid")
    primary = primary[0]
    if {item["segment_length_cm"] for item in primary["ibd_scenarios"]} != {1, 1.5, 2, 3, 5} or len(primary["ibd_scenarios"]) != 5:
        raise IbdSimulationError("d15_scientific_grid_incomplete")
    if d14["domain"]["target_position_bp"] != primary["target"]["position_bp"]:
        raise IbdSimulationError("d15_target_mismatch")
    dense_path = _artifact(d14_root, d14, "DENSE_PHASED_FOUNDER_PANEL")
    observed_path = _artifact(d14_root, d14, "BLINDED_OBSERVED_FOUNDER_PANEL")
    pedigree_path = _artifact(d14_root, d14, "PEDIGREE_TRUTH")
    carriers_path = _artifact(d14_root, d14, "FOUNDER_CARRIER_TRUTH")
    transmission_path = _artifact(d14_root, d14, "MENDELIAN_TRANSMISSION_TRUTH")
    map_path = _artifact(d11a_root, d11, "GENETIC_MAP")
    with _text(pedigree_path) as handle:
        members = list(csv.DictReader(handle, delimiter="\t"))
    with _text(carriers_path) as handle:
        carriers = {row["sample_id"] for row in csv.DictReader(handle, delimiter="\t")}
    eligible = [row["sample_id"] for row in members if row["family_unit_id"].startswith("UNRELATED_") and row["sample_id"] not in carriers]
    with _text(map_path) as handle:
        points = parse_genetic_map(handle, d14["domain"]["chromosome"])
    with _text(observed_path) as handle:
        observed_keys = parse_observed_variant_keys(handle)
    with _text(transmission_path) as handle:
        transmissions = tuple((row["child_id"], row["parent_id"], int(row["parental_haplotype_index"]), int(row["child_haplotype_index"])) for row in csv.DictReader(handle, delimiter="\t"))
    parent_seed = d14["configuration"]["random_seed"]
    seed = int.from_bytes(hashlib.sha256(f"{METHOD_ID}:D15:{parent_seed}:{role}".encode()).digest()[:4], "big")
    injections = build_ibd_injections(
        eligible_samples=eligible, map_points=points,
        chromosome=d14["domain"]["chromosome"], target_bp=primary["target"]["position_bp"],
        scenarios=primary["ibd_scenarios"], random_seed=seed,
    )
    temporary = Path(tempfile.mkdtemp(prefix=f".{run_id}.", dir=workspace_root))
    try:
        private = temporary / "private_truth/d15"
        private.mkdir(parents=True, mode=0o700)
        dense_output, observed_output = private / "ibd_panel.phased.vcf.gz", private / "observed_ibd_panel.vcf.gz"
        # mtime=0 and an empty gzip filename make the payload byte-reproducible.
        with _text(dense_path) as source, dense_output.open("wb") as raw_dense, observed_output.open("wb") as raw_observed:
            with gzip.GzipFile(filename="", fileobj=raw_dense, mode="wb", mtime=0) as zipped_dense, gzip.GzipFile(filename="", fileobj=raw_observed, mode="wb", mtime=0) as zipped_observed:
                with io.TextIOWrapper(zipped_dense, encoding="utf-8", newline="") as dense, io.TextIOWrapper(zipped_observed, encoding="utf-8", newline="") as observed:
                    integrity = transform_ibd_vcf(source, dense, observed, injections=injections, observed_keys=observed_keys, transmissions=transmissions)
        truth = private / "ibd_truth.json"
        truth_records = [{"replication_id": run_id, **asdict(item)} for item in injections]
        truth.write_text(json.dumps(truth_records, indent=2) + "\n", encoding="utf-8")
        artifacts = []
        for path, artifact_role, count in (
            (dense_output, "DENSE_PHASED_IBD_PANEL", integrity["dense_variant_count"]),
            (observed_output, "OBSERVED_IBD_PANEL", integrity["observed_variant_count"]),
            (truth, "IBD_INJECTION_TRUTH", len(injections)),
        ):
            path.chmod(0o600)
            artifacts.append({"role": artifact_role, "path": path.relative_to(temporary).as_posix(), "sha256": _sha(path), "size_bytes": path.stat().st_size, "record_count": count})
        manifest = {
            "schema_version": "1.0.0", "manifest_state": "FINALIZED",
            "ibd_simulation_manifest_id": f"{run_id}.d15", "created_at": created_at,
            "method_id": METHOD_ID, "panel_role": role,
            "scope": "REGIONAL_TECHNICAL_INJECTION", "domain": d14["domain"],
            "inputs": {key: _sha(path) for key, path in paths.items()},
            "seed": {"parent_recombination_seed": parent_seed, "derived_seed": seed,
                     "derivation": "SHA256(method_id:D15:parent_seed:panel_role)[0:4]_BIG_ENDIAN"},
            "truth_policy": "INJECTED_ORIGIN_INTERVAL_BACKGROUND_IBD_UNKNOWN",
            "positive_pair_count": sum(item.is_positive for item in injections),
            "negative_pair_count": sum(not item.is_positive for item in injections),
            "artifacts": artifacts, "integrity": integrity,
            "calibrated": False, "founder_effect_proven": False,
        }
        schema_root = Path(__file__).resolve().parents[3] / "calibration/schemas"
        for name, document in (("ibd_injection_truth", truth_records),
                               ("ibd_simulation_manifest", manifest)):
            schema = _json(schema_root / f"{name}.schema.json")
            jsonschema.Draft202012Validator(
                schema, format_checker=jsonschema.FormatChecker()
            ).validate(document)
        (private / "ibd_simulation_manifest.json").write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")
        (private / "ibd_simulation_manifest.json").chmod(0o600)
        temporary.rename(target_root)
    except Exception:
        shutil.rmtree(temporary)
        raise
    return manifest


def main() -> int:
    """CLI D15 : sorties agrégées, entrées explicitement désignées."""
    parser = argparse.ArgumentParser(description="Injection technique IBD D15")
    for argument in ("d14-root", "d11a-root", "design-root", "replay-root", "workspace-root"):
        parser.add_argument("--" + argument, type=Path, required=True)
    parser.add_argument("--run-id", required=True)
    arguments = vars(parser.parse_args())
    manifest = materialize_ibd(**arguments, created_at=datetime.now(timezone.utc).isoformat())
    print(json.dumps(manifest["integrity"]))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
