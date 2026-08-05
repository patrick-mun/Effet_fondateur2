"""Étape 07 : apparentement KING, pedigree et proposition indépendante."""

from __future__ import annotations

import argparse
import csv
import json
import math
import shutil
import subprocess
import sys
from collections import Counter, defaultdict
from pathlib import Path, PurePosixPath
from time import monotonic
from typing import Any, Sequence

from effet_fondateur.audit import atomic_write_json, read_json, sha256_file
from effet_fondateur.contracts import (
    DocumentValidationError,
    TableValidationError,
    build_file_artifact,
    load_pipeline_config,
    validate_json_document,
    validate_tsv_table,
)
from effet_fondateur.orchestrator.state import utc_now


PAIR_COLUMNS = (
    "PAIR_ID",
    "SAMPLE_ID_1",
    "SAMPLE_ID_2",
    "FID_1",
    "IID_1",
    "FID_2",
    "IID_2",
    "N_SNP",
    "HETHET",
    "IBS0",
    "KINSHIP",
    "INFERRED_DEGREE",
    "INFERRED_RELATION",
    "DECLARED_RELATION",
    "PAIR_ORIGIN",
    "CONCORDANCE_STATUS",
)
DEGREE_SUMMARY_COLUMNS = (
    "INFERRED_DEGREE",
    "PAIR_COUNT",
    "DECLARED_PAIR_COUNT",
    "CRYPTIC_PAIR_COUNT",
    "CONCORDANT_PAIR_COUNT",
    "DISCORDANT_PAIR_COUNT",
)
PEDIGREE_COLUMNS = (
    "RELATION_ID",
    "CHILD_SAMPLE_ID",
    "PARENT_ROLE",
    "DECLARED_PARENT_IID",
    "PARENT_SAMPLE_ID",
    "INFERRED_DEGREE",
    "INFERRED_RELATION",
    "KINSHIP",
    "IBS0",
    "CONCORDANCE_STATUS",
)
PROPOSAL_COLUMNS = (
    "SAMPLE_ID",
    "FID",
    "IID",
    "COMPONENT_ID",
    "PROPOSED_INCLUDED",
    "PROPOSAL_STATUS",
    "REPRESENTATIVE_SAMPLE_ID",
    "RELATION_BASIS",
    "F_MISS",
    "ABS_HETEROZYGOSITY_Z",
    "DECISION_SOURCE",
    "MANUAL_VALIDATION_REQUIRED",
)
DEGREE_ORDER = (
    "DUPLICATE_OR_MZ",
    "FIRST",
    "SECOND",
    "THIRD",
    "FOURTH",
    "UNRELATED",
)
RELATED_DEGREE_RANK = {
    "DUPLICATE_OR_MZ": 0,
    "FIRST": 1,
    "SECOND": 2,
    "THIRD": 3,
    "FOURTH": 4,
    "UNRELATED": 99,
}


class KinshipInputError(ValueError):
    """Signale une entrée ou un paramètre KING invalide."""


class ExternalToolError(RuntimeError):
    """Signale un échec ou une sortie invalide de KING."""


def _build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(prog="infer-kinship")
    parser.add_argument("--stage-inputs", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    return parser


def _resolve_input_path(artifact: dict[str, Any], run_dir: Path) -> Path:
    artifact_path = Path(artifact["path"])
    if artifact_path.is_absolute():
        return artifact_path
    if PurePosixPath(artifact["path"]).parts[:1] == ("stages",):
        return run_dir / artifact_path
    return Path.cwd() / artifact_path


def _artifact_by_id(stage_inputs: dict[str, Any], artifact_id: str) -> dict[str, Any]:
    matches = [
        artifact
        for artifact in stage_inputs["artifacts"]
        if artifact["artifact_id"] == artifact_id
    ]
    if len(matches) != 1:
        raise KinshipInputError(f"{artifact_id}_missing_or_ambiguous")
    return matches[0]


def _validated_artifact_path(artifact: dict[str, Any], run_dir: Path) -> Path:
    physical_path = _resolve_input_path(artifact, run_dir)
    if not physical_path.is_file() or sha256_file(physical_path) != artifact["sha256"]:
        raise KinshipInputError("declared_input_modified")
    return physical_path


def _read_fam(path: Path) -> list[list[str]]:
    rows = [line.split() for line in path.read_text(encoding="utf-8").splitlines() if line]
    if not rows or any(len(row) != 6 for row in rows):
        raise KinshipInputError("malformed_kinship_panel_fam")
    if len({(row[0], row[1]) for row in rows}) != len(rows):
        raise KinshipInputError("duplicate_kinship_panel_sample")
    return rows


def _read_bim(path: Path) -> list[list[str]]:
    rows = [line.split() for line in path.read_text(encoding="utf-8").splitlines() if line]
    if not rows or any(len(row) != 6 for row in rows):
        raise KinshipInputError("malformed_kinship_panel_bim")
    return rows


def _positive_integer(parameters: dict[str, Any], name: str, default: int) -> int:
    value = parameters.get(name, default)
    if isinstance(value, bool) or not isinstance(value, int) or value <= 0:
        raise KinshipInputError(f"invalid_parameter:{name}")
    return value


def _bounded_float(
    parameters: dict[str, Any],
    name: str,
    default: float,
    *,
    minimum: float = 0.0,
    maximum: float = 0.5,
) -> float:
    value = parameters.get(name, default)
    if isinstance(value, bool) or not isinstance(value, (int, float)):
        raise KinshipInputError(f"invalid_parameter:{name}")
    result = float(value)
    if not math.isfinite(result) or not minimum <= result <= maximum:
        raise KinshipInputError(f"invalid_parameter:{name}")
    return result


def _parameters(parameters: dict[str, Any]) -> dict[str, Any]:
    king_degree = _positive_integer(parameters, "king_degree", 4)
    independence_degree = _positive_integer(
        parameters, "relatedness_max_degree_for_independence", 4
    )
    if king_degree > 4 or independence_degree > king_degree:
        raise KinshipInputError("invalid_king_degree_configuration")
    thresholds = {
        "duplicate_or_mz_min": _bounded_float(
            parameters, "duplicate_or_mz_min", 0.354
        ),
        "first_degree_min": _bounded_float(parameters, "first_degree_min", 0.177),
        "second_degree_min": _bounded_float(
            parameters, "second_degree_min", 0.0884
        ),
        "third_degree_min": _bounded_float(parameters, "third_degree_min", 0.0442),
        "fourth_degree_min": _bounded_float(
            parameters, "fourth_degree_min", 0.0221
        ),
    }
    ordered = [
        thresholds["duplicate_or_mz_min"],
        thresholds["first_degree_min"],
        thresholds["second_degree_min"],
        thresholds["third_degree_min"],
        thresholds["fourth_degree_min"],
    ]
    if any(left <= right for left, right in zip(ordered, ordered[1:])):
        raise KinshipInputError("kinship_thresholds_not_strictly_descending")
    return {
        "king_degree": king_degree,
        **thresholds,
        "parent_child_ibs0_max": _bounded_float(
            parameters, "parent_child_ibs0_max", 0.01, maximum=1.0
        ),
        "relatedness_max_degree_for_independence": independence_degree,
        "threads": _positive_integer(parameters, "threads", 4),
        "king_timeout_seconds": _positive_integer(
            parameters, "king_timeout_seconds", 300
        ),
    }


def _resolve_executable(command: str | None) -> str:
    if command is None:
        raise ExternalToolError("king_not_configured")
    executable = shutil.which(command)
    if executable is None:
        raise ExternalToolError("king_not_available")
    return executable


def _king_version(executable: str) -> str | None:
    try:
        completed = subprocess.run(
            [executable, "--version"],
            capture_output=True,
            check=False,
            text=True,
            timeout=10,
        )
    except (OSError, subprocess.TimeoutExpired):
        return None
    output = completed.stdout.strip() or completed.stderr.strip()
    return output.splitlines()[0][:500] if output else None


def _run_king(
    executable: str,
    bed_path: Path,
    output_prefix: Path,
    degree: int,
    threads: int,
    timeout_seconds: int,
) -> None:
    command = [
        executable,
        "-b",
        str(bed_path),
        "--kinship",
        "--degree",
        str(degree),
        "--cpus",
        str(threads),
        "--prefix",
        str(output_prefix),
    ]
    try:
        completed = subprocess.run(
            command,
            capture_output=True,
            check=False,
            text=True,
            timeout=timeout_seconds,
        )
    except (OSError, subprocess.TimeoutExpired) as error:
        raise ExternalToolError("king_execution_failed") from error
    if completed.returncode != 0:
        raise ExternalToolError("king_command_failed")


def _read_king_table(path: Path, *, required: bool) -> list[dict[str, str]]:
    if not path.exists():
        if required:
            raise ExternalToolError(f"missing_king_report:{path.name}")
        return []
    lines = [line.split() for line in path.read_text(encoding="utf-8").splitlines() if line]
    if not lines:
        raise ExternalToolError(f"empty_king_report:{path.name}")
    header = lines[0]
    if any(len(row) != len(header) for row in lines[1:]):
        raise ExternalToolError(f"malformed_king_report:{path.name}")
    return [dict(zip(header, row, strict=True)) for row in lines[1:]]


def _validate_panel_dataset(
    descriptor: dict[str, Any],
    paths: dict[str, Path],
    artifacts: dict[str, dict[str, Any]],
) -> tuple[list[list[str]], list[list[str]]]:
    validate_json_document(descriptor, "plink_dataset.schema.json")
    if descriptor["dataset_id"] != "kinship_panel":
        raise KinshipInputError("unexpected_kinship_panel_dataset_id")
    if descriptor["scope"] != "autosomal_genomewide":
        raise KinshipInputError("kinship_requires_autosomal_panel")
    if paths["kinship_panel_bed"].read_bytes()[:3] != b"\x6c\x1b\x01":
        raise KinshipInputError("invalid_kinship_panel_bed_header")
    for suffix in ("bed", "bim", "fam"):
        artifact_id = f"kinship_panel_{suffix}"
        physical_path = paths[artifact_id]
        declared = descriptor["files"][suffix]
        if (
            declared["name"] != physical_path.name
            or declared["sha256"] != sha256_file(physical_path)
            or artifacts[artifact_id]["sample_set_id"] != descriptor["sample_set_id"]
            or artifacts[artifact_id]["variant_set_id"] != descriptor["variant_set_id"]
        ):
            raise KinshipInputError("kinship_panel_descriptor_mismatch")
    fam_rows = _read_fam(paths["kinship_panel_fam"])
    bim_rows = _read_bim(paths["kinship_panel_bim"])
    if descriptor["sample_count"] != len(fam_rows) or descriptor["variant_count"] != len(
        bim_rows
    ):
        raise KinshipInputError("kinship_panel_count_mismatch")
    return fam_rows, bim_rows


def _sample_context(
    samples_path: Path,
    qc_path: Path,
    fam_rows: list[list[str]],
) -> tuple[
    list[str],
    dict[str, dict[str, Any]],
    dict[tuple[str, str], str],
    dict[str, dict[str, Any]],
]:
    samples_table = validate_tsv_table(samples_path, "samples_master.schema.json")
    qc_table = validate_tsv_table(qc_path, "qc_individual_metrics.schema.json")
    samples_by_plink = {
        (row["FID"], row["IID"]): row for row in samples_table.rows
    }
    qc_by_sample = {row["SAMPLE_ID"]: row for row in qc_table.rows}
    ordered_sample_ids: list[str] = []
    panel_samples: dict[str, dict[str, Any]] = {}
    plink_to_sample: dict[tuple[str, str], str] = {}
    for fam_row in fam_rows:
        plink_id = (fam_row[0], fam_row[1])
        sample = samples_by_plink.get(plink_id)
        if sample is None:
            raise KinshipInputError("kinship_panel_sample_absent_from_registry")
        sample_id = sample["SAMPLE_ID"]
        qc = qc_by_sample.get(sample_id)
        if qc is None or qc["QC_STATUS"] == "EXCLUDED":
            raise KinshipInputError("kinship_panel_sample_absent_from_preliminary_qc")
        ordered_sample_ids.append(sample_id)
        panel_samples[sample_id] = sample
        plink_to_sample[plink_id] = sample_id
    return ordered_sample_ids, panel_samples, plink_to_sample, qc_by_sample


def _float_field(row: dict[str, str], field: str, *, minimum: float | None = None) -> float:
    try:
        value = float(row[field])
    except (KeyError, ValueError) as error:
        raise ExternalToolError(f"invalid_king_field:{field}") from error
    if not math.isfinite(value) or (minimum is not None and value < minimum):
        raise ExternalToolError(f"invalid_king_field:{field}")
    return value


def _integer_field(row: dict[str, str], field: str) -> int:
    try:
        value = int(row[field])
    except (KeyError, ValueError) as error:
        raise ExternalToolError(f"invalid_king_field:{field}") from error
    if value <= 0:
        raise ExternalToolError(f"invalid_king_field:{field}")
    return value


def _king_identifier(row: dict[str, str], *field_names: str) -> str:
    for field_name in field_names:
        value = row.get(field_name)
        if value:
            return value
    raise ExternalToolError(f"missing_king_identifier:{'/'.join(field_names)}")


def _normalized_pair(sample_id_1: str, sample_id_2: str) -> tuple[str, str]:
    if sample_id_1 == sample_id_2:
        raise ExternalToolError("king_self_pair")
    return tuple(sorted((sample_id_1, sample_id_2)))


def _parse_king_pairs(
    within_rows: list[dict[str, str]],
    between_rows: list[dict[str, str]],
    plink_to_sample: dict[tuple[str, str], str],
) -> list[dict[str, Any]]:
    pairs: dict[tuple[str, str], dict[str, Any]] = {}
    for origin, rows in (
        ("WITHIN_FAMILY", within_rows),
        ("BETWEEN_FAMILY", between_rows),
    ):
        for row in rows:
            if origin == "WITHIN_FAMILY":
                fid_1 = fid_2 = _king_identifier(row, "FID")
            else:
                fid_1 = _king_identifier(row, "FID1")
                fid_2 = _king_identifier(row, "FID2")
            iid_1 = _king_identifier(row, "ID1", "IID1")
            iid_2 = _king_identifier(row, "ID2", "IID2")
            try:
                sample_id_1 = plink_to_sample[(fid_1, iid_1)]
                sample_id_2 = plink_to_sample[(fid_2, iid_2)]
            except KeyError as error:
                raise ExternalToolError("king_pair_sample_absent_from_panel") from error
            pair_key = _normalized_pair(sample_id_1, sample_id_2)
            if pair_key in pairs:
                raise ExternalToolError("duplicate_king_pair")
            pairs[pair_key] = {
                "sample_id_1": pair_key[0],
                "sample_id_2": pair_key[1],
                "n_snp": _integer_field(row, "N_SNP"),
                "hethet": _float_field(row, "HetHet", minimum=0.0),
                "ibs0": _float_field(row, "IBS0", minimum=0.0),
                "kinship": _float_field(row, "Kinship"),
                "origin": origin,
            }
    return [pairs[pair_key] for pair_key in sorted(pairs)]


def _classify_relationship(
    kinship: float,
    ibs0: float,
    parameters: dict[str, Any],
) -> tuple[str, str]:
    if kinship >= parameters["duplicate_or_mz_min"]:
        return "DUPLICATE_OR_MZ", "DUPLICATE_OR_MZ"
    if kinship >= parameters["first_degree_min"]:
        relation = (
            "PARENT_CHILD_COMPATIBLE"
            if ibs0 <= parameters["parent_child_ibs0_max"]
            else "FIRST_DEGREE_NON_PARENT_CHILD"
        )
        return "FIRST", relation
    if kinship >= parameters["second_degree_min"]:
        return "SECOND", "SECOND_DEGREE"
    if kinship >= parameters["third_degree_min"]:
        return "THIRD", "THIRD_DEGREE"
    if kinship >= parameters["fourth_degree_min"]:
        return "FOURTH", "FOURTH_DEGREE"
    return "UNRELATED", "UNRELATED"


def _declared_parent_pairs(
    ordered_sample_ids: list[str],
    samples: dict[str, dict[str, Any]],
    plink_to_sample: dict[tuple[str, str], str],
) -> tuple[dict[tuple[str, str], list[dict[str, Any]]], list[dict[str, Any]]]:
    pairs: dict[tuple[str, str], list[dict[str, Any]]] = defaultdict(list)
    relations: list[dict[str, Any]] = []
    for child_sample_id in ordered_sample_ids:
        child = samples[child_sample_id]
        for parent_role, field in (("FATHER", "PID"), ("MOTHER", "MID")):
            parent_iid = child[field]
            if parent_iid == "0":
                continue
            parent_sample_id = plink_to_sample.get((child["FID"], parent_iid))
            relation = {
                "child_sample_id": child_sample_id,
                "parent_role": parent_role,
                "declared_parent_iid": parent_iid,
                "parent_sample_id": parent_sample_id,
            }
            relations.append(relation)
            if parent_sample_id is not None:
                pairs[_normalized_pair(child_sample_id, parent_sample_id)].append(relation)
    return dict(pairs), relations


def _pair_concordance(
    declared_relation: str,
    inferred_degree: str,
    inferred_relation: str,
) -> str:
    if declared_relation == "PARENT_CHILD":
        if inferred_degree != "FIRST":
            return "DISCORDANT_DEGREE"
        if inferred_relation != "PARENT_CHILD_COMPATIBLE":
            return "DISCORDANT_RELATION"
        return "CONCORDANT"
    if declared_relation == "SAME_FAMILY_UNSPECIFIED":
        return (
            "FAMILY_RELATION_UNRESOLVED"
            if inferred_degree != "UNRELATED"
            else "NO_DECLARED_RELATION"
        )
    return (
        "CRYPTIC_RELATEDNESS"
        if inferred_degree != "UNRELATED"
        else "NO_DECLARED_RELATION"
    )


def _decimal(value: float | None) -> str:
    if value is None:
        return ""
    normalized = 0.0 if abs(value) < 5e-13 else value
    return f"{normalized:.8f}".rstrip("0").rstrip(".")


def _pair_rows(
    king_pairs: list[dict[str, Any]],
    samples: dict[str, dict[str, Any]],
    declared_parent_pairs: dict[tuple[str, str], list[dict[str, Any]]],
    parameters: dict[str, Any],
) -> list[dict[str, str]]:
    rows: list[dict[str, str]] = []
    for index, pair in enumerate(king_pairs, start=1):
        sample_id_1 = pair["sample_id_1"]
        sample_id_2 = pair["sample_id_2"]
        sample_1 = samples[sample_id_1]
        sample_2 = samples[sample_id_2]
        pair_key = (sample_id_1, sample_id_2)
        if pair_key in declared_parent_pairs:
            declared_relation = "PARENT_CHILD"
        elif sample_1["FID"] == sample_2["FID"]:
            declared_relation = "SAME_FAMILY_UNSPECIFIED"
        else:
            declared_relation = "NONE"
        inferred_degree, inferred_relation = _classify_relationship(
            pair["kinship"], pair["ibs0"], parameters
        )
        rows.append(
            {
                "PAIR_ID": f"kinship_pair_{index:06d}",
                "SAMPLE_ID_1": sample_id_1,
                "SAMPLE_ID_2": sample_id_2,
                "FID_1": sample_1["FID"],
                "IID_1": sample_1["IID"],
                "FID_2": sample_2["FID"],
                "IID_2": sample_2["IID"],
                "N_SNP": str(pair["n_snp"]),
                "HETHET": _decimal(pair["hethet"]),
                "IBS0": _decimal(pair["ibs0"]),
                "KINSHIP": _decimal(pair["kinship"]),
                "INFERRED_DEGREE": inferred_degree,
                "INFERRED_RELATION": inferred_relation,
                "DECLARED_RELATION": declared_relation,
                "PAIR_ORIGIN": pair["origin"],
                "CONCORDANCE_STATUS": _pair_concordance(
                    declared_relation, inferred_degree, inferred_relation
                ),
            }
        )
    return rows


def _pedigree_rows(
    declared_relations: list[dict[str, Any]],
    pair_rows: list[dict[str, str]],
) -> list[dict[str, str]]:
    pairs_by_key = {
        _normalized_pair(row["SAMPLE_ID_1"], row["SAMPLE_ID_2"]): row
        for row in pair_rows
    }
    rows: list[dict[str, str]] = []
    for index, relation in enumerate(declared_relations, start=1):
        parent_sample_id = relation["parent_sample_id"]
        pair = (
            pairs_by_key.get(
                _normalized_pair(relation["child_sample_id"], parent_sample_id)
            )
            if parent_sample_id is not None
            else None
        )
        if parent_sample_id is None:
            status = "NOT_EVALUATED_PARENT_ABSENT"
        elif pair is None:
            status = "NOT_EVALUATED_KING_PAIR_MISSING"
        else:
            status = pair["CONCORDANCE_STATUS"]
        rows.append(
            {
                "RELATION_ID": f"pedigree_relation_{index:06d}",
                "CHILD_SAMPLE_ID": relation["child_sample_id"],
                "PARENT_ROLE": relation["parent_role"],
                "DECLARED_PARENT_IID": relation["declared_parent_iid"],
                "PARENT_SAMPLE_ID": parent_sample_id or "",
                "INFERRED_DEGREE": pair["INFERRED_DEGREE"] if pair else "",
                "INFERRED_RELATION": pair["INFERRED_RELATION"] if pair else "",
                "KINSHIP": pair["KINSHIP"] if pair else "",
                "IBS0": pair["IBS0"] if pair else "",
                "CONCORDANCE_STATUS": status,
            }
        )
    return rows


def _degree_summary_rows(pair_rows: list[dict[str, str]]) -> list[dict[str, str]]:
    rows: list[dict[str, str]] = []
    for degree in DEGREE_ORDER:
        matching = [row for row in pair_rows if row["INFERRED_DEGREE"] == degree]
        rows.append(
            {
                "INFERRED_DEGREE": degree,
                "PAIR_COUNT": str(len(matching)),
                "DECLARED_PAIR_COUNT": str(
                    sum(row["DECLARED_RELATION"] != "NONE" for row in matching)
                ),
                "CRYPTIC_PAIR_COUNT": str(
                    sum(row["CONCORDANCE_STATUS"] == "CRYPTIC_RELATEDNESS" for row in matching)
                ),
                "CONCORDANT_PAIR_COUNT": str(
                    sum(row["CONCORDANCE_STATUS"] == "CONCORDANT" for row in matching)
                ),
                "DISCORDANT_PAIR_COUNT": str(
                    sum(row["CONCORDANCE_STATUS"].startswith("DISCORDANT") for row in matching)
                ),
            }
        )
    return rows


def _related_edges(
    ordered_sample_ids: list[str],
    pair_rows: list[dict[str, str]],
    declared_parent_pairs: dict[tuple[str, str], list[dict[str, Any]]],
    max_degree: int,
) -> tuple[dict[str, set[str]], dict[tuple[str, str], set[str]]]:
    adjacency = {sample_id: set() for sample_id in ordered_sample_ids}
    edge_sources: dict[tuple[str, str], set[str]] = defaultdict(set)
    for row in pair_rows:
        degree_rank = RELATED_DEGREE_RANK[row["INFERRED_DEGREE"]]
        if degree_rank <= max_degree:
            sample_id_1 = row["SAMPLE_ID_1"]
            sample_id_2 = row["SAMPLE_ID_2"]
            adjacency[sample_id_1].add(sample_id_2)
            adjacency[sample_id_2].add(sample_id_1)
            edge_sources[_normalized_pair(sample_id_1, sample_id_2)].add(
                "KING_RELATEDNESS"
            )
    for sample_id_1, sample_id_2 in declared_parent_pairs:
        adjacency[sample_id_1].add(sample_id_2)
        adjacency[sample_id_2].add(sample_id_1)
        edge_sources[_normalized_pair(sample_id_1, sample_id_2)].add(
            "DECLARED_PARENT_CHILD"
        )
    return adjacency, dict(edge_sources)


def _components(
    ordered_sample_ids: list[str], adjacency: dict[str, set[str]]
) -> list[list[str]]:
    order = {sample_id: index for index, sample_id in enumerate(ordered_sample_ids)}
    unseen = set(ordered_sample_ids)
    components: list[list[str]] = []
    while unseen:
        start = min(unseen, key=order.__getitem__)
        stack = [start]
        component: set[str] = set()
        while stack:
            sample_id = stack.pop()
            if sample_id in component:
                continue
            component.add(sample_id)
            stack.extend(adjacency[sample_id] - component)
        unseen -= component
        components.append(sorted(component, key=order.__getitem__))
    return components


def _quality_key(sample_id: str, qc_by_sample: dict[str, dict[str, Any]]) -> tuple[Any, ...]:
    qc = qc_by_sample[sample_id]
    status_rank = 0 if qc["QC_STATUS"] == "PASS" else 1
    heterozygosity_z = qc["HETEROZYGOSITY_Z"]
    absolute_z = abs(float(heterozygosity_z)) if heterozygosity_z is not None else math.inf
    return status_rank, float(qc["F_MISS"]), absolute_z, sample_id


def _independent_set_rows(
    ordered_sample_ids: list[str],
    samples: dict[str, dict[str, Any]],
    qc_by_sample: dict[str, dict[str, Any]],
    adjacency: dict[str, set[str]],
    edge_sources: dict[tuple[str, str], set[str]],
) -> list[dict[str, str]]:
    selections: dict[str, tuple[bool, str | None, str]] = {}
    for component_index, component in enumerate(
        _components(ordered_sample_ids, adjacency), start=1
    ):
        component_id = f"kinship_component_{component_index:06d}"
        remaining = set(component)
        while remaining:
            representative = min(
                remaining, key=lambda sample_id: _quality_key(sample_id, qc_by_sample)
            )
            selections[representative] = (True, None, component_id)
            excluded_neighbors = adjacency[representative] & remaining
            for neighbor in excluded_neighbors:
                selections[neighbor] = (False, representative, component_id)
            remaining -= excluded_neighbors | {representative}
    rows: list[dict[str, str]] = []
    for sample_id in ordered_sample_ids:
        included, representative, component_id = selections[sample_id]
        qc = qc_by_sample[sample_id]
        heterozygosity_z = qc["HETEROZYGOSITY_Z"]
        relation_basis = ""
        if representative is not None:
            sources = edge_sources[_normalized_pair(sample_id, representative)]
            relation_basis = (
                "KING_AND_DECLARED"
                if len(sources) == 2
                else next(iter(sources))
            )
        rows.append(
            {
                "SAMPLE_ID": sample_id,
                "FID": samples[sample_id]["FID"],
                "IID": samples[sample_id]["IID"],
                "COMPONENT_ID": component_id,
                "PROPOSED_INCLUDED": str(included).lower(),
                "PROPOSAL_STATUS": "RETAIN" if included else "EXCLUDE_RELATED",
                "REPRESENTATIVE_SAMPLE_ID": representative or "",
                "RELATION_BASIS": relation_basis,
                "F_MISS": _decimal(float(qc["F_MISS"])),
                "ABS_HETEROZYGOSITY_Z": (
                    _decimal(abs(float(heterozygosity_z)))
                    if heterozygosity_z is not None
                    else ""
                ),
                "DECISION_SOURCE": "stage_07_quality_aware_greedy_maximal",
                "MANUAL_VALIDATION_REQUIRED": str(not included).lower(),
            }
        )
    return rows


def _write_tsv(
    path: Path, columns: tuple[str, ...], rows: list[dict[str, str]]
) -> None:
    with path.open("w", encoding="utf-8", newline="") as output_file:
        writer = csv.DictWriter(
            output_file, fieldnames=columns, delimiter="\t", lineterminator="\n"
        )
        writer.writeheader()
        writer.writerows({column: row[column] for column in columns} for row in rows)


def _output_artifacts(
    stage_inputs: dict[str, Any],
    descriptor: dict[str, Any],
    contracts: list[tuple[str, str, Path, str | None, str]],
) -> list[dict[str, Any]]:
    producer_stage = f"{stage_inputs['stage_id']}_{stage_inputs['stage_name']}"
    published_dir = PurePosixPath(stage_inputs["published_output_dir"])
    artifacts = []
    for artifact_id, artifact_type, path, schema_name, sensitivity in contracts:
        artifacts.append(
            build_file_artifact(
                physical_path=path,
                published_path=str(published_dir / path.name),
                artifact_id=artifact_id,
                artifact_type=artifact_type,
                media_type=(
                    "application/json"
                    if path.suffix == ".json"
                    else "text/tab-separated-values"
                ),
                producer_stage=producer_stage,
                producer_signature=stage_inputs["signature"],
                schema_name=schema_name,
                schema_version="1.0.0" if schema_name else None,
                assembly=descriptor["assembly"],
                sample_set_id=descriptor["sample_set_id"],
                variant_set_id=descriptor["variant_set_id"],
                sensitivity=sensitivity,
            )
        )
    return artifacts


def execute(stage_inputs_path: Path, output_dir: Path) -> int:
    """Infère l'apparentement et publie une proposition indépendante révisable."""
    started_at = utc_now()
    started_clock = monotonic()
    stage_inputs = read_json(stage_inputs_path)
    validate_json_document(stage_inputs, "stage_inputs.schema.json")
    run_dir = output_dir.parent.parent
    config = load_pipeline_config(run_dir / "config.resolved.yaml")
    parameters = _parameters(stage_inputs["parameters"])
    artifact_ids = (
        "samples_master",
        "qc_individual_metrics",
        "kinship_panel_bed",
        "kinship_panel_bim",
        "kinship_panel_fam",
        "kinship_panel_dataset",
    )
    input_artifacts = {
        artifact_id: _artifact_by_id(stage_inputs, artifact_id)
        for artifact_id in artifact_ids
    }
    paths = {
        artifact_id: _validated_artifact_path(artifact, run_dir)
        for artifact_id, artifact in input_artifacts.items()
    }
    descriptor = read_json(paths["kinship_panel_dataset"])
    fam_rows, _ = _validate_panel_dataset(descriptor, paths, input_artifacts)
    ordered_sample_ids, samples, plink_to_sample, qc_by_sample = _sample_context(
        paths["samples_master"], paths["qc_individual_metrics"], fam_rows
    )
    declared_parent_pairs, declared_relations = _declared_parent_pairs(
        ordered_sample_ids, samples, plink_to_sample
    )

    output_dir.mkdir(parents=True, exist_ok=True)
    king_executable = _resolve_executable(config["tools"]["king"])
    king_version = _king_version(king_executable)
    king_prefix = output_dir / "king_kinship"
    _run_king(
        king_executable,
        paths["kinship_panel_bed"],
        king_prefix,
        parameters["king_degree"],
        parameters["threads"],
        parameters["king_timeout_seconds"],
    )
    within_rows = _read_king_table(king_prefix.with_suffix(".kin"), required=False)
    between_rows = _read_king_table(king_prefix.with_suffix(".kin0"), required=True)
    king_pairs = _parse_king_pairs(within_rows, between_rows, plink_to_sample)
    pair_rows = _pair_rows(king_pairs, samples, declared_parent_pairs, parameters)
    pedigree_rows = _pedigree_rows(declared_relations, pair_rows)
    degree_summary_rows = _degree_summary_rows(pair_rows)
    adjacency, edge_sources = _related_edges(
        ordered_sample_ids,
        pair_rows,
        declared_parent_pairs,
        parameters["relatedness_max_degree_for_independence"],
    )
    proposal_rows = _independent_set_rows(
        ordered_sample_ids, samples, qc_by_sample, adjacency, edge_sources
    )

    pair_path = output_dir / "kinship_pairs.tsv"
    degree_summary_path = output_dir / "kinship_degree_summary.tsv"
    pedigree_path = output_dir / "pedigree_concordance.tsv"
    proposal_path = output_dir / "independent_set_proposal.tsv"
    _write_tsv(pair_path, PAIR_COLUMNS, pair_rows)
    _write_tsv(degree_summary_path, DEGREE_SUMMARY_COLUMNS, degree_summary_rows)
    _write_tsv(pedigree_path, PEDIGREE_COLUMNS, pedigree_rows)
    _write_tsv(proposal_path, PROPOSAL_COLUMNS, proposal_rows)
    validate_tsv_table(pair_path, "kinship_pairs.schema.json")
    validate_tsv_table(degree_summary_path, "kinship_degree_summary.schema.json")
    validate_tsv_table(pedigree_path, "pedigree_concordance.schema.json")
    validate_tsv_table(proposal_path, "independent_set_proposal.schema.json")

    proposed_exclusion_count = sum(
        row["PROPOSAL_STATUS"] == "EXCLUDE_RELATED" for row in proposal_rows
    )
    discordant_pedigree_count = sum(
        row["CONCORDANCE_STATUS"] != "CONCORDANT" for row in pedigree_rows
    )
    manual_validation_required = (
        proposed_exclusion_count > 0 or discordant_pedigree_count > 0
    )
    concordance_counts = Counter(row["CONCORDANCE_STATUS"] for row in pair_rows)
    pedigree_status_counts = Counter(
        row["CONCORDANCE_STATUS"] for row in pedigree_rows
    )
    report_path = output_dir / "kinship_report.json"
    report = {
        "schema_version": "1.0.0",
        "parameters": parameters,
        "sample_count": len(ordered_sample_ids),
        "reported_pair_count": len(pair_rows),
        "related_pair_count": sum(
            row["INFERRED_DEGREE"] != "UNRELATED" for row in pair_rows
        ),
        "cryptic_pair_count": concordance_counts["CRYPTIC_RELATEDNESS"],
        "declared_parent_relation_count": len(pedigree_rows),
        "pedigree_status_counts": dict(sorted(pedigree_status_counts.items())),
        "proposed_inclusion_count": len(proposal_rows) - proposed_exclusion_count,
        "proposed_exclusion_count": proposed_exclusion_count,
        "manual_validation_required": manual_validation_required,
        "selection_method": "quality_aware_greedy_maximal_independent_set_v1",
    }
    atomic_write_json(report_path, report)

    for suffix in ("kin", "kin0", "log"):
        raw_path = king_prefix.with_suffix(f".{suffix}")
        if raw_path.exists():
            raw_path.unlink()

    contracts = [
        (
            "kinship_pairs",
            "kinship_pairs",
            pair_path,
            "kinship_pairs.schema.json",
            "sensitive_genetic",
        ),
        (
            "kinship_degree_summary",
            "kinship_degree_summary",
            degree_summary_path,
            "kinship_degree_summary.schema.json",
            "internal",
        ),
        (
            "pedigree_concordance",
            "pedigree_concordance",
            pedigree_path,
            "pedigree_concordance.schema.json",
            "sensitive_genetic",
        ),
        (
            "independent_set_proposal",
            "independent_set_proposal",
            proposal_path,
            "independent_set_proposal.schema.json",
            "sensitive_genetic",
        ),
        (
            "kinship_report",
            "kinship_report",
            report_path,
            None,
            "internal",
        ),
    ]
    output_artifacts = _output_artifacts(stage_inputs, descriptor, contracts)
    stage_outputs = {
        "schema_version": "1.0.0",
        "run_id": stage_inputs["run_id"],
        "stage_id": stage_inputs["stage_id"],
        "stage_name": stage_inputs["stage_name"],
        "signature": stage_inputs["signature"],
        "artifacts": output_artifacts,
    }
    validate_json_document(stage_outputs, "stage_outputs.schema.json")
    atomic_write_json(output_dir / "stage_outputs.json", stage_outputs)
    audit = {
        "schema_version": "1.0.0",
        "run_id": stage_inputs["run_id"],
        "stage_id": stage_inputs["stage_id"],
        "stage_name": stage_inputs["stage_name"],
        "method_id": "king_robust_kinship_pedigree_concordance_v1",
        "signature": stage_inputs["signature"],
        "started_at": started_at,
        "completed_at": utc_now(),
        "duration_seconds": monotonic() - started_clock,
        "inputs": list(input_artifacts.values()),
        "outputs": output_artifacts,
        "parameters": parameters,
        "tools": [
            {
                "tool": "king",
                "configured": config["tools"]["king"],
                "version": king_version,
            }
        ],
        "counts": {
            "samples": len(ordered_sample_ids),
            "reported_pairs": len(pair_rows),
            "related_pairs": report["related_pair_count"],
            "declared_parent_relations": len(pedigree_rows),
            "proposed_exclusions": proposed_exclusion_count,
        },
        "metrics": {
            "cryptic_pair_count": report["cryptic_pair_count"],
            "pedigree_discordance_count": discordant_pedigree_count,
        },
        "exclusions": [],
        "warnings": (
            [
                {
                    "code": "manual_kinship_review_required",
                    "proposed_exclusion_count": proposed_exclusion_count,
                    "pedigree_issue_count": discordant_pedigree_count,
                }
            ]
            if manual_validation_required
            else []
        ),
        "checks": [
            {"check": "kinship_panel_integrity", "status": "PASS"},
            {"check": "king_outputs", "status": "PASS"},
            {"check": "degree_classification", "status": "PASS"},
            {
                "check": "pedigree_concordance",
                "status": "PASS" if discordant_pedigree_count == 0 else "REVIEW_REQUIRED",
            },
            {
                "check": "independent_set_proposal",
                "status": "REVIEW_REQUIRED" if proposed_exclusion_count else "PASS",
            },
        ],
        "known_limits": [
            "KING ne publie entre familles que les paires jusqu'au degré demandé ; "
            "les autres paires ne sont pas matérialisées comme non apparentées.",
            "La proposition indépendante est gloutonne, déterministe et maximale, "
            "mais ne garantit pas la cardinalité maximum.",
            "Aucune exclusion proposée n'est appliquée automatiquement.",
        ],
        "expected_visualizations": [
            "pseudonymized_kinship_network",
            "kinship_matrix",
            "kinship_distribution",
            "pedigree_concordance",
        ],
        "manual_validation_required": manual_validation_required,
    }
    validate_json_document(audit, "stage_audit.schema.json")
    atomic_write_json(output_dir / "audit.json", audit)
    (output_dir / "checksums.sha256").write_text(
        "".join(
            f"{artifact['sha256']}  {Path(str(artifact['path'])).name}\n"
            for artifact in output_artifacts
        ),
        encoding="utf-8",
    )
    return 0


def main(arguments: Sequence[str] | None = None) -> int:
    """Point d'entrée autonome respectant les codes de retour V2."""
    parser = _build_parser()
    parsed_arguments = parser.parse_args(arguments)
    try:
        return execute(parsed_arguments.stage_inputs, parsed_arguments.output_dir)
    except ExternalToolError as error:
        sys.stderr.write(f"{error}\n")
        return 3
    except (
        KinshipInputError,
        DocumentValidationError,
        TableValidationError,
        OSError,
        ValueError,
        json.JSONDecodeError,
    ) as error:
        parser.error(str(error))
    return 2


if __name__ == "__main__":
    raise SystemExit(main())
