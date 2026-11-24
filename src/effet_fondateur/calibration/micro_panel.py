"""Construction et validation du micro-panel transversal R10.

La release produite est entièrement synthétique et limitée à l'intégrité
technique. Les workers aveugles sont exécutés sous ``sandbox-exec`` sans accès
à ``private_truth/``, ``governance/`` ni aux données réelles du dépôt. Aucun
outil scientifique externe et aucune métrique de performance ne sont lancés.
"""

from __future__ import annotations

import argparse
from collections.abc import Iterable, Mapping, Sequence
from dataclasses import asdict, dataclass
import gzip
from hashlib import sha256
import io
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
from typing import Any

import numpy as np
from numpy.typing import NDArray

from effet_fondateur.calibration.contracts import (
    STAGE_IDS,
    TRUTH_CONSUMERS,
    VIEW_CONSUMERS,
    canonical_sha256,
    file_sha256,
    load_json,
    load_yaml,
    validate_document,
)
from effet_fondateur.calibration.micro_contract import (
    CORE_STAGES,
    STAGE_DEPENDENCIES,
    validate_micro_topology,
)
from effet_fondateur.calibration.readiness import validate_readiness_receipt
from effet_fondateur.calibration.release_guard import (
    partition_pseudonym,
    scan_blinded_document,
    validate_blinded_inventory,
    validate_partition_separation,
    validate_single_factor_pair,
)
from effet_fondateur.calibration.truth_router import (
    allowed_truth_for_stage,
    authorize_truth_access,
)


RELEASE_ID = "panel-master-generic-v1-r10-resource"
RELEASE_DIRECTORY = "micro_panel_generic_v1_r10"
CORE_RUN_ID = "RUN-R10-RESOURCE-CORE"
SENSITIVITY_RUN_ID = "RUN-R10-RESOURCE-SENSITIVITY"
RENDER_RUN_ID = "RUN-R10-RESOURCE-RENDER"
INTERPRET_RUN_ID = "RUN-R10-RESOURCE-INTERPRET"
FIXED_TIME = "2026-11-12T12:00:00Z"
SEED_PREFIX = "Effet_fondateur2|panel_master_generic_v1"
SEED_DOMAINS = (
    "SOURCE", "HAPLOTYPE", "POPULATION", "REFERENCE", "PEDIGREE",
    "MEIOSIS", "TARGET", "FOUNDER", "IBD", "ROH", "MEASUREMENT",
    "VIEW", "PHASING", "NULL_DRAW", "SENSITIVITY",
)


class MicroPanelError(RuntimeError):
    """Erreur bloquante de construction ou de validation R10."""

    def __init__(self, code: str, detail: str | None = None) -> None:
        self.code = code
        self.detail = detail
        super().__init__(code if detail is None else f"{code}: {detail}")


@dataclass(frozen=True)
class Subject:
    master_id: str
    family_id: str | None
    father_id: str | None
    mother_id: str | None
    component_type: str
    role: str
    population_state: str
    isolation_unit_id: str
    batch_id: str


@dataclass(frozen=True)
class Variant:
    chromosome: int
    position_bp: int
    variant_id: str
    ref: str
    alt: str
    position_cm: float
    role: str


@dataclass(frozen=True)
class MicroReleaseAudit:
    release_id: str
    subject_count: int
    reference_count: int
    variant_count: int
    autosome_count: int
    truth_family_count: int
    view_count: int
    stage_count: int
    evaluation_count: int
    markdown_count: int
    run_count: int
    sandbox_private_read_denied: bool
    external_scientific_tool_count: int
    html_pdf_count: int
    logical_bytes: int
    allocated_bytes: int


def _fail(code: str, detail: str | None = None) -> None:
    raise MicroPanelError(code, detail)


def _write_json(path: Path, document: Any) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(
        json.dumps(document, ensure_ascii=False, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )


def _write_tsv(path: Path, columns: Sequence[str], rows: Iterable[Sequence[Any]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    lines = ["\t".join(columns)]
    for row in rows:
        values = ["NA" if value is None else str(value) for value in row]
        if any("\t" in value or "\n" in value for value in values):
            _fail("unsafe_tsv_value", str(path))
        lines.append("\t".join(values))
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def _seed128(namespace: str, domain: str, entity_id: str, replicate_id: str = "0") -> int:
    payload = f"{SEED_PREFIX}|{namespace}|{domain}|{entity_id}|{replicate_id}"
    return int.from_bytes(sha256(payload.encode("utf-8")).digest()[:16], "big")


def _rng(domain: str, entity_id: str) -> np.random.Generator:
    return np.random.Generator(np.random.PCG64DXSM(_seed128("RESOURCE_ONLY", domain, entity_id)))


def _build_subjects() -> list[Subject]:
    subjects: list[Subject] = []
    populations = ("POP-A", "POP-B", "POP-C")
    batch_index = 0

    def add(
        identifier: str,
        family: str | None,
        father: str | None,
        mother: str | None,
        component: str,
        role: str,
        population: str,
    ) -> None:
        nonlocal batch_index
        batch_index += 1
        subjects.append(
            Subject(
                master_id=identifier,
                family_id=family,
                father_id=father,
                mother_id=mother,
                component_type=component,
                role=role,
                population_state=population,
                isolation_unit_id=family or identifier,
                batch_id=f"BATCH-{1 + (batch_index - 1) // 15}",
            )
        )

    for index in range(1, 5):
        family = f"TRIO-{index:02d}"
        population = populations[(index - 1) % 3]
        father, mother = f"{family}-F", f"{family}-M"
        add(father, family, None, None, "TRIO", "FATHER", population)
        add(mother, family, None, None, "TRIO", "MOTHER", population)
        add(f"{family}-C", family, father, mother, "TRIO", "CHILD", population)
    for index in range(1, 5):
        family = f"SIBSHIP-{index:02d}"
        population = populations[index % 3]
        father, mother = f"{family}-F", f"{family}-M"
        add(father, family, None, None, "SIBSHIP_4", "FATHER", population)
        add(mother, family, None, None, "SIBSHIP_4", "MOTHER", population)
        add(f"{family}-C1", family, father, mother, "SIBSHIP_4", "CHILD", population)
        add(f"{family}-C2", family, father, mother, "SIBSHIP_4", "CHILD", population)
    for index in range(1, 5):
        family = f"MULTIGEN-{index:02d}"
        population = populations[(index + 1) % 3]
        grandfather, grandmother = f"{family}-GF", f"{family}-GM"
        parent, spouse = f"{family}-P", f"{family}-S"
        add(grandfather, family, None, None, "MULTIGENERATIONAL_5", "GRANDFATHER", population)
        add(grandmother, family, None, None, "MULTIGENERATIONAL_5", "GRANDMOTHER", population)
        add(parent, family, grandfather, grandmother, "MULTIGENERATIONAL_5", "PARENT", population)
        add(spouse, family, None, None, "MULTIGENERATIONAL_5", "SPOUSE", population)
        add(f"{family}-C", family, parent, spouse, "MULTIGENERATIONAL_5", "CHILD", population)
    for index in range(1, 13):
        identifier = f"UNRELATED-{index:02d}"
        add(identifier, None, None, None, "UNRELATED", "INDEPENDENT", populations[(index - 1) % 3])
    if len(subjects) != 60:
        _fail("subject_count_mismatch", str(len(subjects)))
    return subjects


def _target_catalog(configuration: Mapping[str, Any]) -> list[dict[str, Any]]:
    target = configuration.get("target")
    catalog = target.get("catalog") if isinstance(target, Mapping) else None
    if not isinstance(catalog, list) or len(catalog) != 5:
        _fail("target_catalog_invalid")
    roles = {record.get("role") for record in catalog if isinstance(record, Mapping)}
    if roles != {"STANDARD", "LOW-RECOMB", "HIGH-RECOMB", "LEFT-EDGE", "RIGHT-EDGE"}:
        _fail("target_role_set_mismatch")
    return [dict(record) for record in catalog]


def _build_variants(configuration: Mapping[str, Any]) -> tuple[list[Variant], dict[str, Any]]:
    catalog = _target_catalog(configuration)
    selected_role = configuration["target"]["selected_role"]
    selected = next(record for record in catalog if record["role"] == selected_role)
    design = configuration["variant_design"]
    records: dict[tuple[int, int], dict[str, Any]] = {}
    allele_pairs = (("A", "G"), ("C", "T"), ("G", "A"), ("T", "C"))
    for chromosome in range(1, 23):
        for index, position in enumerate((1_000_000, 10_000_000, 20_000_000), start=1):
            ref, alt = allele_pairs[(chromosome + index) % len(allele_pairs)]
            records[(chromosome, position)] = {
                "chromosome": chromosome,
                "position_bp": position,
                "variant_id": f"SYN-CHR{chromosome:02d}-{position}",
                "ref": ref,
                "alt": alt,
                "role": "AUTOSOMAL_SENTINEL",
            }
    regional_chromosome = int(design["regional_chromosome"])
    for position in range(
        int(design["regional_start_bp"]),
        int(design["regional_end_bp"]) + 1,
        int(design["regional_step_bp"]),
    ):
        ref, alt = allele_pairs[(position // int(design["regional_step_bp"])) % 4]
        records[(regional_chromosome, position)] = {
            "chromosome": regional_chromosome,
            "position_bp": position,
            "variant_id": f"SYN-REG-{regional_chromosome}-{position}",
            "ref": ref,
            "alt": alt,
            "role": "REGIONAL_SENTINEL",
        }
    for target in catalog:
        records[(int(target["chromosome"]), int(target["position_bp"]))] = {
            **target,
            "role": f"TARGET_{target['role']}",
        }
    variants: list[Variant] = []
    for record in sorted(records.values(), key=lambda item: (int(item["chromosome"]), int(item["position_bp"]))):
        chromosome = int(record["chromosome"])
        position = int(record["position_bp"])
        rate_cm_per_mb = 0.7 + chromosome * 0.013
        variants.append(
            Variant(
                chromosome=chromosome,
                position_bp=position,
                variant_id=str(record["variant_id"]),
                ref=str(record["ref"]),
                alt=str(record["alt"]),
                position_cm=round((position / 1_000_000) * rate_cm_per_mb, 6),
                role=str(record["role"]),
            )
        )
    return variants, selected


def _build_haplotypes(
    subjects: Sequence[Subject],
    variants: Sequence[Variant],
    selected_target: Mapping[str, Any],
) -> NDArray[np.int8]:
    subject_index = {subject.master_id: index for index, subject in enumerate(subjects)}
    target_index = next(
        index for index, variant in enumerate(variants)
        if variant.variant_id == selected_target["variant_id"]
    )
    frequency_rng = _rng("HAPLOTYPE", "VARIANT_FREQUENCIES")
    frequencies = frequency_rng.uniform(0.05, 0.45, size=len(variants))
    haplotypes = np.zeros((len(subjects), 2, len(variants)), dtype=np.int8)
    first_parent_by_family: set[str] = set()
    for subject_index_value, subject in enumerate(subjects):
        if subject.father_id is None and subject.mother_id is None:
            generator = _rng("HAPLOTYPE", subject.master_id)
            haplotypes[subject_index_value] = (
                generator.random((2, len(variants))) < frequencies
            ).astype(np.int8)
            if subject.family_id and subject.family_id not in first_parent_by_family:
                haplotypes[subject_index_value, 0, target_index] = 1
                first_parent_by_family.add(subject.family_id)
            continue
        father_index = subject_index[subject.father_id]
        mother_index = subject_index[subject.mother_id]
        meiosis_rng = _rng("MEIOSIS", subject.master_id)
        paternal_copy = int(meiosis_rng.integers(0, 2))
        maternal_copy = int(meiosis_rng.integers(0, 2))
        haplotypes[subject_index_value, 0] = haplotypes[father_index, paternal_copy]
        haplotypes[subject_index_value, 1] = haplotypes[mother_index, maternal_copy]

    regional_indices = [
        index for index, variant in enumerate(variants)
        if variant.chromosome == int(selected_target["chromosome"])
        and 12_000_000 <= variant.position_bp <= 18_000_000
    ]
    carrier_indices = [
        index for index, subject in enumerate(subjects)
        if subject.family_id and haplotypes[index, :, target_index].sum() > 0
    ]
    if carrier_indices:
        founder_template = haplotypes[carrier_indices[0], 0, regional_indices].copy()
        for index in carrier_indices:
            carrier_copy = 0 if haplotypes[index, 0, target_index] else 1
            haplotypes[index, carrier_copy, regional_indices] = founder_template
            haplotypes[index, carrier_copy, target_index] = 1

    unrelated = [index for index, subject in enumerate(subjects) if subject.component_type == "UNRELATED"]
    ibd_indices = [
        index for index, variant in enumerate(variants)
        if variant.chromosome == int(selected_target["chromosome"])
        and 13_000_000 <= variant.position_bp <= 17_000_000
    ]
    for pair_index in range(0, min(10, len(unrelated) - 1), 2):
        donor = unrelated[pair_index]
        recipient = unrelated[pair_index + 1]
        haplotypes[recipient, 0, ibd_indices] = haplotypes[donor, 0, ibd_indices]
    roh_indices = [
        index for index, variant in enumerate(variants)
        if variant.chromosome == int(selected_target["chromosome"])
        and 11_000_000 <= variant.position_bp <= 19_000_000
    ]
    for index in unrelated[:4]:
        haplotypes[index, 1, roh_indices] = haplotypes[index, 0, roh_indices]
    return haplotypes


def _observations(haplotypes: NDArray[np.int8]) -> tuple[NDArray[np.int8], list[dict[str, Any]]]:
    truth = haplotypes.sum(axis=1).astype(np.int8)
    observed = truth.copy()
    events: list[dict[str, Any]] = []
    for event_index, (sample_index, variant_index, event_type) in enumerate(
        ((0, 0, "MISSING_CELL"), (1, 1, "GENOTYPE_FLIP"), (2, 2, "MISSING_SUBJECT_SENTINEL")),
        start=1,
    ):
        before = int(observed[sample_index, variant_index])
        after = -1 if "MISSING" in event_type else 2 - before
        observed[sample_index, variant_index] = after
        events.append(
            {
                "error_event_id": f"MEASURE-{event_index:03d}",
                "sample_index": sample_index,
                "variant_index": variant_index,
                "error_type": event_type,
                "before": before,
                "after": after,
            }
        )
    return observed, events


def _reference_pools() -> list[dict[str, Any]]:
    specifications = (
        ("REF-PHASE-RESOURCE-A", "PHASING", 45),
        ("REF-PHASE-RESOURCE-B", "PHASING", 45),
        ("REF-ANCESTRY-RESOURCE-A", "ANCESTRY", 45),
        ("REF-ANCESTRY-RESOURCE-B", "ANCESTRY", 45),
        ("REF-NULL-RESOURCE-A", "NULL", 30),
        ("REF-NULL-RESOURCE-B", "NULL", 30),
    )
    return [
        {
            "pool_id": pool_id,
            "reference_type": reference_type,
            "member_ids": [f"{pool_id}-S{index:03d}" for index in range(1, count + 1)],
        }
        for pool_id, reference_type, count in specifications
    ]


def _write_plink_dataset(
    root: Path,
    pseudonyms: Sequence[str],
    variants: Sequence[Variant],
    observed: NDArray[np.int8],
) -> list[Path]:
    prefix = root / "observed_autosomal"
    fam = prefix.with_suffix(".fam")
    bim = prefix.with_suffix(".bim")
    bed = prefix.with_suffix(".bed")
    fam.write_text(
        "".join(f"RESOURCE {sample} 0 0 0 -9\n" for sample in pseudonyms),
        encoding="utf-8",
    )
    bim.write_text(
        "".join(
            f"{variant.chromosome}\t{variant.variant_id}\t{variant.position_cm}\t"
            f"{variant.position_bp}\t{variant.alt}\t{variant.ref}\n"
            for variant in variants
        ),
        encoding="utf-8",
    )
    codes = {0: 0b11, 1: 0b10, 2: 0b00, -1: 0b01}
    payload = bytearray((0x6C, 0x1B, 0x01))
    for variant_index in range(len(variants)):
        genotypes = observed[:, variant_index]
        for offset in range(0, len(pseudonyms), 4):
            byte = 0
            for local_index, genotype in enumerate(genotypes[offset : offset + 4]):
                byte |= codes[int(genotype)] << (2 * local_index)
            payload.append(byte)
    bed.write_bytes(bytes(payload))
    return [bed, bim, fam]


def _write_vcf(
    path: Path,
    pseudonyms: Sequence[str],
    variants: Sequence[Variant],
    observed: NDArray[np.int8],
) -> None:
    genotype = {0: "0/0", 1: "0/1", 2: "1/1", -1: "./."}
    with path.open("wb") as raw_handle, gzip.GzipFile(
        filename="",
        mode="wb",
        fileobj=raw_handle,
        mtime=0,
    ) as compressed_handle, io.TextIOWrapper(
        compressed_handle,
        encoding="utf-8",
        newline="",
    ) as handle:
        handle.write("##fileformat=VCFv4.2\n")
        handle.write("##source=panel_master_generic_v1_R10_synthetic\n")
        handle.write("##reference=GRCh38\n")
        for chromosome in range(1, 23):
            handle.write(f"##contig=<ID={chromosome}>\n")
        handle.write("##FORMAT=<ID=GT,Number=1,Type=String,Description=Observed synthetic genotype>\n")
        handle.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t")
        handle.write("\t".join(pseudonyms) + "\n")
        for variant_index, variant in enumerate(variants):
            calls = [genotype[int(value)] for value in observed[:, variant_index]]
            handle.write(
                f"{variant.chromosome}\t{variant.position_bp}\t{variant.variant_id}\t"
                f"{variant.ref}\t{variant.alt}\t.\tPASS\t.\tGT\t"
                + "\t".join(calls)
                + "\n"
            )


def _truth_rows(
    truth_id: str,
    subjects: Sequence[Subject],
    variants: Sequence[Variant],
    pseudonyms: Sequence[str],
    haplotypes: NDArray[np.int8],
    observed: NDArray[np.int8],
    measurement_events: Sequence[Mapping[str, Any]],
    selected_target: Mapping[str, Any],
) -> tuple[tuple[str, ...], list[tuple[Any, ...]]]:
    target_index = next(index for index, variant in enumerate(variants) if variant.variant_id == selected_target["variant_id"])
    if truth_id == "T-SOURCE":
        return ("source_id", "valid", "assembly", "autosome_count"), [("SYNTHETIC-R10", True, "GRCh38", 22)]
    if truth_id == "T-SUBJECT":
        return (
            "master_id", "pseudonym", "family_id", "father_id", "mother_id",
            "component_type", "role", "population_state", "isolation_unit_id", "batch_id",
        ), [tuple(asdict(subject).get(column) for column in ("master_id",)) + (
            pseudonyms[index], subject.family_id, subject.father_id, subject.mother_id,
            subject.component_type, subject.role, subject.population_state,
            subject.isolation_unit_id, subject.batch_id,
        ) for index, subject in enumerate(subjects)]
    if truth_id == "T-GENOTYPE":
        rows = []
        for sample_index, subject in enumerate(subjects):
            for variant_index, variant in enumerate(variants):
                rows.append((subject.master_id, variant.variant_id, variant.chromosome, variant.position_bp, int(haplotypes[sample_index, :, variant_index].sum()), int(observed[sample_index, variant_index])))
        return ("master_id", "variant_id", "chromosome", "position_bp", "true_genotype", "observed_genotype"), rows
    if truth_id == "T-TARGET":
        return ("master_id", "variant_id", "haplotype_1", "haplotype_2", "true_genotype"), [
            (subject.master_id, selected_target["variant_id"], int(haplotypes[index, 0, target_index]), int(haplotypes[index, 1, target_index]), int(haplotypes[index, :, target_index].sum()))
            for index, subject in enumerate(subjects)
        ]
    if truth_id == "T-MEASURE":
        return ("error_event_id", "master_id", "variant_id", "error_type", "before", "after"), [
            (event["error_event_id"], subjects[int(event["sample_index"])].master_id, variants[int(event["variant_index"])].variant_id, event["error_type"], event["before"], event["after"])
            for event in measurement_events
        ]
    if truth_id == "T-KINSHIP":
        return ("pair_id", "master_id_1", "master_id_2", "relationship_class", "kinship_coefficient"), [
            ("PAIR-PARENT-CHILD", subjects[0].master_id, subjects[2].master_id, "FIRST_DEGREE", 0.25),
            ("PAIR-UNRELATED", subjects[-1].master_id, subjects[-2].master_id, "UNRELATED", 0.0),
        ]
    if truth_id == "T-POPULATION":
        return ("master_id", "population_state", "is_atypical"), [
            (subject.master_id, subject.population_state, index in {57, 58, 59})
            for index, subject in enumerate(subjects)
        ]
    if truth_id == "T-PHASE":
        rows = []
        for sample_index, subject in enumerate(subjects):
            for variant_index, variant in enumerate(variants):
                if variant.chromosome == int(selected_target["chromosome"]):
                    rows.append((subject.master_id, variant.variant_id, int(haplotypes[sample_index, 0, variant_index]), int(haplotypes[sample_index, 1, variant_index])))
        return ("master_id", "variant_id", "haplotype_1", "haplotype_2"), rows
    if truth_id == "T-FOUNDER":
        return ("event_id", "family_id", "age_generations", "start_bp", "end_bp", "censoring"), [
            (f"FND-{index:02d}", f"TRIO-{index:02d}", age, 12_000_000, 18_000_000, "NONE")
            for index, age in enumerate((5, 50, 75, 150), start=1)
        ]
    if truth_id == "T-LD":
        regional = [variant for variant in variants if variant.chromosome == int(selected_target["chromosome"])][:2]
        return ("variant_id_1", "variant_id_2", "r2", "d_prime"), [(regional[0].variant_id, regional[1].variant_id, 0.25, 0.5)]
    if truth_id == "T-ROH":
        return ("event_id", "master_id", "chromosome", "start_bp", "end_bp", "class", "contains_target"), [
            (f"ROH-{index:02d}", subjects[-index].master_id, selected_target["chromosome"], 11_000_000, 19_000_000, roh_class, True)
            for index, roh_class in enumerate(("SHORT", "MEDIUM", "LONG", "MULTIPLE_BURDEN"), start=1)
        ]
    if truth_id == "T-IBD":
        lengths = (1.0, 1.5, 2.0, 3.0, 5.0)
        return ("event_id", "master_id_1", "master_id_2", "chromosome", "start_bp", "end_bp", "length_cm", "is_positive"), [
            (f"IBD-{index:02d}", subjects[-(2 * index)].master_id, subjects[-(2 * index - 1)].master_id, selected_target["chromosome"], 13_000_000, 17_000_000, length, True)
            for index, length in enumerate(lengths, start=1)
        ] + [("IBD-NONE", subjects[-1].master_id, subjects[-2].master_id, selected_target["chromosome"], None, None, 0.0, False)]
    if truth_id == "T-SCENARIO":
        return ("pair_id", "domain", "factor_path", "left_value", "right_value", "expected_direction"), [
            (f"PAIR-{domain}", domain, "scenario.factor", "BASELINE", "CHALLENGE", "PRIVATE_EXPECTED")
            for domain in STAGE_DEPENDENCIES["17"]
        ]
    if truth_id == "T-RENDER":
        return ("stage_id", "expected_scope", "html_pdf_allowed", "founder_claim_allowed"), [
            (stage, "TECHNICAL_INTEGRITY", False, False) for stage in STAGE_IDS
        ]
    _fail("unknown_truth_family", truth_id)


def _materialize_truth(
    release_root: Path,
    subjects: Sequence[Subject],
    variants: Sequence[Variant],
    pseudonyms: Sequence[str],
    haplotypes: NDArray[np.int8],
    observed: NDArray[np.int8],
    measurement_events: Sequence[Mapping[str, Any]],
    selected_target: Mapping[str, Any],
) -> dict[str, Any]:
    artifacts = []
    for truth_id, consumers in TRUTH_CONSUMERS.items():
        columns, rows = _truth_rows(
            truth_id, subjects, variants, pseudonyms, haplotypes, observed,
            measurement_events, selected_target,
        )
        relative = Path("private_truth") / "RESOURCE_ONLY" / f"{truth_id.lower()}.tsv"
        path = release_root / relative
        _write_tsv(path, columns, rows)
        artifacts.append(
            {
                "truth_id": truth_id,
                "relative_path": relative.as_posix(),
                "sha256": file_sha256(path),
                "authorized_evaluators": [f"EVAL-{stage}" for stage in consumers],
            }
        )
    manifest = {
        "schema_version": "1.0.0",
        "release_id": RELEASE_ID,
        "disclosure_state": "SEALED_PRIVATE",
        "truth_artifacts": artifacts,
    }
    _write_json(release_root / "private_truth/master_truth_manifest.json", manifest)
    return manifest


def _view_payload(
    view_type: str,
    opaque_view_id: str,
    subjects: Sequence[Subject],
    variants: Sequence[Variant],
    pseudonyms: Sequence[str],
    selected_target: Mapping[str, Any],
) -> dict[str, Any]:
    base: dict[str, Any] = {
        "schema_version": "1.0.0",
        "opaque_view_id": opaque_view_id,
        "partition_id": "RESOURCE_ONLY",
        "synthetic": True,
    }
    if view_type == "V-SOURCE":
        base["observations"] = {"source_format": "SYNTHETIC_VCF_PLINK", "assembly": "GRCh38", "autosome_count": 22}
    elif view_type == "V-TARGET":
        base["observations"] = {"variant_id": selected_target["variant_id"], "chromosome": selected_target["chromosome"], "position_bp": selected_target["position_bp"], "ref": selected_target["ref"], "alt": selected_target["alt"], "call_source": "SYNTHETIC_MOLECULAR"}
    elif view_type == "V-DENSITY":
        base["observations"] = {"density_level": "MICRO_TECHNICAL", "variant_count": len(variants), "autosome_count": 22}
    elif view_type == "V-MEASUREMENT":
        base["observations"] = {"visible_missing_call_count": 2, "visible_nonmissing_call_count": len(subjects) * len(variants) - 2}
    elif view_type == "V-PEDIGREE":
        base["observations"] = {"declared_component_count": 12, "declared_subject_count": len(subjects)}
    elif view_type == "V-REGIONAL":
        base["observations"] = {"chromosome": selected_target["chromosome"], "start_bp": 10_000_000, "end_bp": 20_000_000, "target_variant_id": selected_target["variant_id"]}
    elif view_type == "V-AUTOSOMAL":
        base["observations"] = {"autosomes": list(range(1, 23)), "sample_count": len(subjects), "variant_count": len(variants), "same_subject_order": True}
    elif view_type == "V-PHASING-REFERENCE":
        base["observations"] = {"reference_type": "SYNTHETIC_PHASED", "reference_count": 90, "study_subject_donor_count": 0}
    elif view_type == "V-FOUNDER-NULL":
        base["observations"] = {"internal_unit_count": 12, "external_reference_count": 60, "null_draw_count": 99}
    elif view_type == "V-ROH-IBD":
        base["observations"] = {"genotype_source": "OBSERVED_SYNTHETIC", "phase_source": "ESTIMATED_PLACEHOLDER", "truth_labels_exposed": False}
    elif view_type == "V-SENSITIVITY":
        base["observations"] = {"pair_count": 16, "single_factor_declared": True, "expected_direction_exposed": False}
    return base


def _materialize_blinded_views(
    release_root: Path,
    subjects: Sequence[Subject],
    variants: Sequence[Variant],
    pseudonyms: Sequence[str],
    observed: NDArray[np.int8],
    selected_target: Mapping[str, Any],
    parent_sha256: str,
) -> tuple[dict[str, Any], dict[str, Any]]:
    blind_root = release_root / "blinded_inputs"
    partition_root = blind_root / "RESOURCE_ONLY"
    shared_root = partition_root / "shared"
    shared_root.mkdir(parents=True)
    shared_paths = _write_plink_dataset(shared_root, pseudonyms, variants, observed)
    vcf = shared_root / "observed_autosomal.vcf.gz"
    _write_vcf(vcf, pseudonyms, variants, observed)
    shared_paths.append(vcf)
    reference = shared_root / "synthetic_reference_summary.tsv"
    _write_tsv(
        reference,
        ("reference_type", "subject_count", "phased", "study_overlap_count"),
        (("PHASING", 90, True, 0), ("ANCESTRY", 90, True, 0), ("NULL", 60, True, 0)),
    )
    shared_paths.append(reference)

    view_records = []
    artifact_paths: list[Path] = list(shared_paths)
    for index, (view_type, consumers) in enumerate(VIEW_CONSUMERS.items(), start=1):
        opaque = f"VIEW-{index:016x}"
        view_root = partition_root / opaque
        view_root.mkdir()
        view_path = view_root / "view.json"
        _write_json(
            view_path,
            _view_payload(view_type, opaque, subjects, variants, pseudonyms, selected_target),
        )
        artifact_paths.append(view_path)
        view_records.append(
            {
                "view_type": view_type,
                "opaque_view_id": opaque,
                "parent_sha256": parent_sha256,
                "transformation_id": f"derive-{view_type.lower()}-v1",
                "consumer_stages": list(consumers),
                "inventory_path": "blinded_inputs/blinded_input_inventory.json",
            }
        )
    inventory = {
        "schema_version": "1.0.0",
        "inventory_id": "INVENTORY-R10-RESOURCE-001",
        "release_id": RELEASE_ID,
        "partition_id": "RESOURCE_ONLY",
        "opaque_view_id": "VIEW-ffffffffffffffff",
        "artifacts": [
            {
                "artifact_id": f"BLIND-{index:03d}",
                "relative_path": path.relative_to(blind_root).as_posix(),
                "sha256": file_sha256(path),
                "size_bytes": path.stat().st_size,
                "media_type": (
                    "application/json" if path.suffix == ".json"
                    else "text/tab-separated-values" if path.suffix in {".tsv", ".bim", ".fam"}
                    else "application/octet-stream"
                ),
            }
            for index, path in enumerate(sorted(artifact_paths), start=1)
        ],
    }
    _write_json(blind_root / "blinded_input_inventory.json", inventory)
    view_manifest = {
        "schema_version": "1.0.0",
        "release_id": RELEASE_ID,
        "partition_id": "RESOURCE_ONLY",
        "views": view_records,
    }
    _write_json(release_root / "views/panel_view_manifest.json", view_manifest)
    return inventory, view_manifest


def _materialize_governance(
    repository_root: Path,
    release_root: Path,
    configuration: Mapping[str, Any],
    subjects: Sequence[Subject],
    variants: Sequence[Variant],
    pseudonyms: Sequence[str],
    reference_pools: Sequence[Mapping[str, Any]],
    selected_target: Mapping[str, Any],
    truth_manifest: Mapping[str, Any],
    inventory: Mapping[str, Any],
) -> None:
    governance = release_root / "governance"
    governance.mkdir(exist_ok=True)
    source_catalog = {
        "schema_version": "1.0.0",
        "release_id": RELEASE_ID,
        "sources": [
            {
                "source_id": "SYNTHETIC-GENERATED-R10",
                "source_type": "FULLY_SYNTHETIC",
                "assembly": "GRCh38",
                "real_data": False,
                "license": "INTERNAL_SYNTHETIC",
            }
        ],
    }
    _write_json(governance / "source_catalog.json", source_catalog)
    _write_tsv(
        governance / "seed_ledger.tsv",
        ("partition_namespace", "domain", "entity_id", "replicate_id", "seed128"),
        (("RESOURCE_ONLY", domain, "MASTER", "0", _seed128("RESOURCE_ONLY", domain, "MASTER")) for domain in SEED_DOMAINS),
    )
    _write_tsv(
        governance / "subject_registry.tsv",
        ("master_id", "pseudonym", "partition_id", "isolation_unit_id", "family_id", "role", "batch_id", "activation_rank"),
        (
            (subject.master_id, pseudonyms[index], "RESOURCE_ONLY", subject.isolation_unit_id, subject.family_id, subject.role, subject.batch_id, index + 1)
            for index, subject in enumerate(subjects)
        ),
    )
    subject_contract = {
        "schema_version": "1.0.0",
        "release_id": RELEASE_ID,
        "relative_path": "governance/subject_registry.tsv",
        "sha256": file_sha256(governance / "subject_registry.tsv"),
        "row_count": len(subjects),
        "columns": ["master_id", "pseudonym", "partition_id", "isolation_unit_id", "family_id", "role", "batch_id", "activation_rank"],
        "partition_values": ["CALIBRATION", "VALIDATION_BLIND", "RESOURCE_ONLY"],
    }
    _write_json(governance / "subject_registry.contract.json", subject_contract)
    _write_json(governance / "reference_pools.json", {"release_id": RELEASE_ID, "pools": list(reference_pools)})
    receipt_source = repository_root / "calibration/readiness/R09_READY_TO_RUN.json"
    shutil.copyfile(receipt_source, governance / "readiness_receipt.json")
    release_manifest = {
        "schema_version": "1.0.0",
        "release_id": RELEASE_ID,
        "release_design": "panel_master_generic_v1",
        "assembly": "GRCh38",
        "stage_ids": list(STAGE_IDS),
        "partitions": ["CALIBRATION", "VALIDATION_BLIND", "RESOURCE_ONLY"],
        "scenario_design": {
            "scale": "RESOURCE_ONLY",
            "study_identities": len(subjects),
            "phasing_reference": 90,
            "ancestry_reference": 90,
            "null_reference": 60,
            "variant_count": len(variants),
            "autosomes": 22,
            "performance_decision_allowed": False,
            "configuration_sha256": canonical_sha256(configuration),
            "truth_manifest_sha256": canonical_sha256(truth_manifest),
            "blinded_inventory_sha256": canonical_sha256(inventory),
        },
        "target": {"selection": "PROFILE_CONFIGURED_AUTOSOMAL", "hard_coded": False},
        "qualification_flags": {
            "calibrated": False,
            "founder_effect_proven": False,
            "production_promotion": False,
        },
    }
    _write_json(governance / "panel_release_manifest.json", release_manifest)
    _write_json(
        governance / "target_resolution.json",
        {
            "release_id": RELEASE_ID,
            "selection_rule": "selected_role from micro_panel_v1.yaml",
            "selected_target": dict(selected_target),
            "catalog_size": 5,
            "hard_coded_in_engine": False,
        },
    )


def _separation_inputs(
    subjects: Sequence[Subject],
    reference_pools: Sequence[Mapping[str, Any]],
) -> tuple[list[dict[str, Any]], list[Mapping[str, Any]]]:
    family_ids = sorted({subject.family_id for subject in subjects if subject.family_id})
    assignments = [
        {
            "partition_id": "RESOURCE_ONLY",
            "isolation_unit_ids": sorted({subject.isolation_unit_id for subject in subjects}),
            "master_identity_ids": [subject.master_id for subject in subjects],
            "pedigree_component_ids": family_ids,
            "event_ids": [f"FND-{index:02d}" for index in range(1, 5)] + [f"IBD-{index:02d}" for index in range(1, 6)] + [f"ROH-{index:02d}" for index in range(1, 5)],
            "segment_ids": [f"SEG-{index:02d}" for index in range(1, 14)],
            "duplicate_group_ids": ["DUP-TECH-01"],
            "seed_prefixes": [f"RESOURCE-{domain}" for domain in SEED_DOMAINS],
            "negative_control_ids": [f"NEG-{index:02d}" for index in range(1, 13)],
        }
    ]
    return assignments, list(reference_pools)


def _sandbox_profile(repository_root: Path, release_root: Path) -> str:
    denied_read = (
        release_root / "private_truth",
        release_root / "governance",
        repository_root / "data/input",
    )
    denied_write = (
        release_root / "private_truth",
        release_root / "governance",
        release_root / "blinded_inputs",
        repository_root / "data",
    )
    lines = ["(version 1)", "(allow default)"]
    lines.extend(f'(deny file-read* (subpath "{path}"))' for path in denied_read)
    lines.extend(f'(deny file-write* (subpath "{path}"))' for path in denied_write)
    return "\n".join(lines)


def _run_worker(
    repository_root: Path,
    release_root: Path,
    input_root: str,
    run_root: str,
    stage_id: str,
) -> None:
    profile = _sandbox_profile(repository_root, release_root)
    command = [
        "/usr/bin/sandbox-exec",
        "-p",
        profile,
        sys.executable,
        "-m",
        "effet_fondateur.calibration.micro_stage_worker",
        "--release-root",
        str(release_root),
        "--input-root",
        input_root,
        "--run-root",
        run_root,
        "--stage-id",
        stage_id,
        "--release-id",
        RELEASE_ID,
    ]
    completed = subprocess.run(
        command,
        cwd=repository_root,
        env={"PATH": os.environ.get("PATH", ""), "PYTHONHASHSEED": "0", "TZ": "UTC"},
        text=True,
        capture_output=True,
        check=False,
    )
    if completed.returncode != 0:
        _fail("sandboxed_stage_failed", f"{stage_id}:{completed.stderr.strip()}")


def _prove_private_denial(repository_root: Path, release_root: Path) -> bool:
    profile = _sandbox_profile(repository_root, release_root)
    private_path = release_root / "private_truth/master_truth_manifest.json"
    completed = subprocess.run(
        [
            "/usr/bin/sandbox-exec", "-p", profile, sys.executable, "-m",
            "effet_fondateur.calibration.micro_stage_worker",
            "--probe-private", str(private_path),
        ],
        cwd=repository_root,
        env={"PATH": os.environ.get("PATH", ""), "PYTHONHASHSEED": "0", "TZ": "UTC"},
        text=True,
        capture_output=True,
        check=False,
    )
    return completed.returncode == 0


def _close_run(
    release_root: Path,
    run_id: str,
    stage_ids: Sequence[str],
    inventory_sha256: str,
    dependencies: Mapping[str, Any],
) -> dict[str, Any]:
    run_root = release_root / "runs" / run_id
    outputs: list[dict[str, Any]] = []
    for stage_id in stage_ids:
        stage_manifest = load_json(run_root / stage_id / "stage_outputs.json")
        for artifact in stage_manifest["artifacts"]:
            path = run_root / artifact["relative_path"]
            outputs.append(
                {
                    "stage_id": stage_id,
                    "artifact_id": artifact["artifact_id"],
                    "relative_path": path.relative_to(release_root).as_posix(),
                    "sha256": file_sha256(path),
                }
            )
    manifest: dict[str, Any] = {
        "schema_version": "1.0.0",
        "run_id": run_id,
        "release_id": RELEASE_ID,
        "partition_id": "RESOURCE_ONLY",
        "run_state": "CLOSED",
        "closed_at": FIXED_TIME,
        "manifest_sha256": "0" * 64,
        "consumed_inventory_sha256": inventory_sha256,
        "stage_outputs": outputs,
        "dependencies": dict(dependencies),
    }
    manifest["manifest_sha256"] = canonical_sha256(
        {**manifest, "manifest_sha256": None}
    )
    _write_json(run_root / "pipeline_run_manifest.json", manifest)
    return manifest


def _sensitivity_sources(release_root: Path, core_manifest: Mapping[str, Any]) -> dict[str, Any]:
    source_root = release_root / "runs/sensitivity_sources"
    source_root.mkdir(parents=True)
    rows = []
    dependency_index: dict[str, Any] = {}
    core_outputs_by_stage: dict[str, list[Mapping[str, Any]]] = {}
    for output in core_manifest["stage_outputs"]:
        core_outputs_by_stage.setdefault(output["stage_id"], []).append(output)
    for domain in STAGE_DEPENDENCIES["17"]:
        baseline = {
            "run_id": f"PAIR-{domain}-BASELINE",
            "release_id": RELEASE_ID,
            "run_state": "CLOSED",
            "closed_at": FIXED_TIME,
            "manifest_sha256": sha256(f"{domain}-baseline".encode()).hexdigest(),
            "scenario": {"factor": "BASELINE"},
            "source_stage": domain,
        }
        challenge = {**baseline, "run_id": f"PAIR-{domain}-CHALLENGE", "manifest_sha256": sha256(f"{domain}-challenge".encode()).hexdigest(), "scenario": {"factor": "CHALLENGE"}}
        validate_single_factor_pair(
            baseline,
            challenge,
            allowed_path="scenario.factor",
            ignored_paths=("run_id", "manifest_sha256"),
        )
        baseline_path = source_root / f"{domain}.baseline.json"
        challenge_path = source_root / f"{domain}.challenge.json"
        _write_json(baseline_path, baseline)
        _write_json(challenge_path, challenge)
        rows.append((f"PAIR-{domain}", domain, baseline_path.relative_to(release_root).as_posix(), file_sha256(baseline_path), challenge_path.relative_to(release_root).as_posix(), file_sha256(challenge_path), "scenario.factor", "SIGNED_TECHNICAL"))
        stage_outputs = core_outputs_by_stage.get(domain, [])
        if not stage_outputs:
            _fail("sensitivity_source_stage_missing", domain)
        dependency_index[domain] = {
            "relative_path": stage_outputs[0]["relative_path"],
            "sha256": stage_outputs[0]["sha256"],
        }
    registry = release_root / "governance/sensitivity_run_registry.tsv"
    _write_tsv(
        registry,
        ("pair_id", "domain", "baseline_manifest", "baseline_sha256", "challenge_manifest", "challenge_sha256", "changed_factor", "signature_state"),
        rows,
    )
    sensitivity_input = release_root / "derived_inputs/sensitivity"
    sensitivity_input.mkdir(parents=True)
    _write_json(sensitivity_input / "dependency_index.json", dependency_index)
    _write_json(
        sensitivity_input / "signed_pair_registry_receipt.json",
        {
            "pair_count": 16,
            "domains": list(STAGE_DEPENDENCIES["17"]),
            "registry_sha256": file_sha256(registry),
            "expected_direction_exposed": False,
        },
    )
    _write_derived_inventory(sensitivity_input)
    return dependency_index


def _evaluate_stage(
    repository_root: Path,
    release_root: Path,
    stage_id: str,
    run_manifest: Mapping[str, Any],
    truth_manifest: Mapping[str, Any],
) -> dict[str, Any]:
    truth_ids = sorted(allowed_truth_for_stage(stage_id))
    grant = authorize_truth_access(
        evaluation_id=f"EVAL-{stage_id}-R10",
        stage_id=stage_id,
        requested_truth_ids=truth_ids,
        truth_manifest=truth_manifest,
        run_manifest=run_manifest,
    )
    stage_outputs = [
        output for output in run_manifest["stage_outputs"] if output["stage_id"] == stage_id
    ]
    if not stage_outputs:
        _fail("evaluation_stage_output_missing", stage_id)
    for output in stage_outputs:
        path = release_root / output["relative_path"]
        if file_sha256(path) != output["sha256"]:
            _fail("closed_output_hash_mismatch", stage_id)
    evaluation_root = release_root / "evaluations" / stage_id
    evaluation_root.mkdir(parents=True, exist_ok=True)
    metrics = evaluation_root / f"ETAPE_{stage_id}.metrics.tsv"
    _write_tsv(
        metrics,
        ("metric_id", "value", "unit", "decision_scope", "scientific_performance_evaluated"),
        (("declared_closed_outputs_present", 1, "boolean", "TECHNICAL_INTEGRITY", False),),
    )
    evaluation = {
        "schema_version": "1.0.0",
        "evaluation_id": f"EVAL-{stage_id}-R10",
        "stage_id": stage_id,
        "release_id": RELEASE_ID,
        "closed_run_id": run_manifest["run_id"],
        "truth_grants": truth_ids,
        "metrics_path": metrics.relative_to(release_root).as_posix(),
        "metrics_sha256": file_sha256(metrics),
        "status": "PASS",
        "decision_scope": "TECHNICAL_INTEGRITY",
        "scientific_performance_evaluated": False,
        "technical_checks": [
            "RUN_CLOSED_BEFORE_TRUTH_OPENING",
            "OUTPUTS_DECLARED_AND_HASHED",
            "TRUTH_SCOPE_LEAST_PRIVILEGE",
            "NO_REAL_DATA",
            "NO_EXTERNAL_SCIENTIFIC_TOOL",
        ],
        "counts": {"N_EXPECTED": 1, "N_ELIGIBLE": 1, "N_EVALUABLE": 1, "N_EXCLUDED": 0},
        "criteria_version": "1.0.0",
        "closed": True,
    }
    evaluation_path = evaluation_root / f"ETAPE_{stage_id}.evaluation.json"
    _write_json(evaluation_path, evaluation)
    validate_document(
        evaluation,
        repository_root / "calibration/schemas/stage_evaluation_manifest.schema.json",
    )
    result = release_root / "calibration/resultats" / f"ETAPE_{stage_id}.md"
    truth_text = ", ".join(truth_ids) if truth_ids else "aucune vérité privée"
    result.parent.mkdir(parents=True, exist_ok=True)
    result.write_text(
        f"# Étape {stage_id} — micro-panel R10\n\n"
        f"- Release : `{RELEASE_ID}`\n"
        f"- Run clos : `{run_manifest['run_id']}`\n"
        "- Partition : `RESOURCE_ONLY`\n"
        "- Portée : intégrité technique uniquement\n"
        "- Performance scientifique : non évaluée\n"
        f"- Vérité ouverte après clôture : {truth_text}\n"
        f"- Sorties closes vérifiées : {len(stage_outputs)}\n"
        "- Décision technique : `PASS`\n"
        "- Conclusion scientifique : aucune\n"
        "- Effet fondateur démontré : non\n"
        "- Promotion en production : non\n\n"
        "La décision confirme seulement que les entrées, sorties, empreintes, "
        "gardes et dépendances techniques de cette étape sont traversables. "
        "Les seuils R08 ne sont pas évalués sur le micro-panel.\n",
        encoding="utf-8",
    )
    return evaluation


def _derived_dependency_index(
    release_root: Path,
    destination: Path,
    stage_ids: Sequence[str],
) -> None:
    destination.mkdir(parents=True)
    evidence_root = destination / "evidence"
    evidence_root.mkdir()
    document = {}
    for stage_id in stage_ids:
        evaluation = release_root / "evaluations" / stage_id / f"ETAPE_{stage_id}.evaluation.json"
        metrics = release_root / "evaluations" / stage_id / f"ETAPE_{stage_id}.metrics.tsv"
        if not evaluation.is_file():
            _fail("derived_evaluation_missing", stage_id)
        copied_evaluation = evidence_root / f"ETAPE_{stage_id}.evaluation.json"
        copied_metrics = evidence_root / f"ETAPE_{stage_id}.metrics.tsv"
        shutil.copyfile(evaluation, copied_evaluation)
        shutil.copyfile(metrics, copied_metrics)
        document[stage_id] = {
            "relative_path": copied_evaluation.relative_to(destination).as_posix(),
            "sha256": file_sha256(copied_evaluation),
            "metrics_relative_path": copied_metrics.relative_to(destination).as_posix(),
            "metrics_sha256": file_sha256(copied_metrics),
        }
    _write_json(destination / "dependency_index.json", document)


def _write_derived_inventory(input_root: Path) -> dict[str, Any]:
    """Ferme et empreinte une entrée aval dépourvue de vérité privée."""
    inventory_path = input_root / "input_inventory.json"
    paths = sorted(path for path in input_root.rglob("*") if path.is_file() and path != inventory_path)
    document = {
        "schema_version": "1.0.0",
        "input_id": f"DERIVED-{input_root.name.upper()}",
        "closed": True,
        "artifacts": [
            {
                "relative_path": path.relative_to(input_root).as_posix(),
                "sha256": file_sha256(path),
                "size_bytes": path.stat().st_size,
            }
            for path in paths
        ],
    }
    scan_blinded_document(document)
    _write_json(inventory_path, document)
    return document


def _materialize_stage_chain(
    repository_root: Path,
    release_root: Path,
    inventory: Mapping[str, Any],
    truth_manifest: Mapping[str, Any],
) -> dict[str, dict[str, Any]]:
    inventory_hash = canonical_sha256(inventory)
    runs_root = release_root / "runs"
    runs_root.mkdir()
    core_root = runs_root / CORE_RUN_ID
    core_root.mkdir()
    for stage_id in tuple(stage for stage in CORE_STAGES if stage != "17"):
        _run_worker(
            repository_root,
            release_root,
            "blinded_inputs",
            f"runs/{CORE_RUN_ID}",
            stage_id,
        )
    core_stage_ids = tuple(stage for stage in CORE_STAGES if stage != "17")
    core_manifest = _close_run(
        release_root,
        CORE_RUN_ID,
        core_stage_ids,
        inventory_hash,
        {"16C_ancestry_control": "runs/RUN-R10-RESOURCE-CORE/16A/reference_ancestry_summary.json"},
    )

    _sensitivity_sources(release_root, core_manifest)
    sensitivity_root = runs_root / SENSITIVITY_RUN_ID
    sensitivity_root.mkdir()
    _run_worker(
        repository_root,
        release_root,
        "derived_inputs/sensitivity",
        f"runs/{SENSITIVITY_RUN_ID}",
        "17",
    )
    sensitivity_manifest = _close_run(
        release_root,
        SENSITIVITY_RUN_ID,
        ("17",),
        canonical_sha256(
            load_json(release_root / "derived_inputs/sensitivity/input_inventory.json")
        ),
        {"signed_registry": "governance/sensitivity_run_registry.tsv", "domain_count": 16},
    )

    run_by_stage = {stage: core_manifest for stage in core_stage_ids}
    run_by_stage["17"] = sensitivity_manifest
    for stage_id in tuple(stage for stage in STAGE_IDS if stage not in {"18", "19"}):
        _evaluate_stage(repository_root, release_root, stage_id, run_by_stage[stage_id], truth_manifest)

    render_input = release_root / "derived_inputs/render"
    upstream_for_render = tuple(stage for stage in STAGE_IDS if stage not in {"18", "19"})
    _derived_dependency_index(release_root, render_input, upstream_for_render)
    render_inventory = _write_derived_inventory(render_input)
    render_root = runs_root / RENDER_RUN_ID
    render_root.mkdir()
    _run_worker(repository_root, release_root, "derived_inputs/render", f"runs/{RENDER_RUN_ID}", "18")
    render_manifest = _close_run(
        release_root,
        RENDER_RUN_ID,
        ("18",),
        canonical_sha256(render_inventory),
        {"evaluation_count": len(upstream_for_render), "html_pdf_allowed": False},
    )
    _evaluate_stage(repository_root, release_root, "18", render_manifest, truth_manifest)

    interpretation_input = release_root / "derived_inputs/interpretation"
    upstream_for_interpretation = tuple(stage for stage in STAGE_IDS if stage != "19")
    _derived_dependency_index(release_root, interpretation_input, upstream_for_interpretation)
    markdown_index = {
        stage: {
            "relative_path": f"calibration/resultats/ETAPE_{stage}.md",
            "sha256": file_sha256(release_root / f"calibration/resultats/ETAPE_{stage}.md"),
        }
        for stage in upstream_for_interpretation
    }
    _write_json(interpretation_input / "markdown_index.json", markdown_index)
    markdown_input = interpretation_input / "markdown"
    markdown_input.mkdir()
    for stage in upstream_for_interpretation:
        shutil.copyfile(
            release_root / f"calibration/resultats/ETAPE_{stage}.md",
            markdown_input / f"ETAPE_{stage}.md",
        )
    interpretation_inventory = _write_derived_inventory(interpretation_input)
    interpretation_root = runs_root / INTERPRET_RUN_ID
    interpretation_root.mkdir()
    _run_worker(repository_root, release_root, "derived_inputs/interpretation", f"runs/{INTERPRET_RUN_ID}", "19")
    interpretation_manifest = _close_run(
        release_root,
        INTERPRET_RUN_ID,
        ("19",),
        canonical_sha256(interpretation_inventory),
        {"upstream_markdown_count": 24, "html_pdf_allowed": False},
    )
    _evaluate_stage(repository_root, release_root, "19", interpretation_manifest, truth_manifest)
    return {
        CORE_RUN_ID: core_manifest,
        SENSITIVITY_RUN_ID: sensitivity_manifest,
        RENDER_RUN_ID: render_manifest,
        INTERPRET_RUN_ID: interpretation_manifest,
    }


def build_micro_release(repository_root: Path, release_root: Path) -> MicroReleaseAudit:
    """Construit une nouvelle release R10 dans un répertoire absent."""
    repository_root = repository_root.resolve(strict=True)
    release_root = release_root.resolve(strict=False)
    if release_root.exists():
        _fail("release_target_exists", str(release_root))
    workspace_root = release_root.parent.parent
    if workspace_root != repository_root / "calibration_workspace":
        _fail("release_target_outside_calibration_workspace", str(release_root))
    validate_readiness_receipt(
        repository_root,
        repository_root / "calibration/readiness/R09_READY_TO_RUN.json",
    )
    validate_micro_topology()
    configuration = load_yaml(repository_root / "calibration/configurations/micro_panel_v1.yaml")
    release_root.mkdir(parents=True)
    for directory in (
        "governance", "private_truth/RESOURCE_ONLY", "views", "blinded_inputs",
        "evaluations", "calibration/resultats", "derived_inputs",
    ):
        (release_root / directory).mkdir(parents=True, exist_ok=True)

    subjects = _build_subjects()
    variants, selected_target = _build_variants(configuration)
    secret = sha256(f"{RELEASE_ID}|RESOURCE_ONLY|PSEUDONYM".encode()).digest()
    pseudonyms = [partition_pseudonym(secret, "RESOURCE_ONLY", subject.master_id) for subject in subjects]
    haplotypes = _build_haplotypes(subjects, variants, selected_target)
    observed, measurement_events = _observations(haplotypes)
    reference_pools = _reference_pools()
    assignments, pool_contracts = _separation_inputs(subjects, reference_pools)
    separation = validate_partition_separation(assignments, pool_contracts)

    truth_manifest = _materialize_truth(
        release_root, subjects, variants, pseudonyms, haplotypes, observed,
        measurement_events, selected_target,
    )
    inventory, view_manifest = _materialize_blinded_views(
        release_root, subjects, variants, pseudonyms, observed, selected_target,
        canonical_sha256(configuration),
    )
    _materialize_governance(
        repository_root, release_root, configuration, subjects, variants,
        pseudonyms, reference_pools, selected_target, truth_manifest, inventory,
    )
    _write_json(
        release_root / "governance/separation_audit.json",
        {**asdict(separation), "status": "PASS", "overlap_count": 0},
    )
    _write_json(
        release_root / "governance/separation_contract.json",
        {"assignments": assignments, "reference_pools": list(pool_contracts)},
    )

    schema_root = repository_root / "calibration/schemas"
    validate_document(load_json(release_root / "governance/panel_release_manifest.json"), schema_root / "panel_release_manifest.schema.json")
    validate_document(load_json(release_root / "governance/subject_registry.contract.json"), schema_root / "subject_registry.schema.json")
    validate_document(dict(truth_manifest), schema_root / "master_truth_manifest.schema.json")
    validate_document(dict(view_manifest), schema_root / "panel_view_manifest.schema.json")
    validate_document(dict(inventory), schema_root / "blinded_input_inventory.schema.json")
    validate_blinded_inventory(
        inventory,
        view_root=release_root / "blinded_inputs",
        repository_root=repository_root,
        verify_files=True,
    )
    for path in (release_root / "blinded_inputs").rglob("*.json"):
        scan_blinded_document(load_json(path))
    if not _prove_private_denial(repository_root, release_root):
        _fail("sandbox_private_read_not_denied")

    run_manifests = _materialize_stage_chain(
        repository_root, release_root, inventory, truth_manifest
    )
    for manifest in run_manifests.values():
        validate_document(manifest, schema_root / "pipeline_run_manifest.schema.json")
    audit = validate_micro_release(
        repository_root,
        release_root,
        require_release_checksum_index=False,
    )
    receipt_audit = asdict(audit)
    _write_json(
        release_root / "governance/R10_MICRO_PANEL_RECEIPT.json",
        {
            "schema_version": "1.0.0",
            "release_id": RELEASE_ID,
            "r10_gate": "PASS",
            "scope": "TECHNICAL_INTEGRITY_ONLY",
            "audit_measurement": "BEFORE_RECEIPT_FILE",
            "audit": receipt_audit,
            "qualification_flags": {
                "calibrated": False,
                "founder_effect_proven": False,
                "production_promotion": False,
            },
            "next_authorized_phase": "R11_PILOT_CALIBRATION_AFTER_NEW_REQUEST",
        },
    )
    _write_release_checksum_index(release_root)
    return validate_micro_release(repository_root, release_root)


def _read_tsv_header(path: Path) -> list[str]:
    return path.read_text(encoding="utf-8").splitlines()[0].split("\t")


def _verify_derived_inventory(input_root: Path) -> dict[str, Any]:
    inventory_path = input_root / "input_inventory.json"
    document = load_json(inventory_path)
    scan_blinded_document(document)
    declared = {artifact["relative_path"]: artifact for artifact in document["artifacts"]}
    actual = {
        path.relative_to(input_root).as_posix(): path
        for path in input_root.rglob("*")
        if path.is_file() and path != inventory_path
    }
    if set(declared) != set(actual):
        _fail("derived_inventory_membership_mismatch", input_root.name)
    for relative_path, path in actual.items():
        artifact = declared[relative_path]
        if artifact["sha256"] != file_sha256(path) or artifact["size_bytes"] != path.stat().st_size:
            _fail("derived_inventory_integrity_failed", relative_path)
    return document


def _verify_stage_envelope(stage_root: Path) -> None:
    required = {
        "stage_inputs.json",
        "command.json",
        "stdout.log",
        "stderr.log",
        "stage_outputs.json",
        "audit.json",
        "checksums.sha256",
    }
    actual_names = {path.name for path in stage_root.iterdir() if path.is_file()}
    if not required.issubset(actual_names):
        _fail("stage_envelope_incomplete", stage_root.name)
    checksum_lines = (stage_root / "checksums.sha256").read_text(encoding="utf-8").splitlines()
    declared: dict[str, str] = {}
    for line in checksum_lines:
        digest, separator, name = line.partition("  ")
        if not separator or name in declared:
            _fail("stage_checksum_index_invalid", stage_root.name)
        declared[name] = digest
    expected_files = {
        path.name: path
        for path in stage_root.iterdir()
        if path.is_file() and path.name != "checksums.sha256"
    }
    if set(declared) != set(expected_files):
        _fail("stage_checksum_membership_mismatch", stage_root.name)
    for name, path in expected_files.items():
        if declared[name] != file_sha256(path):
            _fail("stage_checksum_mismatch", f"{stage_root.name}:{name}")


def _write_release_checksum_index(release_root: Path) -> None:
    index_path = release_root / "governance/release_checksums.sha256"
    files = sorted(
        path for path in release_root.rglob("*")
        if path.is_file() and path != index_path
    )
    index_path.write_text(
        "".join(
            f"{file_sha256(path)}  {path.relative_to(release_root).as_posix()}\n"
            for path in files
        ),
        encoding="utf-8",
    )


def _verify_release_checksum_index(release_root: Path) -> None:
    index_path = release_root / "governance/release_checksums.sha256"
    if not index_path.is_file():
        _fail("release_checksum_index_missing")
    declared: dict[str, str] = {}
    for line in index_path.read_text(encoding="utf-8").splitlines():
        digest, separator, relative_path = line.partition("  ")
        if not separator or relative_path in declared:
            _fail("release_checksum_index_invalid")
        declared[relative_path] = digest
    actual = {
        path.relative_to(release_root).as_posix(): path
        for path in release_root.rglob("*")
        if path.is_file() and path != index_path
    }
    if set(declared) != set(actual):
        _fail("release_checksum_membership_mismatch")
    for relative_path, path in actual.items():
        if declared[relative_path] != file_sha256(path):
            _fail("release_checksum_mismatch", relative_path)


def validate_micro_release(
    repository_root: Path,
    release_root: Path,
    *,
    require_release_checksum_index: bool = True,
) -> MicroReleaseAudit:
    """Rouvre une release R10 et vérifie sa porte sans la modifier."""
    repository_root = repository_root.resolve(strict=True)
    release_root = release_root.resolve(strict=True)
    schema_root = repository_root / "calibration/schemas"
    required_roots = (
        "governance/panel_release_manifest.json",
        "governance/source_catalog.json",
        "governance/seed_ledger.tsv",
        "governance/subject_registry.tsv",
        "private_truth/master_truth_manifest.json",
        "views/panel_view_manifest.json",
        "blinded_inputs/blinded_input_inventory.json",
        "governance/readiness_receipt.json",
    )
    for relative_path in required_roots:
        if not (release_root / relative_path).is_file():
            _fail("governance_root_missing", relative_path)
    release_manifest = load_json(release_root / required_roots[0])
    truth_manifest = load_json(release_root / required_roots[4])
    view_manifest = load_json(release_root / required_roots[5])
    inventory = load_json(release_root / required_roots[6])
    validate_document(release_manifest, schema_root / "panel_release_manifest.schema.json")
    validate_document(truth_manifest, schema_root / "master_truth_manifest.schema.json")
    validate_document(view_manifest, schema_root / "panel_view_manifest.schema.json")
    validate_document(inventory, schema_root / "blinded_input_inventory.schema.json")
    validate_blinded_inventory(inventory, view_root=release_root / "blinded_inputs", repository_root=repository_root, verify_files=True)
    separation_contract = load_json(release_root / "governance/separation_contract.json")
    validate_partition_separation(
        separation_contract["assignments"],
        separation_contract["reference_pools"],
    )
    derived_inventories = {
        name: _verify_derived_inventory(release_root / "derived_inputs" / name)
        for name in ("sensitivity", "render", "interpretation")
    }
    for artifact in truth_manifest["truth_artifacts"]:
        path = release_root / artifact["relative_path"]
        if not path.is_file() or file_sha256(path) != artifact["sha256"]:
            _fail("truth_artifact_integrity_failed", artifact["truth_id"])
    if len({artifact["truth_id"] for artifact in truth_manifest["truth_artifacts"]}) != 14:
        _fail("truth_family_count_mismatch")
    if len({view["view_type"] for view in view_manifest["views"]}) != 11:
        _fail("view_count_mismatch")

    bim = release_root / "blinded_inputs/RESOURCE_ONLY/shared/observed_autosomal.bim"
    fam = release_root / "blinded_inputs/RESOURCE_ONLY/shared/observed_autosomal.fam"
    chromosomes = {int(line.split("\t", 1)[0]) for line in bim.read_text(encoding="utf-8").splitlines() if line}
    subject_count = sum(1 for line in fam.read_text(encoding="utf-8").splitlines() if line)
    variant_count = sum(1 for line in bim.read_text(encoding="utf-8").splitlines() if line)
    if chromosomes != set(range(1, 23)):
        _fail("autosome_coverage_mismatch")
    if subject_count != 60:
        _fail("micro_subject_count_mismatch")

    run_manifests = sorted((release_root / "runs").glob("RUN-*/pipeline_run_manifest.json"))
    if len(run_manifests) != 4:
        _fail("closed_run_count_mismatch", str(len(run_manifests)))
    stage_to_run: dict[str, dict[str, Any]] = {}
    for path in run_manifests:
        manifest = load_json(path)
        validate_document(manifest, schema_root / "pipeline_run_manifest.schema.json")
        expected_hash = canonical_sha256({**manifest, "manifest_sha256": None})
        if manifest["manifest_sha256"] != expected_hash:
            _fail("run_manifest_hash_mismatch", manifest["run_id"])
        expected_inventory_hash = {
            CORE_RUN_ID: canonical_sha256(inventory),
            SENSITIVITY_RUN_ID: canonical_sha256(derived_inventories["sensitivity"]),
            RENDER_RUN_ID: canonical_sha256(derived_inventories["render"]),
            INTERPRET_RUN_ID: canonical_sha256(derived_inventories["interpretation"]),
        }[manifest["run_id"]]
        if manifest["consumed_inventory_sha256"] != expected_inventory_hash:
            _fail("run_inventory_hash_mismatch", manifest["run_id"])
        for stage_directory in sorted(
            child for child in path.parent.iterdir() if child.is_dir()
        ):
            _verify_stage_envelope(stage_directory)
        for output in manifest["stage_outputs"]:
            artifact_path = release_root / output["relative_path"]
            if not artifact_path.is_file() or file_sha256(artifact_path) != output["sha256"]:
                _fail("run_output_integrity_failed", output["stage_id"])
            stage_to_run[output["stage_id"]] = manifest
    if set(stage_to_run) != set(STAGE_IDS):
        _fail("stage_chain_incomplete")

    evaluations = sorted((release_root / "evaluations").glob("*/ETAPE_*.evaluation.json"))
    if len(evaluations) != 25:
        _fail("evaluation_count_mismatch", str(len(evaluations)))
    for path in evaluations:
        evaluation = load_json(path)
        validate_document(evaluation, schema_root / "stage_evaluation_manifest.schema.json")
        if evaluation["decision_scope"] != "TECHNICAL_INTEGRITY" or evaluation["scientific_performance_evaluated"] is not False:
            _fail("micro_scientific_decision_forbidden", evaluation["stage_id"])
        authorize_truth_access(
            evaluation_id=evaluation["evaluation_id"],
            stage_id=evaluation["stage_id"],
            requested_truth_ids=evaluation["truth_grants"],
            truth_manifest=truth_manifest,
            run_manifest=stage_to_run[evaluation["stage_id"]],
        )
    markdowns = sorted((release_root / "calibration/resultats").glob("ETAPE_*.md"))
    if len(markdowns) != 25:
        _fail("markdown_count_mismatch", str(len(markdowns)))
    forbidden_reports = list(release_root.rglob("*.html")) + list(release_root.rglob("*.pdf"))
    if forbidden_reports:
        _fail("forbidden_intermediate_report_created")
    for path in release_root.rglob("*"):
        if path.is_symlink():
            _fail("release_symlink_forbidden", str(path))
        if path.is_file() and path.stat(follow_symlinks=False).st_nlink != 1:
            _fail("release_hardlink_forbidden", str(path))
    if require_release_checksum_index:
        _verify_release_checksum_index(release_root)
    files = [path for path in release_root.rglob("*") if path.is_file()]
    logical_bytes = sum(path.stat().st_size for path in files)
    allocated_bytes = sum(path.stat().st_blocks * 512 for path in files)
    if allocated_bytes > 100 * 1024**3:
        _fail("micro_storage_budget_exceeded")
    reference_document = load_json(release_root / "governance/reference_pools.json")
    reference_count = sum(len(pool["member_ids"]) for pool in reference_document["pools"])
    return MicroReleaseAudit(
        release_id=RELEASE_ID,
        subject_count=subject_count,
        reference_count=reference_count,
        variant_count=variant_count,
        autosome_count=len(chromosomes),
        truth_family_count=len(truth_manifest["truth_artifacts"]),
        view_count=len(view_manifest["views"]),
        stage_count=len(stage_to_run),
        evaluation_count=len(evaluations),
        markdown_count=len(markdowns),
        run_count=len(run_manifests),
        sandbox_private_read_denied=True,
        external_scientific_tool_count=0,
        html_pdf_count=0,
        logical_bytes=logical_bytes,
        allocated_bytes=allocated_bytes,
    )


def _parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("action", choices=("build", "validate"))
    parser.add_argument("--repository-root", type=Path, default=Path.cwd())
    parser.add_argument("--release-root", type=Path)
    return parser


def main(argv: list[str] | None = None) -> int:
    arguments = _parser().parse_args(argv)
    repository_root = arguments.repository_root.resolve(strict=True)
    release_root = arguments.release_root or (
        repository_root / "calibration_workspace/releases" / RELEASE_DIRECTORY
    )
    if arguments.action == "build":
        audit = build_micro_release(repository_root, release_root)
    else:
        audit = validate_micro_release(repository_root, release_root)
    print(json.dumps(asdict(audit), ensure_ascii=False, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
