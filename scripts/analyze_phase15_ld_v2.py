"""Exécute la comparaison LD phase 15 v2 dans un dossier derived isolé."""

from __future__ import annotations

import argparse
import csv
import hashlib
import html
import json
import math
import random
import shutil
import statistics
import subprocess
import tempfile
from collections import defaultdict
from pathlib import Path
from typing import Any, Iterable, Mapping, Sequence

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from effet_fondateur.contracts import DocumentValidationError, validate_json_document
from effet_fondateur.ld.matched import (
    MatchedLdError,
    PairMetric,
    PhasedVariant,
    bin_labels,
    centered_pair_metrics,
    count_neighbors_by_bin,
    density_is_matched,
    local_background_metrics,
    median_by_bin,
    parse_phased_genotype,
    r_squared,
)


METHOD_ID = "phased_target_vs_matched_genomewide_ld_v2"
TARGET_VARIANT_ID = "DOCK6_GRCh38_19_11237780_C_A"
DISTANCE_BINS_CM = (0.1, 0.25, 0.5, 1.0)
MINIMUM_NEIGHBOR_MAF = 0.05
ANCHOR_MINOR_COPY_TOLERANCE = 2
DENSITY_MATCH_RATIO = 2.0
MAXIMUM_MATCHED_ANCHORS = 400
MINIMUM_MATCHED_ANCHORS = 100
MINIMUM_CALLED_HAPLOTYPES = 100
RANDOM_SEED = 20260817


class Phase15LdV2Error(RuntimeError):
    """Signale une entrée de run absente ou une sortie non interprétable."""


def _parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(prog="analyze-phase15-ld-v2")
    parser.add_argument("--run-dir", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path)
    parser.add_argument("--bcftools", default="bcftools")
    return parser


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _read_json_object(path: Path, label: str) -> dict[str, Any]:
    """Lit un document JSON sans exposer son contenu dans les erreurs."""

    try:
        document = json.loads(path.read_text(encoding="utf-8"))
    except FileNotFoundError as error:
        raise Phase15LdV2Error(f"{label}_missing:{path}") from error
    except json.JSONDecodeError as error:
        raise Phase15LdV2Error(f"{label}_invalid_json:{path}") from error
    if not isinstance(document, dict):
        raise Phase15LdV2Error(f"{label}_not_object:{path}")
    return document


def _load_verified_stage_artifacts(
    run_dir: Path,
    manifest: Mapping[str, Any],
    stage_name: str,
    stage_outputs_path: Path,
) -> dict[str, dict[str, Any]]:
    """Valide l'identité d'une étape publiée et indexe ses artefacts signés."""

    records = [
        record for record in manifest.get("stages", [])
        if record.get("stage_name") == stage_name
    ]
    if len(records) != 1:
        raise Phase15LdV2Error(f"stage_record_missing_or_ambiguous:{stage_name}")
    record = records[0]
    if record.get("state") not in {"SUCCEEDED", "CACHED"}:
        raise Phase15LdV2Error(f"source_stage_not_successful:{stage_name}")
    audit_path = record.get("audit_path")
    expected_outputs_sha256 = record.get("stage_outputs_sha256")
    signature = record.get("signature")
    if not audit_path or not expected_outputs_sha256 or not signature:
        raise Phase15LdV2Error(f"source_stage_provenance_incomplete:{stage_name}")
    relative_outputs_path = Path(audit_path).parent / "stage_outputs.json"
    if relative_outputs_path.is_absolute() or ".." in relative_outputs_path.parts:
        raise Phase15LdV2Error(f"source_stage_path_unsafe:{stage_name}")
    expected_outputs_path = (run_dir / relative_outputs_path).resolve()
    if not expected_outputs_path.is_relative_to(run_dir):
        raise Phase15LdV2Error(f"source_stage_path_unsafe:{stage_name}")
    if stage_outputs_path.resolve() != expected_outputs_path:
        raise Phase15LdV2Error(f"source_stage_outputs_path_mismatch:{stage_name}")
    if stage_outputs_path.is_symlink() or not stage_outputs_path.is_file():
        raise Phase15LdV2Error(f"source_stage_outputs_missing:{stage_name}")
    if _sha256(stage_outputs_path) != expected_outputs_sha256:
        raise Phase15LdV2Error(f"source_stage_outputs_checksum_mismatch:{stage_name}")

    outputs = _read_json_object(stage_outputs_path, "stage_outputs")
    validate_json_document(outputs, "stage_outputs.schema.json")
    expected_identity = {
        "run_id": manifest.get("run_id"),
        "stage_id": record.get("stage_id"),
        "stage_name": stage_name,
        "signature": signature,
    }
    if any(outputs.get(field) != expected for field, expected in expected_identity.items()):
        raise Phase15LdV2Error(f"source_stage_identity_mismatch:{stage_name}")

    artifacts: dict[str, dict[str, Any]] = {}
    expected_producer_stage = relative_outputs_path.parent.name
    for artifact in outputs["artifacts"]:
        artifact_id = artifact["artifact_id"]
        if artifact_id in artifacts:
            raise Phase15LdV2Error(f"source_artifact_duplicate:{artifact_id}")
        if (
            artifact["producer_stage"] != expected_producer_stage
            or artifact["producer_signature"] != signature
        ):
            raise Phase15LdV2Error(f"source_artifact_producer_mismatch:{artifact_id}")
        artifacts[artifact_id] = artifact
    return artifacts


def _validated_artifact_path(
    run_dir: Path,
    artifacts: Mapping[str, dict[str, Any]],
    artifact_id: str,
    expected_relative_path: str,
    expected_artifact_type: str,
) -> Path:
    """Résout un artefact attendu et vérifie chemin, type et empreinte."""

    artifact = artifacts.get(artifact_id)
    if artifact is None:
        raise Phase15LdV2Error(f"source_artifact_missing:{artifact_id}")
    if (
        artifact["path"] != expected_relative_path
        or artifact["artifact_type"] != expected_artifact_type
    ):
        raise Phase15LdV2Error(f"source_artifact_contract_mismatch:{artifact_id}")
    relative_path = Path(artifact["path"])
    if relative_path.is_absolute() or ".." in relative_path.parts:
        raise Phase15LdV2Error(f"source_artifact_path_unsafe:{artifact_id}")
    unresolved_path = run_dir / relative_path
    path = unresolved_path.resolve()
    if not path.is_relative_to(run_dir):
        raise Phase15LdV2Error(f"source_artifact_path_unsafe:{artifact_id}")
    if unresolved_path.is_symlink() or not path.is_file():
        raise Phase15LdV2Error(f"source_artifact_file_missing:{artifact_id}")
    if _sha256(path) != artifact["sha256"]:
        raise Phase15LdV2Error(f"source_artifact_checksum_mismatch:{artifact_id}")
    return path


def _read_keep(path: Path) -> list[str]:
    samples: list[str] = []
    for line_number, line in enumerate(path.read_text(encoding="utf-8").splitlines(), start=1):
        fields = line.split()
        if len(fields) != 2:
            raise Phase15LdV2Error(f"malformed_keep:{path.name}:{line_number}")
        samples.append(fields[1])
    if len(samples) != len(set(samples)):
        raise Phase15LdV2Error(f"duplicate_sample_in_keep:{path.name}")
    return samples


def _bcf_samples(bcftools: str, path: Path) -> list[str]:
    completed = subprocess.run(
        [bcftools, "query", "-l", str(path)], capture_output=True, check=False,
        text=True, timeout=60,
    )
    if completed.returncode != 0:
        raise Phase15LdV2Error(f"bcftools_sample_query_failed:{path}")
    samples = [line for line in completed.stdout.splitlines() if line]
    if not samples or len(samples) != len(set(samples)):
        raise Phase15LdV2Error(f"invalid_bcf_samples:{path}")
    return samples


def _load_map(path: Path) -> dict[str, tuple[str, int, float]]:
    result: dict[str, tuple[str, int, float]] = {}
    for line_number, line in enumerate(path.read_text(encoding="utf-8").splitlines(), start=1):
        fields = line.split()
        if len(fields) < 4:
            raise Phase15LdV2Error(f"malformed_map:{path.name}:{line_number}")
        chromosome = fields[0].removeprefix("chr")
        variant_id = fields[1]
        try:
            position_cm = float(fields[2])
            position_bp = int(fields[3])
        except ValueError as error:
            raise Phase15LdV2Error(f"malformed_map:{path.name}:{line_number}") from error
        if variant_id in result:
            raise Phase15LdV2Error(f"duplicate_map_variant:{variant_id}")
        result[variant_id] = (chromosome, position_bp, position_cm)
    return result


def _query_phased_variants(
    bcftools: str,
    bcf_path: Path,
    map_rows: dict[str, tuple[str, int, float]],
    *,
    sample_names: Sequence[str] | None = None,
) -> list[PhasedVariant]:
    command = [bcftools, "query"]
    if sample_names is not None:
        command.extend(["-s", ",".join(sample_names)])
    command.extend(["-f", "%CHROM\\t%POS\\t%ID[\\t%GT]\\n", str(bcf_path)])
    process = subprocess.Popen(
        command, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True,
    )
    assert process.stdout is not None
    variants: list[PhasedVariant] = []
    expected_haplotype_count = None if sample_names is None else 2 * len(sample_names)
    try:
        for line_number, line in enumerate(process.stdout, start=1):
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 4:
                raise Phase15LdV2Error(f"malformed_bcf_query:{bcf_path.name}:{line_number}")
            chromosome = fields[0].removeprefix("chr")
            variant_id = fields[2]
            if variant_id not in map_rows:
                raise Phase15LdV2Error(f"variant_missing_from_genetic_map:{variant_id}")
            map_chromosome, map_bp, position_cm = map_rows[variant_id]
            if chromosome != map_chromosome or int(fields[1]) != map_bp:
                raise Phase15LdV2Error(f"variant_map_coordinate_mismatch:{variant_id}")
            haplotypes = tuple(
                allele
                for genotype in fields[3:]
                for allele in parse_phased_genotype(genotype)
            )
            if expected_haplotype_count is not None and len(haplotypes) != expected_haplotype_count:
                raise Phase15LdV2Error(f"unexpected_haplotype_count:{variant_id}")
            variants.append(PhasedVariant(
                chromosome=chromosome,
                variant_id=variant_id,
                position_bp=map_bp,
                position_cm=position_cm,
                haplotypes=haplotypes,
            ))
    except Exception:
        process.kill()
        process.wait()
        raise
    stderr = process.stderr.read() if process.stderr is not None else ""
    return_code = process.wait()
    if return_code != 0:
        raise Phase15LdV2Error(
            f"bcftools_variant_query_failed:{bcf_path.name}:{stderr.strip()[:200]}"
        )
    variants.sort(key=lambda row: (row.position_cm, row.position_bp, row.variant_id))
    return variants


def _subset_variants(
    variants: Sequence[PhasedVariant], sample_indices: Sequence[int],
) -> list[PhasedVariant]:
    haplotype_indices = [index for sample_index in sample_indices for index in (2 * sample_index, 2 * sample_index + 1)]
    return [
        PhasedVariant(
            chromosome=variant.chromosome,
            variant_id=variant.variant_id,
            position_bp=variant.position_bp,
            position_cm=variant.position_cm,
            haplotypes=tuple(variant.haplotypes[index] for index in haplotype_indices),
        )
        for variant in variants
    ]


def _find_variant(variants: Sequence[PhasedVariant], variant_id: str) -> PhasedVariant:
    matches = [variant for variant in variants if variant.variant_id == variant_id]
    if len(matches) != 1:
        raise Phase15LdV2Error(f"target_variant_missing_or_ambiguous:{variant_id}")
    return matches[0]


def _select_matched_anchors(
    chromosome_variants: dict[str, list[PhasedVariant]],
    *,
    target_minor_copy_count: int,
    expected_called_haplotypes: int,
    target_neighbor_counts: dict[str, int],
    seed: int,
    maximum_anchors: int,
) -> list[tuple[PhasedVariant, list[PhasedVariant], dict[str, int]]]:
    randomizer = random.Random(seed)
    eligible_by_chromosome: dict[str, list[tuple[PhasedVariant, list[PhasedVariant], dict[str, int]]]] = {}
    for chromosome, variants in sorted(chromosome_variants.items(), key=lambda item: int(item[0])):
        candidates = []
        for variant in variants:
            if variant.called_haplotype_count != expected_called_haplotypes:
                continue
            maf = variant.minor_allele_frequency
            if maf is None or maf <= 0:
                continue
            minor_copies = round(maf * variant.called_haplotype_count)
            if abs(minor_copies - target_minor_copy_count) > ANCHOR_MINOR_COPY_TOLERANCE:
                continue
            candidates.append(variant)
        randomizer.shuffle(candidates)
        matched = []
        # Le plafond par chromosome limite le coût tout en laissant une large
        # réserve avant l'échantillonnage équilibré final.
        for candidate in candidates[:250]:
            counts = count_neighbors_by_bin(
                candidate, variants, distance_bins_cm=DISTANCE_BINS_CM,
                minimum_neighbor_maf=MINIMUM_NEIGHBOR_MAF,
                minimum_called_haplotypes=MINIMUM_CALLED_HAPLOTYPES,
            )
            if density_is_matched(counts, target_neighbor_counts, ratio=DENSITY_MATCH_RATIO):
                matched.append((candidate, variants, counts))
        eligible_by_chromosome[chromosome] = matched

    selected = []
    chromosome_order = sorted(eligible_by_chromosome, key=int)
    position = 0
    while len(selected) < maximum_anchors:
        added = False
        for chromosome in chromosome_order:
            candidates = eligible_by_chromosome[chromosome]
            if position < len(candidates):
                selected.append(candidates[position])
                added = True
                if len(selected) == maximum_anchors:
                    break
        if not added:
            break
        position += 1
    return selected


def _write_tsv(path: Path, columns: Sequence[str], rows: Iterable[dict[str, object]]) -> None:
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for row in rows:
            writer.writerow({column: row.get(column, "") for column in columns})


def _number(value: float | None) -> str:
    return "" if value is None else format(value, ".10g")


def _quantile(values: Sequence[float], probability: float) -> float | None:
    return None if not values else float(np.quantile(np.asarray(values), probability))


def _plot(
    output_dir: Path,
    summary_rows: Sequence[dict[str, object]],
    target_pair_rows: Sequence[dict[str, object]],
) -> list[Path]:
    labels = [str(row["DISTANCE_BIN"]) for row in summary_rows]
    x = np.arange(len(labels))
    def plot_value(value: object) -> float:
        return float(value) if value not in {"", None} else math.nan

    q25 = np.asarray([plot_value(row["MATCHED_Q25_R2"]) for row in summary_rows])
    q75 = np.asarray([plot_value(row["MATCHED_Q75_R2"]) for row in summary_rows])
    p95 = np.asarray([plot_value(row["MATCHED_P95_R2"]) for row in summary_rows])
    target = np.asarray([plot_value(row["TARGET_MEDIAN_R2_INDEPENDENT"]) for row in summary_rows])
    family = np.asarray([plot_value(row["TARGET_MEDIAN_R2_ALL75"]) for row in summary_rows])
    local = np.asarray([plot_value(row["LOCAL_BACKGROUND_MEDIAN_R2"]) for row in summary_rows])
    visible_values = np.concatenate((q75, p95, target, family, local))
    finite_values = visible_values[np.isfinite(visible_values)]
    display_maximum = min(1.0, max(0.12, math.ceil(float(finite_values.max()) * 20) / 20 + 0.02))

    figure, axes = plt.subplots(1, 2, figsize=(15, 6.2), gridspec_kw={"width_ratios": [1.15, 1]})
    first = axes[0]
    first.fill_between(x, q25, q75, color="#cbd5e1", alpha=0.8, label="Fond génomique apparié — IQR")
    first.plot(x, p95, color="#64748b", linestyle=":", marker="o", label="95e percentile génomique")
    first.plot(x, local, color="#2563eb", marker="o", linewidth=2.2, label="Fond local DOCK6 (65 indépendants)")
    first.plot(x, target, color="#ea580c", marker="o", linewidth=2.8, label="DOCK6 centré — 62+3 indépendants")
    first.plot(x, family, color="#dc2626", marker="x", linestyle="--", linewidth=1.6,
               label="DOCK6 — 75 sujets (sensibilité familiale)")
    first.set_xticks(x, labels)
    first.set_ylim(-0.005, display_maximum)
    first.set_ylabel("r² haplotypique médian")
    first.set_xlabel("Distance génétique à l’ancre")
    first.set_title("LD autour de DOCK6 comparée au fond génomique")
    first.grid(axis="y", alpha=0.18)
    first.legend(fontsize=8.5, loc="upper right")
    for index, row in enumerate(summary_rows):
        if int(row["TARGET_PAIR_COUNT"]) == 0:
            first.annotate(
                "aucune paire cible", (index, display_maximum * 0.16),
                ha="center", va="center", rotation=90, fontsize=8, color="#7c2d12",
            )

    second = axes[1]
    distances = np.asarray([float(row["SIGNED_DISTANCE_CM"]) for row in target_pair_rows])
    independent = np.asarray([float(row["R2_INDEPENDENT_65"]) for row in target_pair_rows])
    all75 = np.asarray([float(row["R2_ALL_75"]) for row in target_pair_rows])
    second.scatter(distances, independent, s=30, alpha=0.8, color="#ea580c", label="65 indépendants")
    second.scatter(distances, all75, s=18, alpha=0.42, color="#dc2626", marker="x", label="75 sujets")
    second.axvline(0, color="#0f172a", linewidth=1.2)
    second.set_xlim(-DISTANCE_BINS_CM[-1], DISTANCE_BINS_CM[-1])
    scatter_maximum = max(float(independent.max()), float(all75.max()))
    second.set_ylim(-0.005, min(1.0, max(0.12, math.ceil(scatter_maximum * 20) / 20 + 0.02)))
    second.set_xlabel("Distance signée à la variation DOCK6 (cM)")
    second.set_ylabel("r² haplotypique")
    second.set_title("Marqueurs voisins centrés sur la variation")
    second.grid(alpha=0.18)
    second.legend(fontsize=9)
    figure.suptitle("Phase 15 v2 — signal LD local et référence génomique appariée", fontsize=16, fontweight="bold")
    figure.text(
        0.5, 0.01,
        "Résultat exploratoire : trois familles indépendantes seulement. La courbe des 75 sujets est une sensibilité confondue par la parenté.",
        ha="center", fontsize=9, color="#7c2d12",
    )
    figure.tight_layout(rect=(0, 0.045, 1, 0.95))
    paths = [output_dir / "phase15_ld_v2_comparison.png", output_dir / "phase15_ld_v2_comparison.svg"]
    figure.savefig(paths[0], dpi=190, bbox_inches="tight", facecolor="white")
    figure.savefig(paths[1], bbox_inches="tight", facecolor="white")
    plt.close(figure)
    return paths


def _markdown(
    summary: dict[str, object], summary_rows: Sequence[dict[str, object]],
) -> str:
    table = "\n".join(
        "| {DISTANCE_BIN} | {TARGET_MEDIAN_R2_INDEPENDENT} | {MATCHED_MEDIAN_R2} | "
        "{MATCHED_P95_R2} | {EMPIRICAL_PERCENTILE} | {BIN_INTERPRETATION} |".format(**row)
        for row in summary_rows
    )
    if summary["evidence_direction"] == "NEUTRAL_EXPLORATORY_SMALL_N":
        run_interpretation = (
            "Dans ce run, aucune classe évaluée ne dépasse le 95e percentile "
            "apparié. Le LD n'ajoute donc pas de soutien autonome à l'hypothèse "
            "fondatrice, sans constituer une preuve contre celle-ci compte tenu "
            "des trois familles seulement."
        )
    else:
        run_interpretation = (
            "Dans ce run, l'excès préspécifié est atteint. Il constitue un soutien "
            "exploratoire qui doit rester secondaire aux résultats IBS et IBD."
        )
    return f"""# Corrections de la phase 15 — LD v2

## But de la correction

La phase 15 historique décrivait le LD local chez les 62 témoins mais ne pouvait
pas calculer un LD centré sur la variation DOCK6, monomorphe dans ce groupe. La
version `{METHOD_ID}` ajoute une comparaison exploratoire sans modifier ni
écraser le run publié.

## Corrections appliquées

- analyse principale sur **65 unités indépendantes** : 62 témoins non apparentés
  et un représentant explicite de chacune des trois familles porteuses ;
- calcul de r² sur les haplotypes phasés, sans demander à PLINK de traiter les
  apparentés comme des fondateurs ;
- voisinage DOCK6 limité à 1 cM et marqueurs voisins de MAF ≥ 5 % ;
- ancrages génomiques appariés sur le nombre de copies de l'allèle mineur
  (tolérance ±{ANCHOR_MINOR_COPY_TOLERANCE}) et la densité de marqueurs dans
  chaque classe de distance (ratio maximal {DENSITY_MATCH_RATIO:g}) ;
- chromosomes génomiques de référence : 1–22 sauf le chromosome 19, pour éviter
  de réutiliser la région cible ;
- graine fixe `{RANDOM_SEED}` et maximum de {MAXIMUM_MATCHED_ANCHORS} fenêtres ;
- une classe sans paire cible reste `NOT_EVALUATED` et ne contraint pas la
  densité des fenêtres de référence ; elle n'est ni fusionnée ni interprétée ;
- analyse des 75 sujets conservée seulement comme **sensibilité confondue par la
  parenté**, jamais comme preuve principale ;
- aucun calcul d'âge et aucune probabilité globale d'effet fondateur dérivés de
  ce résultat.

## Résultat du run

- variation cible : `{TARGET_VARIANT_ID}` ;
- MAF chez les 65 indépendants : {summary['target_maf_independent']:.4f}
  ({summary['target_minor_copy_count_independent']} copies mineures sur
  {summary['independent_haplotype_count']} haplotypes) ;
- fenêtres génomiques appariées retenues : {summary['matched_anchor_count']} ;
- direction qualitative : **{summary['evidence_direction']}** ;
- classes dépassant le 95e percentile génomique :
  {summary['bin_count_above_matched_p95']} sur {len(summary_rows)}.

| Distance | DOCK6 r² médian | Fond médian | Fond P95 | Percentile empirique | Lecture |
|---|---:|---:|---:|---:|---|
{table}

## Interprétation autorisée

{run_interpretation}

Un excès de LD autour de DOCK6 apporte un **soutien secondaire** à la présence
d'un haplotype transmis depuis une origine commune. Il ne démontre pas seul un
effet fondateur. Le nombre d'unités porteuses indépendantes reste égal à trois ;
le résultat doit donc être confronté aux segments IBS, aux appels IBD explicites,
aux ROH et à l'estimation exploratoire de l'âge.

La différence entre la courbe indépendante et celle des 75 sujets quantifie le
risque d'inflation lié aux transmissions familiales. La courbe à 75 sujets ne
doit pas être utilisée pour classer la preuve.

## Fichiers produits

- `phase15_ld_v2_comparison.svg` et `.png` : figure principale ;
- `target_centered_ld_pairs.tsv` : LD de chaque voisin avec la variation ;
- `matched_genomewide_windows.tsv` : résumés par ancre génomique ;
- `ld_v2_summary.tsv` : comparaison par classe de distance ;
- `analysis_summary.json` : paramètres, provenance et conclusion limitée ;
- `checksums.sha256` : empreintes des sorties.
"""


def execute(run_dir: Path, output_dir: Path | None, bcftools: str) -> Path:
    run_dir = run_dir.resolve()
    if not run_dir.is_dir():
        raise Phase15LdV2Error(f"run_dir_missing:{run_dir}")
    if shutil.which(bcftools) is None:
        raise Phase15LdV2Error("bcftools_not_available")
    destination = (output_dir or run_dir / "derived/phase15_ld_v2").resolve()
    if destination.exists():
        raise Phase15LdV2Error(f"output_dir_already_exists:{destination}")

    manifest_path = run_dir / "manifest.json"
    manifest = _read_json_object(manifest_path, "run_manifest")
    validate_json_document(manifest, "run_manifest.schema.json")
    if manifest["run_id"] != run_dir.name:
        raise Phase15LdV2Error("run_manifest_identity_mismatch")

    stage09_outputs = run_dir / "stages/09_freeze_cohorts/stage_outputs.json"
    stage11_outputs = run_dir / "stages/11_prepare_target_region/stage_outputs.json"
    stage12_outputs = run_dir / "stages/12_phase_target_region/stage_outputs.json"
    stage12a_outputs = run_dir / "stages/12A_phase_autosomal_panel/stage_outputs.json"
    stage09_artifacts = _load_verified_stage_artifacts(
        run_dir, manifest, "freeze_cohorts", stage09_outputs,
    )
    stage11_artifacts = _load_verified_stage_artifacts(
        run_dir, manifest, "prepare_target_region", stage11_outputs,
    )
    stage12_artifacts = _load_verified_stage_artifacts(
        run_dir, manifest, "phase_target_region", stage12_outputs,
    )
    stage12a_artifacts = _load_verified_stage_artifacts(
        run_dir, manifest, "phase_autosomal_panel", stage12a_outputs,
    )

    controls_keep = _validated_artifact_path(
        run_dir, stage09_artifacts, "cohort_keep_controls_unrelated",
        "stages/09_freeze_cohorts/keep/controls_unrelated.keep", "plink_keep",
    )
    carriers_keep = _validated_artifact_path(
        run_dir, stage09_artifacts, "cohort_keep_target_carriers_independent",
        "stages/09_freeze_cohorts/keep/target_carriers_independent.keep",
        "plink_keep",
    )
    target_map_path = _validated_artifact_path(
        run_dir, stage11_artifacts, "target_region_bim",
        "stages/11_prepare_target_region/target_region.bim", "plink_bim",
    )
    target_bcf = _validated_artifact_path(
        run_dir, stage12_artifacts, "shapeit5_final_bcf",
        "stages/12_phase_target_region/shapeit5_phasing/target.phased.bcf",
        "shapeit5_final_bcf",
    )
    autosomal_inputs: dict[int, tuple[Path, Path]] = {}
    autosomal_source_paths: list[Path] = []
    for chromosome in range(1, 23):
        if chromosome == 19:
            continue
        bcf_path = _validated_artifact_path(
            run_dir,
            stage12a_artifacts,
            f"autosomal_phased_chr{chromosome}_bcf",
            f"stages/12A_phase_autosomal_panel/chr{chromosome}/study.phased.bcf",
            "autosomal_phased_bcf",
        )
        map_path = _validated_artifact_path(
            run_dir,
            stage12a_artifacts,
            f"autosomal_phased_chr{chromosome}_ibd_map",
            f"stages/12A_phase_autosomal_panel/chr{chromosome}/study.phased.plink.map",
            "autosomal_phased_ibd_map",
        )
        autosomal_inputs[chromosome] = (bcf_path, map_path)
        autosomal_source_paths.extend((bcf_path, map_path))
    source_paths = [
        manifest_path,
        stage09_outputs,
        stage11_outputs,
        stage12_outputs,
        stage12a_outputs,
        controls_keep,
        carriers_keep,
        target_bcf,
        target_map_path,
        *autosomal_source_paths,
    ]

    controls = _read_keep(controls_keep)
    carriers = _read_keep(carriers_keep)
    if len(controls) != 62 or len(carriers) != 3 or set(controls) & set(carriers):
        raise Phase15LdV2Error("unexpected_independent_cohort_composition")
    independent_set = set(controls) | set(carriers)

    all_samples = _bcf_samples(bcftools, target_bcf)
    if len(all_samples) != 75 or not independent_set <= set(all_samples):
        raise Phase15LdV2Error("target_bcf_sample_universe_mismatch")
    independent_samples = [sample for sample in all_samples if sample in independent_set]
    independent_indices = [index for index, sample in enumerate(all_samples) if sample in independent_set]
    if len(independent_samples) != 65:
        raise Phase15LdV2Error("independent_sample_count_mismatch")

    target_map = _load_map(target_map_path)
    target_variants_all = _query_phased_variants(bcftools, target_bcf, target_map)
    target_variants_independent = _subset_variants(target_variants_all, independent_indices)
    target_all = _find_variant(target_variants_all, TARGET_VARIANT_ID)
    target_independent = _find_variant(target_variants_independent, TARGET_VARIANT_ID)
    target_maf = target_independent.minor_allele_frequency
    if target_maf is None or target_maf <= 0:
        raise Phase15LdV2Error("target_not_polymorphic_in_independent_set")
    target_minor_copies = round(target_maf * target_independent.called_haplotype_count)

    target_pairs = centered_pair_metrics(
        target_independent, target_variants_independent,
        distance_bins_cm=DISTANCE_BINS_CM,
        minimum_neighbor_maf=MINIMUM_NEIGHBOR_MAF,
        minimum_called_haplotypes=MINIMUM_CALLED_HAPLOTYPES,
    )
    if not target_pairs:
        raise Phase15LdV2Error("no_target_centered_pairs")
    target_neighbor_counts = count_neighbors_by_bin(
        target_independent,
        target_variants_independent,
        distance_bins_cm=DISTANCE_BINS_CM,
        minimum_neighbor_maf=MINIMUM_NEIGHBOR_MAF,
        minimum_called_haplotypes=MINIMUM_CALLED_HAPLOTYPES,
    )
    evaluated_target_counts = {
        label: sum(metric.distance_bin == label for metric in target_pairs)
        for label in bin_labels(DISTANCE_BINS_CM)
    }
    if target_neighbor_counts != evaluated_target_counts:
        raise Phase15LdV2Error("target_density_evaluation_mismatch")
    local_background = local_background_metrics(
        target_variants_independent, excluded_variant_id=TARGET_VARIANT_ID,
        distance_bins_cm=DISTANCE_BINS_CM,
        minimum_maf=MINIMUM_NEIGHBOR_MAF,
        minimum_called_haplotypes=MINIMUM_CALLED_HAPLOTYPES,
    )

    chromosome_variants: dict[str, list[PhasedVariant]] = {}
    for chromosome, (bcf_path, map_path) in autosomal_inputs.items():
        if _bcf_samples(bcftools, bcf_path) != all_samples:
            raise Phase15LdV2Error(
                f"autosomal_bcf_sample_order_mismatch:chr{chromosome}"
            )
        chromosome_variants[str(chromosome)] = _query_phased_variants(
            bcftools, bcf_path, _load_map(map_path), sample_names=independent_samples,
        )

    selected_anchors = _select_matched_anchors(
        chromosome_variants,
        target_minor_copy_count=target_minor_copies,
        expected_called_haplotypes=target_independent.called_haplotype_count,
        target_neighbor_counts=target_neighbor_counts,
        seed=RANDOM_SEED,
        maximum_anchors=MAXIMUM_MATCHED_ANCHORS,
    )
    if len(selected_anchors) < MINIMUM_MATCHED_ANCHORS:
        raise Phase15LdV2Error(f"insufficient_matched_anchors:{len(selected_anchors)}")

    anchor_rows: list[dict[str, object]] = []
    anchor_bin_values: dict[str, list[float]] = defaultdict(list)
    for anchor, variants, density_counts in selected_anchors:
        metrics = centered_pair_metrics(
            anchor, variants, distance_bins_cm=DISTANCE_BINS_CM,
            minimum_neighbor_maf=MINIMUM_NEIGHBOR_MAF,
            minimum_called_haplotypes=MINIMUM_CALLED_HAPLOTYPES,
        )
        medians = median_by_bin(metrics, DISTANCE_BINS_CM)
        for label in bin_labels(DISTANCE_BINS_CM):
            value = medians[label]
            if value is None:
                continue
            anchor_bin_values[label].append(value)
            anchor_rows.append({
                "ANCHOR_ID": anchor.variant_id,
                "CHROMOSOME": anchor.chromosome,
                "POSITION_BP": anchor.position_bp,
                "POSITION_CM": _number(anchor.position_cm),
                "ANCHOR_MAF": _number(anchor.minor_allele_frequency),
                "MINOR_COPY_COUNT": round((anchor.minor_allele_frequency or 0) * anchor.called_haplotype_count),
                "DISTANCE_BIN": label,
                "ELIGIBLE_NEIGHBOR_COUNT": density_counts[label],
                "EVALUATED_PAIR_COUNT": sum(metric.distance_bin == label for metric in metrics),
                "MEDIAN_R2": _number(value),
            })

    full_by_id = {variant.variant_id: variant for variant in target_variants_all}
    target_pair_rows: list[dict[str, object]] = []
    all75_by_bin: dict[str, list[float]] = defaultdict(list)
    for metric in target_pairs:
        neighbor_all = full_by_id[metric.neighbor_id]
        full_r2, full_called = r_squared(
            target_all.haplotypes, neighbor_all.haplotypes,
            minimum_called=MINIMUM_CALLED_HAPLOTYPES,
        )
        if full_r2 is None:
            raise Phase15LdV2Error(f"all75_pair_not_evaluable:{metric.neighbor_id}")
        all75_by_bin[metric.distance_bin].append(full_r2)
        target_pair_rows.append({
            "CHROMOSOME": metric.chromosome,
            "TARGET_VARIANT_ID": metric.anchor_id,
            "NEIGHBOR_ID": metric.neighbor_id,
            "NEIGHBOR_POSITION_BP": metric.neighbor_position_bp,
            "SIGNED_DISTANCE_CM": _number(metric.signed_distance_cm),
            "DISTANCE_BIN": metric.distance_bin,
            "NEIGHBOR_MAF_INDEPENDENT": _number(metric.neighbor_maf),
            "CALLED_HAPLOTYPES_INDEPENDENT": metric.called_haplotype_count,
            "R2_INDEPENDENT_65": _number(metric.r2),
            "CALLED_HAPLOTYPES_ALL75": full_called,
            "R2_ALL_75": _number(full_r2),
        })

    target_medians = median_by_bin(target_pairs, DISTANCE_BINS_CM)
    summary_rows: list[dict[str, object]] = []
    above_p95_count = 0
    for label in bin_labels(DISTANCE_BINS_CM):
        background = anchor_bin_values[label]
        target_value = target_medians[label]
        full_value = statistics.median(all75_by_bin[label]) if all75_by_bin[label] else None
        local_value = statistics.median(local_background[label]) if local_background[label] else None
        p95 = _quantile(background, 0.95)
        if target_value is None or p95 is None:
            interpretation = "NOT_EVALUATED"
            percentile = None
        else:
            percentile = (sum(value <= target_value for value in background) + 1) / (len(background) + 1)
            interpretation = "ABOVE_MATCHED_P95" if target_value > p95 else "WITHIN_MATCHED_BACKGROUND"
            above_p95_count += interpretation == "ABOVE_MATCHED_P95"
        summary_rows.append({
            "DISTANCE_BIN": label,
            "TARGET_PAIR_COUNT": target_neighbor_counts[label],
            "TARGET_MEDIAN_R2_INDEPENDENT": _number(target_value),
            "TARGET_MEDIAN_R2_ALL75": _number(full_value),
            "LOCAL_BACKGROUND_PAIR_COUNT": len(local_background[label]),
            "LOCAL_BACKGROUND_MEDIAN_R2": _number(local_value),
            "MATCHED_ANCHOR_COUNT": len(background),
            "MATCHED_MEDIAN_R2": _number(_quantile(background, 0.5)),
            "MATCHED_Q25_R2": _number(_quantile(background, 0.25)),
            "MATCHED_Q75_R2": _number(_quantile(background, 0.75)),
            "MATCHED_P95_R2": _number(p95),
            "EMPIRICAL_PERCENTILE": _number(percentile),
            "BIN_INTERPRETATION": interpretation,
        })

    evidence_direction = (
        "SUPPORTIVE_EXPLORATORY_SMALL_N"
        if above_p95_count >= 2
        else "NEUTRAL_EXPLORATORY_SMALL_N"
    )
    summary: dict[str, object] = {
        "schema_version": "1.0.0",
        "method_id": METHOD_ID,
        "run_id": run_dir.name,
        "analysis_role": "SECONDARY_SUPPORTIVE_ONLY",
        "feeds_variant_age": False,
        "primary_sample_set": "62_controls_plus_3_independent_family_representatives",
        "control_count": len(controls),
        "independent_family_count": len(carriers),
        "independent_sample_count": len(independent_samples),
        "independent_haplotype_count": target_independent.called_haplotype_count,
        "family_confounded_sensitivity_sample_count": len(all_samples),
        "target_variant_id": TARGET_VARIANT_ID,
        "target_maf_independent": target_maf,
        "target_minor_copy_count_independent": target_minor_copies,
        "matched_anchor_count": len(selected_anchors),
        "excluded_matched_background_chromosomes": ["19"],
        "distance_bins_cm": list(DISTANCE_BINS_CM),
        "minimum_neighbor_maf": MINIMUM_NEIGHBOR_MAF,
        "anchor_minor_copy_tolerance": ANCHOR_MINOR_COPY_TOLERANCE,
        "density_match_ratio": DENSITY_MATCH_RATIO,
        "random_seed": RANDOM_SEED,
        "verified_autosomal_artifact_count": len(autosomal_source_paths),
        "provenance_verification": "RUN_MANIFEST_STAGE_OUTPUTS_AND_ARTIFACT_SHA256",
        "bin_count_above_matched_p95": above_p95_count,
        "evidence_direction": evidence_direction,
        "interpretation_limit": (
            "Trois familles indépendantes seulement; soutien LD secondaire, "
            "sans probabilité d'effet fondateur ni contribution automatique à l'âge."
        ),
        "source_sha256": {str(path.relative_to(run_dir)): _sha256(path) for path in source_paths},
    }

    destination.parent.mkdir(parents=True, exist_ok=True)
    temporary = Path(tempfile.mkdtemp(prefix=f".{destination.name}.", dir=destination.parent))
    try:
        _write_tsv(
            temporary / "target_centered_ld_pairs.tsv",
            ("CHROMOSOME", "TARGET_VARIANT_ID", "NEIGHBOR_ID", "NEIGHBOR_POSITION_BP",
             "SIGNED_DISTANCE_CM", "DISTANCE_BIN", "NEIGHBOR_MAF_INDEPENDENT",
             "CALLED_HAPLOTYPES_INDEPENDENT", "R2_INDEPENDENT_65",
             "CALLED_HAPLOTYPES_ALL75", "R2_ALL_75"),
            target_pair_rows,
        )
        _write_tsv(
            temporary / "matched_genomewide_windows.tsv",
            ("ANCHOR_ID", "CHROMOSOME", "POSITION_BP", "POSITION_CM", "ANCHOR_MAF",
             "MINOR_COPY_COUNT", "DISTANCE_BIN", "ELIGIBLE_NEIGHBOR_COUNT",
             "EVALUATED_PAIR_COUNT", "MEDIAN_R2"),
            anchor_rows,
        )
        _write_tsv(
            temporary / "ld_v2_summary.tsv",
            ("DISTANCE_BIN", "TARGET_PAIR_COUNT", "TARGET_MEDIAN_R2_INDEPENDENT",
             "TARGET_MEDIAN_R2_ALL75", "LOCAL_BACKGROUND_PAIR_COUNT",
             "LOCAL_BACKGROUND_MEDIAN_R2", "MATCHED_ANCHOR_COUNT", "MATCHED_MEDIAN_R2",
             "MATCHED_Q25_R2", "MATCHED_Q75_R2", "MATCHED_P95_R2",
             "EMPIRICAL_PERCENTILE", "BIN_INTERPRETATION"),
            summary_rows,
        )
        figures = _plot(temporary, summary_rows, target_pair_rows)
        (temporary / "analysis_summary.json").write_text(
            json.dumps(summary, ensure_ascii=False, indent=2, sort_keys=True) + "\n",
            encoding="utf-8",
        )
        (temporary / "CORRECTIONS_PHASE15_LD_V2.md").write_text(
            _markdown(summary, summary_rows), encoding="utf-8",
        )
        cards = "".join(
            f'<figure><img src="{html.escape(path.name)}" alt="Comparaison LD phase 15 v2">'
            '<figcaption>LD DOCK6, fond local et fenêtres génomiques appariées.</figcaption></figure>'
            for path in figures if path.suffix == ".png"
        )
        (temporary / "index.html").write_text(
            "<!doctype html><html lang=\"fr\"><meta charset=\"utf-8\">"
            "<title>Phase 15 LD v2</title><style>body{font-family:system-ui;margin:2rem;"
            "background:#f8fafc;color:#172554;max-width:1500px}figure{background:white;"
            "padding:1rem;border:1px solid #cbd5e1}img{max-width:100%;height:auto}"
            ".warning{background:#fff7ed;color:#7c2d12;padding:1rem;border-left:4px solid #ea580c}"
            "</style><h1>Phase 15 v2 — LD ciblée et fond génomique apparié</h1>"
            f"<p><strong>{html.escape(evidence_direction)}</strong> — {len(selected_anchors)} fenêtres appariées.</p>"
            "<p class=\"warning\">Résultat secondaire et exploratoire : trois familles indépendantes. "
            "La série à 75 sujets est confondue par la parenté.</p>" + cards + "</html>\n",
            encoding="utf-8",
        )
        outputs = sorted(path for path in temporary.iterdir() if path.is_file())
        (temporary / "checksums.sha256").write_text(
            "".join(f"{_sha256(path)}  {path.name}\n" for path in outputs), encoding="utf-8",
        )
        temporary.rename(destination)
    except Exception:
        shutil.rmtree(temporary, ignore_errors=True)
        raise
    return destination


def main() -> int:
    args = _parser().parse_args()
    try:
        destination = execute(args.run_dir, args.output_dir, args.bcftools)
    except (
        DocumentValidationError,
        MatchedLdError,
        OSError,
        Phase15LdV2Error,
        subprocess.SubprocessError,
    ) as error:
        raise SystemExit(f"phase15_ld_v2_failed:{error}") from error
    print(destination)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
