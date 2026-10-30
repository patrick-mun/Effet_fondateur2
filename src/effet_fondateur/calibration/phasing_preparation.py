"""D23a : préparation autosomique liée à la puce, sans exécuter le phasage.

Le modèle de transmission est un processus de Haldane sans interférence.
Les copies fondatrices viennent de groupes D10 distincts ; aucun génotype
du panel de validation n'est ouvert. Les chemins publics sont explicites.
"""

from __future__ import annotations

from collections import Counter, defaultdict
from dataclasses import dataclass
import csv
from datetime import datetime, timezone
import gzip
import hashlib
import json
import logging
from pathlib import Path
import shutil
import subprocess
import tempfile
from typing import Any, Sequence
from urllib.request import urlopen

import numpy as np
import yaml

from effet_fondateur.audit import sha256_file
from effet_fondateur.calibration.export_blinded_inputs import ACPA_COLUMNS, SAMPLE_COLUMNS, pipeline_configuration
from effet_fondateur.calibration.materialize_measurement_errors import deterministic_gzip
from effet_fondateur.calibration.partitioning import parse_population_pedigree, validate_partition_manifest
from effet_fondateur.calibration.simulation.base_haplotypes import parse_population_assignments
from effet_fondateur.calibration.simulation.chip_markers import ChipMarker, ReferenceVariant, harmonize_chip_markers, parse_chip_marker_manifest
from effet_fondateur.calibration.simulation.chip_package import normalize_official_chip_package
from effet_fondateur.calibration.simulation.measurement_errors import read_panel
from effet_fondateur.calibration.technical_scenarios import new_workspace_run
from effet_fondateur.contracts import load_pipeline_config, validate_json_document, validate_tsv_table
from effet_fondateur.contracts.samples import validate_samples_master
from effet_fondateur.references.genetic_maps import ensure_genetic_map_cached, resolve_genetic_map
from effet_fondateur.references.selection import load_reference_selection
from effet_fondateur.references.window import GRCH38_AUTOSOME_LENGTHS


LOGGER = logging.getLogger(__name__)
METHOD_ID = "chip_autosomal_phasing_preparation_v1"
TARGET_SOURCE = "D23A_EXPLICIT_SYNTHETIC_GENOTYPE"


def _validate_document(document: dict[str, Any], schema_name: str) -> None:
    schemas = Path(__file__).resolve().parents[3] / "calibration/schemas"
    validate_json_document(document, schema_name, schemas_dir=schemas)


def write_json(path: Path, document: Any) -> None:
    """Écrit un document déterministe, sans NaN, dans le nouveau run privé."""
    path.write_text(json.dumps(document, indent=2, allow_nan=False) + "\n", encoding="utf-8")


def load_configuration(path: Path) -> dict[str, Any]:
    """Refuse champs inconnus, booléens numériques et paramètres hors domaine."""
    config = json.loads(path.read_text())
    try:
        _validate_document(config, "phasing_preparation_config.schema.json")
    except ValueError:
        raise ValueError("d23a_configuration_invalid") from None
    fixed = dict(schema_version="1.0.0", method_id=METHOD_ID, panel_role="CALIBRATION",
                 target_rule="EXPLICIT_FOUNDER_INDEX_MODULO_3_MENDELIAN_CHILD",
                 reference_rule="EXACT_D10_PHASING_REFERENCE_EXCLUDE_BOTH_DONOR_PARTITIONS")
    integers = ("founder_seed", "meiosis_seed", "phasing_seed", "trio_count", "unrelated_count", "candidate_markers_per_autosome",
                "maximum_markers_per_autosome", "minimum_markers_per_autosome", "minimum_total_markers",
                "target_flank_bp", "source_timeout_seconds")
    if (set(config) != set(fixed) | set(integers) | {"configuration_id", "target_configuration", "target_configuration_sha256", "rare_maf_upper_bound"}
            or any(config.get(key) != value for key, value in fixed.items())
            or any(type(config.get(key)) is not int or config[key] <= 0 for key in integers)
            or not isinstance(config["configuration_id"], str) or not config["configuration_id"]
            or not isinstance(config["target_configuration"], str)
            or config.get("phasing_seed") in {config.get("founder_seed"), config.get("meiosis_seed")}
            or type(config["rare_maf_upper_bound"]) not in (int, float) or not 0 < config["rare_maf_upper_bound"] < 0.5
            or not config["candidate_markers_per_autosome"] >= config["maximum_markers_per_autosome"] >= config["minimum_markers_per_autosome"]):
        raise ValueError("d23a_configuration_invalid")
    return config


def domain_rng(seed: int, *domain: object) -> np.random.Generator:
    """Sous-flux SHA-256 stables, indépendants de l'ordre d'exécution."""
    payload = ":".join((METHOD_ID, str(seed), *(str(value) for value in domain)))
    return np.random.default_rng(int.from_bytes(hashlib.sha256(payload.encode()).digest()[:16], "big"))


def transmission_track(positions_cm: np.ndarray, seed: int, chromosome: int, child_rank: int,
                       parental_side: str) -> tuple[np.ndarray, list[float]]:
    """Retourne la copie parentale 0/1 et les points de recombinaison en cM.

    N suit Poisson(longueur_cM/100), points uniformes conditionnels, pas
    d'interférence ni de carte sexuée. Aucun taux bp/cM arbitraire n'est utilisé.
    """
    if (positions_cm.ndim != 1 or len(positions_cm) < 2 or not np.isfinite(positions_cm).all()
            or np.any(np.diff(positions_cm) < 0) or parental_side not in {"father", "mother"}):
        raise ValueError("d23a_transmission_map_invalid")
    rng = domain_rng(seed, "meiosis", chromosome, child_rank, parental_side)
    event_count = int(rng.poisson((positions_cm[-1] - positions_cm[0]) / 100))
    events = np.sort(rng.uniform(positions_cm[0], positions_cm[-1], event_count))
    initial = int(rng.integers(2))
    track = ((initial + np.searchsorted(events, positions_cm, side="right")) % 2).astype(np.int8)
    return track, events.tolist()


def build_chromosome(founder_alleles: np.ndarray, positions_cm: np.ndarray, chromosome: int,
                     config: dict[str, Any]) -> tuple[np.ndarray, list[dict[str, Any]]]:
    """Assemble fondateurs puis enfants, sans modifier les copies sources."""
    expected = 2 * config["trio_count"] + config["unrelated_count"]
    if founder_alleles.shape != (len(positions_cm), expected, 2) or not np.isin(founder_alleles, [0, 1]).all():
        raise ValueError("d23a_founder_shape_or_allele_invalid")
    result = np.empty((len(positions_cm), expected + config["trio_count"], 2), dtype=np.int8)
    result[:, :expected] = founder_alleles
    transmissions = []
    for child in range(config["trio_count"]):
        for side_index, side in enumerate(("father", "mother")):
            track, events = transmission_track(positions_cm, config["meiosis_seed"], chromosome, child, side)
            parent = child * 2 + side_index
            result[:, expected + child, side_index] = founder_alleles[np.arange(len(positions_cm)), parent, track]
            transmissions.append(dict(child_rank=child, parental_side=side, chromosome=chromosome,
                                      initial_copy=int(track[0]), crossover_positions_cm=events))
    return result, transmissions


def evenly_spaced(items: Sequence[Any], count: int) -> list[Any]:
    """Sélectionne des rangs physiques espacés, sans regarder les génotypes."""
    if count <= 0:
        raise ValueError("d23a_selection_count_invalid")
    if len(items) <= count:
        return list(items)
    return [items[int(index)] for index in np.linspace(0, len(items) - 1, count, dtype=int)]


def _public_path(path: Path) -> None:
    parts = {part.casefold() for part in path.absolute().parts}
    if (not path.exists() or any(item.is_symlink() for item in (path, *path.parents))
            or parts & {"dock6_real", "res_acpa", "temoins", "snp_acpa_accel"}
            or {"data", "input"} <= parts or {"data", "cache"} <= parts):
        raise ValueError("d23a_public_source_path_forbidden")


def _run(command: list[str], timeout: int, ledger: list[dict[str, Any]]) -> str:
    try:
        result = subprocess.run(command, capture_output=True, text=True, timeout=timeout, check=False)
    except (OSError, subprocess.TimeoutExpired) as error:
        raise RuntimeError("d23a_external_command_unavailable_or_timed_out") from error
    ledger.append(dict(command=command, returncode=result.returncode))
    if result.returncode:
        # Ne pas exposer les identifiants que bcftools pourrait citer en erreur.
        raise RuntimeError("d23a_external_command_failed:" + str(result.returncode))
    return result.stdout


def _dual_digest(path: Path) -> tuple[str, str]:
    _public_path(path)
    sha, md5 = hashlib.sha256(), hashlib.md5(usedforsecurity=False)
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            sha.update(block)
            md5.update(block)
    return sha.hexdigest(), md5.hexdigest()


def _metadata_asset(source: dict[str, Any], destination: Path) -> Path:
    url = source["origin"]["canonical_url"]
    if not url.startswith("https://ftp.1000genomes.ebi.ac.uk/"):
        raise ValueError("d23a_metadata_origin_invalid")
    with urlopen(url, timeout=120) as response:
        payload = response.read(10_000_001)
    if len(payload) > 10_000_000 or hashlib.sha256(payload).hexdigest() != source["fingerprint"]["sha256"]:
        raise ValueError("d23a_metadata_sha256_mismatch")
    destination.write_bytes(payload)
    return destination


def _partitions(root: Path, public: dict[str, Any], public_sha: str, metadata: Path) -> tuple[dict[str, Any], dict[str, list[str]], dict[str, str]]:
    if root.parent.name != "calibration_workspace" or any(path.is_symlink() for path in (root, *root.parents)):
        raise ValueError("d23a_partition_root_invalid")
    manifest_path = root / "private_truth/reference_partition_manifest.json"
    manifest = json.loads(manifest_path.read_text())
    # D10 est lié au catalogue 1.0.0 ; la 1.1.0 ajoute la puce. Revalider
    # contre l'instantané exact, jamais réécrire l'engagement historique.
    historical_path = Path(__file__).resolve().parents[3] / "calibration/references/public_reference_catalog_d10_v1.json"
    historical = json.loads(historical_path.read_text())
    historical_sha = sha256_file(historical_path)
    if historical_sha != manifest["source_catalog"]["sha256"]:
        raise ValueError("d23a_historical_partition_catalog_changed")
    old_assets = {row["source_id"]: row for row in historical["sources"]}
    new_assets = {row["source_id"]: row for row in public["sources"]}
    for source_id, asset in old_assets.items():
        if source_id not in new_assets or asset["fingerprint"] != new_assets[source_id]["fingerprint"]:
            raise ValueError("d23a_public_source_version_changed_since_d10")
    memberships = {}
    for row in manifest["partitions"]:
        path = root / row["membership_artifact_path"]
        if path.is_symlink() or not path.is_file() or root.resolve() not in path.resolve().parents or sha256_file(path) != row["membership_sha256"]:
            raise ValueError("d23a_partition_binding_invalid")
        memberships[row["partition_id"]] = path.read_text().splitlines()
    with metadata.open() as handle:
        source_ids, groups = parse_population_pedigree(handle)
    validate_partition_manifest(manifest, memberships_by_partition_id=memberships,
                                source_identifiers=source_ids, independence_group_by_identifier=groups,
                                public_reference_catalog=historical, public_reference_catalog_sha256=historical_sha)
    return manifest, memberships, groups


def choose_founders(members: Sequence[str], groups: dict[str, str], populations: dict[str, str],
                    config: dict[str, Any]) -> list[tuple[tuple[str, int], tuple[str, int]]]:
    """Choisit deux copies par fondateur, sans réutiliser un groupe donneur D10."""
    rng = domain_rng(config["founder_seed"], "founders")
    per_population: dict[str, list[str]] = defaultdict(list)
    used_groups = set()
    # L'ordre ne dépend ni du fichier source ni d'une performance de phasage.
    ranked = sorted(members, key=lambda sample: hashlib.sha256(f'{config["founder_seed"]}:{sample}'.encode()).digest())
    for sample in ranked:
        if groups[sample] not in used_groups:
            per_population[populations[sample]].append(sample)
            used_groups.add(groups[sample])
    count = 2 * config["trio_count"] + config["unrelated_count"]
    result = []
    while len(result) < count:
        progressed = False
        for population in sorted(per_population):
            candidates = per_population[population]
            if len(candidates) >= 2 and len(result) < count:
                first, second = candidates.pop(), candidates.pop()
                result.append(((first, int(rng.integers(2))), (second, int(rng.integers(2)))))
                progressed = True
        if not progressed:
            raise ValueError("d23a_independent_donor_groups_insufficient")
    return result


def _write_vcf(path: Path, variants: list[tuple[str, ...]], alleles: np.ndarray, samples: list[str], *, phased: bool) -> None:
    with deterministic_gzip(path) as handle:
        handle.write('##fileformat=VCFv4.2\n##reference=GRCh38\n')
        for chromosome, length in GRCH38_AUTOSOME_LENGTHS.items():
            handle.write(f"##contig=<ID=chr{chromosome},length={length}>\n")
        handle.write('##FORMAT=<ID=GT,Number=1,Type=String,Description="Synthetic genotype">\n')
        handle.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t" + "\t".join(samples) + "\n")
        for row, calls in zip(variants, alleles, strict=True):
            formatted = [f"{pair[0]}|{pair[1]}" if phased else "/".join(map(str, sorted(pair))) for pair in calls]
            handle.write("\t".join((*row, ".", "PASS", ".", "GT", *formatted)) + "\n")


def _export_observations(root: Path, samples: list[str], variants: list[tuple[str, ...]], alleles: np.ndarray,
                         config: dict[str, Any], target: dict[str, Any], final_root: Path, repository: Path) -> None:
    acpa = root / "acpa_samples"
    acpa.mkdir()
    metadata = root / "metadata"
    metadata.mkdir()
    target_index = next(index for index, row in enumerate(variants) if row[2] == target["project_variant_id"])
    genotypes = []
    for sample_index, sample in enumerate(samples):
        with (acpa / f"{sample}.txt").open("w", newline="") as handle:
            handle.write("# Array Type Name: CytoScan 750K Accel Array\n# UCSC Genomic Version: hg38\n")
            writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
            writer.writerow(ACPA_COLUMNS)
            for row, pair in zip(variants, alleles[:, sample_index], strict=True):
                calls = sorted((row[3], row[4])[int(value)] for value in pair)
                writer.writerow((row[2], "".join(calls), "", row[0].removeprefix("chr"), row[1]))
        genotypes.append("/".join(sorted((target["ref"], target["alt"])[int(value)] for value in alleles[target_index, sample_index])))
    founders = 2 * config["trio_count"] + config["unrelated_count"]
    with (metadata / "samples.tsv").open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(SAMPLE_COLUMNS)
        for index, sample in enumerate(samples):
            if index < 2 * config["trio_count"]:
                family, father, mother, sex = f"F{index // 2:04d}", "0", "0", "MALE" if index % 2 == 0 else "FEMALE"
            elif index >= founders:
                child = index - founders
                family, father, mother, sex = f"F{child:04d}", samples[2 * child], samples[2 * child + 1], "UNKNOWN"
            else:
                family, father, mother, sex = f"U{index:04d}", "0", "0", "UNKNOWN"
            writer.writerow((sample, f"{sample}.txt", family, sample, father, mother, sex, "UNKNOWN", "UNSPECIFIED",
                             genotypes[index], TARGET_SOURCE, "SYNTHETIC", "true", "true", "D23A_SIMULATION"))
    with (root / "target_genotypes.tsv").open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(("SAMPLE_ID", "GENOTYPE", "GENOTYPE_SOURCE", "MEASUREMENT_METHOD"))
        writer.writerows((sample, gt, TARGET_SOURCE, "SIMULATED_NOT_EXPERIMENTAL") for sample, gt in zip(samples, genotypes, strict=True))
    (root / "target_metadata.yaml").write_text(yaml.safe_dump(target, sort_keys=False))
    profile = pipeline_configuration(target)
    profile["project"]["project_id"] = "d23a_inputs_only"
    profile["target"] = {key: target[key] for key in profile["target"]}
    profile["inputs"].update(acpa_samples_dir=str(final_root / "blinded_inputs/acpa_samples"),
                             sample_metadata=str(final_root / "blinded_inputs/metadata/samples.tsv"),
                             target_variant_metadata=str(final_root / "blinded_inputs/target_metadata.yaml"),
                             target_variant_genotypes=str(final_root / "blinded_inputs/target_genotypes.tsv"),
                             reference_panel_catalog=str(repository / "config/references/1000g_high_coverage_grch38_phased.json"),
                             ancestry_reference_catalog=str(repository / "config/references/1000g_ancestry_grch38.json"),
                             genetic_map_catalog=str(repository / "config/references/shapeit4_hapmap_grch38.json"),
                             phasing_reference_selection=str(final_root / "blinded_inputs/reference_selection.json"))
    profile["stages"]["validate_sources"]["parameters"]["required_autosomes"] = list(range(1, 23))
    # Profil volontairement borné à 00–02 : ni revues scientifiques ni succès
    # de cohortes/phasing ne sont fabriqués pendant la préparation.
    (root / "pipeline.yaml").write_text(yaml.safe_dump(profile, sort_keys=False))
    load_pipeline_config(root / "pipeline.yaml")
    (root / "phasing_parameters.yaml").write_text(yaml.safe_dump({
        name: dict(require_reference_sample_selection=True, shapeit5_threads=1, shapeit5_seed=config["phasing_seed"])
        for name in ("phase_target_region", "phase_autosomal_panel")}, sort_keys=False))


def _source_rows(text: str, chromosome: int, donor_count: int) -> tuple[list[ReferenceVariant], dict[int, np.ndarray]]:
    variants, genotypes = [], {}
    duplicated = set()
    for line in text.splitlines():
        fields = line.split("\t")
        if len(fields) != 5 + donor_count:
            raise ValueError("d23a_source_query_width_invalid")
        if fields[0] != f"chr{chromosome}" or len(fields[3]) != 1 or len(fields[4]) != 1 or fields[3] not in "ACGT" or fields[4] not in "ACGT":
            continue
        position = int(fields[1])
        variants.append(ReferenceVariant(chromosome, position, fields[2], fields[3], fields[4]))
        if position in genotypes:
            duplicated.add(position)
        if all(value in {"0|0", "0|1", "1|0", "1|1"} for value in fields[5:]):
            genotypes[position] = np.array([[int(gt[0]), int(gt[2])] for gt in fields[5:]], dtype=np.int8)
    for position in duplicated:
        genotypes.pop(position, None)
    return variants, genotypes



@dataclass(frozen=True)
class ChromosomePreparationContext:
    """Sources explicites communes aux 22 extractions, jamais de vérité validation."""

    source_directory: Path
    technical: dict[str, Any]
    source_specs: dict[int, dict[str, Any]]
    chip_markers: Sequence[ChipMarker]
    map_directory: Path
    target: dict[str, Any]
    config: dict[str, Any]
    origins: list[tuple[tuple[str, int], tuple[str, int]]]
    donors: list[str]
    private: Path
    resources: Path
    ledger: list[dict[str, Any]]
    bcftools: str
    donor_index: dict[str, int]


def _prepare_chromosome(chromosome: int, context: ChromosomePreparationContext) -> tuple[list[tuple[str, ...]], np.ndarray, dict[str, Any], list[dict[str, Any]], dict[str, Any]]:
    """Vérifie puis construit un autosome borné ; bp et cM restent distincts."""
    source_directory = context.source_directory
    technical = context.technical
    source_specs = context.source_specs
    chip_markers = context.chip_markers
    map_directory = context.map_directory
    target = context.target
    config = context.config
    origins = context.origins
    donors = context.donors
    private = context.private
    resources = context.resources
    ledger = context.ledger
    bcftools = context.bcftools
    donor_index = context.donor_index
    LOGGER.info("D23a : autosome %d/22 — empreintes puis extraction bornée", chromosome)
    vcf = source_directory / technical["vcf_filename_template"].format(chromosome=chromosome)
    index = source_directory / technical["index_filename_template"].format(chromosome=chromosome)
    vcf_sha, vcf_md5 = _dual_digest(vcf)
    index_sha, index_md5 = _dual_digest(index)
    if vcf_md5 != source_specs[chromosome]["vcf_md5"] or index_md5 != source_specs[chromosome]["index_md5"]:
        raise ValueError(f"d23a_public_source_md5_invalid:chr{chromosome}")
    fingerprint = dict(chromosome=chromosome, vcf_sha256=vcf_sha, index_sha256=index_sha,
                       vcf_md5=vcf_md5, index_md5=index_md5)
    with (map_directory / f"chr{chromosome}.grch38.tsv").open() as handle:
        map_rows = list(csv.DictReader(handle, delimiter="\t"))
    map_bp = np.array([int(row["POSITION_BP"]) for row in map_rows])
    map_cm = np.array([float(row["POSITION_CM"]) for row in map_rows])
    if len(map_bp) < 2 or np.any(np.diff(map_bp) <= 0) or not np.isfinite(map_cm).all() or np.any(np.diff(map_cm) < 0):
        raise ValueError("d23a_map_invalid")
    chromosome_markers = [marker for marker in chip_markers if marker.chromosome == chromosome]
    position_counts = Counter(marker.position_bp for marker in chromosome_markers)
    eligible = sorted((marker for marker in chromosome_markers if position_counts[marker.position_bp] == 1
                       and map_bp[0] <= marker.position_bp <= map_bp[-1]
                       and not (marker.strand == "UNKNOWN" and {marker.allele_a, marker.allele_b} in ({"A", "T"}, {"C", "G"}))), key=lambda marker: marker.position_bp)
    candidates = {marker.position_bp: marker for marker in evenly_spaced(eligible, config["candidate_markers_per_autosome"])}
    if chromosome == target["chromosome"]:
        candidates.update({marker.position_bp: marker for marker in eligible
                           if abs(marker.position_bp - target["position_bp"]) <= config["target_flank_bp"]})
    position_file = private / f"chr{chromosome}.candidates.tsv"
    position_file.write_text("".join(f"chr{chromosome}\t{position}\t{position}\n" for position in sorted(candidates)))
    queried = _run([bcftools, "query", "--regions-file", str(position_file), "--samples-file", str(private / "donors.ids.txt"),
                    "-f", "%CHROM\t%POS\t%ID\t%REF\t%ALT[\t%GT]\n", str(vcf)], config["source_timeout_seconds"], ledger)
    reference_variants, source_calls = _source_rows(queried, chromosome, len(donors))
    matches, match_audit = harmonize_chip_markers(candidates.values(), reference_variants,
                                                chromosome=chromosome, start_bp=int(map_bp[0]), end_bp=int(map_bp[-1]))
    complete = sorted((match for match in matches if match.marker.position_bp in source_calls), key=lambda match: match.marker.position_bp)
    selected = {match.marker.position_bp: match for match in evenly_spaced(complete, config["maximum_markers_per_autosome"])}
    if chromosome == target["chromosome"]:
        selected.update({match.marker.position_bp: match for match in complete
                         if abs(match.marker.position_bp - target["position_bp"]) <= config["target_flank_bp"]})
        selected.pop(target["position_bp"], None)
    ordered = sorted(selected)
    if len(ordered) < config["minimum_markers_per_autosome"]:
        raise ValueError(f"d23a_insufficient_chip_markers:chr{chromosome}")
    chromosome_rows = [(f"chr{chromosome}", str(position), selected[position].marker.marker_id,
                        selected[position].reference_variant.reference_allele, selected[position].reference_variant.alternate_allele) for position in ordered]
    founder_calls = np.array([[[source_calls[position][donor_index[origin[0]], origin[1]] for origin in pair] for pair in origins] for position in ordered], dtype=np.int8)
    np.savez_compressed(private / f"chr{chromosome}.donor_calls.npz", positions=np.array(ordered),
                        calls=np.array([source_calls[position] for position in ordered]))
    if chromosome == target["chromosome"]:
        if not map_bp[0] <= target["position_bp"] <= map_bp[-1]:
            raise ValueError("d23a_target_outside_map")
        target_rank = int(np.searchsorted(ordered, target["position_bp"]))
        chromosome_rows.insert(target_rank, (f"chr{chromosome}", str(target["position_bp"]), target["project_variant_id"], target["ref"], target["alt"]))
        target_calls = np.array([((0, 0), (0, 1), (1, 1))[index % 3] for index in range(len(origins))], dtype=np.int8)
        founder_calls = np.insert(founder_calls, target_rank, target_calls, axis=0)
    positions_cm = np.interp([int(row[1]) for row in chromosome_rows], map_bp, map_cm)
    values, child_transmissions = build_chromosome(founder_calls, positions_cm, chromosome, config)
    reference_file = resources / f"chr{chromosome}.reference.vcf.gz"
    positions = resources / f"chr{chromosome}.positions.tsv"
    positions.write_text("".join(f"chr{chromosome}\t{position}\t{position}\n" for position in ordered))
    _run([bcftools, "view", "--no-version", "--regions-file", str(positions), "--samples-file", str(private / "reference.ids.txt"),
          "--min-alleles", "2", "--max-alleles", "2", "--types", "snps", "-Oz", "-o", str(reference_file), str(vcf)], config["source_timeout_seconds"], ledger)
    _run([bcftools, "index", "--tbi", str(reference_file)], 120, ledger)
    actual_samples = _run([bcftools, "query", "--list-samples", str(reference_file)], 120, ledger).splitlines()
    if actual_samples != reference:
        raise ValueError("d23a_reference_membership_or_order_invalid")
    reference_count = int(_run([bcftools, "index", "--nrecords", str(reference_file)], 120, ledger).strip())
    if reference_count != len(ordered):
        raise ValueError("d23a_reference_variant_count_invalid")
    query_sites = _run([bcftools, "query", "-f", "%CHROM\t%POS\t%REF\t%ALT\n", str(reference_file)], 120, ledger)
    expected_sites = [(row[0], int(row[1]), row[3], row[4]) for row in chromosome_rows if row[2] != target["project_variant_id"]]
    actual_sites = [(fields[0], int(fields[1]), fields[2], fields[3]) for fields in (line.split("\t") for line in query_sites.splitlines())]
    if actual_sites != expected_sites:
        raise ValueError("d23a_reference_variant_identity_invalid")
    invalid = _run([bcftools, "query", "-i", 'GT="mis" || GT~"/"', "-f", "%CHROM\t%POS\n", str(reference_file)], 120, ledger)
    if invalid.strip():
        raise ValueError("d23a_reference_not_complete_phased")
    maf = np.minimum(values.mean(axis=(1, 2)), 1 - values.mean(axis=(1, 2)))
    hetero = values[:, :, 0] != values[:, :, 1]
    np.savetxt(resources / f"chr{chromosome}.study_map.tsv", np.column_stack(([int(row[1]) for row in chromosome_rows], positions_cm)),
               fmt=("%d", "%.17g"), delimiter="\t", header="POSITION_BP\tPOSITION_CM", comments="")
    summary = dict(chromosome=chromosome, chip_marker_count=len(chromosome_markers), eligible_marker_count=len(eligible),
                          candidate_marker_count=len(candidates), harmonized_marker_count=len(matches),
                          incomplete_donor_marker_count=len(matches) - len(complete), retained_chip_marker_count=len(ordered),
                          observed_marker_count=len(chromosome_rows), reference_variant_count=reference_count,
                          harmonization_exclusions=dict(match_audit.exclusion_counts),
                          rare_variant_count=int(((maf > 0) & (maf <= config["rare_maf_upper_bound"])).sum()),
                          rare_heterozygote_count=int(hetero[(maf > 0) & (maf <= config["rare_maf_upper_bound"])].sum()),
                          common_heterozygote_count=int(hetero[maf > config["rare_maf_upper_bound"]].sum()),
                          monomorphic_variant_count=int((maf == 0).sum()), meiosis_count=len(child_transmissions))
    return chromosome_rows, values, summary, child_transmissions, fingerprint


def prepare_panel(*, configuration: Path, partition_root: Path, source_directory: Path, chip_package: Path,
                  map_archive: Path | None, workspace_root: Path, run_id: str, bcftools_command: str = "bcftools") -> dict[str, Any]:
    """Matérialise et relit un nouveau panel D23a ; ne lance pas SHAPEIT5.

    Les copies donneurs et croisements ne sont écrits que dans private_truth.
    Les profils historiques ne sont ni modifiés ni écrasés. Les sources
    publiques sont intégralement vérifiées ; un échec ne publie pas de PASS.
    """
    destination = new_workspace_run(workspace_root, run_id)
    destination.parent.mkdir(parents=True, exist_ok=True)
    for path in (source_directory, chip_package, *((map_archive,) if map_archive is not None else ())):
        _public_path(path)
    config = load_configuration(configuration)
    repository = Path(__file__).resolve().parents[3]
    public_path = repository / "calibration/references/public_reference_catalog.json"
    public = json.loads(public_path.read_text())
    assets = {row["source_id"]: row for row in public["sources"]}
    technical_path = repository / "config/references/1000g_high_coverage_grch38_phased.json"
    technical = json.loads(technical_path.read_text())["panels"][0]
    binding = assets["igsr_1kg_3202_phased_grch38_20220422"]["catalog_binding"]["sha256"]
    if sha256_file(technical_path) != binding:
        raise ValueError("d23a_technical_catalog_binding_invalid")
    target_path = repository / config["target_configuration"]
    if repository not in target_path.resolve().parents:
        raise ValueError("d23a_target_configuration_path_invalid")
    target = json.loads(target_path.read_text())["target_metadata"]
    if sha256_file(target_path) != config["target_configuration_sha256"]:
        raise ValueError("d23a_public_target_configuration_changed")
    validate_json_document(target, "target_variant_metadata.schema.json")
    bcftools = shutil.which(bcftools_command)
    if bcftools is None:
        raise ValueError("d23a_bcftools_not_found")
    temporary = Path(tempfile.mkdtemp(prefix=f".{run_id}.", dir=destination.parent))
    try:
        temporary.chmod(0o700)
        private, blind, resources = (temporary / name for name in ("private_truth", "blinded_inputs", "reference_assets"))
        for directory in (private, blind, resources):
            directory.mkdir(mode=0o700)
        ledger: list[dict[str, Any]] = []
        version = _run([bcftools, "--version"], 30, ledger).splitlines()[0]
        LOGGER.info("D23a : vérification des métadonnées publiques et des partitions D10")
        metadata = _metadata_asset(assets["igsr_1kg_3202_population_metadata_20130606"], resources / "population_metadata.tsv")
        manifest, memberships, groups = _partitions(partition_root, public, sha256_file(public_path), metadata)
        with metadata.open() as handle:
            populations = parse_population_assignments(handle)
        origins = choose_founders(memberships["donor_calibration"], groups, populations, config)
        donors = sorted({origin[0] for pair in origins for origin in pair})
        if len(donors) != 2 * len(origins) or len({groups[sample] for sample in donors}) != len(donors):
            raise ValueError("d23a_founder_donor_group_reuse")
        reference = memberships["phasing_reference"]
        partition_sha = sha256_file(partition_root / "private_truth/reference_partition_manifest.json")
        selection = dict(schema_version="1.0.0", method_id="explicit_partition_reference_selection_v1", assembly="GRCh38",
                         source_panel_id=technical["panel_id"], partition_manifest_sha256=partition_sha,
                         membership_sha256=hashlib.sha256("".join(sample + "\n" for sample in reference).encode()).hexdigest(),
                         sample_count=len(reference), sample_ids=reference)
        selection_path = blind / "reference_selection.json"
        write_json(selection_path, selection)
        load_reference_selection(selection_path, expected_sha256=sha256_file(selection_path), panel_id=technical["panel_id"])
        (private / "donors.ids.txt").write_text("".join(sample + "\n" for sample in donors))
        (private / "reference.ids.txt").write_text("".join(sample + "\n" for sample in reference))
        write_json(private / "founder_origins.json", origins)
        write_json(private / "configuration.json", config)
        LOGGER.info("D23a : annotation officielle de la puce et cartes GRCh38")
        chip_sha = assets["thermofisher_cytoscan750k_accel_na36_r1"]["fingerprint"]["sha256"]
        marker_path = resources / "chip_markers.tsv"
        chip_audit = normalize_official_chip_package(package_path=chip_package, expected_package_sha256=chip_sha, output_path=marker_path)
        with marker_path.open() as handle:
            chip_markers = parse_chip_marker_manifest(handle)
        map_catalog = repository / "config/references/shapeit4_hapmap_grch38.json"
        resolved_map = resolve_genetic_map(map_catalog, "shapeit4_hapmap_grch38_7a0cab7", "GRCh38")
        if map_archive is not None and sha256_file(map_archive) != resolved_map.archive_sha256:
            raise ValueError("d23a_map_archive_hash_invalid")
        def copy_map(_url: str, destination_path: Path, _timeout: int, _chunk: int) -> None:
            shutil.copyfile(map_archive, destination_path)
        maps = ensure_genetic_map_cached(resolved=resolved_map, chromosome=1, cache_root=resources / "genetic_maps",
                                        **({"downloader": copy_map} if map_archive is not None else {}))
        map_directory = maps.map_path.parent
        rows, matrices, summaries, transmissions, fingerprints = [], [], [], [], []
        source_specs = {row["chromosome"]: row for row in technical["chromosomes"]}
        donor_index = {sample: index for index, sample in enumerate(donors)}
        sample_count = len(origins) + config["trio_count"]
        samples = ["S" + hashlib.sha256(f'{config["founder_seed"]}:sample:{index}'.encode()).hexdigest()[:16] for index in range(sample_count)]
        context = ChromosomePreparationContext(
            source_directory, technical, source_specs, chip_markers, map_directory, target, config, origins, donors, private, resources, ledger, bcftools, donor_index
        )
        for chromosome in range(1, 23):
            chromosome_rows, values, summary, events, fingerprint = _prepare_chromosome(chromosome, context)
            rows.extend(chromosome_rows)
            matrices.append(values)
            summaries.append(summary)
            transmissions.extend(events)
            fingerprints.append(fingerprint)
        alleles = np.concatenate(matrices)
        if len(rows) < config["minimum_total_markers"]:
            raise ValueError("d23a_insufficient_total_markers")
        _write_vcf(private / "truth.phased.vcf.gz", rows, alleles, samples, phased=True)
        _write_vcf(blind / "observed.vcf.gz", rows, alleles, samples, phased=False)
        write_json(private / "transmissions.json", transmissions)
        _export_observations(blind, samples, rows, alleles, config, target, destination, repository)
        LOGGER.info("D23a : relecture indépendante des observations et transmissions exportées")
        audit = audit_exports(temporary, config)
        if any(value != 0 for name, value in audit.items() if name.endswith("mismatch_count")):
            raise ValueError("d23a_export_audit_failed")
        write_json(private / "command_ledger.json", ledger)
        result = dict(schema_version="1.0.0", method_id=METHOD_ID, preparation_state="FINALIZED",
                      created_at=datetime.now(timezone.utc).isoformat(), configuration=config,
                      configuration_sha256=sha256_file(configuration), public_catalog_sha256=sha256_file(public_path),
                      target_configuration_sha256=sha256_file(target_path), partition_manifest_sha256=partition_sha,
                      reference_selection_sha256=sha256_file(selection_path), chip_package_sha256=chip_audit.package_sha256,
                      chip_marker_table_sha256=chip_audit.normalized_marker_sha256,
                      map_archive_sha256=resolved_map.archive_sha256, tool_versions=[version],
                      counts=dict(samples=sample_count, founders=len(origins), trios=config["trio_count"],
                                  independent_simulation_units=config["trio_count"] + config["unrelated_count"],
                                  source_donor_groups=len(donors), reference_samples=len(reference), markers=len(rows),
                                  calibration_donors=len(memberships["donor_calibration"]), validation_donors=len(memberships["donor_validation"])),
                      chromosomes=summaries, source_fingerprints=fingerprints, audit=audit,
                      phasing_executed=False, performance_status="NOT_MEASURED_PREPARATION_ONLY", calibrated=False, founder_effect_proven=False)
        result["artifacts"] = inventory(temporary)
        write_json(temporary / "preparation_manifest.json", result)
        for path in temporary.rglob("*"):
            path.chmod(0o700 if path.is_dir() else 0o600)
        verify_preparation(temporary)
        temporary.rename(destination)
        return result
    finally:
        if temporary.exists():
            shutil.rmtree(temporary)


def inventory(root: Path) -> list[dict[str, Any]]:
    """Inventaire exhaustif avant fermeture ; pas d'identifiants dans le résumé."""
    return [dict(path=path.relative_to(root).as_posix(), sha256=sha256_file(path), size_bytes=path.stat().st_size)
            for path in sorted(root.rglob("*")) if path.is_file() and path.name != "preparation_manifest.json"]


def audit_exports(root: Path, config: dict[str, Any]) -> dict[str, int]:
    """Relit les VCF, appels ACPA, registres et transmissions, sans phasage.

    Les comptes sont exhaustifs sur ce banc, pas des estimations de performance.
    Une phase inversée globalement n'est pas testée ici : elle relève de D23.
    """
    with gzip.open(root / "private_truth/truth.phased.vcf.gz", "rt") as handle:
        truth = read_panel(handle, phased=True)
    with gzip.open(root / "blinded_inputs/observed.vcf.gz", "rt") as handle:
        observed = read_panel(handle, phased=False)
    if truth.samples != observed.samples or truth.variants != observed.variants:
        raise ValueError("d23a_truth_observation_axis_mismatch")
    dosage_mismatch = int(np.count_nonzero(np.sort(truth.alleles, axis=2) != observed.alleles))
    founder_count = 2 * config["trio_count"] + config["unrelated_count"]
    transmissions = json.loads((root / "private_truth/transmissions.json").read_text())
    donors = (root / "private_truth/donors.ids.txt").read_text().splitlines()
    origins = json.loads((root / "private_truth/founder_origins.json").read_text())
    donor_index = {sample: index for index, sample in enumerate(donors)}
    if len(transmissions) != 22 * 2 * config["trio_count"]:
        raise ValueError("d23a_transmission_count_mismatch")
    transmission_mismatch, founder_mismatch = 0, 0
    for chromosome in range(1, 23):
        selected = [index for index, row in enumerate(truth.variants) if row[0] == f"chr{chromosome}"]
        map_values = np.loadtxt(root / f"reference_assets/chr{chromosome}.study_map.tsv", delimiter="\t", skiprows=1)
        if [int(row[0]) for row in map_values] != [int(truth.variants[index][1]) for index in selected]:
            raise ValueError("d23a_export_map_axis_mismatch")
        with np.load(root / f"private_truth/chr{chromosome}.donor_calls.npz", allow_pickle=False) as source:
            by_position = {int(position): index for index, position in enumerate(source["positions"])}
            calls = source["calls"]
            for index in selected:
                position = int(truth.variants[index][1])
                if position in by_position:
                    expected_founders = np.array([[calls[by_position[position], donor_index[origin[0]], origin[1]]
                                                  for origin in pair] for pair in origins], dtype=np.int8)
                    founder_mismatch += int(np.count_nonzero(expected_founders != truth.alleles[index, :founder_count]))
        expected, events = build_chromosome(truth.alleles[selected, :founder_count], map_values[:, 1], chromosome, config)
        actual_events = [row for row in transmissions if row["chromosome"] == chromosome]
        # Dix-sept chiffres conservent exactement les coordonnées float64.
        if [(row["child_rank"], row["parental_side"], row["initial_copy"]) for row in actual_events] != [(row["child_rank"], row["parental_side"], row["initial_copy"]) for row in events]:
            raise ValueError("d23a_transmission_replay_mismatch")
        for first, second in zip(actual_events, events, strict=True):
            if len(first["crossover_positions_cm"]) != len(second["crossover_positions_cm"]) or not np.allclose(first["crossover_positions_cm"], second["crossover_positions_cm"], rtol=1e-12, atol=1e-12):
                raise ValueError("d23a_transmission_events_mismatch")
        transmission_mismatch += int(np.count_nonzero(expected != truth.alleles[selected]))
    metadata_path = root / "blinded_inputs/metadata/samples.tsv"
    metadata = validate_samples_master(metadata_path, root / "blinded_inputs/acpa_samples")
    if tuple(row["SAMPLE_ID"] for row in metadata.rows) != truth.samples:
        raise ValueError("d23a_metadata_sample_order_invalid")
    with metadata_path.open() as handle:
        metadata_rows = list(csv.DictReader(handle, delimiter="\t"))
    target_path = root / "blinded_inputs/target_metadata.yaml"
    target = yaml.safe_load(target_path.read_text())
    target_indices = [index for index, row in enumerate(truth.variants) if row[2] == target["project_variant_id"]]
    if len(target_indices) != 1:
        raise ValueError("d23a_target_missing_or_duplicate")
    target_index = target_indices[0]
    expected_founder_target = np.array([((0, 0), (0, 1), (1, 1))[index % 3] for index in range(founder_count)], dtype=np.int8)
    founder_mismatch += int(np.count_nonzero(expected_founder_target != truth.alleles[target_index, :founder_count]))
    target_table = validate_tsv_table(root / "blinded_inputs/target_genotypes.tsv", "target_variant_genotypes.schema.json")
    if [row["SAMPLE_ID"] for row in target_table.rows] != list(truth.samples):
        raise ValueError("d23a_target_genotype_sample_order_invalid")
    acpa_mismatch, target_mismatch, pedigree_mismatch = 0, 0, 0
    for sample_index, sample in enumerate(truth.samples):
        row = metadata_rows[sample_index]
        if row["CLINICAL_STATUS"] != "UNKNOWN" or row["GROUP_LABEL"] != "UNSPECIFIED" or row["TARGET_GENOTYPE_SOURCE"] != TARGET_SOURCE:
            raise ValueError("d23a_unauthorized_metadata_inference")
        if sample_index >= founder_count:
            child = sample_index - founder_count
            pedigree_mismatch += row["PID"] != truth.samples[2 * child] or row["MID"] != truth.samples[2 * child + 1]
        else:
            pedigree_mismatch += (row["PID"], row["MID"]) != ("0", "0")
        with (root / "blinded_inputs/acpa_samples" / f"{sample}.txt").open() as handle:
            next(handle)
            next(handle)
            records = list(csv.DictReader(handle, delimiter="\t"))
        if len(records) != len(truth.variants):
            raise ValueError("d23a_acpa_variant_count_invalid")
        for index, (record, variant) in enumerate(zip(records, truth.variants, strict=True)):
            expected = sorted((variant[3], variant[4])[int(value)] for value in truth.alleles[index, sample_index])
            acpa_mismatch += (record["Probe Set ID"] != variant[2] or record["Chromosome"] != variant[0].removeprefix("chr")
                              or record["Chromosomal Position"] != variant[1] or sorted(record["Forward Strand Base Calls"]) != expected)
        expected_target = "/".join(sorted((target["ref"], target["alt"])[int(value)] for value in truth.alleles[target_index, sample_index]))
        target_mismatch += row["TARGET_GENOTYPE"] != expected_target or target_table.rows[sample_index]["GENOTYPE"] != expected_target
    return dict(compared_genotype_cells=int(len(truth.variants) * len(truth.samples)),
                genotype_mismatch_count=dosage_mismatch, acpa_mismatch_count=acpa_mismatch, founder_copy_mismatch_count=founder_mismatch,
                transmission_mismatch_count=transmission_mismatch, target_mismatch_count=target_mismatch,
                pedigree_mismatch_count=pedigree_mismatch)


def verify_preparation(root: Path, *, require_blinded_seeds: bool = True) -> dict[str, Any]:
    """Exige l'inventaire fermé exact et les empreintes avant rendu ou usage."""
    if any(path.is_symlink() for path in (root, *root.parents)) or not root.is_dir() or any(path.is_symlink() for path in root.rglob("*")):
        raise ValueError("d23a_preparation_symlink_or_root_invalid")
    manifest = json.loads((root / "preparation_manifest.json").read_text())
    _validate_document(manifest, "phasing_preparation_manifest.schema.json")
    if manifest.get("preparation_state") != "FINALIZED" or manifest.get("phasing_executed") is not False or manifest.get("artifacts") != inventory(root):
        raise ValueError("d23a_preparation_inventory_invalid")
    if require_blinded_seeds:
        config = load_configuration(root / "private_truth/configuration.json")
        parameters = yaml.safe_load((root / "blinded_inputs/phasing_parameters.yaml").read_text())
        if config != manifest["configuration"] or any(
                values.get("shapeit5_seed") != config["phasing_seed"] or values.get("require_reference_sample_selection") is not True
                for values in parameters.values()) or set(parameters) != {"phase_target_region", "phase_autosomal_panel"}:
            raise ValueError("d23a_blinded_phasing_parameters_invalid")
    return manifest


def finalize_preparation(*, source_root: Path, configuration: Path, workspace_root: Path, run_id: str) -> dict[str, Any]:
    """Révise uniquement les paramètres aveugles, sans réextraire de génotype.

    La première préparation conserve ses empreintes et sa graine de simulation.
    Cette clôture ne permet de changer que l'identifiant de configuration et
    la graine du futur phasage. Chaque autre fichier conserve son SHA-256 ;
    les exports et transmissions sont relus avant publication du nouveau run.
    """
    destination = new_workspace_run(workspace_root, run_id)
    if source_root.absolute().parent != destination.parent:
        raise ValueError("d23a_revision_source_root_invalid")
    source = verify_preparation(source_root, require_blinded_seeds=False)
    source_sha = sha256_file(source_root / "preparation_manifest.json")
    config = load_configuration(configuration)
    for key, value in source["configuration"].items():
        if key not in {"configuration_id", "phasing_seed"} and config.get(key) != value:
            raise ValueError("d23a_revision_scientific_configuration_change_forbidden")
    implementation_sha = sha256_file(Path(__file__))
    temporary = Path(tempfile.mkdtemp(prefix=f".{run_id}.", dir=destination.parent))
    try:
        shutil.copytree(source_root, temporary, dirs_exist_ok=True)
        write_json(temporary / "private_truth/configuration.json", config)
        profile_path = temporary / "blinded_inputs/pipeline.yaml"
        profile_path.write_text(profile_path.read_text().replace(str(source_root.absolute()), str(destination)))
        load_pipeline_config(profile_path)
        parameter_path = temporary / "blinded_inputs/phasing_parameters.yaml"
        parameters = yaml.safe_load(parameter_path.read_text())
        if set(parameters) != {"phase_target_region", "phase_autosomal_panel"}:
            raise ValueError("d23a_revision_phasing_parameters_invalid")
        for values in parameters.values():
            values["shapeit5_seed"] = config["phasing_seed"]
        parameter_path.write_text(yaml.safe_dump(parameters, sort_keys=False))
        mutable = {"private_truth/configuration.json", "blinded_inputs/pipeline.yaml", "blinded_inputs/phasing_parameters.yaml"}
        old_files = {row["path"]: row["sha256"] for row in source["artifacts"]}
        new_inventory = inventory(temporary)
        if set(old_files) != {row["path"] for row in new_inventory} or any(
                row["sha256"] != old_files[row["path"]] for row in new_inventory if row["path"] not in mutable):
            raise ValueError("d23a_revision_genetic_artifact_changed")
        audit = audit_exports(temporary, config)
        if audit != source["audit"]:
            raise ValueError("d23a_revision_export_audit_changed")
        manifest = dict(source, configuration=config, configuration_sha256=sha256_file(configuration),
                        created_at=datetime.now(timezone.utc).isoformat(), source_preparation_sha256=source_sha,
                        implementation_sha256=implementation_sha, artifacts=new_inventory)
        write_json(temporary / "preparation_manifest.json", manifest)
        verify_preparation(temporary)
        if sha256_file(source_root / "preparation_manifest.json") != source_sha or sha256_file(Path(__file__)) != implementation_sha:
            raise ValueError("d23a_revision_source_or_implementation_changed")
        temporary.rename(destination)
        return manifest
    finally:
        if temporary.exists():
            shutil.rmtree(temporary)
