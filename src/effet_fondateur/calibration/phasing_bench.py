"""Banc D23 de composants réels, avec amont explicitement simulé.

Les outils ne reçoivent que les observations et les références publiques.
La vérité n'est ouverte que par l'évaluateur, dans un autre module.
"""
from __future__ import annotations

import argparse
import gzip
import json
import logging
from pathlib import Path
import shutil
import subprocess
import sys
from typing import Any, Sequence
from urllib.parse import urlparse
from urllib.request import urlopen

import numpy as np
import yaml

from effet_fondateur.audit import sha256_file
from effet_fondateur.calibration.evaluate_technical import raw_table
from effet_fondateur.calibration.phasing_preparation import verify_preparation
from effet_fondateur.calibration.qc_bench import file_inventory, write_json, write_table
from effet_fondateur.calibration.run_qc_bench import checked_stage
from effet_fondateur.calibration.run_structure_bench import component_inputs
from effet_fondateur.calibration.simulation.measurement_errors import MeasurementPanel, read_panel
from effet_fondateur.calibration.technical_scenarios import new_workspace_run
from effet_fondateur.contracts import build_file_artifact, load_pipeline_config, validate_json_document, validate_tsv_table
from effet_fondateur.references import cache_reference_panel, resolve_reference_panel

LOGGER = logging.getLogger(__name__)
METHOD_ID = "paired_chip_phasing_bench_v1"
PRODUCER = "D23_UPSTREAM_CONTRACT_FIXTURE"
CASES = ("BASELINE", "DENSITY_50", "MISSING_1PCT")


def load_configuration(path: Path) -> dict[str, Any]:
    """Valide le schéma fermé et la séparation des graines de rôle."""
    config = json.loads(path.read_text())
    schemas = Path(__file__).resolve().parents[3] / "calibration/schemas"
    validate_json_document(config, "phasing_bench_config.schema.json", schemas_dir=schemas)
    if len({config[key] for key in ("phasing_seed", "missingness_seed", "bootstrap_seed")}) != 3:
        raise ValueError("d23_seed_collision")
    return config


def observed_profile(panel: MeasurementPanel, case: str, config: dict[str, Any], target: dict[str, Any]) -> MeasurementPanel:
    """Sous-échantillonne en rang ou masque des GT diploïdes, cible protégée."""
    if case not in CASES:
        raise ValueError("d23_profile_unknown")
    chromosomes = np.array([int(row[0].removeprefix("chr")) for row in panel.variants])
    protected = np.array([int(row[0].removeprefix("chr")) == target["chromosome"] and int(row[1]) == target["position_bp"] for row in panel.variants])
    retained = np.ones(len(panel.variants), dtype=bool)
    if case == "DENSITY_50":
        retained[:] = False
        for chromosome in range(1, 23):
            retained[np.flatnonzero(chromosomes == chromosome)[::2]] = True
        retained |= protected
    calls = panel.alleles.copy()
    if case == "MISSING_1PCT":
        mask = np.random.default_rng(config["missingness_seed"]).random(calls.shape[:2]) < config["missingness_fraction"]
        mask[protected] = False
        calls[mask] = -1
    return MeasurementPanel(panel.samples, tuple(row for row, keep in zip(panel.variants, retained, strict=True) if keep), calls[retained])


def write_observation(path: Path, panel: MeasurementPanel) -> None:
    """VCF non phasé neuf ; aucune copie H1/H2 ni champ de vérité exporté."""
    with path.open("w", encoding="utf-8") as handle:
        handle.write('##fileformat=VCFv4.2\n##FORMAT=<ID=GT,Number=1,Type=String,Description="Observed genotype">\n')
        for chromosome in range(1, 23):
            handle.write(f"##contig=<ID=chr{chromosome}>\n")
        handle.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t" + "\t".join(panel.samples) + "\n")
        for row, calls in zip(panel.variants, panel.alleles, strict=True):
            values = ["./." if pair[0] < 0 else "/".join(map(str, sorted(pair))) for pair in calls]
            handle.write("\t".join((*row[:8], "GT", *values)) + "\n")


def _command(command: Sequence[str], root: Path, timeout: int) -> int:
    """Journalise commande et code ; les sorties sensibles restent privées."""
    ledger = root / "commands.jsonl"
    with (root / "commands.log").open("ab") as log:
        result = subprocess.run(list(command), stdout=log, stderr=subprocess.STDOUT, timeout=timeout, check=False)
    with ledger.open("a", encoding="utf-8") as handle:
        handle.write(json.dumps(dict(command=list(command), returncode=result.returncode, timeout_seconds=timeout)) + "\n")
    return result.returncode


def _descriptor(prefix: Path, scope: str, chromosome: int | None) -> Path:
    fam = [row.split() for row in prefix.with_suffix(".fam").read_text().splitlines()]
    bim = [row.split() for row in prefix.with_suffix(".bim").read_text().splitlines()]
    result = dict(schema_version="1.0.0", dataset_id=prefix.name, assembly="GRCh38", scope=scope, target_chromosome=chromosome,
                  sample_set_id="d23_fixture_samples", variant_set_id="d23_fixture_" + sha256_file(prefix.with_suffix(".bim"))[:16],
                  sample_count=len(fam), variant_count=len(bim), marker_mode="intersection", source_format="PLINK_FINAL_QC" if chromosome else "PLINK_PRELIMINARY_QC",
                  allele_orientation="forward_strand_calls", files={suffix: dict(name=prefix.with_suffix("." + suffix).name, sha256=sha256_file(prefix.with_suffix("." + suffix))) for suffix in ("bed", "bim", "fam")})
    validate_json_document(result, "plink_dataset.schema.json")
    output = prefix.with_suffix(".dataset.json")
    write_json(output, result)
    return output


def fixture_inputs(root: Path, preparation: Path, panel: MeasurementPanel, config: dict[str, Any], target: dict[str, Any]) -> list[dict[str, Any]]:
    """Convertit réellement puis déclare les contrats QC comme fixtures, pas PASS amont."""
    fixture = root / "stages" / PRODUCER
    fixture.mkdir(parents=True)
    write_observation(fixture / "observed.vcf", panel)
    timeout = config["command_timeout_seconds"]
    prefix = fixture / "genomewide_pre_qc"
    if _command(["plink", "--vcf", str(fixture / "observed.vcf"), "--double-id", "--keep-allele-order", "--make-bed", "--out", str(prefix)], root, timeout):
        raise RuntimeError("d23_fixture_conversion_failed")
    rows = raw_table(preparation / "blinded_inputs/metadata/samples.tsv")
    by_sample = {row["IID"]: row for row in rows}
    # Les lignes BED suivent toujours le FAM produit par PLINK, jamais un rang supposé.
    fam = [line.split() for line in prefix.with_suffix(".fam").read_text().splitlines()]
    with prefix.with_suffix(".fam").open("w") as handle:
        for line in fam:
            row = by_sample[line[1]]
            handle.write(" ".join((row["FID"], row["IID"], row["PID"], row["MID"], "0", "-9")) + "\n")
    write_table(fixture / "keep.txt", ("FID", "IID"), [])
    (fixture / "keep.txt").write_text("".join(f'{by_sample[line[1]]["FID"]}\t{line[1]}\n' for line in fam))
    target_prefix = fixture / "target_chromosome_all_qc"
    if _command(["plink", "--bfile", str(prefix), "--chr", str(target["chromosome"]), "--keep-allele-order", "--make-bed", "--out", str(target_prefix)], root, timeout):
        raise RuntimeError("d23_fixture_target_extract_failed")
    definitions = []
    descriptors = {}
    for base, scope, chromosome in ((prefix, "autosomal_genomewide", None), (target_prefix, "target_chromosome", target["chromosome"])):
        descriptor = _descriptor(base, scope, chromosome)
        descriptors[base.name] = json.loads(descriptor.read_text())
        definitions.extend((base.name + "_" + suffix, base.with_suffix("." + suffix), None) for suffix in ("bed", "bim", "fam"))
        definitions.append((base.name + "_dataset", descriptor, "plink_dataset.schema.json"))
    registry = fixture / "samples_master.tsv"
    shutil.copyfile(preparation / "blinded_inputs/metadata/samples.tsv", registry)
    genotypes = raw_table(preparation / "blinded_inputs/target_genotypes.tsv")
    audit = [dict(**row, FID=by_sample[row["SAMPLE_ID"]]["FID"], IID=by_sample[row["SAMPLE_ID"]]["IID"], STATUS="ACCEPTED") for row in genotypes]
    audit_columns = ("SAMPLE_ID", "FID", "IID", "GENOTYPE", "GENOTYPE_SOURCE", "MEASUREMENT_METHOD", "STATUS")
    write_table(fixture / "target_genotype_audit.tsv", audit_columns, audit)
    validate_tsv_table(fixture / "target_genotype_audit.tsv", "target_genotype_audit.schema.json")
    validate_tsv_table(registry, "samples_master.schema.json")
    qc = dict(DATASET_ID="target_chromosome_all_qc", COHORT_ID="target_chromosome_all_qc", VARIANT_ID=target["project_variant_id"], CHROMOSOME=str(target["chromosome"]),
              POSITION_BP=str(target["position_bp"]), A1=target["alt"], A2=target["ref"], N_MISS="0", N_GENO=str(len(rows)), F_MISS="0", MAF="", HWE_P="",
              HWE_STATUS="NOT_EVALUATED_POLICY", QC_STATUS="RETAINED_EXCEPTION", FILTER_CODES="D23_FIXTURE", IS_TARGET_VARIANT="true", FILTER_EXCEPTION="true", PRESENT_AFTER_QC="true")
    write_table(fixture / "variant_qc_final.tsv", tuple(qc), [qc])
    validate_tsv_table(fixture / "variant_qc_final.tsv", "variant_qc_final.schema.json")
    definitions += [("samples_master", registry, "samples_master.schema.json"), ("cohort_keep_target_chromosome_all_qc", fixture / "keep.txt", None),
                    ("target_genotype_audit", fixture / "target_genotype_audit.tsv", "target_genotype_audit.schema.json"),
                    ("variant_qc_final", fixture / "variant_qc_final.tsv", "variant_qc_final.schema.json")]
    definitions += [("config_input_" + name, Path(value), None) for name, value in load_pipeline_config(root / "config.resolved.yaml")["inputs"].items()
                    if name in {"target_variant_metadata", "reference_panel_catalog", "ancestry_reference_catalog", "genetic_map_catalog", "phasing_reference_selection"}]
    signature = sha256_file(fixture / "observed.vcf")
    artifacts = []
    for name, path, schema in definitions:
        context = next((descriptor for prefix, descriptor in descriptors.items() if name.startswith(prefix + "_")), {})
        artifacts.append(build_file_artifact(path, str(path), name, name, "application/json" if path.suffix == ".json" else "text/plain", PRODUCER, signature,
                                             schema_name=schema, schema_version="1.0.0" if schema else None, assembly="GRCh38",
                                             sample_set_id=context.get("sample_set_id"), variant_set_id=context.get("variant_set_id")))
    return artifacts


def execute_component(root: Path, stage: str, module: str, parameters: dict[str, Any], artifacts: list[dict[str, Any]], timeout: int) -> dict[str, Any]:
    """Exécute le point d'entrée de production sans simuler les outils."""
    directory = component_inputs(root, stage, module, stage + "_" + module, parameters, artifacts)
    LOGGER.info("D23 composant %s démarré", stage)
    code = _command([sys.executable, "-m", "effet_fondateur.stages." + module, "--stage-inputs", str(directory / "stage_inputs.json"), "--output-dir", str(directory)], root, timeout)
    outputs = checked_stage(root, directory.relative_to(root).as_posix()) if code == 0 else {}
    LOGGER.info("D23 composant %s terminé, code %s", stage, code)
    return dict(stage_id=stage, directory=directory.relative_to(root).as_posix(), returncode=code, outputs=list(outputs.values()))


def prepare_reference_cache(source_dir: Path, cache_root: Path, catalog: Path, chromosome: int, timeout: int) -> None:
    """Copie la vraie source publique ; l'API vérifie les empreintes officielles."""
    resolved = resolve_reference_panel(panel_id="1kg_3202_high_coverage_20220422", assembly="GRCh38", chromosome=chromosome, catalog_path=catalog)
    def downloader(url: str, destination: Path, timeout_seconds: int, chunk_size: int) -> None:
        name = Path(urlparse(url).path).name
        local = source_dir / name
        if local.is_file() and not local.is_symlink():
            shutil.copyfile(local, destination)
        else:
            if name in {resolved.vcf.filename, resolved.index.filename}:
                raise ValueError("d23_declared_local_reference_missing")
            with urlopen(url, timeout=timeout_seconds) as response, destination.open("xb") as handle:
                shutil.copyfileobj(response, handle, chunk_size)
    cache_reference_panel(resolved, cache_root, downloader=downloader, download_timeout_seconds=timeout)


def run_bench(*, preparation_root: Path, configuration_path: Path, workspace_root: Path, run_id: str, reference_source_dir: Path, shapeit_bin_dir: Path) -> dict[str, Any]:
    """Exécute trois profils fermés ; un échec scientifique n'est jamais masqué."""
    destination = new_workspace_run(workspace_root, run_id).absolute()
    preparation = preparation_root.absolute()
    config = load_configuration(configuration_path)
    if preparation.parent != destination.parent or reference_source_dir.is_symlink() or not reference_source_dir.is_dir():
        raise ValueError("d23_source_scope_invalid")
    verify_preparation(preparation)
    if sha256_file(preparation / "preparation_manifest.json") != config["preparation_manifest_sha256"]:
        raise ValueError("d23_preparation_anchor_mismatch")
    for tool in ("plink", "bcftools"):
        if shutil.which(tool) is None:
            raise ValueError("d23_required_tool_missing")
    destination.mkdir(mode=0o700)
    write_json(destination / "configuration.json", config)
    target = yaml.safe_load((preparation / "blinded_inputs/target_metadata.yaml").read_text())
    base = yaml.safe_load((preparation / "blinded_inputs/pipeline.yaml").read_text())
    base["tools"].update(plink="plink", bcftools="bcftools", phasing_adapter=dict(adapter_id="shapeit5_phase_common_rare_v1", phase_common_command=str(shapeit_bin_dir / "SHAPEIT5_phase_common"), phase_rare_command=str(shapeit_bin_dir / "SHAPEIT5_phase_rare"), expected_version="5.1.1"))
    load_pipeline_config(preparation / "blinded_inputs/pipeline.yaml")
    maps = destination / "genetic_maps"
    shutil.copytree(preparation / "reference_assets/genetic_maps", maps)
    cache = destination / "reference_cache"
    prepare_reference_cache(reference_source_dir.absolute(), cache, Path(base["inputs"]["reference_panel_catalog"]), target["chromosome"], config["command_timeout_seconds"])
    with gzip.open(preparation / "blinded_inputs/observed.vcf.gz", "rt") as handle:
        observed = read_panel(handle, phased=False)
    receipts = []
    for case in config["profiles"]:
        LOGGER.info("D23 profil %s", case)
        root = destination / case
        root.mkdir()
        write_json(root / "manifest.json", dict(run_id="d23_" + case.lower(), scope="COMPONENT_BENCH_NOT_FULL_PIPELINE"))
        (root / "config.resolved.yaml").write_text(yaml.safe_dump(base, sort_keys=False))
        load_pipeline_config(root / "config.resolved.yaml")
        profile = observed_profile(observed, case, config, target)
        inputs = fixture_inputs(root, preparation, profile, config, target)
        steps = []
        timeout = config["command_timeout_seconds"]
        steps.append(execute_component(root, "10A", "prepare_autosomal_phasing_panel", dict(minimum_total_variants=config["minimum_total_variants"], minimum_variants_per_autosome=config["minimum_variants_per_autosome"]), inputs, timeout))
        params11 = dict(window_left_bp=config["target_flank_bp"], window_right_bp=config["target_flank_bp"], genetic_map_id="shapeit4_hapmap_grch38_7a0cab7", genetic_map_cache_dir=str(maps), genetic_map_cache_offline=True)
        steps.append(execute_component(root, "11", "prepare_target_region", params11, [row for row in inputs if row["artifact_id"] != "config_input_phasing_reference_selection"], timeout))
        common = dict(reference_panel_id="1kg_3202_high_coverage_20220422", reference_cache_dir=str(cache), require_reference_sample_selection=True,
                      shapeit5_threads=1, shapeit5_seed=config["phasing_seed"])
        if steps[-1]["returncode"] == 0:
            params12 = dict(**common, reference_cache_offline=True, minimum_common_variants=config["minimum_variants_per_autosome"], shapeit5_input_timeout_seconds=timeout, shapeit5_execution_timeout_seconds=timeout)
            steps.append(execute_component(root, "12", "phase_target_region", params12, inputs + steps[-1]["outputs"], timeout))
        if steps[0]["returncode"] == 0:
            params12a = dict(**common, genetic_map_cache_dir=str(maps), genetic_map_id="shapeit4_hapmap_grch38_7a0cab7", reference_source_dir=str(reference_source_dir.absolute()),
                             reference_extract_workers=2, tool_timeout_seconds=timeout, reference_extract_timeout_seconds=timeout, minimum_matched_variants_per_autosome=config["minimum_variants_per_autosome"])
            steps.append(execute_component(root, "12A", "phase_autosomal_panel", params12a, inputs + steps[0]["outputs"], timeout * 22))
        receipts.append(dict(profile=case, observed_markers=len(profile.variants), missing_genotypes=int((profile.alleles[:, :, 0] < 0).sum()), stages=steps))
        write_json(destination / "progress.json", dict(profiles=receipts, state="RUNNING"))
    if sha256_file(preparation / "preparation_manifest.json") != config["preparation_manifest_sha256"]:
        raise ValueError("d23_preparation_changed_during_run")
    receipt = dict(schema_version="1.0.0", method_id=METHOD_ID, closure_state="CLOSED", preparation_manifest_sha256=config["preparation_manifest_sha256"],
                   configuration_sha256=sha256_file(configuration_path), scope="ISOLATED_COMPONENTS_10A_11_12_12A_NOT_04_09", fixture_producer=PRODUCER,
                   implementation_sha256=sha256_file(Path(__file__)), profiles=receipts, artifacts=file_inventory(destination, ("run_index.json",)), calibrated=False, founder_effect_proven=False)
    write_json(destination / "run_index.json", receipt)
    for path in destination.rglob("*"):
        # Les caches de production immuables gardent leurs permissions restrictives.
        if path.is_file() and path.stat().st_mode & 0o222:
            path.chmod(0o600)
    return receipt


def main() -> int:
    """CLI explicite de l'intégration D23, sans action sur l'étude réelle."""
    logging.basicConfig(level=logging.INFO, format="%(asctime)s %(message)s")
    parser = argparse.ArgumentParser(description="Banc borné D23 de phasage")
    for name in ("preparation-root", "configuration-path", "workspace-root", "reference-source-dir", "shapeit-bin-dir"):
        parser.add_argument("--" + name, type=Path, required=True)
    parser.add_argument("--run-id", required=True)
    result = run_bench(**vars(parser.parse_args()))
    print(dict(state=result["closure_state"], profiles=len(result["profiles"])))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
