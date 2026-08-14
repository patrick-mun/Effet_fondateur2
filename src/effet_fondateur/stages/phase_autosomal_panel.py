"""Phase le même panel d'étude sur chacun des 22 autosomes avec SHAPEIT5."""

from __future__ import annotations

import argparse
import bisect
import csv
import hashlib
import json
import math
import shutil
import subprocess
import sys
import tempfile
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path, PurePosixPath
from time import monotonic
from typing import Any, Callable, Sequence

from effet_fondateur.ancestry import cache_ancestry_metadata, cache_reference_extract, load_reference_samples
from effet_fondateur.audit import atomic_write_json, read_json, sha256_file
from effet_fondateur.contracts import DocumentValidationError, build_file_artifact, load_pipeline_config, validate_json_document
from effet_fondateur.orchestrator.state import utc_now
from effet_fondateur.phasing.shapeit5 import parse_shapeit5_adapter_config, probe_shapeit5
from effet_fondateur.references.genetic_maps import ensure_genetic_map_cached, resolve_genetic_map
from effet_fondateur.stages.analyze_reference_ancestry import _extract_remote_reference, _reference_panel


class AutosomalPhasingInputError(ValueError):
    """Signale une entrée/configuration invalide."""


class AutosomalPhasingExternalError(RuntimeError):
    """Signale un échec d'un outil externe ou du réseau de référence."""


class AutosomalPhasingBlockError(RuntimeError):
    """Signale un résultat de phasage scientifiquement incomplet."""


def _artifact(stage_inputs: dict[str, Any], artifact_id: str) -> dict[str, Any]:
    matches = [item for item in stage_inputs["artifacts"] if item["artifact_id"] == artifact_id]
    if len(matches) != 1:
        raise AutosomalPhasingInputError(f"{artifact_id}_missing_or_ambiguous")
    return matches[0]


def _path(artifact: dict[str, Any], run_dir: Path) -> Path:
    candidate = Path(artifact["path"])
    path = run_dir / candidate if PurePosixPath(artifact["path"]).parts[:1] == ("stages",) else candidate
    if not path.is_absolute():
        path = Path.cwd() / path
    if path.is_symlink() or not path.is_file() or sha256_file(path) != artifact["sha256"]:
        raise AutosomalPhasingInputError("declared_input_modified")
    return path


def _integer(parameters: dict[str, Any], name: str, default: int) -> int:
    value = parameters.get(name, default)
    if isinstance(value, bool) or not isinstance(value, int) or value <= 0:
        raise AutosomalPhasingInputError(f"invalid_parameter:{name}")
    return value


def _run(command: Sequence[str], timeout: int, code: str) -> subprocess.CompletedProcess[str]:
    try:
        result = subprocess.run(list(command), capture_output=True, text=True, check=False, timeout=timeout)
    except (OSError, subprocess.TimeoutExpired) as error:
        raise AutosomalPhasingExternalError(code) from error
    if result.returncode != 0:
        detail = (result.stderr or result.stdout or "")[-600:].replace("\n", " ")
        raise AutosomalPhasingExternalError(f"{code}:{result.returncode}:{detail}")
    return result


def _file(path: Path, root: Path) -> dict[str, Any]:
    return {"path": path.relative_to(root).as_posix(), "sha256": sha256_file(path), "size_bytes": path.stat().st_size}


def _write_lines(path: Path, lines: Sequence[str]) -> None:
    path.write_text("".join(f"{line}\n" for line in lines), encoding="utf-8")


def _md5_file(path: Path) -> str:
    digest = hashlib.md5(usedforsecurity=False)
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _validated_local_reference_source(
    source_dir: Path,
    filename: str,
    expected_vcf_md5: str,
    expected_index_md5: str,
) -> Path:
    """Retourne un VCF local uniquement s'il correspond au manifeste officiel."""

    vcf_path = source_dir / filename
    index_path = source_dir / f"{filename}.tbi"
    if not vcf_path.is_file() or not index_path.is_file():
        raise AutosomalPhasingInputError(f"local_reference_source_missing:{filename}")
    if _md5_file(vcf_path) != expected_vcf_md5 or _md5_file(index_path) != expected_index_md5:
        raise AutosomalPhasingInputError(f"local_reference_source_checksum_mismatch:{filename}")
    return vcf_path.resolve()


def _study_variants(bim: Path) -> dict[int, list[tuple[str, int, str, str]]]:
    result = {chromosome: [] for chromosome in range(1, 23)}
    for line in bim.read_text(encoding="utf-8").splitlines():
        fields = line.split()
        if len(fields) != 6:
            raise AutosomalPhasingInputError("autosomal_panel_bim_malformed")
        chromosome, variant_id, position, allele_1, allele_2 = int(fields[0].removeprefix("chr")), fields[1], int(fields[3]), fields[4].upper(), fields[5].upper()
        if chromosome not in result or position <= 0 or allele_1 not in "ACGT" or allele_2 not in "ACGT" or allele_1 == allele_2:
            raise AutosomalPhasingInputError("autosomal_panel_variant_invalid")
        result[chromosome].append((variant_id, position, allele_1, allele_2))
    if any(not rows for rows in result.values()):
        raise AutosomalPhasingBlockError("autosomal_panel_chromosome_empty")
    return result


def _reference_alleles(bcftools: str, vcf: Path, chromosome: int, timeout: int) -> dict[int, tuple[str, str]]:
    query = _run([bcftools, "query", "-f", "%CHROM\t%POS\t%REF\t%ALT\n", str(vcf)], timeout, f"reference_query_failed:chr{chromosome}")
    by_position: dict[int, list[tuple[str, str]]] = {}
    for line in query.stdout.splitlines():
        chrom, raw_position, ref, alt = line.split("\t")
        if int(chrom.removeprefix("chr")) != chromosome or "," in alt:
            continue
        by_position.setdefault(int(raw_position), []).append((ref.upper(), alt.upper()))
    return {position: rows[0] for position, rows in by_position.items() if len(rows) == 1}


def _map_points(path: Path, chromosome: int) -> tuple[list[int], list[float]]:
    with path.open(encoding="utf-8", newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    positions = [int(row["POSITION_BP"]) for row in rows if int(row["CHROMOSOME"]) == chromosome]
    cms = [float(row["POSITION_CM"]) for row in rows if int(row["CHROMOSOME"]) == chromosome]
    if len(positions) < 2 or any(b <= a for a, b in zip(positions, positions[1:])) or any(b < a for a, b in zip(cms, cms[1:])):
        raise AutosomalPhasingBlockError(f"genetic_map_invalid:chr{chromosome}")
    return positions, cms


def _interpolate(position: int, positions: list[int], cms: list[float]) -> float:
    index = bisect.bisect_left(positions, position)
    if index == 0:
        return cms[0]
    if index == len(positions):
        return cms[-1]
    left_bp, right_bp = positions[index - 1], positions[index]
    left_cm, right_cm = cms[index - 1], cms[index]
    return left_cm + (right_cm - left_cm) * (position - left_bp) / (right_bp - left_bp)


def _validated_phased_variants(bcftools: str, bcf: Path, chromosome: int, expected_samples: list[str], timeout: int) -> list[tuple[str, int]]:
    samples = _run([bcftools, "query", "--list-samples", str(bcf)], timeout, f"phased_samples_failed:chr{chromosome}").stdout.splitlines()
    if samples != expected_samples:
        raise AutosomalPhasingBlockError(f"phased_sample_order_mismatch:chr{chromosome}")
    bad = _run([bcftools, "query", "-i", 'GT="mis" || GT~"/"', "-f", "%CHROM\t%POS\n", str(bcf)], timeout, f"phased_gt_check_failed:chr{chromosome}")
    if bad.stdout.strip():
        raise AutosomalPhasingBlockError(f"unphased_or_missing_genotype:chr{chromosome}")
    query = _run([bcftools, "query", "-f", "%ID\t%POS\n", str(bcf)], timeout, f"phased_variant_query_failed:chr{chromosome}")
    rows = [(fields[0], int(fields[1])) for fields in (line.split("\t") for line in query.stdout.splitlines())]
    if len(rows) < 2 or len({variant_id for variant_id, _ in rows}) != len(rows):
        raise AutosomalPhasingBlockError(f"phased_variants_invalid:chr{chromosome}")
    return rows


def execute(stage_inputs_path: Path, output_dir: Path) -> int:
    """Extrait la référence aux marqueurs observés puis phase chromosome par chromosome."""
    started_at, started_clock = utc_now(), monotonic()
    stage_inputs = read_json(stage_inputs_path)
    validate_json_document(stage_inputs, "stage_inputs.schema.json")
    parameters = stage_inputs["parameters"]
    reference_panel_id = parameters.get("reference_panel_id", "1kg_3202_high_coverage_20220422")
    map_id = parameters.get("genetic_map_id", "shapeit4_hapmap_grch38_7a0cab7")
    if not isinstance(reference_panel_id, str) or not isinstance(map_id, str):
        raise AutosomalPhasingInputError("reference_or_map_id_invalid")
    tool_timeout = _integer(parameters, "tool_timeout_seconds", 7200)
    extract_timeout = _integer(parameters, "reference_extract_timeout_seconds", 7200)
    extract_workers = _integer(parameters, "reference_extract_workers", 4)
    extract_chunk = _integer(parameters, "reference_extract_chunk_variants", 1000)
    threads = _integer(parameters, "shapeit5_threads", 2)
    seed = _integer(parameters, "shapeit5_seed", 15052011)
    minimum_matched = _integer(parameters, "minimum_matched_variants_per_autosome", 100)
    cache_root = Path(parameters.get("reference_cache_dir", "data/cache/references"))
    if not cache_root.is_absolute():
        cache_root = Path.cwd() / cache_root
    source_dir_value = parameters.get("reference_source_dir")
    if source_dir_value is not None and (not isinstance(source_dir_value, str) or not source_dir_value):
        raise AutosomalPhasingInputError("invalid_parameter:reference_source_dir")
    source_dir = Path(source_dir_value) if source_dir_value is not None else None
    if source_dir is not None and not source_dir.is_absolute():
        source_dir = Path.cwd() / source_dir
    if source_dir is not None and not source_dir.is_dir():
        raise AutosomalPhasingInputError("local_reference_source_dir_missing")
    offline = parameters.get("reference_cache_offline", False)
    if not isinstance(offline, bool):
        raise AutosomalPhasingInputError("invalid_parameter:reference_cache_offline")
    run_dir = output_dir.parent.parent
    config = load_pipeline_config(run_dir / "config.resolved.yaml")
    ids = (
        "config_input_reference_panel_catalog", "config_input_ancestry_reference_catalog", "config_input_genetic_map_catalog",
        "samples_master", "autosomal_phasing_panel_bed", "autosomal_phasing_panel_bim",
        "autosomal_phasing_panel_fam", "autosomal_phasing_panel_dataset",
    )
    artifacts = {name: _artifact(stage_inputs, name) for name in ids}
    paths = {name: _path(value, run_dir) for name, value in artifacts.items()}
    descriptor = read_json(paths["autosomal_phasing_panel_dataset"])
    validate_json_document(descriptor, "plink_dataset.schema.json")
    if descriptor["assembly"] != "GRCh38" or descriptor["sample_count"] < 1:
        raise AutosomalPhasingInputError("autosomal_panel_descriptor_invalid")
    variants = _study_variants(paths["autosomal_phasing_panel_bim"])
    fam_rows = [line.split() for line in paths["autosomal_phasing_panel_fam"].read_text(encoding="utf-8").splitlines() if line]
    expected_samples = [row[1] for row in fam_rows]
    if len(expected_samples) != descriptor["sample_count"] or len(set(expected_samples)) != len(expected_samples):
        raise AutosomalPhasingInputError("autosomal_panel_samples_invalid")
    plink, bcftools = shutil.which(config["tools"]["plink"]), shutil.which(config["tools"]["bcftools"])
    if plink is None or bcftools is None:
        raise AutosomalPhasingExternalError("plink_or_bcftools_not_found")
    adapter = parse_shapeit5_adapter_config(config["tools"]["phasing_adapter"])
    shapeit5 = probe_shapeit5(adapter)
    metadata_cache = cache_ancestry_metadata(catalog_path=paths["config_input_ancestry_reference_catalog"], cache_root=cache_root, offline=offline, timeout_seconds=extract_timeout)
    reference_samples = [sample.sample_id for sample in load_reference_samples(metadata_cache)]
    panel = _reference_panel(paths["config_input_reference_panel_catalog"], reference_panel_id)
    if panel["assembly"] != descriptor["assembly"]:
        raise AutosomalPhasingInputError("reference_assembly_mismatch")
    sources = {row["chromosome"]: row for row in panel["chromosomes"]}
    resolved_map = resolve_genetic_map(paths["config_input_genetic_map_catalog"], map_id, descriptor["assembly"])
    output_dir.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(prefix=".autosomal_phasing_", dir=output_dir) as temporary_name:
        temporary = Path(temporary_name)
        sample_file = temporary / "reference.samples.txt"
        _write_lines(sample_file, reference_samples)
        position_files: dict[int, Path] = {}
        for chromosome, rows in variants.items():
            path = temporary / f"chr{chromosome}.positions.tsv"
            _write_lines(path, [f"chr{chromosome}\t{position}\t{position}" for _, position, _, _ in rows])
            position_files[chromosome] = path

        def cache_one(chromosome: int):
            source = sources[chromosome]
            filename = panel["vcf_filename_template"].format(chromosome=chromosome)
            source_url = (
                str(_validated_local_reference_source(
                    source_dir,
                    filename,
                    source["vcf_md5"],
                    source["index_md5"],
                ))
                if source_dir is not None
                else f"{panel['base_url']}/{filename}"
            )

            def extractor(url: str, positions_path: Path, samples_path: Path, vcf: Path, index: Path, timeout: int) -> None:
                def runner(command: list[str], command_timeout: int, code: str):
                    return _run(command, command_timeout, code)
                _extract_remote_reference(bcftools=bcftools, source_url=url, positions=positions_path, samples=samples_path, vcf=vcf, index=index, timeout=timeout, chunk_size=extract_chunk, runner=runner)

            return chromosome, cache_reference_extract(
                cache_root=cache_root, panel_id=panel["panel_id"], assembly=panel["assembly"], chromosome=chromosome,
                source_url=source_url, source_vcf_md5=source["vcf_md5"], source_index_md5=source["index_md5"],
                positions_path=position_files[chromosome], samples_path=sample_file, offline=offline,
                timeout_seconds=extract_timeout, extractor=extractor,
            )

        # Consommer les résultats par ordre d'achèvement rend une panne réseau
        # immédiatement visible. `executor.map()` attendrait sinon le premier
        # chromosome dans l'ordre tout en lançant inutilement les 22 tâches.
        executor = ThreadPoolExecutor(max_workers=min(extract_workers, 22))
        futures = {
            executor.submit(cache_one, chromosome): chromosome
            for chromosome in range(1, 23)
        }
        cached_references: dict[int, Any] = {}
        try:
            for future in as_completed(futures):
                chromosome, cached = future.result()
                cached_references[chromosome] = cached
        except Exception:
            for future in futures:
                future.cancel()
            executor.shutdown(wait=True, cancel_futures=True)
            raise
        else:
            executor.shutdown(wait=True)
        chromosome_records: list[dict[str, Any]] = []
        output_artifacts: list[dict[str, Any]] = []
        producer = f"{stage_inputs['stage_id']}_{stage_inputs['stage_name']}"
        published = PurePosixPath(stage_inputs["published_output_dir"])
        total_phased = 0
        for chromosome in range(1, 23):
            chromosome_dir = output_dir / f"chr{chromosome}"
            chromosome_dir.mkdir()
            reference = cached_references[chromosome]
            reference_by_position = _reference_alleles(bcftools, reference.vcf_path, chromosome, tool_timeout)
            selected: list[tuple[str, int, str]] = []
            for variant_id, position, allele_1, allele_2 in variants[chromosome]:
                reference_pair = reference_by_position.get(position)
                if reference_pair is not None and {allele_1, allele_2} == set(reference_pair):
                    selected.append((variant_id, position, reference_pair[0]))
            if len(selected) < minimum_matched:
                raise AutosomalPhasingBlockError(f"reference_matched_density_insufficient:chr{chromosome}")
            extract_ids, reference_alleles = temporary / f"chr{chromosome}.ids", temporary / f"chr{chromosome}.a2.tsv"
            _write_lines(extract_ids, [row[0] for row in selected])
            _write_lines(reference_alleles, [f"{row[0]}\t{row[2]}" for row in selected])
            raw_prefix = temporary / f"chr{chromosome}.study"
            _run([plink, "--bfile", str(paths["autosomal_phasing_panel_bed"].with_suffix("")), "--chr", str(chromosome), "--extract", str(extract_ids), "--a2-allele", str(reference_alleles), "2", "1", "--recode", "vcf-iid", "bgz", "--out", str(raw_prefix)], tool_timeout, f"plink_vcf_export_failed:chr{chromosome}")
            rename = temporary / f"chr{chromosome}.rename.tsv"
            _write_lines(rename, [f"{chromosome}\tchr{chromosome}"])
            study_vcf = chromosome_dir / "study.harmonized.vcf.gz"
            _run([bcftools, "annotate", "--rename-chrs", str(rename), "-Oz", "-o", str(study_vcf), str(raw_prefix.with_suffix(".vcf.gz"))], tool_timeout, f"study_contig_normalization_failed:chr{chromosome}")
            _run([bcftools, "index", "--tbi", str(study_vcf)], tool_timeout, f"study_index_failed:chr{chromosome}")
            study_query = _run([bcftools, "query", "-f", "%POS\t%REF\t%ALT\n", str(study_vcf)], tool_timeout, f"study_allele_query_failed:chr{chromosome}")
            if any(reference_by_position.get(int(fields[0])) != (fields[1], fields[2]) for fields in (line.split("\t") for line in study_query.stdout.splitlines())):
                raise AutosomalPhasingBlockError(f"study_reference_allele_mismatch:chr{chromosome}")
            cached_map = ensure_genetic_map_cached(resolved=resolved_map, chromosome=chromosome, cache_root=Path(parameters.get("genetic_map_cache_dir", "data/cache/references/genetic_maps")), offline=offline)
            source_map = cached_map.map_path.parent.parent / "source" / resolved_map.member_template.format(chromosome=chromosome)
            phased_bcf, log_path = chromosome_dir / "study.phased.bcf", chromosome_dir / "shapeit5.phase_common.log"
            phase = _run([shapeit5.phase_common_path, "--input", str(study_vcf), "--reference", str(reference.vcf_path), "--map", str(source_map), "--region", f"chr{chromosome}", "--output", str(phased_bcf), "--output-format", "bcf", "--log", str(log_path), "--thread", str(threads), "--seed", str(seed + chromosome)], tool_timeout, f"shapeit5_phase_common_failed:chr{chromosome}")
            if not log_path.is_file():
                log_path.write_text(phase.stdout + phase.stderr, encoding="utf-8")
            phased_index = Path(f"{phased_bcf}.csi")
            if not phased_index.is_file():
                _run([bcftools, "index", "--csi", str(phased_bcf)], tool_timeout, f"phased_index_failed:chr{chromosome}")
            phased_rows = _validated_phased_variants(bcftools, phased_bcf, chromosome, expected_samples, tool_timeout)
            if len(phased_rows) != len(selected):
                raise AutosomalPhasingBlockError(f"phased_variant_count_mismatch:chr{chromosome}")
            map_positions, map_cms = _map_points(cached_map.map_path, chromosome)
            ibd_map = chromosome_dir / "study.phased.plink.map"
            with ibd_map.open("w", encoding="utf-8", newline="") as handle:
                previous_cm = -math.inf
                for variant_id, position in phased_rows:
                    cm = _interpolate(position, map_positions, map_cms)
                    if cm < previous_cm:
                        raise AutosomalPhasingBlockError(f"interpolated_map_not_monotonic:chr{chromosome}")
                    # Le libellé doit être identique au contig VCF (`chrN`) :
                    # les deux moteurs refusent silencieusement ou explicitement
                    # une carte `N` associée à un VCF `chrN`.
                    handle.write(f"chr{chromosome}\t{variant_id}\t{cm:.12g}\t{position}\n")
                    previous_cm = cm
            record = {"chromosome": chromosome, "study_input_variants": len(variants[chromosome]), "reference_matched_variants": len(selected), "phased_variants": len(phased_rows), "reference_cache_status": reference.status, "bcf": _file(phased_bcf, output_dir), "index": _file(phased_index, output_dir), "ibd_map": _file(ibd_map, output_dir), "log": _file(log_path, output_dir)}
            chromosome_records.append(record)
            total_phased += len(phased_rows)
            for kind, path, media in (("bcf", phased_bcf, "application/octet-stream"), ("index", phased_index, "application/octet-stream"), ("ibd_map", ibd_map, "text/plain"), ("log", log_path, "text/plain")):
                output_artifacts.append(build_file_artifact(physical_path=path, published_path=str(published / path.relative_to(output_dir).as_posix()), artifact_id=f"autosomal_phased_chr{chromosome}_{kind}", artifact_type=f"autosomal_phased_{kind}", media_type=media, producer_stage=producer, producer_signature=stage_inputs["signature"], assembly=descriptor["assembly"], sample_set_id=descriptor["sample_set_id"], variant_set_id=descriptor["variant_set_id"], sensitivity="sensitive_genetic" if kind == "bcf" else "internal"))
    manifest = {"schema_version": "1.0.0", "method_id": "shapeit5_reference_common_autosomes_v1", "assembly": descriptor["assembly"], "sample_set_id": descriptor["sample_set_id"], "variant_set_id": descriptor["variant_set_id"], "sample_count": len(expected_samples), "variant_count": total_phased, "chromosomes": chromosome_records, "checks": {name: "PASS" for name in ("same_samples_all_chromosomes", "complete_genotypes", "all_genotypes_phased", "reference_alleles", "maps_monotonic", "all_22_autosomes")}}
    validate_json_document(manifest, "autosomal_phasing_manifest.schema.json")
    manifest_path = output_dir / "autosomal_phasing_manifest.json"
    atomic_write_json(manifest_path, manifest)
    output_artifacts.append(build_file_artifact(physical_path=manifest_path, published_path=str(published / manifest_path.name), artifact_id="autosomal_phasing_manifest", artifact_type="autosomal_phasing_manifest", media_type="application/json", producer_stage=producer, producer_signature=stage_inputs["signature"], schema_name="autosomal_phasing_manifest.schema.json", schema_version="1.0.0", assembly=descriptor["assembly"], sample_set_id=descriptor["sample_set_id"], variant_set_id=descriptor["variant_set_id"], sensitivity="internal"))
    stage_outputs = {"schema_version": "1.0.0", "run_id": stage_inputs["run_id"], "stage_id": stage_inputs["stage_id"], "stage_name": stage_inputs["stage_name"], "signature": stage_inputs["signature"], "artifacts": output_artifacts}
    validate_json_document(stage_outputs, "stage_outputs.schema.json")
    atomic_write_json(output_dir / "stage_outputs.json", stage_outputs)
    audit = {"schema_version": "1.0.0", "run_id": stage_inputs["run_id"], "stage_id": stage_inputs["stage_id"], "stage_name": stage_inputs["stage_name"], "method_id": manifest["method_id"], "signature": stage_inputs["signature"], "started_at": started_at, "completed_at": utc_now(), "duration_seconds": monotonic() - started_clock, "inputs": list(artifacts.values()), "outputs": output_artifacts, "parameters": parameters, "tools": [{"tool": "SHAPEIT5_phase_common", "configured": adapter.phase_common_command, "version": shapeit5.phase_common_version}, {"tool": "plink", "configured": config["tools"]["plink"], "version": None}, {"tool": "bcftools", "configured": config["tools"]["bcftools"], "version": None}], "counts": {"samples": len(expected_samples), "phased_variants": total_phased, "autosomes": 22, "reference_samples": len(reference_samples)}, "metrics": {"per_chromosome": [{"chromosome": row["chromosome"], "phased_variants": row["phased_variants"]} for row in chromosome_records]}, "exclusions": [{"code": "not_uniquely_reference_matched", "count": sum(len(variants[row["chromosome"]]) - row["phased_variants"] for row in chromosome_records)}], "warnings": [], "checks": [{"check": name, "status": status} for name, status in manifest["checks"].items()], "known_limits": ["Le phasage autosomal IBD porte sur les SNP bialléliques complets, polymorphes et appariés sans ambiguïté à la référence externe.", "Le variant cible est attribué séparément par le phasage régional avec génotypes explicites."], "expected_visualizations": [], "manual_validation_required": False}
    validate_json_document(audit, "stage_audit.schema.json")
    atomic_write_json(output_dir / "audit.json", audit)
    (output_dir / "checksums.sha256").write_text("".join(f"{item['sha256']}  {item['artifact_id']}\n" for item in output_artifacts), encoding="utf-8")
    return 0


def main(arguments: Sequence[str] | None = None) -> int:
    parser = argparse.ArgumentParser(prog="phase-autosomal-panel")
    parser.add_argument("--stage-inputs", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parsed = parser.parse_args(arguments)
    try:
        return execute(parsed.stage_inputs, parsed.output_dir)
    except AutosomalPhasingBlockError as error:
        sys.stderr.write(f"{error}\n"); return 4
    except AutosomalPhasingExternalError as error:
        sys.stderr.write(f"{error}\n"); return 3
    except (AutosomalPhasingInputError, DocumentValidationError, OSError, ValueError, json.JSONDecodeError) as error:
        sys.stderr.write(f"{error}\n"); return 2


if __name__ == "__main__":
    raise SystemExit(main())
