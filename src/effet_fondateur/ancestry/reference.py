"""Cache immuable des métadonnées de population officielles 1000 Genomes."""

from __future__ import annotations

import fcntl
import hashlib
import json
import os
import shutil
import tempfile
from dataclasses import dataclass
from pathlib import Path
from typing import Callable
from urllib.parse import urlparse
from urllib.request import Request, urlopen

from effet_fondateur.audit import atomic_write_json, sha256_file
from effet_fondateur.contracts import validate_json_document
from effet_fondateur.orchestrator.state import utc_now


Downloader = Callable[[str, Path, int], None]


class AncestryReferenceError(RuntimeError):
    """Signale un catalogue, cache ou fichier de métadonnées incohérent."""


@dataclass(frozen=True)
class CachedAncestryMetadata:
    """Décrit les métadonnées publiques vérifiées et leur statut de cache."""

    status: str
    entry_dir: Path
    population_metadata_path: Path
    unrelated_index_path: Path
    manifest_path: Path


@dataclass(frozen=True)
class ReferenceSample:
    """Métadonnées non sensibles d'un échantillon public de référence."""

    sample_id: str
    population: str
    superpopulation: str


def _download(url: str, destination: Path, timeout_seconds: int) -> None:
    parsed = urlparse(url)
    if parsed.scheme != "https" or parsed.hostname != "ftp.1000genomes.ebi.ac.uk":
        raise AncestryReferenceError("ancestry_metadata_url_not_approved")
    request = Request(url, headers={"User-Agent": "effet-fondateur/0.1"})
    try:
        with urlopen(request, timeout=timeout_seconds) as response:
            final = urlparse(response.geturl())
            if final.scheme != "https" or final.hostname != parsed.hostname:
                raise AncestryReferenceError("ancestry_metadata_redirect_not_approved")
            with destination.open("xb") as handle:
                while chunk := response.read(1024 * 1024):
                    handle.write(chunk)
                handle.flush()
                os.fsync(handle.fileno())
    except AncestryReferenceError:
        raise
    except OSError as error:
        raise AncestryReferenceError("ancestry_metadata_download_failed") from error


def _catalog(path: Path) -> dict[str, object]:
    try:
        document = json.loads(path.read_text(encoding="utf-8"))
        validate_json_document(document, "ancestry_reference_catalog.schema.json")
    except (OSError, ValueError) as error:
        raise AncestryReferenceError("ancestry_reference_catalog_invalid") from error
    return document


def _validate_file(path: Path, expected_sha256: str, code: str) -> None:
    if path.is_symlink() or not path.is_file() or sha256_file(path) != expected_sha256:
        raise AncestryReferenceError(code)


def cache_ancestry_metadata(
    *,
    catalog_path: Path,
    cache_root: Path,
    offline: bool,
    timeout_seconds: int = 300,
    downloader: Downloader = _download,
) -> CachedAncestryMetadata:
    """Télécharge une fois les métadonnées épinglées, puis vérifie chaque réutilisation."""
    catalog = _catalog(catalog_path)
    catalog_sha = sha256_file(catalog_path)
    cache_key = hashlib.sha256(
        json.dumps(
            {"catalog_sha256": catalog_sha, "catalog": catalog},
            sort_keys=True,
            separators=(",", ":"),
        ).encode("utf-8")
    ).hexdigest()
    entry_parent = cache_root / "ancestry_metadata" / str(catalog["catalog_id"])
    entry_parent.mkdir(parents=True, exist_ok=True)
    if entry_parent.is_symlink():
        raise AncestryReferenceError("ancestry_metadata_cache_parent_invalid")
    entry_dir = entry_parent / cache_key
    population_asset = catalog["population_metadata"]
    unrelated_asset = catalog["unrelated_sample_index"]
    assert isinstance(population_asset, dict) and isinstance(unrelated_asset, dict)
    population_name = Path(urlparse(str(population_asset["url"])).path).name
    unrelated_name = Path(urlparse(str(unrelated_asset["url"])).path).name
    manifest_path = entry_dir / "ancestry_metadata_cache_manifest.json"

    def result(status: str) -> CachedAncestryMetadata:
        population_path = entry_dir / population_name
        unrelated_path = entry_dir / unrelated_name
        if entry_dir.stat().st_mode & 0o222:
            raise AncestryReferenceError("ancestry_metadata_cache_writable")
        _validate_file(population_path, str(population_asset["sha256"]), "ancestry_population_metadata_corrupt")
        _validate_file(unrelated_path, str(unrelated_asset["sha256"]), "ancestry_unrelated_index_corrupt")
        if not manifest_path.is_file() or manifest_path.is_symlink():
            raise AncestryReferenceError("ancestry_metadata_manifest_missing")
        return CachedAncestryMetadata(status, entry_dir, population_path, unrelated_path, manifest_path)

    lock_path = entry_parent / f".{cache_key}.lock"
    with lock_path.open("a+b") as lock:
        fcntl.flock(lock.fileno(), fcntl.LOCK_EX)
        if entry_dir.exists():
            return result("HIT")
        if offline:
            raise AncestryReferenceError("ancestry_metadata_cache_offline_miss")
        staging = Path(tempfile.mkdtemp(prefix=f".{cache_key}.", dir=entry_parent))
        try:
            population_path = staging / population_name
            unrelated_path = staging / unrelated_name
            downloader(str(population_asset["url"]), population_path, timeout_seconds)
            downloader(str(unrelated_asset["url"]), unrelated_path, timeout_seconds)
            _validate_file(population_path, str(population_asset["sha256"]), "ancestry_population_metadata_checksum_mismatch")
            _validate_file(unrelated_path, str(unrelated_asset["sha256"]), "ancestry_unrelated_index_checksum_mismatch")
            manifest = {
                "schema_version": "1.0.0",
                "created_at": utc_now(),
                "catalog_id": catalog["catalog_id"],
                "catalog_sha256": catalog_sha,
                "cache_key": cache_key,
                "files": {
                    "population_metadata": {"filename": population_name, "sha256": population_asset["sha256"]},
                    "unrelated_sample_index": {"filename": unrelated_name, "sha256": unrelated_asset["sha256"]},
                },
            }
            atomic_write_json(staging / manifest_path.name, manifest)
            for path in staging.iterdir():
                path.chmod(0o444)
            os.replace(staging, entry_dir)
            entry_dir.chmod(0o555)
        finally:
            if staging.exists():
                shutil.rmtree(staging, ignore_errors=True)
        return result("POPULATED")


def load_reference_samples(cached: CachedAncestryMetadata) -> tuple[ReferenceSample, ...]:
    """Joint les 3 202 métadonnées à la liste officielle des 2 504 non-apparentés."""
    metadata: dict[str, ReferenceSample] = {}
    lines = cached.population_metadata_path.read_text(encoding="utf-8").splitlines()
    if not lines or lines[0].split() != ["FamilyID", "SampleID", "FatherID", "MotherID", "Sex", "Population", "Superpopulation"]:
        raise AncestryReferenceError("ancestry_population_metadata_header_invalid")
    for line in lines[1:]:
        fields = line.split()
        if len(fields) != 7 or fields[1] in metadata:
            raise AncestryReferenceError("ancestry_population_metadata_row_invalid")
        metadata[fields[1]] = ReferenceSample(fields[1], fields[5], fields[6])
    if len(metadata) != 3202:
        raise AncestryReferenceError("ancestry_population_metadata_sample_count_mismatch")
    unrelated: set[str] = set()
    header: list[str] | None = None
    for line in cached.unrelated_index_path.read_text(encoding="utf-8").splitlines():
        if line.startswith("##"):
            continue
        fields = line.split("\t")
        if line.startswith("#ENA_FILE_PATH"):
            header = [field.removeprefix("#") for field in fields]
            continue
        if header is None or len(fields) != len(header):
            raise AncestryReferenceError("ancestry_unrelated_index_row_invalid")
        unrelated.add(fields[header.index("SAMPLE_NAME")])
    if len(unrelated) != 2504 or not unrelated <= metadata.keys():
        raise AncestryReferenceError("ancestry_unrelated_sample_set_mismatch")
    return tuple(metadata[sample_id] for sample_id in sorted(unrelated))
