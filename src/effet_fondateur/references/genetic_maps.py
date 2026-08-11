"""Catalogue et cache immuable des cartes génétiques publiques."""

from __future__ import annotations

import csv
import fcntl
import gzip
import json
import os
import shutil
import tarfile
import tempfile
import time
from contextlib import contextmanager
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Callable, Iterator, Literal
from urllib.parse import urlparse
from urllib.request import Request, urlopen

from effet_fondateur.audit import atomic_write_json, sha256_file
from effet_fondateur.contracts import DocumentValidationError, validate_json_document, validate_tsv_table

Downloader = Callable[[str, Path, int, int], None]
CacheStatus = Literal["HIT", "POPULATED"]
AUTOSOMES = tuple(range(1, 23))
OFFICIAL_HOST = "raw.githubusercontent.com"


class GeneticMapError(RuntimeError):
    """Signale un catalogue, téléchargement ou cache de carte invalide."""


@dataclass(frozen=True)
class ResolvedGeneticMap:
    map_id: str
    catalog_sha256: str
    provider: str
    release_id: str
    assembly: str
    population_scope: str
    method: str
    archive_url: str
    archive_sha256: str
    member_template: str


@dataclass(frozen=True)
class CachedGeneticMap:
    status: CacheStatus
    map_path: Path
    manifest_path: Path
    archive_path: Path
    archive_sha256: str
    map_sha256: str
    resolved: ResolvedGeneticMap


def resolve_genetic_map(catalog_path: Path, map_id: str, assembly: str) -> ResolvedGeneticMap:
    try:
        document = json.loads(catalog_path.read_text(encoding="utf-8"))
        validate_json_document(document, "genetic_map_catalog.schema.json")
    except (OSError, json.JSONDecodeError, DocumentValidationError) as error:
        raise GeneticMapError("genetic_map_catalog_invalid") from error
    matches = [item for item in document["maps"] if item["map_id"] == map_id]
    if len(matches) != 1:
        raise GeneticMapError("genetic_map_not_found_or_ambiguous")
    item = matches[0]
    parsed = urlparse(item["archive_url"])
    if parsed.scheme != "https" or parsed.hostname != OFFICIAL_HOST or parsed.query or parsed.fragment:
        raise GeneticMapError("genetic_map_catalog_unapproved_url")
    if item["assembly"] != assembly or tuple(sorted(item["chromosomes"])) != AUTOSOMES:
        raise GeneticMapError("genetic_map_catalog_scope_mismatch")
    return ResolvedGeneticMap(
        map_id=item["map_id"], catalog_sha256=sha256_file(catalog_path),
        provider=item["provider"], release_id=item["release_id"], assembly=item["assembly"],
        population_scope=item["population_scope"], method=item["method"],
        archive_url=item["archive_url"], archive_sha256=item["archive_sha256"],
        member_template=item["member_template"],
    )


def _download(url: str, destination: Path, timeout: int, chunk_size: int) -> None:
    request = Request(url, headers={"User-Agent": "effet-fondateur/0.1"})
    with urlopen(request, timeout=timeout) as response, destination.open("xb") as output:
        final = urlparse(response.geturl())
        if final.scheme != "https" or final.hostname != OFFICIAL_HOST:
            raise GeneticMapError("genetic_map_download_unapproved_redirect")
        while chunk := response.read(chunk_size):
            output.write(chunk)
        output.flush()
        os.fsync(output.fileno())


@contextmanager
def _lock(path: Path, timeout: float) -> Iterator[None]:
    deadline = time.monotonic() + timeout
    with path.open("a+b") as handle:
        while True:
            try:
                fcntl.flock(handle.fileno(), fcntl.LOCK_EX | fcntl.LOCK_NB)
                break
            except BlockingIOError:
                if time.monotonic() >= deadline:
                    raise GeneticMapError("genetic_map_cache_lock_timeout")
                time.sleep(0.1)
        try:
            yield
        finally:
            fcntl.flock(handle.fileno(), fcntl.LOCK_UN)


def _normalise(source: Path, destination: Path, resolved: ResolvedGeneticMap, chromosome: int) -> None:
    with gzip.open(source, "rt", encoding="utf-8", newline="") as input_file:
        reader = csv.DictReader(input_file, delimiter="\t")
        if reader.fieldnames != ["pos", "chr", "cM"]:
            raise GeneticMapError("genetic_map_source_header_invalid")
        rows = list(reader)
    if len(rows) < 2 or any(int(row["chr"]) != chromosome for row in rows):
        raise GeneticMapError("genetic_map_source_chromosome_invalid")
    with destination.open("x", encoding="utf-8", newline="") as output_file:
        writer = csv.writer(output_file, delimiter="\t", lineterminator="\n")
        writer.writerow(("MAP_ID", "ASSEMBLY", "CHROMOSOME", "POSITION_BP", "POSITION_CM"))
        for row in rows:
            writer.writerow((resolved.map_id, resolved.assembly, chromosome, row["pos"], row["cM"]))
    validate_tsv_table(destination, "genetic_map.schema.json")


def _validated_hit(entry: Path, resolved: ResolvedGeneticMap, chromosome: int) -> CachedGeneticMap | None:
    manifest_path = entry / "manifest.json"
    if not manifest_path.is_file():
        return None
    try:
        manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
        if (
            manifest["map_id"] != resolved.map_id
            or manifest["catalog_sha256"] != resolved.catalog_sha256
        ):
            raise GeneticMapError("genetic_map_cache_identity_mismatch")
        archive = entry / manifest["archive"]["path"]
        map_path = entry / manifest["maps"][str(chromosome)]["path"]
        if manifest["archive"]["sha256"] != resolved.archive_sha256 or sha256_file(archive) != resolved.archive_sha256:
            raise GeneticMapError("genetic_map_cache_archive_mismatch")
        expected_map_sha = manifest["maps"][str(chromosome)]["sha256"]
        if sha256_file(map_path) != expected_map_sha:
            raise GeneticMapError("genetic_map_cache_map_mismatch")
        validate_tsv_table(map_path, "genetic_map.schema.json")
    except (OSError, KeyError, json.JSONDecodeError) as error:
        raise GeneticMapError("genetic_map_cache_manifest_invalid") from error
    return CachedGeneticMap("HIT", map_path, manifest_path, archive, resolved.archive_sha256, expected_map_sha, resolved)


def ensure_genetic_map_cached(*, resolved: ResolvedGeneticMap, chromosome: int, cache_root: Path,
                              offline: bool = False, lock_timeout_seconds: int = 60,
                              download_timeout_seconds: int = 600, download_chunk_size: int = 1_048_576,
                              downloader: Downloader = _download) -> CachedGeneticMap:
    """Télécharge une seule fois l'archive et normalise les 22 autosomes atomiquement."""
    if chromosome not in AUTOSOMES:
        raise GeneticMapError("genetic_map_unsupported_chromosome")
    entry = cache_root / resolved.provider / resolved.map_id / resolved.release_id / resolved.assembly / resolved.archive_sha256
    entry.parent.mkdir(parents=True, exist_ok=True)
    with _lock(entry.parent / f".{resolved.archive_sha256}.lock", lock_timeout_seconds):
        hit = _validated_hit(entry, resolved, chromosome)
        if hit is not None:
            return hit
        if entry.exists():
            raise GeneticMapError("genetic_map_cache_incomplete_entry")
        if offline:
            raise GeneticMapError("genetic_map_cache_offline_miss")
        staging = Path(tempfile.mkdtemp(prefix=".genetic-map-", dir=entry.parent))
        try:
            archive = staging / "genetic_maps.b38.tar.gz"
            downloader(resolved.archive_url, archive, download_timeout_seconds, download_chunk_size)
            if sha256_file(archive) != resolved.archive_sha256:
                raise GeneticMapError("genetic_map_archive_sha256_mismatch")
            source_dir, normalized_dir = staging / "source", staging / "normalized"
            source_dir.mkdir(); normalized_dir.mkdir()
            expected = {resolved.member_template.format(chromosome=value): value for value in AUTOSOMES}
            with tarfile.open(archive, "r:gz") as bundle:
                members = {member.name: member for member in bundle.getmembers()}
                if any(name not in members or not members[name].isfile() for name in expected):
                    raise GeneticMapError("genetic_map_archive_members_missing")
                for name, value in expected.items():
                    source = bundle.extractfile(members[name])
                    if source is None:
                        raise GeneticMapError("genetic_map_archive_member_unreadable")
                    source_path = source_dir / name
                    with source, source_path.open("xb") as output:
                        shutil.copyfileobj(source, output)
                    _normalise(source_path, normalized_dir / f"chr{value}.grch38.tsv", resolved, value)
            maps = {str(value): {"path": f"normalized/chr{value}.grch38.tsv", "sha256": sha256_file(normalized_dir / f"chr{value}.grch38.tsv")} for value in AUTOSOMES}
            manifest = {"schema_version": "1.0.0", "map_id": resolved.map_id,
                        "catalog_sha256": resolved.catalog_sha256,
                        "archive": {"path": archive.name, "source_url": resolved.archive_url, "sha256": resolved.archive_sha256},
                        "maps": maps}
            atomic_write_json(staging / "manifest.json", manifest)
            os.replace(staging, entry)
        except Exception:
            shutil.rmtree(staging, ignore_errors=True)
            raise
        cached = _validated_hit(entry, resolved, chromosome)
        if cached is None:
            raise GeneticMapError("genetic_map_cache_publish_failed")
        return CachedGeneticMap("POPULATED", cached.map_path, cached.manifest_path, cached.archive_path, cached.archive_sha256, cached.map_sha256, resolved)
