"""Cache immuable d'extraits 1000 Genomes limités aux variants utiles."""

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

from effet_fondateur.audit import atomic_write_json, read_json, sha256_file
from effet_fondateur.contracts import validate_json_document
from effet_fondateur.orchestrator.state import utc_now


ExtractRunner = Callable[[str, Path, Path, Path, Path, int], None]


class AncestryExtractCacheError(RuntimeError):
    """Signale une sélection, une extraction ou une entrée de cache invalide."""


@dataclass(frozen=True)
class CachedReferenceExtract:
    """Décrit un extrait public vérifié et son statut de réutilisation."""

    status: str
    entry_dir: Path
    vcf_path: Path
    index_path: Path
    manifest_path: Path


def _approved_source(url: str) -> None:
    parsed = urlparse(url)
    if parsed.scheme != "https" or parsed.hostname != "ftp.1000genomes.ebi.ac.uk":
        raise AncestryExtractCacheError("ancestry_extract_source_url_not_approved")


def cache_reference_extract(
    *,
    cache_root: Path,
    panel_id: str,
    assembly: str,
    chromosome: int,
    source_url: str,
    source_vcf_md5: str,
    source_index_md5: str,
    positions_path: Path,
    samples_path: Path,
    offline: bool,
    timeout_seconds: int,
    extractor: ExtractRunner,
) -> CachedReferenceExtract:
    """Crée une fois un extrait public, puis le vérifie sans accès réseau.

    Le MD5 officiel identifie le VCF complet distant. Le SHA-256 local protège
    l'extrait publié ; il n'est jamais présenté comme le SHA du fichier source.
    """

    _approved_source(source_url)
    if (
        assembly != "GRCh38"
        or not 1 <= chromosome <= 22
        or len(source_vcf_md5) != 32
        or any(character not in "0123456789abcdef" for character in source_vcf_md5)
        or len(source_index_md5) != 32
        or any(character not in "0123456789abcdef" for character in source_index_md5)
        or not positions_path.is_file()
        or not samples_path.is_file()
        or timeout_seconds <= 0
    ):
        raise AncestryExtractCacheError("invalid_ancestry_extract_specification")
    positions_sha = sha256_file(positions_path)
    samples_sha = sha256_file(samples_path)
    identity = {
        "method_id": "bcftools_remote_variant_extract_v1",
        "panel_id": panel_id,
        "assembly": assembly,
        "chromosome": chromosome,
        "source_url": source_url,
        "source_vcf_md5": source_vcf_md5,
        "source_index_md5": source_index_md5,
        "positions_sha256": positions_sha,
        "samples_sha256": samples_sha,
    }
    cache_key = hashlib.sha256(
        json.dumps(identity, sort_keys=True, separators=(",", ":")).encode()
    ).hexdigest()
    entry_parent = cache_root / "ancestry_extracts" / panel_id / f"chr{chromosome}"
    entry_parent.mkdir(parents=True, exist_ok=True)
    if entry_parent.is_symlink():
        raise AncestryExtractCacheError("ancestry_extract_cache_parent_invalid")
    entry_dir = entry_parent / cache_key
    vcf_path = entry_dir / "reference.extract.vcf.gz"
    index_path = entry_dir / "reference.extract.vcf.gz.tbi"
    manifest_path = entry_dir / "ancestry_extract_cache_manifest.json"

    def result(status: str) -> CachedReferenceExtract:
        if entry_dir.is_symlink() or not entry_dir.is_dir() or entry_dir.stat().st_mode & 0o222:
            raise AncestryExtractCacheError("ancestry_extract_cache_permissions_invalid")
        if any(path.is_symlink() or not path.is_file() for path in (vcf_path, index_path, manifest_path)):
            raise AncestryExtractCacheError("ancestry_extract_cache_incomplete")
        manifest = read_json(manifest_path)
        validate_json_document(manifest, "ancestry_extract_cache_manifest.schema.json")
        if manifest["cache_key"] != cache_key or any(manifest[key] != value for key, value in identity.items()):
            raise AncestryExtractCacheError("ancestry_extract_cache_identity_mismatch")
        if (
            sha256_file(vcf_path) != manifest["files"]["vcf"]["sha256"]
            or sha256_file(index_path) != manifest["files"]["index"]["sha256"]
        ):
            raise AncestryExtractCacheError("ancestry_extract_cache_corrupt")
        return CachedReferenceExtract(status, entry_dir, vcf_path, index_path, manifest_path)

    lock_path = entry_parent / f".{cache_key}.lock"
    with lock_path.open("a+b") as lock:
        fcntl.flock(lock.fileno(), fcntl.LOCK_EX)
        if entry_dir.exists():
            return result("HIT")
        if offline:
            raise AncestryExtractCacheError("ancestry_extract_cache_offline_miss")
        staging = Path(tempfile.mkdtemp(prefix=f".{cache_key}.", dir=entry_parent))
        try:
            staging_vcf = staging / vcf_path.name
            staging_index = staging / index_path.name
            extractor(
                source_url,
                positions_path,
                samples_path,
                staging_vcf,
                staging_index,
                timeout_seconds,
            )
            if (
                not staging_vcf.is_file()
                or not staging_index.is_file()
                or staging_vcf.stat().st_size == 0
                or staging_index.stat().st_size == 0
            ):
                raise AncestryExtractCacheError("ancestry_extract_output_missing")
            manifest = {
                "schema_version": "1.0.0",
                "created_at": utc_now(),
                "cache_key": cache_key,
                **identity,
                "files": {
                    "vcf": {"filename": staging_vcf.name, "sha256": sha256_file(staging_vcf), "size_bytes": staging_vcf.stat().st_size},
                    "index": {"filename": staging_index.name, "sha256": sha256_file(staging_index), "size_bytes": staging_index.stat().st_size},
                },
                "checks": {
                    "official_source_identity": "PASS",
                    "selection_identity": "PASS",
                    "local_extract_integrity": "PASS",
                    "study_data_transmitted": False,
                },
            }
            validate_json_document(manifest, "ancestry_extract_cache_manifest.schema.json")
            atomic_write_json(staging / manifest_path.name, manifest)
            for path in staging.iterdir():
                path.chmod(0o444)
            os.replace(staging, entry_dir)
            entry_dir.chmod(0o555)
        finally:
            if staging.exists():
                shutil.rmtree(staging, ignore_errors=True)
        return result("POPULATED")
