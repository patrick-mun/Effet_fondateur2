"""Cache partagé, immuable et vérifié des panels publics de référence."""

from __future__ import annotations

import fcntl
import hashlib
import json
import math
import os
import re
import shutil
import tempfile
import time
from contextlib import contextmanager
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Callable, Iterator, Literal
from urllib.parse import urlparse
from urllib.request import Request, urlopen

from effet_fondateur.audit import atomic_write_json, md5_file, read_json, sha256_file
from effet_fondateur.contracts import DocumentValidationError, validate_json_document
from effet_fondateur.orchestrator.state import utc_now
from effet_fondateur.references.catalog import ResolvedReferencePanel


Downloader = Callable[[str, Path, int, int], None]
CacheStatus = Literal["HIT", "POPULATED"]
SAFE_COMPONENT_PATTERN = re.compile(r"^[A-Za-z0-9][A-Za-z0-9._-]*$")


class ReferenceCacheError(RuntimeError):
    """Signale une erreur de verrouillage, téléchargement ou configuration."""


class ReferenceCacheIntegrityError(ReferenceCacheError):
    """Signale une empreinte ou une entrée publiée incohérente."""


class ReferenceCacheOfflineMiss(ReferenceCacheError):
    """Signale une référence absente alors que le réseau est interdit."""


@dataclass(frozen=True)
class CachedReferencePanel:
    """Décrit les quatre fichiers vérifiés d'une entrée de cache."""

    status: CacheStatus
    cache_entry_key: str
    entry_dir: Path
    manifest_path: Path
    vcf_path: Path
    index_path: Path
    readme_path: Path
    source_manifest_path: Path
    vcf_sha256: str
    index_sha256: str
    readme_sha256: str
    source_manifest_sha256: str


@dataclass(frozen=True)
class _ExpectedFile:
    role: str
    filename: str
    source_url: str
    official_md5: str | None
    expected_sha256: str | None


def _positive_integer(value: int, name: str) -> int:
    if isinstance(value, bool) or not isinstance(value, int) or value <= 0:
        raise ReferenceCacheError(f"invalid_reference_cache_parameter:{name}")
    return value


def _positive_number(value: float, name: str) -> float:
    if (
        isinstance(value, bool)
        or not isinstance(value, (int, float))
        or not math.isfinite(value)
        or value <= 0
    ):
        raise ReferenceCacheError(f"invalid_reference_cache_parameter:{name}")
    return float(value)


def _safe_component(value: str, field: str) -> str:
    if not isinstance(value, str) or SAFE_COMPONENT_PATTERN.fullmatch(value) is None:
        raise ReferenceCacheError(f"unsafe_reference_cache_component:{field}")
    return value


def _url_filename(url: str, field: str) -> str:
    filename = Path(urlparse(url).path).name
    if not filename or filename in {".", ".."} or "/" in filename or "\\" in filename:
        raise ReferenceCacheError(f"unsafe_reference_cache_filename:{field}")
    return filename


def _expected_files(resolved: ResolvedReferencePanel) -> dict[str, _ExpectedFile]:
    expected = {
        "vcf": _ExpectedFile(
            "vcf",
            resolved.vcf.filename,
            resolved.vcf.url,
            resolved.vcf.official_md5,
            None,
        ),
        "index": _ExpectedFile(
            "index",
            resolved.index.filename,
            resolved.index.url,
            resolved.index.official_md5,
            None,
        ),
        "readme": _ExpectedFile(
            "readme",
            _url_filename(resolved.readme_url, "readme"),
            resolved.readme_url,
            None,
            resolved.readme_sha256,
        ),
        "source_manifest": _ExpectedFile(
            "source_manifest",
            _url_filename(resolved.manifest_url, "source_manifest"),
            resolved.manifest_url,
            None,
            resolved.manifest_sha256,
        ),
    }
    filenames = [item.filename for item in expected.values()]
    if len(set(filenames)) != len(filenames):
        raise ReferenceCacheError("reference_cache_filename_collision")
    return expected


def _cache_entry_key(
    resolved: ResolvedReferencePanel,
    expected_files: dict[str, _ExpectedFile],
) -> str:
    identity = {
        "catalog_sha256": resolved.catalog_sha256,
        "panel_id": resolved.panel_id,
        "provider": resolved.provider,
        "release_id": resolved.release_id,
        "assembly": resolved.assembly,
        "chromosome": resolved.chromosome,
        "files": {
            role: {
                "filename": item.filename,
                "source_url": item.source_url,
                "official_md5": item.official_md5,
                "expected_sha256": item.expected_sha256,
            }
            for role, item in sorted(expected_files.items())
        },
    }
    payload = json.dumps(identity, sort_keys=True, separators=(",", ":"))
    return hashlib.sha256(payload.encode("utf-8")).hexdigest()


def _entry_parent(cache_root: Path, resolved: ResolvedReferencePanel) -> Path:
    if (
        isinstance(resolved.chromosome, bool)
        or not isinstance(resolved.chromosome, int)
        or not 1 <= resolved.chromosome <= 22
    ):
        raise ReferenceCacheError("invalid_reference_cache_chromosome")
    components = (
        _safe_component(resolved.provider, "provider"),
        _safe_component(resolved.panel_id, "panel_id"),
        _safe_component(resolved.release_id, "release_id"),
        _safe_component(resolved.assembly, "assembly"),
        f"chr{resolved.chromosome}",
    )
    root = cache_root.expanduser()
    root.mkdir(parents=True, exist_ok=True)
    if root.is_symlink() or not root.is_dir():
        raise ReferenceCacheError("invalid_reference_cache_root")
    current = root
    for component in components:
        current = current / component
        if current.is_symlink():
            raise ReferenceCacheError("invalid_reference_cache_parent")
        current.mkdir(exist_ok=True)
        if not current.is_dir():
            raise ReferenceCacheError("invalid_reference_cache_parent")
    return current


@contextmanager
def _exclusive_lock(lock_path: Path, timeout_seconds: float) -> Iterator[None]:
    deadline = time.monotonic() + timeout_seconds
    with lock_path.open("a+b") as lock_file:
        while True:
            try:
                fcntl.flock(lock_file.fileno(), fcntl.LOCK_EX | fcntl.LOCK_NB)
                break
            except BlockingIOError:
                if time.monotonic() >= deadline:
                    raise ReferenceCacheError("reference_cache_lock_timeout")
                time.sleep(min(0.1, max(0.0, deadline - time.monotonic())))
        try:
            yield
        finally:
            fcntl.flock(lock_file.fileno(), fcntl.LOCK_UN)


def _download_to_path(
    source_url: str,
    destination: Path,
    timeout_seconds: int,
    chunk_size: int,
) -> None:
    request = Request(source_url, headers={"User-Agent": "effet-fondateur/0.1"})
    try:
        with urlopen(request, timeout=timeout_seconds) as response:
            final_url = response.geturl()
            source = urlparse(source_url)
            final = urlparse(final_url)
            if source.scheme != "https" or final.scheme != "https":
                raise ReferenceCacheError("reference_download_requires_https")
            if source.hostname != final.hostname:
                raise ReferenceCacheError("reference_download_cross_host_redirect")
            with destination.open("xb") as output_file:
                while chunk := response.read(chunk_size):
                    output_file.write(chunk)
                output_file.flush()
                os.fsync(output_file.fileno())
    except ReferenceCacheError:
        raise
    except (OSError, ValueError) as error:
        raise ReferenceCacheError(
            f"reference_download_failed:{destination.name}"
        ) from error


def _file_record(path: Path, expected: _ExpectedFile) -> dict[str, Any]:
    if path.is_symlink() or not path.is_file():
        raise ReferenceCacheIntegrityError(
            f"reference_download_not_regular_file:{expected.role}"
        )
    local_sha256 = sha256_file(path)
    if expected.official_md5 is not None:
        local_md5 = md5_file(path)
        if local_md5 != expected.official_md5:
            raise ReferenceCacheIntegrityError(
                f"reference_official_md5_mismatch:{expected.role}"
            )
        return {
            "filename": expected.filename,
            "source_url": expected.source_url,
            "size_bytes": path.stat().st_size,
            "official_md5": expected.official_md5,
            "local_sha256": local_sha256,
        }
    if local_sha256 != expected.expected_sha256:
        raise ReferenceCacheIntegrityError(
            f"reference_expected_sha256_mismatch:{expected.role}"
        )
    return {
        "filename": expected.filename,
        "source_url": expected.source_url,
        "size_bytes": path.stat().st_size,
        "expected_sha256": expected.expected_sha256,
        "local_sha256": local_sha256,
    }


def _manifest_document(
    resolved: ResolvedReferencePanel,
    cache_entry_key: str,
    file_records: dict[str, dict[str, Any]],
) -> dict[str, Any]:
    return {
        "schema_version": "1.0.0",
        "cache_entry_key": cache_entry_key,
        "catalog_id": resolved.catalog_id,
        "catalog_sha256": resolved.catalog_sha256,
        "panel_id": resolved.panel_id,
        "provider": resolved.provider,
        "release_id": resolved.release_id,
        "assembly": resolved.assembly,
        "chromosome": resolved.chromosome,
        "sample_count": resolved.sample_count,
        "sample_scope": resolved.sample_scope,
        "population_scope": resolved.population_scope,
        "created_at": utc_now(),
        "files": file_records,
    }


def _expected_identity(
    resolved: ResolvedReferencePanel,
    cache_entry_key: str,
) -> dict[str, Any]:
    return {
        "cache_entry_key": cache_entry_key,
        "catalog_id": resolved.catalog_id,
        "catalog_sha256": resolved.catalog_sha256,
        "panel_id": resolved.panel_id,
        "provider": resolved.provider,
        "release_id": resolved.release_id,
        "assembly": resolved.assembly,
        "chromosome": resolved.chromosome,
        "sample_count": resolved.sample_count,
        "sample_scope": resolved.sample_scope,
        "population_scope": resolved.population_scope,
    }


def _cached_result(
    entry_dir: Path,
    cache_entry_key: str,
    manifest: dict[str, Any],
    status: CacheStatus,
) -> CachedReferencePanel:
    files = manifest["files"]
    return CachedReferencePanel(
        status=status,
        cache_entry_key=cache_entry_key,
        entry_dir=entry_dir,
        manifest_path=entry_dir / "reference_cache_manifest.json",
        vcf_path=entry_dir / files["vcf"]["filename"],
        index_path=entry_dir / files["index"]["filename"],
        readme_path=entry_dir / files["readme"]["filename"],
        source_manifest_path=entry_dir / files["source_manifest"]["filename"],
        vcf_sha256=files["vcf"]["local_sha256"],
        index_sha256=files["index"]["local_sha256"],
        readme_sha256=files["readme"]["local_sha256"],
        source_manifest_sha256=files["source_manifest"]["local_sha256"],
    )


def _validate_cached_entry(
    entry_dir: Path,
    resolved: ResolvedReferencePanel,
    cache_entry_key: str,
    expected_files: dict[str, _ExpectedFile],
    status: CacheStatus,
) -> CachedReferencePanel:
    if entry_dir.is_symlink() or not entry_dir.is_dir():
        raise ReferenceCacheIntegrityError("reference_cache_entry_not_directory")
    if entry_dir.stat().st_mode & 0o222:
        raise ReferenceCacheIntegrityError("reference_cache_entry_writable")
    manifest_path = entry_dir / "reference_cache_manifest.json"
    if manifest_path.is_symlink() or not manifest_path.is_file():
        raise ReferenceCacheIntegrityError("reference_cache_manifest_missing")
    if manifest_path.stat().st_mode & 0o222:
        raise ReferenceCacheIntegrityError("reference_cache_manifest_writable")
    try:
        manifest = read_json(manifest_path)
        validate_json_document(manifest, "reference_cache_manifest.schema.json")
    except (ValueError, DocumentValidationError) as error:
        raise ReferenceCacheIntegrityError(
            "reference_cache_manifest_invalid"
        ) from error
    for field, expected_value in _expected_identity(
        resolved, cache_entry_key
    ).items():
        if manifest[field] != expected_value:
            raise ReferenceCacheIntegrityError(
                f"reference_cache_identity_mismatch:{field}"
            )
    if set(manifest["files"]) != set(expected_files):
        raise ReferenceCacheIntegrityError("reference_cache_file_roles_mismatch")
    for role, expected in expected_files.items():
        record = manifest["files"][role]
        if record["filename"] != expected.filename:
            raise ReferenceCacheIntegrityError(
                f"reference_cache_filename_mismatch:{role}"
            )
        if record["source_url"] != expected.source_url:
            raise ReferenceCacheIntegrityError(
                f"reference_cache_source_url_mismatch:{role}"
            )
        path = entry_dir / expected.filename
        if path.is_symlink() or not path.is_file():
            raise ReferenceCacheIntegrityError(f"reference_cache_file_missing:{role}")
        if path.stat().st_mode & 0o222:
            raise ReferenceCacheIntegrityError(
                f"reference_cache_file_writable:{role}"
            )
        if path.stat().st_size != record["size_bytes"]:
            raise ReferenceCacheIntegrityError(f"reference_cache_size_mismatch:{role}")
        if sha256_file(path) != record["local_sha256"]:
            raise ReferenceCacheIntegrityError(
                f"reference_cache_sha256_mismatch:{role}"
            )
        if expected.official_md5 is not None:
            if (
                record["official_md5"] != expected.official_md5
                or md5_file(path) != expected.official_md5
            ):
                raise ReferenceCacheIntegrityError(
                    f"reference_cache_md5_mismatch:{role}"
                )
        elif (
            record["expected_sha256"] != expected.expected_sha256
            or record["local_sha256"] != expected.expected_sha256
        ):
            raise ReferenceCacheIntegrityError(
                f"reference_cache_expected_sha256_mismatch:{role}"
            )
    return _cached_result(entry_dir, cache_entry_key, manifest, status)


def _fsync_directory(path: Path) -> None:
    descriptor = os.open(path, os.O_RDONLY)
    try:
        os.fsync(descriptor)
    finally:
        os.close(descriptor)


def _cleanup_stale_staging(
    parent: Path,
    cache_entry_key: str,
    lock_path: Path,
) -> None:
    for candidate in parent.glob(f".{cache_entry_key}.*"):
        if candidate == lock_path:
            continue
        if candidate.is_symlink() or not candidate.is_dir():
            raise ReferenceCacheIntegrityError(
                "reference_cache_unexpected_staging_entry"
            )
        shutil.rmtree(candidate)


def cache_reference_panel(
    resolved: ResolvedReferencePanel,
    cache_root: Path,
    *,
    offline: bool = False,
    lock_timeout_seconds: float = 60,
    download_timeout_seconds: int = 300,
    chunk_size: int = 1024 * 1024,
    downloader: Downloader = _download_to_path,
) -> CachedReferencePanel:
    """Télécharge ou réutilise une entrée immuable après contrôle intégral."""
    if not isinstance(offline, bool):
        raise ReferenceCacheError("invalid_reference_cache_parameter:offline")
    lock_timeout = _positive_number(lock_timeout_seconds, "lock_timeout_seconds")
    download_timeout = _positive_integer(
        download_timeout_seconds, "download_timeout_seconds"
    )
    resolved_chunk_size = _positive_integer(chunk_size, "chunk_size")
    expected_files = _expected_files(resolved)
    cache_entry_key = _cache_entry_key(resolved, expected_files)
    parent = _entry_parent(cache_root, resolved)
    entry_dir = parent / cache_entry_key
    lock_path = parent / f".{cache_entry_key}.lock"
    with _exclusive_lock(lock_path, lock_timeout):
        _cleanup_stale_staging(parent, cache_entry_key, lock_path)
        if entry_dir.exists() or entry_dir.is_symlink():
            return _validate_cached_entry(
                entry_dir,
                resolved,
                cache_entry_key,
                expected_files,
                "HIT",
            )
        if offline:
            raise ReferenceCacheOfflineMiss("reference_cache_offline_miss")
        staging_dir = Path(
            tempfile.mkdtemp(prefix=f".{cache_entry_key}.", dir=parent)
        )
        published = False
        try:
            file_records: dict[str, dict[str, Any]] = {}
            for role, expected in expected_files.items():
                destination = staging_dir / expected.filename
                downloader(
                    expected.source_url,
                    destination,
                    download_timeout,
                    resolved_chunk_size,
                )
                file_records[role] = _file_record(destination, expected)
            manifest = _manifest_document(
                resolved, cache_entry_key, file_records
            )
            validate_json_document(
                manifest, "reference_cache_manifest.schema.json"
            )
            manifest_path = staging_dir / "reference_cache_manifest.json"
            atomic_write_json(manifest_path, manifest)
            for path in staging_dir.iterdir():
                path.chmod(0o444)
            _fsync_directory(staging_dir)
            os.replace(staging_dir, entry_dir)
            published = True
            entry_dir.chmod(0o555)
            _fsync_directory(parent)
        except (DocumentValidationError, OSError) as error:
            raise ReferenceCacheError("reference_cache_publication_failed") from error
        finally:
            if not published:
                shutil.rmtree(staging_dir, ignore_errors=True)
        return _validate_cached_entry(
            entry_dir,
            resolved,
            cache_entry_key,
            expected_files,
            "POPULATED",
        )
