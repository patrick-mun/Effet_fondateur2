import hashlib
import json
import stat
import threading
import time
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import pytest

from effet_fondateur.references import (
    CachedReferencePanel,
    ReferenceAsset,
    ReferenceCacheError,
    ReferenceCacheIntegrityError,
    ReferenceCacheOfflineMiss,
    ResolvedReferencePanel,
    cache_reference_panel,
)


def md5_bytes(content: bytes) -> str:
    return hashlib.md5(content, usedforsecurity=False).hexdigest()


def sha256_bytes(content: bytes) -> str:
    return hashlib.sha256(content).hexdigest()


def reference_fixture() -> tuple[ResolvedReferencePanel, dict[str, bytes]]:
    base_url = (
        "https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/"
        "synthetic_reference"
    )
    filenames = {
        "vcf": "synthetic.chr19.phased.vcf.gz",
        "index": "synthetic.chr19.phased.vcf.gz.tbi",
        "readme": "README.pdf",
        "source_manifest": "manifest.txt",
    }
    contents = {
        f"{base_url}/{filenames['vcf']}": b"synthetic phased vcf\n",
        f"{base_url}/{filenames['index']}": b"synthetic tabix index\n",
        f"{base_url}/{filenames['readme']}": b"synthetic readme\n",
        f"{base_url}/{filenames['source_manifest']}": b"synthetic manifest\n",
    }
    resolved = ResolvedReferencePanel(
        catalog_id="synthetic_catalog_v1",
        catalog_sha256="a" * 64,
        panel_id="synthetic_panel_v1",
        provider="IGSR_1000_GENOMES",
        release_id="synthetic_release_v1",
        assembly="GRCh38",
        chromosome=19,
        sample_count=3,
        sample_scope="ALL_RELEASE_SAMPLES",
        population_scope="ALL_AVAILABLE_1KGP_POPULATIONS",
        readme_url=f"{base_url}/{filenames['readme']}",
        readme_sha256=sha256_bytes(contents[f"{base_url}/{filenames['readme']}"]),
        manifest_url=f"{base_url}/{filenames['source_manifest']}",
        manifest_sha256=sha256_bytes(
            contents[f"{base_url}/{filenames['source_manifest']}"]
        ),
        vcf=ReferenceAsset(
            filenames["vcf"],
            f"{base_url}/{filenames['vcf']}",
            md5_bytes(contents[f"{base_url}/{filenames['vcf']}"]),
        ),
        index=ReferenceAsset(
            filenames["index"],
            f"{base_url}/{filenames['index']}",
            md5_bytes(contents[f"{base_url}/{filenames['index']}"]),
        ),
    )
    return resolved, contents


def downloader_for(
    contents: dict[str, bytes],
    calls: list[str] | None = None,
    delay_seconds: float = 0,
):
    def download(
        source_url: str,
        destination: Path,
        timeout_seconds: int,
        chunk_size: int,
    ) -> None:
        assert timeout_seconds > 0
        assert chunk_size > 0
        if calls is not None:
            calls.append(source_url)
        if delay_seconds:
            time.sleep(delay_seconds)
        destination.write_bytes(contents[source_url])

    return download


def test_cache_populates_atomically_and_reuses_verified_entry(tmp_path: Path) -> None:
    resolved, contents = reference_fixture()
    calls: list[str] = []

    populated = cache_reference_panel(
        resolved,
        tmp_path / "cache",
        downloader=downloader_for(contents, calls),
    )

    assert populated.status == "POPULATED"
    assert len(calls) == 4
    assert populated.vcf_path.read_bytes() == contents[resolved.vcf.url]
    assert populated.index_path.read_bytes() == contents[resolved.index.url]
    assert stat.S_IMODE(populated.entry_dir.stat().st_mode) == 0o555
    assert stat.S_IMODE(populated.vcf_path.stat().st_mode) == 0o444
    manifest = json.loads(populated.manifest_path.read_text(encoding="utf-8"))
    assert manifest["files"]["vcf"]["official_md5"] == resolved.vcf.official_md5
    assert manifest["files"]["vcf"]["local_sha256"] == populated.vcf_sha256

    def forbidden_download(*args, **kwargs) -> None:
        raise AssertionError("network must not be used for a cache hit")

    cached = cache_reference_panel(
        resolved,
        tmp_path / "cache",
        offline=True,
        downloader=forbidden_download,
    )

    assert cached.status == "HIT"
    assert cached.cache_entry_key == populated.cache_entry_key
    assert cached.vcf_sha256 == populated.vcf_sha256


def test_offline_mode_blocks_cache_miss_without_downloading(tmp_path: Path) -> None:
    resolved, _ = reference_fixture()

    def forbidden_download(*args, **kwargs) -> None:
        raise AssertionError("offline mode must not call the downloader")

    with pytest.raises(ReferenceCacheOfflineMiss, match="offline_miss"):
        cache_reference_panel(
            resolved,
            tmp_path / "cache",
            offline=True,
            downloader=forbidden_download,
        )


def test_official_md5_mismatch_never_publishes_entry(tmp_path: Path) -> None:
    resolved, contents = reference_fixture()
    invalid_resolved = ResolvedReferencePanel(
        **{
            **resolved.__dict__,
            "vcf": ReferenceAsset(
                resolved.vcf.filename,
                resolved.vcf.url,
                "0" * 32,
            ),
        }
    )

    with pytest.raises(ReferenceCacheIntegrityError, match="official_md5_mismatch"):
        cache_reference_panel(
            invalid_resolved,
            tmp_path / "cache",
            downloader=downloader_for(contents),
        )

    assert not list((tmp_path / "cache").rglob("reference_cache_manifest.json"))


def test_corrupted_immutable_entry_blocks_without_replacement(tmp_path: Path) -> None:
    resolved, contents = reference_fixture()
    populated = cache_reference_panel(
        resolved,
        tmp_path / "cache",
        downloader=downloader_for(contents),
    )
    populated.vcf_path.chmod(0o644)
    populated.vcf_path.write_bytes(b"tampered\n")

    def forbidden_download(*args, **kwargs) -> None:
        raise AssertionError("a corrupt published entry must not be overwritten")

    with pytest.raises(ReferenceCacheIntegrityError, match="file_writable:vcf"):
        cache_reference_panel(
            resolved,
            tmp_path / "cache",
            downloader=forbidden_download,
        )

    assert populated.vcf_path.read_bytes() == b"tampered\n"


def test_concurrent_requests_download_each_asset_once(tmp_path: Path) -> None:
    resolved, contents = reference_fixture()
    calls: list[str] = []
    calls_lock = threading.Lock()

    def download(
        source_url: str,
        destination: Path,
        timeout_seconds: int,
        chunk_size: int,
    ) -> None:
        with calls_lock:
            calls.append(source_url)
        time.sleep(0.03)
        destination.write_bytes(contents[source_url])

    def populate() -> CachedReferencePanel:
        return cache_reference_panel(
            resolved,
            tmp_path / "cache",
            lock_timeout_seconds=5,
            downloader=download,
        )

    with ThreadPoolExecutor(max_workers=2) as executor:
        results = list(executor.map(lambda _: populate(), range(2)))

    assert sorted(result.status for result in results) == ["HIT", "POPULATED"]
    assert len(calls) == 4
    assert len({result.cache_entry_key for result in results}) == 1


def test_metadata_sha256_mismatch_never_publishes_entry(tmp_path: Path) -> None:
    resolved, contents = reference_fixture()
    invalid_resolved = ResolvedReferencePanel(
        **{
            **resolved.__dict__,
            "readme_sha256": "0" * 64,
        }
    )

    with pytest.raises(
        ReferenceCacheIntegrityError, match="expected_sha256_mismatch:readme"
    ):
        cache_reference_panel(
            invalid_resolved,
            tmp_path / "cache",
            downloader=downloader_for(contents),
        )


def test_retry_removes_only_stale_staging_directory(tmp_path: Path) -> None:
    resolved, contents = reference_fixture()
    populated = cache_reference_panel(
        resolved,
        tmp_path / "first_cache",
        downloader=downloader_for(contents),
    )
    cache_root = tmp_path / "retry_cache"
    parent = (
        cache_root
        / resolved.provider
        / resolved.panel_id
        / resolved.release_id
        / resolved.assembly
        / f"chr{resolved.chromosome}"
    )
    parent.mkdir(parents=True)
    stale = parent / f".{populated.cache_entry_key}.interrupted"
    stale.mkdir()
    (stale / "partial.vcf.gz").write_bytes(b"partial")

    retried = cache_reference_panel(
        resolved,
        cache_root,
        downloader=downloader_for(contents),
    )

    assert retried.status == "POPULATED"
    assert not stale.exists()


def test_cache_rejects_nonfinite_lock_timeout(tmp_path: Path) -> None:
    resolved, contents = reference_fixture()

    with pytest.raises(ReferenceCacheError, match="lock_timeout_seconds"):
        cache_reference_panel(
            resolved,
            tmp_path / "cache",
            lock_timeout_seconds=float("nan"),
            downloader=downloader_for(contents),
        )
