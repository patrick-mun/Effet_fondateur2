import gzip
import json
import shutil
import tarfile
from pathlib import Path

import pytest

from effet_fondateur.audit import sha256_file
from effet_fondateur.references.genetic_maps import (
    GeneticMapError,
    ensure_genetic_map_cached,
    resolve_genetic_map,
)


def _archive(path: Path) -> None:
    sources = path.parent / "map_sources"
    sources.mkdir()
    for chromosome in range(1, 23):
        source = sources / f"chr{chromosome}.b38.gmap.gz"
        with gzip.open(source, "wt", encoding="utf-8") as handle:
            handle.write(f"pos\tchr\tcM\n1\t{chromosome}\t0\n200000\t{chromosome}\t2\n")
    with tarfile.open(path, "w:gz") as bundle:
        for source in sorted(sources.iterdir()):
            bundle.add(source, arcname=source.name)


def _catalog(path: Path, archive_sha256: str) -> None:
    document = {
        "schema_version": "1.0.0", "catalog_id": "synthetic_maps",
        "maps": [{
            "map_id": "synthetic_grch38", "provider": "SHAPEIT4",
            "release_id": "git_synthetic", "assembly": "GRCh38",
            "population_scope": "synthetic", "method": "synthetic",
            "archive_url": "https://raw.githubusercontent.com/example/maps/archive.tar.gz",
            "archive_sha256": archive_sha256,
            "member_template": "chr{chromosome}.b38.gmap.gz",
            "chromosomes": list(range(1, 23)),
        }],
    }
    path.write_text(json.dumps(document), encoding="utf-8")


def test_cache_populates_all_autosomes_then_reuses_archive(tmp_path: Path) -> None:
    archive = tmp_path / "maps.tar.gz"
    _archive(archive)
    catalog = tmp_path / "catalog.json"
    _catalog(catalog, sha256_file(archive))
    resolved = resolve_genetic_map(catalog, "synthetic_grch38", "GRCh38")
    calls = 0

    def copy_download(_url: str, destination: Path, _timeout: int, _chunk: int) -> None:
        nonlocal calls
        calls += 1
        shutil.copyfile(archive, destination)

    first = ensure_genetic_map_cached(
        resolved=resolved, chromosome=19, cache_root=tmp_path / "cache",
        downloader=copy_download,
    )
    second = ensure_genetic_map_cached(
        resolved=resolved, chromosome=19, cache_root=tmp_path / "cache",
        offline=True,
        downloader=lambda *_args: pytest.fail("download must not run on a cache hit"),
    )
    assert first.status == "POPULATED"
    assert second.status == "HIT"
    assert calls == 1
    assert first.archive_path.is_file()
    assert len(list(first.map_path.parent.glob("chr*.grch38.tsv"))) == 22
    assert first.map_sha256 == second.map_sha256


def test_cache_rejects_modified_normalized_map(tmp_path: Path) -> None:
    archive = tmp_path / "maps.tar.gz"
    _archive(archive)
    catalog = tmp_path / "catalog.json"
    _catalog(catalog, sha256_file(archive))
    resolved = resolve_genetic_map(catalog, "synthetic_grch38", "GRCh38")
    cached = ensure_genetic_map_cached(
        resolved=resolved, chromosome=19, cache_root=tmp_path / "cache",
        downloader=lambda _url, destination, _timeout, _chunk: shutil.copyfile(archive, destination),
    )
    cached.map_path.write_text("modified\n", encoding="utf-8")
    with pytest.raises(GeneticMapError, match="genetic_map_cache_map_mismatch"):
        ensure_genetic_map_cached(
            resolved=resolved, chromosome=19, cache_root=tmp_path / "cache", offline=True
        )
