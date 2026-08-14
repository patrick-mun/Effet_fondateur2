import hashlib
import json
from pathlib import Path

import pytest

from effet_fondateur.ancestry import (
    AncestryReferenceError,
    cache_ancestry_metadata,
    load_reference_samples,
)


BASE_URL = "https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage"


def _contents() -> dict[str, bytes]:
    metadata = "FamilyID SampleID FatherID MotherID Sex Population Superpopulation\n"
    metadata += "".join(
        f"F{index:04d} S{index:04d} 0 0 1 POP SUP\n"
        for index in range(3202)
    )
    unrelated = "##FileDate=synthetic\n#ENA_FILE_PATH\tSAMPLE_NAME\n"
    unrelated += "".join(
        f"file{index}\tS{index:04d}{' ' if index < 3 else ''}\n"
        for index in range(2504)
    )
    return {"population.txt": metadata.encode(), "unrelated.index": unrelated.encode()}


def _catalog(path: Path, contents: dict[str, bytes]) -> None:
    document = {
        "schema_version": "1.0.0",
        "catalog_id": "1kg_ancestry_reference_grch38_v1",
        "panel_id": "1kg_3202_high_coverage_20220422",
        "assembly": "GRCh38",
        "population_metadata": {
            "url": f"{BASE_URL}/population.txt",
            "sha256": hashlib.sha256(contents["population.txt"]).hexdigest(),
            "sample_count": 3202,
        },
        "unrelated_sample_index": {
            "url": f"{BASE_URL}/unrelated.index",
            "sha256": hashlib.sha256(contents["unrelated.index"]).hexdigest(),
            "sample_count": 2504,
        },
        "model_scope": "UNRELATED_REFERENCE_SAMPLES_ONLY",
        "projection_scope": "STUDY_SAMPLES_AND_TARGET_HAPLOTYPES",
        "interpretation_policy": "RELATIVE_REFERENCE_POSITIONING_ONLY",
        "prohibited_interpretations": [
            "ETHNIC_IDENTITY_ASSIGNMENT",
            "GENEALOGICAL_ANCESTOR_IDENTIFICATION",
            "LOCAL_ANCESTRY_PROOF_FROM_PCA",
        ],
    }
    path.write_text(json.dumps(document), encoding="utf-8")


def test_ancestry_metadata_is_downloaded_once_then_reused(tmp_path: Path) -> None:
    contents = _contents()
    catalog_path = tmp_path / "catalog.json"
    _catalog(catalog_path, contents)
    calls: list[str] = []

    def downloader(url: str, destination: Path, timeout: int) -> None:
        calls.append(url)
        destination.write_bytes(contents[Path(url).name])

    populated = cache_ancestry_metadata(
        catalog_path=catalog_path,
        cache_root=tmp_path / "cache",
        offline=False,
        downloader=downloader,
    )
    reused = cache_ancestry_metadata(
        catalog_path=catalog_path,
        cache_root=tmp_path / "cache",
        offline=True,
        downloader=lambda *_: pytest.fail("network must not be used on cache hit"),
    )

    assert populated.status == "POPULATED"
    assert reused.status == "HIT"
    assert len(calls) == 2
    samples = load_reference_samples(reused)
    assert len(samples) == 2504
    assert samples[0].sample_id == "S0000"
    assert samples[0].population == "POP"


def test_ancestry_metadata_offline_miss_blocks(tmp_path: Path) -> None:
    contents = _contents()
    catalog_path = tmp_path / "catalog.json"
    _catalog(catalog_path, contents)

    with pytest.raises(AncestryReferenceError, match="offline_miss"):
        cache_ancestry_metadata(
            catalog_path=catalog_path,
            cache_root=tmp_path / "cache",
            offline=True,
        )


def test_ancestry_metadata_corruption_blocks_reuse(tmp_path: Path) -> None:
    contents = _contents()
    catalog_path = tmp_path / "catalog.json"
    _catalog(catalog_path, contents)

    cached = cache_ancestry_metadata(
        catalog_path=catalog_path,
        cache_root=tmp_path / "cache",
        offline=False,
        downloader=lambda url, destination, timeout: destination.write_bytes(
            contents[Path(url).name]
        ),
    )
    cached.entry_dir.chmod(0o755)
    cached.population_metadata_path.chmod(0o644)
    cached.population_metadata_path.write_bytes(b"corrupt")
    cached.population_metadata_path.chmod(0o444)
    cached.entry_dir.chmod(0o555)

    with pytest.raises(AncestryReferenceError, match="corrupt"):
        cache_ancestry_metadata(
            catalog_path=catalog_path,
            cache_root=tmp_path / "cache",
            offline=True,
        )
