"""Catalogues et résolution des références publiques du pipeline V2."""

from effet_fondateur.references.catalog import (
    DEFAULT_REFERENCE_CATALOG_PATH,
    ReferenceAsset,
    ReferenceCatalogError,
    ResolvedReferencePanel,
    load_reference_catalog,
    resolve_reference_panel,
)
from effet_fondateur.references.cache import (
    CachedReferencePanel,
    ReferenceCacheError,
    ReferenceCacheIntegrityError,
    ReferenceCacheOfflineMiss,
    cache_reference_panel,
)
from effet_fondateur.references.window import (
    ExtractedReferenceWindow,
    ReferenceWindowError,
    ReferenceWindowIntegrityError,
    extract_reference_window,
)
from effet_fondateur.references.harmonization import (
    HarmonizedReferenceWindow,
    ReferenceHarmonizationError,
    ReferenceHarmonizationIntegrityError,
    harmonize_reference_window,
)

__all__ = [
    "DEFAULT_REFERENCE_CATALOG_PATH",
    "CachedReferencePanel",
    "ReferenceAsset",
    "ReferenceCacheError",
    "ReferenceCacheIntegrityError",
    "ReferenceCacheOfflineMiss",
    "ReferenceCatalogError",
    "ResolvedReferencePanel",
    "ExtractedReferenceWindow",
    "HarmonizedReferenceWindow",
    "ReferenceHarmonizationError",
    "ReferenceHarmonizationIntegrityError",
    "ReferenceWindowError",
    "ReferenceWindowIntegrityError",
    "load_reference_catalog",
    "cache_reference_panel",
    "extract_reference_window",
    "harmonize_reference_window",
    "resolve_reference_panel",
]
