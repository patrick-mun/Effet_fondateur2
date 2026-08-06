"""Catalogues et résolution des références publiques du pipeline V2."""

from effet_fondateur.references.catalog import (
    DEFAULT_REFERENCE_CATALOG_PATH,
    ReferenceAsset,
    ReferenceCatalogError,
    ResolvedReferencePanel,
    load_reference_catalog,
    resolve_reference_panel,
)

__all__ = [
    "DEFAULT_REFERENCE_CATALOG_PATH",
    "ReferenceAsset",
    "ReferenceCatalogError",
    "ResolvedReferencePanel",
    "load_reference_catalog",
    "resolve_reference_panel",
]
