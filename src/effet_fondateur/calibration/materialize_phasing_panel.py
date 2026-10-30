"""Commande bornée D23a ; aucun lancement de SHAPEIT5 ou de l'étude."""

from __future__ import annotations

import argparse
import logging
from pathlib import Path

from effet_fondateur.calibration.phasing_preparation import finalize_preparation, prepare_panel


def main() -> None:
    """Résout les seules sources publiques explicites et publie un nouveau run."""
    parser = argparse.ArgumentParser(description=__doc__)
    for argument in ("configuration", "workspace-root"):
        parser.add_argument("--" + argument, type=Path, required=True)
    for argument in ("partition-root", "source-directory", "chip-package", "from-preparation"):
        parser.add_argument("--" + argument, type=Path)
    parser.add_argument("--map-archive", type=Path, help="Archive locale épinglée ; sinon acquisition publique vérifiée dans le nouveau workspace")
    parser.add_argument("--run-id", required=True)
    parser.add_argument("--bcftools-command", default="bcftools")
    args = parser.parse_args()
    logging.basicConfig(level=logging.INFO, format="%(message)s")
    if args.from_preparation is not None:
        if any(value is not None for value in (args.partition_root, args.source_directory, args.chip_package, args.map_archive)):
            parser.error("--from-preparation est exclusif des sources d'extraction")
        finalize_preparation(source_root=args.from_preparation, configuration=args.configuration,
                             workspace_root=args.workspace_root, run_id=args.run_id)
    else:
        if any(value is None for value in (args.partition_root, args.source_directory, args.chip_package)):
            parser.error("Les trois sources publiques sont obligatoires sans --from-preparation")
        values = vars(args)
        values.pop("from_preparation")
        prepare_panel(**values)


if __name__ == "__main__":
    main()
