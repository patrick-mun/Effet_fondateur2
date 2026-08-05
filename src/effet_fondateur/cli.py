"""Interface en ligne de commande du pipeline V2."""

from __future__ import annotations

import argparse
from pathlib import Path
from typing import Sequence

from effet_fondateur.contracts import (
    ConfigurationError,
    TableValidationError,
    load_pipeline_config,
    validate_cohorts_frozen,
    validate_samples_master,
    validate_tsv_table,
)
from effet_fondateur.contracts.samples import SAMPLES_SCHEMA_NAME
from effet_fondateur.orchestrator import PipelineError, resume_pipeline, run_pipeline


def _build_parser() -> argparse.ArgumentParser:
    """Déclare les commandes sans exécuter ni charger de donnée scientifique."""
    parser = argparse.ArgumentParser(prog="effet-fondateur")
    subparsers = parser.add_subparsers(dest="command", required=True)

    validate_parser = subparsers.add_parser(
        "validate-config",
        help="Valider une configuration YAML sans lancer le pipeline.",
    )
    validate_parser.add_argument("config", type=Path)

    run_parser = subparsers.add_parser(
        "run",
        help="Créer et exécuter un run V2 dans un dossier immuable.",
    )
    run_parser.add_argument("--config", type=Path, required=True)
    run_parser.add_argument("--runs-dir", type=Path, default=Path("data/runs"))

    resume_parser = subparsers.add_parser(
        "resume",
        help="Reprendre un run V2 après validation de son état.",
    )
    resume_parser.add_argument("--run-dir", type=Path, required=True)

    samples_parser = subparsers.add_parser(
        "validate-samples",
        help="Valider une table samples.master.tsv et ses fichiers sources.",
    )
    samples_parser.add_argument("--table", type=Path, required=True)
    samples_parser.add_argument("--source-root", type=Path, required=True)

    cohorts_parser = subparsers.add_parser(
        "validate-cohorts",
        help="Valider une table cohorts.frozen.tsv contre la table maître.",
    )
    cohorts_parser.add_argument("--table", type=Path, required=True)
    cohorts_parser.add_argument("--samples-master", type=Path, required=True)
    return parser


def main(arguments: Sequence[str] | None = None) -> int:
    """Valide les arguments et exécute la commande V2 demandée."""
    parser = _build_parser()
    parsed_arguments = parser.parse_args(arguments)

    if parsed_arguments.command == "validate-config":
        try:
            load_pipeline_config(parsed_arguments.config)
        except ConfigurationError as error:
            parser.error(str(error))
        print(f"Configuration valide : {parsed_arguments.config}")
        return 0

    if parsed_arguments.command == "run":
        try:
            run_dir = run_pipeline(parsed_arguments.config, parsed_arguments.runs_dir)
        except (ConfigurationError, PipelineError, OSError, ValueError) as error:
            parser.error(str(error))
        print(f"Run V2 créé : {run_dir}")
        return 0

    if parsed_arguments.command == "resume":
        try:
            run_dir = resume_pipeline(parsed_arguments.run_dir)
        except (ConfigurationError, PipelineError, OSError, ValueError) as error:
            parser.error(str(error))
        print(f"Run V2 repris : {run_dir}")
        return 0

    if parsed_arguments.command == "validate-samples":
        try:
            validated_table = validate_samples_master(
                parsed_arguments.table,
                parsed_arguments.source_root,
            )
        except TableValidationError as error:
            parser.error(str(error))
        print(f"Table maître valide : {validated_table.row_count} lignes")
        return 0

    if parsed_arguments.command == "validate-cohorts":
        try:
            samples_table = validate_tsv_table(
                parsed_arguments.samples_master,
                SAMPLES_SCHEMA_NAME,
            )
            known_sample_ids = {row["SAMPLE_ID"] for row in samples_table.rows}
            validated_table = validate_cohorts_frozen(
                parsed_arguments.table,
                known_sample_ids,
            )
        except TableValidationError as error:
            parser.error(str(error))
        print(f"Cohortes valides : {validated_table.row_count} lignes")
        return 0

    parser.error(f"Commande inconnue : {parsed_arguments.command}")
    return 2


if __name__ == "__main__":
    raise SystemExit(main())
