# ld.py
"""
Module LD : calcul du linkage disequilibrium (r² et D)
"""
import os
import subprocess
import logging
import pandas as pd

def run_ld(plink_prefix: str, output_dir: str):
    """
    Calcule le linkage disequilibrium (r² et D) à partir des données PLINK

    :param plink_prefix: chemin du fichier PLINK (sans extension)
    :param output_dir: dossier de sortie
    """
    os.makedirs(output_dir, exist_ok=True)
    output_base = os.path.join(output_dir, "ld")

    cmd = (
        f"plink --bfile {plink_prefix} "
        f"--r2 dprime with-freqs --ld-window-kb 1000 --ld-window 99999 --ld-window-r2 0 "
        f"--out {output_base}"
    )

    try:
        logging.info(f"[LD] Commande exécutée : {cmd}")
        subprocess.run(cmd, shell=True, check=True)
        logging.info("[LD] Calcul r² terminé avec succès.")
    except subprocess.CalledProcessError as e:
        logging.error(f"[LD] Erreur d'exécution PLINK : {e}")
        raise

    ld_file = output_base + ".ld"
    if os.path.exists(ld_file):
        logging.info("[LD] Fichier LD détecté, calcul de D en cours...")
        df = pd.read_csv(ld_file, sep=r'\s+')

        expected = {"R2", "MAF_A", "MAF_B"}
        actual = set(df.columns)
        dprime_col = "Dprime" if "Dprime" in df.columns else ("DP" if "DP" in df.columns else None)

        missing = expected - actual
        if dprime_col is None:
            missing.add("Dprime")

        if not missing:
            def compute_D(row):
                p = row['MAF_A']
                q = row['MAF_B']
                Dmax = min(p * (1 - q), q * (1 - p))
                return row[dprime_col] * Dmax

            df["D"] = df.apply(compute_D, axis=1)
            df.to_csv(output_base + "_with_D.csv", index=False)
            logging.info("[LD] Fichier enrichi avec D sauvegardé.")
        else:
            logging.warning(f"[LD] Colonnes manquantes pour le calcul de D : {missing}")
    else:
        logging.warning("[LD] Aucun fichier .ld généré par PLINK.")