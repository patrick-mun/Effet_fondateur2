# roh.py
"""
Module ROH : Détection des segments d'homozygotie (ROH) avec PLINK
"""
import os
import subprocess
import logging

def run_roh(plink_prefix: str, output_dir: str):
    """
    Exécute la détection des segments ROH via PLINK

    :param plink_prefix: chemin du fichier PLINK (sans extension)
    :param output_dir: dossier de sortie pour les résultats
    """
    os.makedirs(output_dir, exist_ok=True)
    roh_path = os.path.join(output_dir, "roh")

    cmd = f"plink --file {plink_prefix} --homozyg --out {roh_path}"

    try:
        logging.info(f"[ROH] Commande exécutée : {cmd}")
        subprocess.run(cmd, shell=True, check=True)
        logging.info("[ROH] Analyse terminée avec succès.")
    except subprocess.CalledProcessError as e:
        logging.error(f"[ROH] Erreur lors de l'exécution : {e}")
        raise
