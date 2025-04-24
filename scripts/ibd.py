# ibd.py
"""
Module IBD : détection des relations de parenté génétiques via KING ou PLINK
"""
import os
import subprocess
import logging
import shutil

def run_ibd_king(bed_prefix: str, output_dir: str):
    """
    Exécute KING pour détecter les relations IBD

    :param bed_prefix: chemin du fichier binaire PLINK (sans extension)
    :param output_dir: dossier où stocker les résultats
    """
    os.makedirs(output_dir, exist_ok=True)
    output_prefix = os.path.join(output_dir, "ibd_king")

    # Création d'un fichier .fam temporaire compatible KING (statuts 1 ou 2 uniquement)
    fam_path = f"{bed_prefix}.fam"
    fam_temp_path = f"{bed_prefix}_king.fam"
    try:
        with open(fam_path, "r") as f_in, open(fam_temp_path, "w") as f_out:
            for line in f_in:
                parts = line.strip().split()
                if len(parts) >= 6:
                    if parts[5] not in ["1", "2"]:
                        parts[5] = "2"  # Remplace 3 (porteur) ou autres par 2 (atteint)
                    f_out.write(" ".join(parts) + "\n")

        # Remplacer temporairement le .fam
        fam_backup = f"{fam_path}.bak"
        shutil.copy(fam_path, fam_backup)
        shutil.move(fam_temp_path, fam_path)

        cmd = f"king -b {bed_prefix}.bed --kinship --prefix {output_prefix}"
        logging.info(f"[IBD-KING] Commande exécutée : {cmd}")
        subprocess.run(cmd, shell=True, check=True)
        logging.info("[IBD-KING] Analyse KING terminée avec succès.")
    except subprocess.CalledProcessError as e:
        logging.error(f"[IBD-KING] Erreur lors de l'exécution : {e}")
        raise
    finally:
        # Restauration du .fam original
        if os.path.exists(fam_backup):
            shutil.move(fam_backup, fam_path)


def run_ibd_plink(bed_prefix: str, output_dir: str):
    """
    Exécute PLINK pour détecter les relations IBD (pi-hat)

    :param bed_prefix: chemin du fichier binaire PLINK (sans extension)
    :param output_dir: dossier où stocker les résultats
    """
    os.makedirs(output_dir, exist_ok=True)
    output_prefix = os.path.join(output_dir, "ibd_plink")

    cmd = f"plink --bfile {bed_prefix} --genome --out {output_prefix}"

    try:
        logging.info(f"[IBD-PLINK] Commande exécutée : {cmd}")
        subprocess.run(cmd, shell=True, check=True)
        logging.info("[IBD-PLINK] Analyse PLINK terminée avec succès.")
    except subprocess.CalledProcessError as e:
        logging.error(f"[IBD-PLINK] Erreur lors de l'exécution : {e}")
        raise
