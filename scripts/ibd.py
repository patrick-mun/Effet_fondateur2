# ibd.py
"""
Module IBD : détection des relations de parenté génétiques via KING ou PLINK
"""
import os
import subprocess
import logging
import shutil
import pandas as pd
import networkx as nx
import matplotlib.pyplot as plt

def run_ibd_king(bed_prefix: str, output_dir: str):
    os.makedirs(output_dir, exist_ok=True)
    output_prefix = os.path.join(output_dir, "ibd_king")

    fam_path = f"{bed_prefix}.fam"
    fam_temp_path = f"{bed_prefix}_king.fam"
    try:
        with open(fam_path, "r") as f_in, open(fam_temp_path, "w") as f_out:
            for line in f_in:
                parts = line.strip().split()
                if len(parts) >= 6:
                    if parts[5] not in ["1", "2"]:
                        parts[5] = "2"
                    f_out.write(" ".join(parts) + "\n")

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
        if os.path.exists(fam_backup):
            shutil.move(fam_backup, fam_path)

def run_ibd_plink(bed_prefix: str, output_dir: str):
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

def summarize_ibd_king(output_dir: str):
    try:
        kin_file = os.path.join(output_dir, "ibd_king.kin")
        if not os.path.exists(kin_file):
            logging.error(f"[IBD-Summary] Fichier {kin_file} introuvable.")
            return

        df = pd.read_csv(kin_file, delim_whitespace=True)
        df_summary = df[["ID1", "ID2", "InfType", "Kinship"]]

        summary_path = os.path.join(output_dir, "ibd_king_summary.csv")
        df_summary.to_csv(summary_path, index=False)
        logging.info(f"[IBD-Summary] Tableau récapitulatif sauvegardé : {summary_path}")

    except Exception as e:
        logging.error(f"[IBD-Summary] Erreur : {e}")

def plot_ibd_network(output_dir: str, kinship_threshold: float = 0.05):
    try:
        kin_file = os.path.join(output_dir, "ibd_king.kin")
        if not os.path.exists(kin_file):
            logging.error(f"[IBD-Network] Fichier {kin_file} introuvable.")
            return

        df = pd.read_csv(kin_file, delim_whitespace=True)
        df = df[df["Kinship"] >= kinship_threshold]

        G = nx.Graph()
        for _, row in df.iterrows():
            G.add_edge(row["ID1"], row["ID2"], weight=row["Kinship"])

        plt.figure(figsize=(12, 8))
        pos = nx.spring_layout(G, seed=42)
        nx.draw(G, pos, with_labels=True, node_size=700, font_size=10, edge_color='gray')
        plt.title("Réseau de parenté basé sur KING")
        plt.tight_layout()
        plt.savefig(os.path.join(output_dir, "ibd_king_network.png"))
        plt.close()
        logging.info("[IBD-Network] Graphe de parenté sauvegardé.")

    except Exception as e:
        logging.error(f"[IBD-Network] Erreur : {e}")
