# roh_overlap.py
"""
Visualisation comparative des segments ROH entre individus d'un groupe.
Chaque segment est affiché sur une même ligne horizontale, facilitant la détection de recouvrements.
"""

import os
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import logging

def plot_roh_overlap(df_roh, output_path, group_file, target_group="ATTEINT"):
    try:
        group_df = pd.read_csv(group_file, sep="\t", header=None, names=["IID", "Group"])
        selected_iids = group_df[group_df["Group"].str.upper() == target_group.upper()]["IID"].tolist()

        filtered_df = df_roh[df_roh["IID"].isin(selected_iids)]

        if filtered_df.empty:
            logging.warning(f"[ROH Overlap] Aucun segment ROH trouvé pour le groupe {target_group}.")
            return

        # Tri pour affichage ordonné
        filtered_df["IID"] = pd.Categorical(filtered_df["IID"], categories=selected_iids, ordered=True)
        filtered_df = filtered_df.sort_values(["IID", "POS1"])

        # Calcul correct de la zone de recouvrement stricte
        overlap_start = filtered_df["POS1"].max()
        overlap_end = filtered_df["POS2"].min()

        plt.figure(figsize=(12, len(selected_iids)*0.5 + 1))
        for i, (iid, group) in enumerate(filtered_df.groupby("IID")):
            for _, row in group.iterrows():
                plt.hlines(y=i, xmin=row["POS1"], xmax=row["POS2"], color="tab:blue", lw=5)

        plt.yticks(range(len(selected_iids)), selected_iids)
        plt.xlabel("Position sur le chromosome (bp)")
        plt.title(f"Segments ROH - Groupe : {target_group}")

        if overlap_start < overlap_end:
            plt.axvspan(overlap_start, overlap_end, color="lightcoral", alpha=0.5, label="Zone commune")
            plt.legend()
            # Sauvegarde de la zone commune
            overlap_file = os.path.join(os.path.dirname(output_path), "zone_commune.txt")
            with open(overlap_file, "w") as f:
                f.write(f"{overlap_start}\t{overlap_end}\n")
            logging.info(f"[ROH Overlap] Zone commune sauvegardée : {overlap_file}")
        else:
            logging.warning("[ROH Overlap] Pas de véritable recouvrement entre tous les individus.")

        plt.tight_layout()
        os.makedirs(os.path.dirname(output_path), exist_ok=True)
        plt.savefig(output_path)
        plt.close()
        logging.info(f"[ROH Overlap] Graphique sauvegardé : {output_path}")

    except Exception as e:
        logging.error(f"[ROH Overlap] Erreur lors de la génération du graphique : {e}")

