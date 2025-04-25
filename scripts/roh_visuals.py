# roh_visuals.py
"""
Visualisation des segments ROH individuels pour les individus atteints.
Génère des graphiques par individu montrant les régions ROH sur le chromosome.
"""

import pandas as pd
import matplotlib.pyplot as plt
import os
import logging

def load_roh_segments(roh_file):
    try:
        df = pd.read_csv(roh_file, delim_whitespace=True)
        logging.info(f"✔ Fichier ROH chargé : {roh_file} ({df.shape[0]} segments)")
        return df
    except Exception as e:
        logging.error(f"❌ Erreur chargement ROH : {e}")
        return None

def plot_roh_for_individuals(df_roh, output_dir, group_file, only_status="cas"):
    os.makedirs(output_dir, exist_ok=True)

    try:
        group_df = pd.read_csv(group_file, sep="\t", header=None, names=["IID", "Group"])
    except Exception as e:
        logging.error(f"❌ Erreur chargement fichier groupes : {e}")
        return

    selected_iids = group_df[group_df["Group"].str.lower() == only_status.lower()]["IID"].tolist()

    if not selected_iids:
        logging.warning(f"⚠️ Aucun individu trouvé avec le statut : {only_status}")
        return

    filtered_df = df_roh[df_roh["IID"].isin(selected_iids)]

    for iid in selected_iids:
        indiv_df = filtered_df[filtered_df["IID"] == iid]
        if indiv_df.empty:
            logging.info(f"ℹ️ Aucun segment ROH détecté pour {iid}")
            continue

        plt.figure(figsize=(10, 1.5))
        for _, row in indiv_df.iterrows():
            plt.plot([row["POS1"], row["POS2"]], [0.5, 0.5], lw=6)
        plt.title(f"Segments ROH - Individu : {iid}", fontsize=10)
        plt.xlabel("Position sur le chromosome (bp)", fontsize=9)
        plt.yticks([])
        plt.tight_layout()
        output_path = os.path.join(output_dir, f"roh_{iid}.png")
        plt.savefig(output_path)
        plt.close()
        logging.info(f"🧬 ROH plot sauvegardé : {output_path}")
