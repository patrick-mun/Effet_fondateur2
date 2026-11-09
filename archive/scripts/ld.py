# ld.py
"""
Module LD : calcul du linkage disequilibrium (r² et D), et analyse dans les régions ROH.
"""
import os
import subprocess
import logging
import pandas as pd
import matplotlib.pyplot as plt
import json

def run_ld(plink_prefix: str, output_dir: str):
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

def analyse_ld_in_roh(ld_csv_path, roh_file_path, map_file_path, fam_file_path, output_dir, r2_threshold=0.1):
    if not os.path.exists(ld_csv_path) or not os.path.exists(roh_file_path) or not os.path.exists(map_file_path):
        print("[LD vs ROH] Fichier introuvable.")
        return

    df_ld = pd.read_csv(ld_csv_path)
    df_roh = pd.read_csv(roh_file_path, sep=r'\s+')
    df_map = pd.read_csv(map_file_path, sep=r'\s+', header=None, names=["CHR", "SNP", "CM", "BP"])

    snp_pos = df_map.set_index("SNP")["BP"].to_dict()
    df_ld["BP_A"] = df_ld["SNP_A"].map(snp_pos)
    df_ld["BP_B"] = df_ld["SNP_B"].map(snp_pos)

    if df_ld["BP_A"].isnull().all() or df_ld["BP_B"].isnull().all():
        print("[LD vs ROH] Aucune position valide trouvée pour BP_A ou BP_B à partir de la map.")
        return

    if "CHR" not in df_ld.columns:
        if "CHR_A" in df_ld.columns:
            df_ld["CHR"] = df_ld["CHR_A"]
        else:
            df_ld["CHR"] = 1

    def in_same_roh(row):
        if pd.isnull(row["BP_A"]) or pd.isnull(row["BP_B"]):
            return pd.NA
        for _, roh in df_roh.iterrows():
            if (
                roh["CHR"] == row["CHR"]
                and row["BP_A"] >= roh["POS1"]
                and row["BP_B"] <= roh["POS2"]
            ):
                return True
        return False

    try:
        df_ld["in_roh"] = df_ld.apply(in_same_roh, axis=1)

        in_roh_mask = df_ld["in_roh"] == True
        out_roh_mask = df_ld["in_roh"] == False

        n_total = len(df_ld)
        n_in_roh = in_roh_mask.sum()
        n_out_roh = out_roh_mask.sum()
        print(f"[LD vs ROH] Paires totales : {n_total}, dans ROH : {n_in_roh}, hors ROH : {n_out_roh}")

        plt.figure(figsize=(10, 5))
        plt.hist(df_ld[in_roh_mask]["R2"], bins=100, alpha=0.6, label="r² dans ROH", color="mediumseagreen")
        plt.hist(df_ld[out_roh_mask]["R2"], bins=100, alpha=0.6, label="r² hors ROH", color="darkorange")
        plt.xlabel("Valeur de r²")
        plt.ylabel("Nombre de paires de SNPs")
        plt.title("Comparaison de la distribution de r² dans et hors ROH")
        plt.xlim(0, 0.2)
        plt.legend()
        plt.tight_layout()
        os.makedirs(output_dir, exist_ok=True)
        plt.savefig(os.path.join(output_dir, "ld_r2_comparaison_roh.png"))
        plt.close()

        df_ld.to_csv(os.path.join(output_dir, "ld_with_roh_annotation.csv"), index=False)

        mean_r2_in = df_ld[in_roh_mask]["R2"].mean()
        mean_r2_out = df_ld[out_roh_mask]["R2"].mean()

        df_mut = df_map[df_map["SNP"].str.startswith("rsMUT")]
        mutation_position = df_mut.iloc[0]["BP"] if not df_mut.empty else None
        mutation_chr = df_mut.iloc[0]["CHR"] if not df_mut.empty else None
        mutation_in_roh = False
        if mutation_position and mutation_chr:
            for _, roh in df_roh.iterrows():
                if roh["CHR"] == mutation_chr and roh["POS1"] <= mutation_position <= roh["POS2"]:
                    mutation_in_roh = True
                    break

        stats = {
            "Analyse r² ROH": {
                "Moyenne r² dans ROH": round(mean_r2_in, 4),
                "Moyenne r² hors ROH": round(mean_r2_out, 4),
                "Delta": round(mean_r2_in - mean_r2_out, 4),
                "Mutation localisée dans ROH": mutation_in_roh
            }
        }

        with open(os.path.join(output_dir, "ld_r2_stats.json"), "w", encoding="utf-8") as f:
            json.dump(stats, f, indent=4, ensure_ascii=False)

    except Exception as e:
        print(f"[LD vs ROH] Erreur lors de l'annotation in_roh : {e}")
        return
def ld_mutation_zone_analysis(plink_prefix: str, fam_path: str, map_path: str, output_dir: str):
    os.makedirs(output_dir, exist_ok=True)

    cas_file = os.path.join(output_dir, "cas.txt")
    with open(fam_path, "r") as f_in, open(cas_file, "w") as f_out:
        for line in f_in:
            cols = line.strip().split()
            if len(cols) >= 6 and cols[5] == "2":
                f_out.write(f"{cols[0]}\t{cols[1]}\n")

    # Lecture de la position de la mutation depuis .map
    df_map = pd.read_csv(map_path, sep=r'\s+', header=None, names=["CHR", "SNP", "CM", "BP"])
    df_mut = df_map[df_map["SNP"].notna() & df_map["SNP"].str.startswith("rsMUT")]
    if df_mut.empty:
        print("[LD] SNP rsMUT non trouvé dans le fichier .map")
        return
    chr = int(df_mut.iloc[0]["CHR"])
    pos = int(df_mut.iloc[0]["BP"])

    # Crée une liste de SNPs dans la zone cible
    start = max(0, pos - 1_000_000)
    end = pos + 1_000_000
    df_zone = df_map[(df_map["CHR"] == chr) & (df_map["BP"] >= start) & (df_map["BP"] <= end)]
    snp_list_path = os.path.join(output_dir, "snp_list.txt")
    df_zone["SNP"].to_csv(snp_list_path, index=False, header=False)

    mutation_zone_prefix = os.path.join(output_dir, "filtered_mutation_zone")
    subprocess.run([
        "plink", "--bfile", plink_prefix,
        "--extract", snp_list_path,
        "--make-bed", "--out", mutation_zone_prefix
    ], check=True)

    ld_output_prefix = os.path.join(output_dir, "ld_cas_mutation_zone")
    subprocess.run([
        "plink", "--bfile", mutation_zone_prefix,
        "--keep", cas_file,
        "--make-founders",
        "--r2", "dprime", "with-freqs",
        "--ld-window-kb", "1000", "--ld-window", "99999", "--ld-window-r2", "0",
        "--out", ld_output_prefix
    ], check=True)

    # Calcul du D et sauvegarde CSV enrichi
    ld_file = ld_output_prefix + ".ld"
    if os.path.exists(ld_file):
        df = pd.read_csv(ld_file, sep=r'\s+')
        expected = {"R2", "MAF_A", "MAF_B"}
        actual = set(df.columns)
        dprime_col = "Dprime" if "Dprime" in df.columns else ("DP" if "DP" in df.columns else None)

        missing = expected - actual
        if dprime_col is None:
            missing.add("Dprime")

        if not missing:
            def compute_D(row):
                p, q = row['MAF_A'], row['MAF_B']
                Dmax = min(p * (1 - q), q * (1 - p))
                return row[dprime_col] * Dmax
            df["D"] = df.apply(compute_D, axis=1)
            df.to_csv(ld_output_prefix + "_with_D.csv", index=False)
            print("[LD Mutation Zone] Fichier enrichi avec D sauvegardé.")
        else:
            print(f"[LD Mutation Zone] Colonnes manquantes : {missing}")

    print("[LD Mutation Zone] LD recalculé pour la région ciblée chez les cas.")
