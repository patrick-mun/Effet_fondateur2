# gamma.py
import os
import subprocess
import logging
import pandas as pd
import matplotlib.pyplot as plt

def prepare_gamma_input(raw_path: str, map_path: str, output_txt: str, output_fig: str, include_extreme_freqs: bool = True):
    """
    Prépare les fichiers d'entrée pour Gamma à partir des fichiers .raw et .map.

    Les fichiers requis sont :
    - raw_path : un fichier .raw généré par PLINK avec l'option --recodeA (génotypes 0/1/2)
      ex : data/output/complex_simulation/geno/filtered_data.raw
    - map_path : le fichier .map correspondant, filtré par PLINK avec les mêmes SNPs
      ex : data/output/complex_simulation/geno/filtered_data.map
    ⚠️ Ces fichiers doivent provenir de PLINK après filtrage (ex : MAF, GENO, HWE) pour assurer que
    les SNPs présents dans .raw et .map soient cohérents. Ce script s'assure de ne garder que
    les SNPs communs aux deux fichiers.
    """
    try:
        raw_df = pd.read_csv(raw_path, sep=r'\s+', dtype=str)
        raw_df.columns = [col.split('_')[0] if '_' in col else col for col in raw_df.columns]

        for col in raw_df.columns:
            try:
                raw_df[col] = pd.to_numeric(raw_df[col])
            except Exception:
                pass

        map_df = pd.read_csv(map_path, sep=r'\s+', header=None, names=["CHR", "SNP", "CM", "BP"])

        intersection = map_df["SNP"].isin(raw_df.columns)
        if not intersection.any():
            print("⚠️ Aucun SNP de la carte ne correspond aux colonnes du fichier .raw.")
            print("Vérifiez que les noms de SNPs sont identiques dans .map et .raw.")
        else:
            print(f"✔ {intersection.sum()} SNPs communs trouvés entre .map et .raw.")
        map_df = map_df[intersection]

        if map_df["CM"].nunique() <= 1:
            map_df["CM"] = map_df["BP"] / 1_000_000  # 1 cM ≈ 1 Mb (approximation)
            updated_map_path = map_path.replace(".map", "_cm_updated.map")
            map_df.to_csv(updated_map_path, sep="\t", index=False, header=False)
            print(f"\n⚠️ La colonne CM a été recalculée automatiquement. Un fichier mis à jour a été enregistré ici : {updated_map_path}\n")
            confirmation = input("Souhaitez-vous continuer avec ce fichier ? (oui/non) : ").strip().lower()
            if confirmation != "oui":
                print("Opération annulée par l'utilisateur.")
                exit(1)

    except Exception as e:
        logging.error(f"Erreur lecture fichiers .raw ou .map : {e}")
        raise

    snp_cols = [col for col in raw_df.columns if col not in ["FID", "IID", "PAT", "MAT", "SEX", "PHENOTYPE"]]
    data = []

    for snp in snp_cols:
        if snp not in map_df["SNP"].values:
            continue
        cm_pos = map_df.loc[map_df["SNP"] == snp, "CM"].values[0]
        try:
            values = pd.to_numeric(raw_df[snp], errors="coerce").dropna()
            valid = values[values.isin([0, 1, 2])]
            if len(valid) == 0:
                continue
            freq = valid.sum() / (2 * len(valid))
            if include_extreme_freqs or (0 < freq < 1):
                data.append((snp, cm_pos, freq))
        except Exception as e:
            logging.warning(f"Impossible de calculer la fréquence pour {snp} : {e}")
            continue

    if not data:
        logging.error("Aucun SNP valide trouvé pour Gamma.")
        return

    df_out = pd.DataFrame(data, columns=["SNP", "Position_cM", "Freq"])
    df_out = df_out.dropna()
    df_out.to_csv(output_txt, sep="\t", index=False)

    plt.figure(figsize=(10, 4))
    plt.scatter(df_out["Position_cM"], df_out["Freq"], alpha=0.6)
    plt.xlabel("Position (cM)")
    plt.ylabel("Fréquence")
    plt.title("Fréquences alléliques retenues pour Gamma")
    plt.tight_layout()
    plt.savefig(output_fig)
    plt.close()
    logging.info("[Gamma] Fichier d'entrée et graphique générés avec succès.")

def run_gamma(input_txt: str, output_txt: str):
    cmd = f"Gamma -i {input_txt} -o {output_txt}"
    try:
        logging.info(f"[Gamma] Exécution : {cmd}")
        subprocess.run(cmd, shell=True, check=True)
        logging.info("[Gamma] Exécution terminée.")
    except subprocess.CalledProcessError as e:
        logging.error(f"[Gamma] Erreur exécution : {e}")
        raise