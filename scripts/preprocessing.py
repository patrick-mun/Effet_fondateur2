# preprocessing.py
"""
Module de prétraitement :
--------------------------
Ce module inclut les étapes de filtrage de base et la séparation des individus en groupes (cas/témoins).
Utilise PLINK pour effectuer les opérations sur les données génétiques.
Génère également les fichiers binaires (.bed/.bim/.fam), .frq, .hwe, .log, et les sous-ensembles cas/témoins à tous les formats nécessaires pour KING, ADMIXTURE, Gamma, etc.
Vérifie en début la présence et le format des fichiers d'entrée, puis vérifie en sortie la création effective des fichiers.
Génère également des graphiques QC (fréquences alléliques, HWE).
"""
import os
import subprocess
import logging
import csv
import pandas as pd
import matplotlib.pyplot as plt

def check_input_files_exist(input_prefix):
    required = [f"{input_prefix}.ped", f"{input_prefix}.map"]
    missing = [f for f in required if not os.path.exists(f)]
    if missing:
        for f in missing:
            logging.critical(f"❌ Fichier d'entrée manquant : {f}")
        raise FileNotFoundError("Un ou plusieurs fichiers d'entrée sont absents.")
    logging.info("✔ Tous les fichiers d'entrée sont présents.")

def validate_raw_format(raw_path):
    try:
        df = pd.read_csv(raw_path, sep=r'\s+', nrows=10)
    except Exception as e:
        logging.critical(f"Erreur lecture {raw_path} : {e}")
        raise

    required_cols = ["FID", "IID", "SEX", "PHENOTYPE"]
    if not all(col in df.columns for col in required_cols):
        logging.critical("Colonnes obligatoires manquantes dans le .raw")
        raise ValueError("Format .raw invalide")
    snps = [col for col in df.columns if col not in required_cols]
    if len(snps) < 10:
        logging.warning("Moins de 10 SNPs détectés dans le fichier .raw (attention au filtrage)")
    logging.info(f"✔ Format .raw valide ({len(snps)} SNPs détectés)")

def clean_ped_pheno(path):
    with open(path, "r") as infile:
        lines = infile.readlines()

    cleaned_lines = []
    for line in lines:
        parts = line.strip().split()
        if len(parts) >= 6:
            if parts[5] not in ["1", "2", "0"]:
                parts[5] = "0"
            cleaned_lines.append("\t".join(parts) + "\n")

    with open(path, "w") as outfile:
        outfile.writelines(cleaned_lines)
    logging.info("Colonne PHENOTYPE nettoyée (1/2/0 uniquement).")

def create_ped_porteurs(path_in, path_out):
    with open(path_in, "r") as f:
        lines = f.readlines()

    modified = []
    for line in lines:
        parts = line.strip().split()
        if len(parts) >= 6:
            pheno = parts[5]
            if pheno == "1":
                parts[5] = "1"
            elif pheno == "3":
                parts[5] = "2"
            else:
                parts[5] = "0"
            modified.append(" ".join(parts) + "\n")

    with open(path_out, "w") as f:
        f.writelines(modified)
    logging.info(f"✅ Fichier .ped porteurs généré à : {path_out}")

def check_outputs(file_paths, output_csv):
    results = []
    for f in file_paths:
        exists = os.path.exists(f)
        results.append({"fichier": f, "existe": exists})
        if exists:
            logging.info(f"✔ Fichier présent : {f}")
        else:
            logging.warning(f"❌ Fichier manquant : {f}")

    with open(output_csv, "w", newline="") as csvfile:
        writer = csv.DictWriter(csvfile, fieldnames=["fichier", "existe"])
        writer.writeheader()
        writer.writerows(results)

def generate_qc_plots(output_prefix):
    """
    Génère les graphiques de qualité (distribution des MAF et des p-values HWE)
    et produit un résumé statistique des résultats ALL(QT) dans un fichier CSV.
    """
    freq_path = f"{output_prefix}/allele_frequencies.frq"
    if os.path.exists(freq_path):
        df = pd.read_csv(freq_path, sep=r'\s+')
        if "MAF" in df.columns:
            plt.figure()
            plt.hist(df["MAF"], bins=50, color="skyblue", edgecolor="k")
            plt.title("Distribution des MAF (Minor Allele Frequencies)")
            plt.xlabel("Fréquence allélique mineure")
            plt.ylabel("Nombre de SNPs")
            plt.grid(True)
            plt.tight_layout()
            plt.savefig(f"{output_prefix}/qc_maf_distribution.png")
            plt.close()
            logging.info("Graphique MAF sauvegardé.")

    hwe_path = f"{output_prefix}/hwe.hwe"
    if os.path.exists(hwe_path):
        df = pd.read_csv(hwe_path, sep=r'\s+')
        if "P" in df.columns:
            plt.figure()
            df_p = df[df["TEST"] == "ALL"]
            plt.hist(df_p["P"], bins=50, color="salmon", edgecolor="k")
            plt.title("Distribution des p-values HWE (test ALL)")
            plt.xlabel("p-value HWE")
            plt.ylabel("Nombre de SNPs")
            plt.grid(True)
            plt.tight_layout()
            plt.savefig(f"{output_prefix}/qc_hwe_pvalues.png")
            plt.close()
            logging.info("Graphique HWE sauvegardé.")

            # Résumé statistique des p-values HWE (ALL/QT)
            summary = pd.DataFrame({
                "Statistique": [
                    "Nombre total de SNPs testés (ALL)",
                    "SNPs avec p < 0.05",
                    "SNPs avec p < 0.001",
                    "SNPs avec p >= 0.05",
                    "p-value médiane",
                    "p-value minimale",
                    "p-value maximale"
                ],
                "Valeur": [
                    len(df_p),
                    (df_p["P"] < 0.05).sum(),
                    (df_p["P"] < 0.001).sum(),
                    (df_p["P"] >= 0.05).sum(),
                    df_p["P"].median(),
                    df_p["P"].min(),
                    df_p["P"].max()
                ]
            })
            summary.to_csv(f"{output_prefix}/hwe_all_summary.csv", index=False)
            logging.info("Résumé statistique HWE (ALL) sauvegardé dans hwe_all_summary.csv")
          
         


def run_preprocessing(input_prefix: str, output_prefix: str, cas_file: str, temoins_file: str):
    os.makedirs(output_prefix, exist_ok=True)
    check_input_files_exist(input_prefix)

    filtrage_cmd = f"plink --file {input_prefix} --maf 0.01 --geno 0.05 --mind 0.1 --hwe 1e-6 --recode A --out {output_prefix}/filtered_data"
    export_ped_cmd = f"plink --file {input_prefix} --maf 0.01 --geno 0.05 --mind 0.1 --hwe 1e-6 --recode --out {output_prefix}/filtered_data"
    binarisation_cmd = f"plink --file {input_prefix} --maf 0.01 --geno 0.05 --mind 0.1 --hwe 1e-6 --make-bed --out {output_prefix}/filtered_data"
    freq_cmd = f"plink --file {output_prefix}/filtered_data --freq --out {output_prefix}/allele_frequencies"
    hwe_cmd = f"plink --file {output_prefix}/filtered_data --hardy midp --out {output_prefix}/hwe"

    separation_cmds = [
        f"plink --file {output_prefix}/filtered_data --keep {cas_file} --recode A --out {output_prefix}/geno_cas",
        f"plink --file {output_prefix}/filtered_data --keep {temoins_file} --recode A --out {output_prefix}/geno_temoins",
        f"plink --file {output_prefix}/filtered_data --keep {cas_file} --make-bed --out {output_prefix}/geno_cas",
        f"plink --file {output_prefix}/filtered_data --keep {temoins_file} --make-bed --out {output_prefix}/geno_temoins",
        f"plink --file {output_prefix}/filtered_data --keep {cas_file} --recode --out {output_prefix}/geno_cas",
        f"plink --file {output_prefix}/filtered_data --keep {temoins_file} --recode --out {output_prefix}/geno_temoins"
    ]

    try:
        for cmd in [filtrage_cmd, export_ped_cmd, binarisation_cmd, freq_cmd, hwe_cmd] + separation_cmds:
            subprocess.run(cmd, shell=True, check=True)

        clean_ped_pheno(f"{output_prefix}/filtered_data.ped")
        create_ped_porteurs(
            f"{output_prefix}/filtered_data.ped",
            f"{output_prefix}/filtered_data_porteurs.ped"
        )

        vcf_cmd = f"plink --bfile {output_prefix}/filtered_data --recode vcf --out {output_prefix}/filtered_data"
        subprocess.run(vcf_cmd, shell=True, check=True)

        porteurs_cmd = f"plink --file {output_prefix}/filtered_data_porteurs --hardy midp --out {output_prefix}/hwe_porteurs"
        subprocess.run(porteurs_cmd, shell=True, check=True)
        logging.info("[Preprocessing] Test HWE porteurs effectué avec succès.")

        validate_raw_format(f"{output_prefix}/filtered_data.raw")

    except subprocess.CalledProcessError as e:
        logging.error(f"Erreur durant le prétraitement : {e}")
        return

    expected_files = [
        f"{output_prefix}/filtered_data.raw",
        f"{output_prefix}/filtered_data.ped",
        f"{output_prefix}/filtered_data.map",
        f"{output_prefix}/filtered_data.bed",
        f"{output_prefix}/filtered_data.bim",
        f"{output_prefix}/filtered_data.fam",
        f"{output_prefix}/filtered_data.vcf",
        f"{output_prefix}/allele_frequencies.frq",
        f"{output_prefix}/hwe.hwe",
        f"{output_prefix}/hwe_porteurs.hwe",
        f"{output_prefix}/geno_cas.raw",
        f"{output_prefix}/geno_temoins.raw",
        f"{output_prefix}/geno_cas.ped",
        f"{output_prefix}/geno_temoins.ped",
        f"{output_prefix}/geno_cas.map",
        f"{output_prefix}/geno_temoins.map",
        f"{output_prefix}/geno_cas.bed",
        f"{output_prefix}/geno_temoins.bed",
        f"{output_prefix}/geno_cas.bim",
        f"{output_prefix}/geno_temoins.bim",
        f"{output_prefix}/geno_cas.fam",
        f"{output_prefix}/geno_temoins.fam",
        f"{output_prefix}/filtered_data_porteurs.ped"
    ]
    summary_csv = os.path.join(output_prefix, "summary_files.csv")
    check_outputs(expected_files, summary_csv)
    generate_qc_plots(output_prefix)
