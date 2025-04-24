# run_pipeline.py
import os
import logging
import pandas as pd
from scripts.preprocessing import run_preprocessing
from scripts.roh import run_roh
from scripts.ibd import run_ibd_king, run_ibd_plink
from scripts.ld import run_ld
from scripts.gamma import prepare_gamma_input
from scripts.gamma_age_estimation import estimate_mutation_age
from scripts.adegenet import generate_adegenet_script, run_adegenet_script
from scripts.reporting import generate_full_report
from scripts.utils import safe_mkdir

logging.basicConfig(level=logging.INFO)

# Dossiers pour simulation complexe
data_dir = "data/input/complex_simulation"
input_prefix = os.path.join(data_dir, "genotype_data")
output_dir = "data/output/complex_simulation"
raw_path = f"{output_dir}/geno/filtered_data.raw"
map_path = f"{output_dir}/geno/filtered_data.map"
cas_file = os.path.join(data_dir, "cas.txt")
temoins_file = os.path.join(data_dir, "temoins.txt")
group_file = os.path.join(data_dir, "groupes.txt")

def main():
    safe_mkdir(output_dir)

    # Prétraitement
    run_preprocessing(input_prefix, f"{output_dir}/geno", cas_file, temoins_file)

    # ROH
    run_roh(input_prefix, f"{output_dir}/roh")

    # IBD
    run_ibd_king(f"{output_dir}/geno/filtered_data", f"{output_dir}/ibd")
    run_ibd_plink(f"{output_dir}/geno/filtered_data", f"{output_dir}/ibd")

    # LD (corrigé pour utiliser les données filtrées)
    run_ld(f"{output_dir}/geno/filtered_data", f"{output_dir}/ld")

    # Gamma preparation
    gamma_dir = f"{output_dir}/gamma"
    gamma_txt = os.path.join(gamma_dir, "gamma_input.txt")
    gamma_fig = os.path.join(gamma_dir, "gamma_plot.png")
    safe_mkdir(gamma_dir)
    prepare_gamma_input(raw_path, map_path, gamma_txt, gamma_fig, include_extreme_freqs=True)

    # Estimation directe en Python
    df_gamma = pd.read_csv(gamma_txt, sep="\t")
    lengths = df_gamma["Position_cM"].values
    median_freq = df_gamma["Freq"].median()
    result = estimate_mutation_age(lengths[:len(lengths)//2], lengths[len(lengths)//2:], median_freq)

    with open(os.path.join(gamma_dir, "estimation_directe.txt"), "w") as f:
        if result is None or "error" in result:
            logging.warning(f"[Gamma] {result.get('error', 'Estimation échouée') if result else 'Estimation échouée'}")
            f.write(f"[Gamma] Erreur : {result.get('error', 'Estimation échouée') if result else 'Estimation échouée'}\n")
        else:
            f.write(f"Modèle indépendant: τ = {result['independent']['tau']:.2f} "
                    f"({result['independent']['ci'][0]:.2f}–{result['independent']['ci'][1]:.2f})\n")
            tau_corr = result['correlated']['tau']
            ci_corr = result['correlated']['ci']
            rho_corr = result['correlated']['rho']
            if tau_corr is not None and ci_corr[0] is not None:
                f.write(f"Modèle corrélé: τ = {tau_corr:.2f} "
                        f"({ci_corr[0]:.2f}–{ci_corr[1]:.2f}), "
                        f"rho = {rho_corr:.2f}\n")
            else:
                f.write("Modèle corrélé: estimation impossible.\n")

    # Adegenet
    script_path = generate_adegenet_script(raw_path, group_file, f"{output_dir}/adegenet")
    run_adegenet_script(script_path)

    # Rapport
    generate_full_report(output_dir, "rapport_final.pdf", "rapport_final.html")

if __name__ == "__main__":
     main()

