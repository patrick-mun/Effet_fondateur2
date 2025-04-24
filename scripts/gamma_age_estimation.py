# gamma_age_estimation.py
"""
Estimation de l'âge d'une mutation basée sur les longueurs de segments ROH
et les fréquences alléliques (approche via la loi gamma).

Ce module applique deux modèles :
1. Modèle indépendant : les longueurs des deux bras du segment sont considérées comme indépendantes
2. Modèle corrélé : on considère une corrélation entre les deux bras

Ajoute également un résumé PDF des graphiques générés et un fichier texte de synthèse.
"""
import numpy as np
from scipy.stats import gamma
from scipy.optimize import minimize
import logging
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import os

logging.basicConfig(filename="data/output/complex_simulation/gamma/debug_gamma.log", level=logging.DEBUG, 
                    format='%(asctime)s - %(levelname)s - %(message)s')

def compute_correction_factor(n):
    return 1 / (2 * n)

def estimate_tau(lengths, shape=2.0):
    mean_l = np.mean(lengths)
    rate = shape / mean_l
    ci_lower = gamma.ppf(0.025, a=shape, scale=1/rate)
    ci_upper = gamma.ppf(0.975, a=shape, scale=1/rate)
    return shape / rate, (ci_lower, ci_upper)

def negative_log_likelihood(params, lengths1, lengths2):
    tau, rho = params
    if tau <= 0 or not (0 < rho < 1):
        return np.inf
    sigma2 = 1 / tau
    cov = rho * sigma2
    var_matrix = np.array([[sigma2, cov], [cov, sigma2]])
    try:
        inv_var = np.linalg.inv(var_matrix)
    except np.linalg.LinAlgError:
        return np.inf
    log_det = np.linalg.slogdet(var_matrix)[1]
    residuals = np.stack([lengths1 - np.mean(lengths1), lengths2 - np.mean(lengths2)], axis=1)
    quad_form = np.sum(residuals @ inv_var * residuals)
    return 0.5 * (len(lengths1) * log_det + quad_form)

def plot_gamma_diagnostics(l_lengths, r_lengths, total_lengths, tau_ind, tau_corr):
    output_dir = "data/output/complex_simulation/gamma"
    os.makedirs(output_dir, exist_ok=True)

    pdf_path = os.path.join(output_dir, "gamma_diagnostics.pdf")
    with PdfPages(pdf_path) as pdf:
        # Histogrammes des longueurs gauche/droite
        plt.figure()
        plt.hist(l_lengths, bins=50, alpha=0.5, label='Gauche', color='skyblue')
        plt.hist(r_lengths, bins=50, alpha=0.5, label='Droite', color='salmon')
        plt.legend()
        plt.title("Distribution des longueurs ROH gauche/droite")
        plt.xlabel("Longueur (cM)")
        plt.ylabel("Nombre d'individus")
        plt.grid(True)
        pdf.savefig()
        plt.close()

        # Histogramme des longueurs totales
        plt.figure()
        plt.hist(total_lengths, bins=50, color='purple')
        plt.title("Distribution des longueurs totales corrigées")
        plt.xlabel("L_total (cM)")
        plt.ylabel("Nombre d'individus")
        plt.axvline(np.mean(total_lengths), color='red', linestyle='--', label=f"Moyenne")
        plt.legend()
        plt.grid(True)
        pdf.savefig()
        plt.close()

        # Scatterplot des longueurs gauche vs droite
        plt.figure()
        plt.scatter(l_lengths, r_lengths, alpha=0.6, label="Individus")
        plt.plot([0, max(l_lengths.max(), r_lengths.max())], [0, max(l_lengths.max(), r_lengths.max())], 'k--', label="y = x")
        plt.title("Corrélation des longueurs gauche vs droite")
        plt.xlabel("Longueur gauche (cM)")
        plt.ylabel("Longueur droite (cM)")
        plt.legend()
        plt.grid(True)
        pdf.savefig()
        plt.close()

        # Annotation des estimations de tau
        plt.figure(figsize=(6, 2))
        plt.axis('off')
        plt.text(0.1, 0.5, f"Estimation τ (indépendant): {tau_ind:.2f} générations", fontsize=12)
        plt.text(0.1, 0.3, f"Estimation τ (corrélé): {tau_corr:.2f} générations", fontsize=12)
        pdf.savefig()
        plt.close()

def write_gamma_summary(n, l_lengths, r_lengths, total_lengths, tau_i, ci_i, tau_c, rho):
    summary_path = "data/output/complex_simulation/gamma/gamma_summary.txt"
    with open(summary_path, "w") as f:
        f.write("Résumé des longueurs ROH et estimations de l'âge\n")
        f.write("============================================\n")
        f.write(f"Nombre d'individus analysés : {n}\n")
        f.write(f"Moyenne longueurs gauche : {np.mean(l_lengths):.2f} cM\n")
        f.write(f"Moyenne longueurs droite : {np.mean(r_lengths):.2f} cM\n")
        f.write(f"Moyenne longueurs totales corrigées : {np.mean(total_lengths):.2f} cM\n")
        f.write("\n")
        f.write("[Modèle indépendant]\n")
        f.write(f"τ = {tau_i:.2f}\n")
        f.write(f"IC à 95% = ({ci_i[0]:.2f}, {ci_i[1]:.2f})\n")
        f.write("\n")
        f.write("[Modèle corrélé]\n")
        f.write(f"τ = {tau_c:.2f}, ρ = {rho:.2f}\n")

def estimate_mutation_age(l_lengths, r_lengths, median_freq):
    try:
        n = min(len(l_lengths), len(r_lengths))
        l_lengths = np.array(l_lengths[:n])
        r_lengths = np.array(r_lengths[:n])

        logging.debug(f"Taille des vecteurs tronqués : {n}")
        logging.debug(f"Longueurs gauche (extrait) : {l_lengths[:5]}")
        logging.debug(f"Longueurs droite (extrait) : {r_lengths[:5]}")

        cs_correction = compute_correction_factor(n)
        total_lengths = l_lengths + r_lengths - 2 * cs_correction
        logging.debug(f"Correction conditionnelle : {cs_correction}")
        logging.debug(f"Longueurs totales après correction (extrait) : {total_lengths[:5]}")

        # Modèle indépendant
        tau_ind, ci_ind = estimate_tau(total_lengths)
        logging.debug(f"[Indépendant] tau = {tau_ind}, CI = {ci_ind}")

        # Modèle corrélé
        res = minimize(negative_log_likelihood, x0=[tau_ind, 0.5], args=(l_lengths, r_lengths), bounds=[(0.01, None), (0.01, 0.99)])
        if res.success:
            tau_corr, rho = res.x
            logging.debug(f"[Corrélé] tau = {tau_corr}, rho = {rho}")

            plot_gamma_diagnostics(l_lengths, r_lengths, total_lengths, tau_ind, tau_corr)
            write_gamma_summary(n, l_lengths, r_lengths, total_lengths, tau_ind, ci_ind, tau_corr, rho)

            return {
                "independent": {"tau": tau_ind, "ci": ci_ind},
                "correlated": {"tau": tau_corr, "ci": (None, None), "rho": rho}
            }
        else:
            logging.warning("Optimisation échouée")
            return None
    except Exception as e:
        logging.error(f"Erreur d'estimation : {str(e)}")
        return None
