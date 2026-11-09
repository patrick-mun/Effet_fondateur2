# gamma_age_estimation.py
"""
Estimation de l'âge d'une mutation basée sur les longueurs de segments ROH
et les fréquences alléliques, en suivant fidèlement l'algorithme original en R (version anglaise).
"""
import numpy as np
from scipy.stats import gamma
import logging
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import os

logging.basicConfig(filename="data/output/complex_simulation/gamma/debug_gamma.log", level=logging.DEBUG,
                    format='%(asctime)s - %(levelname)s - %(message)s')

def compute_cs_correction(median_freq, chrom_length_cM=100, nb_markers=1000):
    e = 0.01
    p = median_freq**2 + (1 - median_freq)**2
    phi = (chrom_length_cM / 100) / nb_markers
    loci = np.log(e) / np.log(p)
    return loci * phi

def estimate_mutation_age(l_lengths, r_lengths, median_freq, n, cc=0.95, chance_share_correction=True):
    l_lengths = np.array(l_lengths[:n]) / 100
    r_lengths = np.array(r_lengths[:n]) / 100

    if chance_share_correction:
        cs_correction = compute_cs_correction(median_freq)
    else:
        cs_correction = 0

    # Independent genealogy
    i_cs_correction = 0 if n < 10 else cs_correction
    length_correction = (np.sum(l_lengths) + np.sum(r_lengths) - 2 * (n - 1) * i_cs_correction) / (2 * n)
    sum_lengths = np.sum(l_lengths) + np.sum(r_lengths) + 2 * length_correction - 2 * (n - 1) * i_cs_correction
    bc_i = (2 * n - 1) / (2 * n)
    tau_i = (bc_i * 2 * n) / sum_lengths
    g_l_i = gamma.ppf((1 - cc) / 2, a=2 * n, scale=1 / (2 * n * bc_i))
    g_u_i = gamma.ppf(cc + (1 - cc) / 2, a=2 * n, scale=1 / (2 * n * bc_i))
    ci_i = (g_l_i * tau_i, g_u_i * tau_i)

    # Correlated genealogy
    length_correction = (np.sum(l_lengths) + np.sum(r_lengths) - 2 * (n - 1) * cs_correction) / (2 * n)
    l_lengths[np.argmax(l_lengths)] += length_correction + cs_correction
    r_lengths[np.argmax(r_lengths)] += length_correction + cs_correction

    lengths = l_lengths + r_lengths - 2 * cs_correction
    mean_L = np.mean(lengths)
    var_L = np.var(lengths, ddof=1)
    rho = (n * mean_L**2 - var_L * (1 + 2 * n)) / (n * mean_L**2 + var_L * (n - 1))
    n_star = n / (1 + (n - 1) * rho)

    if n_star > n:
        n_star = n
    if n_star < -n:
        n_star = -n
    if rho < -2 / (n - 1):
        n_star = n / (1 + (n - 1) * abs(rho))
    elif -2 / (n - 1) <= rho < -1 / (n - 1):
        n_star = n

    bc_c = (2 * n_star - 1) / (2 * n_star)
    tau_c = (bc_c * 2 * n) / np.sum(lengths)
    g_l_c = gamma.ppf((1 - cc) / 2, a=2 * n_star, scale=1 / (2 * n_star * bc_c))
    g_u_c = gamma.ppf(cc + (1 - cc) / 2, a=2 * n_star, scale=1 / (2 * n_star * bc_c))
    ci_c = (g_l_c * tau_c, g_u_c * tau_c)

    print(f"rho = {rho:.3f} | n_star = {n_star:.2f} | tau_ind = {tau_i:.2f} | tau_corr = {tau_c:.2f}")
    plot_gamma_diagnostics(l_lengths, r_lengths, lengths, tau_i, tau_c)
    write_gamma_summary(n, l_lengths, r_lengths, lengths, tau_i, ci_i, tau_c, rho, ci_c)

    return {
        "independent": {"tau": tau_i, "ci": ci_i},
        "correlated": {"tau": tau_c, "ci": ci_c, "rho": rho}
    }

def plot_gamma_diagnostics(l_lengths, r_lengths, total_lengths, tau_ind, tau_corr):
    output_dir = "data/output/complex_simulation/gamma"
    os.makedirs(output_dir, exist_ok=True)
    pdf_path = os.path.join(output_dir, "gamma_diagnostics.pdf")

    with PdfPages(pdf_path) as pdf:
        plt.figure()
        plt.hist(l_lengths, bins=50, alpha=0.5, label='Left', color='skyblue')
        plt.hist(r_lengths, bins=50, alpha=0.5, label='Right', color='salmon')
        plt.legend()
        plt.title("ROH Segment Lengths (Left/Right)")
        plt.xlabel("Length (cM)")
        plt.ylabel("Count")
        plt.grid(True)
        pdf.savefig()
        plt.close()

        plt.figure()
        plt.hist(total_lengths, bins=50, color='purple')
        plt.title("Corrected Total ROH Lengths")
        plt.xlabel("Total Length (cM)")
        plt.ylabel("Count")
        plt.axvline(np.mean(total_lengths), color='red', linestyle='--', label=f"Mean")
        plt.legend()
        plt.grid(True)
        pdf.savefig()
        plt.close()

        plt.figure()
        plt.scatter(l_lengths, r_lengths, alpha=0.6, label="Individuals")
        max_val = max(l_lengths.max(), r_lengths.max())
        plt.plot([0, max_val], [0, max_val], 'k--', label="y = x")
        plt.title("Left vs Right ROH Lengths")
        plt.xlabel("Left (cM)")
        plt.ylabel("Right (cM)")
        plt.legend()
        plt.grid(True)
        pdf.savefig()
        plt.close()

        plt.figure(figsize=(6, 2))
        plt.axis('off')
        plt.text(0.1, 0.5, f"τ (independent): {tau_ind:.2f} generations", fontsize=12)
        plt.text(0.1, 0.3, f"τ (correlated): {tau_corr:.2f} generations", fontsize=12)
        pdf.savefig()
        plt.close()

def write_gamma_summary(n, l_lengths, r_lengths, total_lengths, tau_i, ci_i, tau_c, rho, ci_c):
    summary_path = "data/output/complex_simulation/gamma/gamma_summary.txt"
    with open(summary_path, "w") as f:
        f.write("Summary of ROH lengths and age estimates\n")
        f.write("========================================\n")
        f.write(f"Sample size: {n}\n")
        f.write(f"Mean left length: {np.mean(l_lengths):.2f} cM\n")
        f.write(f"Mean right length: {np.mean(r_lengths):.2f} cM\n")
        f.write(f"Mean total corrected length: {np.mean(total_lengths):.2f} cM\n\n")
        f.write("[Independent model]\n")
        f.write(f"tau = {tau_i:.2f}\n")
        f.write(f"95% CI = ({ci_i[0]:.2f}, {ci_i[1]:.2f})\n\n")
        f.write("[Correlated model]\n")
        f.write(f"tau = {tau_c:.2f}, rho = {rho:.2f}\n")
        f.write(f"Approx. CI = ({ci_c[0]:.2f}, {ci_c[1]:.2f})\n")
