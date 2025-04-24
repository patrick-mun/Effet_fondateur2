# reporting.py
"""
Génération du rapport final avec figures clés et export PDF et HTML
"""
import os
import pandas as pd
import matplotlib.pyplot as plt
from fpdf import FPDF

def plot_combined_ld(ld_csv_path: str, output_fig: str):
    df = pd.read_csv(ld_csv_path)
    if not all(x in df.columns for x in ["R2", "D"]):
        raise ValueError("Les colonnes R2 et D doivent être présentes")

    fig, ax1 = plt.subplots(figsize=(10, 5))
    ax1.set_xlabel("Index des paires de SNPs")
    ax1.set_ylabel("r²", color="tab:blue")
    ax1.plot(df["R2"].values, color="tab:blue", label="r²")
    ax1.tick_params(axis='y', labelcolor="tab:blue")

    ax2 = ax1.twinx()
    ax2.set_ylabel("D", color="tab:red")
    ax2.plot(df["D"].values, color="tab:red", linestyle="dashed", label="D")
    ax2.tick_params(axis='y', labelcolor="tab:red")

    fig.tight_layout()
    plt.title("Linkage Disequilibrium combiné : r² et D")
    plt.savefig(output_fig)
    plt.close()

def generate_pdf_report(figures: dict, output_pdf: str):
    pdf = FPDF()
    pdf.set_auto_page_break(auto=True, margin=15)
    for title, fig_path in figures.items():
        if not fig_path.lower().endswith(('.png', '.jpg', '.jpeg')):
            continue
        pdf.add_page()
        pdf.set_font("Arial", 'B', 14)
        pdf.cell(0, 10, title, ln=True, align='C')
        pdf.ln(10)
        if os.path.exists(fig_path):
            pdf.image(fig_path, w=180)
        else:
            pdf.set_font("Arial", '', 12)
            pdf.multi_cell(0, 10, f"Figure introuvable : {fig_path}")
    pdf.output(output_pdf)
    print(f"[Reporting] Rapport PDF sauvegardé : {output_pdf}")

def generate_html_report(figures: dict, output_html: str, summary_csv: str = None):
    summary_html = ""
    if summary_csv and os.path.exists(summary_csv):
        try:
            df = pd.read_csv(summary_csv)
            summary_html += "<section id='resume'><h2>Résumé des analyses</h2>"
            summary_html += "<table border='1' cellpadding='5' cellspacing='0' style='width:100%;border-collapse:collapse;'>"
            summary_html += "<tr><th>Module</th><th>Statistique</th><th>Commentaire</th></tr>"
            for _, row in df.iterrows():
                summary_html += f"<tr><td>{row['Module']}</td><td>{row['Statistique']}</td><td>{row['Commentaire']}</td></tr>"
            summary_html += "</table></section>"
        except Exception as e:
            summary_html += f"<p><strong>Erreur lecture summary.csv :</strong> {e}</p>"

    with open(output_html, "w") as f:
        f.write("""
        <html>
        <head>
            <meta charset='utf-8'>
            <title>Rapport d'analyse génétique</title>
            <style>
                body { font-family: Arial, sans-serif; margin: 0; padding: 0; background: #f4f4f4; }
                header { background: #005f73; color: white; padding: 1em; text-align: center; position: sticky; top: 0; }
                nav { background: #0a9396; padding: 0.5em; text-align: center; }
                nav a { color: white; margin: 0 1em; text-decoration: none; font-weight: bold; }
                nav a:hover { text-decoration: underline; }
                section { padding: 2em; background: white; margin: 2em auto; width: 80%; box-shadow: 0 0 10px rgba(0,0,0,0.1); }
                img { max-width: 100%; border: 1px solid #ccc; margin-top: 1em; }
                table { margin-top: 1em; width: 100%; border-collapse: collapse; background: white; }
                th, td { border: 1px solid #999; padding: 0.5em; text-align: center; }
                th { background: #e0f7fa; }
            </style>
        </head>
        <body>
        <header>
            <h1>Rapport de l'analyse génétique</h1>
        </header>
        <nav>
        <a href='#resume'>Résumé</a>
        """)
        for title in figures:
            anchor = title.replace(" ", "_").replace("/", "_").lower()
            f.write(f"<a href='#{anchor}'>{title}</a>")
        f.write("</nav>")

        f.write(summary_html)

        for title, fig_path in figures.items():
            anchor = title.replace(" ", "_").replace("/", "_").lower()
            f.write(f"<section id='{anchor}'>")
            f.write(f"<h2>{title}</h2>")
            if os.path.exists(fig_path):
                f.write(f"<img src='{fig_path}' alt='{title}'>")
            else:
                f.write(f"<p>Figure introuvable : {fig_path}</p>")
            f.write("</section>")

        f.write("</body></html>")
    print(f"[Reporting] Rapport HTML sauvegardé : {output_html}")

def generate_full_report(base_dir: str, output_pdf: str, output_html: str = None, summary_csv: str = None):
    figures = {}

    # ROH
    roh_plot = os.path.join(base_dir, "roh", "roh_plot.png")
    if os.path.exists(roh_plot):
        figures["Analyse ROH"] = roh_plot

    # LD
    ld_csv = os.path.join(base_dir, "ld", "ld_with_D.csv")
    ld_fig = os.path.join(base_dir, "ld", "ld_combined.png")
    if os.path.exists(ld_csv):
        plot_combined_ld(ld_csv, ld_fig)
        figures["Linkage Disequilibrium (r² et D)"] = ld_fig

    # Gamma
    gamma_fig = os.path.join(base_dir, "gamma", "gamma_plot.png")
    if os.path.exists(gamma_fig):
        figures["Analyse Gamma"] = gamma_fig

    # Adegenet
    adegenet_dapc = os.path.join(base_dir, "adegenet", "dapc_plot.png")
    adegenet_tree = os.path.join(base_dir, "adegenet", "dapc_tree.png")
    adegenet_hexp = os.path.join(base_dir, "adegenet", "dapc_hexp.png")
    if os.path.exists(adegenet_dapc):
        figures["Projection DAPC"] = adegenet_dapc
    if os.path.exists(adegenet_tree):
        figures["Arbre phylogénétique"] = adegenet_tree
    if os.path.exists(adegenet_hexp):
        figures["Hobs / Hexp"] = adegenet_hexp

    generate_pdf_report(figures, output_pdf)
    if output_html:
        generate_html_report(figures, output_html, summary_csv)
