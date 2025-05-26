# reporting.py
"""
Génération du rapport final avec figures clés et export PDF et HTML
pour le projet Effet fondateur
"""
import os
import pandas as pd
import matplotlib.pyplot as plt
from fpdf import FPDF
from datetime import datetime
import json
from scripts.ld import analyse_ld_in_roh
from scripts.ld import ld_mutation_zone_analysis


def ascii_safe(text: str) -> str:
    """
    Nettoie les caractères Unicode non supportés par fpdf (latin-1).
    Remplace les caractères typographiques et accents courants.
    """
    replacements = {
        "–": "-", "—": "-", "’": "'", "‘": "'", "“": '"', "”": '"',
        "…": "...", "•": "-", "é": "e", "è": "e", "ê": "e", "ë": "e",
        "à": "a", "â": "a", "ä": "a", "î": "i", "ï": "i",
        "ô": "o", "ö": "o", "ù": "u", "û": "u", "ü": "u", "ç": "c"
    }
    for k, v in replacements.items():
        text = text.replace(k, v)
    return text.encode("latin1", errors="ignore").decode("latin1")

# === EXTRACTION STRUCTURE FAMILLE À PARTIR DU FICHIER .PED ===

def extract_family_structure(ped_path):
    if not os.path.exists(ped_path):
        return pd.DataFrame()
    columns = ["FID", "IID", "PID", "MID", "SEX", "PHENO"]
    df = pd.read_csv(ped_path, sep=r'\s+', header=None, usecols=range(6), names=columns)
    df = df[~df['FID'].str.upper().str.contains('CTRL')]
    return df

# === CHARGER LES INFORMATIONS DU PROJET ===

def load_project_info(info_path):
    if os.path.exists(info_path):
        with open(info_path, "r") as f:
            return json.load(f)
    return {}

# === GÉNÉRATION DU TABLEAU RÉSUMÉ HWE ===

# ------------------------------------ Début Section Figure LD ------------------------------------
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

# ------------------------------------ Fin Section Figure LD ------------------------------------

def generate_hwe_summary_csv(csv_path):
    hwe_path = os.path.join(os.path.dirname(csv_path), "hwe.hwe")
    if not os.path.exists(hwe_path):
        raise FileNotFoundError(f"Fichier HWE introuvable : {hwe_path}")

    df = pd.read_csv(hwe_path, sep=r'\s+')
    total_snps = df.shape[0]
    p_lt_0_05 = (df["P"] < 0.05).sum()
    bonf_threshold = 0.05 / total_snps
    p_lt_bonf = (df["P"] < bonf_threshold).sum()

    summary_df = pd.DataFrame({
        "Module": ["PLINK - HWE"] * 3,
        "Statistique": [
            "SNPs totaux testés",
            "SNPs en déséquilibre (p < 0.05)",
            f"SNPs significatifs après correction Bonferroni (p < {bonf_threshold:.2e})"
        ],
        "Commentaire": [
            f"{total_snps} SNPs analysés pour l'équilibre Hardy-Weinberg.",
            f"{p_lt_0_05} SNPs avec p < 0.05 ({p_lt_0_05 / total_snps:.1%}).",
            f"{p_lt_bonf} SNPs significatifs après correction stricte."
        ]
    })
    os.makedirs(os.path.dirname(csv_path), exist_ok=True)
    summary_df.to_csv(csv_path, index=False)

# === AJOUT: charger le tableau HWE rsMUT ===
def load_hwe_rsMUT_html_table():
    hwe_table_path = os.path.join("data", "output", "complex_simulation", "geno", "hwe_rsMUT_summary.html")
    if os.path.exists(hwe_table_path):
        with open(hwe_table_path, "r") as f:
            return f.read()
    return "<p>Tableau HWE rsMUT non trouvé.</p>"

# === AJOUT: charger le R_sum__des_r_sultats_HWE.csv sous forme HTML ===
def load_hwe_summary_table(summary_csv):
    if os.path.exists(summary_csv):
        try:
            df = pd.read_csv(summary_csv)
            html = "<section><h3>Tableau Résumé des Résultats HWE</h3>"
            html += "<table border='1' cellpadding='5' cellspacing='0' style='width:100%;border-collapse:collapse;'>"
            html += "<tr><th>Module</th><th>Statistique</th><th>Commentaire</th></tr>"
            for _, row in df.iterrows():
                html += f"<tr><td>{row['Module']}</td><td>{row['Statistique']}</td><td>{row['Commentaire']}</td></tr>"
            html += "</table></section>"
            return html
        except Exception as e:
            return f"<p><strong>Erreur lecture summary.csv :</strong> {e}</p>"
    return "<p>Fichier R_sum__des_r_sultats_HWE.csv non trouvé.</p>"

import matplotlib.pyplot as plt
import pandas as pd
import os

def plot_ld_graphics(ld_csv_path: str, output_dir: str):
    """
    Génère deux graphiques à partir d'un fichier LD :
    1. Histogrammes comparatifs de r² et D
    2. Courbes de tendance (binned) de r² et D selon la distance
    """
    if not os.path.exists(ld_csv_path):
        print(f"[LD] Fichier introuvable : {ld_csv_path}")
        return

    df = pd.read_csv(ld_csv_path)

    # Vérification et génération de la colonne "Distance"
    if "Distance" not in df.columns:
        if "BP_A" in df.columns and "BP_B" in df.columns:
            df["Distance"] = abs(df["BP_B"] - df["BP_A"])
        else:
            print("[LD] Impossible de calculer la distance : colonnes BP_A ou BP_B absentes.")
            return

    df = df.dropna(subset=["R2", "D", "Distance"])

    # 1. Histogrammes comparatifs
    plt.figure(figsize=(10, 5))
    plt.hist(df["R2"], bins=100, alpha=0.6, label="r²", color="skyblue")
    plt.hist(df["D"], bins=100, alpha=0.6, label="D", color="tomato")
    plt.xlabel("Valeur du coefficient")
    plt.ylabel("Nombre de paires de SNPs")
    plt.title("Distribution comparée de r² et D")
    plt.legend()
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, "ld_histograms.png"))
    plt.close()

    # 2. Moyenne de r² et D par bin de distance
    df["bin"] = pd.cut(df["Distance"], bins=100)
    grouped = df.groupby("bin", observed=False)[["R2", "D"]].mean()
    bin_centers = [interval.mid for interval in grouped.index.categories]

    plt.figure(figsize=(10, 5))
    plt.plot(bin_centers, grouped["R2"], label="r²", color="blue")
    plt.plot(bin_centers, grouped["D"], label="D", color="red")
    plt.xlabel("Distance entre SNPs (binned)")
    plt.ylabel("Valeur moyenne")
    plt.title("Évolution moyenne de r² et D selon la distance")
    plt.legend()
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, "ld_binned_trend.png"))
    plt.close()


# --------------------------------------------debut mise en page PDF--------------------------------------------------------

# === CLASSE PDF AMÉLIORÉE ===

class CustomPDF(FPDF):
    def header(self):
        pass  

    def section_title(self, title):
        self.set_fill_color(0, 138, 160)
        self.set_text_color(255)
        self.set_font("Arial", "B", 13)
        self.cell(0, 12, title, ln=True, align='C', fill=True)
        self.ln(6)

    def section_paragraph(self, text):
        self.set_font("Arial", "", 11)
        self.set_text_color(50)
        self.multi_cell(0, 8, text)
        self.ln(3)

    def section_table(self, rows):
        self.set_font("Arial", "", 10)
        self.set_text_color(0)
        toggle = True
        for row in rows:
            self.set_fill_color(240, 248, 255) if toggle else self.set_fill_color(255, 255, 255)
            self.set_text_color(0)
            self.cell(0, 8, row, ln=True, fill=True)
            toggle = not toggle
        self.ln(3)

    def hwe_summary_table(self, headers, data):
        self.set_font("Arial", "B", 11)
        self.set_text_color(0)
        self.set_fill_color(224, 247, 250)
        self.set_draw_color(200)
        col_widths = [60, 60, 70]
        table_x = (210 - sum(col_widths)) / 2
        self.set_x(table_x)
        for i, h in enumerate(headers):
            self.cell(col_widths[i], 10, h, border=1, align='C', fill=True)
        self.ln()
        self.set_font("Arial", "", 10)
        toggle = True
        for row in data:
            self.set_fill_color(240, 248, 255) if toggle else self.set_fill_color(255, 255, 255)
            self.set_x(table_x)
            y_start = self.get_y()
            row_height = 20  # Hauteur de ligne doublée
            max_lines = 3
            for i, cell in enumerate(row):
                x = self.get_x()
                y = self.get_y()
                self.multi_cell(col_widths[i], row_height / max_lines, str(cell), border=1, align='C', fill=True)
                self.set_xy(x + col_widths[i], y)
            self.set_y(y_start + row_height)
            toggle = not toggle

    def add_hwe_summary(self, summary_csv):
        if summary_csv and os.path.exists(summary_csv):
            try:
                df = pd.read_csv(summary_csv)
                self.add_page()
                self.section_title("Résumé des analyses HWE")

                commentaire_path = os.path.join("data", "input", "complex_simulation", "resume_commentaire.json")
                commentaire_resume = ""
                if os.path.exists(commentaire_path):
                    with open(commentaire_path, "r", encoding="utf-8") as f:
                        resume_commentaire = json.load(f)
                        commentaire_resume = resume_commentaire.get("Résumé des analyses", "")

                if commentaire_resume:
                    self.section_paragraph(commentaire_resume.encode('latin1', 'replace').decode('latin1'))

                headers = ["Module", "Statistique", "Commentaire"]
                data = df.values.tolist()
                self.hwe_summary_table_html_style(headers, data)

            except Exception as e:
                self.add_page()
                self.set_font("Arial", "", 10)
                self.multi_cell(0, 10, f"Erreur de lecture summary_csv : {e}")

# === GÉNÉRATION DU RAPPORT PDF AMÉLIORÉ ===
def generate_pdf_report(figures: dict, output_pdf: str, summary_csv: str = None, ped_path: str = None):
    pdf = CustomPDF()
    pdf.set_auto_page_break(auto=True, margin=15)

    # Page de garde
    pdf.add_page()
    # Bandeau coloré en haut
    pdf.set_fill_color(0, 95, 115)  # couleur #005f73
    pdf.rect(10, 10, 190, 15, 'F')
    pdf.set_text_color(255)
    pdf.set_xy(10, 10)
    pdf.set_font("Arial", 'B', 14)
    pdf.cell(190, 15, "Projet : DOCK6", ln=True, align='C')

    # Titre centré
    pdf.set_y(40)
    pdf.set_text_color(0)
    pdf.set_font("Arial", 'B', 24)
    pdf.cell(0, 20, "Rapport - Analyse de l'effet fondateur", ln=True, align='C')

    # Date
    pdf.ln(10)
    pdf.set_font("Arial", '', 14)
    pdf.cell(0, 10, ascii_safe(f"Date : {datetime.today().strftime('%d/%m/%Y')}"), ln=True, align='C')
    pdf.ln(20)

    # Section Description médicale
    pdf.add_page()
    pdf.section_title("Description médicale")
    desc_path = os.path.join("data", "input", "complex_simulation", "info_med.json")
    info_med_content = ""
    if os.path.exists(desc_path):
        with open(desc_path, "r", encoding="utf-8") as f:
            try:
                info_med = json.load(f)
                info_med_content = info_med.get("Description", "")
            except:
                info_med_content = ""
    if info_med_content:
        pdf.section_paragraph(str(info_med_content).encode('latin1', 'replace').decode('latin1'))

    # Section Projet
    project_info = load_project_info(os.path.join("data", "input", "complex_simulation", "project_info.json"))
    if project_info:
        pdf.ln(8)
        pdf.section_title("Informations sur la mutation étudiée")
        pdf.set_font("Arial", "", 10)
        pdf.set_text_color(0)
        table_x = (210 - 160) / 2  # centrer le tableau (largeur page A4 - total colonnes)
        col_widths = [80, 80]
        toggle = True
        for k, v in project_info.items():
            pdf.set_x(table_x)
            pdf.set_fill_color(224, 247, 250)
            pdf.set_draw_color(200)
            pdf.set_text_color(0)
            pdf.set_font("Arial", "B", 10)
            pdf.cell(col_widths[0], 10, str(k).encode('latin1', 'replace').decode('latin1'), border=1, align='C', fill=True)
            pdf.set_font("Arial", "", 10)
            pdf.cell(col_widths[1], 10, str(v).encode('latin1', 'replace').decode('latin1'), border=1, align='C', fill=False)
            pdf.ln(10)
            toggle = not toggle
        pdf.ln(3)

    # Familles
    if ped_path and os.path.exists(ped_path):
        # Préparer les couleurs par famille
        df_tmp = extract_family_structure(ped_path)
        fid_order = list(df_tmp['FID'])
        fid_colors = {}
        previous_fid = None
        color_toggle = True
        for fid in fid_order:
            if fid != previous_fid:
                color_toggle = not color_toggle
            fid_colors[fid] = color_toggle
            previous_fid = fid
        pdf.add_page()
        pdf.section_title("Structure des familles")
        df = extract_family_structure(ped_path)
        columns = ["Famille", "Individu", "Père", "Mère", "Sexe", "Phénotype"]
        col_widths = [30, 40, 30, 30, 25, 25]
        pdf.set_font("Arial", "B", 10)
        pdf.set_fill_color(224, 247, 250)
        pdf.set_draw_color(200)
        table_x = (210 - sum(col_widths)) / 2
        pdf.set_x(table_x)
        for i, col in enumerate(columns):
            pdf.set_text_color(0)
            pdf.set_font("Arial", "B", 10)
            pdf.cell(col_widths[i], 8, col, border=1, align='C', fill=True)
        pdf.ln()
        toggle = True
        pdf.set_font("Arial", "", 10)
        current_fid = None
        fid_colors = {}
        current_color = True
        for _, r in df.iterrows():
            table_x = (210 - sum(col_widths)) / 2
            pdf.set_x(table_x)
            table_x = (210 - sum(col_widths)) / 2
            pdf.set_x(table_x)
            pdf.set_text_color(0)
            fill_color = (240, 248, 255) if fid_colors.get(r['FID'], True) else (255, 255, 255)
            pdf.set_fill_color(*fill_color)
            values = [r['FID'], r['IID'], r['PID'], r['MID'], str(r['SEX']), str(r['PHENO'])]
            for i, val in enumerate(values):
                pdf.cell(col_widths[i], 8, str(val), border=1, align='C', fill=True)
            pdf.ln()
            toggle = not toggle
        pdf.ln(3)
        pdf.set_font("Arial", "I", 9)
        pdf.set_text_color(80)
        pdf.multi_cell(0, 7, ascii_safe("Légende : Sexe = 1 (Homme), 2 (Femme) · Phénotype = 1 (non atteint), 2 (atteint), 3 (porteur)"))

        pdf.ln(3)

    # Filtres PLINK
    pdf.add_page()
    pdf.section_title("Critères de filtrage PLINK")
    pdf.set_font("Arial", "", 10)
    pdf.set_text_color(0)
    pdf.set_draw_color(200)
    col_widths = [95, 95]
    table_x = (210 - sum(col_widths)) / 2
    headers = ["Paramètre", "Description"]
    pdf.set_fill_color(224, 247, 250)
    pdf.set_x(table_x)
    for i, col in enumerate(headers):
        pdf.set_font("Arial", "B", 10)
        pdf.cell(col_widths[i], 10, col, border=1, align='C', fill=True)
    pdf.ln()
    filters = [
        ["--maf 0.01", "MAF < 1%"],
        ["--geno 0.05", "geno missing < 5%"],
        ["--mind 0.1", "mind missing < 10%"],
        ["--hwe 1e-6", "HWE p < 1e-6"]
    ]
    for row in filters:
        pdf.set_x(table_x)
        for i, val in enumerate(row):
            pdf.set_font("Arial", "", 10)
            pdf.cell(col_widths[i], 10, val, border=1, align='C')
        pdf.ln()
    pdf.ln(3)
#---------------------------------------début Section résumé des analyses------------------------------------
    def add_hwe_summary(pdf, summary_csv):
        if summary_csv and os.path.exists(summary_csv):
            try:
                df = pd.read_csv(summary_csv)
                pdf.add_page()
                pdf.set_font("Arial", 'B', 14)
                pdf.cell(0, 10, "Résumé des résultats HWE", ln=True, align='C')
                pdf.ln(5)

                commentaire_path = os.path.join("data", "input", "complex_simulation", "resume_commentaire.json")
                resume_commentaire = {}
                if os.path.exists(commentaire_path):
                    with open(commentaire_path, "r", encoding="utf-8") as f:
                        resume_commentaire = json.load(f)
                commentaire_resume = resume_commentaire.get("Résumé des analyses", "")
                if commentaire_resume:
                    pdf.set_font("Arial", "I", 10)
                    pdf.set_text_color(0)
                    pdf.multi_cell(0, 7, ascii_safe(str(commentaire_resume)))
                    pdf.ln(3)

                headers = ["Module", "Statistique", "Commentaire"]
                col_widths = [40, 70, 80]
                table_x = (210 - sum(col_widths)) / 2
                pdf.set_font("Arial", "B", 10)
                pdf.set_fill_color(224, 247, 250)
                pdf.set_draw_color(200)
                pdf.set_text_color(0)
                pdf.set_x(table_x)
                header_height = 14
                for i, col in enumerate(headers):
                    pdf.cell(col_widths[i], header_height, col, border=1, align='C', fill=True)
                pdf.ln(header_height)

                pdf.set_font("Arial", "", 10)
                toggle = True
                for _, r in df.iterrows():
                    y_start = pdf.get_y()
                    row_height = 20
                    values = [r['Module'], r['Statistique'], r['Commentaire']]
                    cell_heights = []
                    for i, val in enumerate(values):
                        text = str(val).encode('latin1', 'replace').decode('latin1')
                        nb_lines = pdf.get_string_width(text) / col_widths[i]
                        nb_lines = int(nb_lines) + 1
                        cell_height = 5 * nb_lines
                        cell_heights.append(cell_height)
                    effective_row_height = max(row_height, max(cell_heights))
                    fill_color = (240, 248, 255) if toggle else (255, 255, 255)
                    pdf.set_fill_color(*fill_color)
                    for i, val in enumerate(values):
                        x_current = table_x + sum(col_widths[:i])
                        pdf.set_xy(x_current, y_start)
                        text = str(val).encode('latin1', 'replace').decode('latin1')
                        pdf.multi_cell(col_widths[i], 5, ascii_safe(text), border=1, align='C', fill=True)
                    pdf.set_y(y_start + effective_row_height)
                    toggle = not toggle
                pdf.ln(3)

            except Exception as e:
                pdf.set_font("Arial", "", 10)
                pdf.multi_cell(0, 10, ascii_safe(f"Erreur de lecture summary_csv : {e}"))



#---------------------------------------Fin section résumé des analyses-----------------------------------------
    # Figures
    for title, path in figures.items():
        if os.path.exists(path):
            pdf.add_page()
            pdf.section_title(title)
            pdf.image(path, w=180)

    pdf.output(output_pdf)
    print(f"[Reporting] Rapport PDF sauvegardé : {output_pdf}")


#----------------------------------------------fin mise en page pdf----------------------------------------------------------------------

#----------------------------------------------debut mise en page html----------------------------------------------------------------------

def generate_html_report(figures: dict, output_html: str, summary_csv: str = None, ped_path: str = None):
    """
    Génère un rapport HTML complet avec :
    - En-tête et navigation
    - Sections : Description médicale, Mutation, Familles, Critères PLINK, Résumé analyses, IBD
    - Affichage des figures
    """
    sections = []
    
    # 1. Description médicale
    info_med = {}
    info_med_path = os.path.join("data", "input", "complex_simulation", "info_med.json")
    if os.path.exists(info_med_path):
        with open(info_med_path) as f:
            info_med = json.load(f)
    desc = info_med.get("Description")
    if desc:
        sections.append(
            f"<section id='description'><h2>Description médicale</h2><p>{desc}</p></section>"
        )

    # 2. Informations sur la mutation
    project_info = load_project_info(os.path.join("data", "input", "complex_simulation", "project_info.json"))
    if project_info:
        rows = ''.join(f"<tr><th>{{k}}</th><td>{{v}}</td></tr>".format(k=k, v=v) for k, v in project_info.items())
        sections.append(
            "<section id='mutation'><h2>Informations sur la mutation étudiée</h2>"
            "<table border='1' style='width:100%;'>{rows}</table></section>".replace("{rows}", rows)
        )

    # 3. Structure des familles
    if ped_path and os.path.exists(ped_path):
        df_ped = extract_family_structure(ped_path)
        rows = []
        current_fid = None
        toggle = False
        for _, r in df_ped.iterrows():
            if r['FID'] != current_fid:
                toggle = not toggle
                current_fid = r['FID']
            bg = '#eef' if toggle else '#fff'
            rows.append(
                f"<tr style='background:{bg};'><td>{r.FID}</td><td>{r.IID}</td><td>{r.PID}</td>"
                f"<td>{r.MID}</td><td>{r.SEX}</td><td>{r.PHENO}</td></tr>"
            )
        sections.append(
            "<section id='familles'><h2>Structure des familles</h2>"
            "<table border='1' style='width:100%;'><tr>"
            "<th>Famille</th><th>IID</th><th>PID</th><th>MID</th><th>Sexe</th><th>Phéno</th></tr>"
            + ''.join(rows) +
            "</table></section>"
        )

    # 4. Critères PLINK 
    if summary_csv and os.path.exists(summary_csv):
        try:
            df_sum = pd.read_csv(summary_csv)
            # Critères
            criteria = [
                ("--maf 0.01", "MAF < 1%"),
                ("--geno 0.05", "geno missing < 5%"),
                ("--mind 0.1", "mind missing < 10%"),
                ("--hwe 1e-6", "HWE p < 1e-6")
            ]
            crit_rows = ''.join(f"<tr><td>{c[0]}</td><td>{c[1]}</td></tr>" for c in criteria)
            sections.append(
                f"<section id='critere'><h2>Critères de filtrage PLINK</h2>"
                f"<table border='1' style='width:100%;'>{crit_rows}</table></section>"
            )
            
        except Exception as e:
            sections.append(f"<p>Erreur lecture summary_csv: {e}</p>")


    # Résumé des analyses HWE 
    if summary_csv and os.path.exists(summary_csv):
        try:
            df = pd.read_csv(summary_csv)
            
            # Lire le commentaire général
            commentaire_path = os.path.join("data", "input", "complex_simulation", "resume_commentaire.json")
            resume_commentaire = {}
            if os.path.exists(commentaire_path):
                with open(commentaire_path, "r", encoding="utf-8") as f:
                    resume_commentaire = json.load(f)
            commentaire_resume = resume_commentaire.get("Résumé des analyses", "")

            # Construire les lignes du tableau
            table_rows = ''.join(
                f"<tr style='background-color:#{'f0f8ff' if i%2==0 else 'ffffff'};'>"
                f"<td>{r.Module}</td><td>{r.Statistique}</td><td>{r.Commentaire}</td></tr>"
                for i, r in enumerate(df.itertuples())
            )

            sections.append(
                "<section id='resume'><h2>Résumé des analyses HWE</h2>" +
                (f"<p><em>{commentaire_resume}</em></p>" if commentaire_resume else "") +
                "<table border='1' style='width:100%;border-collapse:collapse;margin-top:1em;'>"
                "<tr style='background:#e0f7fa;'><th>Module</th><th>Statistique</th><th>Commentaire</th></tr>"
                f"{table_rows}</table></section>"
            )
        except Exception as e:
            sections.append(f"<p><strong>Erreur lecture summary_csv :</strong> {e}</p>")



    # Résumé IBD
    try:
        plink_file = os.path.join("data", "output", "complex_simulation", "ibd", "ibd_plink.genome")
        king_file = os.path.join("data", "output", "complex_simulation", "ibd", "ibd_king.kin0")
        network_img_rel = "ibd/ibd_king_network.png"
        network_img_abs = os.path.join(os.path.dirname(output_html), network_img_rel)

        if os.path.exists(plink_file) and os.path.exists(king_file):
            df_pl = pd.read_csv(plink_file, sep=r'\s+')
            df_kg = pd.read_csv(king_file, sep=r'\s+')

            # PLINK
            if 'RT' in df_pl.columns:
                ps = df_pl['RT'].value_counts().rename_axis('Relation').reset_index(name='PLINK')
            else:
                print("[DEBUG] Colonne 'RT' absente dans ibd_plink.genome")
                ps = pd.DataFrame(columns=['Relation', 'PLINK'])
        
            # KING
            if 'Kinship' in df_kg.columns:
                def cls(x):
                    if x > 0.354: return 'Jumeaux'
                    if x > 0.177: return '1er'
                    if x > 0.0884: return '2e'
                    if x > 0.0442: return '3e'
                    if x > 0.0221: return '4e'
                    return 'Non'
                df_kg['Relation'] = df_kg['Kinship'].apply(cls)
                ks = df_kg['Relation'].value_counts().rename_axis('Relation').reset_index(name='KING')
            else:
                print("[DEBUG] Colonne 'Kinship' absente dans ibd_king.kin0")
                ks = pd.DataFrame(columns=['Relation', 'KING'])

            ibd = ps.merge(ks, on='Relation', how='outer').fillna(0)

            ibd_rows = ''.join(
                f"<tr><td>{r.Relation}</td><td>{int(r.PLINK)}</td><td>{int(r.KING)}</td></tr>"
                for r in ibd.itertuples()
            )

            # Création de la section HTML
            section_ibd_html = (
                "<section id='ibd'><h2>Résumé des analyses IBD</h2>"
                "<table border='1' style='width:100%; border-collapse: collapse;'>"
                "<tr><th>Relation</th><th>PLINK</th><th>KING</th></tr>"
                f"{ibd_rows}</table>"
                "<p style='margin-top: 1em;'><strong>Légende :</strong></p>"
                "<ul>"
                "<li><strong>PO</strong> : Parent–Enfant</li>"
                "<li><strong>FS</strong> : Frères/Sœurs</li>"
                "<li><strong>2e</strong> : Degré 2 (grand-parent, oncle/tante)</li>"
                "<li><strong>3e</strong> : Degré 3 (cousin germain)</li>"
                "<li><strong>4e</strong> : Degré 4 (arrière-cousin)</li>"
                "<li><strong>OT</strong> : Autres / hors seuil</li>"
                "<li><strong>UN</strong> : Non apparentés</li>"
                "<li><strong>Non</strong> : Non classé (en dessous du seuil de KING)</li>"
                "<li><strong>Jumeaux</strong> : Relation gémellaire (KING uniquement)</li>"
                "</ul>"
            )

            # Ajout de l'image si elle existe (à l'intérieur du bloc !)
            if os.path.exists(network_img_abs):
                section_ibd_html += f"<img src='{network_img_rel}' alt='Graphe IBD KING' style='margin-top:1em; max-width:100%;'>"

            section_ibd_html += "</section>"
            sections.append(section_ibd_html)

            # Section d'interprétation du graphe IBD (lecture JSON dédiée)
            ibd_comment_path = os.path.join("data", "input", "complex_simulation", "resume_commentaire_idb.json")
            if os.path.exists(ibd_comment_path):
                print("[DEBUG] Chargement des commentaires IBD depuis JSON")
                try:
                    with open(ibd_comment_path, "r", encoding="utf-8") as f:
                        commentaire_ibd = json.load(f)

                    ibd_comment_html = "<section id='ibd_commentaire'><h2>Interprétation du graphe IBD KING</h2>"
                    for titre, texte in commentaire_ibd.items():
                        ibd_comment_html += f"<h3>{titre}</h3><p>{texte}</p>"
                    ibd_comment_html += "</section>"

                    sections.append(ibd_comment_html)
                except Exception as e:
                    sections.append(f"<p><strong>Erreur lecture commentaire IBD :</strong> {e}</p>")
            else:
                print(f"[DEBUG] Fichier JSON IBD non trouvé : {ibd_comment_path}")


    except Exception as e:
        sections.append(f"<p><strong>Erreur IBD :</strong> {e}</p>")


   # Résumé visuel LD (r² et D)
    ld_hist_rel = "ld/ld_histograms.png"
    ld_trend_rel = "ld/ld_binned_trend.png"
    ld_roh_rel = "ld/ld_r2_comparaison_roh.png"

    ld_hist_abs = os.path.join(os.path.dirname(output_html), ld_hist_rel)
    ld_trend_abs = os.path.join(os.path.dirname(output_html), ld_trend_rel)
    ld_roh_abs = os.path.join(os.path.dirname(output_html), ld_roh_rel)

    if os.path.exists(ld_hist_abs) and os.path.exists(ld_trend_abs):
        # Commentaire LD1 (histogrammes)
        ld_comment_html = ""
        ld1_path = os.path.join("data", "input", "complex_simulation", "resume_commentaire_ld1.json")
        if os.path.exists(ld1_path):
            try:
                with open(ld1_path, "r", encoding="utf-8") as f:
                    ld1_data = json.load(f)
                    if "Analyse LD" in ld1_data:
                        ld_comment_html = "<div style='margin-top:1em;'><h3>Interprétation des résultats LD</h3><ul>"
                        for k, v in ld1_data["Analyse LD"].items():
                            ld_comment_html += f"<li><strong>{k} :</strong> {v}</li>"
                        ld_comment_html += "</ul></div>"
            except Exception as e:
                print(f"[DEBUG] Erreur JSON LD1 : {e}")

        # Commentaire LD2 (binned trend)
        ld_comment_html2 = ""
        ld2_path = os.path.join("data", "input", "complex_simulation", "resume_commentaire_ld2.json")
        if os.path.exists(ld2_path):
            try:
                with open(ld2_path, "r", encoding="utf-8") as f:
                    ld2_data = json.load(f)
                    if "Analyse LD binned" in ld2_data:
                        ld_comment_html2 = "<div style='margin-top:1em;'><h3>Interprétation du profil LD par distance</h3><ul>"
                        for k, v in ld2_data["Analyse LD binned"].items():
                            ld_comment_html2 += f"<li><strong>{k} :</strong> {v}</li>"
                        ld_comment_html2 += "</ul></div>"
            except Exception as e:
                print(f"[DEBUG] Erreur JSON LD2 : {e}")

        # Commentaire LD vs ROH (ld_r2_comparaison_roh.png)
        ld_roh_comment_html = ""
        roh_ld_json = os.path.join("data", "input", "complex_simulation", "resume_commentaire_ld_roh.json")
        if os.path.exists(roh_ld_json):
            try:
                with open(roh_ld_json, "r", encoding="utf-8") as f:
                    roh_ld_data = json.load(f)
                    if "Analyse LD vs ROH" in roh_ld_data:
                        ld_roh_comment_html = "<div style='margin-top:1em;'><h3>Interprétation LD vs ROH</h3><ul>"
                        for k, v in roh_ld_data["Analyse LD vs ROH"].items():
                            ld_roh_comment_html += f"<li><strong>{k} :</strong> {v}</li>"
                        ld_roh_comment_html += "</ul></div>"
            except Exception as e:
                print(f"[DEBUG] Erreur JSON LD ROH : {e}")

        # Construction de la section HTML
        html_ld_section = (
            "<section id='ld_summary'><h2>Analyse du déséquilibre de liaison (LD)</h2>"
            "<p>Trois visualisations complémentaires comparent les coefficients <strong>r²</strong> et <strong>D</strong> :</p>"
            "<ul>"
            "<li><strong>Histogramme :</strong> distribution globale de r² et D</li>"
            "<li><strong>Courbes binned :</strong> évolution moyenne de r² et D selon la distance entre SNPs</li>"
            "<li><strong>Comparaison dans/hors ROH :</strong> distribution de r² selon les zones d'homozygotie</li>"
            "</ul>"
            f"<h3>Distribution des coefficients r² et D</h3><img src='{ld_hist_rel}' alt='Histogrammes LD'>"
            f"{ld_comment_html}"
            f"<h3>Tendance moyenne par distance</h3><img src='{ld_trend_rel}' alt='Tendance LD binned'>"
            f"{ld_comment_html2}"
        )

        if os.path.exists(ld_roh_abs):
            html_ld_section += (
                f"<h3>Distribution de r² dans et hors ROH</h3><img src='{ld_roh_rel}' alt='Comparaison LD ROH'>"
                f"{ld_roh_comment_html}"
            )

        html_ld_section += "</section>"
        sections.append(html_ld_section)

    else:
        print(f"[DEBUG] Graphiques LD manquants : {ld_hist_abs}, {ld_trend_abs}")

    # Section : Comparaison r² global vs r² dans la zone de la mutation chez les cas
    r2_global_path = os.path.join("data", "output", "complex_simulation", "ld", "ld_with_D.csv")
    r2_cas_path = os.path.join("data", "output", "complex_simulation", "ld", "ld_cas_mutation_zone_with_D.csv")
    compare_hist = os.path.join("data", "output", "complex_simulation", "ld", "compare_r2_mutation_zone.png")
    compare_box = os.path.join("data", "output", "complex_simulation", "ld", "boxplot_r2_mutation_zone.png")

    if os.path.exists(r2_global_path) and os.path.exists(r2_cas_path):
        try:
            df_global = pd.read_csv(r2_global_path)
            df_cas = pd.read_csv(r2_cas_path)

            r2_global = df_global["R2"].dropna()
            r2_cas = df_cas["R2"].dropna()

            # Histogramme
            import matplotlib.pyplot as plt
            plt.figure(figsize=(10, 5))
            plt.hist(r2_global, bins=50, alpha=0.6, label="r² global", color="skyblue", density=True)
            plt.hist(r2_cas, bins=50, alpha=0.6, label="r² atteints (±1Mb mutation)", color="salmon", density=True)
            plt.xlabel("Valeur de r²")
            plt.ylabel("Densité")
            plt.title("Comparaison de la distribution du r²")
            plt.legend()
            plt.tight_layout()
            plt.savefig(compare_hist)
            plt.close()

            # Boxplot
            plt.figure(figsize=(8, 5))
            plt.boxplot([r2_global, r2_cas], labels=["r² global", "r² atteints (mutation zone)"])
            plt.ylabel("Valeur de r²")
            plt.title("Comparaison des distributions de r²")
            plt.grid(axis="y")
            plt.tight_layout()
            plt.savefig(compare_box)
            plt.close()
            
           # Section HTML
            html_r2_section = (
                "<section id='compare_r2'><h2>Analyse comparative du r² global et local chez les atteints</h2>"
                "<p>Cette section compare le déséquilibre de liaison (r²) dans l'ensemble du génome avec celui observé dans la région ±1Mb autour de la mutation chez les individus atteints.</p>"
                f"<h3>Histogramme comparatif</h3><img src='{os.path.relpath(compare_hist, os.path.dirname(output_html))}' alt='Histogramme r² global vs atteints'>"
            )

            # 🔽 Commentaire sous histogramme
            hist_comment_path = os.path.join("data", "input", "complex_simulation", "resume_commentaire_ld_histogramme.json")
            if os.path.exists(hist_comment_path):
                try:
                    with open(hist_comment_path, "r", encoding="utf-8") as f:
                        hist_data = json.load(f)
                        if "Analyse LD histogramme" in hist_data:
                            html_r2_section += "<div style='margin-top:1em;'><h3>Interprétation de l'histogramme</h3><ul>"
                            for k, v in hist_data["Analyse LD histogramme"].items():
                                html_r2_section += f"<li><strong>{k} :</strong> {v}</li>"
                            html_r2_section += "</ul></div>"
                except Exception as e:
                    print(f"[DEBUG] Erreur lecture JSON commentaire histogramme : {e}")

            # 🔽 Ajout du boxplot
            html_r2_section += (
                f"<h3>Boxplot comparatif</h3><img src='{os.path.relpath(compare_box, os.path.dirname(output_html))}' alt='Boxplot r² global vs mutation zone'>"
            )

            # 🔽 Commentaire sous boxplot
            boxplot_comment_path = os.path.join("data", "input", "complex_simulation", "resume_commentaire_ld_boxplot.json")
            if os.path.exists(boxplot_comment_path):
                try:
                    with open(boxplot_comment_path, "r", encoding="utf-8") as f:
                        boxplot_data = json.load(f)
                        if "Analyse LD boxplot" in boxplot_data:
                            html_r2_section += "<div style='margin-top:1em;'><h3>Interprétation du boxplot</h3><ul>"
                            for k, v in boxplot_data["Analyse LD boxplot"].items():
                                html_r2_section += f"<li><strong>{k} :</strong> {v}</li>"
                            html_r2_section += "</ul></div>"
                except Exception as e:
                    print(f"[DEBUG] Erreur lecture JSON commentaire boxplot : {e}")

            html_r2_section += "</section>"
            sections.append(html_r2_section)
            
        except Exception as e:
            print(f"[DEBUG] Erreur comparaison r² global vs atteints : {e}")
            
        # 🔽 Conclusion et synthèse LD
        ld_synthese_path = os.path.join("data", "input", "complex_simulation", "resume_commentaire_ld_synthese.json")
        if os.path.exists(ld_synthese_path):
            try:
                with open(ld_synthese_path, "r", encoding="utf-8") as f:
                    ld_synthese_data = json.load(f)
                    synthese = ld_synthese_data.get("Conclusion et synthèse de la LD", {})
                    
                    html_synthese_section = "<section id='conclusion_ld'><h2>Conclusion et synthèse de la LD</h2>"
                    for titre, texte in synthese.items():
                        html_synthese_section += f"<h3>{titre}</h3><p>{texte}</p>"
                    html_synthese_section += "</section>"

                    sections.append(html_synthese_section)
            except Exception as e:
                print(f"[DEBUG] Erreur lecture JSON synthèse LD : {e}")

    # 🔽 Ajout graphique de superposition ROH chez les atteints
    roh_overlap_path = os.path.join("data", "output", "complex_simulation", "roh","figures","roh_overlap.png")
    if os.path.exists(roh_overlap_path):
        roh_rel_path = os.path.relpath(roh_overlap_path, os.path.dirname(output_html))
        roh_overlap_html = (
            "<section id='roh_overlap'>"
            "<h2>Superposition des segments ROH chez les atteints</h2>"
            "<p>Ce graphique visualise les segments d’homozygotie (ROH) individuels chez les personnes atteintes, "
            "ainsi que la région commune partagée par l’ensemble d’entre eux (zone surlignée en rose). "
            "La cohérence et le chevauchement de ces segments soutiennent l’hypothèse d’un héritage commun d’origine fondatrice.</p>"
            f"<img src='{roh_rel_path}' alt='Segments ROH atteints' style='max-width:100%; margin-top:1em;'>"
            "</section>"
        )
        sections.append(roh_overlap_html)
    else:
        print("[DEBUG] roh_overlap.png non trouvé – section ROH overlap non ajoutée.")


    # 7. Résumé ROH depuis roh.hom
    roh_file = os.path.join("data", "output", "complex_simulation", "roh", "roh.hom")
    if os.path.exists(roh_file):
        try:
            df_roh = pd.read_csv(roh_file, sep=r'\s+')
            display_cols = ["FID", "IID", "CHR", "KB", "NSNP"]
            if all(col in df_roh.columns for col in display_cols):
                rows = ''.join(
                    f"<tr><td>{r.FID}</td><td>{r.IID}</td><td>{r.CHR}</td><td>{r.KB}</td><td>{r.NSNP}</td></tr>"
                    for r in df_roh.itertuples()
                )
                sections.append(
                    "<section id='roh_table'><h2>Résumé des segments ROH</h2>"
                    "<table border='1' style='width:100%;'>"
                    "<tr><th>FID</th><th>IID</th><th>CHR</th><th>Longueur (Kb)</th><th>Nb SNPs</th></tr>"
                    + rows +
                    "</table></section>"
                )
        except Exception as e:
            sections.append(f"<p><strong>Erreur lecture roh.hom :</strong> {e}</p>")

        # 🔽 Section GAMMA – Estimation de l'âge de la mutation
        gamma_summary_path = os.path.join("data", "output", "complex_simulation", "gamma", "gamma_summary.txt")
        if os.path.exists(gamma_summary_path):
            try:
                with open(gamma_summary_path, "r", encoding="utf-8") as f:
                    lines = f.readlines()

                section_gamma = "<section id='gamma'><h2>Analyse Gamma – Estimation de l’âge de la mutation</h2>"
                section_gamma += (
                    "<p>L’analyse Gamma permet d’estimer l’âge d’une mutation en générations, en se basant sur la longueur des segments d’homozygotie (ROH) flanquant la mutation chez les individus atteints.</p>"
                    "<ul>"
                )

                for line in lines:
                    if line.strip():
                        section_gamma += f"<li>{line.strip()}</li>"
                section_gamma += "</ul>"

                # Ajout commentaire d'interprétation
                section_gamma += (
                    "<p>Ces résultats suggèrent que la mutation pourrait être apparue dans les <strong>3 à 5 dernières générations</strong>, selon le modèle retenu. "
                    "La forte corrélation intra-segmentaire (ρ ≈ 0.95) justifie l'utilisation du modèle corrélé dans un contexte d'homogénéité génétique.</p>"
                    "<p>Cette estimation appuie l’hypothèse d’un effet fondateur récent dans la population étudiée.</p>"
                    "</section>"
                )

                sections.append(section_gamma)

            except Exception as e:
                print(f"[DEBUG] Erreur lecture Gamma summary : {e}")


        # 🔽 Intégration graphique DAPC (Adegenet)
        dapc_img = os.path.join("data", "output", "complex_simulation", "adegenet", "dapc_plot.png")
        if os.path.exists(dapc_img):
            dapc_rel = os.path.relpath(dapc_img, os.path.dirname(output_html))
            html_dapc_section = (
                "<section id='dapc'><h2>Analyse DAPC – Structure génétique des individus</h2>"
                "<p>La DAPC (Discriminant Analysis of Principal Components) permet de visualiser la structuration des individus "
                "selon leurs profils génotypiques. Elle met en évidence les regroupements par familles, porteurs, atteints ou témoins, "
                "et aide à détecter une différenciation potentielle au sein de la population analysée.</p>"
                f"<img src='{dapc_rel}' alt='Projection DAPC' style='max-width:100%; margin-top:1em;'>"
            )

            # 🔽 Ajout commentaire DAPC depuis JSON
            dapc_comment_path = os.path.join("data", "input", "complex_simulation", "resume_commentaire_dapc.json")
            if os.path.exists(dapc_comment_path):
                try:
                    with open(dapc_comment_path, "r", encoding="utf-8") as f:
                        dapc_data = json.load(f)
                        dapc_html = "<div style='margin-top:1em;'><h3>Interprétation de la DAPC</h3><ul>"
                        for titre, texte in dapc_data.get("Interprétation de la DAPC", {}).items():
                            dapc_html += f"<li><strong>{titre} :</strong> {texte}</li>"
                        dapc_html += "</ul></div>"
                        html_dapc_section += dapc_html
                except Exception as e:
                    print(f"[DEBUG] Erreur lecture JSON DAPC : {e}")

            html_dapc_section += "</section>"
            sections.append(html_dapc_section)
        else:
            print("[DEBUG] Image DAPC introuvable – section DAPC non ajoutée.")

    # 🔽 Conclusion générale intégrée
    conclusion_path = os.path.join("data", "input", "complex_simulation", "resume_commentaire_conclusion_generale.json")
    if os.path.exists(conclusion_path):
        try:
            with open(conclusion_path, "r", encoding="utf-8") as f:
                conclusion_data = json.load(f)
                general_conclusion = conclusion_data.get("Conclusion générale", {})

                html_conclusion = "<section id='conclusion_generale'><h2>Conclusion générale</h2>"
                for sous_titre, texte in general_conclusion.items():
                    html_conclusion += f"<h3>{sous_titre}</h3><p>{texte}</p>"
                html_conclusion += "</section>"

                sections.append(html_conclusion)
        except Exception as e:
            print(f"[DEBUG] Erreur lecture JSON conclusion générale : {e}")

        # 🔽 Signature de l'auteur
        html_signature = (
            "<section id='signature' style='margin-top:2em; text-align:right;'>"
            "<h2 style='text-align:left;'>Auteur</h2>"
            "<p><strong>Patrick MUNIER</strong><br>"
            "Laboratoire de génétique – Projet DOCK6<br>"
            "Version : 1.0 ALPHA<br>"
            "Date : 03/04/2025</p>"
            "</section>"
        )
        sections.append(html_signature)

        #----------------------------------------------fin mise en page html----------------------------------------------------------------------

        html = [
        "<!DOCTYPE html>",
        "<html lang='fr'>",
        "<head><meta charset='utf-8'><title>Rapport Effet fondateur</title>"
        "<style>"
        "body { font-family: Arial, sans-serif; background:#f8f8f8; color:#222; margin:0; padding:0; }"
        "header { background:#005f73; color:white; padding:1em; text-align:center; }"
        "nav { background:#0a9396; padding:0.5em; text-align:center; }"
        "nav a { color:white; margin:0 1em; text-decoration:none; font-weight:bold; }"
        "nav a:hover { text-decoration:underline; }"
        "section { padding:2em; margin:2em auto; background:white; width:80%; box-shadow:0 0 5px rgba(0,0,0,0.1); border-radius: 8px; }"
        "h2 { font-size:1.8em; margin-bottom:0.4em; color:#003049; border-bottom:2px solid #ccc; padding-bottom:0.2em; }"
        "h3 { font-size:1.3em; margin-top:1.5em; color:#0a9396; }"
        "p, li { font-size:1.05em; line-height:1.6; text-align: justify; }"
        "table { width:100%; border-collapse:collapse; margin-top:1em; }"
        "th, td { border:1px solid #ccc; padding:0.5em; text-align:center; }"
        "th { background:#e0f7fa; }"
        "img { max-width:100%; margin-top:1em; border:1px solid #ddd; border-radius: 6px; }"
        "</style></head>",
        "<body>",
        # 🔽 Avertissement juste après <body> et avant <header>
        "<section id='avertissement' style='background:#fff3cd; color:#856404; border:1px solid #ffeeba; "
        "padding:1.5em; margin:2em auto; width:80%; border-radius:8px;'>"
        "<h2 style='color:#856404;'>Avertissement – Données simulées</h2>"
        "<p>Cette analyse de l’effet fondateur a été réalisée à partir de données générées informatiquement. "
        "Elle a pour unique objectif de démontrer la faisabilité et la validité du traitement bioinformatique mis en œuvre. "
        "Les résultats présentés dans ce rapport sont issus d’un pipeline expérimental testé en conditions simulées.<br><br>"
        "<strong>Aucune donnée ni résultat contenu dans ce rapport ne doit être utilisé à des fins de publication scientifique.</strong>"
        "</p></section>",
        "<header><h1>Rapport Effet fondateur</h1><nav>" +
        ''.join(f"<a href='#{k}'>{k}</a>" for k in [
            'description','mutation','familles','critere','resume','ibd'
        ]) +
        "</nav></header>",
    ]
    html += sections




    # 7. Ajout des figures
    for title, path in figures.items():
        anchor = title.replace(' ', '_')
        html.append(f"<section id='{anchor}'><h2>{title}</h2><img src='{path}' alt='{title}'></section>")

    html.append("</body></html>")

    with open(output_html, 'w') as f:
        f.write("\n".join(html))

    if not sections:
        print("[ERREUR] Aucune section générée !")
    else:
        print(f"[DEBUG] {len(sections)} sections prêtes à être écrites.")

    with open(output_html, 'w') as f:
        f.write("\n".join(html))
    print(f"[Reporting] Rapport HTML sauvegardé dans : {output_html}")


# === RAPPORT COMPLET ===
def generate_full_report(base_dir: str, output_pdf: str, output_html: str = None, summary_csv: str = None):
    """
    Génère un rapport PDF et HTML complet contenant :
    - Figures de qualité (MAF, HWE)
    - Résultats ROH, LD, Gamma, DAPC
    - Résumé des analyses depuis un fichier CSV
    - Intégration HTML de la structure familiale, mutation et IBD si disponible
    """
    figures = {}

    # Fichier CSV de résumé des analyses
    summary_csv_default = os.path.join(base_dir, "geno", "R_sum__des_r_sultats_HWE.csv")
    if summary_csv is None:
        summary_csv = summary_csv_default
    if not os.path.exists(summary_csv):
       

    # Graphiques de qualité (MAF, HWE)
        maf_plot = os.path.join(base_dir, "geno", "qc_maf_distribution.png")
        hwe_plot = os.path.join(base_dir, "geno", "qc_hwe_pvalues.png")
        if os.path.exists(maf_plot):
            figures["Distribution des MAF"] = maf_plot
        if os.path.exists(hwe_plot):
            figures["Distribution HWE (p-values)"] = hwe_plot

    # Génération des graphiques LD
    ld_csv = os.path.join(base_dir, "ld", "ld_with_D.csv")
    ld_dir = os.path.join(base_dir, "ld")
    plot_ld_graphics(ld_csv, ld_dir)

    # Génération des graphiques LD
    ld_csv = os.path.join(base_dir, "ld", "ld_with_D.csv")
    ld_fig = os.path.join(base_dir, "ld", "ld_combined.png")
    if os.path.exists(ld_csv):
        plot_combined_ld(ld_csv, ld_fig)
        analyse_ld_in_roh(
        ld_csv_path=os.path.join(base_dir, "ld", "ld_with_D.csv"),
        roh_file_path=os.path.join(base_dir, "roh", "roh.hom"),
        map_file_path=os.path.join("data", "input", "complex_simulation", "genotype_data.map"),
        fam_file_path=os.path.join(base_dir, "geno", "filtered_data.fam"),  # ← chemin ajouté
        output_dir=os.path.join(base_dir, "ld")
    )
    
       
    
    # Analyse LD dans la région ±1Mb autour de la mutation rsMUT (ex: position 11087401 sur chr19)
    ld_mutation_zone_analysis(
        plink_prefix=os.path.join(base_dir, "geno", "filtered_data"),
        fam_path=os.path.join(base_dir, "geno", "filtered_data.fam"),
        map_path=os.path.join(base_dir, "geno", "filtered_data.map"),
        output_dir=os.path.join(base_dir, "ld")
    )

    # ROH
    roh_plot = os.path.join(base_dir, "roh", "roh_plot.png")
    if os.path.exists(roh_plot):
        figures["Analyse ROH"] = roh_plot
    roh_overlap = os.path.join(base_dir, "roh", "roh_overlap.png")
    if os.path.exists(roh_overlap):
        figures["Segments ROH – Zone commune"] = roh_overlap


    # Génération des rapports (HTML enrichi + PDF simple)
    ped_path = os.path.join(base_dir, "geno", "filtered_data.ped")
    generate_pdf_report(figures, output_pdf)
    if output_html:
        generate_html_report(figures, output_html, summary_csv, ped_path=ped_path)