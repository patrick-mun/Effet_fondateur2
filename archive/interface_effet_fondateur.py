import streamlit as st
import subprocess
import sys
import pandas as pd
import streamlit.components.v1 as components
from pathlib import Path


# Répertoires par défaut
DEFAULT_PED = "data/input/complex_simulation/genotype_data.ped"
DEFAULT_MAP = "data/input/complex_simulation/genotype_data.map"
REPORT_PATH = "data/output/complex_simulation/rapport_final.html"
SUMMARY_CSV = "data/output/complex_simulation/summary.csv"
WIKI_DIR = "Module_WIKI"
README_FILE = "README.md"

st.set_page_config(page_title="Effet Fondateur - Pipeline", layout="wide")

st.sidebar.title("Navigation")
choice = st.sidebar.radio("Menu", ["🔁 Chargement & Analyse", "📊 Résultats & Rapport", "📚 Wiki"])

def launch_pipeline():
    with st.spinner("Exécution du pipeline en cours..."):
        try:
            subprocess.run([sys.executable, "run_pipeline.py"], check=True)
            st.success("Analyse terminée avec succès.")
        except subprocess.CalledProcessError as e:
            st.error(f"Erreur lors de l'exécution : {e}")

if choice == "🔁 Chargement & Analyse":
    st.header("🔁 Chargement des fichiers et lancement du pipeline")

    ped_file = st.file_uploader("Sélectionner le fichier .ped", type="ped")
    map_file = st.file_uploader("Sélectionner le fichier .map", type="map")

    if ped_file and map_file:
        ped_path = Path("data/input/complex_simulation/user_input.ped")
        map_path = Path("data/input/complex_simulation/user_input.map")
        ped_path.write_bytes(ped_file.getbuffer())
        map_path.write_bytes(map_file.getbuffer())
        st.success("Fichiers .ped et .map chargés.")
    else:
        ped_path = Path(DEFAULT_PED)
        map_path = Path(DEFAULT_MAP)
        st.info("Utilisation des fichiers par défaut.")

    if st.button("🚀 Lancer l’analyse complète"):
        launch_pipeline()

elif choice == "📊 Résultats & Rapport":
    st.header("📊 Résultats de l'analyse et rapport")

    if Path(REPORT_PATH).exists():
        st.subheader("Aperçu rapide du rapport HTML")
        components.html(Path(REPORT_PATH).read_text(), height=600, scrolling=True)
        st.markdown(f"[📄 Ouvrir dans un nouvel onglet]({REPORT_PATH})", unsafe_allow_html=True)
    else:
        st.warning("Le rapport HTML n'a pas encore été généré.")

    if Path(SUMMARY_CSV).exists():
        df_summary = pd.read_csv(SUMMARY_CSV)
        st.subheader("Résumé des résultats")
        st.dataframe(df_summary)
    else:
        st.info("Aucun fichier de résumé (summary.csv) disponible.")

elif choice == "📚 Wiki":
    st.header("📚 Documentation Wiki")

    if Path(README_FILE).exists():
        st.subheader("README du projet")
        st.markdown(Path(README_FILE).read_text(), unsafe_allow_html=True)
    else:
        st.warning("README.md introuvable.")

    if Path(WIKI_DIR).exists():
        md_files = list(Path(WIKI_DIR).glob("*.md"))
        selected_md = st.selectbox("Choisir un module à afficher :", [f.name for f in md_files])
        if selected_md:
            content = Path(WIKI_DIR, selected_md).read_text()
            st.markdown(f"### {selected_md}", unsafe_allow_html=True)
            st.markdown(content, unsafe_allow_html=True)
    else:
        st.warning("Dossier 'Module_WIKI' non trouvé.")
