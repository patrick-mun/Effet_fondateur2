# generate_rsid_list.py
import pandas as pd

# Chemin vers le fichier CytoScan de la puce ACPA
acpa_file = "data/input/complex_simulation/chr19_rs_82206960_(CytoScan750K_Array).cy750K.txt"

# Chargement et extraction des rsID uniques
df = pd.read_csv(acpa_file, sep="\t", comment="#")
rsids = df["dbSNP RS ID"].dropna().unique()

# Écriture dans un fichier texte
with open("data/input/complex_simulation/acpa_rsids.txt", "w") as f:
    for rsid in rsids:
        f.write(f"{rsid}\n")

print("✅ Fichier acpa_rsids.txt généré.")
