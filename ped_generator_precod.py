import pandas as pd
import numpy as np
import os
import subprocess
from random import choice
from collections import defaultdict

# Étape 1 : Charger le VCF comme source de SNPs
vcf_path = "data/output/complex_simulation/acpa_chr19_valid.vcf"
df_vcf = pd.read_csv(vcf_path, sep="\t", comment="#", header=None,
                     names=["CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO"])

SNP_IDS = df_vcf["ID"].tolist()
SNP_ALLELES = list(zip(df_vcf["REF"], df_vcf["ALT"]))
POSITIONS = df_vcf["POS"].tolist()
N_SNP = len(SNP_IDS)

# Étape 2 : Demander une mutation cible à positionner dans les SNPs
mutation_pos = int(input(f"Entrez la position exacte (base) de la mutation à dater (entre {POSITIONS[0]} et {POSITIONS[-1]}) : "))
mut_allele1 = input("Allèle 1 (référence) : ").strip().upper()
mut_allele2 = input("Allèle 2 (muté) : ").strip().upper()

if mutation_pos not in POSITIONS:
    # Trouver la position correcte d'insertion
    insertion_index = next((i for i, pos in enumerate(POSITIONS) if pos > mutation_pos), len(POSITIONS))
    mutation_id = f"rsMUT_{mutation_pos}"
    mutation_alleles = (mut_allele1, mut_allele2)
    SNP_IDS.insert(insertion_index, mutation_id)
    SNP_ALLELES.insert(insertion_index, mutation_alleles)
    POSITIONS.insert(insertion_index, mutation_pos)
    print(f"✅ Mutation ajoutée à la position {mutation_pos} avec ID {mutation_id} à l'index {insertion_index}.")
else:
    insertion_index = POSITIONS.index(mutation_pos)

N_SNP = len(SNP_IDS)

# Étape 3 : Génération du fichier .map à partir du VCF
os.makedirs("data/input/complex_simulation", exist_ok=True)
with open("data/input/complex_simulation/genotype_data.map", "w") as f:
    for i in range(N_SNP):
        f.write(f"19 {SNP_IDS[i]} 0 {POSITIONS[i]}\n")

# Étape 4 : Fonctions de génération de génotypes
def simulate_genotype():
    return [(choice(alleles), choice(alleles)) for alleles in SNP_ALLELES]

def inherit_genotype(parent1, parent2):
    return [(choice(g1), choice(g2)) for g1, g2 in zip(parent1, parent2)]

def generate_parent_genotype(gp1, gp2):
    if gp1 and gp2:
        return inherit_genotype(gp1, gp2)
    else:
        return simulate_genotype()

# Étape 5 : Définir la structure des familles dynamiquement
nb_familles = int(input("Combien de familles voulez-vous simuler ? "))

familles = []
for f in range(nb_familles):
    print(f"\nFamille {f+1} :")
    print("Définir la structure des grands-parents :")

    def input_gp_count(role):
        n = int(input(f"Combien de grands-parents {role} ? "))
        statuses = []
        for i in range(n):
            statut = int(input(f"  Statut du grand-parent {role} #{i+1} (1 = sain, 2 = atteint, 3 = porteur) : "))
            statuses.append(statut)
        return n, statuses

    nb_gp_h_paternel, _ = input_gp_count("hommes du côté paternel")
    nb_gp_h_maternel, _ = input_gp_count("hommes du côté maternel")
    nb_gp_f_paternel, _ = input_gp_count("femmes du côté paternel")
    nb_gp_f_maternel, _ = input_gp_count("femmes du côté maternel")

    print("Définir la structure des parents :")
    status_pere = int(input("Statut du père (1 = sain, 2 = atteint, 3 = porteur) : "))
    status_mere = int(input("Statut de la mère (1 = saine, 2 = atteinte, 3 = porteuse) : "))

    print("Définir les enfants :")
    n_enfants = int(input("Combien d'enfants ? "))
    enfants_statuts = []
    for i in range(n_enfants):
        statut = int(input(f"Statut de l'enfant {i+1} (1 = sain, 2 = atteint, 3 = porteur) : "))
        enfants_statuts.append(statut)

    familles.append({
        "gp_h_paternel": nb_gp_h_paternel,
        "gp_h_maternel": nb_gp_h_maternel,
        "gp_f_paternel": nb_gp_f_paternel,
        "gp_f_maternel": nb_gp_f_maternel,
        "status_pere": status_pere,
        "status_mere": status_mere,
        "enfants_statuts": enfants_statuts
    })

# Étape 6 : Définir les témoins
nb_temoins = int(input("Combien de témoins souhaitez-vous générer ? "))
temoins = []
for i in range(nb_temoins):
    geno = simulate_genotype()
    sex = 1 if i % 2 == 0 else 2
    temoins.append(("CTRL", f"CTRL_{i+1}", '0', '0', sex, 1, geno))

# Étape 7 : Sélection validée de la région homozygote
def get_validated_homozygous_region():
    while True:
        start = int(input(f"Position de début de la région homozygote (0 à {N_SNP - 1}) : "))
        end = int(input(f"Position de fin (entre {start + 1} et {N_SNP}) : "))
        if end - start >= 100 and (POSITIONS[end - 1] - POSITIONS[start]) >= 1000000:
            if insertion_index < start or insertion_index >= end:
                print("⚠️ Attention : la mutation ciblée ne se trouve pas dans la région homozygote sélectionnée.")
            return start, end
        print("❌ La région sélectionnée ne respecte pas les critères PLINK (≥100 SNPs, ≥1Mb).")
        for i in range(N_SNP - 100):
            if POSITIONS[i + 99] - POSITIONS[i] >= 1000000:
                print(f"💡 Exemple suggéré : entre {i} et {i + 100} (≈ {(POSITIONS[i + 99] - POSITIONS[i]) / 1000:.1f} kb)")
                break
        print("Veuillez réessayer.")

# Étape 8 : Construction des génotypes individuels
individus = []

for idx, famille in enumerate(familles):
    fid = f"F{idx+1}"

    gph1 = simulate_genotype() if famille['gp_h_paternel'] else []
    gpf1 = simulate_genotype() if famille['gp_f_paternel'] else []
    gph2 = simulate_genotype() if famille['gp_h_maternel'] else []
    gpf2 = simulate_genotype() if famille['gp_f_maternel'] else []

    individus.extend([
        (fid, f"{fid}_GP_HP1", '0', '0', 1, 1, gph1) if gph1 else None,
        (fid, f"{fid}_GP_FP1", '0', '0', 2, 1, gpf1) if gpf1 else None,
        (fid, f"{fid}_GP_HM1", '0', '0', 1, 1, gph2) if gph2 else None,
        (fid, f"{fid}_GP_FM1", '0', '0', 2, 1, gpf2) if gpf2 else None,
    ])

    pere = generate_parent_genotype(gph1, gpf1)
    mere = generate_parent_genotype(gph2, gpf2)

    individus.extend([
        (fid, f"{fid}_PERE", f"{fid}_GP_HP1" if gph1 else '0', f"{fid}_GP_FP1" if gpf1 else '0', 1, famille['status_pere'], pere),
        (fid, f"{fid}_MERE", f"{fid}_GP_HM1" if gph2 else '0', f"{fid}_GP_FM1" if gpf2 else '0', 2, famille['status_mere'], mere)
    ])

    start, end = get_validated_homozygous_region()

    for i, statut in enumerate(famille['enfants_statuts']):
        geno = inherit_genotype(pere, mere)
        if statut == 2:
            for j in range(start, end):
                allele = geno[j][0]
                geno[j] = (allele, allele)
            geno[insertion_index] = (mut_allele2, mut_allele2)
        individus.append((fid, f"{fid}_ENF{i+1}", f"{fid}_PERE", f"{fid}_MERE", 1 if i % 2 == 0 else 2, statut, geno))

individus.extend(temoins)
individus = [i for i in individus if i]

# Étape 9 : Écriture du fichier .ped
with open("data/input/complex_simulation/genotype_data.ped", "w") as f:
    for ind in individus:
        geno_str = ' '.join(f"{a} {b}" for a, b in ind[6])
        f.write(f"{ind[0]} {ind[1]} {ind[2]} {ind[3]} {ind[4]} {ind[5]} {geno_str}\n")

# Étape 10 : Génération du fichier groupes.txt
with open("data/input/complex_simulation/groupes.txt", "w") as f:
    for ind in individus:
        if ind[0].startswith("CTRL"):
            label = "TEMOIN"
        elif ind[5] == 2:
            label = "ATTEINT"
        elif ind[5] == 3:
            label = "PORTEUR"
        else:
            label = "FAMILLE"
        f.write(f"{ind[1]}\t{label}\n")

# Étape 11 : Conversion .ped/.map → .bed/.bim/.fam avec PLINK
print("\n✅ Conversion vers le format binaire PLINK...")
subprocess.run([
    "plink",
    "--file", "data/input/complex_simulation/genotype_data",
    "--make-bed",
    "--out", "data/input/complex_simulation/genotype_data"
])
print("✅ Fichiers .bed/.bim/.fam générés avec succès.")
