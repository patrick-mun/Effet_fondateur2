# ped_generator_precod.py
"""
Pré-code pour la génération d'un fichier .ped et .map simulé
Objectif : simuler une structure familiale complète avec des SNPs diploïdes, une mutation délétère
homozygote chez certains enfants, et intégrer des erreurs ou données manquantes pour tester un pipeline.
"""

# Étape 1 : Définir les constantes générales dynamiquement
BASES = ['A', 'C', 'G', 'T']  # Les 4 bases possibles
N_SNP = int(input("Combien de SNPs souhaitez-vous simuler ? "))

print("Souhaitez-vous définir manuellement les fréquences alléliques ?")
manual_freq = input("(oui/non) : ").strip().lower() == "oui"

from random import choice, choices
import numpy as np
import os
import subprocess

# Étape 2 : Génération de positions irrégulières croissantes pour le fichier .map
positions = np.cumsum(np.random.randint(500, 1500, size=N_SNP))
while positions[-1] < 10_000_000:
    positions = np.cumsum(np.random.randint(500, 1500, size=N_SNP))
os.makedirs("data/input/complex_simulation", exist_ok=True)
with open("data/input/complex_simulation/genotype_data.map", "w") as f:
    for i in range(N_SNP):
        f.write(f"1 rs{i+1} 0 {positions[i]}\n")

# Étape 3 : Fonctions de génération de génotypes
def simulate_genotype():
    return [(choice(alleles), choice(alleles)) for alleles in SNP_ALLELES]

def inherit_genotype(parent1, parent2):
    return [(choice(g1), choice(g2)) for g1, g2 in zip(parent1, parent2)]

def generate_parent_genotype(gp1, gp2):
    if gp1 and gp2:
        return inherit_genotype(gp1, gp2)
    else:
        return simulate_genotype()

if manual_freq:
    SNP_ALLELES = []
    for i in range(N_SNP):
        print(f"SNP {i+1}:")
        a1 = input("  Allèle 1 : ").strip().upper()
        a2 = input("  Allèle 2 : ").strip().upper()
        if a1 not in BASES or a2 not in BASES:
            raise ValueError(f"Allèles non valides : {a1}, {a2}. Choisissez parmi {BASES}.")
        SNP_ALLELES.append((a1, a2))
else:
    SNP_ALLELES = [(choice(BASES), choice(BASES)) for _ in range(N_SNP)]

for i, (a1, a2) in enumerate(SNP_ALLELES):
    if a1 not in BASES or a2 not in BASES:
        raise ValueError(f"Allèle invalide détecté à la position {i+1} : ({a1}, {a2})")

# Étape 4 : Définir la structure des familles dynamiquement
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

    nb_gp_h_paternel, st_gp_h_paternel = input_gp_count("hommes du côté paternel")
    nb_gp_h_maternel, st_gp_h_maternel = input_gp_count("hommes du côté maternel")
    nb_gp_f_paternel, st_gp_f_paternel = input_gp_count("femmes du côté paternel")
    nb_gp_f_maternel, st_gp_f_maternel = input_gp_count("femmes du côté maternel")

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

# Étape 5 : Définir les témoins
nb_temoins = int(input("Combien de témoins souhaitez-vous générer ? "))
temoins = []
for i in range(nb_temoins):
    geno = simulate_genotype()
    sex = 1 if i % 2 == 0 else 2
    temoins.append(("CTRL", f"CTRL_{i+1}", '0', '0', sex, 1, geno))

# Étape 6 : Sélection validée de la région homozygote
def get_validated_homozygous_region():
    while True:
        start = int(input(f"Position de début de la région homozygote (0 à {N_SNP - 1}) : "))
        end = int(input(f"Position de fin (entre {start + 1} et {N_SNP}) : "))
        if end - start >= 100 and (positions[end - 1] - positions[start]) >= 1000000:
            return start, end
        print("❌ La région sélectionnée ne respecte pas les critères PLINK (≥100 SNPs, ≥1Mb).")
        for i in range(N_SNP - 100):
            if positions[i + 99] - positions[i] >= 1000000:
                print(f"💡 Exemple suggéré : entre {i} et {i + 100} (≈ {(positions[i + 99] - positions[i]) / 1000:.1f} kb)")
                break
        print("Veuillez réessayer.")

# Étape 7 : Construction des génotypes individuels
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
        individus.append((fid, f"{fid}_ENF{i+1}", f"{fid}_PERE", f"{fid}_MERE", 1 if i % 2 == 0 else 2, statut, geno))

individus.extend(temoins)
individus = [i for i in individus if i]

# Étape 8 : Écriture du fichier .ped
with open("data/input/complex_simulation/genotype_data.ped", "w") as f:
    for ind in individus:
        geno_str = ' '.join(f"{a} {b}" for a, b in ind[6])
        f.write(f"{ind[0]} {ind[1]} {ind[2]} {ind[3]} {ind[4]} {ind[5]} {geno_str}\n")

# Étape 9 : Génération du fichier groupes.txt
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

# Étape 10 : Conversion .ped/.map → .bed/.bim/.fam avec PLINK
print("\n✅ Conversion vers le format binaire PLINK...")
subprocess.run([
    "plink",
    "--file", "data/input/complex_simulation/genotype_data",
    "--make-bed",
    "--out", "data/input/complex_simulation/genotype_data"
])
print("✅ Fichiers .bed/.bim/.fam générés avec succès.")
