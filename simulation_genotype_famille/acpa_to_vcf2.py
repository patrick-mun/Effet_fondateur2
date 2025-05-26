import pandas as pd
import os
from glob import glob
from collections import defaultdict

# Répertoire contenant les fichiers ACPA
input_dir = "data/input/complex_simulation/acpa_samples"
acpa_files = sorted(glob(os.path.join(input_dir, "*.txt")))

# Dictionnaires pour collecter les informations
snp_positions = {}
snp_alleles = defaultdict(set)

# Lecture brute des fichiers ACPA sans entête structuré
for file_path in acpa_files:
    try:
        df = pd.read_csv(
            file_path,
            sep=None,
            engine="python",
            comment='#',
            header=None
        )
        print(f"✅ Fichier chargé : {file_path} ({df.shape[0]} lignes)")

    except Exception as e:
        print(f"❌ Erreur lors de la lecture de {file_path} : {e}")
        continue

    for _, row in df.iterrows():
        try:
            rsid = str(row[6]).strip()
            chrom = str(row[7]).strip()
            pos = int(float(str(row[8]).strip()))
            genotype = str(row[5]).strip().upper()
        except Exception:
            continue

        if not rsid or not genotype:
            continue

        if chrom == "19" and pos:
            snp_positions[rsid] = pos

        for base in genotype:
            if base in {"A", "C", "G", "T"}:
                snp_alleles[rsid].add(base)

# Créer deux fichiers :
# 1. VCF avec les positions valides
# 2. Liste des rsID sans position valide à corriger manuellement

vcf_valid_records = []
vcf_missing_positions = []

for rsid, alleles in snp_alleles.items():
    alleles = sorted(list(alleles))
    if len(alleles) == 0:
        continue
    elif len(alleles) == 1:
        ref = alleles[0]
        alt = ref  # Duplique l'allèle observé pour compatibilité VCF
    else:
        ref, alt = alleles[:2]

    if rsid in snp_positions:
        vcf_valid_records.append({
            "CHROM": "19",
            "POS": snp_positions[rsid],
            "ID": rsid,
            "REF": ref,
            "ALT": alt,
            "QUAL": ".",
            "FILTER": "PASS",
            "INFO": "."
        })
    else:
        vcf_missing_positions.append({
            "ID": rsid,
            "REF": ref,
            "ALT": alt
        })

# Écriture du VCF avec positions valides
os.makedirs("data/output/complex_simulation", exist_ok=True)

if vcf_valid_records:
    vcf_valid_df = pd.DataFrame(vcf_valid_records).sort_values(by="POS")
    vcf_valid_path = "data/output/complex_simulation/acpa_chr19_valid.vcf"
    with open(vcf_valid_path, "w") as f:
        f.write("##fileformat=VCFv4.2\n")
        f.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")
        vcf_valid_df.to_csv(f, sep='\t', index=False, header=False)
    print(f"✅ VCF écrit : {vcf_valid_path}")
else:
    print("❌ Aucun SNP avec position valide trouvé.")

# Écriture du fichier avec les rsID manquants
vcf_missing_df = pd.DataFrame(vcf_missing_positions)
missing_path = "data/output/complex_simulation/acpa_chr19_missing_positions.tsv"
vcf_missing_df.to_csv(missing_path, sep='\t', index=False)
print(f"📄 RSID sans position : {missing_path}")