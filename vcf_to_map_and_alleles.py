# vcf_to_map_and_alleles.py
import pandas as pd
import os
import gzip

vcf_path = "data/input/complex_simulation/acpa_chr19_filtered.vcf.gz"
map_path = "data/input/complex_simulation/genotype_data.map"
alleles_out = "data/input/complex_simulation/snp_alleles.py"

snp_ids = []
positions = []
alleles = []

with gzip.open(vcf_path, "rt") as f:
    for line in f:
        if line.startswith("#"):
            continue
        fields = line.strip().split("\t")
        chrom, pos, rsid, ref, alt = fields[:5]
        if rsid == "." or alt == ".":
            continue
        alt_allele = alt.split(",")[0]  # gérer les multi-alleliques
        snp_ids.append(rsid)
        positions.append(pos)
        alleles.append((ref, alt_allele))

# Écriture du .map
with open(map_path, "w") as f:
    for rsid, pos in zip(snp_ids, positions):
        f.write(f"19 {rsid} 0 {pos}\n")

# Générer un fichier Python contenant SNP_ALLELES
with open(alleles_out, "w") as f:
    f.write("# Auto-generated SNP_ALLELES list from filtered VCF\n")
    f.write("SNP_ALLELES = [\n")
    for a1, a2 in alleles:
        f.write(f"    ('{a1}', '{a2}'),\n")
    f.write("]\n")

print(f"✅ {len(alleles)} SNPs exportés vers :")
print(f"- .map : {map_path}")
print(f"- SNP_ALLELES (Python) : {alleles_out}")
