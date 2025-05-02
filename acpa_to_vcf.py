# acpa_to_vcf.py
import pandas as pd

infile = "data/input/complex_simulation/chr19_rs_82206960_(CytoScan750K_Array).cy750K.txt"
outfile = "data/input/complex_simulation/acpa_chr19.vcf"

df = pd.read_csv(infile, sep="\t", comment="#")
df = df[["dbSNP RS ID", "Chromosomal Position", "Forward Strand Base Calls"]].dropna()

with open(outfile, "w") as f:
    # En-tête VCF standard
    f.write("##fileformat=VCFv4.2\n")
    f.write("##source=ACPA_CytoScan750K_Converted\n")
    f.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")

    for _, row in df.iterrows():
        rsid = row["dbSNP RS ID"]
        pos = int(row["Chromosomal Position"])
        alleles = row["Forward Strand Base Calls"]
        ref = alleles[0]
        alt = alleles[1] if len(alleles) > 1 and alleles[1] != alleles[0] else "N"
        f.write(f"19\t{pos}\t{rsid}\t{ref}\t{alt}\t.\t.\t.\n")

print("✅ Fichier VCF généré :", outfile)
