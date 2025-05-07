#!/bin/bash

# === CONFIGURATION ===

# Fichier VCF source de gnomAD v3.1.2 (chr19)
VCF_URL="https://storage.googleapis.com/gnomad-public/release/3.1.2/genomes/gnomad.genomes.v3.1.2.sites.chr19.vcf.bgz"
TBI_URL="${VCF_URL}.tbi"

# Fichiers locaux
VCF_FILE="gnomad_chr19.vcf.bgz"
TBI_FILE="${VCF_FILE}.tbi"
META_URL="https://storage.googleapis.com/gnomad-public/release/3.1/ht/genomes/sample_annotations.ht/metadata.tsv.bgz"
META_FILE="metadata.tsv"
SELECTED_IDS="selected_ids.txt"
ACPA_RSIDS="acpa_rsids.txt"  # <=== fournis ce fichier toi-même
OUTPUT="acpa_chr19_filtered.vcf.gz"

# === ÉTAPE 1 : Télécharger VCF + index ===
echo "🔽 Téléchargement du VCF chr19..."
curl -C - -O "$VCF_URL"
curl -C - -O "$TBI_URL"

# === ÉTAPE 2 : Télécharger et décompresser les métadonnées ===
echo "🔽 Téléchargement des métadonnées d’échantillons..."
curl -O "$META_URL"
gunzip -f "$(basename $META_URL)"

# === ÉTAPE 3 : Extraire les IDs des populations d’intérêt ===
echo "🔎 Extraction des IDs pour les populations ciblées..."
awk -F'\t' 'tolower($0) ~ /french|orcadian|tuscan|sardinian|han|dai|dinka|maasai|malagasy|hadza|somali|telugu|bengali|gujarati|brahmin|mala|xibo|hezhen/' "$META_FILE" | cut -f1 > "$SELECTED_IDS"

echo "✔️ $(wc -l < $SELECTED_IDS) individus sélectionnés."

# === ÉTAPE 4 : Extraire les variantes ACPA pour ces individus ===
if [ ! -f "$ACPA_RSIDS" ]; then
  echo "❌ Fichier $ACPA_RSIDS introuvable. Fournis ta liste (~3500 rsID), un par ligne."
  exit 1
fi

echo "🚀 Extraction des variants ACPA pour les individus sélectionnés..."
bcftools view -i "ID=@$ACPA_RSIDS" -S "$SELECTED_IDS" -Oz -o "$OUTPUT" "$VCF_FILE"

echo "✅ Fichier final : $OUTPUT"
