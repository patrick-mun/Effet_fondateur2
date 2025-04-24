#!/bin/bash

# Étape 1 : ROH avec PLINK (paramètres ajustés)
plink --file ../data/input/genotype_data \
      --homozyg \
      --homozyg-snp 30 \
      --homozyg-kb 300 \
      --homozyg-density 100 \
      --homozyg-gap 1000 \
      --out ../data/output/roh_results/roh

# Étape 2 : IBD avec KING
king -b ../data/input/genotype_data.bed --kinship --prefix ../data/output/ibd_results/ibd

# Étape 4 : Datation avec Gamma
Gamma -i ../data/input/genotype_data.gamma_input.txt -o ../data/output/gamma_results/gamma_output.txt
