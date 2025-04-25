
# Pipeline – Analyse de l’effet fondateur

Ce pipeline exécute une analyse complète pour démontrer un **effet fondateur** autour d'une mutation génétique, en utilisant des données SNP issues de puces ACPA.

## 🧪 Objectif

- Identifier les segments d'autozygotie (ROH)
- Détecter les liens familiaux cachés (IBD)
- Analyser la structure populationnelle via les SNPs
- Dater l'apparition de la mutation via déséquilibre de liaison (Gamma)
- Générer des graphiques et un rapport PDF

---

## ⚙️ Dépendances

- Python 3.x
- Modules Python : `pandas`, `matplotlib`, `fpdf`
- Logiciels : `PLINK`, `KING`, `Gamma`, `R` (avec les packages `adegenet`, `poppr`, `ape`, `ggtree`)

---

## 📁 Structure des fichiers attendus

Placer vos fichiers d'entrée dans `data/input/` :

| Fichier                   | Description                             | Comment l’obtenir                                |
|--------------------------|-----------------------------------------|---------------------------------------------------|
| `genotype_data.ped`      | Données SNP au format texte             | Export depuis Affymetrix ou conversion VCF → PED  |
| `genotype_data.map`      | Carte des SNPs                          | Généré avec `.ped` via PLINK                      |
| `genotype_data.bed`      | Fichier binaire compressé               | `plink --make-bed`                                |
| `genotype_data.bim`      | Informations SNPs en binaire            | Généré avec `.bed`                                |
| `genotype_data.fam`      | Métadonnées individus                   | Généré avec `.bed`                                |
| `genotype_data.raw`      | Tableau de génotypes 0/1/2              | `plink --recodeA`                                 |
| `genotype_data.map`      | Carte génétique (positions en cM)       | Nécessaire pour Gamma                             |

---

## 🚀 Lancer le pipeline

### Exécution complète
```bash
python pipeline_dock6.py
```

### Options CLI

| Option             | Description                                                             |
|--------------------|-------------------------------------------------------------------------|
| `--version`        | Affiche la version du pipeline                                          |
| `--step`           | Exécute une étape spécifique (`roh`, `ibd`, `adegenet`, `gamma`, `all`) |
| `--debug`          | Affiche les messages de débogage détaillés                              |
| `--resume`         | Réservé pour relancer un run interrompu (à venir)                       |

Exemples :
```bash
python pipeline_dock6.py --step roh
python pipeline_dock6.py --step gamma --debug
```

---

## 📦 Résultats générés

| Dossier              | Contenu                                           |
|----------------------|---------------------------------------------------|
| `output/roh_results` | Segments ROH détectés par PLINK                   |
| `output/ibd_results` | Liens familiaux détectés par KING                 |
| `output/adegenet_results` | Arbres phylogénétiques, matrices de distance |
| `output/gamma_results` | `gamma_input.txt`, graphique PNG, rapport PDF   |

---

## 📄 Rapport final

Un fichier PDF est généré automatiquement avec :
- Un graphique des fréquences alléliques des SNPs informatifs
- Un tableau statistique : moyenne, min, max, nombre total de SNPs

---

## 📌 Bonnes pratiques

- Vérifiez la cohérence entre `.raw` et `.map`
- Supprimez les SNPs non informatifs (invariants)
- Travaillez sur un sous-ensemble de données avant d'appliquer au jeu complet

---

## 🧠 Auteur

Patrick MUNIER  
https://github.com/patrick-mun/Effet_fondateur2/tree/main
Laboratoire de génétique, Projet DOCK6  
Version : 1.0 ALPHA
date : 03/04/2025

---

