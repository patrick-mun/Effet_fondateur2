# 🧬 DOCK6 Pipeline – Analyse de l’effet fondateur

## 📌 Objectif du pipeline

Le pipeline **DOCK6** automatise l’analyse génétique d’un **effet fondateur** à partir de données SNP issues de puces de génotypage. Il combine des outils éprouvés (PLINK, KING, Gamma, R/Adegenet) pour :

- Détecter les segments d’autozygotie (ROH)
- Identifier les liens familiaux cachés (IBD)
- Représenter la structure populationnelle
- Estimer la date d’apparition d’une mutation
- Générer graphiques et rapports

---

## 📦 Dépendances et environnement

### Langages et outils requis
- Python ≥ 3.8
- R (avec packages : `adegenet`, `poppr`, `ape`)
- Binaries : `PLINK`, `KING`, `Gamma`

### Modules Python
```bash
pip install pandas matplotlib fpdf
```

---

## 📁 Structure du projet

```
DOCK6_PROJECT/
├── data/
│   ├── input/         ← Fichiers PED, MAP, RAW, BED, etc.
│   └── output/        ← Résultats des différentes étapes
├── scripts/           ← Scripts R et commandes shell
├── pipeline_dock6.py  ← Script principal Python
├── dock6_app.py       ← Interface Streamlit (optionnelle)
└── README.md          ← Ce document
```

---

## 🔄 Étapes du pipeline

| Étape        | Outil       | Description principale                             |
|--------------|-------------|----------------------------------------------------|
| ROH          | PLINK       | Détection des segments d’autozygotie               |
| IBD          | KING        | Identification des relations de parenté            |
| Structure    | Adegenet    | Analyse des clusters et arbres phylogénétiques     |
| Gamma        | Gamma       | Estimation de la date d’apparition de la mutation  |

---

## 🚀 Exécution du pipeline

### Lancer toutes les étapes automatiquement
```bash
python pipeline_dock6.py
```

### Lancer une étape spécifique
```bash
python pipeline_dock6.py --step roh
python pipeline_dock6.py --step gamma --debug
```

Options disponibles :
- `--step` : `roh`, `ibd`, `adegenet`, `gamma`, `all`
- `--debug` : active le mode verbeux (logs détaillés)
- `--resume` : option future pour relancer un run interrompu

---

## 📤 Résultats générés

| Dossier              | Contenu                                               |
|----------------------|-------------------------------------------------------|
| `roh_results/`       | Fichiers `.hom` (segments homozygotes PLINK)          |
| `ibd_results/`       | Table `.kin0` (parenté KING)                          |
| `adegenet_results/`  | Arbres, graphiques DAPC, matrices de distance         |
| `gamma_results/`     | Fichier `.txt`, graphique de fréquence, rapport PDF   |

---

## 🧾 Rapport final

Le pipeline peut générer un **rapport PDF** intégrant :
- Le graphique Gamma (`frequence_allelique.png`)
- Les statistiques de fréquence (moyenne, min, max)
- Les arbres phylogénétiques (via Adegenet)

---

## 📊 Interface Web (optionnelle)

L’interface Streamlit (`dock6_app.py`) permet de :
- Uploader les fichiers `.raw` et `.map`
- Choisir une étape à exécuter
- Visualiser les résultats par onglet
- Télécharger le rapport final

Lancement local :
```bash
streamlit run dock6_app.py
```

---

## ✅ Bonnes pratiques

- Vérifier la correspondance entre `.raw` et `.map`
- Supprimer les SNPs non informatifs avec PLINK (`--maf`, `--geno`)
- Tester d’abord le pipeline sur un sous-échantillon réduit

---

## 🧠 Auteurs et informations

**Auteur principal** : Patrick MUNIER  
Laboratoire de génétique – Projet DOCK6  
Version : 1.0 ALPHA  
Date : 03/04/2025

---

📚 **Besoin de plus de détails sur chaque étape ?** Consultez les modules dédiés pour une explication complète des outils utilisés dans le pipeline :
👉 [PLINK](./module/Module_PLINK.md) | [KING](./module/Module_king.md) | [Gamma](./module/Module_gamma.md) | [R/Adegenet](./module/Module_adgenet.md)

