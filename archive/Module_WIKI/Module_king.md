# 🧬 Wiki DOCK6 – Module : KING

## 📌 Contexte d'utilisation

Ce module présente **KING**, un logiciel essentiel dans le pipeline DOCK6, utilisé pour identifier les relations de parenté (IBD) entre individus à partir de données SNP. Il est particulièrement utile pour :

- Détecter des individus apparentés sans lien familial documenté
- Valider ou découvrir des regroupements autour d’un ancêtre commun
- Renforcer les résultats obtenus avec les segments ROH (PLINK) et les fréquences alléliques (Gamma)

---

## 📖 Récapitulatif du document

- **🔍 Présentation de KING**
  - Objectifs et cas d’usage
- **📁 Formats de fichiers supportés**
  - Fichiers binaires PLINK requis
- **🔄 Fichiers utilisés dans le pipeline DOCK6**
  - Fichiers d’entrée et sorties générées
- **⚙️ Commande utilisée (DOCK6)**
  - Explication détaillée de chaque paramètre
- **📊 Exemple de fichier de sortie**
  - Interprétation des colonnes du fichier `.kin0`
- **🔍 Interprétation des résultats**
  - Cas d’étude et seuils de parenté
- **🧪 Intégration dans DOCK6**
  - Emplacement dans le pipeline et automatisation
- **✅ Avantages et ⚠️ Limites**
  - Ce qu’il faut savoir
- **💡 Bonnes pratiques**
  - Recommandations d’usage

---

## 🔍 Présentation de KING

**Nom complet** : KING (Kinship-based INference for GWAS)  
**Site officiel** : [https://www.kingrelatedness.com/](https://www.kingrelatedness.com/)

### Objectifs et utilisations principales
- Détection des relations de parenté jusqu’au 5ᵉ degré
- Vérification des liens familiaux dans les données de génotypage
- Mise en évidence de structures familiales non déclarées

---

## 📁 Formats de fichiers supportés

KING utilise le trio de fichiers binaires générés par PLINK :

- `.bed` : génotypes binaires
- `.bim` : liste et positions des SNPs
- `.fam` : métadonnées des individus

---

## 🔄 Fichiers utilisés dans le pipeline DOCK6

### 📥 Fichier d’entrée :
```
Contenu : fichiers binaires PLINK
- genotype_data.bed
- genotype_data.bim
- genotype_data.fam
```
Ces fichiers sont générés automatiquement à partir des données `.ped/.map` grâce à PLINK avec l’option `--make-bed`.

### 📤 Fichier de sortie :
```
Nom : ibd.kin0
Contenu : tableau listant toutes les paires d’individus avec leur coefficient de parenté estimé.
Utilité : identifier les relations de type frères/sœurs, parent-enfant, cousins, etc.
```
---

## ⚙️ Commande utilisée (DOCK6)

### 💡 Cette commande lance le calcul de parenté à partir des fichiers binaires générés par PLINK.
```bash
king -b genotype_data.bed --kinship --prefix ibd
```
🔗 [Documentation officielle KING](https://www.kingrelatedness.com/manual.html)

### Explication des options :
- `-b` : spécifie le fichier `.bed` (le programme retrouve automatiquement les fichiers `.bim` et `.fam`)
- `--kinship` : calcule les coefficients de parenté entre toutes les paires
- `--prefix` : définit le préfixe des fichiers de sortie (`ibd.kin0`, `ibd.log`, etc.)

---

## 📊 Exemple de fichier `ibd.kin0`
```
FID1 IID1 FID2 IID2 RT    HetHet  IBS0  Kinship
FAM1 IND001 FAM1 IND002 FS    1284    43    0.2510
FAM2 IND003 FAM2 IND004 UN    1203    78    0.0801
FAM2 IND003 FAM2 IND005 UN    1220    65    0.1213
```

**Colonnes :**
- `FID1`, `FID2` : identifiants familiaux des deux individus comparés
- `IID1`, `IID2` : identifiants individuels
- `RT` : type de relation estimée (UN = inconnue, FS = frères/soeurs, PO = parent-enfant...)
- `HetHet` : nombre de SNPs hétérozygotes chez les deux individus
- `IBS0` : nombre de SNPs pour lesquels aucun allèle n’est partagé (Identical By State = 0)
- `Kinship` : coefficient de parenté estimé (ex. 0.25 = frères/soeurs, 0.125 = demi-frères/sœurs, etc.)
---

## 🔍 Interprétation des résultats

Dans une analyse d’effet fondateur, KING permet de :

- Détecter une parenté insoupçonnée entre individus d’un même cluster ROH
- Quantifier la proximité génétique entre porteurs d’une mutation
- Appuyer l’hypothèse d’un fondateur commun à travers des coefficients > 0.1

*(Insérez ici une heatmap ou matrice de parenté illustrant les relations entre individus)*

---

## 🧪 Intégration dans DOCK6

KING est exécuté automatiquement par `run_pipeline.py` via le module
`scripts/ibd.py`. Les résultats sont enregistrés dans le sous-dossier `ibd/`
du dossier de sortie de l'analyse, puis intégrés aux rapports.

---

## ✅ Avantages
- Rapide, efficace, adapté à de grands jeux de données
- Prend en charge les fichiers PLINK directement
- Estimation précise jusqu’au 5ᵉ degré de parenté

## ⚠️ Limites
- Sous-estime parfois la parenté dans les populations très homogènes
- Sensible aux SNPs de mauvaise qualité ou aux taux de génotypage faibles

---

## 💡 Bonnes pratiques
- Filtrer les SNPs avec PLINK avant d’utiliser KING (`--geno`, `--maf`)
- Utiliser au moins 10000 SNPs informatifs avec MAF > 0.01
- Vérifier visuellement les matrices de parenté avec Python ou R (heatmap, dendrogramme)
- Toujours croiser les résultats avec ceux de PLINK et Gamma


