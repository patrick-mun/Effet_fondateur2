# 🧬 Wiki EFFET FONDATEUR – Module : R / Adegenet

## 📌 Contexte d'utilisation

Ce module décrit **Adegenet** et les outils R associés (poppr, ape), utilisés dans le pipeline EFFET FONDATEUR pour :

- Réaliser une analyse de structure populationnelle (clustering DAPC)
- Générer un arbre phylogénétique basé sur la distance génétique
- Obtenir une visualisation claire des relations entre individus

---

## 📖 Récapitulatif du document

- **🔍 Présentation d’Adegenet et des packages R associés**
- **📁 Formats de fichiers supportés**
- **🔄 Fichiers utilisés dans le pipeline DOCK6**
- **⚙️ Script R utilisé**
- **📊 Exemples de fichiers de sortie**
- **🔍 Interprétation des résultats**
- **🧪 Intégration dans DOCK6**
- **✅ Avantages et ⚠️ Limites**
- **💡 Bonnes pratiques**

---

## 🔍 Présentation d’Adegenet

**Nom complet** : Adegenet (Analysis of Genetic Data)

**Langage** : R

**Packages associés** :
- `adegenet` : clustering et visualisation
- `poppr` : distances génétiques
- `ape` : arbres phylogénétiques

### Objectifs et utilisations principales

- Clustering par DAPC (Discriminant Analysis of Principal Components)
- Calcul de distances génétiques
- Génération d’arbres phylogénétiques (type Neighbour Joining)

---

## 📁 Formats de fichiers supportés

Adegenet accepte :
- `.raw` : format PLINK recodé (0/1/2)
- `.gen` : format Genepop (optionnel)

Fichier requis : `genotype_data.raw`

---

## 🔄 Fichiers utilisés dans le pipeline EFFET FONDATEUR

### 📥 Fichier d’entrée :
```
Nom : genotype_data.raw
Format : fichier texte tabulé (.raw), issu de la commande PLINK --recodeA
Contenu :
- FID, IID, sexe, phénotype
- Colonnes de génotypes codés 0/1/2 pour chaque SNP
Exemple (extrait) :
FID IID SEX rs123_A rs124_G rs125_C rs126_T
FAM001 IND001 1 2 1 0 2
FAM001 IND002 2 1 2 1 1
```
Ce fichier sert de base à l’import dans R pour l’analyse DAPC, les matrices de distance, et l’arbre phylogénétique.

### 📤 Fichiers de sortie :
- `dapc.png` : graphique de l’analyse DAPC (structure)
- `phylo_tree.png` : arbre phylogénétique
- `tree.newick` : arbre au format Newick
- `distance_matrix.csv` : matrice des distances génétiques

---

## ⚙️ Script R utilisé

Fichier : `scripts/adegenet_analysis.R`

Voici le script avec des commentaires détaillés pour chaque étape :

```r
# Chargement des bibliothèques nécessaires
library(adegenet)   # Pour l'analyse DAPC et la lecture des données PLINK
library(poppr)      # Pour le calcul des distances génétiques
library(ape)        # Pour la génération de l'arbre phylogénétique

# Lecture des données au format PLINK .raw
# Ce fichier doit être généré par PLINK avec l'option --recodeA
data <- read.PLINK("../data/input/genotype_data.raw")

# Analyse de structure populationnelle par DAPC
dapc1 <- dapc(data)  # Calcul du clustering DAPC
png("../data/output/adegenet_results/dapc.png")
scatter(dapc1)       # Tracé des clusters avec les composantes discriminantes
dev.off()

# Calcul de la distance génétique pair-à-pair
# Transformation des données génétiques pour en extraire une matrice
# de distance Euclidienne entre individus

# Conversion des données pour obtenir une matrice numérique
matrix_data <- tab(data)
distgen <- dist(matrix_data)

# Export de la matrice de distance en CSV
write.csv(as.matrix(distgen), "../data/output/adegenet_results/distance_matrix.csv")

# Construction d'un arbre phylogénétique Neighbour-Joining à partir des distances
tree <- nj(distgen)

# Visualisation de l'arbre en format radial
png("../data/output/adegenet_results/phylo_tree.png")
plot(tree, typ='fan', cex=0.7)
dev.off()

# Export de l'arbre au format Newick (utilisable par d'autres logiciels)
write.tree(tree, file="../data/output/adegenet_results/tree.newick")
```library(adegenet)
library(poppr)
library(ape)

data <- read.PLINK("../data/input/genotype_data.raw")

dapc1 <- dapc(data)
png("../data/output/adegenet_results/dapc.png")
scatter(dapc1)
dev.off()

distgen <- dist(tab(data))
write.csv(as.matrix(distgen), "../data/output/adegenet_results/distance_matrix.csv")

tree <- nj(distgen)
png("../data/output/adegenet_results/phylo_tree.png")
plot(tree, typ='fan', cex=0.7)
dev.off()

write.tree(tree, file="../data/output/adegenet_results/tree.newick")
```

---

## 📊 Exemples de fichiers de sortie

### 🎨 Visualisation des résultats attendus

#### 1. Structure DAPC
Représentation typique attendue dans `dapc.png` :

- Axe X : composantes discriminantes (ex : discriminant 1)
- Axe Y : composantes secondaires ou variation intra-groupe
- Couleurs : groupes génétiques détectés
- Points : individus

*(exemple visuel à insérer ici : `dapc.png`)*

#### 2. Arbre phylogénétique
L’arbre `phylo_tree.png` représente la proximité génétique entre individus. 
Les branches plus courtes indiquent une proximité forte.

- Forme : cercle (type 'fan')
- Noeuds : regroupements d’individus similaires
- Étiquettes : noms des individus (`IID`)

*(exemple visuel à insérer ici : `phylo_tree.png`)*


### `dapc.png`
- Visualisation de la structure populationnelle
- Chaque point représente un individu, les couleurs définissent les clusters

### `phylo_tree.png`
- Arbre généalogique basé sur la distance génétique entre individus

### `distance_matrix.csv`
```
,IND001,IND002,IND003
IND001,0.000,0.251,0.342
IND002,0.251,0.000,0.187
IND003,0.342,0.187,0.000
```

### `tree.newick`
- Format arbre : `((IND001:0.2,IND002:0.2):0.3,IND003:0.5);`

---

## 🔍 Interprétation des résultats

- Le DAPC met en évidence des **groupes génétiques distincts**, révélateurs de l’effet fondateur
- L’arbre phylogénétique permet d’identifier visuellement les **proches génétiques**
- Les matrices de distances complètent les relations visuelles

*(Insérez ici les figures `dapc.png` et `phylo_tree.png` pour une interprétation visuelle)*

---

## 🧪 Intégration dans DOCK6

Adegenet est exécuté via :
```bash
Rscript scripts/adegenet_analysis.R
```

Les résultats sont placés dans `output/adegenet_results/` et intégrés dans :
- Le rapport PDF global
- L’interface Streamlit (onglet Structure)

---

## ✅ Avantages

- Visualisation intuitive de la structure génétique
- Compatible avec les formats issus de PLINK
- Repose sur des méthodes statistiques robustes (DAPC, NJ)

## ⚠️ Limites

- Difficile à interpréter avec des jeux de données trop petits
- Sensible aux données manquantes ou aux erreurs de format

---

## 💡 Bonnes pratiques

- Toujours générer le fichier `.raw` via `plink --recodeA`
- Vérifier la cohérence des identifiants dans `.raw`
- Ne pas forcer un nombre de clusters trop élevé dans DAPC
- Inspecter visuellement les arbres pour détecter les anomalies

