# 🧜‍♂️ Wiki DOCK6 – Module : PLINK

## 📌 Contexte d'utilisation

Ce module décrit **PLINK**, un logiciel clé dans le pipeline DOCK6 pour la détection d’effet fondateur. Il est principalement utilisé pour :

- Convertir les données de génotypage en formats compatibles avec les autres outils
- Détecter les segments d’homozygotie (ROH)
- Générer les fichiers nécessaires pour KING, Gamma, et l’analyse R/adegenet

---

## 📖 Récapitulatif du document

- **🔍 Présentation de PLINK**
  - Objectifs et utilisations principales
- **📁 Formats de fichiers supportés**
  - Format texte PED/MAP
  - Format binaire BED/BIM/FAM
  - Autres formats compatibles
- **🔄 Formats utilisés dans le pipeline DOCK6**
  - Fichiers d'entrée
  - Fichiers de sortie générés
- **📊 Exemples détaillés de fichiers**
  - PED, MAP, FAM, BIM, RAW
- **⚙️ Commandes utiles (DOCK6)**
  - Exemples pratiques et paramètres détaillés
- **🔍 Visualisation de fichiers intermédiaires**
  - Exemples réels de fichiers générés
- **📊 Comparaison avant/après filtrage**
  - Effets du nettoyage sur les données
- **🧼 Stratégie de nettoyage et pré-traitement**
  - Méthodes pour assurer la qualité des données
- **🧬 Interprétation des résultats**
  - ROH et effet fondateur
  - Gamma (fréquences alléliques)
  - KING (détection de parenté)
- **🧪 Intégration dans DOCK6**
  - Étapes automatiques via `run_pipeline.py`
- **✅ Avantages et ⚠️ Limites**
  - Forces et précautions d'utilisation
- **💡 Astuce**
  - Conseils pratiques pour une utilisation optimale

---

## 🔍 Présentation de PLINK

**Nom complet** : PLINK ("Whole genome association analysis toolset")\
**Site officiel** : [https://www.cog-genomics.org/plink/](https://www.cog-genomics.org/plink/)\
**Version recommandée** : PLINK 1.9 ou PLINK 2.0

### Objectifs et utilisations principales

- Analyses d'association génétique (GWAS)
- Étude de l'homozygotie (segments ROH)
- Filtrage et manipulation des SNPs
- Calcul de fréquences alléliques
- Conversion de formats (PED ↔ BED, RAW, VCF, etc.)

---

## 📁 Formats de fichiers supportés

### ➔ Format texte PED/MAP

- `.ped` : génotypes individuels (SNPs en paires alléliques)
- `.map` : carte des SNPs (position génomique)

### ➔ Format binaire BED/BIM/FAM

- `.bed` : génotypes compressés
- `.bim` : liste et positions des SNPs
- `.fam` : métadonnées des individus

### ➔ Autres formats

- `.raw` : table SNPs codés 0/1/2, pour R et Gamma
- `.vcf`, `.tped`, `.tfam` : formats compatibles via conversion

---

## 🔄 Formats utilisés dans le pipeline DOCK6

### 📅 Fichiers d’entrée :

- `genotype_data.ped` / `genotype_data.map`
- `genotype_data.bed` / `.bim` / `.fam`
- `genotype_data.raw` : pour R et Gamma

### 📄 Fichiers de sortie générés :

- `roh.hom` : segments ROH détectés
- `genotype_data.raw` : table pour analyse
- `.bed/.bim/.fam` : pour KING
- `gamma_input.txt` : dérivé à partir de `.raw` et `.map`

Ces fichiers assurent la transition vers les autres logiciels du pipeline : KING, R/adegenet et Gamma.

---

## 📊 Exemples de fichiers

### 1. Fichier PED

```
FID IID PAT MAT SEX PHENOTYPE SNP1 SNP2 SNP3 SNP4
FAM001 IND001 0 0 1 1 A A G G C C T T
FAM001 IND002 0 0 2 1 A G G G C T T T
FAM002 IND003 0 0 1 2 A A G G T T T T
```

**Colonnes :**

- `FID` : Identifiant de la famille (ex : FAM001)
- `IID` : Identifiant individuel (ex : IND001)
- `PAT`, `MAT` : Identifiants des parents (0 si inconnu)
- `SEX` : Sexe (1 = mâle, 2 = femelle)
- `PHENOTYPE` : Valeur du phénotype (1 = sain, 2 = malade, -9 = valeur manquante)
- `SNPx` : Allèles pour chaque SNP (codés par lettres nucléotidiques, ex : A/G)èles pour chaque SNP (codés par lettres nucléotidiques)

### 2. Fichier MAP

```
CHR SNP_ID GEN_DIST PHYS_POS
1 rs123 0.0 10100
1 rs124 0.0 10300
1 rs125 0.0 10500
1 rs126 0.0 10800
```

**Colonnes :**

- `CHR` : Numéro du chromosome
- `SNP_ID` : Identifiant unique du SNP
- `GEN_DIST` : Distance génétique (en centiMorgans, souvent 0)
- `PHYS_POS` : Position physique sur le chromosome (en paires de bases)

### 3. Fichier FAM

```
FID IID PAT MAT SEX PHENOTYPE
FAM001 IND001 0 0 1 1
FAM001 IND002 0 0 2 1
FAM002 IND003 0 0 1 2
```

**Colonnes :**

- `FID` : Identifiant de la famille
- `IID` : Identifiant individuel
- `PAT`, `MAT` : Identifiants des parents
- `SEX` : Sexe biologique (1 = mâle, 2 = femelle)
- `PHENOTYPE` : Code du phénotype (état de santé ou groupe d’étude)

### 4. Fichier BIM

```
CHR SNP_ID CM POS ALLELE1 ALLELE2
1 rs123 0 10100 A G
1 rs124 0 10300 G T
1 rs125 0 10500 C T
1 rs126 0 10800 T A
```

**Colonnes :**

- `CHR` : Chromosome
- `SNP_ID` : Identifiant du SNP
- `CM` : Distance génétique cumulée (en cM)
- `POS` : Position physique du SNP sur le chromosome
- `ALLELE1`, `ALLELE2` : Les deux allèles possibles pour ce SNP

### 5. Fichier RAW

```
FID IID PAT MAT SEX PHENOTYPE rs123_A rs124_G rs125_C rs126_T
FAM001 IND001 0 0 1 1 2 1 0 2
FAM001 IND002 0 0 2 1 1 2 1 1
FAM002 IND003 0 0 1 2 2 0 2 2
```

**Colonnes :**

- `FID` : Identifiant de la famille (ex : FAM001)
- `IID` : Identifiant de l'individu (ex : IND001)
- `PAT`, `MAT` : Identifiants des parents (0 si inconnus)
- `SEX` : Sexe biologique (1 = mâle, 2 = femelle)
- `PHENOTYPE` : Statut de l'individu (1 = sain, 2 = malade, -9 = inconnu)
- `rsXXX_Y` : Valeurs codées 0/1/2 pour chaque SNP, représentant respectivement 0, 1 ou 2 copies de l’allèle mineur : Identifiants et métadonnées (voir fichier PED)
- `rsXXX_Y` : Valeurs génotypiques codées 0/1/2 pour chaque SNP, indiquant le nombre d’allèles mineurs

---

## ⚙️ Commandes utiles (DOCK6)

### 💡 Exemple de commande complète pour détection ROH personnalisée

Voici une commande combinant plusieurs paramètres recommandés pour une analyse fine des segments d’homozygotie :

```bash
plink --bfile input_data \
      --homozyg \
      --homozyg-snp 50 \
      --homozyg-kb 1000 \
      --homozyg-density 50 \
      --homozyg-gap 1000 \
      --homozyg-window-snp 50 \
      --homozyg-window-het 1 \
      --homozyg-window-missing 5 \
      --out roh_custom
```

Cette configuration permet de détecter uniquement les segments d’homozygotie fiables, en limitant les erreurs liées aux hétérozygotes et aux données manquantes. Le fichier généré `roh_custom.hom` contiendra les segments filtrés selon ces critères.

---

Chaque commande PLINK ci-dessous permet d'exécuter une tâche spécifique dans le cadre de l'analyse de l'effet fondateur. Elles servent à transformer les formats de fichiers, filtrer les données génétiques, détecter l'homozygotie ou préparer les données pour d'autres outils du pipeline DOCK6. Des explications sont fournies pour chaque paramètre, ainsi que des cas d’usage fréquents dans ce type d’analyse.

Chaque commande PLINK ci-dessous permet d'exécuter une tâche spécifique dans le cadre de l'analyse de l'effet fondateur. Elles servent à transformer les formats de fichiers, filtrer les données génétiques, détecter l'homozygotie ou préparer les données pour d'autres outils du pipeline DOCK6.



**Conversion PED → BED** :
Cette commande convertit les fichiers de génotypage au format texte (.ped/.map) en un format binaire (.bed/.bim/.fam) plus rapide à manipuler pour les étapes suivantes du pipeline.

```bash
plink --file input_data --make-bed --out output_data
```
🔗 [Documentation officielle – make-bed](https://www.cog-genomics.org/plink/1.9/data#make-bed)

**Explication** :

- `--file` : spécifie la base du nom de fichier pour `.ped` et `.map`

- `--make-bed` : génère les fichiers `.bed`, `.bim`, `.fam`

- `--out` : définit le préfixe des fichiers de sortie



**Détection ROH** :
Cette commande détecte les segments d'homozygotie (ROH) chez les individus, utiles pour identifier des signatures d'effet fondateur.

```bash
plink --file input_data --homozyg --out roh_output
```
🔗 [Documentation officielle – homozyg](https://www.cog-genomics.org/plink/1.9/ibc#homozyg)

**Explication** :

- `--homozyg` : lance la détection des segments d’homozygotie
- Produit un fichier `.hom` avec les segments ROH par individu

**Paramètres avancés** (optionnels mais recommandés pour affiner l’analyse ROH) :

- `--homozyg-snp` : nombre minimal de SNPs dans un segment ROH (ex : `--homozyg-snp 50`)
- `--homozyg-kb` : longueur minimale en kb d’un segment ROH (ex : `--homozyg-kb 1000`)
- `--homozyg-density` : densité maximale entre SNPs (ex : `--homozyg-density 50`)
- `--homozyg-gap` : distance max en kb entre deux SNPs consécutifs (ex : `--homozyg-gap 1000`)
- `--homozyg-window-snp` : taille de la fenêtre glissante (ex : `--homozyg-window-snp 50`)
- `--homozyg-window-het` : nombre maximal d’hétérozygotes autorisés par fenêtre (ex : `--homozyg-window-het 1`)
- `--homozyg-window-missing` : nombre de génotypes manquants permis par fenêtre (ex : `--homozyg-window-missing 5`)



**Générer .raw pour R/Adegenet/Gamma** :
Cette commande exporte les données génétiques au format texte (.raw), codées en 0/1/2, nécessaires pour les analyses dans R (adegenet) ou Gamma.

```bash
plink --bfile input_data --recodeA --out output_data
```
🔗 [Documentation officielle – recodeA](https://www.cog-genomics.org/plink/1.9/formats#recodeA)

**Explication** :

- `--bfile` : utilise les fichiers binaires `.bed/.bim/.fam`

- `--recodeA` : exporte un fichier `.raw` avec codage 0/1/2



**Filtrage SNPs** :
Cette commande nettoie les données en supprimant les SNPs de mauvaise qualité ou peu informatifs selon leur taux de données manquantes et leur fréquence.

```bash
plink --bfile input_data --geno 0.05 --maf 0.01 --make-bed --out filtered_data
```
🔗 [Documentation officielle – --geno / --maf](https://www.cog-genomics.org/plink/1.9/filtering#geno)

**Explication** :

- `--geno 0.05` : supprime les SNPs avec plus de 5% de données manquantes

- `--maf 0.01` : supprime les SNPs avec une fréquence allèlique mineure < 1%

- `--make-bed` : génère un nouveau jeu de fichiers binaires filtrés



**Extraire une région chromosomique** :
Cette commande permet d’extraire uniquement les SNPs situés dans une région spécifique d’un chromosome, utile pour cibler une mutation donnée.

```bash
plink --bfile data --chr 1 --from-bp 10000 --to-bp 50000 --make-bed --out chr1_region
```
🔗 [Documentation officielle – extraire une région](https://www.cog-genomics.org/plink/1.9/filtering#chr)

**Explication** :

- Extraction des SNPs sur le chromosome 1 entre les positions 10 000 et 50 000



**Calcul des fréquences alléliques** :
Cette commande calcule la fréquence des allèles dans l’ensemble de la population étudiée.

```bash
plink --bfile data --freq --out freq_output
```
🔗 [Documentation officielle – fréquences alléliques](https://www.cog-genomics.org/plink/1.9/summary#freq)

**Explication** :

- Génère un fichier `.frq` contenant la fréquence de chaque allèle



**Créer une liste de SNPs informatifs (bivariés)** :
Cette commande permet de conserver uniquement les SNPs polymorphes (non invariants) pour les analyses statistiques.

```bash
plink --bfile data --exclude snps_invariants.txt --make-bed --out informative_data
```
🔗 [Documentation officielle – exclure SNPs](https://www.cog-genomics.org/plink/1.9/data#extract)

**Explication** :

- Supprime les SNPs invariants (non polymorphes) listés dans `snps_invariants.txt`



**Obtenir un sous-échantillon d’individus** :
Cette commande permet de restreindre l’analyse à un groupe spécifique d’individus (par exemple une famille ou un sous-groupe).

```bash
plink --bfile data --keep subset_individuals.txt --make-bed --out subset_data
```
🔗 [Documentation officielle – sous-échantillonnage](https://www.cog-genomics.org/plink/1.9/data#keep)

**Explication** :

- Garde uniquement les individus listés dans `subset_individuals.txt`



---

## 🔍 Visualisation de fichiers intermédiaires

Voici quelques exemples de contenu réel que vous pouvez retrouver dans les fichiers générés par PLINK, utiles pour interpréter rapidement les résultats de traitement :

### 🧾 Extrait de `roh_output.hom`

```
FID     IID     CHR     SNP1    SNP2    KB      NSNP    HOM_TYPE
FAM001  IND001  1       rs123   rs126   800     25      1
FAM001  IND002  2       rs200   rs220   1000    30      1
```

- `CHR` : numéro du chromosome
- `SNP1`, `SNP2` : bornes du segment ROH
- `KB` : taille du segment en kilobases
- `NSNP` : nombre de SNPs dans le segment

### 📊 Extrait de `freq_output.frq`

```
CHR     SNP     A1      A2      MAF     NCHROBS
1       rs123   G       A       0.35    198
1       rs124   T       G       0.48    200
```

- `A1` : allèle mineur
- `A2` : allèle majeur
- `MAF` : fréquence allèlique mineure
- `NCHROBS` : nombre de chromosomes observés (2 × nombre d’individus)

### 📑 Extrait de `genotype_data.raw`

```
FID IID SEX rs123_A rs124_G rs125_C rs126_T
FAM001 IND001 1 2 1 0 2
FAM001 IND002 2 1 2 1 1
```

- `FID` : Identifiant de la famille
- `IID` : Identifiant individuel
- `SEX` : Sexe (1 = mâle, 2 = femelle)
- `rsXXX_Y` : Valeurs génotypiques codées 0/1/2 représentant le nombre d’allèles mineurs (0 = homozygote majeur, 1 = hétérozygote, 2 = homozygote mineur)

### 🧬 Extrait de `genotype_data.bim`

```
1 rs123 0 10100 A G
1 rs124 0 10300 G T
1 rs125 0 10500 C T
1 rs126 0 10800 T A
```

- `CHR` : Numéro de chromosome
- `SNP` : Identifiant du SNP
- `CM` : Position génétique (souvent 0 si inconnue)
- `BP` : Position physique en paires de bases
- `ALLELE1`, `ALLELE2` : Les deux allèles possibles

### 🗺️ Extrait de `genotype_data.map`

```
1 rs123 0.0 10100
1 rs124 0.0 10300
1 rs125 0.0 10500
1 rs126 0.0 10800
```

- `CHR` : Chromosome
- `SNP` : Identifiant du SNP
- `GEN_DIST` : Distance génétique (souvent 0)
- `PHYS_POS` : Position physique du SNP (en pb)

---

Pour mieux comprendre l’impact des étapes de filtrage ou de transformation, il peut être utile de visualiser rapidement les fichiers générés par PLINK :

- **Comparer le nombre de SNPs avant/après filtrage** :

```bash
wc -l input_data.bim
wc -l filtered_data.bim
```

- **Visualiser les fréquences alléliques** :

```bash
head freq_output.frq
```

- **Vérifier les segments ROH** :

```bash
less roh_output.hom
```

---

## 📊 Comparaison avant/après filtrage

L’étape de filtrage permet de nettoyer les données en supprimant les SNPs ou individus de mauvaise qualité. Voici un exemple illustratif des effets d’un tel filtrage sur un jeu de données :

| Critère                              | Avant filtrage | Après filtrage |
|--------------------------------------|----------------|----------------|
| Nombre total de SNPs                | 500,000        | 450,000        |
| Nombre total d’individus            | 1,000          | 980            |
| Taux moyen de données manquantes    | 2.5%           | 0.5%           |
| SNPs avec MAF < 1%                  | 25,000         | 0              |
| SNPs invariants (aucune variation)  | 10,000         | 0              |

Ces résultats sont obtenus à l’aide des commandes précédentes et permettent de garantir que seuls les SNPs informatifs et de qualité sont conservés.

---

## 🧼 Stratégie de nettoyage et pré-traitement

Avant toute analyse, il est recommandé d’appliquer une série de filtres pour garantir la qualité des données :

1. **Filtrage des individus et SNPs avec trop de données manquantes** :

```bash
plink --bfile data --mind 0.1 --geno 0.05 --make-bed --out step1_filtered
```

2. **Filtrage par fréquence allèlique mineure** :

```bash
plink --bfile step1_filtered --maf 0.01 --make-bed --out step2_filtered
```

3. **Élimination des SNPs invariants (optionnel)** :

```bash
plink --bfile step2_filtered --exclude snps_invariants.txt --make-bed --out final_data
```

Ces étapes sont essentielles pour éviter les biais et les faux positifs dans les résultats statistiques.

---

## 🧬 Interprétation des résultats ROH et effet fondateur

L’identification de segments ROH (homozygotie par descendance) est une méthode clé pour repérer des événements d’effet fondateur. Dans les populations issues d’un petit nombre d’ancêtres communs, on observe souvent :

- Une **augmentation du nombre et de la longueur des segments ROH**
- Une **distribution non aléatoire de ces segments**, souvent centrée autour d’une mutation fondatrice
- Des individus porteurs de la mutation avec des ROH partagés dans une même région génomique

*(Insérez ici un graphique illustrant clairement la distribution des segments ROH autour d'une mutation fondatrice.)*

Ces observations peuvent être renforcées par des analyses complémentaires (KING pour la parenté, Gamma pour la datation) dans le pipeline DOCK6.

---

## 🔗 Interprétation des fréquences alléliques (Gamma) et parenté (KING)

### 📈 Gamma – Fréquences alléliques informatives

L’analyse via **Gamma** repose sur l’évolution des fréquences alléliques autour de la mutation fondatrice. Des fréquences intermédiaires (autour de 0.3–0.7) indiquent un **déséquilibre de liaison**, témoin de la présence d’un haplotype ancestral commun.

- Un pic de fréquence allélique cohérent autour d’un locus cible peut suggérer un effet fondateur récent
- La largeur du signal permet d’estimer l’ancienneté de l’événement (plus étroit = plus ancien)

*(Insérez ici un graphique montrant clairement un pic typique de fréquence allélique indiquant un effet fondateur.)*

### 👪 KING – Détection de parenté cachée

L’outil **KING** détecte les liens familiaux au sein d’un groupe d’individus génotypés. En contexte fondateur, on retrouve souvent :

- Une **parenté élevée** entre individus non apparentés sur le papier (IBD élevé)
- Des **groupes familiaux disjoints** partageant les mêmes segments ROH ou allèles rares

*(Insérez ici un exemple visuel, tel qu'une heatmap ou un dendrogramme, illustrant clairement la parenté élevée entre individus supposément non apparentés.)*

Ces éléments renforcent l’hypothèse d’un ancêtre commun récent et viennent compléter les résultats de PLINK.

---

L’identification de segments ROH (homozygotie par descendance) est une méthode clé pour repérer des événements d’effet fondateur. Dans les populations issues d’un petit nombre d’ancêtres communs, on observe souvent :

- Une **augmentation du nombre et de la longueur des segments ROH**
- Une **distribution non aléatoire de ces segments**, souvent centrée autour d’une mutation fondatrice
- Des individus porteurs de la mutation avec des ROH partagés dans une même région génomique

Ces observations peuvent être renforcées par des analyses complémentaires (KING pour la parenté, Gamma pour la datation) dans le pipeline DOCK6.

---

## 🧪 Intégration dans DOCK6

PLINK est exécuté automatiquement dans `run_pipeline.py` selon les étapes
suivantes :

1. **Préparation** : conversion des données (`--make-bed`, `--recodeA`)
2. **ROH** : détection d’autozygotie (`--homozyg`)
3. **Sorties** : fichiers requis pour KING, Gamma et R

---

## ✅ Avantages

- Très rapide même avec de grands jeux de données
- Format de sortie largement standardisé
- Documentation et communauté abondante

## ⚠️ Limites

- Interface CLI uniquement
- Certains formats ou analyses nécessitent PLINK 2.0
- PLINK est sensible aux erreurs de formatage de fichiers ; une vérification préalable rigoureuse est nécessaire
- Les résultats ROH et de fréquences alléliques dépendent fortement du choix des paramètres ; des tests préliminaires sont recommandés

---

### 💡 Astuce

Il est recommandé d’utiliser PLINK dans un environnement conda ou via une image Docker pour assurer la reproductibilité sur différentes machines.

