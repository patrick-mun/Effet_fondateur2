# Pipeline d'analyse de l'effet fondateur — DOCK6

Ce projet analyse un possible effet fondateur autour d'une mutation du gène
`DOCK6` à partir de données SNP au format PLINK. Il combine prétraitement,
détection des ROH, analyses de parenté IBD, déséquilibre de liaison, estimation
de l'âge de la mutation, structure populationnelle et génération de rapports.

## Environnement requis

### Python

- Python 3.10 ou version plus récente (Python 3.12 recommandé sur macOS Intel)
- dépendances déclarées dans `requirements.txt`

Installation recommandée :

```bash
python3.12 -m venv .venv
source .venv/bin/activate
python -m pip install --upgrade pip
python -m pip install -r requirements.txt
```

### Activer l'environnement du projet

L'environnement doit être activé dans chaque nouveau terminal avant de lancer
les scripts Python :

```bash
cd /Users/utilisateur/Documents/python_programme/Effet_fondateur2
source .venv/bin/activate
```

Le terminal affiche alors généralement `(.venv)` au début de la ligne. Pour
vérifier que le bon interpréteur est utilisé :

```bash
which python
python --version
```

Le chemin affiché par `which python` doit se terminer par `.venv/bin/python`.
Une fois l'environnement activé, les commandes peuvent être lancées simplement :

```bash
python run_pipeline.py
streamlit run interface_effet_fondateur.py
pytest test/
```

Pour quitter l'environnement :

```bash
deactivate
```

Le projet importe directement les bibliothèques suivantes :

- `pandas`, `numpy` et `scipy` pour les calculs ;
- `matplotlib` pour les graphiques ;
- `networkx` pour le réseau de parenté ;
- `fpdf2` pour le rapport PDF ;
- `streamlit` pour l'interface ;
- `pytest` et `pytest-html` pour les tests.

### Programmes externes

Les exécutables suivants doivent être accessibles depuis le `PATH` :

- `plink` 1.9 ou compatible : filtrage, ROH, IBD, HWE et LD ;
- `king` : estimation des relations de parenté ;
- `Rscript` : exécution de l'analyse Adegenet.

Le binaire `Gamma` est uniquement nécessaire pour utiliser la fonction
`run_gamma()` ou l'ancien script shell. Le pipeline principal utilise
actuellement l'estimation Python de `scripts/gamma_age_estimation.py`.

Vérification rapide :

```bash
python3 --version
plink --version
king --version
Rscript --version
```

### Packages R

Le script Adegenet nécessite :

```r
install.packages(c("adegenet", "ggplot2", "ape", "ade4", "poppr"))
```

## Données attendues

Le pipeline principal utilise actuellement ces chemins fixes :

```text
data/input/complex_simulation/genotype_data.ped
data/input/complex_simulation/genotype_data.map
data/input/complex_simulation/groupes.txt
data/input/complex_simulation/cas.txt
data/input/complex_simulation/temoins.txt
```

Les fichiers `cas.txt` et `temoins.txt` sont des fichiers PLINK `--keep` à deux
colonnes (`FID IID`). Ils doivent être présents avant l'exécution complète.

### Protocole de conversion ACPA vers PLINK

Ce protocole extrait un chromosome, attribue les rsID et produit les fichiers
PLINK sans modifier les exports ACPA originaux. Les exports peuvent provenir
des puces CytoScan 750K Array et Accel, mais doivent utiliser le build `hg38`.

#### 1. Activer l'environnement

Depuis la racine du projet :

```bash
source .venv/bin/activate
```

#### 2. Placer et vérifier les exports ACPA

Placer un export ChAS `.txt` par individu dans :

```text
data/input/complex_simulation/acpa_samples/
```

Vérifier les fichiers détectés :

```bash
find data/input/complex_simulation/acpa_samples \
  -maxdepth 1 -type f -name '*.txt' -print
```

#### 3. Créer le fichier `samples.tsv`

Créer un modèle à partir des noms de fichiers ACPA :

```bash
python -m simulation_genotype_famille.acpa_to_plink \
  --input-dir data/input/complex_simulation/acpa_samples \
  --create-metadata-template data/input/complex_simulation/samples.tsv
```

Si le fichier existe déjà, le convertisseur refuse de l'écraser. L'examiner
avant d'envisager l'option `--force`.

#### 4. Modifier `samples.tsv`

Ouvrir le tableau, par exemple avec TextEdit :

```bash
open -a TextEdit data/input/complex_simulation/samples.tsv
```

Conserver la ligne d'en-tête et renseigner une ligne par individu à analyser :

| Colonne | Modification attendue |
| --- | --- |
| `FILE` | Conserver exactement le nom du fichier ACPA présent dans `acpa_samples`. |
| `FID` | Indiquer l'identifiant familial. Les apparentés doivent partager le même `FID`. |
| `IID` | Indiquer un identifiant individuel unique, sans espace. |
| `PID` | Indiquer l'`IID` du père s'il est présent, sinon `0`. |
| `MID` | Indiquer l'`IID` de la mère si elle est présente, sinon `0`. |
| `SEX` | Utiliser `1` pour un homme, `2` pour une femme ou `0` si inconnu. |
| `PHENOTYPE` | Utiliser `2` pour atteint, `1` pour témoin, `-9` ou `0` si inconnu. |
| `GROUP` | Utiliser notamment `ATTEINT` ou `TEMOIN`; choisir un nom explicite pour les autres groupes. |

Exemple fictif :

```text
FILE	FID	IID	PID	MID	SEX	PHENOTYPE	GROUP
patient_01.txt	F1	PATIENT_01	0	0	1	2	ATTEINT
patient_02.txt	F1	PATIENT_02	PATIENT_01	0	2	2	ATTEINT
temoin_01.txt	CTRL	TEMOIN_01	0	0	1	1	TEMOIN
```

Supprimer du tableau les lignes correspondant à des fichiers qui ne doivent pas
être analysés. Le fichier donné avec `--rsid-reference-acpa` n'a pas besoin de
figurer dans `samples.tsv` s'il sert uniquement de référence d'annotation. Un
autre exemple fictif est disponible dans
`simulation_genotype_famille/samples.example.tsv`.

#### 5. Lancer l'annotation du chromosome 19

Le fichier `SNP_82404457_(CytoScan750K_Array).txt` contient déjà des rsID. Il
sert de dictionnaire local `Probe Set ID + position hg38 → rsID`; ce n'est pas
le fichier qui détermine les individus analysés.

```bash
python -m simulation_genotype_famille.acpa_to_plink \
  --input-dir data/input/complex_simulation/acpa_samples \
  --metadata data/input/complex_simulation/samples.tsv \
  --output-prefix data/output/acpa_chr19/genotype_data \
  --chromosome 19 \
  --rsid-reference-acpa \
  'data/input/complex_simulation/acpa_samples/SNP_82404457_(CytoScan750K_Array).txt'
```

Le traitement sélectionne d'abord le chromosome demandé, puis attribue les rsID
aux sondes conservées. Il propage en priorité les correspondances fiables de
l'export ACPA annoté donné avec `--rsid-reference-acpa`. Pour les positions
restantes, `bcftools` interroge uniquement les blocs utiles du VCF officiel
NCBI dbSNP build 157 sur GRCh38.p14. Le fichier complet de 28 Go n'est pas
téléchargé. La source exacte utilisée est enregistrée dans le rapport JSON.

À une position comportant plusieurs entrées dbSNP, seuls les SNV dont les
allèles `REF/ALT` sont compatibles avec les `Forward Strand Base Calls` sont
considérés. Une correspondance ambiguë ou absente conserve temporairement le
`Probe Set ID` et est signalée dans les rapports. Ajouter `--require-rsid` pour
exclure ces sondes au lieu de conserver leur identifiant de sonde.

Le convertisseur utilise l'intersection des sondes entre les échantillons et
refuse d'écraser une sortie existante sans `--force`. Il ne modifie jamais les
exports ACPA originaux : les rsID sont écrits dans les sorties PLINK et le
fichier d'audit.

Pour conserver une analyse précédente, modifier le dossier indiqué dans
`--output-prefix`, par exemple `data/output/acpa_chr19_analyse2/genotype_data`.
N'utiliser `--force` que pour remplacer volontairement tous les résultats du
dossier choisi.

#### 6. Examiner les résultats d'annotation

Afficher le rapport :

```bash
cat data/output/acpa_chr19/acpa_conversion_report.json
```

Vérifier en priorité `sample_count`, `output_marker_count`,
`rsid_marker_count`, `unresolved_rsid_count` et `excluded_marker_count`.

Ouvrir le tableau détaillé :

```bash
open data/output/acpa_chr19/acpa_marker_audit.tsv
```

Les fichiers générés sont :

```text
data/output/acpa_chr19/genotype_data.ped
data/output/acpa_chr19/genotype_data.map
data/output/acpa_chr19/groupes.txt
data/output/acpa_chr19/cas.txt
data/output/acpa_chr19/temoins.txt
data/output/acpa_chr19/acpa_conversion_report.json
data/output/acpa_chr19/acpa_excluded_markers.tsv
data/output/acpa_chr19/acpa_marker_audit.tsv
```

Le rapport JSON résume les échantillons, appels manquants, marqueurs conservés,
rsID résolus et exclusions. Le TSV d'audit conserve la correspondance entre
sonde, rsID, position hg38, allèles observés et statut d'annotation.

#### 7. Valider les fichiers avec PLINK

```bash
plink \
  --file data/output/acpa_chr19/genotype_data \
  --make-bed \
  --out data/output/acpa_chr19/genotype_data_checked
```

PLINK doit terminer sans erreur et indiquer le nombre de variants, d'individus
et le taux global de génotypage.

#### 8. Options complémentaires

Les rsID sont utilisés par défaut dans le fichier MAP. Pour conserver l'union
des sondes, désactiver l'interrogation NCBI ou fournir un VCF dbSNP local
indexé :

```bash
--marker-mode union
--no-dbsnp-annotation
--dbsnp-vcf /chemin/vers/dbsnp_grch38.vcf.gz
```

Les sondes absentes d'un échantillon sont codées `0 0` en mode `union`.

## Exécution

### Pipeline complet

Depuis la racine du projet :

```bash
.venv/bin/python run_pipeline.py
```

Le script exécute successivement :

1. le filtrage et le contrôle qualité PLINK ;
2. la détection et la visualisation des ROH ;
3. les analyses IBD avec KING et PLINK ;
4. le calcul du LD global et autour de la mutation ;
5. la préparation des fréquences et l'estimation de l'âge de la mutation ;
6. l'analyse DAPC avec R/Adegenet ;
7. la génération des rapports HTML et PDF.

Les résultats sont écrits dans `data/output/complex_simulation/`.

Le script ne possède pas encore d'options pour exécuter une seule étape.

### Interface Streamlit

```bash
.venv/bin/streamlit run interface_effet_fondateur.py
```

L'interface permet de consulter le wiki et les résultats et de lancer le
pipeline. Le raccordement des fichiers uploadés au pipeline doit encore être
corrigé : le pipeline continue actuellement à lire `genotype_data.ped/map`.

### Tests

```bash
pytest test/
```

Pour produire le rapport HTML et l'ouvrir sur macOS :

```bash
./run_all_tests.sh
```

## Organisation

```text
run_pipeline.py                 orchestration principale
interface_effet_fondateur.py    interface Streamlit
scripts/                        modules d'analyse
simulation_genotype_famille/    préparation ACPA et simulation PED/MAP
Module_WIKI/                    documentation scientifique des outils
data/input/                     données sources et métadonnées
data/output/                    résultats générés
test/                           tests automatisés
```

La description détaillée de PLINK, KING, Gamma et Adegenet se trouve dans
`Module_WIKI/`.

## État du projet

Le pipeline correspond à une version alpha de recherche. Les résultats doivent
être contrôlés avant toute interprétation biologique ou médicale. Plusieurs
points de raccordement et de validation restent à fiabiliser avant une nouvelle
analyse complète.

Auteur principal : Patrick MUNIER — Laboratoire de génétique, projet DOCK6.
