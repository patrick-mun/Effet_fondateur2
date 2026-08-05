# Pipeline d'analyse de l'effet fondateur — DOCK6

Ce projet analyse un possible effet fondateur autour d'une mutation du gène
`DOCK6` à partir de données SNP au format PLINK. Il combine prétraitement,
détection des ROH, analyses de parenté IBD, déséquilibre de liaison, estimation
de l'âge de la mutation, structure populationnelle et génération de rapports.

## Environnement requis

### Python

- Python 3.12
- dépendances déclarées dans `requirements.txt`

Installation recommandée :

```bash
python3.12 -m venv .venv
source .venv/bin/activate
python -m pip install --upgrade pip
python -m pip install -r requirements.txt
python -m pip install -e .
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
- `PyYAML` et `jsonschema` pour la configuration V2 ;
- `matplotlib` pour les graphiques ;
- `networkx` pour le réseau de parenté ;
- `fpdf2` pour le rapport PDF ;
- `streamlit` pour l'interface ;
- `pytest` et `pytest-html` pour les tests.

### Configuration V2

La V2 est développée parallèlement au pipeline historique. Une configuration
peut être validée sans lire les données ni lancer d'analyse :

```bash
effet-fondateur validate-config config/pipeline.example.yaml
effet-fondateur validate-config config/studies/dock6.example.yaml
```

Le profil DOCK6 conserve volontairement les coordonnées et allèles non confirmés
à `null`. Leur absence bloquera les futures étapes scientifiques concernées au
lieu de provoquer l'utilisation de valeurs supposées.

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

### Simuler des témoins non apparentés

Le script `simulate_unrelated_controls.py` remplace l'ancienne simulation à
fréquence fixe. Il estime, pour chaque SNP, les fréquences alléliques et le taux
de données manquantes dans le PED source, puis tire indépendamment les deux
allèles de chaque témoin sous Hardy–Weinberg.

Chaque témoin simulé reçoit un `FID` et un `IID` uniques, des parents `0/0`, le
phénotype témoin `1` et le groupe `TEMOIN_SIMULE`. La graine aléatoire et toutes
les hypothèses sont enregistrées dans un rapport JSON. Les fichiers source ne
sont jamais modifiés.

#### Simulation interactive

Cette commande demande le nombre de témoins à créer :

```bash
python -m simulation_genotype_famille.simulate_unrelated_controls \
  --source-prefix data/output/acpa_chr19/genotype_data \
  --output-prefix \
    data/output/acpa_chr19_avec_temoins_simules/genotype_data \
  --seed 20260804
```

Pour fournir directement le nombre, par exemple 50 témoins :

```bash
python -m simulation_genotype_famille.simulate_unrelated_controls \
  --source-prefix data/output/acpa_chr19/genotype_data \
  --output-prefix \
    data/output/acpa_chr19_avec_temoins_simules/genotype_data \
  --n-controls 50 \
  --seed 20260804
```

Par défaut, le PED produit contient les individus ACPA source suivis des témoins
simulés. Ajouter `--controls-only` pour produire uniquement les témoins. Ajouter
`--no-simulate-missingness` pour ne pas reproduire les taux de données
manquantes observés. Comme pour le convertisseur ACPA, `--force` est requis pour
remplacer volontairement des sorties existantes.

Sorties principales :

```text
data/output/acpa_chr19_avec_temoins_simules/genotype_data.ped
data/output/acpa_chr19_avec_temoins_simules/genotype_data.map
data/output/acpa_chr19_avec_temoins_simules/groupes.txt
data/output/acpa_chr19_avec_temoins_simules/cas.txt
data/output/acpa_chr19_avec_temoins_simules/temoins.txt
data/output/acpa_chr19_avec_temoins_simules/simulation_marker_frequencies.tsv
data/output/acpa_chr19_avec_temoins_simules/control_simulation_report.json
```

Valider le jeu combiné avec PLINK :

```bash
plink \
  --file data/output/acpa_chr19_avec_temoins_simules/genotype_data \
  --make-bed \
  --out data/output/acpa_chr19_avec_temoins_simules/genotype_data_checked
```

> **Limite scientifique importante :** ces témoins sont indépendants par
> construction, mais les SNP sont simulés indépendamment. Ils ne reproduisent
> donc ni le déséquilibre de liaison, ni les haplotypes, ni les ROH d'une
> population réelle. Les fréquences sont en outre estimées sur la petite cohorte
> familiale source, enrichie en individus atteints. Utiliser ces témoins pour
> tester le fonctionnement du pipeline ou pour une exploration clairement
> étiquetée, jamais comme preuve biologique d'un effet fondateur. Avec un nombre
> limité de SNP, certaines paires indépendantes peuvent aussi dépasser
> légèrement un seuil KING faible par fluctuation aléatoire.

### Injecter et documenter la mutation

L'injection se fait **après** la conversion ACPA et, pour les essais actuels,
après l'ajout des témoins simulés. Le module `inject_mutation.py` crée une copie
du PED/MAP, insère un marqueur `rsMUT...` à sa position physique et produit les
métadonnées destinées au rapport. Il ne modifie jamais les fichiers source.

Les génotypes de mutation ne sont pas déduits du statut clinique. Chaque
individu doit recevoir explicitement `REF/REF`, `REF/ALT`, `ALT/ALT` ou
`MISSING`, par règle de groupe ou par fichier individuel. Cette séparation
évite de transformer une hypothèse clinique en donnée génétique silencieuse.

#### 1. Créer la fiche de mutation

Copier le modèle dans les données de travail :

```bash
cp simulation_genotype_famille/mutation.example.json \
  data/input/complex_simulation/mutation_info.json
```

Champs requis pour l'injection :

- `gene` : symbole officiel du gène, par exemple `DOCK6` ;
- `mutation_name` : nom lisible du variant ;
- `chromosome` et `position_bp` : coordonnées génomiques ;
- `reference_allele` et `alternate_allele` : bases du brin de référence ;
- `assembly` : assemblage, `GRCh38` par défaut ;
- `marker_id` : facultatif, sinon `rsMUT_<position>` est créé.

Les champs `hgvs_c`, `hgvs_p`, `transcript`, `disease`, `inheritance`,
`variant_type`, `dbsnp_id`, `clinvar_id`, `laboratory_method`, `source` et
`notes` ne modifient pas l'analyse PLINK, mais documentent le rapport. Pour un
variant d'épissage, `hgvs_p` peut rester vide si l'effet protéique n'est pas
établi.

> **Contrôle indispensable :** `c.1833-1G>T` décrit le variant par rapport à un
> transcrit. Les allèles écrits dans le PED doivent, eux, correspondre au brin
> génomique de référence de l'assemblage choisi. Il faut donc confirmer la
> position, le transcrit et les allèles avant l'injection.

#### 2. Définir les génotypes

Pour une règle commune à un groupe, répéter `--group-genotype`. Exemple pour
une hypothèse récessive à confirmer avec les résultats moléculaires :

```text
ATTEINT=ALT/ALT
HTZ=REF/ALT
SAINS=REF/REF
TEMOIN_SIMULE=REF/REF
```

Pour corriger ou définir un individu séparément, copier puis compléter le
modèle TSV :

```bash
cp simulation_genotype_famille/mutation_genotypes.example.tsv \
  data/input/complex_simulation/mutation_genotypes.tsv
```

Le fichier comporte les colonnes `FID`, `IID` et `GENOTYPE`. Une affectation
individuelle est prioritaire sur la règle de groupe correspondante.

#### 3. Lancer l'injection

Avec le jeu comprenant les 50 témoins simulés :

```bash
python -m simulation_genotype_famille.inject_mutation \
  --source-prefix \
    data/output/acpa_chr19_avec_temoins_simules/genotype_data \
  --output-prefix data/output/acpa_chr19_mutation/genotype_data \
  --mutation-config data/input/complex_simulation/mutation_info.json \
  --groups \
    data/output/acpa_chr19_avec_temoins_simules/groupes.txt \
  --group-genotype ATTEINT=ALT/ALT \
  --group-genotype HTZ=REF/ALT \
  --group-genotype SAINS=REF/REF \
  --group-genotype TEMOIN_SIMULE=REF/REF
```

Ajouter si nécessaire :

```bash
--individual-genotypes data/input/complex_simulation/mutation_genotypes.tsv
```

L'injecteur refuse un individu sans génotype explicite, un MAP non trié, un
chromosome discordant et toute sortie existante. Utiliser `--force` uniquement
pour remplacer volontairement une injection précédente.

#### 4. Contrôler les sorties

```text
data/output/acpa_chr19_mutation/genotype_data.ped
data/output/acpa_chr19_mutation/genotype_data.map
data/output/acpa_chr19_mutation/mutation_info.json
data/output/acpa_chr19_mutation/project_info.json
data/output/acpa_chr19_mutation/mutation_genotype_audit.tsv
data/output/acpa_chr19_mutation/mutation_injection_report.json
```

Le rapport JSON contient les empreintes SHA-256 des PED/MAP source et produits,
les effectifs par génotype et l'index d'insertion. `project_info.json` reprend
les libellés affichables par le rapport final. Vérifier ensuite :

```bash
plink \
  --file data/output/acpa_chr19_mutation/genotype_data \
  --make-bed \
  --out data/output/acpa_chr19_mutation/genotype_data_checked

plink \
  --file data/output/acpa_chr19_mutation/genotype_data \
  --mendel \
  --out data/output/acpa_chr19_mutation/mendel_check
```

Toute erreur mendélienne au marqueur `rsMUT...` doit être résolue ou documentée
avant le pipeline. Pour l'analyse actuelle, copier ensuite de façon contrôlée
le PED, le MAP, `groupes.txt`, `cas.txt`, `temoins.txt` et `project_info.json`
dans `data/input/complex_simulation/`, car `run_pipeline.py` lit encore ce
dossier fixe.

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
.venv/bin/python -m pytest test tests
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
simulation_genotype_famille/    préparation ACPA, simulation et injection
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
