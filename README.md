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
