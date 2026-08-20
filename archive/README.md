# Archive — pipeline V1

Ce dossier conserve le pipeline historique (« V1 »), remplacé par la V2
(`src/effet_fondateur/`). Il est conservé pour référence et vérification des
dépendances, mais **ne doit plus être exécuté pour produire un résultat
scientifique** (voir `SESSION.md`, section « Problèmes connus »).

## Contenu

- `run_pipeline.py` : orchestration V1 (PLINK, ROH, IBD, LD, Gamma, adegenet,
  rapport). Utilise des chemins codés en dur
  (`data/input/complex_simulation/...`) et mélange calcul, visualisation et
  interprétation.
- `scripts/` : fonctions appelées par `run_pipeline.py` (prétraitement, ROH,
  IBD, LD, Gamma, adegenet, reporting, utilitaires).
- `interface_effet_fondateur.py` : interface Streamlit pilotant
  `run_pipeline.py`. Incohérence connue : écrit `user_input.ped/map` alors que
  le pipeline lit `genotype_data.ped/map`.
- `Module_WIKI/` : documentation du pipeline V1 (PLINK, adegenet, Gamma, KING),
  consommée uniquement par l'interface Streamlit ci-dessus.

## Avant suppression définitive

Conformément à `SESSION.md` : archiver un script V1 seulement après
cartographie de ses dépendances et validation de son remplacement V2. Ce
dossier peut être supprimé une fois confirmé que chaque fonctionnalité listée
ci-dessus a un équivalent V2 validé (`src/effet_fondateur/stages/` et
`docs/modules/`).
