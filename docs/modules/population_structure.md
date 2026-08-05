# Étape 08 — Structure populationnelle

## Responsabilité

`08_analyze_population_structure` estime les axes de variation genome-wide sur
une référence aussi indépendante que possible, puis projette tous les individus
du panel sans laisser les apparentés déterminer les axes.

L'étape est exploratoire. Elle ne transforme aucune proposition d'outlier en
exclusion effective et n'utilise ni statut clinique, ni groupe descriptif pour
ajuster la PCA.

## Entrées

- `kinship_panel.bed/.bim/.fam` et son descripteur de l'étape `06` ;
- `independent_set_proposal.tsv` de l'étape `07` ;
- `samples.master.tsv` pour les identifiants, familles, groupes et lots.

Les empreintes, les effectifs et les ensembles d'individus sont contrôlés avant
le calcul. La référence contient exactement les individus dont
`PROPOSED_INCLUDED` vaut `true`. Cette proposition peut rester en attente de
validation humaine ; cette limite est inscrite dans l'audit et devra être
résolue avant le gel des cohortes.

## Modèle PCA

PLINK exporte temporairement le nombre d'allèles A1 sous forme de dosages
`0`, `1` ou `2`. Pour chaque variant, la fréquence `p` est estimée uniquement
sur la référence indépendante. Le dosage est standardisé par :

```text
(dosage - 2p) / sqrt(2p(1-p))
```

Un génotype manquant est imputé à `2p`, donc à zéro après standardisation. Un
variant est exclu du modèle s'il est monomorphe dans la référence ou si son taux
d'appel y est inférieur à `min_reference_call_rate`.

La matrice est divisée par la racine du nombre de variants informatifs puis
décomposée par SVD. Le signe de chaque axe est rendu déterministe en imposant un
loading positif au variant de poids absolu maximal. Les apparentés sont ensuite
projetés avec les mêmes fréquences, écarts-types et loadings ; modifier leurs
dosages ne peut donc pas modifier les axes ni les scores de référence.

Le nombre réellement calculé est limité par `requested_components`, le rang de
la référence et dix colonnes contractuelles `PC1` à `PC10`. Les colonnes au-delà
du rang disponible restent vides.

## Outliers

Sur les premiers axes configurés, l'étape calcule :

```text
distance² = somme(PC_k² / valeur_propre_k)
```

La p-value associée utilise une approximation du chi-deux. Elle sert uniquement
à signaler des individus pour revue. Le seuil nominal `outlier_alpha` n'est pas
une preuve d'ascendance, une règle clinique ou un filtre automatique. Toute
proposition ajoute `population_structure_exclusion_approval` au manifest.

Aucune référence populationnelle externe n'est utilisée dans cette première
version ; les axes décrivent uniquement les individus de l'étude. La DAPC reste
hors périmètre tant qu'une question supervisée et une validation du
surapprentissage ne sont pas définies.

## Sorties

- `population_scores.tsv` : coordonnées, statut référence/projeté et contexte ;
- `population_eigenvalues.tsv` : valeurs propres et variance expliquée ;
- `population_outliers.tsv` : distance, p-value et proposition révisable ;
- `population_variant_loadings.tsv` : fréquences, échelles et loadings du modèle ;
- `population_structure_report.json` : résumé non individuel ;
- `stage_outputs.json`, `audit.json` et `checksums.sha256`.

Les scores, outliers et loadings sont classés `sensitive_genetic`. Le fichier
PLINK `.raw` temporaire n'est jamais publié.

## Codes de retour

- `0` : calcul et contrats valides, même si une revue manuelle est requise ;
- `2` : configuration ou cohérence des entrées invalide ;
- `3` : échec PLINK ou dosage exporté incohérent ;
- `4` : référence, nombre de variants informatifs ou rang PCA insuffisant.
