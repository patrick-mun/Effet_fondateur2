# LD local secondaire — contrat scientifique V2

## Question et portée

L'étape `15_analyze_local_ld` décrit la corrélation locale entre variants dans
des cohortes indépendantes déjà figées. Elle n'identifie pas l'apparentement,
ne démontre pas l'IBD, ne définit pas l'haplotype fondateur et ne date pas le
variant.

La méthode versionnée est `plink19_local_ld_secondary_v1`. PLINK 1.9 est appelé
deux fois afin de ne pas confondre deux définitions :

- `--r2` : carré de la corrélation des dosages alléliques, mesure primaire ;
- `--r2 dprime` : `D′` absolu estimé par maximum de vraisemblance haplotypique,
  avec le `r²` haplotypique correspondant conservé séparément.

Python ne recalcule aucune statistique de LD. Un export dosage temporaire PLINK
sert seulement au décompte des individus appelés pour chaque paire.

## Cohortes

`controls_unrelated` est la cohorte descriptive principale. Les cohortes
`target_carriers_independent` et `family_noncarriers` sont secondaires et
restent séparées. Elles sont lues depuis `cohorts.frozen.tsv` de l'étape `09` ;
aucune sélection individuelle n'est reconstruite à l'étape `15`.

Moins de 5 individus conduit à `NOT_EVALUATED`. Entre 5 et 19 individus, les
valeurs sont `EXPLORATORY_SMALL_N`. À partir de 20 individus, elles peuvent
recevoir le statut `DESCRIPTIVE_PRIMARY`. Ces seuils configurables sont des
gardes de projet contre des estimations manifestement instables, pas une
garantie de précision ni une frontière biologique.

Une paire exige par défaut cinq individus appelés et une MAF d'au moins 0,05
pour chacun de ses variants. Cette garde configurable vise la stabilité
descriptive ; elle ne constitue ni une règle biologique ni un filtre causal.

## Univers de variants et comparabilité

L'univers commun est l'intersection du jeu chromosome cible final de l'étape
`10` et de la région/carte validée de l'étape `11`. Les variants doivent être
autosomiques, bialléliques, identifiés sans ambiguïté et conserver le même ordre
REF/ALT. Les mêmes limites bp/cM et le même univers sont appliqués à toutes les
cohortes évaluées. Un variant monomorphe dans une cohorte reste visible dans
l'audit mais ses paires ne sont pas interprétées.

Le variant cible est conservé. Toute paire qui l'implique porte
`INVOLVES_TARGET=true`, car la séparation préalable entre porteurs et
non-porteurs conditionne directement sa fréquence et peut rendre la comparaison
circulaire.

## Sorties et limites

`local_ld_pairs.tsv` conserve les deux `r²`, `D′`, les distances bp/cM,
l'effectif appelé et un statut contrôlé. `local_ld_summary.tsv` agrège par
cohorte et classe de distance uniquement lorsqu'au moins le nombre configuré de
paires est évaluable. Les sorties PLINK natives sont conservées pour audit.

Il n'y a ni blocs LD, ni correction de tests multiples, ni p-value par paire,
ni comparaison inférentielle naïve entre matrices. Les différences visibles
peuvent refléter la taille d'échantillon, les fréquences alléliques, la structure
populationnelle, la sélection sur le génotype cible, la densité de la puce ou
la résolution de la carte.

## Articulation avec les étapes 13 et 14

L'étape `15` ne consomme aucune sortie des étapes `13–14` et ne leur fournit
aucune entrée. Le partage IBS centré cible reste défini à l'étape `13`; les bras
en cM reçus par Gamma restent définis à l'étape `14`. Le LD local peut seulement
être présenté comme contexte descriptif secondaire lors de la synthèse finale.
