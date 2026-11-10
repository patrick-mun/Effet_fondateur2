# D28 — Préspecification fermée de l’évaluateur d’enrichissement 16B

Statut : `FINALIZED_BEFORE_RUN`
Autorisation : demande utilisateur du 10 novembre 2026
Rapport : `CAL-REP-16B`

## 1. Question et portée

D28 mesure la calibration sous le nul et la puissance descriptive de la
primitive réelle d’enrichissement haplotypique de l’étape 16B. Le banc
réutilise les dix scénarios biologiques synthétiques D24, chacun projeté sur
les profils `MAXIMUM`, `BASELINE`, `DENSITY_50` et `MISSING_1PCT`.

Le primaire comporte cinq scénarios positifs et cinq scénarios nuls au profil
`MAXIMUM`. Les autres profils sont des sensibilités appariées et ne créent pas
40 réplications biologiques. Les scénarios ont trois familles indépendantes ;
ce banc ne mesure donc pas l’effet du nombre de familles.

Le calcul n’utilise ni `dock6_real`, ni `data/input`, ni cache interne. Il
n’exécute pas l’orchestrateur complet, bcftools ou la référence 1000 Genomes de
production. Il appelle directement `evaluate_exact_sharing`,
`sample_external_null`, `summarize_null` et la règle de statut de 16B sur les
VCF synthétiques fermés de D24.

## 2. Séparation appels–vérité

La phase d’appel reçoit uniquement les sorties closes de l’étape 13, le VCF
phasé candidat, la carte, les cohortes et les haplotypes porteurs. Elle ne lit
pas `founder_haplotype_case.json`, ne connaît pas le rôle positif/nul et publie
les cas sous le SHA-256 de leur identifiant. Les probabilités, statuts et
non-évaluations sont scellés avant l’ouverture de la vérité D24.

L’évaluation vérifie ensuite les empreintes du manifeste et des appels, ouvre
la vérité et agrège par scénario indépendant. Aucune exclusion postérieure à
la vérité n’est autorisée.

## 3. Statistique et distribution nulle

Pour trois représentants familiaux présélectionnés, chaque bras du partage IBS
exact s’étend depuis la cible jusqu’au premier désaccord ou génotype manquant.
La cible sert d’ancre mais est exclue de la signature. Avec les bornes communes
`L` et `R` et la position cible `T`, la statistique est :

`S_obs = [cM(T) − cM(L)] + [cM(R) − cM(T)]`.

Chaque tirage nul choisit sans remise trois individus du fond synthétique non
apparenté, puis une copie haplotypique par individu. La même primitive calcule
`S_b`. Pour `N` tirages évaluables et `k = Σ 1(S_b ≥ S_obs)`, la probabilité
empirique corrigée est :

`p_emp = (1 + k) / (1 + N)`.

Dix mille tirages évaluables sont demandés par cas, avec au plus cent mille
tentatives. Les tentatives non évaluables restent comptées séparément. Un cas
sans sortie 13 compatible ou sans deux marqueurs partagés sur chaque flanc est
`NOT_EVALUATED` ; il n’est jamais remplacé par une probabilité nulle.

## 4. Seuil du banc et absence de classification de production

Le seuil de production `empirical_classification_threshold` reste
obligatoirement `null`. La sortie auditée de 16B doit donc rester
`NOT_CLASSIFIED` lorsqu’une probabilité est disponible, même si elle est très
faible. `CAL-ENRICH-001` compte toute émission indue de
`UNUSUAL_TARGET_CENTERED_SHARING` ou `NO_UNUSUAL_SHARING_DETECTED`.

Pour calculer uniquement les caractéristiques opératoires du banc, D28 fixe
avant résultats la règle descriptive `appel_banc = 1[p_emp ≤ 0,05]`. Ce seuil
n’est ni un seuil de production, ni un seuil d’acceptation scientifique. Les
seuils d’acceptation de `CAL-ENRICH-002` et `CAL-ENRICH-003` restent `null`.

## 5. Métriques et incertitude

Sur les cinq scénarios nuls indépendants du profil primaire :

`FPR_nul = nombre d’appels_banc nuls / 5`.

Sur les cinq scénarios positifs indépendants :

`puissance = nombre d’appels_banc positifs / 5`.

Une non-évaluation reste au dénominateur et ne compte pas comme appel, ce qui
rend la puissance conservatrice. Le taux conditionnel parmi les seuls cas
évaluables est également publié, mais n’est pas le critère primaire. Les
intervalles descriptifs à 95 % utilisent Wilson, avec le scénario biologique
comme unité indépendante. Les distributions de probabilités, nombres de
tirages, dépassements et motifs de non-évaluation sont publiés par profil.

## 6. Verdicts autorisés

Le cœur peut recevoir `PASS_CORE_EVALUATOR` si les sources et appels sont
empreintés, les probabilités suivent la correction `+1`, les non-évaluations
sont conservées, la cible est exclue et aucune classification de production
n’est émise sans seuil.

`CAL-ENRICH-001` peut réussir son contrat fixe. `CAL-ENRICH-002` et
`CAL-ENRICH-003` restent toujours `REVIEW_REQUIRED` car leurs seuils
`pilot_then_seal` ne sont pas scellés. La qualification globale est donc
`REVIEW_REQUIRED_UNSEALED`, avec `production_promotion=false`,
`calibrated=false` et `founder_effect_proven=false`.

Le faible effectif, les trois familles fixes, le fond nul partagé, l’usage du
phasage-vérité D24 et l’absence de référence externe 1000 Genomes interdisent
une généralisation à l’étude réelle ou à l’étape 16B bout-en-bout.

## 7. Rapport indépendant

Le rapport autonome comporte quatorze sections, au moins trois graphiques
quantitatifs, les équations, numérateurs, dénominateurs, intervalles, taux de
non-évaluation, sensibilités appariées, empreintes et limites. Le rendu ne
recalcule aucune métrique et ne contient aucun identifiant individuel.
