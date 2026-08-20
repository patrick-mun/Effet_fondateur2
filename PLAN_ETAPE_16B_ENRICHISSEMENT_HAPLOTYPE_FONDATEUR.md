# Plan de développement — étape 16B, enrichissement du partage haplotypique fondateur

## 1. Objet et décision d'architecture

L'étape proposée est nommée :

`16B_evaluate_founder_haplotype_enrichment`

Elle est insérée après `16A_analyze_reference_ancestry` et avant
`17_run_sensitivity_analyses`, sans renuméroter les étapes `17–19`.

L'étape `15_analyze_local_ld` ne doit pas être transformée. Elle répond déjà à
une question distincte et valide : décrire le déséquilibre de liaison de fond
dans la région, principalement chez les témoins indépendants. L'étape 16B doit
répondre à une autre question :

> Le segment haplotypique partagé autour de la variation cible par des familles
> porteuses indépendantes est-il plus long ou plus rare que le partage attendu
> fortuitement au même locus dans les haplotypes de fond ?

Cette séparation évite de présenter comme un LD populationnel une statistique
conditionnée par le recrutement des porteurs, leur génotype et leur structure
familiale.

## 2. Place dans le pipeline V2

Le graphe cible devient :

```text
12 phasage -> 13 haplotype fondateur/IBS -> 14 datation
09 + 10 + 11 -> 15 LD local de fond
09 + 10 + 11 -> 16 ROH
06 + 12 + 16 -> 16A positionnement 1000G
02 + 09 + 11 + 12 + 13 + 16A -> 16B rareté du partage haplotypique
13 + 14 + 15 + 16 + 16A + 16B -> 17 sensibilités
17 -> 18 visualisations -> 19 rapport et revue finale
```

### 2.1 Dépendances directes recommandées

| Producteur | Artefacts requis | Usage dans 16B |
|---|---|---|
| `02` | `samples_master` | identités maître et familles, sans afficher les identifiants |
| `09` | `cohorts_frozen` | unités indépendantes gelées et rôles analytiques |
| `11` | `target_genetic_map` | positions bp/cM et ordre des marqueurs |
| `12` | `shapeit5_final_bcf`, index, `carrier_haplotypes`, référence harmonisée et son index | haplotypes de l'étude, copies mutantes, haplotypes de fond et référence locale |
| `13` | `founder_segments`, `founder_consensus`, matrice pairwise et résumé | définition primaire déjà figée du candidat IBS et sélection des unités |
| `16A` | résumé, audit des variants et provenance de référence | confirmation des 2 504 individus/5 008 haplotypes et du même univers harmonisé |

L'étape 15 n'est pas une entrée de calcul de 16B. Son résumé peut être présenté
à côté du résultat 16B par l'étape 18, mais il ne doit ni définir les bornes du
segment, ni sélectionner les marqueurs, ni modifier la statistique primaire.

L'étape 16 n'est pas une entrée de calcul non plus. Le ROH reste une mesure
d'autozygotie individuelle séparée du partage inter-familial.

### 2.2 Criticité

16B est secondaire et exploratoire par défaut. Une non-évaluation ne bloque ni
la datation de 14, ni les résultats déjà publiés de 13, 15, 16 ou 16A. Si 16B
est activée, son intégrité devient toutefois obligatoire pour ses figures, ses
sensibilités et sa section du rapport.

## 3. Hypothèse et unité d'analyse

### 3.1 Hypothèse primaire

Sous une origine ancestrale commune récente, les chromosomes portant la même
variation rare peuvent conserver autour de la cible un haplotype partagé plus
long que celui obtenu par hasard entre chromosomes de fond au même locus.

Le résultat attendu n'est pas une preuve automatique d'identité par descendance
ni d'effet fondateur. Il constitue un argument empirique supplémentaire,
conditionnel au phasage, à la carte, à la densité ACPA, à la sélection des
familles et à la représentativité des références.

### 3.2 Unités indépendantes

L'unité primaire est la famille indépendante gelée, jamais l'individu ou la
copie chromosomique brute. Dans le premier profil DOCK6 actuellement observé :

- 9 porteurs explicites ;
- 5 génotypes `A/A` et 4 génotypes `C/A` ;
- 14 copies chromosomiques portant `A` ;
- mais seulement 3 familles indépendantes.

L'analyse principale utilise donc trois unités, pas neuf individus ni quatorze
observations indépendantes.

Les copies mutantes supplémentaires servent uniquement à :

- vérifier la ségrégation et la cohérence intrafamiliales ;
- préciser les événements de recombinaison ;
- définir le segment familial minimal ;
- détecter plusieurs arrière-plans haplotypiques au sein d'une famille ;
- alimenter des sensibilités explicitement étiquetées.

Elles ne réduisent pas artificiellement l'incertitude inter-familiale.

## 4. Statistique primaire préspécifiée

### 4.1 Définition du segment observé

La règle d'extension doit être exactement celle de la méthode primaire de 13,
`target_centered_exact_ibs_v1`, sans redéfinition après lecture des résultats :

1. partir de la variation cible configurée ;
2. exclure la cible de la signature comparée au fond afin d'éviter une
   concordance triviale déterminée par le génotype ;
3. avancer séparément à gauche et à droite ;
4. arrêter un bras au premier marqueur manquant ou discordant selon le contrat
   strict de 13 ;
5. exprimer les deux bras en nombre de marqueurs, bp et cM ;
6. conserver les limites physiques et génétiques exactes.

La sélection du chromosome mutant représentatif de chaque famille doit être
déterministe et décidée avant le test. Le choix recommandé est de réutiliser
exactement la copie porteuse déjà retenue par 13. Il est interdit de choisir a
posteriori, dans chaque famille, la copie donnant le segment le plus long.

### 4.2 Statistique principale

La statistique principale proposée est :

```text
T_TOTAL_CM = LEFT_SHARED_CM + RIGHT_SHARED_CM
```

où les deux longueurs sont celles du segment exact commun aux unités familiales
indépendantes.

Les statistiques secondaires, publiées sans multiplier les conclusions, sont :

- `LEFT_SHARED_CM` et `RIGHT_SHARED_CM` séparément ;
- `MIN_ARM_CM = min(LEFT_SHARED_CM, RIGHT_SHARED_CM)` ;
- longueur totale en bp ;
- nombre de marqueurs informatifs partagés ;
- fréquence exacte de la signature observée dans le fond interne et externe.

Une autre statistique primaire ne pourra être adoptée qu'avant implémentation
et avec justification documentée. Il ne faut pas sélectionner après coup celle
qui donne la plus petite probabilité empirique.

## 5. Distributions nulles

### 5.1 Fond interne de l'étude

Le fond interne est construit à partir des haplotypes non porteurs, fiables et
appelés dans le BCF final de 12. Les chromosomes portant la variation cible et
les unités familiales utilisées pour l'observation sont exclus du fond.

Chaque tirage nul doit :

1. contenir le même nombre d'unités que l'observation primaire ;
2. sélectionner des haplotypes appartenant à des individus distincts ;
3. ne jamais prendre `H1` et `H2` du même individu dans un même tirage ;
4. utiliser le même locus, les mêmes marqueurs et la même carte ;
5. appliquer exactement la même règle d'arrêt que pour l'observation ;
6. conserver les non-évaluations dues aux données manquantes, sans les convertir
   en longueurs nulles.

Avec environ 120 haplotypes de fond, l'énumération exhaustive des triplets est
préférable si elle reste compatible avec les contraintes d'individus distincts.
Elle évite une erreur Monte-Carlo inutile et produit le dénominateur exact du
fond interne disponible.

### 5.2 Fond externe 1000 Genomes

Le fond externe utilise les 5 008 haplotypes locaux des 2 504 individus 1000G
non apparentés confirmés par 16A. Les haplotypes doivent provenir de la
référence harmonisée de 12 et être limités au même univers de variants validé
par 16A.

Chaque tirage doit sélectionner des individus distincts. Pour trois unités,
l'espace complet est trop grand ; un Monte-Carlo reproductible est donc admis
avec :

- graine obligatoire et enregistrée ;
- nombre de tirages configuré avant exécution ;
- minimum recommandé de 100 000 tirages évaluables ;
- poursuite jusqu'au nombre évaluable demandé, avec plafond d'essais audité ;
- intervalle binomial Monte-Carlo autour de la probabilité estimée ;
- aucun arrêt anticipé fondé sur la valeur observée.

Le primaire externe utilise toutes les populations de référence. Les analyses
par superpopulation sont des sensibilités descriptives, car 1000G ne représente
pas exhaustivement les histoires réunionnaise, malgache ou comorienne. Elles ne
doivent pas dépendre d'une attribution automatique d'ascendance aux individus
de l'étude.

### 5.3 Probabilité empirique

Pour un Monte-Carlo de `N` tirages évaluables :

```text
p_empirical = (1 + count(T_null >= T_observed)) / (1 + N)
```

Pour une énumération exhaustive, l'étape publie à la fois la fraction exacte et
une version conservatrice avec correction `+1`, clairement distinguées.

Une valeur empirique n'est calculée que si :

- le segment observé satisfait les gardes de 13 ;
- trois unités indépendantes au minimum sont évaluables ;
- le fond contient assez de tirages valides ;
- les marqueurs et règles sont identiques entre observation et nul.

Le seuil éventuel de classification doit être préspécifié. À défaut, la valeur
reste numérique avec le statut `NOT_CLASSIFIED`, conformément à la philosophie
de 17.

## 6. Analyses secondaires et sensibilités internes

Les analyses suivantes sont utiles mais ne remplacent pas le primaire :

1. **Tous les chromosomes mutants explicites** : décrire les 14 copies portant
   `A`, regroupées par famille, et vérifier leur compatibilité avec le consensus.
2. **Leave-one-family-out** : répéter la description sur les paires de familles ;
   avec deux familles, aucune probabilité primaire ne doit être annoncée.
3. **Une autre copie intrafamiliale préspécifiée** : vérifier que le résultat ne
   dépend pas du représentant choisi, sans retenir le maximum.
4. **Masquage des marqueurs peu appelés** : seuils définis avant exécution.
5. **Variantes des règles de terminaison** : exact IBS strict, tolérance limitée
   ou exclusion de marqueurs incertains, chacune versionnée.
6. **Sous-populations 1000G** : uniquement comme sensibilité externe.
7. **Fenêtres physiques/génétiques** : la fenêtre ne peut être élargie au-delà
   de la région préparée et cartographiée sans produire un nouveau run.

Les scénarios modifiant la carte, le phasage, les cohortes ou les seuils de QC
restent des runs V2 distincts consolidés par 17. Ils ne sont pas recalculés
silencieusement dans le dossier primaire de 16B.

## 7. Interprétation scientifique autorisée

### 7.1 Formulations possibles

- `NOT_EVALUATED` : effectif, phase, marqueurs ou fond insuffisants.
- `NO_UNUSUAL_SHARING_DETECTED` : le partage observé n'est pas inhabituel dans
  le fond selon la règle préspécifiée.
- `UNUSUAL_TARGET_CENTERED_SHARING` : partage plus rare que le seuil
  préspécifié, compatible avec un haplotype ancestral commun.
- `MULTIPLE_CARRIER_BACKGROUNDS` : chromosomes mutants incompatibles avec une
  origine haplotypique unique dans la résolution disponible.
- `NOT_CLASSIFIED` : métriques publiées sans seuil préalable.

Même en cas de partage inhabituel, le langage obligatoire reste :

> résultat compatible avec un haplotype ancestral commun et apportant un
> argument supplémentaire à l'hypothèse d'effet fondateur.

### 7.2 Conclusions interdites

16B ne peut pas, isolément :

- prouver un effet fondateur ou un ancêtre unique ;
- transformer IBS en IBD sans méthode IBD validée ;
- estimer l'âge de la mutation en dehors de 14 ;
- inférer une origine ethnique ou généalogique ;
- annoncer une association causale ;
- additionner les résultats de 13, 15, 16, 16A et 16B dans un score composite ;
- traiter les membres apparentés comme observations indépendantes ;
- interpréter l'absence de LD classique chez les témoins comme preuve positive.

## 8. Articulation avec le LD général de 15

La comparaison souhaitée doit être formulée comme deux observations parallèles :

1. **Étape 15** : décroissance du LD de fond entre variants polymorphes chez les
   témoins indépendants, avec `r²` et `D′` séparés.
2. **Étape 16B** : rareté de la longueur d'un haplotype exact partagé autour de
   la cible entre familles porteuses indépendantes.

Ces mesures n'ont ni la même unité, ni la même distribution nulle. Elles ne
doivent donc pas être soustraites, divisées ou soumises à un test naïf de
différence de moyennes.

La présentation conjointe autorisée est :

> Dans une région où le LD de fond décroît rapidement, les familles porteuses
> partagent autour de la cible un segment dont la longueur est située à tel
> quantile de la distribution de partage fortuit au même locus.

## 9. Contrats de configuration

Ajouter au schéma de configuration :

```yaml
stages:
  evaluate_founder_haplotype_enrichment:
    enabled: false
    parameters:
      method: target_centered_empirical_haplotype_sharing_v1
      primary_statistic: total_shared_cm
      minimum_independent_units: 3
      minimum_flank_markers: 2
      internal_null_mode: exhaustive
      external_null_draws: 100000
      external_null_max_attempts: 1000000
      random_seed: 161602026
      minimum_evaluable_null_draws: 10000
      empirical_classification_threshold: null
      run_superpopulation_sensitivities: true
      bcftools_timeout_seconds: 300
```

Règles de validation :

- types stricts et bornes positives ;
- graine obligatoire si Monte-Carlo ;
- seuil facultatif mais fixé avant le run ;
- cohérence entre nombre de tirages et minimum évaluable ;
- méthode et statistique appartenant à des énumérations versionnées ;
- aucun gène, chromosome, allèle ou emplacement codé en dur.

## 10. Artefacts de sortie proposés

### 10.1 Tables sensibles

- `founder_haplotype_units.tsv` : unités familiales et copies utilisées,
  pseudonymisables mais classées `sensitive_genetic` ;
- `founder_haplotype_boundaries.tsv` : limites observées par bras et unité ;
- `founder_haplotype_family_consistency.tsv` : cohérence intrafamiliale des
  copies mutantes ;
- `founder_haplotype_null_draws.tsv.gz` : tirages ou histogramme détaillé selon
  la taille, sans identifiant d'étude ;
- `founder_haplotype_variant_audit.tsv` : marqueurs inclus, rejetés et raisons.

### 10.2 Résumés partageables en interne

- `founder_haplotype_enrichment_summary.tsv` : observation, fonds, quantiles,
  probabilités et intervalles ;
- `founder_haplotype_enrichment_summary.json` : méthode, statuts, effectifs,
  statistique primaire, limites et interdictions ;
- `audit.json`, `stage_inputs.json`, `stage_outputs.json` et
  `checksums.sha256` conformes à l'orchestrateur.

Le résumé ne contient aucun identifiant individuel. Chaque artefact enregistre
assemblage, variation cible, ensemble de variants, unités bp/cM, signatures des
producteurs et SHA-256.

## 11. Schémas à créer ou modifier

Créer au minimum :

- `founder_haplotype_units.schema.json` ;
- `founder_haplotype_boundaries.schema.json` ;
- `founder_haplotype_family_consistency.schema.json` ;
- `founder_haplotype_null_draws.schema.json` ;
- `founder_haplotype_variant_audit.schema.json` ;
- `founder_haplotype_enrichment_summary.schema.json`.

Modifier :

- `pipeline_config.schema.json` ;
- les schémas de sensibilité pour le domaine
  `FOUNDER_HAPLOTYPE_ENRICHMENT` ;
- les schémas de figures, de complétude et de rapport pour le neuvième domaine ;
- les exemples de configuration générique et DOCK6.

Tous les schémas restent stricts (`additionalProperties: false` lorsque le
contrat l'impose) et sont testés sur exemples valides et invalides.

## 12. Modules et changements de code prévus

### 12.1 Nouveau domaine scientifique

Créer :

```text
src/effet_fondateur/founder_enrichment/
  __init__.py
  model.py            # types et résultats immuables
  observed.py         # statistique observée, sans I/O orchestrateur
  null_internal.py    # énumération exhaustive
  null_external.py    # Monte-Carlo reproductible 1000G
  statistics.py       # p empirique, quantiles et intervalles
  publication.py      # tables et résumés

src/effet_fondateur/stages/
  evaluate_founder_haplotype_enrichment.py
```

Le cœur scientifique doit être composé de fonctions pures testables. Le module
de stage se limite à résoudre et valider les artefacts, lancer le calcul puis
publier atomiquement.

### 12.2 Orchestrateur

Modifier `orchestrator/pipeline.py` pour :

- enregistrer `16B` après `16A` et avant `17` ;
- déclarer les dépendances et artefacts requis ;
- conserver la criticité secondaire ;
- intégrer la signature du module et des paramètres ;
- appliquer les codes V2 : entrée/configuration `2`, outil externe `3`, blocage
  scientifique ou intégrité `4`, erreur interne `5` ;
- préserver les tentatives échouées et publier seulement après validation.

Une reprise doit réutiliser 16B uniquement si sa signature, ses entrées et tous
ses SHA-256 sont inchangés.

### 12.3 Référence 1000G

16B ne doit ni télécharger une nouvelle référence non cataloguée ni transmettre
de donnée d'étude. Elle réutilise :

- la référence locale harmonisée de 12 ;
- la liste officielle des 2 504 non-apparentés et la provenance confirmée par
  16A ;
- les mêmes conventions REF/ALT et le même ordre des haplotypes.

Si les sorties actuelles de 16A ne suffisent pas à prouver l'identité exacte du
panel haplotypique consommé, ajouter à 16A un petit descripteur versionné de la
référence locale. Ne pas publier une duplication massive des génotypes si le
BCF harmonisé de 12 et ses empreintes suffisent.

## 13. Intégration aux étapes 17–19

### 13.1 Étape 17 — sensibilités

Ajouter le domaine indépendant :

`FOUNDER_HAPLOTYPE_ENRICHMENT`

Les métriques comparées sont : statut, nombre d'unités, bras gauche/droit,
longueur totale, quantile nul, probabilité empirique et nombre de tirages
évaluables. Aucun vote avec `FOUNDER_IBS`, `LOCAL_LD`, `ROH` ou
`REFERENCE_ANCESTRY` n'est autorisé.

### 13.2 Étape 18 — visualisations

Ajouter un neuvième domaine :

`FOUNDER_HAPLOTYPE_ENRICHMENT`

Figure recommandée :

- histogramme ou fonction de survie de `T_null` interne et externe ;
- ligne verticale pour `T_observed` ;
- bras gauche/droit présentés séparément dans un panneau secondaire ;
- effectif indépendant explicite (`3 familles`) ;
- tirages évaluables et non évaluables ;
- probabilités empiriques et intervalles ;
- mention visible « partage IBS centré cible, pas preuve IBD ».

Une planche comparative peut juxtaposer la courbe de LD de 15 et la distribution
16B, avec deux axes et deux légendes clairement séparés. Elle ne calcule aucune
statistique nouvelle.

Le nombre minimal de domaines de 18 passe de huit à neuf lorsque 16B est
activée. Si 16B est désactivée, elle apparaît comme `SKIPPED`, jamais comme une
ancienne figure réutilisée.

### 13.3 Étape 19 — rapport

Ajouter une section distincte contenant :

- la question scientifique et la statistique préspécifiée ;
- les trois familles comme effectif indépendant ;
- les résultats interne et 1000G séparés ;
- les sensibilités ;
- la relation contextuelle avec le LD de fond de 15 ;
- les limites de phase, carte, densité, démographie et petit effectif ;
- une formulation prudente compatible/incompatible/non concluante.

Le paquet de faits transmis à l'assistance rédactionnelle reste agrégé. Les
haplotypes individuels, génotypes et identifiants ne doivent jamais entrer dans
le prompt ou le rapport partageable.

## 14. Tests obligatoires

### 14.1 Tests unitaires scientifiques

- segment observé symétrique et asymétrique ;
- arrêt au premier manque et à la première discordance ;
- exclusion correcte de la cible de la signature de fond ;
- refus de sélectionner a posteriori la meilleure copie familiale ;
- regroupement de plusieurs porteurs d'une même famille en une unité ;
- tirages internes avec individus distincts ;
- interdiction de sélectionner `H1/H2` du même individu dans un tirage ;
- énumération exhaustive vérifiée sur petit exemple calculable à la main ;
- Monte-Carlo déterministe sous une graine fixe ;
- formule empirique avec correction `+1` ;
- gestion des tirages non évaluables ;
- multiples arrière-plans mutants donnant `MULTIPLE_CARRIER_BACKGROUNDS` ;
- absence de seuil donnant `NOT_CLASSIFIED` ;
- unités et conversions bp/cM.

### 14.2 Tests de contrats et sécurité

- schémas JSON/TSV stricts ;
- assemblage, cible, allèles, ordre et empreintes discordants bloquants ;
- aucune déduction de génotype depuis le statut clinique ou le groupe ;
- aucune donnée de l'étude envoyée au réseau ;
- référence 1000G limitée aux 2 504 individus attendus ;
- chemins confinés au run et refus des liens symboliques dangereux ;
- aucun identifiant individuel dans les résumés, figures ou faits du rapport ;
- cache corrompu ou absent en mode hors ligne bloquant.

### 14.3 Tests orchestrateur et intégration

- ordre `16A -> 16B -> 17` ;
- étape désactivée donnant `SKIPPED` ;
- `NOT_EVALUATED` publié sans échec de l'étape ;
- reprise sans recalcul si les empreintes sont identiques ;
- nouvelle tentative si le code, la graine ou un artefact change ;
- tentative échouée conservée ;
- publication atomique ;
- test synthétique de bout en bout avec trois familles, fond interne et fausse
  référence externe locale ;
- aucun téléchargement ou calcul réel dans les tests unitaires.

### 14.4 Smoke tests externes

Les smoke tests utilisent uniquement des fixtures temporaires et de petits BCF
phasés. `bcftools` peut être exécuté réellement sur ces fixtures. Aucun smoke
test ne doit écrire dans `data/input/`, `data/output/` ou un run réel.

## 15. Critères d'acceptation scientifique

16B est considérée validée uniquement si :

1. la statistique primaire et le modèle nul sont documentés avant analyse ;
2. les unités indépendantes correspondent exactement au gel de 09 et à 13 ;
3. l'observation et chaque tirage nul utilisent le même algorithme ;
4. la cible n'est pas comptée dans la signature de fond ;
5. les apparentés n'augmentent pas artificiellement l'effectif indépendant ;
6. les fonds interne et externe sont séparés ;
7. le Monte-Carlo est reproductible ;
8. les échecs de phase, carte ou intégrité bloquent avant le calcul ;
9. `NOT_EVALUATED` et `NOT_CLASSIFIED` restent visibles ;
10. aucune preuve IBD, origine ethnique ou conclusion fondatrice automatique
    n'est produite ;
11. les figures et le rapport ne recalculent rien ;
12. une revue scientifique humaine reste obligatoire.

## 16. Ordre de développement recommandé

1. Figer ce document après revue scientifique.
2. Écrire les schémas de configuration et de sorties.
3. Créer des fixtures phasées synthétiques avec résultat manuel connu.
4. Implémenter les fonctions pures de segment et statistique observée.
5. Implémenter l'énumération interne exhaustive.
6. Implémenter le Monte-Carlo externe reproductible.
7. Ajouter publication, audit et gardes d'intégrité.
8. Raccorder 16B à l'orchestrateur.
9. Ajouter le domaine à 17.
10. Ajouter le neuvième domaine à 18.
11. Ajouter la section et les faits agrégés à 19.
12. Mettre à jour README, `PIPELINE_V2_PRECODE.md`, les modules documentaires et
    `SESSION.md`.
13. Lancer les tests ciblés, puis la suite moderne complète.
14. Réaliser un audit scientifique transversal `13–16B` avant toute donnée
    réelle.
15. Préparer une nouvelle configuration de run et attendre l'autorisation de
    l'utilisateur avant exécution réelle.

## 17. Migration et run réel

Le run existant
`2026-08-12T072902Z_dock6_reunion_founder_effect_b2075d52` doit rester immuable.
Sa configuration résolue et son manifeste ne doivent pas être modifiés pour y
injecter 16B.

Après implémentation, 16B devra être exécutée dans un nouveau run V2 possédant
une configuration signée qui l'active. Les caches publics vérifiés pourront être
réutilisés, mais aucun ancien dossier d'étape ne sera copié ou déclaré réussi
sans passer par les mécanismes d'intégrité prévus.

Il est interdit de développer, tester ou lancer 16B dans le dossier temporaire
du run 16A actuellement actif.

## 18. Limites anticipées pour le profil DOCK6

- Trois familles donnent une puissance limitée et un résultat exploratoire.
- Les quatorze chromosomes mutants ne sont pas quatorze réplications
  indépendantes.
- La densité ACPA peut surestimer la continuité entre marqueurs non observés.
- Les erreurs ou incertitudes de phase peuvent déplacer les limites.
- La carte génétique et ses intervalles d'interpolation conditionnent les cM.
- Le fond interne est local et de taille modeste.
- Le fond 1000G est large mais démographiquement imparfait pour La Réunion.
- Un haplotype rare peut être compatible avec un fondateur sans démontrer un
  événement fondateur unique.
- Plusieurs haplotypes porteurs peuvent refléter recombinaison ancienne,
  homoplasie, erreur de phase ou origines multiples ; ils doivent être audités,
  pas forcés dans un consensus.

## 19. Références méthodologiques de départ

- Gandolfo LC, Bahlo M, Speed TP. *Dating rare mutations from small samples
  with dense marker data* — partage haplotypique continu et petits effectifs :
  <https://pmc.ncbi.nlm.nih.gov/articles/PMC4125402/>.
- Henden et al. *Identifying individuals with rare disease variants by
  inferring shared ancestral haplotypes from SNP array data* :
  <https://pmc.ncbi.nlm.nih.gov/articles/PMC11970371/>.
- Identity-by-descent analysis appliquée à des événements fondateurs SOD1 :
  <https://pmc.ncbi.nlm.nih.gov/articles/PMC7414871/>.
- Documentation PLINK 1.9 sur le LD, pour maintenir la séparation entre 15 et
  16B : <https://www.cog-genomics.org/plink/1.9/ld>.

Ces références guident la méthode mais ne remplacent pas la validation sur
fixtures synthétiques, la revue statistique et la préspécification du test.
