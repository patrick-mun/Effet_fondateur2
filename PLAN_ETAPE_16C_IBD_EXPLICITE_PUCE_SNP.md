# Plan de mise en place — appel IBD explicite sur données de puce SNP

## 1. Objectif

Ajouter au pipeline V2 une confirmation explicite de l'identité par descendance
(`IBD`) autour du variant cible, sans modifier ni remplacer :

- l'IBS exact conservateur de l'étape `13` ;
- la datation exploratoire de l'étape `14` ;
- l'enrichissement empirique de l'étape `16B`.

La nouvelle branche doit répondre séparément aux questions suivantes :

1. Hap-IBD détecte-t-il un segment contenant le variant cible entre les trois
   familles indépendantes ?
2. Refined IBD confirme-t-il ce segment avec une méthode différente ?
3. Existe-t-il une intersection IBD cohérente entre toutes les paires de
   familles ?
4. Ce signal est-il absent ou rare chez les témoins indépendants comparables ?
5. Les données de puce sont-elles assez denses pour rendre l'analyse
   interprétable ?

Une extension ultérieure utilisera les variants autosomiques de la puce pour
estimer une éventuelle parenté lointaine avec ERSA. Aucun séquençage du génome
entier (`WGS`) n'est requis dans ce plan.

## 2. Périmètre et non-objectifs

### 2.1 Périmètre immédiat

La première livraison ajoute une étape ciblée :

```text
16C_call_target_centered_ibd
```

Elle utilise les données ACPA/CytoScan déjà disponibles, le phasage SHAPEIT5,
la carte génétique GRCh38 et les cohortes gelées du run.

### 2.2 Extension ultérieure

La seconde livraison ajoute une branche genome-wide sur les variants de puce :

```text
16D_prepare_genomewide_ibd_panel
16E_estimate_distant_relatedness
```

Elle produit les segments IBD sur les 22 autosomes, puis alimente ERSA.

### 2.3 Non-objectifs

Cette branche ne doit pas :

- convertir automatiquement un partage IBS en preuve IBD ;
- conclure à un effet fondateur sur la seule présence d'un segment ;
- choisir les seuils après observation du résultat DOCK6 ;
- imputer silencieusement des génotypes et les publier comme observés ;
- compter plusieurs personnes ou chromosomes d'une même famille comme unités
  familiales indépendantes ;
- attribuer une origine ethnique ou géographique ;
- envoyer les données de l'étude vers un service externe.

## 3. Données déjà disponibles

Le run réel validé fournit notamment :

- `75` individus après QC ;
- `62` témoins indépendants ;
- `3` familles porteuses indépendantes ;
- environ `162 183` variants autosomiques communs après QC final chez les
  témoins ;
- `2 418` variants sur le chromosome cible dans le jeu final ;
- `365` variants dans la région préparée ;
- `362` variants dans le BCF phasé après trois exclusions mendéliennes ;
- un segment IBS strict commun de `1,43649051839 cM`, environ `903 kb` et
  `18` marqueurs flanquants ;
- un panel de référence 1000 Genomes GRCh38 disponible dans un cache local
  vérifié.

Le BCF actuellement publié par l'étape `12` couvre seulement la région cible.
Il n'existe pas encore de jeu de puce phasé sur les 22 autosomes pour ERSA.

## 4. Contraintes de faisabilité

Hap-IBD utilise par défaut une longueur minimale de sortie de `2 cM` et une
graine contenant au moins `100` marqueurs. Le segment IBS strict actuel,
`1,44 cM` et `18` marqueurs, risque donc de ne pas être détecté sous les
paramètres standards.

Cette absence éventuelle sera interprétée comme une limite de résolution de la
puce, pas automatiquement comme une réfutation de l'IBS observé.

L'étape doit pouvoir publier :

- `SUPPORTED_EXPLICIT_IBD` lorsque les critères forts sont satisfaits ;
- `PARTIAL_IBD_SUPPORT` lorsque la confirmation est incomplète ;
- `SENSITIVITY_ONLY` lorsque le signal dépend uniquement de seuils
  exploratoires ;
- `NOT_DETECTED` lorsque les données sont suffisantes mais les appels sont
  négatifs ;
- `NOT_EVALUABLE` lorsque la densité ou la qualité ne permettent pas un appel
  fiable.

## 5. Décisions scientifiques préalables

Avant toute exécution réelle, figer dans la configuration :

1. les paramètres primaires de chaque outil ;
2. la grille des sensibilités ;
3. les critères de chevauchement entre outils ;
4. la tolérance maximale entre limites ;
5. le minimum de marqueurs par segment ;
6. la règle de fusion des petits écarts ;
7. les critères de rareté chez les témoins ;
8. le traitement des haplotypes `BOTH` ;
9. le nombre minimal de familles positives ;
10. les règles de non-conclusion.

Aucun seuil ne doit être choisi parce qu'il rend le résultat réel positif.

## 6. Outils retenus

### 6.1 Première livraison

- Hap-IBD : appel principal sur VCF phasé ;
- Refined IBD : confirmation indépendante ;
- `bcftools` : préparation et validation des VCF ;
- Java 17 : environnement d'exécution reproductible.

IBDseq n'est pas prioritaire, car il est principalement conçu pour des données
de séquençage non phasées. RaPID est destiné aux cohortes beaucoup plus grandes
et n'apporte pas d'avantage évident pour `75` individus.

### 6.2 Acquisition et cache

Les outils externes doivent être stockés dans un cache local immuable avec :

- URL officielle ;
- nom et version exacte ;
- licence ;
- SHA-256 attendu et observé ;
- date d'acquisition ;
- version Java ;
- commande de sonde ;
- test synthétique non sensible.

La configuration ne doit pas réutiliser le champ historique
`local_ibd_adapter`, actuellement réservé et explicitement refusé par la
méthode IBS de l'étape `13`. Ajouter un bloc distinct, par exemple :

```yaml
tools:
  explicit_ibd_adapters:
    java_command: java
    hap_ibd_jar: /chemin/epingle/hap-ibd.jar
    hap_ibd_sha256: null
    refined_ibd_jar: /chemin/epingle/refined-ibd.jar
    refined_ibd_sha256: null
```

Le schéma final devra exiger des empreintes réelles et des versions exactes.

## 7. Lot 0 — rapport de faisabilité

Avant d'intégrer les outils au pipeline :

1. inventorier les `2 418` variants du chromosome 19 ;
2. compter les variants possédant une carte cM valide ;
3. compter les variants polymorphes ;
4. mesurer les appels manquants par variant et par individu ;
5. déterminer le nombre de variants complets utilisables par Hap-IBD ;
6. mesurer la densité en marqueurs par cM autour de DOCK6 ;
7. estimer combien de marqueurs couvrent `1`, `1,5`, `2`, `3` et `5 cM` ;
8. vérifier que le chromosome 19 complet peut être phasé avec SHAPEIT5 et la
   référence locale ;
9. publier uniquement des comptes agrégés et un statut de faisabilité.

Livrable proposé :

```text
explicit_ibd_feasibility.json
explicit_ibd_marker_density.tsv
```

Si les paramètres standards ne peuvent matériellement pas être satisfaits,
l'analyse primaire reste autorisée à produire `NOT_EVALUABLE` et les
sensibilités restent clairement secondaires.

## 8. Préparation du chromosome 19

### 8.1 Extension du panel

Préparer un panel chromosome 19 depuis le jeu `target_chromosome_all_qc` de
l'étape `10`, au lieu de se limiter à la fenêtre actuelle de `362` variants.

Le panel doit :

- conserver les `75` individus dans l'ordre maître ;
- conserver le variant cible même s'il est absent de 1000 Genomes ;
- utiliser les coordonnées et allèles GRCh38 canoniques ;
- disposer d'une position cM validée pour chaque variant publié ;
- conserver un audit exhaustif des exclusions ;
- ne jamais extrapoler une carte hors des ancres autorisées.

### 8.2 Phasage

Réutiliser l'adaptateur SHAPEIT5 `phase_common`/`phase_rare` et les pedigrees
existants. Le nouveau phasage doit rester distinct du BCF régional historique
et publier ses propres :

- entrées ;
- paramètres ;
- logs ;
- versions ;
- contrôles Mendel ;
- scores de confiance ;
- empreintes.

### 8.3 Univers commun aux outils

Hap-IBD et Refined IBD doivent consommer exactement le même ensemble de
variants. L'audit doit classer chaque variant :

- `INCLUDED` ;
- `TRUE_MISSING_CALL` ;
- `MONOMORPHIC` ;
- `MENDEL_EXCLUDED` ;
- `ALLELE_MISMATCH` ;
- `MAP_MISSING` ;
- `TOOL_FILTERED_MINOR_ALLELE_COUNT` ;
- autre motif explicitement contracté.

Hap-IBD interdit les allèles manquants. La stratégie primaire doit donc retirer
les marqueurs incomplets, sans compléter leurs génotypes. Une stratégie
d'imputation ne pourra être ajoutée que dans une branche distincte et ne devra
jamais requalifier les dosages imputés en observations.

## 9. Définition des chromosomes mutants

Les attributions de l'étape `12` sont interprétées ainsi :

- `H1` : H1 est le chromosome mutant candidat ;
- `H2` : H2 est le chromosome mutant candidat ;
- `BOTH` : H1 et H2 portent tous deux la mutation et restent deux haplotypes
  distincts.

Pour chaque famille, plusieurs chromosomes mutants peuvent confirmer un
segment, mais la famille compte toujours pour une seule unité indépendante.

Pour conclure à une intersection inter-familiale cohérente, il faut pouvoir
assigner un haplotype mutant unique et continu par famille sur le segment
commun. Il est interdit d'utiliser H1 pour une paire et H2 pour une autre si
aucune assignation globale cohérente n'existe.

Lorsque plusieurs assignations sont possibles :

- toutes les solutions sont auditées ;
- la règle de sélection est appliquée également aux simulations et témoins ;
- le statut ne dépasse pas `PARTIAL_IBD_SUPPORT` si l'ambiguïté n'est pas
  résolue par le pedigree ou les transmissions.

## 10. Appels IBD primaires

### 10.1 Hap-IBD

Exécuter d'abord les paramètres officiels standards, sans adaptation au
résultat DOCK6 :

- `min-seed = 2 cM` ;
- `min-output = 2 cM` ;
- `min-markers = 100` ;
- `min-mac` et tolérances d'écart documentés ;
- nombre de threads et mémoire Java bornés.

### 10.2 Refined IBD

Utiliser les paramètres standards documentés. Les fusions de segments ne sont
autorisées que dans une sortie séparée et selon une règle préspécifiée, par
exemple :

- écart inférieur à `0,6 cM` ;
- au maximum un homozygote discordant ;
- limites avant/après fusion conservées.

### 10.3 Isolation des segments cibles

Pour chaque outil :

1. conserver toutes les sorties brutes dans le run ;
2. extraire les segments contenant la coordonnée du variant cible ;
3. vérifier que les haplotypes impliqués portent explicitement l'allèle cible ;
4. produire les résultats pour les trois paires de familles ;
5. calculer l'intersection commune en bp et cM ;
6. comparer ces limites à l'IBS strict de l'étape `13`, sans l'utiliser pour
   modifier les appels.

## 11. Grille de sensibilité

La grille exacte sera figée après le lot de faisabilité et avant le run réel.
Base proposée :

| Scénario | Longueur minimale | Marqueurs minimaux | Rôle |
|---|---:|---:|---|
| Primaire | `2 cM` | `100` | appel à haute spécificité |
| S1 | `1,5 cM` | calibré | proche du segment IBS observé |
| S2 | `1 cM` | calibré | exploratoire |
| S3 | `3 cM` | standard renforcé | contrôle de spécificité |
| S4 | identique | identique | sans fusion des écarts |
| S5 | identique | identique | fusion limitée préspécifiée |

Les seuils de marqueurs des scénarios S1/S2 doivent provenir des simulations,
pas du nombre de marqueurs observé dans le segment réel.

## 12. Calibration synthétique

Construire des fixtures reproduisant :

- `75` individus ;
- la densité réelle du chromosome 19 ;
- les fréquences alléliques observées ;
- le taux de génotypes manquants ;
- la structure des trois familles ;
- les erreurs et incertitudes de phase plausibles ;
- des segments IBD injectés de `1`, `1,5`, `2`, `3` et `5 cM` ;
- des simulations sans IBD.

Pour chaque outil et scénario, mesurer :

- sensibilité ;
- faux positifs ;
- biais de longueur ;
- précision des limites ;
- fréquence de détection des trois paires ;
- stabilité de l'assignation haplotypique ;
- influence des marqueurs manquants.

Les graines sont fixes et publiées. Les simulations ne doivent contenir aucune
donnée individuelle réelle.

## 13. Analyse des trois familles

Les trois comparaisons indépendantes sont :

```text
FAM001–FAM002
FAM001–FAM003
FAM002–FAM003
```

Pour chaque paire et chaque outil, publier :

- statut de l'appel ;
- chromosomes H1/H2 concernés ;
- limites bp/cM ;
- longueur bp/cM ;
- nombre de marqueurs ;
- nombre d'écarts ou discordances ;
- présence du variant cible ;
- paramètres ayant produit l'appel ;
- rôle primaire ou sensibilité.

Construire ensuite un graphe familial et rechercher une assignation cohérente
d'un haplotype mutant par famille. La seule présence de trois arêtes obtenues
avec des haplotypes incompatibles ne constitue pas un segment fondateur commun.

## 14. Témoins et fréquence de fond

### 14.1 Témoins internes

Utiliser les `62` témoins indépendants gelés, avec :

- le même univers de variants ;
- le même phasage ;
- les mêmes paramètres ;
- la même règle de fusion ;
- la même logique de sélection haplotypique.

Mesurer :

- la fréquence d'un segment IBD couvrant DOCK6 ;
- la fréquence d'un segment au moins aussi long ;
- la fréquence du même haplotype ;
- la fréquence d'une intersection analogue entre trois individus indépendants.

### 14.2 Référence 1000 Genomes

Conserver séparément :

- l'ensemble de la référence ;
- les cinq superpopulations ;
- des témoins appariés descriptivement sur la PCA locale.

Aucune population ne doit être retenue après coup parce qu'elle produit la
probabilité la plus faible.

## 15. Concordance entre outils

Calculer pour Hap-IBD et Refined IBD :

- présence/absence par paire ;
- intersection en bp et cM ;
- différence des limites gauche/droite ;
- proportion de chevauchement ;
- concordance des haplotypes ;
- présence du variant cible ;
- statut sous paramètres primaires ;
- stabilité dans les sensibilités.

La tolérance maximale entre limites doit être fixée avant l'analyse réelle à
partir de la calibration synthétique.

## 16. Règles de classification

### 16.1 `SUPPORTED_EXPLICIT_IBD`

Exiger au minimum :

- Hap-IBD et Refined IBD positifs ;
- les trois paires inter-familiales positives ;
- une assignation haplotypique globale cohérente ;
- une intersection commune contenant le variant cible ;
- des limites concordantes selon la tolérance préspécifiée ;
- un scénario primaire ou validé par calibration ;
- une fréquence de fond sous le seuil défini avant le run.

### 16.2 `PARTIAL_IBD_SUPPORT`

Utiliser lorsque :

- un seul outil est positif ;
- seulement deux paires sur trois sont positives ;
- les limites sont insuffisamment concordantes ;
- l'assignation haplotypique reste ambiguë ;
- la rareté chez les témoins est insuffisamment établie.

### 16.3 `SENSITIVITY_ONLY`

Utiliser lorsque le signal apparaît uniquement avec des paramètres sous les
seuils primaires et dont la spécificité est limitée.

### 16.4 `NOT_DETECTED`

Utiliser lorsque la calibration montre une puissance suffisante mais qu'aucun
appel conforme n'est observé.

### 16.5 `NOT_EVALUABLE`

Utiliser lorsque la densité, les génotypes complets, le phasage ou la carte ne
permettent pas d'atteindre les exigences minimales.

Même `SUPPORTED_EXPLICIT_IBD` ne doit pas être traduit automatiquement en
« effet fondateur prouvé ». La fréquence populationnelle, la généalogie et les
explications concurrentes restent distinctes.

## 17. Artefacts et contrats

Arborescence proposée :

```text
explicit_ibd/
├── explicit_ibd_feasibility.json
├── ibd_marker_audit.tsv
├── ibd_analysis_units.tsv
├── hap_ibd_segments.tsv.gz
├── refined_ibd_segments.tsv.gz
├── target_ibd_pair_matrix.tsv
├── target_ibd_consensus.tsv
├── ibd_tool_concordance.tsv
├── ibd_control_frequency.tsv
├── ibd_parameter_sensitivity.tsv
├── explicit_ibd_summary.json
└── tool_logs/
```

Créer des schémas JSON/TSV stricts pour tous les artefacts publiés. Les tables
contenant des identifiants ou des haplotypes sont `sensitive_genetic`. Le
rapport et les figures n'utilisent que des pseudonymes ou agrégats.

Le résumé doit conserver explicitement :

- `ibs_only_from_step13` ;
- `explicit_ibd_supported` ;
- `founder_effect_proven: false` ;
- `geographic_origin_inferred: false` ;
- statut primaire ;
- statuts des sensibilités ;
- outils, versions, empreintes et paramètres.

## 18. Intégration à l'orchestrateur

Définition proposée :

```text
stage_id: 16C
stage_name: call_target_centered_ibd
critical: false
```

Dépendances minimales :

- `build_sample_registry` ;
- `freeze_cohorts` ;
- `qc_final` ;
- `prepare_target_region` ;
- `phase_target_region` ;
- `infer_founder_haplotype` ;
- `analyze_reference_ancestry`.

L'étape `16C` doit être indépendante de la statistique de `16B` pour l'appel
primaire, mais pourra publier une comparaison descriptive avec le segment IBS
et l'enrichissement déjà calculés.

Après intégration :

- placer `16C` avant `17` ;
- ajouter le domaine `EXPLICIT_IBD` aux sensibilités ;
- ajouter une figure distincte à l'étape `18` ;
- ajouter une section prudente au rapport `19` ;
- incrémenter les contrats de figure et de complétude sans modifier les runs
  historiques.

## 19. Visualisation

Ajouter une figure montrant sur le même axe :

- position du variant cible ;
- segment IBS strict de l'étape `13` ;
- segments Hap-IBD par paire ;
- segments Refined IBD par paire ;
- intersection commune ;
- paramètres primaires versus sensibilités ;
- fréquence agrégée chez les témoins.

La figure doit distinguer visuellement :

```text
IBS observé
→ enrichissement empirique
→ appel IBD explicite
→ hypothèse fondatrice
```

Elle ne doit pas produire de score composite de preuve.

## 20. Tests

### 20.1 Tests unitaires

- parsing Hap-IBD ;
- parsing Refined IBD ;
- contrôle des versions et empreintes ;
- préparation VCF sans allèle manquant ;
- vérification de la carte ;
- sélection des haplotypes mutants ;
- cohérence globale H1/H2 ;
- intersection avec le variant cible ;
- intersection multi-familles ;
- fréquence chez les témoins ;
- concordance entre outils ;
- classification et non-conclusion ;
- blocage des paramètres non préspécifiés.

### 20.2 Tests synthétiques

- aucun segment ;
- segment partagé par deux familles ;
- segment partagé par les trois ;
- outils discordants ;
- segment ne contenant pas la cible ;
- segment fréquent chez les témoins ;
- segment détecté uniquement en sensibilité ;
- haplotypes pairwise incompatibles globalement ;
- densité insuffisante ;
- appels manquants incompatibles avec Hap-IBD.

### 20.3 Tests d'intégration

Les tests ordinaires simulent les outils externes. Un test d'intégration séparé
utilise les vrais JAR sur un petit VCF synthétique non sensible. Aucune donnée
génétique réelle ne doit entrer dans les fixtures versionnées.

## 21. Validation avant run réel

Ordre obligatoire :

1. `git diff --check` ;
2. import du nouveau module ;
3. validation des nouveaux schémas ;
4. tests unitaires ciblés ;
5. smoke synthétique avec adaptateurs simulés ;
6. smoke synthétique avec les vrais outils ;
7. suite moderne complète ;
8. vérification des licences, versions et SHA-256 ;
9. revue scientifique des seuils ;
10. validation d'une nouvelle configuration de run ;
11. lancement réel uniquement après autorisation explicite.

Le run réel terminé `2026-08-13T120222Z_dock6_reunion_founder_effect_8333bf0d`
reste immuable. L'étape `16C` nécessitera un nouveau run ou un mécanisme
d'import inter-run explicitement contracté ; elle ne doit jamais être injectée
manuellement dans ce run historique.

## 22. Extension genome-wide et ERSA

ERSA nécessite des segments IBD sur l'ensemble des autosomes, pas uniquement la
région DOCK6. Cette extension utilise les variants de puce existants.

### 22.1 Préparation

1. construire un univers autosomique commun après QC ;
2. conserver une densité suffisante par chromosome ;
3. phaser les `22` autosomes avec SHAPEIT5 et les cartes GRCh38 ;
4. publier les génotypes observés et les complétions internes séparément ;
5. exécuter Hap-IBD et Refined IBD genome-wide ;
6. masquer les régions connues comme problématiques ou excessivement riches en
   IBD selon une règle sourcée ;
7. calibrer le fond de segments chez les témoins indépendants.

### 22.2 ERSA

Fournir à ERSA :

- nombre de segments partagés par paire ;
- longueurs en cM ;
- distribution du fond chez les témoins ;
- seuil minimal de segment ;
- paramètres et version du modèle.

Publier pour chaque paire inter-familiale :

- degré de parenté le plus vraisemblable ;
- vraisemblances des degrés concurrents ;
- intervalle ou ensemble compatible ;
- statut `RELATED`, `UNRELATED_WITHIN_POWER` ou `NOT_EVALUABLE` ;
- limites de détection.

Une absence de parenté détectée par ERSA ne réfute pas un ancêtre trop ancien
pour la résolution de la puce.

## 23. Risques principaux et réponses

| Risque | Conséquence | Réponse prévue |
|---|---|---|
| Segment de `1,44 cM` sous les seuils standards | absence d'appel primaire | `NOT_EVALUABLE` ou sensibilité clairement séparée |
| Seulement `18` marqueurs dans l'IBS strict | faible résolution | calibration sur densité réelle, chromosome 19 complet |
| Variants manquants incompatibles avec Hap-IBD | perte de marqueurs | filtre complet auditée, aucune imputation cachée |
| Ambiguïté `BOTH` | sélection post hoc des chromosomes | assignation globale cohérente et calibration identique |
| Faux positifs aux seuils courts | surinterprétation | simulations nulles, deux outils, témoins |
| Refined IBD ancien | maintenance limitée | version figée, adaptateur isolé, logs complets |
| ERSA sans WGS | puissance réduite pour parentés très lointaines | modèle de puissance et statut de non-évaluation |
| Populations 1000G imparfaites pour La Réunion | fréquence de fond biaisée | témoins internes prioritaires, analyses par strate séparées |

## 24. Découpage opérationnel

### Jalon A — faisabilité

- inventaire des marqueurs chromosome 19 ;
- densité par cM ;
- compatibilité Hap-IBD/Refined IBD ;
- définition finale des sensibilités ;
- décision `GO`, `GO_SENSITIVITY_ONLY` ou `NOT_EVALUABLE`.

### Jalon B — socle logiciel

- Java et outils épinglés ;
- adaptateurs ;
- parsers ;
- contrats ;
- fixtures synthétiques.

### Jalon C — étape 16C

- préparation/phasing chromosome 19 ;
- appels des deux outils ;
- consensus inter-familial ;
- comparaison aux témoins ;
- résumé, audit et tests.

### Jalon D — restitution

- sensibilité `EXPLICIT_IBD` ;
- figure dédiée ;
- section du rapport ;
- revue scientifique manuelle.

### Jalon E — extension ERSA

- phasage autosomique de la puce ;
- IBD genome-wide ;
- calibration du fond ;
- estimation de parenté distante.

## 25. Critères de fin

La première livraison est terminée lorsque :

- les deux outils s'exécutent depuis des versions et empreintes figées ;
- les données d'entrée sont validées et entièrement auditées ;
- les paramètres primaires et sensibilités sont préspécifiés ;
- les simulations quantifient puissance et faux positifs ;
- les trois paires familiales et les témoins sont analysés sans pseudo-réplication ;
- l'intersection cible et la concordance entre outils sont publiées ;
- tous les contrats et tests passent ;
- le rapport sépare clairement IBS, IBD et effet fondateur ;
- aucun run réel n'a été modifié ou lancé sans autorisation.
