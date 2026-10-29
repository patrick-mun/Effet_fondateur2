# Plan de développement du panel de simulation et de calibration

**Projet :** Effet_fondateur2
**Version du plan :** 0.3
**Date :** 29 octobre 2026
**Statut :** plan de développement — aucune simulation volumineuse autorisée à
ce stade
**Périmètre :** qualification par simulation des étapes `00` à `19`, avec une
priorité scientifique pour `12`, `13`, `14` et `16C`

## Résumé exécutif

- Le panel de calibration sera développé dans un sous-système isolé. Il ne
  lira aucun génotype, identifiant ou résultat individuel de `dock6_real`.
- Deux panels indépendants seront produits : un panel de calibration et un
  panel de validation aveugle. La vérité simulée ne sera jamais transmise au
  pipeline analysé.
- Chaque étape du pipeline fera l'objet d'un rapport autonome en français,
  autoportant en HTML, imprimable en PDF et accompagné de données de rapport
  structurées et d'empreintes SHA-256.
- Le développement commencera par un micro-panel, puis un pilote de
  `500 + 500` individus. Le passage à `2 500 + 2 500` individus dépendra d'un
  contrôle explicite du temps, de la mémoire, du stockage et de la qualité des
  simulations.
- Aucun résultat ne passera automatiquement `16C` à `calibrated: true` et
  aucune correction de datation ne sera appliquée silencieusement. Toute
  qualification exigera un rapport audité, une validation aveugle et une
  décision humaine.

## 1. Objectifs et définition de la réussite

### 1.1 Objectif principal

Construire un environnement reproductible capable de générer des données
génétiques synthétiques avec une vérité connue, de les faire analyser à
l'aveugle par le pipeline réel, puis de mesurer séparément les performances de
chaque étape.

### 1.2 Objectifs secondaires

1. Mesurer les erreurs logicielles, les erreurs de format et les ruptures de
   contrat entre étapes.
2. Quantifier la sensibilité, la spécificité, les biais et les limites de
   résolution des analyses scientifiques.
3. Calibrer les scénarios IBD de `16C` sans utiliser le résultat DOCK6 pour
   choisir les seuils.
4. Évaluer le biais et la couverture des intervalles de datation de `14` sur
   des âges connus.
5. Améliorer le rendu en publiant, pour chaque analyse, un rapport clair qui
   distingue vérité, observation, erreur, incertitude et conclusion.
6. Constituer une suite de non-régression réduite pour les évolutions futures
   du pipeline.

### 1.3 Ce que la calibration ne devra pas prétendre

- Elle ne démontrera pas à elle seule un effet fondateur réel pour DOCK6.
- Elle ne remplacera pas entièrement une validation biologique externe.
- Elle ne garantira pas des intervalles plus étroits. Une calibration honnête
  pourra conduire à élargir les incertitudes.
- Elle ne permettra pas de modifier rétroactivement les runs existants.
- Elle ne devra pas être utilisée pour choisir après coup des paramètres qui
  rendent le résultat DOCK6 positif.

### 1.4 Critère global de réussite

Le panel sera qualifié uniquement si les conditions suivantes sont réunies :

1. séparation vérifiée entre données aveugles et vérité ;
2. absence vérifiée de lecture de `data/input/dock6_real/` ;
3. reproductibilité exacte avec une même configuration et une même graine ;
4. validation de tous les schémas et de toutes les empreintes ;
5. rapports autonomes produits pour chaque étape activée ;
6. critères scientifiques figés avant l'ouverture du panel aveugle ;
7. résultats du panel aveugle conformes aux critères préspécifiés ;
8. rapport transversal bout en bout sans incohérence non résolue ;
9. décision humaine explicite de qualification.

## 2. Sources d'autorité

Les sources sont hiérarchisées afin d'éviter les décisions implicites.

*Tableau 1 — Hiérarchie des sources d'autorité du projet de calibration.*

| Priorité | Source | Autorité |
|---:|---|---|
| 1 | `AGENTS.md` | Sécurité, modification, tests et reproductibilité |
| 2 | `docs/modules/*.md` | Méthodes scientifiques et limites des étapes |
| 3 | `schemas/*.json` | Contrats des artefacts et des configurations |
| 4 | `config/studies/dock6.16b.next.yaml` | Paramètres du profil actuel à reproduire |
| 5 | Code sous `src/effet_fondateur/` | Comportement réellement exécuté |
| 6 | Tests sous `tests/` | Comportements observables déjà garantis |
| 7 | Rapports de calibration futurs | Performances mesurées sur vérité simulée |

Le run DOCK6 du 27 octobre 2026 pourra servir uniquement à inventorier les
dimensions agrégées du problème. Ses résultats individuels et son issue ne
seront jamais utilisés pour régler le simulateur ou les seuils.

## 3. Principes non négociables

### 3.1 Isolation des données

- Interdire au processus de simulation l'accès à `data/input/dock6_real/`.
- N'utiliser que des références publiques, des cartes documentées et des
  paramètres agrégés autorisés.
- Utiliser des identifiants synthétiques sans correspondance avec une personne.
- Ne jamais copier de donnée réelle dans `calibration/` ni dans
  `calibration_workspace/`.

### 3.2 Séparation calibration–validation

- Créer les deux panels avec des graines, des individus et des transmissions
  distincts.
- Autoriser le réglage uniquement sur le panel de calibration.
- Sceller les critères d'acceptation avant d'ouvrir les résultats du panel de
  validation.
- Interdire toute modification de seuil motivée par le panel aveugle.

### 3.3 Séparation vérité–analyse

- Stocker la vérité dans `private_truth/`.
- Stocker les entrées visibles par le pipeline dans `blinded_inputs/`.
- Refuser toute référence à `private_truth/` dans un manifest de pipeline.
- Réaliser la comparaison à la vérité dans un évaluateur séparé, uniquement
  après la fin du run analysé.

### 3.4 Reproductibilité

- Enregistrer les graines de chaque générateur.
- Enregistrer les paramètres, versions et commandes.
- Calculer les SHA-256 de toutes les entrées, vérités et sorties publiées.
- Produire un identifiant immuable pour chaque run de calibration.
- Rendre impossible l'écrasement silencieux d'un run existant.

### 3.5 Conclusion scientifique

- Distinguer `PASS`, `FAIL`, `NOT_EVALUABLE` et `REVIEW_REQUIRED`.
- Ne jamais transformer une absence de détection en preuve d'absence.
- Ne jamais présenter IBS, ROH, LD ou proximité PCA comme une preuve IBD.
- Conserver `founder_effect_proven: false` dans la calibration analytique.
- Marquer toute conclusion issue de ce système comme
  `SIMULATION_CALIBRATED`.

### 3.6 Conventions de développement obligatoires

Tout changement sous `calibration/`, `src/effet_fondateur/calibration/` ou
`tests/calibration/` doit respecter
`calibration/CONVENTIONS_DEVELOPPEMENT.md`. Ce document décline les règles de
`AGENTS.md`, qui reste l'autorité supérieure en cas de divergence.

Avant validation d'une étape de développement :

1. les noms, responsabilités, types et unités doivent être explicites ;
2. toute fonction publique ou scientifique complexe doit documenter ses
   entrées, sorties, hypothèses et erreurs ;
3. les graines, versions, paramètres et empreintes nécessaires à la
   reproductibilité doivent être enregistrés ;
4. les erreurs et codes de retour externes ne doivent jamais être ignorés ;
5. un test ciblé doit couvrir chaque nouveau comportement, cas limite et erreur
   attendue ;
6. les outils externes doivent être simulés en test unitaire et réservés aux
   tests d'intégration explicitement identifiés ;
7. aucune donnée réelle, vérité générée, sortie ou rapport ne doit entrer dans
   Git ;
8. `git diff --check`, les imports et les tests proportionnés au changement
   doivent réussir.

Une porte de passage ne peut pas être déclarée réussie si ces contrôles ne sont
pas documentés. Un changement scientifique doit en plus publier sa référence,
ses hypothèses, une comparaison à une vérité connue et un identifiant de
méthode compatible avec sa portée.

## 4. Architecture cible

```text
calibration/
├── PLAN_DEVELOPPEMENT_PANEL_SIMULATION.md
├── README.md
├── references/
│   ├── public_reference_catalog.json
│   └── CYTOSCAN_750K_ACCEL_MARKERS.md
├── protocoles/
│   ├── PROTOCOLE_GENERAL.md
│   ├── CRITERES_ACCEPTATION.yaml
│   └── REGISTRE_DECISIONS.md
├── configurations/
│   ├── micro_panel.yaml
│   ├── window_design.yaml
│   ├── pilot_calibration.yaml
│   ├── pilot_validation.yaml
│   ├── qualification_calibration.yaml
│   └── qualification_validation.yaml
├── schemas/
│   ├── calibration_config.schema.json
│   ├── truth_manifest.schema.json
│   ├── calibration_run_manifest.schema.json
│   ├── public_reference_catalog.schema.json
│   ├── reference_partition_manifest.schema.json
│   ├── base_haplotype_manifest.schema.json
│   ├── chip_marker_harmonization_manifest.schema.json
│   ├── population_structure_config.schema.json
│   ├── population_structure_manifest.schema.json
│   ├── pedigree_config.schema.json
│   ├── pedigree_manifest.schema.json
│   ├── founder_simulation_config.schema.json
│   ├── founder_simulation_manifest.schema.json
│   ├── simulation_window_config.schema.json
│   ├── simulation_window_manifest.schema.json
│   ├── simulation_window_replay.schema.json
│   ├── stage_evaluation.schema.json
│   └── calibration_decision.schema.json
└── templates/
    └── rapport_scientifique_calibration.html

calibration_workspace/           # généré, entièrement ignoré par Git
└── <calibration_run_id>/
    ├── blinded_inputs/
    ├── private_truth/
    ├── pipeline_runs/
    ├── evaluations/
    ├── rapports/
    ├── figures/
    ├── cache/
    ├── manifest.json
    └── checksums.sha256

src/effet_fondateur/calibration/
├── configuration.py
├── materialize_partitions.py
├── partitioning.py
├── materialize_base_haplotypes.py
├── materialize_chip_observation.py
├── materialize_population_structure.py
├── materialize_pedigrees.py
├── materialize_founder_recombination.py
├── materialize_simulation_windows.py
├── materialize_window_replay_receipt.py
├── simulation/
│   ├── base_haplotypes.py
│   ├── chip_markers.py
│   ├── population_structure.py
│   ├── pedigrees.py
│   ├── founder_recombination.py
│   └── windows.py
├── truth/
├── evaluation/
├── reporting/
└── cli.py

tests/calibration/
├── fixtures/
├── test_configuration.py
├── test_simulation.py
├── test_blinding.py
├── test_evaluation.py
├── test_reporting.py
└── test_end_to_end.py
```

Le code, les schémas, les protocoles et les configurations seront versionnés
sous `calibration/`, `src/` et `tests/`. Les VCF, BCF, vérités, sorties d'outils
et rapports générés seront écrits sous `calibration_workspace/`, entièrement
ignoré par Git. Les chemins de travail seront fournis au lancement et ne seront
pas codés en dur dans une configuration versionnée.

## 5. Dimensionnement progressif

### 5.1 Micro-panel technique

- 24 à 60 individus synthétiques ;
- deux petites strates de population ;
- plusieurs degrés de parenté ;
- au moins un trio ;
- un segment IBD positif et un scénario négatif ;
- une mutation synthétique sans lien avec une personne réelle ;
- un chromosome réduit pour les tests ordinaires ;
- 22 autosomes uniquement pour le test d'intégration identifié.

Ce panel sert à contrôler les formats, les erreurs attendues et les contrats. Il
ne produit aucune conclusion de performance.

### 5.2 Pilote

- 500 individus pour le panel de calibration ;
- 500 individus distincts pour le panel de validation ;
- combinaison de sujets indépendants, individus avec ascendance mélangée,
  trios et pedigrees ;
- scénarios de segments IBD, ROH et haplotypes fondateurs à vérité connue.

Le pilote doit mesurer les ressources réelles avant toute montée en charge.

### 5.3 Qualification

- 2 500 individus pour la calibration ;
- 2 500 individus distincts pour la validation aveugle ;
- cinq strates de référence et une composante mélangée documentée ;
- au moins 100 trios par panel ;
- au moins 200 paires indépendantes par classe de parenté évaluée ;
- 3 000 paires négatives indépendantes pour l'IBD par panel ;
- 400 à 500 répétitions positives par longueur IBD évaluée ;
- 500 répétitions légères par âge simulé pour la datation.

La qualification à 5 000 individus n'est autorisée qu'après le point de
contrôle de ressources du pilote.

### 5.4 Grille scientifique minimale

La grille initiale doit couvrir :

- âges : 5, 10, 15, 20, 30, 40, 50, 75, 100 et 150 générations ;
- familles indépendantes : 3, 5 et 10 ;
- longueurs IBD : 1, 1,5, 2, 3 et 5 cM ;
- scénarios négatifs sans segment injecté ;
- taux de données manquantes et erreurs de génotypage à plusieurs niveaux ;
- erreurs de phase et densités de marqueurs à plusieurs niveaux ;
- positions autosomiques variées, sans calibration limitée à DOCK6.

### 5.5 Dimensionnement des fenêtres de simulation

Les fenêtres utilisées pour mesurer les performances ne seront pas définies par
une largeur physique arbitraire. Pour un âge vrai `g`, la distance `X` au
premier événement sur un bras suit le modèle D14
`X ~ exponentielle(g / 100)` lorsque `X` est exprimée en `cM`. La probabilité
qu'un bras soit censuré par une demi-fenêtre de largeur `W` vaut donc
`exp(-g × W / 100)`.

La série principale utilisera une borne préspecifiée de `0,05` pour la
probabilité théorique de censure de chaque bras. La largeur génétique minimale
de chaque côté de la cible sera par conséquent :

`W_min(g) = -100 × ln(0,05) / g = 299,573 / g cM`.

Cette convention limite à `0,0975` la probabilité théorique qu'au moins un des
deux bras soit censuré, sous l'hypothèse d'indépendance du modèle D14. Elle ne
permet ni de rejeter ni de retirer un tirage effectivement censuré.

*Tableau 2 — Demi-largeur génétique minimale de la série principale.*

| Âge vrai (générations) | Minimum par bras (cM) | Minimum total (cM) |
|---:|---:|---:|
| 5 | 59,915 | 119,829 |
| 10 | 29,957 | 59,915 |
| 15 | 19,972 | 39,943 |
| 20 | 14,979 | 29,957 |
| 30 | 9,986 | 19,972 |
| 40 | 7,489 | 14,979 |
| 50 | 5,991 | 11,983 |
| 75 | 3,994 | 7,988 |
| 100 | 2,996 | 5,991 |
| 150 | 1,997 | 3,994 |

Les limites `bp` seront obtenues par interpolation inverse de la carte publique
figée en D09. L'extrapolation de carte, l'approximation `1 cM = 1 Mb` et le
rabotage silencieux aux bornes d'un chromosome sont interdits. Si une cible ne
dispose pas de la marge exigée sur ses deux côtés, le couple cible–âge est
déclaré `INFEASIBLE_WINDOW` avant tout tirage ; une autre position autosomique
préspecifiée peut être utilisée, mais la fenêtre ne peut pas être choisie après
observation des recombinaisons ou des résultats du pipeline.

Deux séries resteront séparées :

- la série principale, dimensionnée par la formule ci-dessus et utilisée pour
  les métriques de performance et de datation ;
- une série de défi, volontairement tronquée, utilisée uniquement pour vérifier
  la détection et le traitement de la censure. Elle ne contribuera ni au biais,
  ni à la RMSE, ni à une courbe de correction non explicitement stratifiée sur
  la censure.

Pour D15, le domaine devra en plus contenir intégralement chaque segment IBD
injecté de `1`, `1,5`, `2`, `3` ou `5 cM`. Chaque borne vraie devra rester à
l'intérieur de la fenêtre et posséder, à l'extérieur du segment, le nombre de
marqueurs flanquants exigé par le scénario analysé. Le segment observable devra
atteindre le seuil de marqueurs du scénario `16C` correspondant. Un scénario
qui ne satisfait pas ces conditions sera `NOT_EVALUABLE`; sa longueur ou son
seuil ne seront pas modifiés pour le faire passer.

Le micro-panel D33 matérialisera au moins une fenêtre principale techniquement
faisable et une fenêtre de défi censurée. Le pilote et la qualification
couvriront toute la grille d'âges sur plusieurs positions autosomiques
présélectionnées. L'appel des 22 autosomes exigé par `16C` reste un contrôle de
l'univers génomique complet et ne signifie pas que chaque âge doit être simulé
à chaque position.

## 6. Contrat commun des rapports indépendants

### 6.1 Format

Chaque validation produira :

1. un rapport HTML français, autoportant et imprimable en A4 ;
2. un résumé JSON conforme à un schéma ;
3. les tables TSV agrégées nécessaires ;
4. deux à quatre figures SVG ou CSS pertinentes ;
5. un fichier d'empreintes SHA-256 ;
6. un statut final parmi `PASS`, `FAIL`, `NOT_EVALUABLE` et
   `REVIEW_REQUIRED`.

Le rapport HTML sera dérivé du modèle scientifique du skill
`report-formater`, avec CSS linéaire, pagination Paged.js, sommaire cliquable,
pied de page et numérotation `Page n / total`.

### 6.2 Sections obligatoires

Chaque rapport contiendra exactement les sections suivantes :

1. titre, identifiant, date, version et étape ;
2. résumé exécutif en trois à cinq points ;
3. question évaluée ;
4. sources et autorité ;
5. configuration, graines et versions ;
6. vérité simulée ;
7. entrées aveugles réellement consommées ;
8. méthode d'évaluation ;
9. résultats et intervalles d'incertitude ;
10. comparaison aux critères préspécifiés ;
11. anomalies, exclusions et écarts au protocole ;
12. limites ;
13. décision et action suivante ;
14. annexes de provenance et empreintes.

### 6.3 Règles de rendu

- Ne publier que des résultats agrégés.
- Numéroter sans rupture toutes les figures et tables.
- Citer la source de chaque valeur et de chaque figure.
- Employer au maximum trois niveaux de titres.
- Utiliser des cartes de chiffres clés pour trois ou quatre métriques au plus.
- Utiliser une courbe de calibration pour la datation.
- Utiliser une matrice ou des barres pour les classifications.
- Utiliser un diagramme de flux pour la provenance et l'aveuglement.
- Ne jamais dépendre d'une bibliothèque graphique réseau au moment du rendu.

### 6.4 Contrôle qualité éditorial

Avant publication, chaque rapport devra réussir :

- vérification de toutes les sections promises ;
- absence de texte provisoire ou de champ non résolu ;
- continuité des numéros de figures et tables ;
- validité des renvois internes ;
- accents, unités et dates cohérents ;
- absence de débordement en A4 ;
- figures et tables non coupées ;
- pied de page présent sur toutes les pages de contenu ;
- couverture et sommaire sans pied de page ;
- contrôle du PDF exporté sur une seconde lecture.

### 6.5 Convention de nommage

Chaque rapport sera rangé dans un dossier autonome :

```text
calibration_workspace/<calibration_run_id>/rapports/<report_id>_<nom_court>/
├── index.html
├── report_data.json
├── tables/
├── figures/
└── checksums.sha256
```

Le fichier `index.html` restera autoportant. Les autres fichiers permettront de
recalculer son contenu et d'en vérifier la provenance sans extraire le HTML.

## 7. Registre des rapports à produire

Chaque étape possède son propre rapport. Un rapport transversal ne remplace
jamais un rapport d'étape.

*Tableau 3 — Registre contractuel des rapports indépendants.*

| ID | Étape évaluée | Rapport indépendant |
|---|---|---|
| CAL-REP-00 | `00_initialize_run` | Initialisation, environnement et immutabilité |
| CAL-REP-01 | `01_validate_sources` | Sources, formats et empreintes |
| CAL-REP-02 | `02_build_sample_registry` | Registre synthétique et alignement |
| CAL-REP-03 | `03_convert_acpa` | Fidélité de conversion génotypique |
| CAL-REP-04 | `04_prepare_target_variant_dataset` | Mutation synthétique et génotypes cibles |
| CAL-REP-05 | `05_qc_preliminary` | Contrôle qualité préliminaire |
| CAL-REP-06 | `06_build_kinship_panel` | Panel de parenté et sélection des marqueurs |
| CAL-REP-07 | `07_infer_kinship` | Classification des degrés de parenté |
| CAL-REP-08 | `08_analyze_population_structure` | Structure et atypiques de population |
| CAL-REP-09 | `09_freeze_cohorts` | Gel des cohortes et indépendance familiale |
| CAL-REP-10 | `10_qc_final` | Contrôle qualité final |
| CAL-REP-10A | `10A_prepare_autosomal_phasing_panel` | Univers autosomique commun |
| CAL-REP-11 | `11_prepare_target_region` | Région cible et carte génétique |
| CAL-REP-12 | `12_phase_target_region` | Phase régionale et orientation de la mutation |
| CAL-REP-12A | `12A_phase_autosomal_panel` | Phase autosomique et erreurs de commutation |
| CAL-REP-13 | `13_infer_founder_haplotype` | Segment fondateur IBS et bornes |
| CAL-REP-14 | `14_estimate_variant_age` | Biais et couverture de la datation Gamma |
| CAL-REP-15 | `15_analyze_local_ld` | Fidélité du LD local |
| CAL-REP-16 | `16_analyze_roh` | Détection et longueur des ROH |
| CAL-REP-16A | `16A_analyze_reference_ancestry` | Positionnement global et local |
| CAL-REP-16B | `16B_evaluate_founder_haplotype_enrichment` | Calibration de l'enrichissement IBS |
| CAL-REP-16C | `16C_call_explicit_ibd` | Sensibilité, spécificité et concordance IBD |
| CAL-REP-17 | `17_run_sensitivity_analyses` | Stabilité inter-runs et scénarios |
| CAL-REP-18 | `18_build_visualizations` | Fidélité et complétude des figures |
| CAL-REP-19 | `19_build_report` | Fidélité du rapport scientifique final |
| CAL-REP-E2E | Pipeline complet | Validation bout en bout et absence de fuite |
| CAL-REP-DEC | Décision | Qualification, restrictions et durée de validité |

## 8. Points de contrôle par domaine

### 8.1 Intégrité technique — étapes 00 à 04

Le domaine réussit si :

- tous les fichiers attendus sont présents ;
- les schémas sont valides ;
- les SHA-256 correspondent ;
- chaque individu synthétique est aligné sans perte ni duplication ;
- les génotypes après conversion sont identiques à la vérité ;
- aucune donnée réelle n'est ouverte ;
- une même graine reproduit exactement les mêmes artefacts.

### 8.2 Contrôle qualité — étapes 05 et 10

Le domaine réussit si :

- les défauts injectés sont détectés et classés selon le contrat ;
- les exclusions attendues correspondent à la vérité ;
- aucun défaut non injecté n'est fabriqué par la conversion ;
- les cas non évaluables restent explicitement non évaluables ;
- les seuils ne sont pas modifiés après ouverture du panel aveugle.

Le contrôle de discordance de sexe restera `NOT_EVALUABLE` tant qu'un univers
de chromosomes sexuels avec vérité adaptée n'aura pas été ajouté.

### 8.3 Parenté — étapes 06 et 07

Le domaine réussit si :

- la vérité couvre les degrés testés et les sujets non apparentés ;
- la matrice de confusion et les intervalles sont publiés ;
- les erreurs sont stratifiées par ascendance, données manquantes et densité ;
- la proposition d'ensemble indépendant ne contient aucune paire interdite ;
- les performances du panel aveugle respectent les seuils figés après pilote.

### 8.4 Structure de population et cohortes — étapes 08 et 09

Le domaine réussit si :

- les groupes simulés et les mélanges connus sont traçables ;
- les atypiques injectés sont retrouvés sans règle circulaire ;
- les décisions de cohorte sont reproductibles ;
- les familles indépendantes restent indépendantes après gel ;
- aucune proximité PCA n'est interprétée comme une origine géographique.

### 8.5 Univers de phasage — étapes 10A et 11

Le domaine réussit si :

- les mêmes variants et le même ordre d'échantillons sont utilisés ;
- les 22 autosomes sont présents dans le test d'intégration complet ;
- les cartes génétiques et positions `bp`/`cM` sont cohérentes ;
- la cible et ses marqueurs flanquants sont conservés ;
- toute exclusion possède un motif audité.

### 8.6 Phasage — étapes 12 et 12A

Le domaine réussit si :

- le taux d'erreur de commutation est calculé contre la phase vraie ;
- les erreurs ponctuelles de phase de variants rares sont séparées des
  commutations de blocs ;
- l'orientation de la mutation synthétique est correctement classée ;
- les résultats sont stratifiés par fréquence allélique, densité et manque ;
- les deux panels aveugles restent entièrement séparés.

### 8.7 Haplotype fondateur — étape 13

Le domaine réussit si :

- la mutation vraie reste dans le segment récupéré lorsqu'il est détectable ;
- la sensibilité et la spécificité sont mesurées sur scénarios positifs et
  négatifs ;
- les erreurs de bornes gauche et droite sont publiées en `bp` et `cM` ;
- les familles absentes, discordantes ou mal phasées sont exclues avec motif ;
- un segment IBS n'est jamais requalifié automatiquement en IBD.

### 8.8 Datation — étape 14

Le domaine réussit si :

- les estimations sont comparées à l'âge vrai pour chaque point de la grille ;
- le biais, l'erreur absolue, la RMSE et la couverture à 95 % sont publiés ;
- la calibration est séparée pour 3, 5 et 10 familles indépendantes ;
- la stabilité leave-one-family-out est publiée ;
- l'effet des erreurs de phase, de carte et de bornes est quantifié ;
- toute correction proposée reçoit un nouvel identifiant de méthode ;
- le modèle Gamma historique reste visible comme référence non modifiée.

La cible proposée pour la couverture à 95 % sera figée après le pilote et avant
la validation aveugle. Elle devra être accompagnée de son intervalle binomial,
pas seulement d'une valeur ponctuelle.

### 8.9 LD — étape 15

Le domaine réussit si :

- les `r²` et distances sont comparés aux valeurs calculées sur la vérité ;
- les erreurs sont stratifiées par fréquence allélique et densité ;
- le fond génomique est construit sans chromosome cible lorsque le contrat
  l'exige ;
- l'influence de la parenté est montrée séparément ;
- le LD ne nourrit pas silencieusement la datation.

### 8.10 ROH — étape 16

Le domaine réussit si :

- les ROH injectées sont comparées aux appels PLINK ;
- sensibilité, faux positifs et erreurs de longueur sont publiés ;
- les résultats sont stratifiés par longueur et densité ;
- les ROH dues à l'autozygotie sont distinguées du segment fondateur partagé ;
- les ROH ne sont pas utilisées directement comme longueurs Gamma.

### 8.11 Ascendance de référence — étape 16A

Le domaine réussit si :

- les composantes simulées sont connues et séparées du panel de référence ;
- les positions globales et locales sont stables entre répétitions ;
- les mélanges simulés sont placés entre références selon la vérité ;
- les limites dues aux références larges sont explicitement publiées ;
- aucune étiquette de population n'est convertie en origine individuelle.

### 8.12 Enrichissement haplotypique — étape 16B

Le domaine réussit si :

- les probabilités empiriques sont calibrées sous scénarios nuls ;
- la puissance est publiée pour plusieurs longueurs et nombres de familles ;
- les tirages non évaluables sont comptés séparément ;
- la cible est exclue de la signature comme dans le contrat réel ;
- aucune classe positive n'est introduite sans seuil préspécifié.

Tant que `empirical_classification_threshold` reste nul, le rapport doit
conserver `NOT_CLASSIFIED` même si la simulation montre une forte puissance.

### 8.13 IBD explicite — étape 16C

Le domaine réussit si :

- Hap-IBD et Refined IBD reçoivent le même univers et la même carte ;
- les 22 autosomes sont appelés avant interrogation de la cible ;
- les paires négatives mesurent la fréquence de faux support concordant ;
- les segments positifs mesurent sensibilité et erreur de bornes ;
- les scénarios 2 cM/100, 1,5 cM/50 et 1 cM/30 sont séparés ;
- le nombre de paires et les intervalles de confiance sont publiés ;
- la borne supérieure de fréquence de faux positif est comparée au maximum
  préspécifié de `0,05` ;
- le panel aveugle confirme la calibration sans ajustement post hoc.

`calibrated: true` ne pourra être proposé que si `CAL-REP-16C` et
`CAL-REP-E2E` sont qualifiés et si le rapport de calibration est référencé par
chemin et SHA-256 dans une nouvelle configuration.

### 8.14 Sensibilités, figures et rapport final — étapes 17 à 19

Le domaine réussit si :

- les comparaisons inter-runs de `17` retrouvent les différences injectées ;
- aucune différence absente n'est inventée ;
- `18` ne recalcule aucune statistique scientifique ;
- chaque figure possède une provenance et une empreinte ;
- `19` reprend exactement les faits contrôlés des rapports sources ;
- les statuts non évaluables et exploratoires restent visibles ;
- le rapport final passe le contrôle éditorial et l'export PDF.

## 9. Étapes de développement, une à une

### D01 — Geler le périmètre initial

**Action :** confirmer que la première tranche couvre `12`, `13`, `14` et
`16C`, sans exclure la qualification future des autres étapes.
**Livrable :** entrée signée dans `REGISTRE_DECISIONS.md`.
**Réussite :** aucune ambiguïté sur les étapes prioritaires ou différées.

### D02 — Créer l'arborescence isolée

**Action :** créer les dossiers suivis sous `calibration/` et la zone locale
`calibration_workspace/` entièrement ignorée. Conserver aussi les protections
pour les anciens chemins `calibration/runs/` et `calibration/local/`.
**Livrable :** structure décrite à la section 4.
**Réussite :** aucun fichier généré n'apparaît dans `git status`.

### D03 — Rédiger le protocole général

**Action :** figer objectifs, unités statistiques, scénarios, aveuglement et
règles de décision.
**Livrable :** `PROTOCOLE_GENERAL.md`.
**Réussite :** revue humaine du protocole avant développement scientifique.

### D04 — Rédiger les critères d'acceptation

**Action :** définir les métriques, intervalles et portes de passage par étape.
**Livrable :** `CRITERES_ACCEPTATION.yaml`.
**Réussite :** tous les critères possèdent une unité, une direction et une
règle pour les valeurs manquantes.

### D05 — Définir les schémas de configuration

**Action :** créer et tester le schéma des panels, graines et scénarios.
**Livrable :** `calibration_config.schema.json`.
**Réussite :** refus des graines absentes, effectifs incohérents et chemins
interdits.

### D06 — Définir le manifest de vérité

**Action :** formaliser parenté, ascendance, phase, mutation, segments IBD,
ROH, âge et erreurs injectées.
**Livrable :** `truth_manifest.schema.json`.
**Réussite :** toute vérité consommée par un évaluateur est déclarée et typée.

### D07 — Définir le manifest de run de calibration

**Action :** enregistrer sources, versions, commandes, graines, fichiers et
empreintes.
**Livrable :** `calibration_run_manifest.schema.json`.
**Réussite :** un run peut être reproduit sans connaissance orale.

### D08 — Construire les gardes d'isolation

**Action :** refuser les chemins réels, les liens symboliques non autorisés et
les collisions d'identifiants.
**Livrable :** module de sécurité et tests dédiés.
**Réussite :** les tests démontrent que `dock6_real` ne peut pas être ouvert.

### D09 — Inventorier et figer les références publiques

**Action :** vérifier licences, versions, assemblage, cartes et empreintes.
**Livrable :** catalogue versionné des références autorisées.
**Réussite :** chaque source possède identifiant, version, origine et SHA-256.

Le catalogue D09 distingue l'empreinte de la source épinglée de celle des
charges utiles acquises. Lorsque le fournisseur ne publie que des MD5 par
fichier, le SHA-256 du manifeste officiel fige la release, puis chaque fichier
téléchargé doit être contrôlé par l'empreinte fournisseur et recevoir un
SHA-256 local dans le manifest D07 avant utilisation. Une autorisation d'usage
local ne vaut pas autorisation de redistribution.

### D10 — Séparer donneurs simulés et référence de phasage

**Action :** créer des partitions disjointes pour éviter une validation
optimiste.
**Livrable :** manifest des partitions.
**Réussite :** intersection nulle des identifiants entre partitions interdites.

Le manifest D10 reste privé dans le workspace : il publie seulement les rôles,
effectifs et empreintes des listes placées sous `private_truth/partitions/`.
La séparation est contrôlée par composante connexe de pedigree ou de parenté,
pas seulement par identifiant. Un groupe contenant à la fois un candidat à la
référence de phasage et un autre donneur est retiré en bloc de la référence.
La référence de phasage peut être commune aux panels de calibration et de
validation, mais les groupes donneurs restent distincts entre eux et de cette
référence.

### D11 — Générer les haplotypes de base

**Action :** produire les deux chromosomes de chaque individu synthétique en
préservant fréquences, LD et carte.
**Livrable :** vérité haplotypique dense D11a, génotypes dérivés et univers
observable CytoScan 750K Accel D11b.
**Réussite :** concordance génotype–haplotype à 100 %, liaisons puce–1000G
non ambiguës et nombre de variants observable égal au nombre de marqueurs liés.

La méthode D11 réapparie, au sein d'une même strate publique, les deux copies
de donneurs distincts. Chaque haplotype donneur est utilisé exactement une fois
et n'est pas modifié : les fréquences alléliques et le LD haplotypique du pool
sélectionné sont donc invariants par permutation. La paire diploïde d'un
échantillon public n'est jamais reproduite. La provenance donneur–copie reste
privée, tandis que le manifest ne publie que des effectifs et empreintes.

À ce stade, « cohérence mendélienne » désigne uniquement la dérivation exacte
du génotype non phasé à partir de ses deux haplotypes. Aucune transmission
parent–enfant n'existe encore en D11 ; son évaluation est explicitement différée
à D13 pour ne pas revendiquer un contrôle biologique qui n'est pas applicable.
Les positions `cM` sont interpolées linéairement entre points de carte, sans
extrapolation hors du domaine couvert.

D11 est séparée en deux couches. D11a conserve la vérité 1000G dense nécessaire
aux transmissions, recombinaisons et mesures de référence. D11b construit une
observation distincte limitée aux SNP réellement présents sur la
`CytoScan750K_Accel`, version `NA36.r1`, assemblage `GRCh38`. Le paquet
constructeur doit posséder une entrée Thermo Fisher revue dans le catalogue
D09. D11b lit directement sa base d'annotations, valide son contenu et produit
dans le workspace privé la table normalisée ; le paquet et cette table sont
tous deux liés par SHA-256. Aucun cache ni export d'échantillon réel ne peut
définir l'univers de marqueurs.

L'harmonisation D11b exige une coordonnée exacte et des allèles compatibles
avec l'orientation ou le brin complémentaire. Elle exclut et compte les
positions dupliquées, les variants 1000G multiples, les allèles incompatibles
et les SNP palindromiques de brin inconnu. Le VCF observable conserve
l'orientation REF/ALT de 1000G et enregistre la relation de sonde séparément.
Le test D11a de 100 kb et son run D11b réel constituent une intégration
technique bornée : le paquet officiel est lié et le prérequis de D12 est
satisfait pour ce domaine. Avec un seul marqueur de puce dans la fenêtre, ce
diagnostic ne remplace ni le futur micro-panel D33 ni son contrôle de densité.

### D12 — Générer la structure de population

**Action :** produire strates et mélanges avec proportions connues.
**Livrable :** vérité d'ascendance et panel aveugle.
**Réussite :** proportions valides, total égal à 1 et séparation des panels.

**Prérequis bloquant :** un manifest D11a et un manifest D11b valides, liés au
même panel et au même domaine. D12 ne peut pas utiliser le seul VCF dense comme
entrée observable du pipeline.

**Méthode retenue :** D12 conserve les diploïdes D11 non sélectionnés comme
strates pures. Pour un mélange entre deux strates, deux diploïdes sources
produisent les paires réciproques `A1|B2` et `B1|A2`. Chaque copie haplotypique
est ainsi utilisée exactement une fois et le mélange individuel vaut exactement
50/50 sur le domaine. Une proportion inégale ou plus de deux composantes
nécessiterait une mosaïque chromosomique ; elle est interdite à D12 et différée
à D14.

La configuration D12 est finalisée, empreintée et explicite sur les strates,
les proportions, l'effectif pair de chaque mélange et la graine
`haplotypes`. La sortie contient séparément le VCF phasé dense privé, la vérité
d'ascendance, la provenance limitée aux identifiants synthétiques D11 et un VCF
observable réduit aux marqueurs D11b. Ce dernier ne contient aucune étiquette
d'ascendance. Les rôles `CALIBRATION` et `VALIDATION` utilisent des préfixes
distincts et ne peuvent être mélangés dans un même run.

### D13 — Générer pedigrees et parentés

**Action :** simuler trios et relations jusqu'au troisième degré.
**Livrable :** pedigrees vrais et génotypes descendants.
**Réussite :** transmissions mendéliennes sans incohérence non injectée.

**Prérequis bloquant :** un manifest D12 finalisé dont le panel dense, le panel
observable, la vérité d'ascendance et la provenance haplotypique sont tous liés
par SHA-256. D13 refuse de déduire l'ascendance d'un descendant à partir de la
seule étiquette diploïde de ses parents.

**Méthode retenue :** une unité familiale de huit individus part de quatre
fondateurs D12 distincts et appartenant à la même strate. Un couple de
grands-parents produit deux germains ; chacun forme ensuite un couple avec un
conjoint fondateur distinct et produit un enfant. Une unité couvre ainsi quatre
trios, les premier et deuxième degrés, puis une paire de cousins germains au
troisième degré. Les coefficients de parenté théoriques `0,25`, `0,125` et
`0,0625` sont calculés récursivement depuis le pedigree ; ils ne sont ni
estimés à partir des génotypes, ni assimilés à une quantité IBD réalisée sur la
fenêtre.

D13 tire, avec la graine `pedigrees`, une copie haplotypique entière de chaque
parent. Aucune recombinaison n'est autorisée avant D14. La provenance D13→D12
rend chaque transmission vérifiable et l'ascendance du descendant est calculée
depuis les deux strates haplotypiques effectivement transmises. Un descendant
de parents 50/50 peut donc être pur pour une composante ou rester 50/50 sur le
domaine ; une moyenne parentale systématique est interdite.

Les vérités privées comprennent le pedigree, les transmissions, la parenté
théorique, la provenance haplotypique et l'ascendance réalisée. Le VCF aveugle
conserve seulement l'univers de marqueurs D12 et ne révèle ni famille, ni
ascendance, ni identifiant D12. La porte D13 exige une bijection exacte des
fondateurs, des familles disjointes, zéro erreur mendélienne, les quatre classes
de relation attendues, des empreintes valides et une concordance intégrale
entre provenance et ascendance.

### D14 — Simuler le fondateur et les recombinaisons

**Action :** transmettre un haplotype fondateur pendant un âge connu.
**Livrable :** bornes vraies gauche/droite par famille et génération.
**Réussite :** distances `bp`, `cM` et générations cohérentes.

**Prérequis bloquant :** les manifests D11a, D12 et D13 doivent être finalisés,
liés par identifiant et SHA-256, et décrire exactement le même panel, chromosome
et domaine. La carte interpolée D11a est l'unique autorité pour les conversions
`bp`↔`cM` ; aucune vitesse moyenne codée en dur n'est autorisée.

**Méthode retenue :** pour un âge vrai `g`, la distance au premier événement de
recombinaison de chaque côté du variant suit indépendamment une loi
exponentielle de taux `g` par Morgan, soit `g / 100` par cM. Si le premier
événement tombe hors du domaine simulé, le bras est censuré au bord. Une limite
censurée n'est jamais comptée comme un événement de recombinaison et sa longueur
ne doit pas être interprétée comme la longueur ancestrale complète.

Un haplotype D13 extérieur aux familles sert de matrice fondatrice commune. Une
origine haplotypique distincte est choisie dans chaque famille à partir d'un
descendant terminal ; le segment fondateur est mosaïqué sur cette origine, puis
propagé à tous les membres qui l'ont reçue. Cette propagation par origine D12
préserve exactement les transmissions mendéliennes D13. Les origines ne peuvent
pas traverser deux familles et la matrice fondatrice ne peut appartenir à une
famille simulée.

Le variant cible synthétique doit correspondre à un SNV biallélique déjà
présent et correctement orienté dans le VCF. Son fond allélique est remis à la
référence, puis l'allèle alternatif est posé uniquement sur les copies
porteuses. Cette règle évite d'inventer une base GRCh38 sans référence FASTA.
La vérité privée distingue segments, censures, porteurs, mutation et mosaïque
d'ascendance locale. Le VCF aveugle conserve uniquement les marqueurs
observables et ne reçoit aucun libellé de vérité.

Le diagnostic de 100 kb à 150 générations contrôle la chaîne, l'injection et la
censure, mais il ne suffit pas à valider la distribution des recombinaisons :
sa longueur génétique totale est seulement d'environ `0,136 cM`. Une fenêtre
plus large reste obligatoire pour le micro-panel et la qualification.

### D14b — Construire les fenêtres génétiques de simulation

**Action :** sélectionner puis matérialiser, avant tout tirage D15, des domaines
génomiques dimensionnés en `cM` selon la section 5.5.
**Livrable :** configuration versionnée, manifest privé des fenêtres, tables de
coordonnées `bp`/`cM`, inventaires des variants denses et observables,
empreintes SHA-256 et reçu privé du rejeu D11a–D14.
**Réussite :** chaque fenêtre principale respecte sa demi-largeur théorique,
chaque fenêtre de défi est explicitement étiquetée, et aucune cible n'est
choisie selon un résultat simulé ou un résultat DOCK6.

La configuration doit fixer avant matérialisation : le rôle du panel, l'âge ou
le scénario IBD, le chromosome, la cible synthétique, la carte D09, la source
VCF D11, le profil de puce D11b, la probabilité maximale de censure par bras et
les seuils de marqueurs. Les positions candidates seront déterminées uniquement
à partir des références publiques et des règles de densité. Leur ordre de
sélection sera déterministe et lié aux graines déjà réservées ; aucune nouvelle
graine implicite n'est autorisée.

Le manifest doit enregistrer pour chaque fenêtre : les limites physiques et
génétiques, les distances gauche/droite, la probabilité théorique de censure,
le statut `PRIMARY` ou `CENSORING_CHALLENGE`, les nombres de variants denses et
observables de chaque côté, les seuils IBD applicables et les empreintes des
sources. Les identifiants de donneurs et les chemins locaux absolus y sont
interdits.

La matérialisation échoue sans produire de run accepté si l'une des conditions
suivantes survient : carte non couvrante, extrapolation nécessaire, cible hors
VCF, cible non biallélique, marge génétique insuffisante, marqueurs flanquants
insuffisants, seuil `16C` impossible à atteindre, discordance d'assemblage,
source non figée, rôle de panel ambigu ou collision de chemin de sortie. Une
fenêtre calibration et une fenêtre validation peuvent partager la même
définition publique et la même coordonnée pour permettre une comparaison
appariée, mais leurs identifiants de scénario, graines, haplotypes, individus,
transmissions et artefacts générés doivent rester distincts.

Après acceptation d'une définition de fenêtre, la chaîne D11a–D14 est rejouée
sur ce domaine sous de nouveaux identifiants de run. Un reçu D14b distinct lie
par SHA-256 le manifest de fenêtre et les cinq manifests D11a, D11b, D12, D13
et D14 ; aucun manifest finalisé n'est modifié rétroactivement. Un tirage
principal effectivement censuré reste dans la vérité et n'est jamais rejoué
dans le but d'obtenir un événement interne.

Les tests ciblés devront au minimum vérifier la formule, les dix âges de la
grille, l'interpolation inverse, les bornes chromosomiques, les cartes non
monotones, les densités insuffisantes, la séparation des deux séries, le refus
de `dock6_real`, le déterminisme et la séparation calibration–validation. La
porte D14b exige aussi un rejeu octet pour octet identique avec la même
configuration et la même graine.

### D15 — Injecter les scénarios IBD

**Action :** générer segments positifs et négatifs aux longueurs prévues.
**Livrable :** vérité IBD par paire, chromosome et scénario.
**Réussite :** aucune cible positive dans les scénarios négatifs.

**Prérequis bloquant :** un manifest D14b finalisé et valide. D15 ne peut ni
élargir une fenêtre à la volée, ni déplacer une cible après un tirage, ni
reclasser une fenêtre de défi comme fenêtre principale.

**Méthode D15 v1 :** `disjoint_unrelated_haplotype_copy_v1`. Le manifest de
fenêtre D14b constitue la configuration scientifique scellée des longueurs et
seuils. Le reçu D14b lie D11a et D14 ; les artefacts consommés sont vérifiés
par SHA-256 avant lecture. Une dérivation déterministe SHA-256 de la graine
`recombination`, du rôle de panel et de l'identifiant de méthode fournit la
graine D15, sans nouvelle graine maîtresse implicite.

Chaque longueur reçoit deux sujets positifs, deux sujets négatifs et une
matrice distincte, sélectionnés sans remise parmi les individus hors pedigree
et hors porteurs D14. Tous les sujets et matrices sont disjoints entre
scénarios. H1 de la matrice est copié sur H1 des deux positifs dans un intervalle
centré sur la cible ; H2 et les génotypes hors intervalle restent inchangés.
Les transmissions D13/D14 sont contrôlées sur tous les variants. Les seuils
observables sont revérifiés après conversion des bornes en bp ; une densité
insuffisante bloque l'acceptation du run.

La vérité distingue longueur demandée et bornes réalisées en bp/cM. Elle
décrit l'origine commune injectée, pas les limites maximales d'un appel.
Les négatifs ne reçoivent aucune injection et possèdent des bornes nulles :
le fond 1000G peut conserver un partage ancestral inconnu ou de l'IBS. Il est
interdit de filtrer les négatifs en fonction des sorties des outils IBD.
La fréquence d'appels sur ces paires mesure donc le fond sans injection ; elle
ne doit pas être étiquetée comme une spécificité contre une absence absolue
d'IBD connue. Ce point devra rester explicite dans l'évaluateur D29.

La matrice et les bornes inscrites dans la vérité D15 constituent une surcharge
locale d'origine sur le fond D14. Si matrice et receveur sont de strates
différentes, les évaluateurs d'ascendance devront appliquer cette surcharge
avant comparaison ; les étiquettes d'ascendance antérieures ne sont pas une
vérité inchangée pour les copies modifiées.

Le run régional technique contient une paire positive et une paire négative
par longueur. Ces sujets ne deviennent pas automatiquement porteurs de la
mutation ou représentants familiaux. Il ne satisfait ni le nombre de familles
requis pour le consensus 16C, ni l'appel des 22 autosomes, ni les effectifs de
qualification ; `calibrated` reste `false`. Les copies contrôlent le générateur
et les futurs évaluateurs par paire. Une origine synthétique distincte ne
démontre pas à elle seule l'indépendance du fond de population.

Le rapport indépendant `CAL-DEV-D15` publie les agrégats, limites et empreintes
de cette validation de développement. Le rapport de performance `CAL-REP-16C`
sera produit par l'évaluateur ultérieur ; il ne peut être remplacé par ce
contrôle technique.

### D16 — Injecter les ROH

**Action :** créer des segments autozygotes de longueurs connues, distincts ou
superposés au locus fondateur selon le scénario.
**Livrable :** vérité ROH.
**Réussite :** distinction explicite entre ROH et segment IBD partagé.

**Prérequis bloquant :** un run D15 finalisé, lié par SHA-256 aux mêmes
manifests D14 et D11a, et une configuration D16 finalisée de même rôle de
panel. Les entrées sont vérifiées avant lecture ; la fenêtre D14b et les
seuils de marqueurs ne sont pas modifiés.

**Méthode D16 v1 :** `within_sample_h1_to_h2_roh_copy_v1`. Une copie locale
de l'origine H1 sur H2 impose une autozygotie synthétique à un sujet sans
reconstruire de généalogie consanguine globale. La graine est dérivée par
SHA-256 de la graine de recombinaison D15, de la méthode et du rôle de panel,
dans un domaine D16 distinct. Aucun tirage ne dépend d'un scan PLINK.

La configuration pré-micro-panel comporte quatre scénarios :

- un sujet sans ROH injectée, avec bornes nulles ;
- une ROH excluant la position cible (`TARGET_DISTINCT`) ;
- une ROH couvrant la cible, sur un sujet sans injection IBD D15 ;
- une ROH couvrant la cible, sur le premier sujet d'une paire positive D15
  explicitement nommée, afin d'évaluer la coexistence ROH/IBD.

Les trois premiers sujets sont distincts, hors pedigree, hors porteurs D14 et
hors paires ou matrices D15. Le quatrième reprend un positif D15 hors pedigree.
Toutes les copies H1 sont conservées : la vérité IBD H1 demeure valide tandis
que H2 devient localement autozygote. Les négatifs et matrices D15 ne peuvent
pas être modifiés. H2 reste inchangée hors intervalle et toutes les
transmissions sont contrôlées. Une incompatibilité d'une copie IBD déjà
injectée bloque également le run.

Les intervalles sont fournis en bp fermés, avec longueur
`end_bp - start_bp + 1`, et convertis en cM sur la carte D11a. La configuration
d'intégration fixe trois intervalles de 2 Mb ; les contrôles exigent au moins
1 500 kb, 50 sondes, au plus 50 kb par sonde et un trou interne maximal de
1 000 kb, avec des marqueurs flanquants. Une insuffisance bloque la
matérialisation sans déplacement des bornes ni abaissement des seuils.

La vérité ROH, son audit de densité et la vérité IBD conservée sont des
artefacts privés distincts, empreintés. Le cas couplé crée également un
partage H2–H1 dans la paire D15 sur l'intersection ROH/IBD. Cet intervalle
est inscrit dans une vérité IBD induite séparée et contrôlé dans le VCF ;
il ne doit être compté ni comme faux positif, ni comme répétition indépendante.
La matrice D15 reste une origine de simulation, hors des paires désignées
par le protocole d'évaluation. L'origine H1 copiée constitue une
surcharge d'ascendance sur H2, à appliquer par les futurs évaluateurs.
`TARGET_DISTINCT` signifie uniquement que la cible est hors ROH ; les
coordonnées peuvent chevaucher le segment fondateur d'une autre famille.
Les bornes injectées ne sont pas les limites maximales d'un appel PLINK, et
le témoin sans injection peut conserver une ROH de fond inconnue.

Le rapport indépendant `CAL-DEV-D16` valide la génération et la conservation
des contrats. Les performances ROH et leurs intervalles restent
`NOT_EVALUABLE` jusqu'aux appels réels et à l'évaluateur D26. Le fardeau
genome-wide, les autres classes de longueur et le panel aveugle restent à
qualifier. Aucune longueur ROH ne remplace une entrée de datation Gamma.

### D17 — Injecter les imperfections de mesure

**Action :** appliquer manque, erreurs génotypiques, densité réduite et erreurs
de phase selon des profils versionnés.
**Livrable :** registre des erreurs injectées.
**Réussite :** aucune erreur non enregistrée dans la vérité.

### D18 — Exporter les entrées aveugles

**Action :** produire les formats attendus par le pipeline sans y inclure la
vérité.
**Livrable :** `blinded_inputs/` et inventaire.
**Réussite :** validation des formats et absence de champs de vérité.

### D19 — Implémenter l'évaluateur technique

**Action :** comparer fichiers, effectifs, ordres, génotypes et empreintes.
**Livrable :** évaluateurs `00–04`.
**Réussite :** détection de toute perte, duplication ou substitution injectée.

### D20 — Implémenter l'évaluateur QC

**Action :** comparer alertes et exclusions aux défauts injectés.
**Livrable :** évaluateurs `05` et `10`.
**Réussite :** matrice de confusion et motifs complets.

### D21 — Implémenter l'évaluateur de parenté

**Action :** comparer degrés KING et ensembles indépendants à la vérité.
**Livrable :** évaluateurs `06` et `07`.
**Réussite :** métriques par degré, ascendance et niveau de manque.

### D22 — Implémenter l'évaluateur de structure

**Action :** mesurer séparation, mélanges et atypiques.
**Livrable :** évaluateurs `08` et `09`.
**Réussite :** décisions reproductibles sans fuite d'étiquette vraie.

### D23 — Implémenter l'évaluateur de phasage

**Action :** calculer erreurs de commutation, erreurs ponctuelles et phase de la
mutation.
**Livrable :** évaluateurs `10A`, `11`, `12` et `12A`.
**Réussite :** métriques stratifiées et bornes d'incertitude.

### D24 — Implémenter l'évaluateur d'haplotype fondateur

**Action :** comparer segment récupéré et segment vrai.
**Livrable :** évaluateur `13`.
**Réussite :** sensibilité, spécificité et erreurs de bornes publiées.

### D25 — Implémenter l'évaluateur de datation

**Action :** comparer âge Gamma et âge vrai sur toute la grille.
**Livrable :** évaluateur `14` et courbes de calibration.
**Réussite :** biais, RMSE, couverture et sensibilités familiales publiés.

### D26 — Implémenter les évaluateurs LD et ROH

**Action :** comparer `r²`, distances et segments ROH à la vérité.
**Livrable :** évaluateurs `15` et `16`.
**Réussite :** erreurs de détection et de longueur quantifiées.

### D27 — Implémenter l'évaluateur d'ascendance

**Action :** comparer les positions globales/locales aux mélanges simulés.
**Livrable :** évaluateur `16A`.
**Réussite :** stabilité et limites publiées sans attribution géographique.

### D28 — Implémenter l'évaluateur d'enrichissement

**Action :** mesurer calibration sous le nul et puissance sous alternatives.
**Livrable :** évaluateur `16B`.
**Réussite :** probabilités, puissance et non-évaluations vérifiées.

### D29 — Implémenter l'évaluateur IBD explicite

**Action :** comparer appels, consensus et bornes à la vérité IBD.
**Livrable :** évaluateur `16C`.
**Réussite :** fréquence de faux positifs, sensibilité et concordance avec
intervalles de confiance.

### D30 — Implémenter les évaluateurs 17 à 19

**Action :** contrôler sensibilités, figures et rapport consolidé.
**Livrable :** évaluateurs `17`, `18` et `19`.
**Réussite :** aucune statistique inventée ou recalculée par le rendu.

### D31 — Construire le moteur de rapports

**Action :** adapter le modèle `recherche.html`, le CSS linéaire et Paged.js du
skill `report-formater`.
**Livrable :** générateur de rapports français autoportants.
**Réussite :** un même JSON produit un HTML déterministe et imprimable.

### D32 — Construire chaque rapport indépendant

**Action :** implémenter les rapports listés à la section 7.
**Livrable :** un HTML, un JSON et des empreintes par rapport.
**Réussite :** aucun rapport manquant et aucune conclusion transmise seulement
par un rapport transversal.

### D33 — Exécuter le micro-panel

**Action :** lancer les tests unitaires, les outils simulés puis l'intégration
réelle identifiée.
**Livrable :** premier jeu complet de rapports techniques.
**Réussite :** zéro fuite, zéro incohérence de schéma et zéro artefact non
traçable.

Le micro-panel doit inclure au moins une fenêtre principale D14b, une fenêtre
de défi censurée et les contrôles D15 positif/négatif. Le chemin de code d'un
événement interne doit être couvert par une fixture technique déterministe,
sans filtrer les tirages du panel principal. Ce passage vérifie la chaîne et
les contrats ; son faible effectif ne permet pas d'estimer une performance
scientifique.

### D34 — Exécuter le pilote de calibration

**Action :** lancer le panel de 500 individus destiné au réglage.
**Livrable :** métriques, ressources et rapports de calibration.
**Réussite :** résultats exploitables et ressources compatibles avec la montée
en charge.

### D35 — Figer les seuils scientifiques

**Action :** décider les seuils à partir du pilote et les signer avant la
validation.
**Livrable :** version scellée de `CRITERES_ACCEPTATION.yaml`.
**Réussite :** empreinte enregistrée et fichier non modifié pendant le panel
aveugle.

### D36 — Exécuter le pilote de validation aveugle

**Action :** lancer les 500 individus distincts sans accès à leur vérité.
**Livrable :** rapports aveugles, puis évaluation après clôture.
**Réussite :** aucune modification de code ou seuil entre run et ouverture de
la vérité.

### D37 — Autoriser ou refuser la qualification à 5 000 individus

**Action :** examiner temps, mémoire, stockage, coûts et qualité.
**Livrable :** décision `GO`, `REVISE` ou `STOP`.
**Réussite :** décision motivée et signée avant toute génération volumineuse.

### D38 — Exécuter les panels de qualification

**Action :** lancer `2 500 + 2 500` avec critères scellés.
**Livrable :** rapports complets `CAL-REP-00` à `CAL-REP-19`.
**Réussite :** fin sans violation du protocole ou rapport manquant.

### D39 — Produire la validation bout en bout

**Action :** consolider uniquement les décisions des rapports indépendants.
**Livrable :** `CAL-REP-E2E`.
**Réussite :** cohérence des effectifs, statuts, versions et empreintes.

### D40 — Produire la décision de qualification

**Action :** statuer sur chaque méthode et chaque restriction d'usage.
**Livrable :** `CAL-REP-DEC`.
**Réussite :** décision humaine explicite, datée et liée aux empreintes.

### D41 — Intégrer les améliorations au pipeline

**Action :** proposer, dans des changements séparés, les corrections de biais,
contrats de provenance et améliorations de rendu validées.
**Livrable :** tests, documentation et nouvelle configuration candidate.
**Réussite :** aucune modification de méthode sans identifiant nouveau,
comparaison de référence et rapport de qualification.

### D42 — Maintenir la qualification

**Action :** définir les changements qui invalident la calibration.
**Livrable :** matrice de compatibilité des versions.
**Réussite :** recalibration obligatoire lors d'un changement de carte, panel,
outil, seuil, algorithme ou contrat scientifique.

## 10. Portes de passage

*Tableau 4 — Portes de passage et autorisations successives.*

| Porte | Autorise | Conditions minimales |
|---|---|---|
| G0 | Développement | Plan et périmètre approuvés |
| G1 | Micro-panel | Schémas, isolation, conventions, manifest D14b et tests unitaires réussis |
| G2 | Pilote calibration | Micro-panel et rapports techniques `PASS` |
| G3 | Pilote aveugle | Critères scellés et calibration documentée |
| G4 | Qualification 5 000 | Pilote aveugle, ressources et QA acceptables |
| G5 | Décision scientifique | Tous les rapports indépendants disponibles |
| G6 | Intégration pipeline | `CAL-REP-E2E` et `CAL-REP-DEC` approuvés |

Une porte échouée bloque la phase suivante. Elle ne peut pas être contournée par
une modification manuelle d'un statut.

## 11. Tests à prévoir

### 11.1 Tests unitaires

- déterminisme des graines ;
- transmissions mendéliennes ;
- conversions `bp`/`cM` ;
- injection et absence de segments ;
- séparation vérité–aveugle ;
- refus des chemins réels ;
- schémas et erreurs attendues ;
- calcul des métriques et intervalles ;
- génération déterministe des rapports.

### 11.2 Tests d'intégration simulés

- PLINK, KING, bcftools, SHAPEIT5, Hap-IBD et Refined IBD simulés ;
- erreurs, timeouts et sorties incomplètes ;
- reprise d'un run sans écrasement ;
- contrôle des signatures et producteurs.

### 11.3 Tests d'intégration réels identifiés

- petit jeu sans donnée réelle ;
- outils externes réellement exécutés ;
- dossier temporaire uniquement ;
- aucune écriture dans `data/input/` ou les runs scientifiques existants.

### 11.4 Qualification complète

- suite ciblée après chaque lot ;
- suite complète Python avant une porte majeure ;
- comparaison à un résultat de référence ;
- validation visuelle et PDF de chaque famille de rapports.

## 12. Gestion des risques

*Tableau 5 — Risques principaux et mesures de maîtrise.*

| Risque | Conséquence | Mesure de maîtrise |
|---|---|---|
| Fuite de vérité | Performance artificiellement élevée | Séparation physique et tests d'accès |
| Suradaptation | Mauvaise généralisation | Panel aveugle distinct et critères scellés |
| Simulation irréaliste | Calibration optimiste | Plusieurs modèles et imperfections injectées |
| Coût trop élevé | Qualification inachevée | Micro-panel, pilote et porte G4 |
| Réutilisation d'un résultat réel | Circularité | Interdiction de `dock6_real` et audit des chemins |
| Faux gain de précision | Intervalle trompeur | Couverture et biais obligatoires |
| Fenêtre trop courte | Censure confondue avec une recombinaison | Dimensionnement D14b en cM et série de défi séparée |
| Rapport incomplet | Décision non traçable | Contrat commun et QA `report-formater` |
| Changement d'outil | Calibration périmée | Matrice de compatibilité et recalibration |

## 13. Améliorations autorisées après qualification

Les résultats pourront conduire à :

- corriger un bug démontré par la vérité simulée ;
- ajuster un seuil uniquement selon le panel de calibration ;
- ajouter une provenance de calibration obligatoire à `16C` ;
- publier une courbe de calibration de la datation ;
- ajouter un modèle corrigé sous un nouvel identifiant ;
- améliorer les intervalles lorsqu'une couverture correcte est démontrée ;
- élargir les intervalles si la simulation montre une sous-couverture ;
- améliorer les figures et textes de prudence ;
- renforcer les tests de non-régression.

Toute modification scientifique devra être séparée d'un simple changement de
rendu et accompagnée d'un test ciblé.

## 14. Critère de clôture du projet de calibration

Le projet sera considéré comme terminé lorsque :

1. les 25 rapports d'étapes et les deux rapports transversaux existent ;
2. leurs données, figures et empreintes sont vérifiées ;
3. le panel aveugle est évalué sans modification postérieure des seuils ;
4. les limites et domaines non évaluables sont publiés ;
5. la décision de qualification est signée ;
6. une configuration candidate séparée référence le rapport et son SHA-256 ;
7. la petite suite de non-régression est intégrée aux tests ;
8. les runs et rapports générés restent hors Git ;
9. les données et runs DOCK6 existants sont restés immuables.

## 15. Références méthodologiques initiales

1. Gandolfo LC, Bahlo M, Speed TP. *Dating Rare Mutations from Small Samples
   with Dense Marker Data*. Genetics. 2014;197:1315–1327.
   DOI : `10.1534/genetics.114.164616`.
2. Zhou Y, Browning SR, Browning BL. *A Fast and Simple Method for Detecting
   Identity-by-Descent Segments in Large-Scale Data*. American Journal of Human
   Genetics. 2020;106:426–437. DOI : `10.1016/j.ajhg.2020.02.010`.
3. Hofmeister RJ et al. *Accurate rare variant phasing of whole-genome and
   whole-exome sequencing data in the UK Biobank*. Nature Genetics.
   2023;55:1243–1249. DOI : `10.1038/s41588-023-01415-w`.
4. Talts S, Betancourt M, Simpson D, Vehtari A, Gelman A. *Validating Bayesian
   Inference Algorithms with Simulation-Based Calibration*. 2018.
   `arXiv:1804.06788`.
5. Documentation interne : `docs/modules/variant_age.md`.
6. Documentation interne : `docs/modules/explicit_ibd.md`.

## 16. Prochaine action autorisée

Après approbation de ce plan, la prochaine action est limitée aux étapes D01 à
D08 : protocole, critères, arborescence, schémas et gardes d'isolation. Aucun
panel de 500 ou 2 500 individus et aucun outil externe réel ne sera lancé avant
une nouvelle revue du point de contrôle G1.
