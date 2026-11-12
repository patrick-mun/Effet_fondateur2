# Plan unifié du panel de simulation et de la calibration

**Projet :** Effet_fondateur2

**Version :** 1.0

**Date :** 12 novembre 2026

**Statut :** plan directeur approuvé pour la remise à zéro ; aucune nouvelle
analyse scientifique n'est autorisée avant validation du panel maître

## 1. Décision et périmètre

La campagne de calibration conduite jusqu'à D28 est abandonnée comme chaîne de
qualification. Elle a utilisé des bancs utiles séparément, mais construits avec
des individus, vérités, densités, régions et paramètres qui ne composent pas
une release scientifique commune.

La nouvelle campagne partira d'un panel maître conçu à partir des besoins
complets des étapes du pipeline `00` à `19`. Toutes les vues analysées seront
dérivées de cette même release et resteront reliées aux mêmes identités,
pedigrees, graines et vérités.

La remise à zéro concerne uniquement les contrats, configurations, protocoles
et résultats de calibration. Le code sous `src/`, les tests sous `tests/`, la
documentation fonctionnelle, les configurations du logiciel et les données du
projet sont conservés.

Après la dernière étape seulement, un rapport HTML/PDF unique rassemblera les
résultats, courbes, équations, limites et modifications logicielles proposées.
Pendant la campagne, chaque étape produira uniquement un fichier Markdown de
résultats et ses données sources.

## 2. Objectif et limites scientifiques

### 2.1 Objectif

Construire une validation reproductible capable de :

1. produire une vérité génétique synthétique cohérente ;
2. dériver des observations réalistes sans révéler cette vérité ;
3. exécuter le pipeline actif dans ses conditions réelles ;
4. comparer chaque sortie à une vérité indépendante ;
5. mesurer performances, biais, incertitudes et domaines non évaluables ;
6. qualifier ou refuser séparément chaque méthode ;
7. identifier les modifications nécessaires au logiciel sans les appliquer
   silencieusement pendant la validation.

### 2.2 Ce que la calibration ne démontrera pas

- Elle ne prouvera pas à elle seule un effet fondateur réel autour de DOCK6.
- Elle ne remplacera pas une validation biologique ou externe.
- Elle ne transformera pas IBS, ROH, LD ou proximité PCA en preuve IBD.
- Elle ne garantira pas des intervalles plus étroits.
- Elle ne permettra pas de régler un seuil à partir du résultat réel DOCK6.
- Elle ne permettra pas de déclarer un outil qualifié sur une réussite
  purement technique.

Jusqu'à une décision scientifique finale explicite :

- `calibrated=false` ;
- `founder_effect_proven=false` ;
- `production_promotion=false`.

## 3. Règles non négociables

### 3.1 Isolation

- Ne jamais lire ni copier `data/input/dock6_real/` pour construire ou régler
  le panel.
- Ne pas utiliser le cache interne comme source de simulation.
- Utiliser uniquement des identifiants synthétiques.
- Conserver les données générées, vérités privées et sorties volumineuses sous
  `calibration_workspace/`, ignoré par Git.
- Refuser les chemins interdits, les liens symboliques sortant du périmètre et
  les collisions d'identifiants.

### 3.2 Indépendance

- Les sujets de calibration et de validation auront des graines, transmissions
  et composantes de pedigree distinctes.
- Les références de phasage seront disjointes des sujets évalués.
- Les seuils seront réglés uniquement sur la calibration.
- Les critères seront scellés avant l'ouverture de la validation aveugle.
- Aucun seuil ni algorithme ne sera modifié après observation de cette
  validation.

### 3.3 Séparation vérité–analyse

- `private_truth/` contiendra la vérité inaccessible au pipeline.
- `blinded_inputs/` contiendra uniquement ce que le pipeline peut observer.
- Les manifests du pipeline ne référenceront jamais `private_truth/`.
- L'évaluateur ouvrira la vérité seulement après clôture et empreinte du run.
- La personne ou le processus qui exécute l'analyse ne recevra pas les classes
  positives, négatives ni les âges vrais.

### 3.4 Reproductibilité

Chaque release enregistrera :

- toutes les graines aléatoires ;
- les sources, licences, versions et assemblages ;
- les versions de Python, PLINK, KING, bcftools, SHAPEIT5, Hap-IBD,
  Refined IBD et Java lorsqu'ils sont utilisés ;
- les cartes génétiques et leurs unités ;
- les commandes et paramètres effectifs ;
- les SHA-256 des entrées, vérités, vues et sorties ;
- les relations producteur–artefact–consommateur ;
- les écarts au protocole et les reprises.

Un dossier finalisé ne pourra pas être écrasé ou complété silencieusement.

## 4. Sources d'autorité

En cas de divergence, l'ordre suivant s'applique :

1. `AGENTS.md` pour la sécurité, les modifications et la reproductibilité ;
2. le présent plan pour l'organisation de la calibration ;
3. `docs/modules/*.md` pour les contrats scientifiques du pipeline ;
4. les schémas et configurations recréés après conception du panel ;
5. le code réellement exécuté sous `src/effet_fondateur/` ;
6. les tests sous `tests/` ;
7. les résultats clos de la nouvelle campagne.

Le comportement actuel du logiciel est un objet à évaluer, pas une vérité
scientifique. Une divergence entre documentation, code et sortie sera notée
comme anomalie avant toute proposition de correction.

## 5. Produit central : une release de panel maître

### 5.1 Définition

Une release est un ensemble immuable qui relie :

- les sujets synthétiques et leurs rôles ;
- les haplotypes d'origine ;
- la structure de population et les mélanges ;
- les pedigrees et transmissions ;
- les mutations synthétiques ;
- les recombinaisons et segments fondateurs ;
- les segments IBD et ROH indépendants du fondateur ;
- les erreurs de mesure ;
- les vérités de génotype, phase, parenté, ascendance et segments ;
- toutes les observations dérivées ;
- les partitions de calibration et de validation ;
- les paramètres et empreintes.

Une étape ne pourra pas créer son propre panel parallèle. Si un besoin manque,
la release sera refusée avant exécution et son contrat sera corrigé de manière
prospective.

### 5.2 Couches obligatoires

La release contiendra quatre couches séparées :

1. **vérité biologique simulée** : haplotypes, transmissions, phases,
   recombinaisons, segments et âges ;
2. **vérité des événements de mesure** : marqueurs supprimés, génotypes
   perturbés, données manquantes, erreurs de phase et duplications ;
3. **observations aveugles** : fichiers réellement consommables par le
   pipeline ;
4. **sorties du pipeline** : artefacts produits sans accès aux vérités.

Les évaluations constitueront une cinquième couche créée seulement après
clôture du run.

### 5.3 Partitions

Les partitions minimales seront :

- `CALIBRATION` : réglage des seuils et analyse des erreurs ;
- `VALIDATION_BLIND` : individus, familles et transmissions indépendants ;
- `PHASING_REFERENCE` : référence distincte des sujets évalués ;
- `NEGATIVE_CONTROLS` : absence connue de l'événement testé ;
- `RESOURCE_ONLY` : petit sous-ensemble utilisé uniquement pour mesurer les
  coûts techniques.

L'indépendance sera contrôlée par composante connexe de pedigree ou de parenté,
pas uniquement par identifiant :

`CALIBRATION ∩ VALIDATION_BLIND = ∅`

`(CALIBRATION ∪ VALIDATION_BLIND) ∩ PHASING_REFERENCE = ∅`

### 5.4 Vues dérivées

Les vues suivantes seront prévues avant la génération :

- génotypes denses et vérité haplotypique ;
- observation à la densité maximale réellement exploitable de la puce ;
- sous-densités contrôlées dérivées des mêmes marqueurs ;
- vue autosomale sur les 22 autosomes ;
- vues régionales autour de plusieurs cibles synthétiques ;
- vue avec pedigrees et vue sans pedigrees, toutes choses égales par ailleurs ;
- niveaux contrôlés de données manquantes et d'erreurs de génotypage ;
- phase vraie, génotypes non phasés et phase estimée ;
- scénarios fondateurs positifs et scénarios nuls ;
- segments IBD positifs et paires négatives ;
- ROH injectées et fonds sans ROH injectée ;
- vues d'ascendance globale et locale.

Chaque vue héritera des mêmes identifiants, rôles, vérités et empreintes de
release. Une comparaison de densité ou de pedigree utilisera exactement les
mêmes sujets, références et paramètres hors facteur étudié.

## 6. Besoins couverts pour chaque étape du pipeline

Le tableau suivant est le contrat initial de couverture. Il sera complété par
les formats exacts et les effectifs lors de la conception détaillée, avant
toute génération.

| Étape | Données ou vue requise | Vérité indépendante | Mesures principales |
|---|---|---|---|
| 00 — initialisation | release miniature, environnement et manifests | liste attendue des artefacts et empreintes | complétude, immutabilité, reprise |
| 01 — validation des sources | sources publiques autorisées et cas volontairement invalides | catalogue scellé, versions, assemblage, empreintes | acceptations et rejets exacts |
| 02 — registre d'échantillons | métadonnées synthétiques, doublons et discordances injectés | identité, rôle, famille et correspondances | pertes, duplications, alignement |
| 03 — conversion ACPA | génotypes synthétiques multi-échantillons | génotypes avant conversion | concordance allèles, individus et variants |
| 04 — variant cible | porteurs, non-porteurs et cas non évaluables synthétiques | génotype cible et origine de son attribution | classification et exclusions |
| 05 — QC préliminaire | vues avec manque, hétérozygotie, duplications et erreurs contrôlés | événements de mesure injectés | sensibilité, spécificité, motifs d'exclusion |
| 06 — panel de parenté | 22 autosomes, densités et LD contrôlés | ensemble de marqueurs admissibles | fidélité de sélection et indépendance |
| 07 — parenté | pedigrees, trios, degrés variés et non apparentés | coefficients et degrés issus des pedigrees | matrice de confusion, biais, intervalles |
| 08 — structure | strates et mélanges synthétiques | proportions d'ascendance simulées | positionnement, atypiques, stabilité |
| 09 — gel des cohortes | groupes synthétiques et composantes familiales | affectations et indépendances attendues | fidélité du gel, absence de fuite familiale |
| 10 — QC final | défauts stratifiés par cohorte | défauts et exclusions attendus | décisions finales et faux rejets |
| 10A — univers autosomal | mêmes individus sur 22 autosomes | variants, ordre et cartes attendus | complétude et cohérence de l'univers |
| 11 — région cible | plusieurs cibles autosomiques et cartes génétiques | positions, allèles, fenêtres et marqueurs | conservation de la cible et conversions bp/cM |
| 12 — phasage régional | génotypes non phasés, pedigrees factoriels et références | phase vraie et chromosome porteur | erreurs de commutation, erreurs ponctuelles |
| 12A — phasage autosomal | mêmes sujets, densités et références sur 22 autosomes | phase vraie autosomale | erreurs par chromosome, densité et fréquence |
| 13 — haplotype fondateur | fondateurs vrais, partages fortuits et scénarios nuls | segment fondateur et bornes vraies | sensibilité, spécificité, erreurs de bornes |
| 14 — datation | âges et familles indépendantes connus, censure contrôlée | âge, recombinaisons et longueurs vraies | biais, MAE, RMSE, couverture, stabilité |
| 15 — LD local | haplotypes vrais, parenté et densités factorielles | fréquences et LD sur vérité | erreur de r² et effet de la parenté |
| 16 — ROH | ROH de longueurs connues et scénarios sans injection | segments ROH vrais | sensibilité, faux positifs, erreurs de longueur |
| 16A — ascendance | références distinctes et mélanges connus | ascendance globale et locale | stabilité et erreur de positionnement |
| 16B — enrichissement | porteurs, non-porteurs, familles et scénarios nuls | partage fondateur vrai et classes nulles | calibration nulle, puissance, non-évaluables |
| 16C — IBD explicite | phase estimée, 22 autosomes, segments positifs et paires négatives | segments IBD et bornes vraies | sensibilité, spécificité, concordance, bornes |
| 17 — sensibilités | vues appariées ne différant que par un facteur | facteur modifié et résultat de référence | variation attribuable, stabilité inter-runs |
| 18 — visualisations | résultats clos et données graphiques empreintées | registre des valeurs et unités autorisées | fidélité, complétude, absence de recalcul |
| 19 — synthèse finale | tous les Markdown et artefacts clos | registre de décisions et provenance | fidélité des faits, limites et cohérence |

Une étape marquée `NOT_EVALUABLE` pour un besoin absent ne pourra pas être
requalifiée en `PASS`.

## 7. Grille scientifique et dimensionnement

### 7.1 Facteurs à croiser prospectivement

La conception devra couvrir au minimum :

- plusieurs positions autosomiques, sans se limiter à DOCK6 ni au chromosome 5 ;
- âges fondateurs de 5, 10, 15, 20, 30, 40, 50, 75, 100 et 150 générations ;
- 3, 5 et 10 familles indépendantes ;
- segments IBD de 1, 1,5, 2, 3 et 5 cM ;
- scénarios sans fondateur et sans segment injecté ;
- plusieurs fréquences du variant synthétique ;
- densité maximale exploitable de la puce et sous-densités contrôlées ;
- présence et absence de pedigrees ;
- plusieurs niveaux de données manquantes, erreurs de génotypage et erreurs de
  phase ;
- strates d'ascendance et mélanges connus ;
- longueurs de ROH et niveaux d'autozygotie variés ;
- censure gauche, droite, double et absence de censure.

Les plans factoriels complets ou fractionnaires seront décidés avant génération.
Un facteur ne pourra pas être ajouté après observation d'un échec sans créer
une nouvelle version de release.

### 7.2 Échelles

**Micro-panel transversal**

- environ 24 à 60 sujets ;
- au moins deux strates, un trio et plusieurs degrés de parenté ;
- un cas positif et un cas négatif pour chaque mécanisme majeur ;
- vue réduite pour les contrôles rapides et vue 22 autosomes pour
  l'intégration réelle ;
- aucune estimation de performance scientifique.

**Pilote**

- point de départ : 500 sujets de calibration et 500 sujets de validation
  distincts ;
- mesure du temps, de la mémoire, du stockage et des taux de non-évaluables ;
- estimation préliminaire des variances et ajustement prospectif du
  dimensionnement final.

**Qualification**

- point de départ : 2 500 sujets de calibration et 2 500 sujets de validation ;
- au moins 100 trios par partition si le calcul de puissance confirme ce
  minimum ;
- au moins 200 paires indépendantes par classe de parenté évaluée ;
- au moins 3 000 paires négatives indépendantes pour l'IBD par partition ;
- environ 400 à 500 répétitions positives par longueur IBD ;
- environ 500 répétitions légères par âge pour la datation.

Ces valeurs sont des minima de conception hérités des besoins scientifiques.
Elles seront confirmées ou augmentées par une justification de précision et de
puissance à R07. Elles ne pourront pas être réduites pour des raisons de
ressources sans changement explicite de portée.

### 7.3 Fenêtres génétiques

Pour un âge vrai `g`, la distance `X` au premier événement sur un bras suit le
modèle :

`X ~ exponentielle(g / 100)`, avec `X` en cM.

Pour une demi-fenêtre de largeur `W`, la probabilité théorique de censure
d'un bras vaut :

`P(censure) = exp(-g × W / 100)`.

Avec une limite préspecifiée de 0,05 par bras :

`W_min(g) = -100 × ln(0,05) / g = 299,573 / g cM`.

Les limites physiques seront obtenues par interpolation de la carte publique.
L'extrapolation, l'approximation silencieuse `1 cM = 1 Mb` et le rabotage aux
bornes chromosomiques sont interdits. Un couple cible–âge sans marge suffisante
sera `INFEASIBLE_WINDOW` avant tout tirage.

Deux séries seront séparées :

- série principale dimensionnée pour les métriques ;
- série de défi volontairement tronquée pour tester la censure.

La série de défi ne contribuera pas silencieusement au biais, à la RMSE ou à
une courbe de correction.

## 8. Contrats à créer après la conception

Les schémas ne seront recréés qu'après validation des besoins. Le socle minimal
comprendra :

1. `panel_release_manifest` : identité, version, sources, graines et couches ;
2. `subject_registry` : rôles, partitions et composantes familiales ;
3. `master_truth_manifest` : inventaire typé de toutes les vérités ;
4. `panel_view_manifest` : filtre, transformation et héritage de chaque vue ;
5. `blinded_input_inventory` : seules entrées accessibles au pipeline ;
6. `pipeline_run_manifest` : commandes, versions, sorties et clôture ;
7. `stage_evaluation_manifest` : vérité autorisée et métriques par étape ;
8. `readiness_receipt` : preuve que tous les prérequis sont satisfaits ;
9. `calibration_decision` : qualification, restrictions et compatibilité.

La porte d'exécution sera :

`READY_TO_RUN = contrat ∧ données ∧ effectifs ∧ indépendance ∧ outils ∧ aveuglement ∧ ressources ∧ résultats_attendus`.

Chaque terme devra être démontré par un artefact empreinté. Une valeur
implicite, un chemin absent ou un effectif insuffisant donnera `false`.

## 9. Résultat obligatoire de chaque étape

Chaque étape exécutée produira :

`calibration/resultats/ETAPE_<identifiant>.md`

Ce fichier contiendra :

1. identifiant, date et version de méthode ;
2. question et hypothèse testées ;
3. identifiant et empreinte de la release ;
4. vues, partitions et effectifs réellement consommés ;
5. commandes, paramètres, graines et versions ;
6. vérité ouverte par l'évaluateur ;
7. métriques, unités et intervalles d'incertitude ;
8. comparaison aux critères préspecifiés ;
9. décision `PASS`, `FAIL`, `NOT_EVALUABLE` ou `REVIEW_REQUIRED` ;
10. exclusions, anomalies et écarts au protocole ;
11. limites d'interprétation ;
12. modifications du logiciel à envisager, sans mise en œuvre implicite ;
13. chemins et SHA-256 des tables nécessaires aux courbes finales.

Le Markdown ne contiendra aucune donnée génétique individuelle ni vérité
sensible. Les tables volumineuses, données graphiques et sorties resteront
dans le workspace de la release.

Aucun rapport HTML/PDF intermédiaire ne sera créé.

## 10. Déroulement définitif de la nouvelle calibration

### R00 — Fixer le périmètre

**État :** réalisé.

Inventorier les anciens espaces, distinguer calibration et logiciel, et
valider la remise à zéro sans archive.

### R01 — Nettoyer l'ancienne campagne

Supprimer définitivement tout le contenu de `calibration/`, sauf le présent
plan, et tout le contenu de `calibration_workspace/`. Ne rien supprimer sous
`src/`, `tests/`, `docs/`, `config/`, `data/` ou
`simulation_genotype_famille/`.

**Porte :** seul ce plan subsiste sous `calibration/`, le workspace est vide
et aucun fichier logiciel n'a été touché.

### R02 — Vérifier l'intégrité du logiciel

Exécuter les imports, les tests ciblés puis la suite de référence adaptée. Les
échecs préexistants seront distingués des régressions.

**Porte :** aucune régression imputable au nettoyage.

### R03 — Compléter la matrice des besoins

Pour chaque étape de la section 6, préciser formats, unités, vérités, tailles,
outils, ressources, erreurs attendues et données graphiques.

**Porte :** chaque entrée possède un producteur et chaque sortie un
consommateur ou une justification.

### R04 — Établir le graphe des dépendances

Décrire le chemin producteur–artefact–consommateur de la génération jusqu'au
résultat final.

**Porte :** aucun artefact implicite, aucune vérité orpheline et aucune
dépendance circulaire.

### R05 — Concevoir le panel maître

Figer sujets, références, haplotypes, pedigrees, événements fondateurs, IBD,
ROH, structure, erreurs et graines.

**Porte :** tous les besoins R03–R04 sont couverts avant génération.

### R06 — Définir les partitions et les vues

Figer les partitions indépendantes et toutes les transformations dérivées du
panel maître.

**Porte :** disjonctions démontrées, vues appariées et vérité inaccessible au
pipeline.

### R07 — Dimensionner les scénarios

Justifier les effectifs, répétitions, précisions cibles, ressources et plans
factoriels.

**Porte :** dimensionnement signé avant observation d'une performance.

### R08 — Préspecifier les critères par étape

Fixer métriques, seuils, intervalles, règles de valeurs manquantes et domaines
non évaluables.

**Porte :** aucune décision ne dépend d'un résultat encore inconnu.

### R09 — Recréer les contrats et configurations

Créer seulement les schémas, configurations et gardes nécessaires au panel
unifié. Réutiliser le code existant uniquement après vérification de son
adéquation.

**Porte :** validation statique, tests ciblés et reçu `READY_TO_RUN` du
micro-panel.

### R10 — Construire le micro-panel transversal

Générer une petite release capable de traverser réellement les étapes 00–19,
y compris les 22 autosomes lorsque l'étape l'exige.

**Porte :** chaîne techniquement complète, vérités séparées, artefacts
traçables ; aucune conclusion de performance.

### R11 — Exécuter le pilote de calibration

Exécuter les outils réels sur la partition de calibration et produire les
Markdown par étape.

**Porte :** résultats exploitables, ressources mesurées, causes des
non-évaluables comprises.

### R12 — Corriger puis figer

Les défauts démontrés peuvent conduire à des propositions de correction du
logiciel. Chaque correction autorisée sera testée et obligera à reprendre les
étapes affectées. Ensuite, méthodes et seuils seront empreintés et scellés.

**Porte :** code candidat, critères et compatibilités figés avant validation.

### R13 — Exécuter la validation aveugle

Exécuter la partition indépendante sans accès à la vérité, clore le run, puis
ouvrir la vérité pour l'évaluation.

**Porte :** aucun changement de code, méthode ou seuil entre clôture et
évaluation.

### R14 — Décider la qualification et produire le rapport final

Décider séparément la qualification de chaque méthode, ses restrictions et sa
durée de validité. Produire ensuite un unique rapport HTML/PDF avec les courbes,
équations, commentaires d'étapes, limites et modifications logicielles.

**Porte :** décision humaine explicite reliée aux versions et empreintes.

## 11. Portes globales

| Porte | Autorise | Conditions |
|---|---|---|
| G0 | conception | plan unifié approuvé |
| G1 | création des contrats | besoins, dépendances, panel et vues validés |
| G2 | micro-panel | schémas, isolation et dimensionnement validés |
| G3 | pilote | micro-panel transversal techniquement complet |
| G4 | gel | pilote interprétable et corrections closes |
| G5 | validation aveugle | code, méthodes et seuils scellés |
| G6 | qualification | validation close sans violation |
| G7 | rapport et promotion éventuelle | décision humaine et traçabilité complète |

Une porte échouée bloque la suivante. Elle ne peut pas être contournée par la
modification manuelle d'un statut.

## 12. Tests et contrôles

### 12.1 Tests unitaires

- déterminisme des graines ;
- transmissions mendéliennes ;
- conversions bp/kb/cM ;
- injection et absence des événements ;
- séparation des partitions et de la vérité ;
- refus des chemins réels ;
- validation des schémas ;
- calcul des métriques et intervalles ;
- censure et cas non évaluables ;
- fidélité du résultat Markdown.

PLINK, KING, bcftools, SHAPEIT5, Hap-IBD et Refined IBD seront simulés dans les
tests unitaires.

### 12.2 Intégration

- petit jeu entièrement synthétique ;
- outils externes réels uniquement dans les tests identifiés ;
- dossiers temporaires ou workspace de calibration ;
- reprise sans écrasement ;
- contrôle des producteurs et signatures ;
- exécution sur 22 autosomes lorsque le contrat l'exige.

### 12.3 Contrôles avant une porte majeure

1. `git diff --check` ;
2. imports du package et de l'orchestrateur ;
3. tests ciblés ;
4. suite complète avant pilote, gel et validation ;
5. vérification des SHA-256 ;
6. audit de fuite de vérité ;
7. contrôle des effectifs indépendants ;
8. revue humaine des critères.

## 13. Analyse des risques

| Risque | Conséquence | Prévention |
|---|---|---|
| panel différent par étape | résultats non composables | release maître et vues dérivées |
| fuite de vérité | performances artificiellement élevées | séparation physique et clôture avant ouverture |
| suradaptation | échec en validation | partitions indépendantes et seuils scellés |
| densité confondue avec pedigree | causalité erronée | plan factoriel sur mêmes sujets |
| simulation trop simple | calibration optimiste | erreurs, structure, nulls et défis prospectifs |
| effectifs dépendants | intervalles trop étroits | unité indépendante déclarée par métrique |
| fenêtre trop courte | censure prise pour recombinaison | dimensionnement en cM et série de défi séparée |
| phase vraie donnée à l'outil | évaluation conditionnelle trompeuse | entrée non phasée et phase vraie privée |
| correction post hoc | biais de validation | version nouvelle et reprise des étapes affectées |
| ressources insuffisantes | qualification inachevée | micro-panel, pilote et porte de ressources |
| rapport final infidèle | conclusion non traçable | valeurs issues uniquement des Markdown clos |

## 14. Modifications du logiciel

Les résultats d'étape peuvent proposer :

- correction d'un bug démontré ;
- clarification d'un contrat ;
- nouveau seuil issu uniquement de la calibration ;
- nouvel identifiant de méthode pour une modification statistique ;
- élargissement d'un intervalle sous-couvrant ;
- meilleure gestion des non-évaluables ;
- amélioration de provenance ou de rendu ;
- test de non-régression supplémentaire.

Une proposition sera inscrite dans le Markdown de l'étape. Elle ne sera
appliquée que dans un lot séparé, avec test ciblé et reprise explicite des
résultats qu'elle invalide. Une amélioration de rendu ne modifiera jamais une
statistique close.

## 15. Rapport final unique

Après R14, le rapport final contiendra :

- architecture et provenance du panel ;
- séparation calibration–validation ;
- résultats de chaque étape ;
- effectifs indépendants et incertitudes ;
- courbes et tables dérivées des sources empreintées ;
- équations réellement utilisées, variables, unités et hypothèses ;
- critères, décisions et domaines non évaluables ;
- analyse des échecs et sensibilités ;
- liste des modifications logicielles réalisées ou encore proposées ;
- restrictions de qualification et conditions de recalibration.

Le rapport sera produit avec `report-formater`, contrôlé visuellement et
exporté en HTML autoportant et PDF. Il ne recalculera aucune statistique.

## 16. Clôture et invalidation

La calibration sera close uniquement si :

1. toutes les étapes requises possèdent un Markdown clos ;
2. la validation aveugle a été exécutée sans modification post hoc ;
3. les décisions sont reliées aux versions et empreintes ;
4. les limites et non-évaluables sont visibles ;
5. le rapport final est fidèle aux sources ;
6. aucune donnée réelle n'a été utilisée pour régler le système ;
7. une décision humaine autorise ou refuse explicitement la promotion.

Un changement de carte, panel, densité, outil, seuil, algorithme, format,
contrat scientifique ou domaine d'usage invalidera tout ou partie de la
qualification selon une matrice de compatibilité à produire en clôture.

## 17. Prochaine action

La prochaine action est R01 :

1. inspecter de nouveau `git status` ;
2. supprimer tout le contenu de `calibration/` sauf ce plan ;
3. supprimer tout le contenu de `calibration_workspace/` sans archive ;
4. vérifier qu'aucun fichier logiciel n'a été touché ;
5. exécuter R02 avant toute nouvelle conception ou analyse.

Aucun panel, outil externe ou rapport ne doit être lancé avant la réussite de
R02.
