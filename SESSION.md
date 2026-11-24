# Continuité de session — Effet_fondateur2

**Date :** 15 novembre 2026

**Branche :** `agent/calibration-d24-d28`

**Décision active :** R00 à R10 sont réalisés ; R11 est arrêté à la porte
`RESOURCE_BLOCKED` pour la qualification scientifique. La campagne D01–D42 est
abandonnée et ne doit pas être restaurée. Le candidat R12 initial est conservé
mais rejeté. La release indépendante R12b a franchi G5 et R13 est terminé avec
25/25 étapes `PASS` dans une portée strictement technique. R14 est clos par un
rapport technique de 13 pages et une décision explicite de non-qualification
scientifique. R07 a ensuite été révisé en v2 sans réduire les effectifs ; la
confirmation de ressources R11b est terminée avec le verdict
`RESOURCE_READY_FOR_NEW_SCIENTIFIC_RELEASE`. La release R12c a ensuite franchi
G5. R13 est maintenant clos dans une portée amendée limitée à l'intégration
des outils. R14b est terminé par un rapport final d'intégration de 16 pages ;
la performance R07/R08 complète reste non évaluée.

## État décisionnel

- Le plan directeur est
  `calibration/PLAN_DEVELOPPEMENT_PANEL_SIMULATION.md`, version unifiée 1.0 du
  12 novembre 2026, organisée en R00–R14. Sa section 5 contient la conception
  R05 du panel maître et sa section 6 la matrice R03 et le graphe R04 complets.
- Le micro-panel R10 entièrement synthétique existe dans le workspace ignoré.
  Aucun outil scientifique externe, pipeline réel ou donnée réelle n'a été
  lancé ou utilisé. Les rendus produits par les tests sont restés dans leurs
  dossiers temporaires.
- Le logiciel actif conserve les libellés génériques des rapports et
  visualisations. R09 ajoute un noyau contractuel de calibration séparé ; les
  anciens bancs exploratoires ne sont pas des producteurs autorisés.
- Maintenir `calibrated=false`, `founder_effect_proven=false` et
  `production_promotion=false`.
- R07 v2 ferme désormais le dimensionnement scientifique, mais n'annule pas
  rétroactivement la décision R14 et n'ouvre pas G5 scientifique sans R11b.

## R01 — remise à zéro terminée

- `calibration/` ne contient plus que le plan unifié.
- 119 anciens fichiers suivis ont été supprimés : 2 fichiers à la racine,
  31 configurations, 16 protocoles, 4 références, 65 schémas et 1 modèle.
- Le fichier non suivi `calibration/AUDIT_PREVISIONNEL_D29_D42.md` a également
  été supprimé. Au total, R01 a supprimé 120 fichiers et 5 sous-répertoires de
  `calibration/`.
- Tout le contenu de `calibration_workspace/` a été supprimé définitivement :
  29 623 fichiers, 1 809 sous-répertoires et 146 entrées de premier niveau.
  Le volume alloué avant suppression était de 9 199 776 Kio, soit environ
  8,8 Gio ; le workspace occupe maintenant 0 Kio et est vide.
- Les reliquats en lecture seule de deux caches D23 ont été rendus inscriptibles
  uniquement dans `calibration_workspace/`, puis supprimés.
- Aucune archive, copie, branche, stash, tag ou réécriture d'historique n'a été
  créé. Aucun chemin sous `src/`, `tests/`, `docs/`, `config/`, `data/` ou
  `simulation_genotype_famille/` n'a été supprimé ou modifié.

## R02 — intégrité du logiciel vérifiée

- `git diff --check` : réussi.
- Imports `effet_fondateur` et `effet_fondateur.orchestrator.pipeline` :
  réussis.
- Tests ciblés de configuration, démarrage et orchestration : 20 réussis en
  6,21 s.
- Suite logicielle de référence, exécutée avec
  `--ignore=tests/calibration` : 334 tests réussis en 723,04 s.
- Les 62 modules de `tests/calibration/` ont été exclus de la référence car ils
  appartiennent à l'ancien contrat de calibration ; 41 référencent directement
  les configurations, protocoles, références, schémas ou modèles supprimés.
  Ils sont classés obsolètes par dépendance, sans restauration ni modification,
  et seront réévalués/adaptés à R09 après conception du nouveau panel.
- Aucun échec logiciel ou régression imputable au nettoyage n'a été observé.

## R03 — matrice détaillée terminée

- Les 25 identifiants contractuels (`00`–`19`, avec `10A`, `12A` et
  `16A`–`16C`) possèdent leurs entrées et formats, vérités privées, unités,
  tailles structurelles, outils, classes de ressources, erreurs attendues et
  données graphiques obligatoires.
- Les effectifs finaux et seuils ne sont pas anticipés : ils restent réservés
  respectivement à R07 et R08.
- Une fermeture producteur–sorties–consommateurs existe pour chaque étape ;
  R04 l'a développée au niveau des artefacts et a vérifié les cycles.
- Les écarts prospectifs sont consignés sans modification du logiciel :
  contrats de release absents en 00, absence d'estimateur de proportions en
  16A, sensibilités 16B/16C à raccorder en 17, et séparation des rendus
  HTML/PDF de 18–19 jusqu'au rapport final R14.
- Validation : `git diff --check`, imports V2 et contrôle automatique des
  25 lignes réussis ; 20 tests ciblés de configuration et d'orchestration
  réussis en 9,60 s.
- Aucun panel, outil scientifique, rapport HTML/PDF ou donnée réelle n'a été
  utilisé. Le code, les tests et `calibration_workspace/` restent inchangés.

## Remise en cohérence générique avant R04

- Le moteur est désormais cadré comme un logiciel générique pour toute
  variation cible autosomique configurée. DOCK6 reste uniquement un profil
  d'étude historique et le nom d'un chemin de données réelles protégé.
- `AGENTS.md`, le README, le plan unifié et les contrats documentaires 16A/16C
  ne présentent plus DOCK6 ni le chromosome 19 comme périmètre du logiciel.
  Le tutoriel ACPA sur le chromosome 19 est explicitement identifié comme un
  exemple historique paramétrable.
- Les libellés 16C et de convergence sont génériques. La couverture et le pied
  du rapport 19 affichent désormais le gène ou l'identifiant du variant résolu
  depuis la configuration, avec repli sur « variant cible ».
- Deux régressions sont couvertes : rapport construit pour `GENE_EXAMPLE` sans
  mention DOCK6, et titres génériques des domaines `EXPLICIT_IBD` et
  `POPULATION_CONVERGENCE`.
- Les anciens bancs sous `src/effet_fondateur/calibration/` ne sont ni restaurés
  ni adaptés : ils appartiennent à la campagne abandonnée et restent réservés
  à la reconstruction R09.
- Validation : `git diff --check`, imports V2 et compilation ciblée réussis ;
  20 tests ciblés réussis en 8,25 s ; suite de référence hors
  `tests/calibration/` : 335 tests réussis en 721,08 s.
- Aucun panel, outil scientifique, rapport HTML/PDF persistant ou donnée réelle
  n'a été utilisé. Les méthodes et résultats scientifiques ne sont pas
  modifiés. Maintenir `calibrated=false`, `founder_effect_proven=false` et
  `production_promotion=false`.

## R04 — graphe des dépendances terminé

- Le graphe normatif relie la génération, les vues aveugles, les 25 étapes,
  la clôture, les évaluateurs, les SVG, la synthèse et le rapport terminal R14.
- Huit racines de gouvernance, quatorze familles de vérité privées et onze
  vues aveugles possèdent maintenant un producteur et des consommateurs
  autorisés. La vérité n'est accessible qu'après clôture aux évaluateurs
  explicitement autorisés.
- Chaque étape possède une enveloppe commune, des sorties déclarées et une
  justification terminale. Les sorties dynamiques doivent être fermées dans
  leur manifest ; la découverte postérieure par glob est interdite.
- `16C` doit consommer un contrôle explicite de `16A`; `17` doit utiliser un
  registre signé de runs clos couvrant aussi `16B` et `16C`; `18` ne consomme
  que des tables d'évaluation closes; `19` assemble les 24 Markdown amont.
- Huit écarts d'implémentation sont enregistrés pour R09. Ils bloquent toute
  calibration, mais aucun code ou test n'a été modifié pendant R04.
- Contrôles réussis : `git diff --check`, imports V2, couverture statique des
  25 étapes, 14 vérités et 11 vues, contrôle d'acyclicité (`cycle=false`) et
  20 tests ciblés de configuration/orchestration réussis en 5,63 s.
- Aucun panel, outil scientifique externe, rapport HTML/PDF ou donnée réelle
  n'a été utilisé. `calibration_workspace/` reste vide.

## R05 — conception du panel maître terminée

- La conception immuable est `panel_master_generic_v1`, entièrement
  synthétique, sur les 22 autosomes `GRCh38` et sans cible codée en dur.
- Treize composants couvrent sujets, haplotypes, cibles, pedigrees, trois
  populations latentes, trois références synthétiques disjointes, événements
  fondateurs, IBD, ROH, erreurs, graines et défis techniques.
- Le catalogue fixe huit types de pedigrees, cinq rôles de cibles sur des
  autosomes distincts, six classes fondatrices, les âges 5–150 générations,
  les unités familiales 3/5/10 et les segments IBD 1–5 cM.
- Dix familles d'erreurs sont appliquées après la vérité biologique. Quinze
  domaines de graines utilisent une dérivation SHA-256 déterministe et des
  flux `PCG64DXSM` indépendants.
- La couverture statique confirme 14 familles de vérité, 11 vues et les 25
  étapes. Les effectifs, répétitions et intensités restent réservés à R07 ;
  aucune valeur de performance ou seuil R08 n'a été anticipé.
- Deux écarts sont transmis à R09 : accepter des références entièrement
  synthétiques dans le profil de calibration sans altérer le profil public de
  production, et ne réutiliser aucun ancien générateur sans démonstration
  d'adéquation au panel R05.
- Validation : `git diff --check`, imports V2, contrôles statiques R05 et 20
  tests ciblés de configuration/orchestration réussis en 6,83 s.
- R05 n'a modifié ni code ni test. Aucun panel, outil scientifique externe,
  rapport HTML/PDF ou donnée réelle n'a été utilisé ; le workspace reste vide.

## R06 — partitions et vues terminées

- Les sujets d'étude sont affectés avant génération à `CALIBRATION`,
  `VALIDATION_BLIND` ou `RESOURCE_ONLY`. L'unité indivisible réunit toute
  identité, duplicat, composante de pedigree, réplication fondatrice, événement
  IBD et base biologique appariée qui partagent une relation.
- Les trois familles de références synthétiques possèdent chacune un sous-pool
  propre à la calibration et un autre à la validation. Les contrôles négatifs
  sont des strates générées indépendamment dans ces deux partitions, jamais un
  pool partagé.
- Les onze vues ont maintenant leurs filtres, transformations, interdits et
  consommateurs exacts. Leur ordre de dérivation est déterministe ; les
  génotypes d'étude restent non phasés et aucune vue ne peut modifier la vérité
  maître.
- Les paires de sensibilité restent dans une partition et ne diffèrent que par
  un facteur signé. L'étape 17 ne consomme que des runs clos et n'accède jamais
  à la direction attendue privée.
- La vérité, les graines et les identités maîtres sont physiquement séparées
  des entrées aveugles. Le pipeline doit fonctionner sur une liste blanche de
  chemins résolus et ne peut recevoir la vérité ; l'évaluateur n'y accède
  qu'après clôture et seulement par famille autorisée.
- Trois réalisations sont réservées à R09 : isolation de processus et gardes de
  chemins, projections pseudonymisées avec références propres aux partitions,
  et registre signé avec diff canonique à facteur unique.
- Validation : `git diff --check`, imports V2, contrôle statique des onze vues
  et de l'état R06/R07 réussis ; 20 tests ciblés de configuration et
  d'orchestration réussis (18 en 5,36 s et 2 en 0,13 s).
- R06 n'a modifié ni code ni test. Aucun panel, outil scientifique externe,
  rapport HTML/PDF ou donnée réelle n'a été utilisé ; le workspace reste vide.

## R07 — dimensionnement des scénarios terminé

- Le dimensionnement distingue sujets actifs, unités indépendantes et mondes
  de scénario. Tous appartiennent à la même release ; aucune identité ne passe
  entre mondes ou partitions et les comparaisons appariées ne gonflent pas les
  effectifs.
- La qualification nominale comprend 4 000 identités, 200 trios et trois
  références synthétiques de 3 000, 1 800 et 3 000 sujets par partition. Les
  réserves scellées peuvent doubler ces nombres uniquement par la règle
  prospective estimée sur le pilote de calibration.
- Les cellules primaires utilisent 500 unités indépendantes ; les contrôles
  négatifs 3 000 unités ; les outils externes stochastiques 50 graines. Ces
  valeurs correspondent notamment à une demi-largeur Wilson maximale de
  4,4 points et à une borne nulle unilatérale d'environ 0,001.
- Les niveaux sont fixés pour cinq cibles, dix âges, 3/5/10 familles, trois
  fréquences ALT, trois densités de 600 000/300 000/100 000 variants, neuf
  états populationnels, cinq longueurs IBD, les erreurs de mesure et quatre
  classes ROH numériques.
- Les plans factoriels sont fermés par domaine. La datation utilise 500
  événements par âge avec sous-ensembles familiaux emboîtés ; l'IBD utilise
  500 paires positives par longueur et 3 000 paires négatives ; les tirages
  nuls 16B sont fixés à 99/999/9 999 selon l'échelle.
- Les plafonds sont de 100 Gio/12 h pour le micro, 1 Tio/7 jours pour le
  pilote et 8 Tio/45 jours pour la qualification, avec plafonds CPU/RAM et
  timeouts par classe. Un dépassement bloque la cellule au lieu de réduire un
  effectif ou une densité.
- Validation : `git diff --check`, imports V2 et contrôles statiques des
  effectifs, formules, fenêtres, état R07 et prochaine action réussis ; 20 tests
  ciblés de configuration et d'orchestration réussis en 5,46 s.
- R07 n'a modifié ni code ni test. Aucun panel, outil scientifique externe,
  rapport HTML/PDF ou donnée réelle n'a été utilisé ; le workspace reste vide.

## R08 — critères par étape préspecifiés

- Les 25 identifiants possèdent désormais leurs métriques primaires, unités
  indépendantes, seuils, intervalles, règles de multiplicité et domaines
  `NOT_EVALUABLE`. Un `PASS` exige toutes les métriques primaires ; une
  métrique secondaire ne peut jamais compenser un échec.
- Les contrôles déterministes `00`–`04`, `09`, `10A`, `11`, `18` et `19`
  exigent une conformité exacte lorsque la règle est binaire. Une sortie ou un
  outil manquant après préparation, un timeout ou une exclusion excessive est
  un `FAIL`, pas un non-évaluable.
- Les proportions utilisent Wilson ou la borne binomiale exacte, les mesures
  continues 10 000 bootstraps par grappe et les comparaisons déclenchées par
  plusieurs contrastes la correction de Holm. Les décisions sont calculées
  sur la précision complète avant arrondi.
- Les règles scientifiques couvrent notamment QC, parenté, stabilité PCA,
  phase régionale/autosomale, segments fondateurs, datation, LD, ROH,
  positionnement de référence, enrichissement, IBD explicite et sensibilités.
  Les proportions d'ascendance 16A restent `NOT_EVALUABLE` faute de méthode
  dédiée dans le logiciel.
- Les absences portent toujours un code et les dénominateurs attendus,
  éligibles, évaluables et exclus. Les impossibilités structurelles sont
  distinguées des défaillances du logiciel ; une précision insuffisante après
  toute la réserve R07 donne `NOT_EVALUABLE`, jamais un `PASS` forcé.
- Les critères pourront seulement être amendés avant validation par une
  nouvelle version signée, avec invalidation des résultats de calibration
  affectés. Aucun changement n'est permis après ouverture de la validation.
- Validation : `git diff --check`, imports V2 et contrôle statique des 25
  étapes, quatre statuts, seuils R08, sections du plan et prochaine action
  réussis ; 20 tests ciblés de configuration et d'orchestration réussis en
  5,35 s.
- R08 n'a modifié ni code ni test. Aucun panel, outil scientifique externe,
  rapport HTML/PDF ou donnée réelle n'a été utilisé ; le workspace reste vide.

## R09 — contrats et configurations recréés

- Neuf schémas JSON autonomes couvrent la release, le registre des sujets, la
  vérité maître, les vues, l'inventaire aveugle, les runs clos, les évaluations,
  la préparation et la décision terminale.
- Cinq configurations YAML ferment les 25 étapes, huit racines, quatorze
  vérités, onze vues, critères R08, dimensions R07 du micro-panel, politique
  d'isolation et harnais simulé des outils. La section R08 exacte est liée par
  SHA-256 au bundle contractuel.
- Le nouveau noyau `contracts`, `release_guard`, `truth_router`, `evaluation`
  et `readiness` impose les chemins résolus, la séparation des partitions, les
  pseudonymes propres à chaque partition, le diff à facteur unique, l'accès à
  la vérité après clôture seulement et l'ordre de décision R08.
- Les écarts R04/R06 sont fermés statiquement : décisions humaines signées,
  dépendance `16C` vers le contrôle `16A`, registre `17` sur seize domaines,
  sorties dynamiques exhaustives et aucun HTML/PDF aux étapes 18–19.
- Les 62 anciens modules de test, leur fixture et leur README — 64 fichiers au
  total — ont été supprimés, sans restaurer leurs contrats D. Ils sont remplacés
  par six modules R09 et 45 tests ciblés. Les anciens modules sources restent
  consultables mais sont exclus des exports et producteurs autorisés.
- Le reçu `calibration/readiness/R09_READY_TO_RUN.json` est positif seulement
  pour `R10_MICRO_PANEL_CONSTRUCTION`. Il ne certifie aucune donnée générée et
  n'autorise aucun verdict scientifique ; la vérification physique commence à
  R10.
- Validation : `git diff --check`, imports V2 et calibration, compilation du
  paquet, reçu et empreintes réussis ; 45 tests R09 réussis en 0,64 s ; 20
  tests ciblés de configuration/orchestration réussis en 5,46 s ; suite
  complète : 380 tests réussis en 711,50 s.
- Aucun panel, outil scientifique externe, donnée réelle ou rapport HTML/PDF
  n'a été utilisé. `calibration_workspace/` reste vide.

## R10 — micro-panel transversal construit

- La release `panel-master-generic-v1-r10-resource` est matérialisée sous
  `calibration_workspace/releases/micro_panel_generic_v1_r10/`. Le
  constructeur refuse l'écrasement et la cible est résolue depuis la
  configuration générique, sans cible moléculaire codée en dur.
- Elle contient 60 sujets synthétiques, six sous-références disjointes pour
  240 individus, 169 variants sur les 22 autosomes, 14 familles de vérité
  privée et 11 vues aveugles. Les formats PLINK BED/BIM/FAM et VCF gzip sont
  produits de manière déterministe en Python.
- Les 25 étapes ont traversé quatre runs clos ; `17` attend le registre signé
  des seize domaines, `18` les évaluations amont et `19` les 24 Markdown. Les
  25 évaluations et les 25 Markdown portent exclusivement sur
  `TECHNICAL_INTEGRITY_ONLY`, sans mesure de performance scientifique.
- L'isolation de processus refuse effectivement la lecture directe de la
  vérité privée par le travailleur. Les inventaires d'entrée sont fermés, les
  sorties et la release sont empreintées, et aucun lien sortant n'est admis.
- Un rejeu intégral temporaire a reproduit les 484 chemins et empreintes, avec
  l'agrégat SHA-256
  `5ff90ef966cbc7bdd29bc46414b6375a95c946b1667ee7e7cb53fae7fdecd01f` ;
  sa copie temporaire a ensuite été supprimée.
- La validation indépendante finale compte 60 sujets, 240 références, 169
  variants, 22 autosomes, 14 vérités, 11 vues, 25 étapes, 25 évaluations, 25
  Markdown et 4 runs. Elle confirme zéro outil scientifique externe et zéro
  HTML/PDF. La release occupe 2 740 224 octets alloués, très sous le plafond
  de 100 Gio.
- Validations : `git diff --check`, imports V2/calibration, reçu R09 (15
  preuves, 12 contrôles) et validation physique R10 réussis ; 56 tests de
  calibration réussis en 0,90 s ; 20 tests ciblés de configuration et
  d'orchestration réussis en 5,42 s ; suite complète : 391 tests réussis en
  680,90 s.
- R10 n'a lu aucune donnée réelle et n'a exécuté ni PLINK, ni KING, ni
  bcftools, ni Java, ni outil IBD. Le pipeline logiciel de production et ses
  méthodes scientifiques n'ont pas été modifiés pendant R10 ; seuls le
  sous-système de calibration, ses contrats, tests et documentations ont été
  étendus. Aucun rapport HTML/PDF n'a été créé.

## R11 — pilote de calibration exécuté, porte non franchie

- Le toolchain isolé R11 contient SHAPEIT5 5.1.1, OpenJDK 17.0.21.0.2,
  Hap-IBD 1.0-15Jun23.92f et Refined IBD 17Jan20.102. Les adaptateurs du dépôt
  acceptent leurs versions et empreintes. PLINK 1.9, KING 2.3.2 et bcftools
  1.21 sont ceux du `PATH` ; aucun outil système n'a été remplacé.
- La release fermée
  `pilot_panel_generic_v1_r11_calibration_attempt02` contient 500 identités
  d'étude, 500 références de phase, 300 références d'ascendance et 500
  références nulles synthétiques. Les vues PLINK emboîtées comptent exactement
  600 000, 300 000 et 100 000 variants pour 500 sujets ; la vue régionale
  compte 1 200 variants.
- Le premier essai R11 est conservé : les six commandes ont été rejetées avant
  exécution car le `/usr/bin/time` macOS n'accepte pas `-o`. L'unique nouvel
  essai autorisé a corrigé seulement ce wrapper de mesure et n'a changé aucun
  paramètre scientifique.
- Dans le second essai, les trois appels PLINK, KING et les quatre appels
  bcftools ont réussi. KING a évalué 124 359 paires et en a placé 50 au-dessus
  de 0,354. Les maxima observés sont 52 305 920 octets RSS pour PLINK,
  99 725 312 pour KING et 2 785 280 pour bcftools ; leurs temps cumulés sont
  respectivement 3,97 s, 5,84 s et 0,54 s.
- Les dix graines SHAPEIT5 ont échoué de manière identique avant phasage :
  `AC field is needed in file`. Le harnais R11 a omis les tags INFO `AC/AN`
  que la préparation du pipeline de production ajoute déjà avec bcftools. Les
  outils Hap-IBD et Refined IBD ont été vérifiés mais leur analyse scientifique
  n'a pas été lancée puisque la phase amont n'existe pas.
- Les 25 Markdown R11 sont clos : 3 statuts pilotes `PASS`, 20
  `NOT_EVALUABLE` et 2 `FAIL` (12 et 16C). La release contient 18 commandes
  mesurées, 10 échecs attendus de SHAPEIT5, zéro HTML/PDF et occupe
  318 509 056 octets alloués. Sa validation physique et ses empreintes passent.
- Validations : `git diff --check` et imports V2/calibration réussis ; 87 tests
  ciblés R09–R11, SHAPEIT5, IBD et orchestrateur réussis en 92,36 s ; suite
  complète de référence : 396 tests réussis en 746,90 s.
- Aucune donnée réelle n'a été lue. Le pipeline de production et ses méthodes
  scientifiques n'ont pas été modifiés ; R11 ajoute uniquement une
  configuration, un harnais et des tests de calibration. Les trois drapeaux de
  qualification restent faux.

## R12 — correction ciblée et candidat figé

- La release R11 `pilot_panel_generic_v1_r11_calibration_attempt02` est restée
  immuable. R12 a dérivé une nouvelle release non écrasable,
  `panel_master_generic_v1_r12_frozen_candidate`, sans ouvrir
  `VALIDATION_BLIND`.
- Les VCF synthétiques d'étude et de référence ont reçu `INFO/AC` et `INFO/AN`
  par `bcftools +fill-tags`, puis ont été réindexés. Les effectifs restent 500
  sujets et 500 références, la région reste à 1 200 variants et les critères
  R08 sont byte-identiques à R11.
- Les dix graines SHAPEIT5 5.1.1 ont toutes produit 1 200 variants pour 500
  sujets avec une concordance génotypique minimale de 1,0. Le SER commun moyen
  est 0,490925 : cette mesure ne qualifie pas la méthode et n'est pas masquée.
- Hap-IBD a appelé 178 segments sur 101 paires et retrouvé 23/25 duplicats
  désignés ; Refined IBD a appelé 122 segments sur 44 paires et retrouvé
  22/25 duplicats. Le Jaccard des paires vaut 0,435644. Les grilles de longueur
  et les contrôles requis restent incomplets.
- Les étapes affectées 12, 12A et 16C sont donc closes en `NOT_EVALUABLE`, avec
  motif `OPTIONAL_CELL_ABSENT_BY_DESIGN`. R12 a exécuté 27 commandes mesurées,
  sans échec, avec un maximum de 122 339 328 octets RSS et 1 535,64 s murales
  cumulées. La release occupe 5 599 232 octets alloués.
- La matrice `compatibility_matrix_v1.yaml`, le harnais R12, les critères, la
  configuration, les modules de production et les outils sont empreintés. Le
  reçu initial `calibration/readiness/R12_FROZEN_CANDIDATE.json` a ensuite été
  invalidé par l'audit pré-R13 ; toute dérive impose un nouveau gel explicite.
- Validations : `git diff --check`, imports V2/orchestrateur/calibration,
  validation physique et SHA-256 de la release réussis ; 65 tests de
  calibration réussis en 1,15 s ; suite complète : 400 tests réussis en
  739,67 s. Aucun HTML/PDF ni donnée réelle n'a été utilisé.
- La sonde temporaire AC/AN (856 Kio alloués) a été supprimée après succès. Le
  pipeline et les méthodes scientifiques de production n'ont pas été modifiés
  pendant R12 ; les trois drapeaux de qualification restent faux.

### Audit pré-R13 et réouverture

- L'audit a confirmé l'absence totale de `blinded_inputs/VALIDATION_BLIND` et
  de `private_truth/VALIDATION_BLIND` dans le candidat, ainsi que l'absence
  d'un générateur, exécuteur et évaluateur R13 dans les empreintes gelées.
- Les seuls producteurs existants sont spécifiques à `RESOURCE_ONLY` (R10) ou
  `CALIBRATION` (R11). Réutiliser la calibration comme validation ou écrire un
  exécuteur après gel violerait l'aveuglement.
- R11 demeure sous sa porte avec vingt étapes `NOT_EVALUABLE`. Les mesures ne
  permettent pas les projections R07 obligatoires pour 4 000 sujets, 3 000
  références de phase, 22 autosomes, 50 graines et les mondes IBD ; une
  extrapolation linéaire non vérifiée est explicitement interdite.
- Le reçu R12 est corrigé en `ready_to_run=false`, `gate_g5_passed=false` et
  `candidate_state=REJECTED_INCOMPLETE`. Le candidat reste immuable pour audit,
  aucune vérité de validation n'a été créée ou ouverte et R13 n'a pas démarré.

### R11 — reprise et mesures de ressources du 13 novembre 2026

- L'inventaire a confirmé que les 25 étapes contractuelles 00–19 sont bien
  présentes dans le pipeline V2 actif. L'insuffisance concernait le harnais
  transversal R11, pas la disparition des outils du logiciel.
- Une sonde synthétique de ressources SHAPEIT5, générique et sans donnée
  réelle, a été lancée sur 500 sujets d'étude, 500 références et 27 000
  variants. Elle a été arrêtée explicitement par l'utilisateur pendant la
  quatrième itération principale, après 13 itérations achevées sur 15.
- L'essai partiel
  `resource_probes/r11_shapeit_median_chr_attempt01` est conservé et marqué
  `INTERRUPTED_BY_USER`. Il est exclu des preuves scientifiques, des
  projections de ressources et de toute qualification ; sa sortie incomplète
  ne sera pas reprise.
- Le nouvel essai `r11_shapeit_median_chr_attempt02` a terminé les 15
  itérations SHAPEIT5 sur 500 sujets, 500 références et 26 000 variants :
  3 902,81 s murales, 572 559 360 octets de RSS maximal et 260,83 s de durée
  HMM médiane. Le BCF final contient 500 sujets et 26 000 variants.
- Deux points d'échelle supplémentaires, à densité et paramètres constants,
  ont réussi : `125+125` en 336,90 s avec 201 592 832 octets RSS, et
  `250+250` en 1 236,99 s avec 326 062 080 octets RSS. Tous les BCF et index
  ont été relus et empreintés.
- Le plan prospectif `r11_shapeit_sample_scaling_v1` fixe cinq effectifs
  égaux, un contrôle du ratio étude/référence `4:3`, un contrôle de densité,
  une borne prédictive unilatérale à 95 % et le facteur de sécurité R07 de 2.
- Le point `375+375` `attempt01` a été interrompu par une coupure de session
  après 12 itérations sur 15, sans erreur SHAPEIT5. Aucun BCF final n'existe ;
  l'essai est marqué `INTERRUPTED_SESSION` et interdit pour toute projection.
- Le nouvel essai `375+375` `attempt02` a terminé ses 15 itérations en
  2 660,57 s, avec 478 208 000 octets de RSS maximal et 181,18 s de durée HMM
  médiane. Son BCF final contient 375 sujets et 26 000 variants ; aucune sortie
  de l'essai interrompu n'a été réutilisée.
- Le cinquième point égal `187+187` a terminé en 725,80 s, avec 258 396 160
  octets de RSS maximal. Les cinq points égaux `125+125`, `187+187`,
  `250+250`, `375+375` et `500+500` sont donc tous complets, relus et
  empreintés.
- Le contrôle d'allocation `286+214`, à total constant avec `250+250`, a
  terminé en 1 434,45 s et 350 228 480 octets de RSS. Les rapports appariés
  valent 1,159629 pour le temps et 1,074116 pour le RSS : l'allocation doit
  être corrigée explicitement dans les projections.
- Le contrôle `D_MID` préspecifié avant exécution conserve strictement un
  variant commun sur deux sur les mêmes 250 sujets et 250 références. Il a
  terminé sur 13 000 variants en 607,42 s et 210 345 984 octets de RSS, soit
  0,491047 du temps et 0,645110 du RSS du contrôle pleine densité.
- La projection `r11_shapeit_qualification_projection_v1` ajuste séparément
  temps, RSS et CPU par régression puissance log–log sur les cinq points
  égaux. À 4 000 sujets + 3 000 références, la borne prédictive unilatérale à
  95 %, la correction d'allocation `4:3` et le facteur de sécurité 2 donnent
  126,70 h murales, 6,67 Gio de RSS et 121,82 heures-cœur par appel médian.
- Le plafond R2 de 72 h est dépassé. Le seul phasage commun projeté sur
  `22 autosomes × 50 graines` atteint 134 000 heures-cœur, au-dessus des
  120 000 heures-cœur de qualification avant tout coût IBD ou autre étape.
  La décision R11 est `RESOURCE_BLOCKED`; aucun crédit de parallélisation non
  mesuré, aucune réduction d'effectif et aucune baisse de densité n'a été
  appliqué.
- Aucun processus SHAPEIT5 ne demeure actif. R12 reste rejeté pour la
  qualification scientifique. Les trois drapeaux de qualification restent
  faux.

## R12b — validation aveugle technique figée, G5 franchie

- Une nouvelle release non écrasable,
  `panel_master_generic_v1_r12b_validation_blind`, a été construite sans
  réutiliser les entrées du candidat R12 rejeté. Elle utilise l'espace de noms
  indépendant `VALIDATION_BLIND_R12B_V1`, 60 identités pseudonymisées et 169
  variants synthétiques génériques.
- Les onze vues `VALIDATION_BLIND` sont fermées par un inventaire de 11
  artefacts dont tailles et SHA-256 sont vérifiés. Les quatorze familles de
  vérité sont présentes sous `private_truth/VALIDATION_BLIND`, scellées et
  routées seulement vers les évaluateurs autorisés après clôture du run.
- Les composants `r13_runner.py` et `r13_evaluator.py`, le worker technique,
  la topologie et les configurations applicables sont figés dans le manifeste
  R12b par neuf empreintes. Le runner ne copie dans son workspace que la liste
  blanche aveugle ; l'évaluateur exige d'abord un manifeste `CLOSED` valide.
- R12b ne contient aucun dossier `runs/` ou `evaluations/` : R13 n'a pas été
  lancé et aucune vérité privée n'a été ouverte. Aucun outil scientifique
  externe, aucune donnée réelle et aucun rapport HTML/PDF n'a été utilisé.
- Le reçu `calibration/readiness/R12B_G5_READY.json` porte
  `gate_g5_passed=true` et autorise uniquement
  `R13_TECHNICAL_BLIND_VALIDATION`. Il interdit toute conclusion de
  performance scientifique, d'effet fondateur ou de promotion.
- R07 est explicitement `REVISION_REQUIRED` pour les effectifs de
  qualification, l'enveloppe SHAPEIT5, le budget autosomes × graines et les
  mondes IBD. Cette révision n'empêche pas le contrôle d'intégrité technique
  R13, mais reste obligatoire avant une qualification scientifique complète.
- Validation calibration : 71 tests R09–R12b réussis, dont un dry-run R13 jetable
  de bout en bout sur une release temporaire ; compilation des trois
  nouveaux composants, audit physique de la release et empreintes réussis.

## R13 — validation aveugle technique terminée

- Le runner figé a créé le workspace ignoré
  `calibration_workspace/r13_validation_blind` et clos le run
  `RUN-R13-VALIDATION-BLIND-001` le 13 novembre 2026 à 08:02:07 UTC.
- Les 25 étapes ont produit 26 sorties déclarées et empreintées. Le manifeste
  canonique porte SHA-256
  `f04ba79119ead5078764a27586903fe28a9c0cac757de319ac08e6cb22ac3a6b`.
  Pendant le run, `truth_opened=false`, `scientific_tools_executed=false` et
  `real_data_accessed=false`.
- Après validation du manifeste `CLOSED`, l'évaluateur figé a ouvert seulement
  les familles de vérité autorisées pour chaque étape. Les 25 évaluations sont
  `PASS` avec la portée `TECHNICAL_INTEGRITY` et
  `scientific_performance_evaluated=false`.
- Le bundle contient 25 métriques, 25 évaluations, 266 fichiers et 436 910
  octets logiques. Aucun HTML/PDF n'a été créé. Les validateurs du run et des
  évaluations repassent après clôture.
- Le reçu suivi est
  `calibration/readiness/R13_TECHNICAL_VALIDATION.json`. Il maintient R07 à
  `REVISION_REQUIRED`, G6 scientifique à `false` et les trois drapeaux de
  qualification à `false`.

## R14 — rapport final technique et décision terminés

- Le rapport final autorisé à R14 est produit dans
  `calibration_workspace/r14_technical_closure_attempt02` : HTML autonome et
  PDF A4 de 13 pages, quatre figures, huit tables, équations, limites et
  inventaire de dix sources empreintées.
- La décision `R14-TECHNICAL-CLOSURE-001` porte
  `technical_integrity=PASS`, `scientific_qualification=NOT_GRANTED`,
  `founder_effect=NOT_PROVEN` et `production_promotion=NOT_AUTHORIZED`.
- Le reçu suivi `calibration/readiness/R14_TECHNICAL_CLOSURE.json` conserve les
  empreintes des quatre artefacts livrés, le contrôle éditorial et les trois
  drapeaux de qualification à `false`.
- La première exportation `calibration_workspace/r14_technical_closure`, partie
  avant la fin de la pagination Paged.js, ne contenait que deux pages. Elle est
  conservée pour audit mais rejetée et ne doit pas être livrée.
- La recette de la seconde exportation confirme 13 pages, les figures 1–4, les
  tables 1–8, zéro lien interne manquant, zéro identifiant dupliqué, zéro
  ressource externe et une revue visuelle positive. Aucun outil scientifique,
  aucune donnée réelle et aucun recalcul scientifique n'ont été introduits.
- Validation finale : 74 tests de calibration et 409 tests du dépôt réussis ;
  `git diff --check`, les imports V2 et le validateur R14 réussissent également.

## R07 v2 — révision scientifique terminée après R14

- La configuration normative est
  `calibration/configurations/scientific_dimensioning_v2.yaml` et son reçu
  déterministe `calibration/readiness/R07_SCIENTIFIC_DIMENSIONING_V2.json`.
- Les minima scientifiques sont conservés : 4 000 identités d'étude, 3 000
  références de phase, 500 unités primaires, 3 000 négatifs, 500 événements
  fondateurs par âge et 50 graines d'outil. Aucun effectif, densité ou seuil R08
  n'a été abaissé pour contourner le blocage R11.
- Le calendrier SHAPEIT5 est maintenant exact : 10 graines sur les 22
  autosomes à `D_MAX`, puis 50 graines sur trois chromosomes sentinelles pour
  chacune des trois densités, soit 640 appels. Les sensibilités `D_MID` et
  `D_LOW` sont qualifiables uniquement dans ce périmètre sentinelle déclaré.
- À partir des mesures R11, la projection vaut 53 399,586235 heures-cœur et
  66 749,482793 avec une réserve additionnelle de 25 %, sous l'allocation de
  75 000 heures-cœur. Le timeout R2 passe prospectivement à 168 h, le
  calendrier à 120 jours et la concurrence à 24 workers, soit au plus
  160,062 Gio selon la borne sûre de 6,669253 Gio par worker.
- L'IBD est fermé à 90 mondes indépendants et 540 appels outil × densité : 25
  mondes positifs, 30 négatifs, 30 de localisation secondaire et 5 composés.
  Son allocation est 30 000 heures-cœur, mais sa projection mesurée reste à
  produire ; le reçu conserve donc `scientific_g5_authorized=false`.
- R11b doit mesurer les ponts d'échelle SHAPEIT5 2 000 + 1 500 et
  4 000 + 3 000, `D_LOW`, les lots IBD de 100 paires, le stockage et la
  faisabilité à 24 workers. R07 v2 autorise ces sondes synthétiques, pas une
  campagne scientifique complète ni l'accès à des données réelles.
- Validation : reçu suivi identique au calcul déterministe, 77 tests de
  calibration et 412 tests du dépôt réussis ; `git diff --check`, compilation
  et imports V2/calibration réussissent.

## R11b — confirmation de ressources terminée

- Le workspace reprenable est
  `calibration_workspace/resource_probes/r11b_scientific_resources_v1`. Son
  dessin `RESOURCE_ONLY` est lié aux empreintes de R07 v2, de la projection
  R11, du harnais et des trois outils externes ; aucune donnée réelle ni vérité
  de validation n'est utilisée.
- Les entrées synthétiques emboîtées sont préparées : cible 4 000 + 3 000,
  pont 2 000 + 1 500, contrôle `D_LOW` 250 + 250 et lot IBD de 100 paires sous
  trois densités. L'inventaire ferme 31 artefacts.
- Les 18 sondes IBD sont terminées : trois répétitions pour chaque combinaison
  de Hap-IBD/Refined IBD et `D_MAX`/`D_MID`/`D_LOW`. Les temps muraux unitaires
  vont de 0,7285 à 1,9724 s et les RSS de 104 177 664 à 350 093 312 octets.
- La sonde SHAPEIT5 `D_LOW` est close sur 250 + 250 et 4 334 variants :
  195,453559 s murales, 195,19 s CPU et 118 890 496 octets RSS.
- SHAPEIT5 avait créé lui-même l'index `.csi`; le premier post-traitement a
  donc refusé de le recréer après un calcul réussi. L'amendement
  `R11B-HARNESS-AMENDMENT-001` corrige seulement l'idempotence et la reprise du
  reçu, sans changer les entrées ni paramètres et sans relancer les HMM.
- Le pont 2 000 + 1 500 a été interrompu proprement le 13 novembre 2026 à
  10:33 UTC, à la demande de l'utilisateur avant le transport de la machine.
  Aucun reçu de pont n'a été émis. SHAPEIT5, son enveloppe de mesure et la
  continuation initiale ont été arrêtés ; la cible 4 000 + 3 000 n'a pas
  démarré. Les journaux et le reçu `INTERRUPTED_BY_USER` sont conservés dans
  `probes/shapeit_bridge_attempt01` et exclus de toute projection.
- Le pont a repris dans `shapeit_bridge_attempt02`, sans réutiliser la mesure
  interrompue. L'amendement `R11B-HARNESS-AMENDMENT-002` ne change ni les
  entrées ni les paramètres scientifiques ; il route seulement le pont vers
  la nouvelle tentative. Le pont a terminé en 14 873,458 s avec
  1 465 393 152 octets de RSS maximal ; la cible 4 000 + 3 000 a terminé en
  28 074,519 s avec 2 538 668 032 octets de RSS maximal. Les deux sorties
  contiennent 26 000 variants et les effectifs attendus.
- La décision `R11B-SCIENTIFIC-RESOURCE-CONFIRMATION-V1`, produite le
  14 novembre 2026 à 07:41 (La Réunion), vaut
  `RESOURCE_READY_FOR_NEW_SCIENTIFIC_RELEASE`. Les cinq contrôles passent :
  17 291,871 heures-cœur SHAPEIT, 1,604 heures-cœur IBD, 122,402 Gio de RSS
  agrégé à 24 workers, 30,267 jours projetés et 29,633 Gio de stockage avec
  facteur 2.
- La session `r11b_continuation` est close. Cette décision reste une preuve de
  ressources seulement : `scientific_g5_authorized=false`, aucune performance
  scientifique n'est mesurée et les trois drapeaux de qualification restent à
  `false`.
- Le reçu suivi `calibration/readiness/R11B_RESOURCE_CONFIRMATION.json` lie
  R07 v2, la décision, les trois mesures SHAPEIT5 et les 18 mesures IBD. Il
  autorise uniquement `R12C_BUILD_AND_FREEZE_NEW_SCIENTIFIC_RELEASE` et garde
  `scientific_g5_passed=false`.
- Le contrat `scientific_validation_v1.yaml` et le préflight suivi
  `R12C_SCIENTIFIC_PREFLIGHT.json` lient les minima R07, les 25 étapes R08 et
  le reçu R11b. La projection R11b est réordonnancée sans réduction sur huit
  workers : 90,802 jours projetés, sous la limite R07 de 120 jours.
- Le runner et l'évaluateur scientifiques sont gelés en états
  `FULL_EXECUTION_READY` et `FULL_EVALUATION_READY`. Le runner exécute sans
  shell des lots de huit, conserve des checkpoints reprenables et exige le
  reçu G5 ; l'évaluateur n'ouvre la vérité qu'après clôture.
- Le manifeste `scientific_metric_rules_v1.yaml` fige maintenant, par les
  empreintes de deux fragments, 135 règles dont 21 contextuelles couvrant les
  25 étapes dans l'ordre R08. Le validateur exige `LCB95` pour un minimum,
  `UCB95` pour un maximum, `CI95` pour une équivalence et interdit d'inventer
  un intervalle sur les contrôles déterministes. Le préflight ne classe donc
  plus les règles parmi les composants manquants, mais elles n'autorisent
  encore aucune évaluation.
- Le runner construit et vérifie maintenant l'inventaire exact des 640 appels
  SHAPEIT5 (`340/150/150`) et des 90 mondes IBD déclinés en 540 appels
  outil–densité. L'évaluateur applique l'ordre décisionnel R08 aux métriques
  non arrondies, distingue échec ponctuel et imprécision, traite les contextes
  composites et agrège les 25 étapes seulement après clôture. Ces moteurs
  sont désormais autorisés uniquement par le reçu G5 R12c.
- Le socle `scientific_statistics.py` et les producteurs génériques figent les
  méthodes d'inférence des 135 règles : 65 contrôles déterministes, 39
  proportions Wilson unilatérales, 4 FPR rares Clopper–Pearson exactes et 27
  statistiques continues avec 10 000 bootstraps par grappe. La dérivation des
  graines utilise HMAC-SHA-256 vers PCG64DXSM ; Holm et la correction
  Monte-Carlo `(b+1)/(N+1)` sont implémentés. L'audit est désormais
  `METHODS_AND_STAGE_BINDINGS_FROZEN`.
- Les producteurs spécifiques couvrent les 25 étapes `00` à `19` depuis des
  preuves closes. L'audit est `FULL_STAGE_PRODUCERS_FROZEN` avec `25/25` et
  aucun producteur manquant.
- La release `panel_master_generic_v1_r12c_scientific_validation` est
  matérialisée dans `calibration_workspace/releases/` : 4 000 identités
  pseudonymisées, 14 familles de vérité privées scellées, zéro run ou
  évaluation R13. Les constructions intermédiaires ont été déplacées dans les
  sauvegardes temporaires `/tmp/effet-fondateur-r12c-*-backup.*` avant chaque
  reconstruction.
- Le reçu suivi `R12C_G5_READY.json` porte `gate_g5_passed=true`, autorise R13
  scientifique sur huit workers sans réduction de volume et garde les trois
  drapeaux de qualification à `false`.
- Le fractionnement est maintenant granulaire :
  `r13_checkpoint_ledger.py` suit séparément les 1 180 appels externes, limite
  chaque réservation à 1–8 appels et produit un JSON autoritatif ainsi que des
  tableaux TSV et Markdown avec `Fait=OUI/NON`. Les sorties terminées sont
  empreintées ; une sortie absente est reprise, et l'action `interrupt` rend
  les appels encore actifs reprenables après un arrêt de la machine.
- Le registre du futur run est initialisé sous
  `calibration_workspace/r13_scientific/RUN-R13-SCIENTIFIC-001` : 1 180
  `PENDING`, zéro `RUNNING`, zéro `COMPLETED`, zéro `FAILED`. La vérification
  exhaustive des empreintes est disponible avant chaque reprise ou clôture.
- Le préflight du premier appel `SHAPEIT-DMAX-CHR01-SEED01` a révélé que R12c
  ne contient aucun VCF/BCF, aucune carte génétique et aucun worker gelé pour
  matérialiser les shards pourtant déclarés
  `DETERMINISTIC_SHARDS_AT_R13_WORK_ITEM_START`. L'appel reste donc `PENDING`,
  tentative 0, et aucun outil externe n'a été lancé. Le reçu d'audit est
  `governance/FIRST_CHECKPOINT_PREFLIGHT.json` dans le run R13 scientifique.
  Les entrées R11b `RESOURCE_ONLY` ne doivent pas être substituées : leur plus
  proche mesure (4 000 + 3 000, 26 000 variants, 1 thread) vaut 28 074,519 s,
  mais ne constitue pas une preuve R13.
- Ce blocage est corrigé par `R12C-EXECUTABLE-REFREEZE-002`. Le matérialiseur
  `r13_scientific_materializer.py` construit à la demande les shards SHAPEIT5
  et les 90 mondes IBD sous trois densités ; le worker
  `r13_scientific_call_worker.py` exécute un seul appel sans shell, mesure ses
  ressources et ferme le checkpoint. La vérité privée est écrite hors du run
  aveugle et son chemin n'est jamais transmis au worker.
- Le modèle figé utilise 1 200 000 rangs autosomaux et des sélections emboîtées
  1/2, 1/4 et 1/12 pour obtenir exactement 600 000, 300 000 et 100 000
  variants. Il comprend trois populations Markov et 200 trios mendéliens. Les
  cartes sont synthétiques linéaires dans un cadre GRCh38 et ne doivent pas
  être interprétées comme des cartes humaines empiriques.
- Les tests réduits simulent les outils externes ; trois intégrations réelles
  ont aussi réussi : SHAPEIT5 5.1.1 sur 60 sujets + 60 références et 20 trios
  (zéro erreur mendélienne, 60/60 variants produits), puis Hap-IBD et Refined
  IBD sur 200 sujets et 6 000 variants. Avant gel, 106 tests ciblés passent ;
  après gel, 24 tests R12c/R13 passent et l'autorisation de matérialisation est
  valide. La validation finale regroupe 107 tests ciblés réussis.
- Le reçu G5 révision 2 porte `gate_g5_passed=true`, lie 1 180 appels
  matérialisables, remplace par son empreinte l'ancien reçu incomplet et garde
  les trois drapeaux de qualification à `false`.
- Le shard réel `D_MAX` du chromosome 1 est matérialisé dans la vue aveugle :
  46 524 variants, 4 000 sujets, 3 000 références et 200 trios. La vérité
  privée séparée contient 4 000 sujets et les mêmes 46 524 variants ; les
  trois effectifs ont été contrôlés directement avec `bcftools`. Le manifeste
  a été publié le 14 novembre 2026 à 22:23:46 (La Réunion).
- Le préflight du worker vaut `READY_FOR_SINGLE_CALL_EXECUTION`, lie SHAPEIT5
  5.1.1 et confirme que la vérité privée n'est pas accessible au worker. Le
  checkpoint `SHAPEIT-DMAX-CHR01-SEED01`, ordinal 1, a démarré en tentative 1
  le 14 novembre 2026 à 22:24:16 (La Réunion), graine 15052001 et un thread.
  Il a chargé 46 524 variants, 4 000 cibles, 3 000 références et 200 trios
  avec zéro erreur mendélienne, puis atteint la première itération burn-in.
- À la demande de l'utilisateur, cette tentative a été interrompue proprement
  à 18:27:09 UTC pour laisser refroidir la machine. Le reçu de tentative
  conserve `KeyboardInterrupt`, tandis que le ledger autoritatif conserve
  `INTERRUPTED_BY_USER`.
- La reprise du 15 novembre 2026 a passé le préflight, puis réservé uniquement
  le même checkpoint avec `attempts=2`. Cette tentative a exécuté trois
  burn-in complets puis atteint le burn-in 4/5, mais sa session d'exécution a
  disparu à 11:40:20 sans processus, reçu ni erreur SHAPEIT. Le ledger a été
  fermé avec `EXTERNAL_EXEC_SESSION_TERMINATED_NO_RECEIPT` : cet arrêt externe
  après 1 h 19 min 41 s n'est pas un échec scientifique.
- La tentative 3 a démarré le 15 novembre 2026 à 11:54:38 (La Réunion), dans
  la session macOS `screen` détachée
  `r13_shapeit_dmax_chr01_seed01_a03`, réadoptée par le PID 1. Elle a terminé
  avec succès à 17:53:45 : code retour 0, 21 548,10 s murales (5 h 59 min
  08 s), 21 203,21 s CPU utilisateur, 293,28 s CPU système et
  2 248 622 080 octets de RSS maximal. Le BCF final contient exactement
  4 000 échantillons et 46 524 variants ; son index CSI, le journal et les
  empreintes du reçu sont valides. Les journaux d'erreur sont vides et aucun
  processus ne reste actif. Le registre contient 1 `COMPLETED`, 1 179
  `PENDING`, zéro `RUNNING`, zéro `FAILED` et `resume_safe=true` après
  vérification des empreintes.
- Décision méthodologique du 15 novembre 2026 : R13 ne doit pas redémontrer la
  validité scientifique intrinsèque de SHAPEIT5 ni des autres outils déjà
  établie par leurs publications. La qualification locale doit démontrer leur
  utilisation correcte dans le pipeline selon les critères publiés : identité
  et version de l'outil, entrées et paramètres conformes, contrôles familiaux
  ou synthétiques appropriés, sorties intègres, métriques prévues et
  provenance auditable. Le checkpoint en cours reste une preuve d'intégration
  grandeur réelle ; aucun autre appel ne doit être réclamé avant un amendement
  explicite de R13 et du gel concerné réduisant les répétitions redondantes.
- L'amendement `R13-LIGHTWEIGHT-INTEGRATION-001` a été audité et clos le
  15 novembre 2026 à 14:37:22 UTC. Sa portée est strictement
  `R13_TOOL_INTEGRATION_VERIFICATION`. Le checkpoint
  `SHAPEIT-DMAX-CHR01-SEED01` tentative 3 est la preuve principale : SHAPEIT5
  5.1.1, code retour 0, 4 000 sujets, 3 000 références, 46 524 variants,
  200 trios, zéro erreur mendélienne et BCF final 4 000 × 46 524 vérifié avec
  bcftools. Les sorties, le reçu et le journal sont liés par SHA-256 ; le
  worker n'a pas eu accès à la vérité.
- Les contrôles positifs R12 soutiennent l'intégration Hap-IBD (23/25
  duplicats, 178 segments) et Refined IBD (22/25, 122 segments). Les 189
  fichiers indexés de cette release ont été revérifiés. Les publications de
  référence sont liées au contrat pour la validité intrinsèque des méthodes,
  sans transformer ces contrôles locaux en qualification de performance.
- Le registre historique reste inchangé à 1 `COMPLETED` et 1 179 `PENDING`.
  La vue amendée classe ces 1 179 appels comme retirés non exécutés et en
  autorise zéro nouveau. L'overlay
  `governance/r13_execution_amendment.json` bloque désormais toute nouvelle
  réservation ; le reçu suivi est
  `calibration/readiness/R13_INTEGRATION_AMENDMENT.json`. La vérité privée
  reste fermée, la performance complète reste non évaluée et les trois
  drapeaux de qualification restent faux.
- Validation précédant le premier lancement : 102 tests ciblés
  R07/R11b/R12c/R13 réussis en 11,90 s après le gel final ; imports V2 et
  `git diff --check` réussissent.
- Validation avant publication du 15 novembre 2026 : les 192 tests ciblés de
  calibration et de rendu réussissent en 25,92 s ; la suite complète réussit
  avec 516 tests en 723,23 s. Les imports V2, `git diff --check`, l'audit des
  fichiers candidats et le contrôle de motifs de secrets réussissent aussi.
- Validation de l'amendement R13 : l'audit réel, le contrôle bcftools du BCF,
  les empreintes du reçu et le verrou d'exécution réussissent. Les 180 tests de
  calibration hors fichier R09 passent en 15,98 s ; la suite calibration
  complète donne 184 succès et un échec préexistant en 16,45 s. Les imports V2
  et `git diff --check` réussissent. Le seul test R09 défaillant est antérieur à
  cet amendement : son reçu suivi conserve l'ancienne empreinte de
  `CONVENTIONS_DEVELOPPEMENT.md` alors que le fichier présent dans `HEAD` a une
  autre empreinte. Cette incohérence indépendante n'a pas été corrigée ici.

## R14b — rapport final d'intégration terminé

- Le rapport `R14B-INTEGRATION-CLOSURE-001` est produit dans
  `calibration_workspace/r14b_integration_closure` avec le modèle autonome
  `report-formater`. Il contient un HTML autonome et un PDF A4 de 16 pages,
  quatre figures, onze tables, trois liens DOI et aucun espace réservé ni
  ressource externe de rendu.
- Les graphiques montrent la hiérarchie des preuves, le temps et la mémoire
  SHAPEIT5 observés, les détections IBD avec IC95 de Wilson descriptifs, et les
  temps IBD moyens avec étendues sur trois répétitions par cellule. Les
  équations d'intégration, de Wilson, de résumé des répétitions et de la
  borne conservatrice R11b sont explicitées avec unités et hypothèses.
- La conclusion est `COMPLETED_WITH_TOOL_INTEGRATION_VERIFIED` : SHAPEIT5 est
  `VERIFIED`, Hap-IBD et Refined IBD sont `SUPPORTED`, tandis que la performance
  R07/R08 complète reste `NOT_EVALUATED`, l'effet fondateur `NOT_PROVEN` et la
  promotion `NOT_AUTHORIZED`. Aucun nouvel appel n'est autorisé.
- Le reçu suivi est
  `calibration/readiness/R14B_INTEGRATION_CLOSURE.json`. Il lie par SHA-256 la
  décision, l'inventaire des sources, le HTML, le PDF et le manifeste. Les
  drapeaux `calibrated`, `founder_effect_proven` et `production_promotion`
  restent tous à `false`.
- Le contrôle visuel a couvert la couverture, la traçabilité, les résultats
  centraux et la dernière page. Une coupure des indices 13–16 de l'inventaire
  a été corrigée avant le rendu accepté.
- Validation finale : 15 tests ciblés R13/R14b puis 186 tests de calibration
  hors fichier R09 réussissent ; `git diff --check`, les imports V2 et le
  validateur R14b réussissent. Le fichier R09 isolé conserve ses quatre succès
  et son unique échec préexistant `receipt_evidence_mismatch`, déjà documenté
  ci-dessus et indépendant de R14b.

## État Git

- R01–R02 ont été commités et poussés dans `65f1cbc`.
- R03–R14, la clarification du cadrage générique, le noyau de calibration,
  les contrats, tests et gels R12/R12b/R12c ont été commités et poussés dans
  `d02e81a` sur `agent/calibration-d24-d28`.
- `calibration_workspace/` contient la release R10, les essais et projections
  R11, le toolchain R11, le candidat rejeté R12, la release R12b et le
  workspace R13 technique clos, la release R12c scientifique, ainsi que le
  rapport R14 accepté, son premier rendu rejeté et le rapport R14b final ; il
  reste ignoré par Git.
- Aucun commit, push ou merge n'est autorisé sans nouvelle demande.

## Prochaine action

R12c exécutable, G5 révision 2, l'amendement R13 allégé et R14b sont terminés.
Ne réclamer aucun des 1 179 appels retirés : le verrou d'exécution doit rester
actif. La prochaine décision éventuelle doit porter séparément sur l'archivage
du cycle ou sur la nécessité d'une qualification de performance R07/R08
complète ; cette dernière exige un nouveau plan explicite et ne peut pas être
déduite de la seule preuve d'intégration.
Ne toucher à aucune donnée réelle sans protocole autorisé. Les
drapeaux `calibrated`, `founder_effect_proven` et `production_promotion`
restent tous à `false`.
