# Continuité de session — Effet_fondateur2

**Date :** 26 novembre 2026

**Branche :** `codex/refactor-external-adapters`

**Base :** `origin/main` à `9dfc0bc`

## État

La campagne de calibration R00–R14b reste close et archivée sous
`archive/calibration/`. Elle ne fait pas partie du pipeline actif et ne doit pas
être exécutée ou restaurée dans le paquet V2.

Un audit initial du code actif a été réalisé. Le refactoring SHAPEIT5 est
terminé pour ce lot. `build_shapeit5_inputs` est désormais un orchestrateur
court qui délègue la validation des sources, la préparation des chemins,
l'export VCF, la politique mendélienne, les contrôles scientifiques, le
nettoyage et la construction du manifeste à des fonctions séparées. Des objets
internes typés transportent le plan, les chemins de staging et le résultat
validé sans modifier le contrat publié.

La centralisation progressive des commandes externes a commencé dans
`external_tools.py`. Ce module fournit le runner sans shell, la résolution des
exécutables, les sondes non bloquantes de version et la normalisation des
échecs techniques. Les modules de préparation et d'exécution SHAPEIT5
l'utilisent via de fins adaptateurs qui préservent leurs exceptions et codes
d'erreur métier historiques.

Le lot QC/kinship utilise maintenant la même couche commune dans
`qc_preliminary`, `qc_final`, `build_kinship_panel` et `infer_kinship`. Chaque
module conserve la traduction exacte de ses erreurs PLINK ou KING historiques.
Les quatre modules ne contiennent plus d'appel direct à `subprocess.run` ou
`shutil.which`.

Le lot conversion/préparation utilise également cette couche dans
`convert_acpa`, `prepare_target_variant_dataset` et `prepare_target_region`.
Les commandes PLINK et leurs codes d'erreur historiques sont inchangés. Aucun
appel direct à `subprocess.run` ou `shutil.which` ne subsiste dans ces trois
modules ; `convert_acpa` conserve seulement `shutil.rmtree` pour le nettoyage
de son spool temporaire.

Le lot analyses de population utilise maintenant la couche commune dans
`analyze_population_structure` et `analyze_reference_ancestry`. La couche sait
désormais écrire une sortie standard textuelle directement dans un nouveau
fichier et conserve le stderr ainsi que l'état de timeout dans son exception
technique. Les adaptateurs locaux reconstruisent toujours exactement les codes
d'erreur PLINK et bcftools historiques.

Le lot de phasage autosomique utilise maintenant la couche commune dans
`prepare_autosomal_phasing_panel` et `phase_autosomal_panel`. L'exception
technique conserve aussi stdout afin que les adaptateurs locaux puissent
reproduire les détails combinés stdout/stderr historiques. Les commandes,
cartes, contrôles d'échantillons et validations de phasage restent locaux aux
étapes.

Le lot d'analyses secondaires utilise maintenant la couche commune dans
`analyze_roh` avec `roh/analysis.py`, puis `analyze_local_ld` avec
`ld/local.py`. La résolution et la contrainte de version PLINK 1.9, ainsi que
les exécutions ROH/LD, sont centralisées. Les copies de sorties natives gardent
leurs usages locaux de `shutil`.

Le lot IBD et haplotype fondateur est également centralisé. Les requêtes
bcftools de `call_explicit_ibd`, la sonde Java et les lancements Hap-IBD et
Refined IBD passent par la couche commune. La sonde de version et les requêtes
bcftools de `infer_founder_haplotype`, ainsi que les requêtes de
`evaluate_founder_haplotype_enrichment`, utilisent le même adaptateur. Les
exceptions métier, codes retour, sorties capturées et journaux historiques
sont conservés.

Le lot d'infrastructure des références est terminé dans `references/window.py`
et `references/harmonization.py`. La résolution de bcftools et toutes les
commandes d'extraction, de contrôle et de normalisation passent maintenant par
la couche commune. Les runners injectables, les erreurs métier historiques et
le nettoyage des répertoires de staging sont conservés.

La centralisation des adaptateurs externes est désormais achevée dans les
derniers composants autonomes. `phasing/shapeit5.py` délègue la résolution et
les sondes de version, `phasing/execution.py` délègue la résolution résiduelle
de bcftools, et `orchestrator/environment.py` utilise les sondes communes sans
rendre l'inventaire bloquant. Les seuls appels directs restants sont les deux
frontières intentionnelles : `external_tools.py`, qui implémente la couche
commune, et `orchestrator/runner.py`, qui lance les modules Python des étapes.

Un audit séparé du code mort a ensuite couvert tous les fichiers Python actifs
du dépôt, y compris `scripts/`, les simulations et les tests, en excluant
explicitement `.venv`, `archive/`, les caches, builds et données. Les candidats
ont été vérifiés par AST, comptage lexical global et recherches ciblées. Ont été
retirés : l'exception jamais référencée `PhaseTargetRegionExternalError`, la
propriété interne jamais lue `_ReferenceVariant.canonical_key`, la constante
jamais lue `INDEPENDENT_COHORT_IDS`, le wrapper privé `_mendel_errors` utilisé
uniquement par deux assertions de test, et cinq imports inutilisés. Les tests
appellent désormais directement `_mendel_check().error_count`. Les façades de
réexport, symboles publics et usages par les scripts racine ont été conservés.

La configuration pytest limite désormais la collecte à `tests/`. Une invocation
simple de `pytest` ne collecte donc plus les tests historiques archivés, qui
produisaient 37 erreurs de collecte hors du périmètre actif.

Le premier lot prioritaire de distribution est terminé. Les schémas, les trois
catalogues de référence et les ressources du rapport sont désormais embarqués
dans le paquet et résolus avec `importlib.resources`. Les anciens chemins
`schemas/` et `config/references/` restent disponibles comme liens de
compatibilité. La wheel déclare aussi NumPy, SciPy et FPDF2, puis un test la
construit et l'installe hors du dépôt avant d'exercer ces ressources.

Le deuxième lot remplace les quatre assertions de production par des erreurs
métier explicites. Le catalogue d'ascendance refuse maintenant des actifs mal
typés, la lecture des échantillons bcftools refuse un résultat absent, la
statistique d'enrichissement exige toutes ses valeurs évaluées et la résolution
PLINK n'effectue plus qu'un seul appel à `shutil.which`.

## Décisions

- Ne modifier aucune méthode scientifique, formule, donnée ou schéma dans ce
  premier lot de refactoring.
- Conserver les codes d'erreur et le comportement observable de la préparation
  SHAPEIT5.
- Conserver une seule copie versionnée des ressources distribuées dans
  `src/effet_fondateur/resources/` ; les chemins racine ne sont que des liens de
  compatibilité.
- Migrer les commandes externes module par module en conservant une traduction
  locale des erreurs métier ; ne pas créer une exception générique visible dans
  les contrats d'étape.
- Aucun commit ni push sans demande explicite.

## Audit — priorités restantes

- Les fonctions monolithiques restantes se concentrent encore dans plusieurs
  étapes scientifiques et de préparation ; elles doivent être traitées sans
  mêler de changement scientifique.
- La centralisation progressive des outils scientifiques externes est
  terminée. Les imports `subprocess` ou `CompletedProcess` qui subsistent hors
  des deux frontières intentionnelles servent uniquement aux annotations de
  types.
- Ruff contrôle désormais les erreurs Python et le code manifestement inutilisé
  sur le paquet actif et les tests. Son premier passage a confirmé puis permis
  de retirer une variable LD, une affectation de rapport et un import de test
  inutilisés, avec 24 tests ciblés réussis avant et après.
- Mypy est configuré sur un périmètre initial explicite de six modules
  centralisés (adaptateurs externes, environnement, phasage et références).
  Étendre ce périmètre module par module plutôt que masquer globalement la dette
  de typage restante.
- La couverture lignes/branches inclut les sous-processus Python des tests
  d'intégration et exclut la copie temporaire de la wheel afin d'éviter un
  double comptage. La référence actuelle est de 77,4 %, sans seuil bloquant.
- Le ré-audit après nettoyage ne trouve plus de définition, méthode ou
  constante limitée à sa seule déclaration, ni d'import inutilisé hors façades
  de réexport explicites. Cette preuve reste statique et conservatrice face aux
  usages dynamiques ou externes ; une ligne non couverte n'est pas, à elle
  seule, une preuve de code mort.

Points positifs vérifiés : aucune cible génétique propre à l'étude n'est codée
en dur dans `src/`, et les usages aléatoires repérés utilisent des graines
explicites.

## Validation

- `git diff --check` réussi ;
- imports `effet_fondateur` et `effet_fondateur.orchestrator.pipeline` réussis ;
- compilation de `src/` et `simulation_genotype_famille/` réussie ;
- `pip check` réussi ;
- construction et installation isolée de la wheel réussies ;
- présence et utilisation des schémas, catalogues et ressources de rapport
  vérifiées depuis la wheel installée hors du dépôt ;
- 32 tests ciblés ascendance/LD/enrichissement réussis ;
- 32 tests ciblés pour les entrées, l'exécution et l'intégration SHAPEIT5,
  ainsi que le nouvel adaptateur externe, réussis ;
- validations séquentielles réussies après chaque migration : 21 tests pour le
  QC préliminaire, 18 pour le QC final, 23 pour le panel kinship et 16 pour
  l'inférence KING ;
- validations séquentielles réussies pour la conversion et la préparation : 16
  tests pour `convert_acpa`, 16 pour le jeu du variant cible et 18 pour la région
  cible ;
- validations séquentielles réussies pour les analyses de population : 16 tests
  pour la structure de population et 22 pour l'ascendance de référence ;
- validations séquentielles réussies pour le phasage autosomique : 21 tests
  après la préparation du panel et 23 après l'exécution du phasage ;
- validations séquentielles réussies pour les analyses secondaires : 16 tests
  pour ROH et 17 pour le LD local ;
- validations séquentielles réussies pour le dernier lot : 12 tests pour IBD
  explicite, 7 pour l'inférence du haplotype fondateur et 19 pour
  l'enrichissement ;
- validation transversale du lot et de l'intégration de la phase cible : 49
  tests réussis ;
- validations séquentielles réussies pour l'infrastructure des références : 10
  tests pour l'extraction de fenêtre et 10 pour l'harmonisation ; validation
  transversale avec phasage et ascendance : 65 tests réussis ;
- validations séquentielles réussies pour le dernier lot : 10 tests pour la
  sonde SHAPEIT5, 12 pour l'exécution du phasage et 19 pour l'inventaire
  d'environnement avec ses consommateurs ; validation transversale références,
  phasage et environnement : 81 tests réussis ;
- audit de code mort validé par petits lots : 16 tests pour les définitions de
  phasage, 20 pour les références et 32 pour les imports/constante restants ;
- dernier lot de code inutilisé détecté par Ruff : 24 tests ciblés réussis avant
  et après suppression ;
- `ruff check src/effet_fondateur tests` réussi ;
- `mypy` réussi sur les six modules du périmètre initial ;
- suite active complète avec couverture des sous-processus : 407 tests réussis
  en 1 759,43 s, couverture lignes/branches de 77,4 % ;
- aucun nœud `assert` dans le paquet actif ;
- `pip check`, imports du paquet et `git diff --check` réussis.

## Prochaine action

Commit/push du lot combiné uniquement sur demande explicite. Pour poursuivre
l'audit, étendre mypy et examiner les zones faiblement couvertes un module à la
fois ; confirmer chaque absence d'usage par recherche statique et tests ciblés
avant toute suppression supplémentaire.
