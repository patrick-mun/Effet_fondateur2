# Continuité de session — Effet_fondateur2

**Date :** 24 novembre 2026

**Branche :** `agent/archive-calibration`

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
- Les trois derniers appels directs aux outils externes se trouvent dans IBD
  explicite, l'inférence du haplotype fondateur et l'enrichissement.
- Aucun outil de lint, typage statique ou couverture n'est configuré dans le
  dépôt actif.

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
- aucun nœud `assert` dans le paquet actif ;
- suite active complète : 383 tests réussis en 769,81 s.

## Prochaine action

Terminer la centralisation dans `call_explicit_ibd`, puis
`infer_founder_haplotype` et `evaluate_founder_haplotype_enrichment`, un module
à la fois avec tests ciblés avant chaque migration.
