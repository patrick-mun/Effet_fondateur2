# Suivi de session

Dernière mise à jour : 7 août 2026

## État du dépôt

- Dépôt GitHub : `patrick-mun/Effet_fondateur2`.
- La PR `#4` porte la livraison de l'étape `14` depuis la branche
  `feature/v2-variant-age`, au commit fonctionnel `c844823`.
- PR `#3` fusionnée dans `main` le 6 août 2026 ; elle valide les étapes `08–13`.
- PR V2 `#2` fusionnée dans `main` le 5 août 2026.
- PR `#1` fusionnée dans `main` le 4 août 2026.
- Environnement : Python 3.12.13 dans `.venv`.
- Outils disponibles et fonctionnels : PLINK 1.9, KING, bcftools, Rscript et
  packages R du projet. SHAPEIT5 5.1.1 est installé dans l'environnement Conda
  isolé `effet-fondateur-shapeit5` ; `SHAPEIT5_phase_common` et
  `SHAPEIT5_phase_rare` démarrent correctement.
  Les commandes `Gamma`/`gamma` présentes dans le `PATH` pointent actuellement
  vers un document HTML invalide et ne sont pas utilisables.

## Bilan des avancées

- Inventaire du dépôt réalisé et artefacts manifestement obsolètes supprimés.
- Environnement Python reproductible installé et documenté dans `README.md`.
- Instructions de codage, tests, maintenance et sécurité ajoutées à
  `AGENTS.md`.
- Convertisseur multi-échantillons `acpa_to_plink.py` créé :
  - lecture des exports ACPA ;
  - sélection du chromosome ;
  - annotation dbSNP hg38 ;
  - génération de PED, MAP, groupes, cas, témoins et rapports QC ;
  - refus d'écrasement sans `--force`.
- Protocole ACPA complet et préparation de `samples.tsv` documentés.
- Simulateur `simulate_unrelated_controls.py` créé pour générer des témoins
  techniques reproductibles sous HWE, avec graine, audit et limites explicites.
- Injecteur `inject_mutation.py` créé :
  - insertion triée d'un marqueur `rsMUT...` dans une copie PED/MAP ;
  - métadonnées du gène, variant, HGVS c./p., transcrit, position et allèles ;
  - génotypes obligatoirement affectés par règle explicite de groupe ou par
    individu ;
  - audit TSV, rapport JSON, empreintes SHA-256 et données pour le rapport final.
- README mis à jour avec les commandes de conversion, simulation, injection et
  validation PLINK/Mendel.
- Revue scientifique et technique du pipeline historique réalisée : les
  analyses KING sur le seul chromosome 19, Gamma, DAPC, ROH, LD et le reporting
  ne doivent pas encore être utilisés pour une interprétation biologique.
- `PIPELINE_V2_PRECODE.md` définit l'architecture cible en 20 étapes :
  orchestrateur sans logique scientifique, contrats d'artefacts, audits
  immuables, visualisations séparées, reprise sur échec, tests, migration et
  structure optimale du dépôt. Cette architecture est générique pour l'étude
  d'un variant autosomique cible ; DOCK6 sur le chromosome 19 est son premier
  profil d'étude, pas une spécialisation du pipeline.
- Premier socle V2 créé : paquet Python sous `src/`, schéma JSON strict de la
  configuration, exemples générique et DOCK6, commande `validate-config` et
  tests ciblés. Les champs moléculaires DOCK6 inconnus restent explicitement à
  `null` et les politiques interdisent l'inférence clinique des génotypes.
- Orchestrateur minimal V2 créé avec bootstrap `00`, manifest atomique, journal
  append-only, contrats d'entrées/sorties et d'audit, empreintes SHA-256,
  publication atomique, conservation des tentatives échouées et reprise après
  contrôle d'intégrité. Le profil `config/testing/synthetic.example.yaml`
  permet de le tester sans donnée génétique ni outil externe.
- Contrats TSV V2 créés pour `samples.master.tsv` et `cohorts.frozen.tsv` :
  colonnes ordonnées, valeurs canoniques, clés uniques, provenance obligatoire
  des génotypes cibles, chemins sources confinés, cohérence parentale et motifs
  d'exclusion. Les commandes `validate-samples` et `validate-cohorts` effectuent
  ces contrôles sans afficher de donnée individuelle.
- Revue de maintenabilité V2 réalisée après lecture du code : séparation du
  runner en modules de catalogue, signature, intégrité et environnement ;
  centralisation des descripteurs d'artefacts ; commentaires de décision et
  guides d'extension ajoutés dans `docs/modules/`. Les limites de concurrence,
  timeout, cache inter-run et migration de manifest restent explicitement
  documentées avant production.
- Étape V2 `01_validate_sources` implémentée et validée sur fixtures synthétiques :
  - inventaire non destructif des exports ACPA et empreintes SHA-256 ;
  - contrôle UTF-8, métadonnées ChAS, colonnes, lignes, assemblage, autosomes et
    chromosome cible ;
  - détection des noms dupliqués, sources modifiées et liens symboliques ;
  - contrats TSV versionnés pour `source_inventory.tsv` et `source_qc.tsv` ;
  - déclaration des sources externes dans `stage_inputs.json` et intégration de
    leurs empreintes à la signature de reprise ;
  - aucune lecture de génotype à des fins d'interprétation et aucune donnée
    génétique dérivée produite.
- Étape V2 `02_build_sample_registry` implémentée et validée sur fixtures
  synthétiques :
  - consommation explicite de `source_inventory.tsv`, `source_qc.tsv` et de la
    table utilisateur déclarée par `inputs.sample_metadata` ;
  - validation du contrat `samples.master.tsv`, des identifiants, filiations,
    sources uniques et provenances de génotypes cibles ;
  - refus explicite des génotypes provenant du statut clinique ou du groupe ;
  - pseudonymes déterministes et table `sample_registry_review.tsv` sensible ;
  - approbation humaine versionnée facultative dans la configuration ; sans
    approbation, le run conserve `sample_registry_approval` comme décision en
    attente ;
  - blocage de la reprise si la table utilisateur ou un artefact dépendant a été
    modifié ; aucune nouvelle dépendance logicielle ajoutée.
- Étape V2 `03_convert_acpa` implémentée et validée sur fixtures synthétiques et
  avec un smoke test PLINK 1.9 réel dans un dossier temporaire :
  - blocage tant que `sample_registry_approval` reste en attente ;
  - lecture unique de chaque export ACPA référencé, puis alignement depuis des
    spools temporaires supprimés après validation ;
  - production distincte de `genomewide_base.bed/.bim/.fam` sur les 22
    autosomes et de `target_chromosome_base.bed/.bim/.fam` ;
  - politiques `intersection` et `union`, génotypes invalides convertis en
    manquants, conflits de coordonnées et variants multialléliques exclus ;
  - rsID absents ou conflictuels audités sans exclusion automatique ;
  - descripteurs de jeux PLINK, alignement des échantillons, audit des variants,
    QC par échantillon et chromosome, empreintes et rapport de conversion ;
  - échecs, timeouts ou sorties PLINK invalides rapportés avec le code V2 `3` ;
  - aucune nouvelle dépendance : PLINK était déjà requis et disponible.
- Étape V2 `04_prepare_target_variant_dataset` implémentée et validée sur
  fixtures synthétiques et avec PLINK 1.9 réel dans un dossier temporaire :
  - dépendances directes explicites vers le registre approuvé de l'étape `02`
    et le jeu chromosome cible de l'étape `03` ;
  - contrats versionnés pour la définition moléculaire confirmée, les génotypes
    individuels, leur audit et le résumé Mendel ;
  - couverture exacte de tous les échantillons, cohérence avec le registre et
    refus des sources dérivées du statut clinique ou du groupe ;
  - injection dans une copie, contrôle de la coordonnée et des allèles, puis
    comparaison des empreintes BED/BIM/FAM avant et après ;
  - contrôle PLINK Mendel limité au variant cible, bloquant avec le code V2 `4`
    en cas d'erreur et noté `NOT_APPLICABLE` lorsqu'aucun trio n'est évaluable ;
  - échec externe PLINK rapporté avec le code V2 `3` et tentative non publiée
    conservée pour diagnostic ;
  - aucune nouvelle dépendance : PyYAML, jsonschema et PLINK étaient déjà requis.
- Audit qualité du groupe `00–04` réalisé : architecture, contrats, sécurité des
  génotypes, codes d'échec, conservation des tentatives et contrôle scientifique
  sont cohérents. Limites conservées : confirmations moléculaires externes
  obligatoires et statut Mendel non applicable sans trio complet.
- Étape V2 `05_qc_preliminary` implémentée et validée sur fixtures synthétiques
  et avec PLINK 1.9 réel dans un dossier temporaire :
  - dépendances directes vers le registre approuvé de l'étape `02` et le jeu
    genome-wide de l'étape `03` ;
  - métriques PLINK individuelles, variants, fréquences et hétérozygotie, plus
    missingness stratifié par lot ACPA ;
  - exclusions automatiques limitées aux seuils techniques de données
    manquantes, matérialisées par des listes `--remove` et `--exclude` ;
  - MAF faible, hétérozygotie extrême, différentiel entre lots et duplicats
    potentiels conservés comme alertes sans exclusion ;
  - absence volontaire de HWE avant le gel des cohortes et discordance de sexe
    notée `NOT_EVALUATED` sur le jeu strictement autosomique ;
  - scan préliminaire des duplicats conditionné par un minimum de variants
    polymorphes, la décision finale restant réservée à l'étape `07` ;
  - aucune nouvelle dépendance : PLINK, PyYAML et jsonschema étaient déjà requis.
- Étape V2 `06_build_kinship_panel` implémentée et validée sur fixtures
  synthétiques et avec PLINK 1.9 réel dans un dossier temporaire :
  - consommation exclusive du jeu `genomewide_pre_qc` et des métriques variants
    de l'étape `05`, avec contrôle des empreintes et des ensembles ;
  - recalcul de la MAF après les exclusions individuelles techniques, sans
    appliquer aveuglément les alertes préliminaires ;
  - exclusion auditée des faibles MAF et des régions complexes explicitement
    configurées, puis pruning PLINK `--indep-pairwise` paramétré ;
  - publication d'un BED/BIM/FAM distinct pour KING et la structure, sans
    modification de l'ordre ni de l'ensemble des individus ;
  - audit de chaque variant, couverture des 22 autosomes, distances entre
    marqueurs et distribution agrégée du LD résiduel ;
  - blocage avec le code V2 `4` si le panel est trop petit, déséquilibré ou ne
    couvre pas les autosomes requis ; échec PLINK rapporté avec le code `3` ;
  - clarification que le nombre effectif opérationnel est le nombre de marqueurs
    après pruning, pas une dimension effective statistique ;
  - aucune nouvelle dépendance : PLINK, PyYAML et jsonschema étaient déjà requis.
- Étape V2 `07_infer_kinship` implémentée et validée sur fixtures synthétiques
  et avec KING 2.3.2 réel dans un dossier temporaire :
  - consommation explicite du panel indépendant de l'étape `06`, du registre
    maître et des métriques individuelles de l'étape `05` ;
  - exécution KING `--kinship` jusqu'au degré configuré et normalisation des
    sorties intra-famille `.kin` et inter-familles `.kin0` ;
  - seuils versionnés pour duplicats/jumeaux et premier à quatrième degrés,
    avec IBS0 pour distinguer une relation compatible parent-enfant ;
  - conservation des parents déclarés absents et des paires KING manquantes
    comme statuts `NOT_EVALUATED`, sans les ignorer silencieusement ;
  - séparation des relations concordantes, discordantes, familiales non
    résolues et cryptiques ;
  - proposition déterministe d'un ensemble indépendant maximal orienté qualité,
    sans prétendre garantir la cardinalité maximum ;
  - provenance de chaque exclusion proposée : KING, pedigree déclaré ou les
    deux, sans aucune application automatique ;
  - ajout de la décision `kinship_exclusion_approval` au manifest dès qu'une
    exclusion ou une anomalie de pedigree requiert une revue humaine ;
  - aucune nouvelle dépendance : KING, PyYAML et jsonschema étaient déjà requis.
- Étape V2 `08_analyze_population_structure` implémentée et validée sur fixtures
  synthétiques et avec PLINK 1.9/KING 2.3.2 réels dans un dossier temporaire :
  - consommation du panel `06`, du registre maître et de la proposition
    indépendante `07`, avec contrôle des empreintes et ensembles d'individus ;
  - export temporaire des dosages A1 par PLINK, sans publication du `.raw` ;
  - estimation des fréquences, centrage et réduction exclusivement sur les
    individus proposés indépendants, avec imputation des manquants à `2p` ;
  - exclusion auditée du modèle des variants monomorphes ou insuffisamment
    appelés dans la référence ;
  - PCA par SVD avec signe déterministe et projection séparée de tous les
    apparentés, qui ne peuvent pas modifier les axes ni les scores de référence ;
  - publication des scores, valeurs propres, métriques exploratoires d'outlier,
    fréquences, échelles et loadings nécessaires à la reproduction ;
  - distinction explicite entre l'ensemble complet projeté et l'identifiant de
    la référence indépendante ayant servi à ajuster le modèle ;
  - approximation du chi-deux réservée au signalement exploratoire, sans
    exclusion automatique ni interprétation clinique ou d'ascendance ;
  - ajout de `population_structure_exclusion_approval` au manifest uniquement
    lorsqu'un outlier est proposé pour revue humaine ;
  - aucune référence populationnelle externe et aucune DAPC dans cette version.
- Étape V2 `09_freeze_cohorts` implémentée et validée sur fixtures synthétiques
  et sur une chaîne temporaire `01` à `09` avec PLINK 1.9/KING 2.3.2 réels :
  - consommation explicite du registre, des génotypes cibles acceptés, du QC,
    de la proposition indépendante et des résultats populationnels ;
  - publication d'une matrice complète de sept cohortes avec une ligne par
    individu, y compris pour chaque exclusion et sa provenance ;
  - classification porteur/non-porteur fondée exclusivement sur le génotype
    individuel accepté de l'étape `04`, jamais sur le statut ou le groupe ;
  - groupes témoins déclarés explicitement par configuration, sans utilisation
    comme source de génotype ;
  - représentants indépendants déterministes par famille selon le QC technique,
    avec unités explicites pour contrôles, familles porteuses et composants KING ;
  - approbations inter-run liées au SHA-256 exact des propositions `07`/`08`,
    avec liste complète des exclusions et date horodatée ;
  - résolution des décisions manuelles par l'orchestrateur seulement après
    validation et publication atomique de l'étape ; un échec les conserve ;
  - identifiants approuvés confinés à `cohort_decisions.tsv`, classé sensible,
    tandis que l'audit et le rapport ne publient que des comptes ;
  - génération de sept fichiers PLINK `--keep`, dont les cohortes optionnelles
    peuvent être vides sans ambiguïté.
- Étape V2 `10_qc_final` implémentée et validée sur fixtures synthétiques et sur
  une chaîne temporaire `01` à `10` avec PLINK 1.9/KING 2.3.2 réels :
  - publication séparée de jeux finaux pour les témoins indépendants, les
    porteurs indépendants et tous les individus du chromosome cible ;
  - exclusion individuelle sur missingness avant recalcul des métriques variant ;
  - HWE et MAF appliqués uniquement aux témoins, sans HWE ni filtre MAF chez
    les porteurs indépendants ;
  - MAF locale et missingness variant sur le chromosome cible, avec conservation
    forcée et auditée du variant cible lorsqu'un seuil standard échoue ;
  - alertes de lot descriptives sans exclusion automatique ;
  - nouveaux identifiants déterministes d'échantillons et variants pour chaque
    jeu BED/BIM/FAM final, avec descripteurs et empreintes.
- Étape V2 `11_prepare_target_region` implémentée et validée sur fixtures
  synthétiques et sur une chaîne temporaire `01` à `11` avec PLINK 1.9/KING
  2.3.2 réels :
  - extraction PLINK d'une fenêtre physique gauche/droite configurable depuis
    le jeu local final de l'étape `10`, sans modifier les individus ;
  - carte TSV à identifiant unique, assemblage explicite, positions bp
    strictement croissantes et positions cM monotones ;
  - interpolation décimale uniquement entre deux ancres dont l'écart respecte
    le maximum configuré, sans extrapolation ni approximation `1 Mb = 1 cM` ;
  - contrôle des coordonnées, allèles, doublons, ordre et présence du variant
    cible avant publication ;
  - écriture des cM validés dans le BIM régional et publication de
    `target_genetic_map.tsv` ;
  - manifest générique `PLINK_BED_BIM_FAM_WITH_CM`, sans imposer prématurément
    l'adaptateur de phasage de l'étape `12`.
- Sous-étape V2 `12.0` implémentée et validée sans donnée génétique :
  - contrat `shapeit5_phase_common_rare_v1` figé sur SHAPEIT5 `5.1.1`, licence
    MIT, autosomes GRCh38, sorties BCF et carte génétique obligatoire ;
  - `phase_common` retenu pour le scaffold partagé avec la référence et
    `phase_rare` pour le variant cible rare lorsqu'il est absent du panel ;
  - pedigree enfant-père-mère, graine, threads et taille efficace explicitement
    représentés dans les commandes reproductibles ;
  - schéma de configuration strict, sonde de version exacte pour les deux
    exécutables et inventaire technique du run adaptés au contrat ;
  - aucune référence téléchargée et aucun phasage exécuté à cette sous-étape.
- Sous-étape V2 `12.1` implémentée et validée sans téléchargement de panel :
  - catalogue JSON versionné et schéma strict pour le panel haute couverture
    phasé de `3 202` échantillons 1000 Genomes sur GRCh38 ;
  - couverture explicite des 22 autosomes, de tous les échantillons et de toutes
    les populations disponibles, sans sélection implicite ;
  - résolution hors ligne par panel, assemblage et chromosome, limitée au
    domaine officiel `ftp.1000genomes.ebi.ac.uk` ;
  - URL et SHA-256 du README et du manifeste épinglés, avec MD5 fournisseur du
    VCF et de son index pour chaque autosome ;
  - aucun VCF de référence téléchargé et aucune donnée privée transmise.
- Sous-étape V2 `12.2` implémentée et validée sans téléchargement du VCF 1000G :
  - cache partagé adressé par catalogue, panel, release, assemblage, chromosome,
    URL et empreintes attendues ;
  - verrou `flock`, écriture temporaire sur le même système de fichiers et
    publication par renommage atomique ;
  - MD5 fournisseur vérifié pour VCF/index, SHA-256 épinglé pour README/manifeste
    et SHA-256 local recalculé pour les quatre fichiers ;
  - entrées publiées en lecture seule et intégralement revérifiées avant chaque
    réutilisation, sans remplacement automatique d'une corruption ;
  - reprise des seuls dossiers temporaires interrompus et mode hors ligne
    bloquant un cache absent sans appel réseau.
- Sous-étape V2 `12.3` implémentée et validée sans téléchargement du VCF 1000G :
  - extraction régionale bgzip avec `bcftools` et création d'un index tabix ;
  - conservation obligatoire de tous les échantillons et de leur ordre source,
    sans sélection de population ;
  - contrôle des longueurs canoniques des contigs GRCh38, du champ `GT`, de
    l'effectif catalogue, du nombre de variants et de l'absence de génotypes
    appelés non phasés ;
  - revalidation des empreintes du cache avant lecture et publication atomique
    d'un manifeste liant la fenêtre à la release et aux fichiers sources.
- Sous-étape V2 `12.4` implémentée et validée sans lancement de SHAPEIT5 :
  - décomposition biallélique des SNV/indels de la fenêtre de référence ;
  - jointure des marqueurs ACPA par assemblage, chromosome, position, REF et ALT,
    sans utiliser rsID ni identifiant de sonde comme clé scientifique ;
  - audit des correspondances directes, inversions A1/A2, variants propres à
    l'étude, discordances alléliques et collisions de sondes ;
  - absence de correction automatique de brin et blocage des ambiguïtés touchant
    le variant cible ;
  - conservation autorisée du variant cible absent du panel pour `phase_rare` ;
  - publication atomique d'une référence harmonisée, de son index, de la table
    d'audit et d'un manifeste lié par SHA-256 aux entrées `11` et `12.3`.
  - intégration à l'étape orchestrée `12_phase_target_region`, avec artefacts
    `12.3` et `12.4` séparés, audit V2, codes de retour et reprise stricte ;
  - décision de ne pas ajouter de FASTA : les représentations GRCh38 doivent
    déjà être minimales et toute évolution de cette méthode sera distincte ;
  - complétion autorisée d'un allèle PLINK `0` uniquement depuis une
    correspondance de référence unique, sans inférer de génotype individuel.
- Sous-étape V2 `12.5` implémentée et validée sans lancement de SHAPEIT5 :
  - sélection auditée des variants communs et conservation obligatoire du
    variant cible comme `COMMON_TARGET` ou `RARE_TARGET` ;
  - export PLINK du VCF d'étude avec les `SAMPLE_ID` maître, REF canonique imposé
    par `--a2-allele`, puis revalidation bcftools de l'ordre et des REF/ALT ;
  - carte SHAPEIT `pos/chr/cM`, pedigree trio/duo sans en-tête avec `NA` pour un
    parent absent et fichier vide autorisé lorsqu'aucun parent n'est présent ;
  - VCF bgzip/indexé, tables de sélection et de correspondance des individus,
    manifeste lié aux empreintes et publication atomique ;
  - intégration orchestrée, audit, codes de retour, reprise stricte et délai
    externe configurable, sans exécuter `phase_common` ni `phase_rare`.
- Sous-étape V2 `12.6` implémentée et validée :
  - ajout des tags `AC/AN` requis par `phase_common` dans le VCF étude `12.5` ;
  - correction du contrat `phase_rare 5.1.1`, qui ne prend pas `--log`, avec
    capture stdout/stderr dans un journal dédié ;
  - exécution séquentielle `phase_common`/`phase_rare`, BCF/index CSI, versions,
    graine, threads, `Ne`, régions et délais enregistrés ;
  - contrôle des empreintes, individus, variants, génotypes et erreurs
    mendéliennes avant/après, puis concordance avec le génotype cible explicite ;
  - attribution `H1/H2/BOTH`, transmissions trio/duo, score PP des singletons,
    zones non fiables et revue manuelle lorsque la confiance manque ;
  - publication atomique des sorties, tables versionnées, manifeste, audit,
    codes d'erreur et reprise stricte.
- Sous-étape V2 `12.7` implémentée et validée :
  - consolidation de douze contrôles couvrant intégrité, `AC/AN`, variants,
    individus, génotypes, cible, Mendel, transmissions et confiance ;
  - publication d'un résumé agrégé `NONE/H1/H2/BOTH` sans identifiant
    individuel et d'un manifeste final lié aux sorties `12.5–12.6` ;
  - warning de confiance conservé sans exclusion automatique, avec revue
    manuelle obligatoire lorsqu'un porteur reste non fiable ;
  - publication atomique, reprise stricte et déclaration des trois
    visualisations pseudonymisées attendues pour l'étape `18` ;
  - clôture de l'étape `12` après fixtures synthétiques, suite complète et smoke
    réel temporaire.

## Jeux de données actuellement disponibles

### Témoins ACPA réels

- 71 exports témoins ont été analysés avec KING sur les 22 autosomes.
- Après QC et pruning LD, 16 233 marqueurs ont été utilisés pour l'analyse de
  parenté.
- Aucun apparentement du premier au troisième degré n'a été détecté ; 8 paires
  impliquant 11 individus étaient compatibles avec un quatrième degré.
- Une couverture minimale du réseau a conduit à exclure 5 témoins, en priorité
  parmi ceux ayant le plus de génotypes manquants, afin de conserver la cohorte
  indépendante maximale.
- 66 témoins restent dans
  `data/input/complex_simulation/acpa_samples/controls/`.
- Les 5 témoins exclus ont été déplacés sans suppression dans
  `data/input/complex_simulation/acpa_samples/controls_exclus_king_degree4/`.
- Une seconde analyse KING des 66 témoins ne détecte plus aucune paire jusqu'au
  quatrième degré.
- Le jeu chromosome 19 combinant les 9 échantillons initiaux et ces 66 témoins
  réels n'a pas encore été construit.

### Conversion ACPA réelle

- Dossier conseillé pour la dernière relance : `data/output/acpa_chr19_relance/`.
- 9 individus réels.
- 2 725 marqueurs du chromosome 19.
- Aucun rsID non résolu dans le rapport de conversion.
- PED/MAP et fichiers `groupes.txt`, `cas.txt`, `temoins.txt` présents.

### Jeu de test avec témoins simulés

- Dossier : `data/output/acpa_chr19_avec_temoins_simules/`.
- 9 individus ACPA et 50 témoins synthétiques, soit 59 individus.
- 2 725 marqueurs avant injection de la mutation.
- Groupes : 5 `ATTEINT`, 3 `HTZ`, 1 `SAINS`, 50 `TEMOIN_SIMULE`.
- Graine de simulation : `20260804`.
- Réservé aux tests techniques et exploratoires, sans interprétation biologique
  du LD, des haplotypes ou des ROH.

### Entrée actuellement lue par le pipeline

- `data/input/complex_simulation/genotype_data.ped` contient 127 individus.
- `data/input/complex_simulation/genotype_data.map` contient 2 726 marqueurs,
  dont `rsMUT_11237780`.
- Ce jeu est historique et ne correspond pas au nouveau jeu ACPA de 9 individus
  avec 50 témoins simulés.
- `cas.txt` et `temoins.txt` sont absents de ce dossier, donc le pipeline actuel
  ne doit pas encore être lancé sur cette entrée.
- Aucune injection définitive n'a été enregistrée dans
  `data/output/acpa_chr19_mutation/` ; la validation de l'injecteur a été faite
  uniquement dans un dossier temporaire.

## Validations effectuées

- `pip check` : aucune dépendance cassée.
- Imports Python principaux et `run_pipeline` réussis.
- Packages R `adegenet`, `ggplot2`, `ape`, `ade4` et `poppr` disponibles.
- Conversion ACPA réelle validée par PLINK : 9 individus, 2 725 marqueurs et
  taux de génotypage de 98,04 %.
- Jeu avec témoins simulés validé par PLINK et KING : 59 individus et 2 725
  marqueurs.
- Injection technique temporaire validée par PLINK : 59 individus, 2 726
  marqueurs et taux de génotypage de 98,25 %.
- Analyse KING genome-wide des témoins réels : 71 individus analysés, 16 233
  marqueurs après QC/pruning, puis 66 individus conservés sans paire détectée
  jusqu'au quatrième degré.
- Contrôle mendélien temporaire : 10 erreurs, dont une au marqueur injecté avec
  les règles de groupe d'exemple.
- Tests ciblés de préparation des données : 16 réussis.
- Tests ciblés de l'étape `04` : 6 réussis.
- Tests ciblés de l'étape `05` : 6 réussis.
- Tests ciblés de l'étape `06` : 7 réussis.
- Tests ciblés de l'étape `07` : 8 réussis.
- Tests ciblés de l'étape `08` : 8 réussis.
- Tests ciblés de l'étape `09` : 7 réussis.
- Les validations séparées des branches `15` et `16` totalisaient respectivement
  171 et 173 tests réussis avant leur intégration locale pour l'étape `17`.
- Contrôle ciblé final de la fenêtre, du cache, du catalogue, de la
  configuration, de l'orchestrateur et de SHAPEIT5 : 45 réussis.
- Dernière exécution combinée antérieure des suites modernes et historiques :
  129 réussis et 4 échecs historiques sans rapport avec les étapes `05` à `11` ;
  la suite historique n'a pas été relancée pour `12.0`.
- Smoke test temporaire avec PLINK 1.9 réel : 2 individus, 23 marqueurs
  autosomiques et 2 marqueurs du chromosome cible, conversion réussie.
- Smoke test temporaire de l'étape `04` avec PLINK 1.9 réel : 2 individus,
  passage de 1 à 2 variants et Mendel `NOT_APPLICABLE` faute de trio.
- Smoke test temporaire de l'étape `05` avec PLINK 1.9 réel : 3 individus et 22
  variants conservés, rapports individu/variant/lot valides, aucun HWE appliqué
  et scan de duplicats `NOT_EVALUATED` faute de 100 variants polymorphes.
- Smoke test temporaire de l'étape `06` avec PLINK 1.9 réel : 3 individus, 22
  marqueurs conservés, 22 autosomes couverts et LD résiduel `NOT_EVALUATED`
  faute de paire dans la fenêtre sur ce jeu minimal.
- Smoke test temporaire de l'étape `07` avec KING 2.3.2 réel : 3 individus,
  une relation parent-enfant déclarée auditée, une exclusion proposée sur la
  seule base du pedigree et validation manuelle correctement requise.
- Smoke test temporaire de l'étape `08` avec PLINK 1.9 et KING 2.3.2 réels :
  3 individus, 22 marqueurs, 11 variants informatifs, 2 composantes calculées,
  aucun outlier et aucun fichier de dosage temporaire publié.
- Smoke test temporaire de la chaîne `01` à `09` avec PLINK 1.9 et KING 2.3.2
  réels : 3 individus, 7 cohortes, 2 contrôles indépendants, 1 porteur total et
  indépendant, 3 individus dans la référence de structure et sur le chromosome
  cible, sans décision manuelle résiduelle.
- Smoke test temporaire de la chaîne `01` à `10` avec PLINK 1.9 et KING 2.3.2
  réels : jeux finaux de 2 individus/11 variants chez les témoins, 1/22 chez le
  porteur indépendant et 3/2 sur le chromosome cible ; aucun individu exclu et
  variant cible conservé sans exception sur cette fixture.
- Smoke test temporaire de la chaîne `01` à `11` avec PLINK 1.9 et KING 2.3.2
  réels : région de 3 individus et 2 variants, une position de carte exacte,
  une position interpolée, variant cible en deuxième position et manifest de
  phasage valide.
- Smoke technique `12.0` : les exécutables installés
  `SHAPEIT5_phase_common` et `SHAPEIT5_phase_rare` répondent tous deux avec la
  version exacte `5.1.1` attendue par le contrat.
- Smoke métadonnées `12.1` contre le serveur EBI : les 44 MD5 des VCF/index des
  22 autosomes et les deux SHA-256 du manifeste et du README correspondent au
  catalogue local.
- Smoke réseau `12.2` : le téléchargeur Python HTTPS a acquis le manifeste EBI
  de 5 092 octets et validé son SHA-256 épinglé, sans télécharger aucun VCF.
- Smoke réel `12.3` avec bcftools 1.21 : VCF bgzip temporaire de 2 échantillons
  et 2 variants phasés, extraction régionale, index tabix et manifeste valides.
- Smoke réel `12.4` avec bcftools 1.21 : 2 variants de référence phasés,
  3 variants d'étude, 2 correspondances canoniques, décomposition, filtrage,
  index tabix et publication temporaire valides.
- Smoke réel `12.5` avec PLINK 1.9 et bcftools 1.21 : 3 individus, dont un trio,
  2 variants, renommage vers les identifiants maître, REF/ALT canoniques et
  index tabix valides ; la complétion d'un REF PLINK `0` par `--a2-allele` est
  également validée sur une sortie temporaire.
- Smoke réel `12.6` avec SHAPEIT5 5.1.1 et bcftools 1.21 : 60 individus,
  200 variants communs et une cible rare ; les gardes ont d'abord bloqué une
  fixture avec erreurs mendéliennes puis un génotype explicite aux mauvais
  allèles. Après correction de la fixture, les deux passes, 201 variants finaux,
  l'indexation et l'attribution `H1` ont réussi ; le PP `0,542` a été classé
  `UNRELIABLE` sous le seuil `0,9`, sans masquer le résultat.
- Smoke réel `12.7` sur les sorties SHAPEIT5 précédentes : douze contrôles QC,
  200 variants communs, 201 variants finaux, 60 individus et un porteur `H1`
  correctement comptabilisés ; le warning PP faible est conservé et le résumé
  final exige une revue manuelle.
- Étape `13` implémentée avec la méthode conservatrice
  `target_centered_exact_ibs_v1` : sélection des unités porteuses indépendantes,
  limites exactes en bp/cM, consensus et matrice pairwise, fréquence de fond
  hors mutation, exclusions explicites et distinction IBS/IBD. Les tests
  synthétiques couvrent un candidat partagé, l'effectif insuffisant,
  l'orchestration et la reprise ; aucune analyse des données réelles n'a été
  lancée.
- Étape `14` implémentée avec `gamma_gandolfo_2014_v1` : correction des bras
  individuels de l'étape 13, reproduction de l'exemple Gamma officiel, modèle
  corrélé primaire, modèle indépendant et leave-one-family-out en sensibilités,
  conversions 25/28/30 ans par génération et règles de non-conclusion sous
  trois unités. Aucun exécutable Gamma externe ni donnée réelle n'est utilisé.
- Étape `15` implémentée avec `plink19_local_ld_secondary_v1` comme analyse
  descriptive secondaire : cohorte principale de témoins indépendants,
  porteurs indépendants et non-porteurs familiaux séparés, univers local commun,
  `r²` génotypique distinct du `D′` haplotypique ML, effectifs et paires
  insuffisants explicitement non évalués, et aucune dépendance vers les sorties
  `13–14`. Les tests et le smoke PLINK utilisent uniquement des données
  synthétiques temporaires ; aucune analyse réelle n'a été lancée.
- Étape `16` implémentée avec `plink19_array_roh_secondary_v1` : fardeau
  genome-wide sur l'intersection exacte des variants chez témoins et porteurs
  indépendants, analyse séparée du chromosome cible, intersection individuelle
  avec le variant, conservation des individus sans segment, avertissement des
  cibles hétérozygotes dans un ROH et `F_ROH` interdit sans dénominateur sourcé.
  Les seuils insuffisants publient `NOT_EVALUATED` sans relâchement automatique.
  Aucun calcul n'a été lancé sur les données réelles.

## Suivi des étapes du pipeline V2

Statuts autorisés : `NON_FAIT`, `EN_COURS`, `FAIT_A_VALIDER`, `VALIDE` et
`BLOQUE`. Une étape est `VALIDE` uniquement lorsque son script, ses contrats,
ses tests, son audit et sa documentation sont cohérents.

| Étape | Nom | Statut | Tests ciblés | Audit d'étape | Remarque |
|---:|---|---|---|---|---|
| `00` | Initialisation du run | `VALIDE` | Oui | Oui | Bootstrap et reprise synthétique validés |
| `01` | Validation des sources | `VALIDE` | Oui | Oui | Fixtures synthétiques et reprise validées |
| `02` | Métadonnées maître | `VALIDE` | Oui | Oui | Revue manuelle et reprise synthétiques validées |
| `03` | Conversion et harmonisation ACPA | `VALIDE` | Oui | Oui | Double sortie PLINK et smoke réel validés |
| `04` | Variant cible | `VALIDE` | Oui | Oui | Injection explicite et smoke PLINK réel validés |
| `05` | QC préliminaire genome-wide | `VALIDE` | Oui | Oui | Missingness technique, alertes et smoke PLINK réel validés |
| `06` | Panel indépendant | `VALIDE` | Oui | Oui | MAF post-QC, pruning LD et smoke PLINK réel validés |
| `07` | Apparentement et duplicats | `VALIDE` | Oui | Oui | KING, pedigree, proposition indépendante et smoke réel validés |
| `08` | Structure populationnelle | `VALIDE` | Oui | Oui | PCA indépendante, projection, outliers et smoke réel validés |
| `09` | Gel des cohortes | `VALIDE` | Oui | Oui | Génotypes explicites, revues SHA, unités et fichiers keep validés |
| `10` | QC final | `VALIDE` | Oui | Oui | Trois politiques, exception cible et smoke réel validés |
| `11` | Région cible | `VALIDE` | Oui | Oui | Carte GRCh38 bornée, cM monotones et smoke réel validés |
| `12` | Référence phasée et phasage | `VALIDE` | Oui | Oui | `12.0–12.7` validées, QC final et smoke réel inclus |
| `13` | Haplotype fondateur et IBD local | `VALIDE` | Oui | Oui | IBS exact centré cible, aucune revendication IBD |
| `14` | Datation du variant | `VALIDE` | Oui | Oui | Gamma corrélé primaire, seuils petits effectifs |
| `15` | LD local secondaire | `FAIT_A_VALIDER` | Oui | Oui | PR brouillon `#5`, intégrée localement pour préparer `17` |
| `16` | ROH secondaire | `VALIDE` | Oui | Oui | Deux périmètres, autozygotie distincte de l'IBS/IBD |
| `17` | Analyses de sensibilité | `NON_FAIT` | Non | Non | — |
| `18` | Visualisations consolidées | `NON_FAIT` | Non | Non | — |
| `19` | Rapport et revue finale | `NON_FAIT` | Non | Non | — |

### Découpage opérationnel de l'étape 12

L'étape reste générique pour tout variant autosomique. DOCK6 et le chromosome
19 ne sont que la première configuration d'étude. Le panel principal utilise
par défaut tous les échantillons et toutes les populations de la référence
phasée ; les sous-populations sont réservées aux analyses de sensibilité.

| Sous-étape | Responsabilité | Statut | Validation attendue |
|---:|---|---|---|
| `12.0` | Figer le contrat et choisir l'adaptateur de phasage | `VALIDE` | SHAPEIT5 5.1.1, MIT, GRCh38, pedigree, formats, paramètres et sonde réelle validés |
| `12.1` | Résoudre le panel depuis un catalogue versionné | `VALIDE` | 1000G GRCh38, 3 202 échantillons, 22 autosomes, README, manifeste et MD5 validés |
| `12.2` | Télécharger ou réutiliser le cache immuable | `VALIDE` | Verrou, publication atomique, reprise, MD5/SHA-256, lecture seule et hors ligne validés |
| `12.3` | Extraire la fenêtre de référence avec tous les échantillons | `VALIDE` | VCF régional bgzip/indexé, génotypes phasés et assemblage concordant |
| `12.4` | Normaliser et harmoniser référence/ACPA | `VALIDE` | Étape orchestrée, schémas, reprise, codes V2, smoke bcftools réel et décision sans FASTA validés |
| `12.5` | Construire les entrées étude, pedigree et référence | `VALIDE` | VCF/index, carte, pedigree, cible rare, ordre/identités, schémas, reprise et smoke réel validés |
| `12.6` | Exécuter l'adaptateur et attribuer le chromosome porteur | `VALIDE` | Deux passes réelles, transmissions, PP, Mendel, génotypes explicites, zones non fiables et reprise validés |
| `12.7` | Publier le QC, les haplotypes et le smoke test | `VALIDE` | Douze contrôles, agrégats, warning/revue, reprise, audit complet et smoke réel validés |

Les sous-étapes `12.1–12.5` doivent pouvoir être testées sans lancer le logiciel
de phasage. Un échec d'empreinte, d'assemblage ou d'allèle doit donc bloquer
avant l'opération externe coûteuse. Le cache de référence est partagé entre les
runs, mais chaque run enregistre exactement la release, les paramètres et les
empreintes qu'il consomme.

### Audits qualité par groupe de cinq étapes

| Groupe | Audit requis après | Statut | Contenu minimal |
|---|---:|---|---|
| `00–04` | étape `04` | `VALIDE` | architecture, contrats, sécurité, tests, audit scientifique |
| `05–09` | étape `09` | `VALIDE` | QC, KING, structure, indépendance et cohortes |
| `10–14` | étape `14` | `VALIDE` | QC final, phasage, haplotypes et datation |
| `15–19` | étape `19` | `NON_FAIT` | analyses secondaires, figures, rapport et reproductibilité |

La revue de maintenabilité du socle réalisée avant l'étape `01` est un contrôle
préliminaire ; elle ne remplace pas l'audit complet obligatoire après l'étape
`04`.

L'audit groupé `05–09` confirme que les filtres automatiques restent limités au
QC technique, que KING et la PCA ne produisent que des propositions, que les
revues sont liées aux artefacts par empreinte, que les décisions ne sont
résolues qu'après succès, que les données individuelles sont classées sensibles
et que chaque cohorte/raison/unité reste reproductible. Aucun blocage nouveau
n'a été identifié ; les limites scientifiques de la PCA exploratoire, de la
sélection gloutonne et des groupes témoins configurés restent documentées.

L'audit groupé `10–14` confirme que le variant cible reste conservé par le QC,
que la carte génétique interdit extrapolation et approximation physique, que le
phasage est confronté aux génotypes moléculaires explicites, que l'étape 13 ne
transforme jamais IBS en preuve IBD et que Gamma ne reçoit que des bras
individuels positifs en cM provenant d'une unité par famille. L'exemple officiel
Gamma est reproduit à l'arrondi, les modes corrélé/indépendant et les
sensibilités sont séparés, et les petits effectifs produisent une non-conclusion
ou un résultat exploratoire. La suite V2 complète passe désormais avec 171
tests. Les limites restantes — résolution des marqueurs, phase, carte,
partage IBS et correction fortuite désactivée — imposent une revue manuelle.

Le contrat de l'étape `15` confirme que le LD reste une branche secondaire
indépendante des étapes `13–14`. Les cohortes sont lues depuis le gel `09`, une
unité familiale dupliquée rend la cohorte non évaluable, les mêmes variants et
distances sont utilisés entre cohortes, et les paires impliquant la cible sont
signalées comme conditionnées par la sélection porteur/non-porteur.

## Problèmes connus

1. Les allèles génomiques REF/ALT de la mutation DOCK6 doivent être confirmés
   sur le brin de référence GRCh38 ; la notation HGVS du transcrit ne suffit pas.
2. Les génotypes individuels réels de la mutation restent à renseigner. Les
   règles de groupe d'exemple créent au moins une incompatibilité mendélienne.
3. `run_pipeline.py` utilise encore des chemins fixes dans
   `data/input/complex_simulation/` et ouvre automatiquement le rapport HTML.
4. L'interface Streamlit écrit `user_input.ped/map`, tandis que le pipeline lit
   `genotype_data.ped/map`.
5. La datation Gamma du pipeline V1 sépare les positions en deux moitiés au
   lieu d'utiliser des longueurs gauche/droite pertinentes autour de la
   mutation ; l'étape V2 `14` ne reprend pas cette approximation.
6. Les exécutables `Gamma` et `gamma` installés dans `/usr/local/bin` contiennent
   une page HTML au lieu d'un programme exécutable. Gamma reste optionnel et ne
   doit pas être déclaré disponible avant réinstallation et validation.
7. Quatre tests historiques échouent encore :
   - confirmation interactive non neutralisée dans un test Gamma ;
   - argument `n` absent d'un test d'estimation d'âge ;
   - deux mocks de prétraitement ne créent pas `filtered_data.ped`.
8. Le simulateur familial historique et les anciens convertisseurs sont encore
   conservés tant que leur dernière logique utile n'a pas été vérifiée.
9. Le pipeline V1 mélange calcul, visualisation et interprétation et peut lire
   des sorties historiques après un échec ; il ne doit pas être relancé pour
   produire un résultat scientifique avant la migration V2.
10. La signature d'étape lie le module, la configuration et les artefacts, mais
    pas encore l'empreinte du binaire externe ni celle des modules partagés ; ce
    durcissement transversal reste requis avant un usage de production.

## Priorités de la prochaine session

1. Définir puis implémenter l'étape `17` des analyses de sensibilité sur la
   branche d'intégration locale des étapes `15–16`.
2. Préparer la table maître réelle et son approbation humaine sans déduire les
   génotypes cibles du statut clinique ou du groupe.
3. Préparer les revues réelles `kinship_exclusion_approval` et
   `population_structure_exclusion_approval` sur un premier run, puis les lier
   par SHA-256 dans la configuration d'un nouveau run destiné au gel.
4. Intégrer les 66 témoins réels dans le QC genome-wide, KING, la structure et le
   gel des cohortes avant toute analyse locale.
5. Confirmer les données moléculaires et génotypes individuels du variant DOCK6,
   puis valider son intégration et Mendel sur le jeu chromosome 19 définitif.
6. Valider ensuite la datation avec des jeux synthétiques de référence avant
   toute analyse réelle.

## Décisions à conserver

- Ne jamais déduire silencieusement un génotype de mutation du statut clinique.
- Utiliser les témoins synthétiques uniquement pour tester le fonctionnement.
- Ne pas écraser `data/input/` ou `data/output/` sans demande explicite et
  `--force` lorsque l'outil le prévoit.
- Considérer `data/output/` comme reproductible, mais pas comme source primaire.
- Vérifier PLINK, l'audit de mutation et les erreurs mendéliennes avant le
  pipeline complet.
- Ne pas lancer le pipeline complet sans prévenir : il écrit de nombreuses
  sorties et ouvre le navigateur.
- Séparer les jeux genome-wide et chromosome 19 ainsi que calcul scientifique,
  visualisation et interprétation.
- Identifier scientifiquement chaque variant par assemblage, chromosome,
  position, REF et ALT ; conserver les identifiants de sonde et rsID comme
  annotations, sans exclure un marqueur uniquement parce que son rsID manque.
- Générer chaque figure uniquement depuis des sorties et audits du run courant.
- Utiliser par défaut toutes les populations du panel phasé de référence ; les
  sélections de populations ne servent qu'aux analyses de sensibilité.
- Ne jamais envoyer les données privées de l'étude vers un service externe ; le
  réseau ne sert qu'à acquérir les références publiques dans un cache immuable.
- Conserver le variant cible dans les données de l'étude même s'il est absent du
  panel, et aligner les variants d'abord par assemblage, chromosome, position,
  REF et ALT plutôt que par rsID ou identifiant de sonde.
- Archiver un script V1 seulement après cartographie de ses dépendances et
  validation de son remplacement V2.

## Reprise rapide

Au début de la prochaine session :

```bash
cd /Users/utilisateur/Documents/python_programme/Effet_fondateur2
source .venv/bin/activate
git checkout main
git pull --ff-only origin main
git status -sb
```

Lire ensuite `AGENTS.md`, ce fichier et `PIPELINE_V2_PRECODE.md`. La prochaine
action de développement attendue est d'implémenter l'étape `17` sur la branche
d'intégration locale de `15–16`, avant toute analyse des données réelles.
