# Suivi de session

Dernière mise à jour : 5 août 2026

## État du dépôt

- Dépôt GitHub : `patrick-mun/Effet_fondateur2`.
- Branche active : `feature/v2-orchestrator-foundation`, synchronisée avec sa
  branche distante.
- PR V2 en brouillon : `#2`, ciblant `main`.
- PR `#1` fusionnée dans `main` le 4 août 2026.
- Commit de fusion : `f2fec8f`.
- Environnement : Python 3.12.13 dans `.venv`.
- Outils disponibles et fonctionnels : PLINK 1.9, KING, bcftools, Rscript et
  packages R du projet. Les commandes `Gamma`/`gamma` présentes dans le `PATH`
  pointent actuellement vers un document HTML invalide et ne sont pas utilisables.

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
- Suite moderne `tests/`, incluant les étapes `01` à `03` et les contrôles de
  maintenabilité V2 : 61 réussis.
- Ensemble des suites Python actuelles : 73 réussis et 4 échecs historiques.
- Smoke test temporaire avec PLINK 1.9 réel : 2 individus, 23 marqueurs
  autosomiques et 2 marqueurs du chromosome cible, conversion réussie.

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
| `04` | Variant cible | `NON_FAIT` | Non | Non | — |
| `05` | QC préliminaire genome-wide | `NON_FAIT` | Non | Non | — |
| `06` | Panel indépendant | `NON_FAIT` | Non | Non | — |
| `07` | Apparentement et duplicats | `NON_FAIT` | Non | Non | — |
| `08` | Structure populationnelle | `NON_FAIT` | Non | Non | — |
| `09` | Gel des cohortes | `NON_FAIT` | Non | Non | Schéma préparé, script non implémenté |
| `10` | QC final | `NON_FAIT` | Non | Non | — |
| `11` | Région cible | `NON_FAIT` | Non | Non | — |
| `12` | Phasage | `NON_FAIT` | Non | Non | Méthode à décider |
| `13` | Haplotype fondateur et IBD local | `NON_FAIT` | Non | Non | Méthode à décider |
| `14` | Datation du variant | `NON_FAIT` | Non | Non | Méthodes à valider |
| `15` | LD local secondaire | `NON_FAIT` | Non | Non | Exploratoire par défaut |
| `16` | ROH secondaire | `NON_FAIT` | Non | Non | Exploratoire par défaut |
| `17` | Analyses de sensibilité | `NON_FAIT` | Non | Non | — |
| `18` | Visualisations consolidées | `NON_FAIT` | Non | Non | — |
| `19` | Rapport et revue finale | `NON_FAIT` | Non | Non | — |

### Audits qualité par groupe de cinq étapes

| Groupe | Audit requis après | Statut | Contenu minimal |
|---|---:|---|---|
| `00–04` | étape `04` | `NON_FAIT` | architecture, contrats, sécurité, tests, audit scientifique |
| `05–09` | étape `09` | `NON_FAIT` | QC, KING, structure, indépendance et cohortes |
| `10–14` | étape `14` | `NON_FAIT` | QC final, phasage, haplotypes et datation |
| `15–19` | étape `19` | `NON_FAIT` | analyses secondaires, figures, rapport et reproductibilité |

La revue de maintenabilité du socle réalisée avant l'étape `01` est un contrôle
préliminaire ; elle ne remplace pas l'audit complet obligatoire après l'étape
`04`.

## Problèmes connus

1. Les allèles génomiques REF/ALT de la mutation DOCK6 doivent être confirmés
   sur le brin de référence GRCh38 ; la notation HGVS du transcrit ne suffit pas.
2. Les génotypes individuels réels de la mutation restent à renseigner. Les
   règles de groupe d'exemple créent au moins une incompatibilité mendélienne.
3. `run_pipeline.py` utilise encore des chemins fixes dans
   `data/input/complex_simulation/` et ouvre automatiquement le rapport HTML.
4. L'interface Streamlit écrit `user_input.ped/map`, tandis que le pipeline lit
   `genotype_data.ped/map`.
5. La datation Gamma sépare actuellement les positions en deux moitiés au lieu
   d'utiliser des longueurs gauche/droite pertinentes autour de la mutation.
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

## Priorités de la prochaine session

1. Implémenter `04_prepare_target_variant_dataset` sur des copies synthétiques,
   avec validation moléculaire et contrôle mendélien simulé.
2. Préparer la table maître réelle et son approbation humaine sans déduire les
   génotypes cibles du statut clinique ou du groupe.
3. Migrer ensuite le QC préliminaire genome-wide et le panel indépendant de
   marqueurs.
4. Intégrer les 66 témoins réels dans le QC genome-wide, KING, la structure et le
   gel des cohortes avant toute analyse locale.
5. Confirmer les données moléculaires et génotypes individuels du variant DOCK6,
   puis valider son intégration et Mendel sur le jeu chromosome 19 définitif.
6. Implémenter ensuite phasage, haplotype fondateur et datation avec jeux
   synthétiques de référence avant analyse réelle.

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
- Archiver un script V1 seulement après cartographie de ses dépendances et
  validation de son remplacement V2.

## Reprise rapide

Au début de la prochaine session :

```bash
cd /Users/utilisateur/Documents/python_programme/Effet_fondateur2
source .venv/bin/activate
git status -sb
```

Lire ensuite `AGENTS.md`, ce fichier et `PIPELINE_V2_PRECODE.md`. La prochaine
action attendue est d'implémenter `04_prepare_target_variant_dataset` en
consommant le jeu du chromosome cible de l'étape `03`, les métadonnées
moléculaires confirmées et uniquement des génotypes individuels explicitement
documentés.
