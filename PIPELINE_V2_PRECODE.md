# Pipeline V2 — architecture, contrats et pré-code

## 1. Statut du document

Ce document décrit l'architecture cible d'un pipeline générique d'étude d'un
possible effet fondateur autour d'un variant cible. Le projet `DOCK6` sur le
chromosome 19 constitue le premier cas d'utilisation, pas une spécialisation de
l'architecture. Il s'agit d'un pré-code : il fixe les responsabilités, les
contrats de fichiers, l'ordre des analyses, les contrôles, les audits et les
critères d'arrêt avant toute réécriture du pipeline.

Il ne constitue ni une implémentation, ni une validation scientifique des
méthodes. Chaque seuil et chaque outil scientifique devront être confirmés,
documentés et testés avant activation dans le pipeline de production.

Principes non négociables :

- ne jamais déduire le génotype du variant cible du statut clinique, familial ou
  du groupe ;
- conserver l'origine explicite de chaque génotype du variant cible ;
- ne jamais modifier les exports ACPA sources ;
- séparer les analyses genome-wide des analyses du chromosome et de la région
  cibles ;
- ne pas traiter des apparentés comme des observations indépendantes ;
- ne pas dater le variant avant d'avoir défini les haplotypes porteurs et leurs
  limites de recombinaison en unités génétiques ;
- ne jamais produire une conclusion biologique codée en dur ;
- rendre chaque résultat traçable jusqu'aux entrées, paramètres, outils et
  exclusions qui l'ont produit.

## 2. Objectifs de l'architecture V2

Le pipeline V2 doit :

1. valider les données et métadonnées avant tout calcul ;
2. construire un jeu autosomique genome-wide pour le QC, l'apparentement et la
   structure populationnelle ;
3. construire un jeu spécifique au chromosome cible pour l'étude locale ;
4. figer des cohortes analytiques explicites et reproductibles ;
5. identifier les chromosomes portant réellement le variant cible ;
6. rechercher et mesurer un haplotype fondateur partagé autour du variant cible ;
7. réserver le LD et les ROH à des analyses secondaires correctement définies ;
8. dater le variant uniquement à partir de longueurs haplotypiques valides ;
9. produire des visualisations à partir de résultats audités ;
10. générer un rapport factuel qui distingue observation, limitation,
    interprétation et conclusion validée manuellement.

### 2.1 Portée initiale et extensions

La première implémentation générique cible les variants autosomiques diploïdes
bialléliques, notamment les SNV et petits indels dont les allèles et coordonnées
sont définis sans ambiguïté sur un assemblage de référence.

Les chromosomes X et Y, l'ADN mitochondrial, les CNV, les variants structuraux,
les mosaïques et les analyses multi-variants nécessitent des contrats et méthodes
spécifiques. Ils doivent être ajoutés par adaptateurs explicites, jamais traités
comme un simple changement de numéro de chromosome.

## 3. Principe général d'orchestration

L'orchestrateur ne contient aucune formule scientifique et ne transforme aucun
génotype. Il pilote des scripts indépendants, vérifie leurs contrats et maintient
l'état du run.

Chaque étape est un processus séparé :

```text
orchestrateur
  -> valide les entrées déclarées de l'étape
  -> calcule la signature de l'étape
  -> lance le script d'analyse dans un sous-processus
  -> vérifie le code de retour
  -> valide les fichiers de sortie et leurs schémas
  -> enregistre les empreintes et métriques
  -> lance le script de visualisation associé si l'étape le prévoit
  -> valide les figures et leurs fichiers de provenance
  -> met à jour l'audit global de façon atomique
  -> autorise ou bloque les étapes dépendantes
```

L'orchestrateur doit appeler les scripts avec une liste d'arguments, sans
`shell=True`, et enregistrer la commande structurée sans exposer de donnée
sensible.

### 3.1 États possibles d'une étape

Chaque étape possède exactement un état parmi :

- `PENDING` : non commencée ;
- `RUNNING` : processus actif ;
- `SUCCEEDED` : calcul et validation des sorties réussis ;
- `FAILED` : erreur technique ou sortie invalide ;
- `BLOCKED` : dépendance ou décision scientifique manquante ;
- `SKIPPED` : étape explicitement désactivée avec justification ;
- `CACHED` : résultats réutilisés après validation stricte de la signature.

Une étape critique en `FAILED` ou `BLOCKED` interdit la génération d'un rapport
présenté comme complet.

### 3.2 Codes de retour standardisés

Tous les scripts d'étape doivent utiliser les mêmes codes :

| Code | Signification |
|---:|---|
| `0` | succès et sorties valides |
| `2` | erreur de configuration ou de schéma d'entrée |
| `3` | échec d'un outil externe |
| `4` | critère scientifique bloquant non satisfait |
| `5` | erreur interne inattendue |

Les avertissements non bloquants sont inscrits dans l'audit et ne modifient pas
le code de retour.

### 3.3 Reprise et cache

La signature d'une étape est calculée à partir de :

- la version du script ;
- la version de son schéma ;
- les paramètres utiles de configuration ;
- les versions des outils externes ;
- les empreintes SHA-256 de toutes les entrées ;
- l'identifiant de la méthode scientifique sélectionnée.

Une étape n'est réutilisable que si sa signature, ses sorties, leurs empreintes
et son audit sont tous valides. La simple présence d'un fichier ne suffit jamais.

## 4. Organisation des runs et immutabilité

Chaque exécution possède un identifiant unique, par exemple :

```text
2026-08-05T083015Z_study_a1b2c3d4
```

Le dossier du run est créé une seule fois. Les étapes écrivent d'abord dans un
dossier temporaire du run, puis publient leurs sorties par renommage atomique.
Une étape réussie ne réécrit jamais ses sorties.

```text
data/runs/<run_id>/
  manifest.json
  events.jsonl
  config.resolved.yaml
  environment.json
  stages/
  report/
```

Les sources restent dans `data/input/`. Les résultats reproductibles restent
dans `data/runs/`. Aucun script scientifique ne lit directement un ancien
dossier de résultats.

## 5. Configuration du projet

Le fichier de configuration principal est validé avant la création du run.
Il référence les autres documents sans recopier silencieusement leur contenu.

### 5.1 Sections minimales

```yaml
project:
  project_id: founder_effect_study
  assembly: GRCh38

target:
  gene: GENE_EXAMPLE
  chromosome: 1
  position_bp: 100000
  ref: A
  alt: G
  transcript: TRANSCRIPT_EXAMPLE
  project_variant_id: target_GRCh38_1_100000_A_G
  rsid: null

inputs:
  acpa_samples_dir: chemin
  sample_metadata: chemin
  target_variant_metadata: chemin
  target_variant_genotypes: chemin
  genetic_map: chemin

tools:
  plink: plink
  king: king
  bcftools: bcftools
  rscript: Rscript
  phasing_adapter: non_configure
  local_ibd_adapter: non_configure

policies:
  overwrite_inputs: false
  allow_clinical_genotype_inference: false
  pseudonymize_figures: true
  stop_on_critical_warning: true

stages:
  stage_name:
    enabled: true
    parameters: {}
```

### 5.2 Profil d'étude DOCK6

DOCK6 est fourni dans un profil d'étude distinct. Les champs moléculaires non
confirmés restent explicitement manquants et bloquent les étapes 04 et suivantes
au lieu d'être remplacés par des valeurs supposées.

```yaml
project:
  project_id: dock6_reunion_founder_effect
  assembly: GRCh38

target:
  gene: DOCK6
  chromosome: 19
  position_bp: null
  ref: null
  alt: null
  transcript: null
  project_variant_id: null
  rsid: null
```

### 5.3 Paramètres scientifiques

Les seuils ne doivent pas être dispersés dans les scripts. Chaque seuil doit
posséder :

- un nom explicite ;
- une unité ;
- une valeur ;
- une justification ;
- une référence ou une décision de projet ;
- une indication `exploratoire` ou `validé` ;
- la date et l'auteur de la décision.

Un seuil sans justification bloque l'étape qui en dépend en mode production.

## 6. Contrats de fichiers communs

### 6.1 Règles générales

- texte en UTF-8 ;
- TSV pour les tables humaines et scientifiques ;
- JSON pour les audits et manifests ;
- JSON Lines pour le journal append-only ;
- colonnes et types validés par schéma versionné ;
- valeurs manquantes représentées selon le schéma, jamais par mélange de
  chaînes arbitraires ;
- coordonnées génomiques toujours accompagnées de l'assemblage ;
- unités incluses dans le nom de colonne lorsque nécessaire : `_bp`, `_kb`,
  `_cm`, `_generations` ;
- identifiants d'échantillons sans espace et uniques dans le projet ;
- aucun identifiant réel dans une figure destinée au partage par défaut.

### 6.2 Table maître des échantillons

Fichier cible : `metadata/samples.master.tsv`.

Colonnes minimales :

| Colonne | Rôle |
|---|---|
| `SAMPLE_ID` | identifiant stable interne |
| `SOURCE_FILE` | chemin relatif de l'export source |
| `FID` | famille PLINK |
| `IID` | individu PLINK |
| `PID` | père déclaré ou `0` |
| `MID` | mère déclarée ou `0` |
| `SEX` | sexe déclaré selon convention documentée |
| `CLINICAL_STATUS` | statut clinique explicite |
| `GROUP_LABEL` | groupe descriptif, non utilisé comme génotype |
| `TARGET_GENOTYPE` | génotype explicite ou valeur manquante |
| `TARGET_GENOTYPE_SOURCE` | origine du génotype |
| `ARRAY_BATCH` | lot ou série ACPA |
| `INCLUDE_GENOMEWIDE` | éligibilité technique initiale |
| `INCLUDE_TARGET_CHROMOSOME` | éligibilité technique initiale |
| `NOTES_CODE` | référence vers une note contrôlée |

Contraintes :

- unicité de `SAMPLE_ID` et de `FID + IID` ;
- existence de `SOURCE_FILE` ;
- cohérence parent-enfant sans inférer un lien absent ;
- `TARGET_GENOTYPE` obligatoirement associé à `TARGET_GENOTYPE_SOURCE` ;
- impossibilité de remplir `TARGET_GENOTYPE` depuis `CLINICAL_STATUS` ou
  `GROUP_LABEL`.

Encodage canonique V2 :

- champ manquant : cellule TSV vide, jamais `NA`, `N/A`, `None` ou `-` ;
- booléens : `true` et `false` en minuscules ;
- `SEX` : `UNKNOWN`, `MALE`, `FEMALE` ou `NOT_APPLICABLE` ;
- `CLINICAL_STATUS` : `AFFECTED`, `UNAFFECTED`, `UNKNOWN` ou `NOT_REPORTED` ;
- génotype cible diploïde : `REF/ALT` non phasé ou `REF|ALT` phasé avec les
  séquences alléliques explicites ;
- `SOURCE_FILE` : chemin POSIX relatif au répertoire source autorisé ;
- identifiants : caractères ASCII alphanumériques, point, tiret et soulignement,
  sans espace.

### 6.3 Table des cohortes figées

Fichier cible : `cohorts/cohorts.frozen.tsv`.

Une ligne par échantillon et par cohorte analytique :

| Colonne | Rôle |
|---|---|
| `COHORT_ID` | nom versionné de la cohorte |
| `SAMPLE_ID` | individu concerné |
| `ROLE` | porteur, non-porteur, témoin, apparenté, autre |
| `INCLUDED` | booléen |
| `EXCLUSION_CODE` | motif contrôlé si exclu |
| `DECISION_SOURCE` | étape ou validation manuelle |
| `INDEPENDENT_UNIT_ID` | unité indépendante retenue |
| `FAMILY_REPRESENTATIVE` | indique une sélection par famille |

Une cohorte figée est immuable dans le run. Toute nouvelle décision crée un
nouveau run ou une nouvelle version explicitement liée.

Les rôles canoniques sont `CARRIER`, `NON_CARRIER`, `CONTROL`, `RELATED` et
`OTHER`. Une ligne exclue (`INCLUDED=false`) exige un `EXCLUSION_CODE`; une
ligne incluse laisse ce champ vide. `FAMILY_REPRESENTATIVE=true` exige une unité
`INDEPENDENT_UNIT_ID` explicite.

### 6.4 Formats génétiques canoniques

Pour compatibilité avec PLINK 1.9 et KING, le format binaire canonique est le
triplet BED/BIM/FAM. PED/MAP est réservé aux imports et aux exports nécessitant
les allèles textuels.

Chaque triplet doit être accompagné de :

- `dataset.json` : assemblage, population, filtres, effectifs et provenance ;
- `samples.tsv` : ordre exact des individus ;
- `variants.tsv` : ordre exact des variants, coordonnées et statut QC ;
- `checksums.sha256` : empreintes de tous les artefacts.

### 6.5 Descripteurs d'artefacts entre étapes

Un script ne reçoit pas une collection de chemins libres. L'orchestrateur lui
fournit un `stage_inputs.json` validé qui décrit chaque artefact autorisé.

Champs minimaux d'un artefact :

| Champ | Description |
|---|---|
| `artifact_id` | identifiant unique dans le run |
| `artifact_type` | type contrôlé : `plink_bed`, `sample_table`, etc. |
| `path` | chemin relatif au run |
| `media_type` | format physique du fichier |
| `schema_name` | schéma applicable |
| `schema_version` | version exacte du schéma |
| `sha256` | empreinte du fichier |
| `producer_stage` | étape ayant produit l'artefact |
| `producer_signature` | signature de cette étape |
| `assembly` | assemblage si donnée génomique |
| `sample_set_id` | identifiant de l'ensemble d'individus |
| `variant_set_id` | identifiant de l'ensemble de variants |
| `sensitivity` | `public`, `internal`, `sensitive_genetic` |

Pour un triplet BED/BIM/FAM, le descripteur représente un artefact composé et
contient l'empreinte de chaque fichier. Un script doit refuser deux artefacts
génétiques dont les `sample_set_id`, `variant_set_id` ou assemblages sont
incompatibles avec son contrat.

Chaque étape publie symétriquement `stage_outputs.json`. L'orchestrateur le
valide avant d'ajouter les nouveaux artefacts au registre du run.

### 6.6 Schémas des tables scientifiques critiques

#### Métriques QC des individus

`sample_qc.tsv` contient au minimum :

```text
SAMPLE_ID  SAMPLE_SET_ID  N_VARIANTS_CALLED  CALL_RATE
HET_RATE  F_MISS  ARRAY_BATCH  QC_STATUS  EXCLUSION_CODE
```

Les colonnes supplémentaires doivent être versionnées. `QC_STATUS` ne peut
prendre que `PASS`, `WARN`, `FAIL` ou `NOT_EVALUATED`.

#### Métriques QC des variants

`variant_qc.tsv` contient au minimum :

```text
VARIANT_ID  PROBE_ID  RSID  CHR  POSITION_BP  ASSEMBLY  REF  ALT
CALL_RATE  MAF  HWE_CONTROL_P  QC_STATUS  EXCLUSION_CODE
```

`HWE_CONTROL_P` reste manquant avant le gel des témoins indépendants. Le variant
cible est identifié par un champ dédié dans la version complète du schéma et ne
peut pas être supprimé silencieusement par un filtre général.

#### Paires d'apparentement

`kinship_pairs.tsv` contient au minimum :

```text
SAMPLE_ID_1  SAMPLE_ID_2  N_SNP  KINSHIP  IBS0
INFERRED_DEGREE  DECLARED_RELATION  CONCORDANCE_STATUS
```

Le fichier est classé `sensitive_genetic`. Les résumés partageables ne
contiennent que des comptes agrégés.

#### Coordonnées de structure populationnelle

`population_scores.tsv` contient au minimum :

```text
SAMPLE_ID  SAMPLE_SET_ID  PC1  PC2  PC3
PROJECTED  OUTLIER_STATUS  OUTLIER_REASON
```

Le nombre total de composantes et leur variance expliquée sont enregistrés dans
un artefact séparé `population_eigenvalues.tsv`.

#### Haplotypes porteurs et segments partagés

`founder_segments.tsv` contient une ligne par unité porteuse indépendante :

```text
INDEPENDENT_UNIT_ID  SAMPLE_ID  FAMILY_ID  CARRIER_HAPLOTYPE_ID
TARGET_VARIANT_ID  LEFT_BOUND_BP  TARGET_BP  RIGHT_BOUND_BP
LEFT_LENGTH_CM  RIGHT_LENGTH_CM  PHASING_CONFIDENCE
SEGMENT_METHOD  SEGMENT_STATUS  EXCLUSION_CODE
```

Les longueurs représentent des distances entre le variant cible et chaque limite,
jamais des positions absolues. Les trois positions en bp doivent vérifier :

```text
LEFT_BOUND_BP <= TARGET_BP <= RIGHT_BOUND_BP
```

#### Entrée de datation

`variant_age_input.tsv` est dérivé de `founder_segments.tsv` et contient :

```text
INDEPENDENT_UNIT_ID  FAMILY_ID  LEFT_LENGTH_CM  RIGHT_LENGTH_CM
INCLUDED  EXCLUSION_CODE  MAP_VERSION  SEGMENT_METHOD
```

Une unité exclue reste visible avec `INCLUDED=false`. Aucun script de datation
ne reçoit directement un MAP, un RAW ou une liste de positions de SNP.

### 6.7 Politique d'identification des variants et des rsID

Chaque variant possède un `VARIANT_ID` interne unique et stable dans le run. Le
variant cible possède en plus un `PROJECT_VARIANT_ID` défini par le projet. Leur
identité scientifique repose d'abord sur :

```text
assemblage + chromosome + position + REF + ALT
```

Les identifiants suivants sont des annotations distinctes :

- `PROBE_ID` : identifiant de sonde de la puce ;
- `RSID` : identifiant dbSNP lorsqu'il existe et est résolu sans ambiguïté ;
- `HGVS_C` et `HGVS_P` : descriptions liées à un transcrit déclaré ;
- identifiant interne par coordonnées, utilisé comme repli reproductible.

Un rsID n'est pas obligatoire pour PLINK, KING, le phasage familial ou la
comparaison d'haplotypes si les variants sont uniques, ordonnés, alignés sur le
même assemblage et orientés avec des allèles cohérents. Le pipeline ne doit donc
pas exclure automatiquement un marqueur uniquement parce qu'aucun rsID n'a été
résolu.

Le rsID reste utile pour :

- joindre dbSNP et des annotations externes ;
- comparer les résultats à la littérature ;
- harmoniser avec certains panels de référence ;
- faciliter la revue humaine et l'échange avec d'autres projets.

Toute attribution de rsID conserve sa source, sa version, sa méthode et son
statut : `resolved`, `ambiguous`, `conflicting` ou `unresolved`. Un rsID ne
remplace jamais la vérification des coordonnées et des allèles.

### 6.8 Classification et protection des données

- les exports ACPA, génotypes, paires KING et haplotypes individuels sont
  `sensitive_genetic` ;
- les audits globaux utilisent des métriques agrégées ;
- les identifiants affichés sont pseudonymisés par défaut ;
- les rapports partageables n'incluent pas les tables individuelles sensibles ;
- les journaux et messages d'erreur ne recopient jamais une ligne de génotype ;
- les chemins absolus et noms de fichiers sources peuvent être masqués dans les
  rapports exportés ;
- les politiques de conservation et d'accès doivent être définies avant tout
  partage hors de l'environnement local.

## 7. Architecture d'audit

### 7.1 Manifest global

`manifest.json` représente l'état courant du run. Il est réécrit de façon
atomique par l'orchestrateur uniquement.

Le manifest est un état de contrôle, pas un artefact qui se référence lui-même.
Il conserve les empreintes de `audit.json` et `stage_outputs.json` de chaque
étape ; sa propre intégrité repose sur l'écriture atomique et le journal
append-only, ce qui évite toute dépendance circulaire d'empreinte.

Contenu minimal :

- identifiant du run ;
- configuration résolue et son empreinte ;
- environnement et versions ;
- liste ordonnée des étapes ;
- état, signature, dates et durée de chaque étape ;
- chemins relatifs des audits d'étape ;
- statut global : incomplet, valide techniquement, prêt pour revue scientifique
  ou bloqué ;
- liste des décisions manuelles encore requises.

### 7.2 Journal append-only

`events.jsonl` reçoit des événements structurés :

```json
{
  "timestamp": "...",
  "stage": "07_kinship",
  "event": "stage_succeeded",
  "severity": "INFO",
  "details": {}
}
```

Le journal ne contient ni génotypes, ni secrets, ni commandes exposant des
identifiants sensibles.

### 7.3 Audit d'étape

Chaque étape réussie publie `audit.json` avec :

- `schema_version` ;
- `stage_id` et `method_id` ;
- signature de l'étape ;
- entrées et SHA-256 ;
- sorties et SHA-256 ;
- paramètres et unités ;
- outils et versions ;
- effectifs avant et après ;
- métriques principales ;
- exclusions avec codes de motif ;
- avertissements ;
- contrôles passés et échoués ;
- limites connues ;
- visualisations attendues ;
- validation manuelle éventuellement requise.

Les exclusions individuelles détaillées restent dans un TSV local audité. Le
manifest global n'enregistre que les comptes agrégés.

### 7.4 Provenance des figures

Chaque figure `figure_name.png` ou `figure_name.pdf` possède un fichier
`figure_name.figure.json` contenant :

- script et version ;
- fichiers sources et empreintes ;
- paramètres graphiques utiles ;
- filtres appliqués ;
- effectif réellement représenté ;
- groupes et unités ;
- légende complète ;
- avertissements d'interprétation ;
- statut de pseudonymisation.

## 8. Graphe des étapes

```text
00 initialisation -> 01 validation des sources -> 02 métadonnées maître
  -> 03 conversion et harmonisation ACPA

03 -> 04 validation et injection du variant sur le chromosome cible
03 -> 05 QC préliminaire genome-wide -> 06 panel indépendant
   -> 07 apparentement et duplicats -> 08 structure populationnelle

04 + 07 + 08 -> 09 résolution et gel des cohortes -> 10 QC final
04 + 10 -> 11 préparation de la région cible -> 12 phasage
        -> 13 haplotype fondateur et IBD local -> 14 datation
10 + 11 -> 15 LD local secondaire
10 + 11 -> 16 ROH secondaire
13 + 14 si activée + étapes secondaires activées -> 17 analyses de sensibilité
17 -> 18 visualisations consolidées -> 19 rapport et revue finale
```

Les branches `15` et `16` ne sont pas des prérequis de la datation. La datation
dépend du phasage et des longueurs haplotypiques validées, pas des ROH globaux.

### 8.1 Matrice de dépendances et criticité

| Étape | Dépendances directes | Criticité par défaut |
|---|---|---|
| `00` | configuration | critique |
| `01` | `00` | critique |
| `02` | `01` | critique et validation manuelle |
| `03` | `02` | critique |
| `04` | `02`, sortie du chromosome cible de `03` | critique pour les analyses du variant cible |
| `05` | sortie genome-wide de `03` | critique |
| `06` | `05` | critique |
| `07` | `02`, `05`, `06` | critique et validation des exclusions |
| `08` | `07` | critique pour figer une cohorte comparable |
| `09` | `04`, `07`, `08` | critique et validation manuelle |
| `10` | `03`, `04`, `09` | critique |
| `11` | `04`, `10` | critique |
| `12` | `09`, `11` | critique pour haplotype et datation |
| `13` | `12` | critique pour conclure sur un haplotype fondateur |
| `14` | `13` | critique uniquement si une datation est demandée |
| `15` | `09`, `10`, `11` | secondaire, exploratoire par défaut |
| `16` | `09`, `10`, `11` | secondaire, exploratoire par défaut |
| `17` | `13`, `14` si activée, `15` et `16` si activées | critique pour les résultats concernés |
| `18` | toutes les étapes activées dans un état terminal | critique pour le rapport |
| `19` | `18` et tous les audits requis | critique |

Une étape secondaire désactivée est `SKIPPED`, pas `SUCCEEDED`. Son absence doit
être visible dans le rapport et ne doit jamais être remplacée par une sortie
d'un run précédent.

## 9. Pré-code détaillé des étapes

### Étape 00 — Initialisation du run

**Script cible** : `stages/initialize_run.py`.

**Responsabilité** : créer un run vide et reproductible.

Cette étape est un bootstrap appelé une seule fois par la CLI avant la boucle
générique du DAG. Elle est auditée comme les autres étapes dans
`00_initialize_run/`, mais n'est pas relancée par l'orchestrateur après création
du run.

**Entrées** :

- configuration principale ;
- schémas de configuration ;
- environnement local.

**Traitements** :

1. valider le YAML/JSON ;
2. résoudre tous les chemins sans suivre de chemin hors projet non autorisé ;
3. vérifier que les sources existent sans les ouvrir en écriture ;
4. capturer les versions Python et outils ;
5. créer l'identifiant de run ;
6. écrire la configuration résolue ;
7. initialiser `manifest.json` et `events.jsonl`.

**Sorties de contrôle du run** : environnement, configuration résolue et
manifest initial. Ces documents racine sont administrés par l'orchestrateur.
L'artefact propre à l'étape est un résumé d'initialisation non sensible.

**Critères bloquants** : configuration invalide, source absente, assemblage non
déclaré ou outil requis par une étape activée indisponible.

### Étape 01 — Validation des sources

**Script cible** : `stages/validate_sources.py`.

**Responsabilité** : inventorier les fichiers sources sans produire de données
génétiques dérivées.

**Contrôles ACPA** :

- nombre de fichiers attendu et observé ;
- extension, encodage et en-tête ;
- colonnes obligatoires ;
- assemblage déclaré ;
- type de puce ;
- présence des chromosomes autosomiques et du chromosome cible ;
- unicité des noms ;
- taille et empreinte SHA-256 ;
- cohérence des versions d'annotation ;
- détection des fichiers tronqués ou vides.

**Sorties** :

- `source_inventory.tsv` ;
- `source_qc.tsv` ;
- `audit.json`.

**Visualisation** : aucune donnée génétique ; diagramme facultatif des comptes
de fichiers par lot et par statut de validation.

### Étape 02 — Construction des métadonnées maître

**Script cible** : `stages/build_sample_registry.py`.

**Responsabilité** : créer le registre unique des individus et familles.

**Traitements** :

1. charger la table fournie par l'utilisateur ;
2. valider les identifiants PLINK ;
3. vérifier les doublons de fichier et d'individu ;
4. valider la structure parentale déclarée ;
5. joindre les fichiers ACPA sans inférence clinique ;
6. valider la provenance des génotypes du variant cible ;
7. attribuer un pseudonyme graphique stable si configuré ;
8. publier `samples.master.tsv`.

**Critères bloquants** : génotype du variant cible sans source, fichier assigné à
plusieurs individus, identifiant dupliqué ou incohérence de pedigree non résolue.

**Validation manuelle** : la table maître doit être approuvée avant conversion.

### Étape 03 — Conversion et harmonisation ACPA

**Script cible** : `stages/convert_acpa.py`.

Le script publie deux artefacts distincts au cours de la même étape afin de lire
chaque export ACPA une seule fois autant que possible. Les adaptateurs internes
peuvent être séparés, mais l'orchestrateur ne considère la conversion comme
réussie que si les deux jeux attendus et leur table de correspondance sont
valides.

**Responsabilité** : lire chaque export ACPA une fois autant que possible et
produire deux jeux cohérents.

**Jeu A — genome-wide autosomique** : chromosomes 1 à 22, destiné au QC,
apparentement et structure.

**Jeu B — chromosome cible** : densité maximale compatible, destiné à l'analyse
de la région cible déclarée dans la configuration.

**Contrôles** :

- coordonnées GRCh38 ;
- orientation des allèles ;
- variants bialléliques ;
- conflits de position et d'identifiant ;
- génotypes invalides convertis explicitement en manquants ;
- intersection ou union de sondes selon une politique déclarée ;
- taux de génotypage par fichier et par chromosome ;
- ordre identique et documenté des individus.

**Sorties** : BED/BIM/FAM, tables d'audit de variants, rapport de conversion et
empreintes pour chacun des deux jeux.

**Visualisations** :

- taux de génotypage par échantillon ;
- nombre de marqueurs par chromosome ;
- distribution des marqueurs manquants par lot ;
- densité des sondes sur le chromosome cible.

### Étape 04 — Variant cible : validation et injection contrôlée

**Script cible** : `stages/prepare_target_variant_dataset.py`.

**Responsabilité** : vérifier la définition moléculaire puis injecter, si
nécessaire, le variant cible dans une copie du jeu du chromosome cible.

**Entrées obligatoires** : assemblage, chromosome, position, REF, ALT,
transcrit, HGVS c., HGVS p. si connue, méthode de mesure et génotype individuel
explicite.

**Contrôles** :

- REF/ALT compatibles avec GRCh38 ;
- position unique ;
- cohérence HGVS/transcrit ;
- génotype cible fourni pour chaque individu inclus ;
- aucune règle dérivée du statut clinique ;
- contrôle mendélien après injection ;
- comparaison des empreintes avant/après.

**Sorties** : jeu du chromosome cible enrichi, audit individuel local, rapport
Mendel, rapport d'injection et liste des décisions bloquantes.

**Critère bloquant** : toute incompatibilité moléculaire ou génotype obligatoire
non renseigné.

### Étape 05 — QC préliminaire genome-wide

**Script cible** : `stages/qc_preliminary.py`.

**Responsabilité** : produire un jeu techniquement fiable pour l'apparentement,
sans HWE dépendant d'une cohorte encore non résolue.

**Analyses** :

- taux de données manquantes par individu et variant ;
- hétérozygotie autosomique ;
- discordance de sexe si les données et métadonnées le permettent ;
- variants non bialléliques ;
- duplicats techniques potentiels ;
- fréquence allélique minimale adaptée au panel d'apparentement ;
- comparaison des métriques par lot ACPA ;
- différentiel de données manquantes entre lots, sans interprétation clinique.

**Important** : aucun filtre HWE définitif à cette étape, car les témoins
indépendants ne sont pas encore figés.

**Sorties** : jeu pré-QC, métriques individuelles, métriques variants et liste
d'alertes. Les exclusions automatiques doivent être limitées aux erreurs
techniques clairement définies.

**Visualisations** : missingness, hétérozygotie, lot, MAF et matrice de données
manquantes agrégée.

### Étape 06 — Construction du panel indépendant de marqueurs

**Script cible** : `stages/build_kinship_panel.py`.

**Responsabilité** : construire un panel autosomique genome-wide adapté à KING
et à la structure populationnelle.

**Traitements** :

1. exclure chromosomes sexuels et mitochondrie ;
2. appliquer les critères techniques préliminaires ;
3. recalculer la MAF après les exclusions individuelles techniques ;
4. exclure les régions complexes configurées si la méthode l'exige ;
5. appliquer un pruning LD documenté ;
6. vérifier la couverture des 22 autosomes ;
7. calculer le nombre opérationnel de marqueurs et leur répartition ;
8. publier un BED/BIM/FAM dédié, distinct du jeu du chromosome cible.

**Critère bloquant** : couverture insuffisante ou concentration excessive des
marqueurs dans une faible partie du génome.

**Visualisations** : marqueurs conservés par chromosome, distances entre
marqueurs et distribution agrégée du LD résiduel. L'absence de paire évaluable
est notée `NOT_EVALUATED`.

### Étape 07 — Apparentement, duplicats et concordance du pedigree

**Script cible** : `stages/infer_kinship.py`.

**Responsabilité** : estimer l'apparentement genome-wide et le comparer aux
relations déclarées.

**Traitements** :

- KING sur le panel indépendant ;
- classification des degrés avec seuils versionnés ;
- recherche de duplicats ou jumeaux ;
- comparaison relation observée/relation déclarée ;
- identification séparée des relations attendues et cryptiques ;
- construction d'un graphe de parenté ;
- proposition, sans application silencieuse, d'un ensemble indépendant maximal ;
- stratégie de sélection tenant compte de la qualité des génotypes.

La première implémentation utilise une sélection gloutonne déterministe orientée
qualité. Elle garantit un ensemble maximal, pas nécessairement un ensemble de
cardinalité maximum ; cette limite est inscrite dans l'audit.

**Sorties** :

- paires et coefficients dans un fichier sensible local ;
- résumé agrégé par degré ;
- concordance pedigree/KING ;
- proposition d'exclusions ;
- audit sans identifiants dans sa partie partageable.

**Validation manuelle** : toute exclusion d'individu doit être approuvée ou
provenir d'une politique pré-validée explicitement activée.

**Visualisations** : réseau pseudonymisé, matrice de kinship, distribution des
coefficients et tableau de concordance.

### Étape 08 — Structure populationnelle

**Script cible** : `stages/analyze_population_structure.py`.

**Responsabilité** : détecter les axes de structure et les individus aberrants
sur le panel genome-wide indépendant.

**Stratégie** :

- estimer les axes sur une base aussi indépendante que possible ;
- ne pas laisser plusieurs membres d'une même famille déterminer les axes ;
- projeter ensuite les apparentés si nécessaire ;
- documenter les références externes si elles sont introduites ;
- distinguer PCA exploratoire et décision d'exclusion ;
- réserver la DAPC à une question définie, avec génotypes correctement encodés,
  nombre de composantes validé et contrôle du surapprentissage.

**Sorties** : coordonnées, variance expliquée, métriques d'outliers et décisions
proposées.

**Visualisations** : scree plot, PCA pseudonymisée par cohorte, lot et famille,
ainsi que projection des individus exclus potentiels.

### Étape 09 — Résolution et gel des cohortes

**Script cible** : `stages/freeze_cohorts.py`.

**Responsabilité** : transformer les résultats QC, KING, structure, pedigree et
génotype cible en cohortes analytiques explicites.

**Cohortes minimales** :

- `controls_unrelated` ;
- `target_carriers_all` ;
- `target_carriers_independent` ;
- `family_noncarriers` ;
- `affected_noncarriers` si présents ;
- `target_chromosome_all_qc` ;
- `genomewide_structure_reference`.

**Règles** :

- un individu peut appartenir à plusieurs cohortes avec des rôles distincts ;
- l'indépendance est définie par une unité explicite ;
- pour la datation, la stratégie doit éviter de compter plusieurs fois un même
  haplotype ancestral familial ;
- aucune cohorte ne se fonde implicitement sur le seul statut `ATTEINT` ;
- chaque exclusion possède un code et une source.

**Sorties** : `cohorts.frozen.tsv`, fichiers PLINK `--keep` et audit de décision.

**Critère bloquant** : absence de génotypes cibles fiables ou impossibilité
de définir une unité porteuse indépendante.

### Étape 10 — QC final par cohorte

**Script cible** : `stages/qc_final.py`.

**Responsabilité** : appliquer les filtres propres à chaque analyse après gel
des cohortes.

**Jeu témoins indépendants** : HWE, MAF, missingness et lot.

**Jeu porteurs** : contrôle de qualité sans utiliser HWE comme critère
d'exclusion du variant cible ou d'un signal potentiellement causal.

**Jeu local cible** : conservation forcée et auditée du variant cible, même s'il
échoue un filtre exploratoire standard ; toute exception est enregistrée.

**Sorties** : jeux finaux versionnés par cohorte et matrices de concordance des
variants.

**Visualisations** : HWE chez témoins, MAF par cohorte, missingness comparative
et densité locale finale.

### Étape 11 — Préparation de la région cible

**Script cible** : `stages/prepare_target_region.py`.

**Responsabilité** : préparer les données autour du variant cible pour le phasage
et les analyses locales.

**Traitements** :

- validation finale des coordonnées ;
- sélection d'une fenêtre physique large configurable ;
- jointure avec une carte de recombinaison GRCh38 ;
- interpolation contrôlée des positions en cM ;
- vérification de monotonie de la carte ;
- contrôle des allèles et doublons ;
- conservation de l'ordre des variants ;
- création des formats requis par l'adaptateur de phasage.

**Interdiction** : remplacer une carte génétique absente par `1 Mb = 1 cM` en
mode production.

**Visualisations** : densité physique et génétique, taux de recombinaison local
et position du variant cible.

### Étape 12 — Phasage et identification du chromosome porteur

**Script cible** : `stages/phase_target_region.py` avec adaptateur d'outil externe.

**Responsabilité** : produire des haplotypes et identifier celui qui porte le
variant cible chez chaque porteur.

**Prérequis à décider avant implémentation** : outil de phasage, version,
référence éventuelle, compatibilité GRCh38 et prise en charge du pedigree.

**Traitements** :

- phasage pedigree-aware lorsque disponible ;
- utilisation des génotypes cibles explicites ;
- contrôle des inversions d'allèles ;
- contrôle Mendel avant et après ;
- score ou niveau de confiance du phasage ;
- attribution du chromosome porteur ;
- marquage des régions non fiables.

**Critère bloquant** : variant cible non phasable chez un nombre suffisant d'unités
indépendantes selon la méthode de datation choisie.

**Visualisations** : haplotypes familiaux, transmissions compatibles et zones de
phase incertaine, avec pseudonymes.

### Étape 13 — Haplotype fondateur et IBD local

**Script cible** : `stages/infer_founder_haplotype.py`.

**Responsabilité** : tester l'existence d'un segment ancestral partagé autour de
du variant cible et mesurer ses limites.

**Traitements** :

1. sélectionner uniquement les chromosomes porteurs validés ;
2. comparer les allèles phasés autour du variant cible ;
3. détecter les segments partagés avec une méthode locale choisie et versionnée ;
4. distinguer IBS et IBD selon les informations disponibles ;
5. mesurer pour chaque unité indépendante la limite gauche et droite ;
6. convertir les limites en bp et cM ;
7. comparer la fréquence du segment chez témoins et non-porteurs ;
8. calculer un intervalle partagé central et ses variantes de sensibilité ;
9. signaler les marqueurs informatifs insuffisants ;
10. ne pas conclure à un fondateur unique si plusieurs haplotypes porteurs
    incompatibles sont observés.

**Sorties** : table par chromosome porteur, limites gauche/droite, matrice de
partage, haplotype consensus éventuel et audit des discordances.

**Visualisations** :

- piste d'haplotypes alignés autour du variant cible ;
- segments partagés par unité indépendante ;
- matrice de concordance allélique ;
- intervalle consensus en bp et cM ;
- fréquence du segment selon les cohortes.

### Étape 14 — Datation du variant cible

**Script cible** : `stages/estimate_variant_age.py`.

**Responsabilité** : estimer l'âge uniquement depuis les longueurs génétiques
gauche/droite des haplotypes porteurs indépendants.

**Entrées obligatoires** :

- une ligne par unité porteuse indépendante ;
- longueur gauche en cM ;
- longueur droite en cM ;
- méthode de définition des limites ;
- confiance du phasage ;
- fréquence pertinente à la méthode ;
- paramètres démographiques requis ;
- version de la carte génétique.

**Traitements** :

- validation des longueurs positives et des unités ;
- exécution de Gamma avec paramètres explicites si la méthode est validée ;
- méthode alternative ou analyse de comparaison lorsque disponible ;
- intervalles de confiance ;
- bootstrap ou sensibilité leave-one-family-out ;
- comparaison avec et sans unités incertaines ;
- signalement des hypothèses non vérifiables.

**Interdictions** :

- utiliser des positions absolues de marqueurs comme longueurs ;
- couper une liste de marqueurs en deux moitiés ;
- utiliser le nombre total d'individus comme nombre de chromosomes porteurs ;
- annoncer un nombre de générations codé en dur.

**Visualisations** : longueurs gauche/droite, distributions bootstrap,
estimations par scénario et forêt des intervalles de confiance.

### Étape 15 — LD local secondaire

**Script cible** : `stages/analyze_local_ld.py`.

**Responsabilité** : décrire le LD autour du variant cible sans le confondre avec
l'apparentement ou l'haplotype IBD.

**Cohortes** : témoins indépendants, porteurs indépendants et éventuellement
non-porteurs familiaux dans des analyses séparées.

**Traitements** :

- calcul direct de `r²` et `D′` par l'outil ;
- mêmes fenêtres, variants et conventions entre cohortes comparées ;
- prise en compte de l'effectif très faible des porteurs ;
- exclusion ou pondération familiale définie avant calcul ;
- aucune reconstruction manuelle de `D` depuis `D′` et les MAF sans modèle
  allélique complet ;
- analyses descriptives, sans test naïf sur des millions de paires corrélées.

**Visualisations** : matrices de LD comparables, décroissance selon distance en
cM et zoom centré sur le variant cible.

### Étape 16 — ROH secondaire

**Script cible** : `stages/analyze_roh.py`.

**Responsabilité** : étudier l'autozygotie comme phénomène distinct de
l'haplotype fondateur.

**Traitements** :

- paramètres adaptés à la densité ACPA ;
- analyse genome-wide du fardeau ROH si les données le permettent ;
- analyse locale du chromosome cible séparée ;
- calcul d'intersections par chromosome et par segment ;
- comparaison descriptive entre cohortes indépendantes ;
- recherche du variant cible dans un ROH par individu, jamais uniquement au niveau
  global ;
- distinction explicite entre homozygotie, IBS et IBD.

**Visualisations** : fardeau ROH, longueurs, pistes individuelles par chromosome
et chevauchements locaux correctement calculés.

### Étape 17 — Analyses de sensibilité

**Script cible** : `stages/run_sensitivity_analyses.py`.

**Responsabilité** : mesurer la stabilité des résultats principaux.

Scénarios minimaux :

- retrait d'une famille à la fois ;
- retrait d'un porteur à faible confiance de phase ;
- plusieurs seuils de définition de l'haplotype ;
- plusieurs fenêtres locales ;
- cartes génétiques ou méthodes d'interpolation autorisées ;
- seuils de pruning LD pour les analyses genome-wide ;
- comparaison avec et sans apparentés éloignés selon la question.

Chaque scénario possède sa propre signature et ne remplace jamais le résultat
principal.

**Visualisations** : tableau de stabilité, forêt des estimations et carte des
scénarios qui changent la conclusion technique.

### Étape 18 — Visualisations consolidées

**Orchestrateur cible** : `visualization/build_figure_index.py`.

Les scripts de visualisation par domaine sont indépendants des scripts de
calcul :

- `visualization/plot_source_qc.py` ;
- `visualization/plot_genotype_qc.py` ;
- `visualization/plot_target_variant_qc.py` ;
- `visualization/plot_kinship_panel.py` ;
- `visualization/plot_kinship.py` ;
- `visualization/plot_population_structure.py` ;
- `visualization/plot_cohorts.py` ;
- `visualization/plot_target_map.py` ;
- `visualization/plot_phasing.py` ;
- `visualization/plot_founder_haplotypes.py` ;
- `visualization/plot_local_ld.py` ;
- `visualization/plot_roh.py` ;
- `visualization/plot_variant_age.py` ;
- `visualization/plot_sensitivity.py`.

Règles :

- un script graphique ne lit jamais directement les exports ACPA ;
- il lit uniquement des sorties validées et leurs audits ;
- il échoue si l'effectif ou les unités attendues ne correspondent pas ;
- il ne modifie aucun résultat scientifique ;
- il ne masque pas les valeurs manquantes ou exclusions ;
- il produit une figure de remplacement explicite « résultat indisponible »
  uniquement pour un rapport incomplet, jamais une ancienne figure ;
- toutes les figures utilisent des légendes générées depuis les paramètres
  réels ;
- les couleurs et catégories sont centralisées dans une configuration graphique.

`figure_index.json` référence toutes les figures, leur provenance, leur niveau
de sensibilité et la section de rapport autorisée.

### Étape 19 — Rapport, interprétation et validation finale

**Scripts cibles** :

- `reporting/build_interpretation_facts.py` ;
- `reporting/render_report.py` ;
- `reporting/validate_report.py`.

**Principe** : l'interprétation automatique se limite à des faits calculés et à
des limitations structurées.

Exemples autorisés :

- « KING ne détecte aucune paire au-dessus du seuil configuré dans la cohorte
  analysée » ;
- « cinq unités porteuses indépendantes possèdent des longueurs gauche et droite
  valides » ;
- « l'estimation est sensible au retrait d'une famille ».

Exemples interdits sans validation humaine :

- « le variant cible provient certainement d'un fondateur unique » ;
- « le variant cible date de trois à cinq générations » si cette phrase n'est pas
  produite depuis le résultat courant ;
- « un ROH commun prouve un effet fondateur ».

Le rapport distingue :

1. données et provenance ;
2. contrôles qualité ;
3. décisions de cohorte ;
4. résultats principaux ;
5. analyses secondaires ;
6. sensibilité ;
7. limites ;
8. interprétation proposée ;
9. validation scientifique manuelle.

Chaque section indique son état : `VALID`, `EXPLORATORY`, `BLOCKED`, `SKIPPED`
ou `FAILED`.

Deux modes de rapport sont distingués :

- `TECHNICAL_INCOMPLETE` : état du run, audits disponibles, échecs et artefacts
  manquants, sans conclusion scientifique ;
- `READY_FOR_SCIENTIFIC_REVIEW` : rapport complet respectant tous les critères
  ci-dessous.

Le rapport final ne peut être marqué `READY_FOR_SCIENTIFIC_REVIEW` que si :

- toutes les étapes critiques ont réussi ;
- toutes les entrées possèdent une provenance ;
- les cohortes sont figées ;
- les figures correspondent au run courant ;
- aucune conclusion codée en dur n'est présente ;
- les limitations obligatoires sont affichées ;
- le validateur de rapport confirme les empreintes et effectifs.

## 10. Contrat CLI commun des scripts

Chaque script d'étape expose au minimum :

```text
python -m effet_fondateur.stages.<stage> \
  --config data/runs/<run_id>/config.resolved.yaml \
  --run-dir data/runs/<run_id> \
  --input-manifest data/runs/<run_id>/stages/<stage_id>/stage_inputs.json \
  --output-dir <temporary_stage_output>
```

Options communes :

- `--validate-only` : valider sans calcul ;
- `--log-level` : niveau de journalisation ;
- `--threads` : ressources autorisées ;
- `--seed` : obligatoire pour toute étape aléatoire ;
- `--tool-path` : fourni par configuration, jamais codé en dur ;
- `--force` : interdit dans un run publié, réservé à un dossier temporaire.

Chaque script écrit :

- journaux techniques sur stderr ou fichier de log ;
- résumé utilisateur non sensible sur stdout ;
- résultats uniquement dans `--output-dir` ;
- `audit.json` avant publication atomique.

### 10.1 Matrice contractuelle des scripts

Les noms ci-dessous sont les artefacts principaux. Chaque dossier d'étape
contient en plus `stage_inputs.json`, `stage_outputs.json`, `audit.json`,
`checksums.sha256` et les journaux techniques.

| Étape | Entrées principales | Sorties principales | Visualisation |
|---|---|---|---|
| `00` | configuration utilisateur, schémas | `config.resolved.yaml`, `environment.json`, `manifest.json` | aucune |
| `01` | répertoires sources déclarés | `source_inventory.tsv`, `source_qc.tsv` | `plot_source_qc.py` |
| `02` | inventaire validé, métadonnées utilisateur, variant cible déclaré | `samples.master.tsv`, `sample_registry_review.tsv` | aucune |
| `03` | exports ACPA validés, `samples.master.tsv` | `genomewide_base.*`, `target_chromosome_base.*`, `sample_alignment.tsv`, audits variants | `plot_genotype_qc.py` |
| `04` | `target_chromosome_base.*`, métadonnées et génotypes cibles | `target_variant.*`, `target_genotype_audit.tsv`, `mendel.tsv` | `plot_target_variant_qc.py` |
| `05` | `genomewide_base.*`, lots et registre | `genomewide_pre_qc.*`, `sample_qc.tsv`, `variant_qc_preliminary.tsv` | `plot_genotype_qc.py` |
| `06` | `genomewide_pre_qc.*`, paramètres pruning | `kinship_panel.*`, `pruned_variants.tsv`, `panel_coverage.tsv`, `panel_residual_ld_bins.tsv` | `plot_kinship_panel.py` |
| `07` | `kinship_panel.*`, pedigree déclaré | `kinship_pairs.tsv`, `kinship_degree_summary.tsv`, `pedigree_concordance.tsv`, `independent_set_proposal.tsv` | `plot_kinship.py` |
| `08` | panel genome-wide, résultat KING, lots et groupes descriptifs | `population_scores.tsv`, `population_eigenvalues.tsv`, `population_outliers.tsv` | `plot_population_structure.py` |
| `09` | registre, variant cible, KING, structure, QC | `cohorts.frozen.tsv`, fichiers `keep`, `cohort_decisions.tsv` | `plot_cohorts.py` |
| `10` | jeux de base, cohortes figées | jeux QC par cohorte, `sample_qc_final.tsv`, `variant_qc.tsv` | `plot_genotype_qc.py` |
| `11` | jeu du chromosome cible QC, variant cible, carte génétique | `target_region.*`, `target_genetic_map.tsv`, formats de phasage | `plot_target_map.py` |
| `12` | région cible, pedigree, génotypes cibles | haplotypes phasés, `carrier_haplotypes.tsv`, `phasing_qc.tsv` | `plot_phasing.py` |
| `13` | haplotypes phasés, porteurs indépendants, témoins | `founder_segments.tsv`, `founder_consensus.tsv`, matrice de partage | `plot_founder_haplotypes.py` |
| `14` | `variant_age_input.tsv`, paramètres de méthode | `variant_age_estimates.tsv`, `variant_age_scenarios.tsv` | `plot_variant_age.py` |
| `15` | jeux de la région cible comparables par cohorte | tables LD natives, `local_ld_summary.tsv` | `plot_local_ld.py` |
| `16` | jeux QC et cohortes | `roh_segments.tsv`, `roh_burden.tsv`, `target_in_roh.tsv` | `plot_roh.py` |
| `17` | résultats principaux et secondaires activés | `sensitivity_scenarios.tsv`, `sensitivity_results.tsv` | `plot_sensitivity.py` |
| `18` | audits et figures validées | `figure_index.json`, contrôle de complétude | index uniquement |
| `19` | manifest, audits, faits, index des figures | rapport HTML/PDF, `report_validation.json` | assemblage uniquement |

Le suffixe `.*` d'un jeu PLINK désigne un artefact composé BED/BIM/FAM et ses
descripteurs, pas un glob résolu par les scripts.

### 10.2 Validation d'entrée obligatoire dans chaque script

Avant tout calcul, le script :

1. valide `stage_inputs.json` ;
2. recalcule les empreintes ;
3. refuse tout chemin non déclaré ;
4. valide le schéma de chaque table ;
5. vérifie assemblage, `sample_set_id` et `variant_set_id` ;
6. vérifie que les effectifs concordent avec l'audit producteur ;
7. vérifie que les dépendances sont `SUCCEEDED` ou `CACHED` ;
8. enregistre les validations dans son audit ;
9. s'arrête avant de créer une sortie si un contrôle critique échoue.

Les scripts graphiques appliquent exactement la même validation aux artefacts
scientifiques qu'ils consomment.

## 11. Pseudo-code de l'orchestrateur

```text
fonction run_pipeline(config_path):
    config = valider_configuration(config_path)
    run = exécuter_bootstrap_etape_00(config)
    dag = construire_dag_etapes_01_a_19(config)

    pour chaque etape dans ordre_topologique(dag):
        si etape désactivée:
            enregistrer SKIPPED avec justification
            continuer

        si dépendance FAILED ou BLOCKED:
            enregistrer BLOCKED
            continuer

        entrees = résoudre_artefacts(etape)
        valider_schemas_entree(etape, entrees)
        signature = calculer_signature(etape, entrees, config, outils)

        si cache_strictement_valide(signature):
            publier référence CACHED
            continuer

        enregistrer événement stage_started
        créer dossier temporaire dédié

        résultat = lancer_sous_processus_sans_shell(etape, arguments)

        si résultat.code_retour != 0:
            collecter audit d'échec
            publier FAILED ou BLOCKED selon le code
            supprimer uniquement le dossier temporaire validé
            continuer selon politique d'arrêt

        valider_schemas_sortie(etape)
        valider_cohérence_effectifs_et_identifiants(etape)
        calculer_empreintes_sortie(etape)
        valider_audit_etape(etape)
        publier sorties atomiquement

        si visualisation associée:
            lancer script graphique sur sorties publiées
            valider figure et figure.json

        mettre à jour manifest atomiquement
        enregistrer événement stage_succeeded

    calculer statut global

    si étapes critiques valides:
        construire faits d'interprétation
        rendre rapport
        valider rapport
    sinon:
        rendre uniquement rapport technique incomplet si demandé

    retourner statut global du run
```

## 12. Tests et critères d'acceptation

### 12.1 Tests unitaires

- validation de chaque schéma ;
- conversion d'unités bp/kb/cM ;
- calcul des signatures ;
- transitions d'état ;
- sélection d'un ensemble indépendant ;
- classification des degrés KING ;
- intersection correcte de segments ;
- séparation positions absolues/longueurs ;
- impossibilité d'inférer un génotype depuis un groupe ;
- génération des faits d'interprétation depuis des audits synthétiques.

### 12.2 Tests d'intégration

- outils externes simulés avec sorties réalistes minimales ;
- reprise après échec ;
- détection d'une sortie ancienne ou modifiée ;
- publication atomique ;
- run complet sur données synthétiques non sensibles ;
- absence de lecture depuis un autre run ;
- cohérence des effectifs entre étapes.

### 12.3 Tests scientifiques de référence

- pedigree synthétique avec degrés connus ;
- cohorte avec haplotype fondateur simulé et limites connues ;
- cohorte sans effet fondateur ;
- variant cible sur plusieurs haplotypes incompatibles ;
- cas de phasage incertain ;
- ROH multiples par individu ;
- carte génétique non linéaire ;
- datation vérifiée contre un résultat de référence documenté.

### 12.4 Tests des visualisations

- figure créée seulement depuis les sorties du run ;
- fichier `.figure.json` obligatoire ;
- effectif affiché égal à l'audit ;
- pseudonymisation active ;
- absence de figure obsolète après échec ;
- rendu minimal vérifié sur données synthétiques.

## 13. Migration et archivage de l'ancien pipeline

L'archivage est une phase séparée de la réécriture.

### 13.1 Inventaire préalable

Pour chaque script actuel, enregistrer :

- fonctions publiques ;
- appels depuis le dépôt ;
- tests existants ;
- sorties produites ;
- comportement scientifique utile ;
- défauts connus ;
- remplacement V2 prévu ;
- décision : conserver, refactoriser, remplacer ou archiver.

### 13.2 Règles d'archivage

- ne jamais archiver un script encore importé ;
- ne jamais déplacer une logique scientifique sans test de référence ;
- ne pas supprimer les scripts archivés dans la même opération ;
- placer les scripts historiques dans un dossier non importable ;
- ajouter un manifest d'archive avec date, raison et remplacement ;
- interdire à l'orchestrateur V2 de lire les sorties historiques ;
- conserver l'historique Git comme source principale de traçabilité.

### 13.3 Contenu du manifest d'archive

```text
ancien chemin
nouveau chemin d'archive
date d'archivage
raison
fonctionnalité remplacée par
tests de remplacement
limitations historiques
```

## 14. Décisions encore requises avant implémentation

1. outil et stratégie de phasage ;
2. méthode de détection IBD locale ;
3. carte génétique GRCh38 et licence de redistribution ;
4. seuils QC adaptés à la taille finale des cohortes ;
5. seuil d'apparentement retenu pour chaque analyse ;
6. stratégie d'unité indépendante par famille et chromosome porteur ;
7. méthode Gamma exacte et implémentation de référence ;
8. méthode alternative de datation ;
9. politique de pseudonymisation et niveau de partage des rapports ;
10. format final du rapport scientifique et procédure de validation humaine.

Ces décisions doivent apparaître dans un registre de décisions versionné. Elles
ne doivent pas être résolues implicitement pendant le codage.

## 15. Définition de terminé pour la V2

La V2 est considérée prête pour une première analyse réelle lorsque :

- le run est entièrement configurable ;
- les sources sont immuables ;
- les jeux genome-wide et du chromosome cible sont distincts et audités ;
- KING utilise un panel autosomique indépendant ;
- les cohortes reposent sur le génotype cible explicite ;
- le pedigree et l'indépendance sont pris en compte ;
- le phasage et l'haplotype fondateur disposent de tests scientifiques ;
- Gamma ne reçoit que des longueurs gauche/droite en cM ;
- toutes les figures possèdent une provenance ;
- le rapport ne peut pas utiliser une ancienne sortie ;
- aucun échec critique n'est silencieux ;
- un run synthétique de référence est reproductible ;
- les décisions scientifiques en attente sont visibles et bloquantes.

## 16. Structure cible optimale du dépôt

Cette structure est une cible de migration. Elle ne doit pas être créée en bloc
avant que les premiers modules V2 soient définis et testés.

```text
Effet_fondateur2/
├── AGENTS.md
├── README.md
├── SESSION.md
├── PIPELINE_V2_PRECODE.md
├── pyproject.toml
├── requirements.txt
├── config/
│   ├── pipeline.example.yaml
│   ├── plotting.yaml
│   ├── studies/
│   │   └── dock6.example.yaml
│   └── decisions/
│       ├── scientific_decisions.yaml
│       └── thresholds.yaml
├── docs/
│   ├── architecture.md
│   ├── data_contracts.md
│   ├── scientific_methods.md
│   ├── audit_and_reproducibility.md
│   ├── migration_plan.md
│   └── modules/
│       ├── acpa.md
│       ├── plink.md
│       ├── king.md
│       ├── phasing.md
│       ├── haplotypes.md
│       ├── ld.md
│       ├── roh.md
│       └── dating.md
├── schemas/
│   ├── pipeline_config.schema.json
│   ├── run_manifest.schema.json
│   ├── stage_audit.schema.json
│   ├── figure_provenance.schema.json
│   ├── samples_master.schema.json
│   └── cohorts_frozen.schema.json
├── src/
│   └── effet_fondateur/
│       ├── __init__.py
│       ├── cli.py
│       ├── orchestrator/
│       │   ├── __init__.py
│       │   ├── pipeline.py
│       │   ├── dag.py
│       │   ├── catalog.py
│       │   ├── runner.py
│       │   ├── signatures.py
│       │   ├── integrity.py
│       │   ├── environment.py
│       │   ├── state.py
│       │   ├── cache.py
│       │   └── recovery.py
│       ├── audit/
│       │   ├── __init__.py
│       │   ├── manifest.py
│       │   ├── events.py
│       │   ├── checksums.py
│       │   └── provenance.py
│       ├── contracts/
│       │   ├── __init__.py
│       │   ├── configuration.py
│       │   ├── documents.py
│       │   ├── tables.py
│       │   ├── artifacts.py
│       │   ├── samples.py
│       │   ├── cohorts.py
│       │   ├── genetics.py
│       │   └── stage_io.py
│       ├── external/
│       │   ├── __init__.py
│       │   ├── plink.py
│       │   ├── king.py
│       │   ├── bcftools.py
│       │   ├── rscript.py
│       │   ├── phasing.py
│       │   └── local_ibd.py
│       ├── stages/
│       │   ├── __init__.py
│       │   ├── initialize_run.py
│       │   ├── validate_sources.py
│       │   ├── build_sample_registry.py
│       │   ├── convert_acpa.py
│       │   ├── prepare_target_variant_dataset.py
│       │   ├── qc_preliminary.py
│       │   ├── build_kinship_panel.py
│       │   ├── infer_kinship.py
│       │   ├── analyze_population_structure.py
│       │   ├── freeze_cohorts.py
│       │   ├── qc_final.py
│       │   ├── prepare_target_region.py
│       │   ├── phase_target_region.py
│       │   ├── infer_founder_haplotype.py
│       │   ├── estimate_variant_age.py
│       │   ├── analyze_local_ld.py
│       │   ├── analyze_roh.py
│       │   └── run_sensitivity_analyses.py
│       ├── visualization/
│       │   ├── __init__.py
│       │   ├── common.py
│       │   ├── plot_source_qc.py
│       │   ├── plot_genotype_qc.py
│       │   ├── plot_target_variant_qc.py
│       │   ├── plot_kinship_panel.py
│       │   ├── plot_kinship.py
│       │   ├── plot_population_structure.py
│       │   ├── plot_cohorts.py
│       │   ├── plot_target_map.py
│       │   ├── plot_phasing.py
│       │   ├── plot_founder_haplotypes.py
│       │   ├── plot_local_ld.py
│       │   ├── plot_roh.py
│       │   ├── plot_variant_age.py
│       │   ├── plot_sensitivity.py
│       │   └── build_figure_index.py
│       └── reporting/
│           ├── __init__.py
│           ├── build_interpretation_facts.py
│           ├── render_report.py
│           ├── validate_report.py
│           └── templates/
│               ├── report.html.j2
│               └── report.css
├── scripts/
│   ├── run_pipeline_v2.py
│   ├── validate_run.py
│   └── inspect_audit.py
├── tests/
│   ├── unit/
│   │   ├── orchestrator/
│   │   ├── audit/
│   │   ├── contracts/
│   │   ├── stages/
│   │   ├── visualization/
│   │   └── reporting/
│   ├── integration/
│   │   ├── test_pipeline_recovery.py
│   │   ├── test_external_tool_failures.py
│   │   └── test_complete_synthetic_run.py
│   ├── scientific/
│   │   ├── test_known_pedigree.py
│   │   ├── test_founder_haplotype_reference.py
│   │   ├── test_no_founder_effect.py
│   │   ├── test_roh_intersections.py
│   │   └── test_variant_age_reference.py
│   └── fixtures/
│       ├── synthetic_acpa/
│       ├── synthetic_pedigrees/
│       ├── synthetic_haplotypes/
│       └── tool_outputs/
├── data/
│   ├── input/
│   │   ├── acpa/
│   │   │   ├── samples/
│   │   │   └── controls/
│   │   ├── metadata/
│   │   ├── target_variant/
│   │   ├── references/
│   │   │   ├── dbsnp/
│   │   │   └── genetic_maps/
│   │   └── README.md
│   ├── runs/
│   │   └── <run_id>/
│   │       ├── manifest.json
│   │       ├── events.jsonl
│   │       ├── config.resolved.yaml
│   │       ├── environment.json
│   │       ├── stages/
│   │       │   ├── 00_initialize_run/
│   │       │   ├── 01_validate_sources/
│   │       │   ├── 02_sample_registry/
│   │       │   ├── 03_acpa_conversion/
│   │       │   ├── 04_target_variant/
│   │       │   ├── 05_qc_preliminary/
│   │       │   ├── 06_kinship_panel/
│   │       │   ├── 07_kinship/
│   │       │   ├── 08_population_structure/
│   │       │   ├── 09_cohorts/
│   │       │   ├── 10_qc_final/
│   │       │   ├── 11_target_region/
│   │       │   ├── 12_phasing/
│   │       │   ├── 13_founder_haplotype/
│   │       │   ├── 14_variant_age/
│   │       │   ├── 15_local_ld/
│   │       │   ├── 16_roh/
│   │       │   ├── 17_sensitivity/
│   │       │   ├── 18_visualizations/
│   │       │   └── 19_report/
│   │       └── report/
│   │           ├── figures/
│   │           ├── figure_index.json
│   │           ├── report.html
│   │           ├── report.pdf
│   │           └── report_validation.json
│   └── README.md
├── archive/
│   ├── README.md
│   ├── archive_manifest.tsv
│   └── legacy_pipeline_v1/
│       ├── scripts/
│       ├── documentation/
│       └── tests/
└── tools/
    └── README.md
```

### 16.1 Règles associées à cette structure

- `src/` contient la logique importable et testable ;
- `scripts/` ne contient que de minces points d'entrée utilisateur ;
- `schemas/` est versionné avec le code ;
- `config/` contient des exemples et décisions, jamais des secrets ;
- `data/input/` est en lecture seule pour le pipeline ;
- `data/runs/` contient uniquement des résultats reproductibles ;
- `archive/` n'est jamais ajouté au chemin d'import Python ;
- un seul dossier `tests/` remplace progressivement la coexistence actuelle de
  `test/` et `tests/` ;
- les anciens résultats ne sont pas déplacés automatiquement dans la nouvelle
  structure ;
- aucune donnée génétique réelle ne doit être commitée.

## 17. Première séquence de mise en œuvre recommandée

1. valider ce document et enregistrer les décisions ouvertes ;
2. créer les schémas du manifest, de l'audit et des métadonnées ;
3. implémenter l'orchestrateur minimal avec une étape factice ;
4. migrer la validation ACPA et la table maître ;
5. construire la conversion genome-wide en une lecture efficace des puces ;
6. migrer le QC préliminaire et KING avec tests scientifiques ;
7. figer la logique de cohortes ;
8. migrer le jeu du chromosome cible et l'injection contrôlée du variant ;
9. choisir puis prototyper le phasage et l'IBD local ;
10. implémenter les visualisations auditées ;
11. implémenter la datation seulement après validation des longueurs ;
12. construire le nouveau rapport ;
13. établir la carte de remplacement des scripts V1 ;
14. archiver les scripts inutiles après validation du run synthétique V2.

L'étape technique `T00_synthetic_stage` peut être utilisée pendant la mise en
œuvre pour tester les contrats, les échecs et la reprise sans donnée génétique.
Elle n'appartient pas au DAG scientifique `00`–`19` et doit rester désactivée
dans les profils d'étude réels.
