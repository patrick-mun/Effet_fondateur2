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

**Dépendances directes** : étapes `02`, `06` et `07`.

**Responsabilité** : détecter les axes de structure et les individus aberrants
sur le panel genome-wide indépendant.

**Stratégie** :

- exporter temporairement les dosages du panel `06` avec PLINK ;
- estimer les axes uniquement sur les individus `PROPOSED_INCLUDED=true` de la
  proposition indépendante `07` ;
- centrer chaque dosage sur `2p` et le réduire par
  `sqrt(2p(1-p))`, avec `p` estimé dans cette référence ;
- imputer les génotypes manquants par `2p`, soit zéro après standardisation ;
- exclure du modèle les variants monomorphes ou insuffisamment appelés dans la
  référence, avec un code explicite ;
- ajuster la PCA par SVD, orienter chaque axe de façon déterministe et projeter
  ensuite tous les individus apparentés sans recalculer les axes ;
- publier fréquences, écarts-types et loadings nécessaires à la reproduction de
  la projection ;
- documenter qu'aucune référence externe n'est introduite dans cette version ;
- distinguer PCA exploratoire et décision d'exclusion ;
- réserver la DAPC à une question définie, avec génotypes correctement encodés,
  nombre de composantes validé et contrôle du surapprentissage.

La métrique d'outlier est la somme des carrés des scores divisés par leur valeur
propre sur les premiers axes configurés. Sa p-value du chi-deux est une
approximation exploratoire, pas un critère clinique ou une exclusion automatique.

**Sorties** : `population_scores.tsv`, `population_eigenvalues.tsv`,
`population_outliers.tsv`, `population_variant_loadings.tsv`, rapport et audit.

**Critère bloquant** : référence indépendante trop petite, nombre de variants
informatifs insuffisant ou rang PCA nul.

**Visualisations** : scree plot, PCA pseudonymisée par cohorte, lot et famille,
ainsi que projection des individus exclus potentiels.

### Étape 09 — Résolution et gel des cohortes

**Script cible** : `stages/freeze_cohorts.py`.

**Dépendances directes** : étapes `02`, `04`, `05`, `07` et `08`.

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

Le gel publie une ligne par individu et par cohorte, y compris pour les individus
non retenus. Les cohortes indépendantes utilisent au plus un représentant par
famille, sélectionné avec les métriques techniques de l'étape `05`. Les
génotypes porteurs/non-porteurs proviennent exclusivement de l'audit accepté de
l'étape `04`. Les groupes témoins sont une liste explicite de `GROUP_LABEL`
configurée ; ils ne servent jamais à déduire le génotype cible.

Une proposition d'exclusion de `07` ou `08` bloque le gel tant qu'une revue ne
fournit pas l'identifiant de revue, le rôle du réviseur, une date horodatée,
l'empreinte SHA-256 exacte de l'artefact relu et la liste complète des individus
proposés. Cette revue est déclarée dans la configuration d'un nouveau run. Une
approbation accepte les exclusions proposées ; changer les représentants exige
un nouveau calcul amont. Les identifiants individuels restent dans la table
sensible `cohort_decisions.tsv` et ne sont pas recopiés dans l'audit partageable.

**Sorties** : `cohorts.frozen.tsv`, `cohort_summary.tsv`,
`cohort_decisions.tsv`, sept fichiers PLINK `--keep`, rapport et audit.

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

**Ordre des traitements** : extraction par le fichier `--keep` figé, calcul du
missingness individuel, exclusion explicite des individus, recalcul des
métriques variant, exclusion explicite des variants puis publication. Les
différences de lot restent descriptives et ne provoquent aucune exclusion.

**Sorties** : trois jeux BED/BIM/FAM finaux avec descripteurs versionnés,
`sample_qc_final.tsv`, `variant_qc_final.tsv`, `batch_qc_final.tsv`, rapport et
audit. L'exception du variant cible conserve les filtres échoués sans transformer
le génotype en mesure fiable.

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

**Contrat implémenté** : carte TSV à identifiant unique et assemblage explicite,
positions physiques strictement croissantes, cM monotones, interpolation
uniquement entre deux ancres séparées d'au plus le seuil configuré, sans
extrapolation. Le jeu PLINK régional reçoit les cM validés et un manifest
`PLINK_BED_BIM_FAM_WITH_CM` indépendant du futur outil de phasage.

**Sorties** : `target_region.bed/.bim/.fam`, descripteur versionné,
`target_genetic_map.tsv`, `phasing_input_manifest.json`, rapport et audit.

**Visualisations** : densité physique et génétique, taux de recombinaison local
et position du variant cible.

### Étape 12 — Phasage et identification du chromosome porteur

**Script cible** : `stages/phase_target_region.py` avec adaptateur d'outil externe.

**Responsabilité** : produire des haplotypes et identifier celui qui porte le
variant cible chez chaque porteur.

**Contrat retenu en 12.0** : adaptateur `shapeit5_phase_common_rare_v1`, SHAPEIT5
`5.1.1` sous licence MIT, référence phasée obligatoire, autosomes GRCh38, carte
génétique obligatoire et pedigree enfant-père-mère pris en charge. Le scaffold
commun est produit par `phase_common`, puis les variants rares ou absents du
panel sont traités par `phase_rare` afin de ne pas perdre le variant cible.

**Catalogue retenu en 12.1** : panel
`1kg_3202_high_coverage_20220422`, release phasée SNV/INDEL/SV de `3 202`
échantillons 1000 Genomes sur GRCh38, toutes populations et tous échantillons
par défaut. Le catalogue local versionné couvre les 22 autosomes, épingle le
README et le manifeste fournisseur, puis expose les MD5 officiels du VCF et de
son index sans effectuer de téléchargement.

**Cache retenu en 12.2** : entrée adressée par les identifiants et empreintes de
la référence, verrou inter-processus, téléchargement dans un dossier temporaire,
contrôles MD5/SHA-256 puis renommage atomique. Une entrée publiée est en lecture
seule et intégralement revérifiée avant réutilisation. Le mode hors ligne bloque
un cache absent sans accès réseau ; une corruption bloque sans remplacement
automatique.

**Extraction retenue en 12.3** : fenêtre `chrN:START-END` produite par
`bcftools` depuis le cache immuable, en conservant exactement tous les
échantillons dans leur ordre source. Les longueurs de contig GRCh38, le champ
`GT`, l'effectif, l'index tabix, le nombre de variants et l'absence de génotypes
appelés non phasés sont contrôlés avant une publication atomique. Le manifeste
lie la sortie à la release et aux empreintes précises du cache.

**Contrat retenu en 12.4** : décomposition biallélique des SNV/indels de la
référence, conservation stricte des échantillons et harmonisation de chaque
marqueur ACPA par assemblage, chromosome, position, REF et ALT. Les inversions
A1/A2 sont auditées, les corrections de brin ne sont jamais automatiques, les
collisions entre sondes restent explicites et le variant cible peut rester
spécifique à l'étude pour `phase_rare`. Une représentation non minimale ou une
ambiguïté touchant le variant cible bloque avant SHAPEIT5.
Les entrées GRCh38 doivent déjà utiliser une représentation minimale : aucune
FASTA ni correction automatique de gauche-alignement n'est ajoutée à ce
contrat. Un allèle PLINK `0` peut seulement être complété depuis une
correspondance publique unique, sans modifier ni inférer les génotypes.

**Contrat retenu en 12.5** : VCF d'étude bgzip/indexé utilisant les `SAMPLE_ID`
maître, REF canonique imposé avec `--a2-allele`, puis contrôle bcftools de
l'ordre des individus et de chaque identité `CHROM/POS/ID/REF/ALT`. Les variants
communs alimentent `phase_common`; le variant cible est conservé comme
`RARE_TARGET` lorsqu'il est absent du panel. La carte SHAPEIT contient
`pos/chr/cM`. Le pedigree sans en-tête contient `enfant père mère`, emploie les
mêmes identifiants que le VCF et utilise `NA` pour un parent absent. Un pedigree
vide est valide et devra conduire `12.6` à omettre son argument. Les fichiers,
versions d'outils, effectifs, régions et empreintes sont liés dans un manifeste
publié atomiquement. SHAPEIT5 n'est pas lancé à cette sous-étape.

**Contrat retenu en 12.6** : exécution séquentielle de `phase_common` puis
`phase_rare` 5.1.1 avec graine, threads, `Ne`, régions et délais enregistrés.
Les individus, variants, génotypes et erreurs mendéliennes sont comparés avant
et après. Le chromosome porteur `H1/H2/BOTH` est dérivé du GT cible phasé et
doit concorder avec le génotype moléculaire explicite. `--score-singletons`
fournit le PP lorsqu'il existe ; un PP absent ou inférieur au seuil ne disparaît
pas, mais rend l'attribution non fiable, crée une zone d'alerte et impose une
revue manuelle. Les transmissions trio/duo restent explicites sans inventer une
origine parentale ambiguë.

**Contrat retenu en 12.7** : consolidation sans recalcul de douze contrôles
versionnés couvrant intégrité, `AC/AN`, nombres de variants, ordre des individus,
préservation des génotypes, cible, Mendel, transmissions et confiance. Les
effectifs `NONE/H1/H2/BOTH` sont publiés sans identifiant individuel. Un warning
de confiance est conservé dans le résumé final et impose une revue manuelle.
Les deux manifestes sources et toutes les sorties sont liés par SHA-256, puis le
dossier QC est publié atomiquement. Les visualisations pseudonymisées sont
déclarées mais restent produites à l'étape 18.

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

**Contrat retenu** : la méthode primaire `target_centered_exact_ibs_v1` compare
les haplotypes porteurs fiables d'unités indépendantes dans le BCF final phasé.
Elle part du variant cible et s'étend de chaque côté tant que tous les chromosomes
porteurs ont un allèle appelé identique ; le premier manque ou la première
discordance arrête le segment. Les limites sont reportées en bp et en distances
cM depuis la cible. Par défaut, trois unités indépendantes et deux marqueurs
informatifs sur chaque flanc sont requis. Le résultat est nommé
`IBS_SHARED_CANDIDATE` : en l'absence d'un adaptateur IBD local séparément validé
sur des données assez denses, il ne constitue jamais une preuve IBD. Le variant
cible est retiré de la signature mesurée chez les témoins et non-porteurs afin
d'éviter une comparaison trivialement déterminée par leur statut mutationnel.
Les porteurs de phase non fiable sont exclus mais conservés dans l'audit. Un
homozygote cible compte comme une unité familiale et n'est informatif qu'aux
marqueurs homozygotes identiques, conformément au cas récessif de Gamma. Toute
insuffisance d'effectif ou de marqueurs produit `NO_FOUNDER_CONCLUSION` et impose
une revue.

**Traitements** :

1. sélectionner uniquement les chromosomes porteurs validés ;
2. comparer les allèles phasés autour du variant cible ;
3. détecter les segments partagés avec une méthode locale choisie et versionnée ;
4. distinguer IBS et IBD selon les informations disponibles ;
5. mesurer pour chaque unité indépendante la limite gauche et droite ;
6. convertir les limites en bp et cM ;
7. comparer la fréquence du segment chez témoins et non-porteurs ;
8. calculer un intervalle partagé central et les segments pairwise nécessaires
   aux variantes de sensibilité ;
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

**Contrat retenu** : l'estimateur primaire est
`gamma_gandolfo_2014_v1`, réimplémentation auditée des formules et corrections
de bord de Gandolfo, Bahlo et Speed (2014). Le modèle corrélé est primaire afin
de tenir compte de l'apparentement cryptique plausible dans une population
insulaire ; le modèle indépendant est toujours publié comme sensibilité. Trois
unités familiales indépendantes permettent seulement une estimation
exploratoire et cinq sont requises pour une estimation primaire. En dessous de
trois unités, ou si un bras individuel est nul, l'étape publie explicitement
`NOT_ESTIMATED`. La correction de partage fortuit est désactivée par défaut :
elle ne peut être activée qu'avec fréquence allélique médiane, nombre de
marqueurs chromosomiques et longueur chromosomique en cM explicitement fournis.
Les sensibilités comprennent le modèle indépendant, le leave-one-family-out et
les conversions à 25, 28 et 30 ans par génération. EstiAge reste différé tant
qu'un adaptateur local reproductible n'est pas validé.

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
- calcul de Gamma corrélé et indépendant avec paramètres explicites ;
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
l'apparentement, l'IBD local, l'haplotype fondateur ou la datation.

**Contrat retenu** : la méthode `plink19_local_ld_secondary_v1` est une analyse
descriptive secondaire. La cohorte principale est `controls_unrelated`. Les
cohortes `target_carriers_independent` et `family_noncarriers` sont analysées
séparément lorsqu'elles existent et satisfont les seuils d'effectif ; elles ne
sont jamais fusionnées. Une unité indépendante ne peut apparaître qu'une fois
par cohorte. Le même univers de variants bialléliques, le même ordre allélique,
la même fenêtre et les mêmes limites de distance sont utilisés pour toute
comparaison. Cet univers est l'intersection du jeu local QC de l'étape `10` et
de la région cartographiée de l'étape `11`, sans sélection depuis les résultats
des étapes `13–14`.

PLINK 1.9 exécute deux calculs natifs distincts. `--r2` fournit le `r²` primaire
de corrélation des dosages alléliques. `--r2 dprime` fournit `D′` absolu par
estimation haplotypique de maximum de vraisemblance et son propre `r²`, conservé
comme mesure secondaire explicitement nommée. L'étape ne mélange jamais ces
deux définitions. Un export PLINK temporaire sert uniquement à compter les
génotypes appelés par paire ; aucune statistique de LD n'est recalculée en
Python.

Les seuils par défaut sont une politique exploratoire de projet, pas des
constantes biologiques : moins de `5` individus produit `NOT_EVALUATED`, de `5`
à `19` individus produit `EXPLORATORY_SMALL_N`, et au moins `20` individus est
requis pour le statut descriptif principal. Chaque paire exige au moins `5`
individus appelés et une MAF d'au moins `0,05` pour chaque variant dans la
cohorte. Les paires qui ne
satisfont pas ces gardes restent auditables avec une raison contrôlée ; une
valeur absente n'est jamais remplacée par zéro. Les classes de distance sans au
moins `10` paires évaluables ne reçoivent pas de résumé numérique.

**Traitements** :

- calcul natif séparé du `r²` génotypique et de `D′` par PLINK 1.9 ;
- restriction configurable en bp et cM, sans extrapolation de carte ;
- univers de variants commun et conventions alléliques identiques entre les
  cohortes effectivement comparées ;
- comptage des génotypes appelés par paire et signalement des monomorphes ;
- statut explicite selon l'effectif et le nombre de paires informatives ;
- signalement séparé des paires impliquant le variant cible, car le recrutement
  par génotype conditionne mécaniquement leur distribution chez les porteurs et
  non-porteurs ;
- aucune pondération familiale a posteriori : seules les cohortes indépendantes
  figées à l'étape `09` sont admises ;
- aucune reconstruction manuelle de `D` depuis `D′` et les MAF ;
- aucune définition de blocs haplotypiques et aucun test naïf sur des millions
  de paires corrélées.

**Articulation avec `13–14`** : l'étape `15` dépend exclusivement des étapes
`09–11` et peut être exécutée en parallèle conceptuel de `13–14`. Elle ne lit ni
`founder_segments.tsv` ni une estimation d'âge, ne modifie aucune limite de
segment et n'alimente jamais Gamma. Elle peut seulement contextualiser, après
coup, la structure de corrélation locale dans laquelle un partage IBS a été
observé. Une concordance ou une différence de LD ne prouve ni ne réfute une
origine fondatrice unique et ne valide pas l'âge estimé.

**Sorties** : tables PLINK natives par cohorte et méthode,
`local_ld_variants.tsv`, `local_ld_pairs.tsv`, `local_ld_summary.tsv`, résumé
JSON agrégé et audit. Les tables par paire sont classées `sensitive_genetic`.

**Critère non bloquant** : une cohorte trop petite, monomorphe ou sans assez de
paires informatives publie `NOT_EVALUATED` avec ses raisons. Puisque l'étape est
secondaire, cette insuffisance ne bloque ni les étapes `13–14` ni le rapport
technique ; elle interdit seulement une interprétation LD pour cette cohorte.

**Visualisations** : matrices de LD comparables, décroissance selon distance en
cM et zoom centré sur le variant cible.

### Étape 16 — ROH secondaire

**Script cible** : `stages/analyze_roh.py`.

**Responsabilité** : étudier l'autozygotie individuelle comme phénomène distinct
de l'haplotype fondateur partagé entre individus.

**Contrat retenu** : la méthode `plink19_array_roh_secondary_v1` utilise PLINK
1.9 avec tous les paramètres de scan explicités. Elle possède deux périmètres
qui ne sont jamais agrégés. Le périmètre `GENOMEWIDE_BURDEN` compare
descriptivement `controls_unrelated` et `target_carriers_independent` sur
l'intersection exacte de leurs variants genome-wide QC de l'étape `10`. Le
périmètre `TARGET_CHROMOSOME` recherche les ROH sur le jeu chromosome cible de
tous les individus QC et détermine, individu par individu, si la coordonnée du
variant cible appartient à un segment. Les `family_noncarriers` peuvent être
décrits dans ce second périmètre, mais pas introduits dans le fardeau genome-wide
faute de jeu autosomique final comparable publié par l'étape `10`.

Le profil exploratoire par défaut requiert un ROH d'au moins `1 500 kb` et `50`
variants, une densité moyenne d'au moins un variant par `50 kb`, aucun trou de
plus de `1 000 kb`, une fenêtre de `50` variants, au plus un hétérozygote et
cinq appels manquants par fenêtre, et un taux minimal de fenêtres concordantes
de `0,05`. Ces valeurs sont configurables et doivent être reprises comme
scénarios à l'étape `17`; elles ne sont pas des constantes biologiques. PLINK
1.9 autorisant par défaut un nombre total illimité d'hétérozygotes dans un
segment, le nombre observé reste audité via `PHET` et toute interprétation est
limitée par le profil de fenêtre. Aucun seuil n'est relâché automatiquement si
la densité ACPA est insuffisante : le périmètre est `NOT_EVALUATED`.

**Traitements** :

- contrôler l'intégrité, l'assemblage, les individus et les variants des trois
  jeux QC de l'étape `10` ;
- construire l'univers genome-wide commun sans fusionner les cohortes ;
- exiger une couverture autosomique et un effectif minimal configurés avant le
  scan genome-wide ;
- exécuter PLINK séparément par cohorte genome-wide avec les mêmes variants et
  paramètres ;
- conserver chaque individu, y compris lorsqu'aucun ROH n'est détecté ;
- publier nombre de segments, longueur totale, longueur maximale et classes de
  longueur par individu ;
- ne calculer un véritable `F_ROH` que si un dénominateur autosomique en kb,
  accompagné de sa source, est explicitement configuré ;
- exécuter séparément le scan du chromosome cible sur tous les individus QC ;
- calculer l'intersection exacte `POS1 <= TARGET_BP <= POS2` pour chaque
  individu et conserver le génotype cible moléculaire accepté ;
- signaler un génotype cible hétérozygote à l'intérieur d'un ROH comme une
  discordance d'interprétation compatible avec la tolérance de scan ou une
  incertitude technique ;
- ne jamais construire une « zone commune » en mélangeant plusieurs segments
  ou chromosomes avant d'avoir défini une règle d'appariement explicite ;
- réserver les comparaisons de groupe aux unités indépendantes et aux effectifs
  suffisants.

**Articulation avec `13–15`** : l'étape `16` ne lit ni les segments IBS de
`13`, ni les estimations Gamma de `14`, ni les paires LD de `15`. Un ROH décrit
une homozygotie continue au sein d'un individu ; l'étape `13` cherche un
haplotype partagé entre individus. L'absence de ROH chez un hétérozygote
n'affaiblit donc pas automatiquement l'hypothèse fondatrice. Les analyses
croisées ROH–LD ou les changements de profil sont réservés à `17`.

**Sorties** : sorties natives PLINK par périmètre, `roh_segments.tsv`,
`roh_burden.tsv`, `target_in_roh.tsv`, `roh_cohort_summary.tsv`, résumé agrégé
et audit. Les tables individuelles sont classées `sensitive_genetic`.

**Critère non bloquant** : effectif, couverture ou densité insuffisants publient
`NOT_EVALUATED` avec une raison contrôlée. L'étape reste secondaire et ne bloque
pas `13–15`; elle interdit seulement une interprétation ROH pour le périmètre
concerné.

**Visualisations** : fardeau ROH, longueurs, pistes individuelles par chromosome
et chevauchements locaux correctement calculés.

### Étape 17 — Analyses de sensibilité

**Script cible** : `stages/run_sensitivity_analyses.py`.

**Responsabilité** : mesurer la stabilité des résultats principaux sans créer
une nouvelle preuve ni remplacer le run primaire.

**Contrat scientifique figé avant implémentation** :

- le run primaire reste la seule analyse principale ; un scénario de
  sensibilité est une perturbation préspécifiée et étiquetée exploratoire ;
- les scénarios qui changent une cohorte, une fenêtre, une carte, un seuil ou
  un pruning sont des runs V2 distincts et immuables, jamais des recalculs
  cachés dans le dossier de l'étape `17` ;
- chaque scénario conserve son `run_id`, son empreinte de configuration, les
  signatures des étapes comparées et les SHA-256 des résumés consommés ;
- un seul axe peut différer du primaire par scénario. Les analyses à plusieurs
  facteurs doivent être déclarées séparément comme exploratoires et ne sont pas
  utilisées pour attribuer une instabilité à un facteur précis ;
- l'assemblage, le variant cible moléculaire, la provenance des génotypes et la
  définition primaire porteur/non-porteur doivent rester identiques ; le statut
  clinique ou le groupe ne peuvent jamais devenir une source de génotype ;
- les unités indépendantes restent la référence. Un scénario avec apparentés
  éloignés est explicitement secondaire et ne remplace aucune estimation
  obtenue sur unités indépendantes ;
- les conclusions sont comparées séparément pour l'IBS local de `13`, la
  datation de `14`, le LD secondaire de `15` et les ROH secondaires de `16` ;
  aucun score composite et aucun vote entre domaines ne sont autorisés ;
- une absence de résultat, un petit effectif ou une étape non exécutée donnent
  `NOT_EVALUATED` et ne sont jamais comptés comme une confirmation ;
- la stabilité catégorielle signifie uniquement que le statut technique du
  domaine ne change pas entre le primaire et tous les scénarios évaluables ;
  elle ne prouve ni IBD, ni origine unique, ni causalité ;
- la variation quantitative est toujours publiée. Elle n'est classée stable ou
  instable que si une tolérance a été préspécifiée dans la configuration ; en
  l'absence de tolérance, le statut est `NOT_CLASSIFIED` ;
- toute inversion de statut, franchissement d'une tolérance ou impossibilité
  systématique d'évaluer un scénario impose une revue manuelle.

**Cohortes et comparabilité** : le primaire utilise les porteurs indépendants
et témoins indépendants gelés à `09`. Le retrait d'une famille ou d'un porteur
s'applique à l'unité indépendante correspondante. Les scénarios avec apparentés
éloignés doivent conserver une colonne d'unité familiale et ne peuvent pas être
interprétés comme une augmentation équivalente de l'effectif indépendant. Les
analyses LD et ROH conservent leurs cohortes et rôles propres définis à `15–16`.

Scénarios minimaux :

- retrait d'une famille à la fois ;
- retrait d'un porteur à faible confiance de phase ;
- plusieurs seuils de définition de l'haplotype ;
- plusieurs fenêtres locales ;
- cartes génétiques ou méthodes d'interpolation autorisées ;
- seuils de pruning LD pour les analyses genome-wide ;
- comparaison avec et sans apparentés éloignés selon la question ;
- profils ROH préspécifiés autour des paramètres exploratoires de `16` ;
- sélection de populations du panel phasé uniquement comme sensibilité, jamais
  comme modification silencieuse du panel primaire.

Un registre externe versionné décrit le run primaire, les runs de sensibilité,
le facteur modifié, les domaines attendus et les empreintes. L'étape contrôle
les manifestes, configurations, audits et résumés, puis publie une ligne par
scénario et domaine. Un scénario absent du registre ne peut pas être ajouté a
posteriori au tableau confirmatoire ; il doit être déclaré exploratoire.

Chaque scénario possède sa propre signature et ne remplace jamais le résultat
principal.

**Sorties** : registre résolu, comparaison des statuts par scénario et domaine,
plages quantitatives disponibles, matrice de stabilité, résumé agrégé et audit.
Le résumé agrégé reste sans identifiant individuel ; le registre résolu est
`sensitive_genetic` s'il contient l'identifiant pseudonymisé d'une unité retirée.

**Limites** : la sensibilité ne corrige ni un biais partagé par tous les runs,
ni une erreur de génotype, ni une carte inadéquate, ni la faible densité ACPA.
Un résultat robuste aux scénarios testés peut rester faux sous une hypothèse non
testée. L'étape ne fournit pas de validation externe sur une autre population
ou une autre technologie.

**Visualisations** : tableau de stabilité, forêt des estimations et carte des
scénarios qui changent la conclusion technique.

### Étape 18 — Visualisations consolidées

**Orchestrateur cible** : `visualization/build_figure_index.py`.

**Responsabilité** : transformer uniquement les artefacts scientifiques déjà
validés du run courant en vues descriptives séparées, sans recalcul, classement
causal ni synthèse quantitative entre domaines. La méthode de consolidation est
versionnée `validated_current_run_figures_v1`.

**Entrées minimales** :

- étape `08` : `population_scores.tsv`, `population_eigenvalues.tsv` et
  `population_outliers.tsv` ;
- étape `13` : `founder_segments.tsv` et `founder_analysis_summary.json` ;
- étape `14` : `variant_age_estimates.tsv` et `variant_age_scenarios.tsv` ;
- étape `15` : `local_ld_summary.tsv` ;
- étape `16` : `roh_cohort_summary.tsv` ;
- étape `17` : `sensitivity_comparisons.tsv` et
  `sensitivity_stability.tsv` ;
- descripteurs de ces artefacts, signatures productrices, SHA-256, schémas et
  audits producteurs déclarés par le run courant.

Les exports ACPA, BED/BIM/FAM, BCF/VCF, sorties PLINK natives et tout chemin
`data/output/` historique sont hors contrat. L'étape ne peut pas découvrir un
fichier par parcours de répertoire ou par glob : chaque entrée doit être
déclarée dans `stage_inputs.json`, appartenir au même `run_id` et provenir de
son étape attendue dans un état `SUCCEEDED` ou `CACHED`.

**Sorties versionnées** :

- `population_structure.svg`, `founder_ibs.svg`, `variant_age.svg`,
  `local_ld.svg`, `roh.svg` et `sensitivity.svg` ;
- `visualization_gallery.html`, rendu de consultation prioritaire assemblé
  depuis `figure_index.json` sans table individuelle ;
- `visualization_gallery.pdf`, rendu secondaire paginé depuis le même index et
  les mêmes SVG avec `fpdf2` ;
- `visualization_render_manifest.json`, qui lie l'index, le HTML et le PDF par
  SHA-256 et atteste l'absence de recalcul scientifique ;
- un document `<figure>.figure.json` par figure, conforme au schéma
  `figure_provenance.schema.json` ;
- `figure_index.json`, conforme au schéma `figure_index.schema.json` ;
- `visualization_completeness.json`, conforme au schéma
  `visualization_completeness.schema.json` ;
- `stage_outputs.json`, `audit.json` et `checksums.sha256` usuels.

Chaque entrée d'index possède un domaine unique parmi `POPULATION_STRUCTURE`,
`FOUNDER_IBS`, `VARIANT_AGE`, `LOCAL_LD`, `ROH` et `SENSITIVITY`, un statut
`RENDERED`, `NOT_EVALUATED` ou `BLOCKED`, son niveau de sensibilité et la liste
de ses sources exactes. Les figures et provenances sont classées au moins
`sensitive_genetic`, même si elles n'affichent que des pseudonymes, afin de ne
pas déclassifier implicitement une vue dérivée de données génétiques.

Le HTML est le format de consultation principal : navigation par domaine,
statuts, limites et liens vers les provenances. Le PDF est produit ensuite et
conserve le même ordre, les mêmes statuts et les mêmes figures. Les deux formats
sont `sensitive_genetic`, ne chargent aucune ressource réseau et ne lisent aucun
artefact scientifique supplémentaire.

**Contrôles d'intégrité bloquants par figure** :

1. chemin déclaré, confiné au run courant et distinct de `data/output/` ;
2. SHA-256, signature productrice, étape productrice et version de schéma ;
3. validation complète des JSON/TSV consommés ;
4. cohérence de `sample_set_id`, `variant_set_id`, assemblage et unités lorsque
   ces champs sont applicables ;
5. égalité des effectifs et statuts redondants entre résumé, table et audit ;
6. aucune valeur numérique présente lorsqu'un statut est `NOT_EVALUATED` ou
   `NOT_ESTIMATED`, et aucune unité bp/kb/cM/générations interchangeable ;
7. absence d'identifiant individuel autre qu'un pseudonyme graphique autorisé.

Une incohérence bloque la figure concernée avant écriture de son SVG et produit
une entrée `BLOCKED` avec un code non sensible. Elle ne doit pas être remplacée
par une figure ancienne. Les autres domaines intègres peuvent être rendus, mais
l'étape 18 reste incomplète et l'étape 19 ne peut pas présenter le rapport comme
prêt pour revue scientifique.

**Comportement `NOT_EVALUATED`** : une non-évaluation valide n'est pas une
erreur d'intégrité. Elle produit un panneau neuf et audité indiquant le domaine,
le statut, les effectifs observés, les valeurs manquantes, les exclusions ou la
raison contrôlée et les limites pertinentes. Une valeur manquante n'est jamais
convertie en zéro et une étape `SKIPPED` reste distincte d'un résultat négatif.

**Visualisations minimales** :

- structure populationnelle : scree plot et projection PC1–PC2 dans un SVG
  commun, avec pseudonymes, distinction de la référence indépendante, des
  apparentés projetés et des outliers exploratoires ; effectifs, composantes
  manquantes et limites sont affichés sans attribution automatique d'ascendance ;
- IBS : longueurs gauche/droite en cM par unité pseudonymisée, inclusions,
  exclusions, confiance manquante, effectif et rappel explicite « IBS, pas
  preuve IBD » ;
- datation : estimation corrélée primaire distincte du modèle indépendant et
  des scénarios exploratoires, intervalles en générations, petits effectifs,
  `NOT_ESTIMATED` et limites Gamma ;
- LD : médianes `r²` et `D′` séparées par cohorte et classe de distance,
  effectifs, paires évaluables et `NOT_EVALUATED`, sans fusion de cohortes ;
- ROH : fardeau genome-wide séparé de l'intersection au chromosome cible,
  effectifs évalués, valeurs manquantes et avertissement autozygotie ≠ IBS/IBD ;
- sensibilités : résultat primaire visuellement distinct des scénarios
  exploratoires, domaines en lignes séparées, scénarios non évalués visibles et
  aucune moyenne, vote ou score composite.

**Articulation avec `08` et `13–17`** : l'étape 18 ne renvoie aucune donnée vers
les étapes scientifiques et ne change aucun de leurs statuts. Elle représente
la PCA exploratoire de `08` sans attribuer d'ascendance, le partage IBS de `13`
sans le renommer IBD, la datation de `14` dans ses unités publiées, le LD de
`15` et les ROH de `16` comme analyses secondaires séparées, et les sensibilités
de `17` comme comparaisons au primaire. Une concordance visuelle entre domaines
ne constitue ni preuve causale, ni preuve d'un fondateur unique, ni validation
de l'âge estimé.

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
- les petits effectifs, exclusions, valeurs manquantes, `NOT_EVALUATED` et
  limites sont inscrits dans la figure et dans sa provenance ;
- aucun script ne produit de score composite d'effet fondateur ni de langage
  causal ;
- le scénario primaire et les scénarios exploratoires utilisent des styles et
  libellés distincts.

`figure_index.json` référence toutes les figures, leur provenance, leur niveau
de sensibilité et la section de rapport autorisée.

### Étape 19 — Rapport, interprétation et validation finale

**Scripts cibles** :

- `reporting/build_interpretation_facts.py` ;
- `reporting/render_report.py` ;
- `reporting/validate_report.py`.

**Principe** : l'interprétation automatique se limite à des faits calculés et à
des limitations structurées.

La méthode de rapport est versionnée `reviewable_scientific_report_v1`. Elle
consomme uniquement les artefacts validés du run courant : les figures, index et
manifest de rendu de `18`, ainsi que les résumés KING, paires pseudonymisées et
contrôles de concordance de `07`. La configuration résolue et les audits du run
sont lus uniquement pour présenter les paramètres effectivement utilisés ; les
chemins d'entrée, secrets, identifiants sources et données génétiques brutes ne
sont jamais injectés dans le rapport ou dans un prompt.

**Sorties de brouillon** :

- `interpretation_facts.json` : faits, effectifs, unités, statuts et limites
  contrôlés, sans conclusion causale ;
- `interpretation_prompt.txt` : prompt versionné construit depuis ce paquet de
  faits et destiné à une IA configurée explicitement ;
- `interpretation_draft.json` : proposition structurée, toujours étiquetée
  `AI_DRAFT` ou `DETERMINISTIC_DRAFT` et jamais considérée comme validée ;
- `report_draft.html` et `report_editor.js` : HTML local éditable, sans ressource
  réseau, qui permet de télécharger un `report_review.json` ;
- `report_validation.json` : état `AWAITING_HUMAN_REVIEW`,
  `READY_FOR_SCIENTIFIC_REVIEW` ou `TECHNICAL_INCOMPLETE` ;
- après relecture : `report_final.html` verrouillé, puis `report_final.pdf`
  produit depuis le même contenu approuvé.

L'appel à une IA externe est désactivé par défaut. Un fournisseur ne peut être
activé que par une configuration explicite et ne reçoit que
`interpretation_facts.json`. Le modèle, sa version, le SHA-256 du prompt et le
mode de génération sont enregistrés ; aucun jeton n'est écrit dans les sorties.
Une génération IA ne modifie jamais un résultat scientifique.

**Cycle de validation humaine** : `AI_DRAFT` → `HUMAN_EDITED` → `APPROVED` →
`FINAL`. Chaque zone de texte conserve son identifiant de section et ses faits
sources. Le bouton « Valider le rapport » contrôle les champs obligatoires et
télécharge la décision de relecture ; la finalisation côté pipeline vérifie son
schéma, les empreintes du brouillon et des faits, les formulations interdites et
la complétude avant de verrouiller le HTML. Une correction ultérieure crée une
nouvelle révision au lieu de réécrire l'approbation précédente.

Le rapport commence par une fiche d'étude issue des paramètres résolus : run,
assemblage, cible, effectifs, seuils QC/KING, paramètres de structure, fenêtres
IBS/LD, profil ROH, datation et sensibilités. Il présente ensuite un tableau KING
pseudonymisé et un réseau de parenté lorsque les artefacts validés le permettent.
KING reste un apparentement genome-wide et n'est jamais décrit comme une preuve
d'IBD locale. La PCA est primaire ; une DAPC ne pourra apparaître que comme
analyse exploratoire si une étape amont publie ultérieurement un artefact DAPC
validé. L'étape `19` ne reprend jamais directement le PNG ou le code historique.

**Garde-fous rédactionnels** : tous les nombres et unités du texte doivent être
présents dans le paquet de faits ; `NOT_EVALUATED` interdit une conclusion ; le
primaire reste distinct de l'exploratoire ; les petits effectifs, exclusions,
valeurs manquantes et limites sont obligatoires. Les termes affirmant une preuve,
une causalité, une origine unique ou une confirmation automatique déclenchent
`REVIEW_REQUIRED`. Aucun score composite d'effet fondateur n'est produit.

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
| `08` | panel genome-wide, résultat KING, lots et groupes descriptifs | `population_scores.tsv`, `population_eigenvalues.tsv`, `population_outliers.tsv`, `population_variant_loadings.tsv` | `plot_population_structure.py` |
| `09` | registre, variant cible, KING, structure, QC | `cohorts.frozen.tsv`, `cohort_summary.tsv`, fichiers `keep`, `cohort_decisions.tsv` | `plot_cohorts.py` |
| `10` | jeux de base, cohortes figées | jeux QC par cohorte, `sample_qc_final.tsv`, `variant_qc_final.tsv` | `plot_genotype_qc.py` |
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

1. méthode de détection IBD locale ;
2. carte génétique GRCh38 et licence de redistribution ;
3. seuils QC adaptés à la taille finale des cohortes ;
4. seuil d'apparentement retenu pour chaque analyse ;
5. stratégie d'unité indépendante par famille et chromosome porteur ;
6. méthode Gamma exacte et implémentation de référence ;
7. méthode alternative de datation ;
8. politique de pseudonymisation et niveau de partage des rapports ;
9. format final du rapport scientifique et procédure de validation humaine.

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
