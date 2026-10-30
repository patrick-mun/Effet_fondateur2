# Calibration par simulation

Ce dossier contient uniquement les éléments reproductibles qui doivent être
versionnés avec le pipeline : plan, protocoles, configurations, schémas et
modèles de rapports.

Les développements suivent obligatoirement
`CONVENTIONS_DEVELOPPEMENT.md`, déclinaison locale de `AGENTS.md` pour la
simulation, l'aveuglement, les tests et le calcul scientifique.

## Séparation des espaces

- `calibration/` : documents et contrats suivis par Git ;
- `calibration/references/` : liste fermée des références publiques D09 ;
- `src/effet_fondateur/calibration/` : futur code de simulation et
  d'évaluation suivi par Git ;
- `tests/calibration/` : futurs tests et minuscules fixtures synthétiques ;
- `calibration_workspace/` : panels, vérités, runs, figures et rapports
  générés, entièrement ignorés par Git.

Les chemins historiques éventuels `calibration/runs/` et
`calibration/local/` sont également ignorés par précaution.

## Règles de sécurité

- Ne jamais copier de donnée de `data/input/dock6_real/` dans ces espaces.
- Ne jamais utiliser un résultat individuel DOCK6 pour régler un seuil.
- Ne jamais committer un VCF, BCF, fichier de vérité ou rapport généré.
- Fournir le chemin du workspace au lancement ; ne pas enregistrer de chemin
  absolu local dans une configuration versionnée.
- Vérifier les graines, versions et SHA-256 avant toute qualification.

## Matérialisation D10

La commande suivante construit un nouveau jeu de partitions à partir des deux
métadonnées publiques vérifiées par D09 :

```bash
.venv/bin/python -m effet_fondateur.calibration.materialize_partitions \
  --public-catalog calibration/references/public_reference_catalog.json \
  --population-metadata <fichier_public_population_3202> \
  --unrelated-index <index_public_2504> \
  --workspace-root calibration_workspace \
  --partition-run-id <identifiant_nouveau>
```

Le dossier cible doit être neuf. La commande n'affiche aucun identifiant et
écrit les listes avec des permissions privées sous
`private_truth/partitions/`.

## Matérialisation D11

La D11a lit les VCF 1000G depuis un chemin public fourni explicitement. Elle ne
copie pas ces références dans le dépôt et n'enregistre aucun chemin absolu dans
le manifest : seuls le nom de fichier et les SHA-256 sont conservés.

```bash
.venv/bin/python -m effet_fondateur.calibration.materialize_base_haplotypes \
  --public-catalog calibration/references/public_reference_catalog.json \
  --technical-catalog config/references/1000g_high_coverage_grch38_phased.json \
  --source-vcf <dossier_public_1000G>/1kGP_high_coverage_Illumina.chr19.filtered.SNV_INDEL_SV_phased_panel.vcf.gz \
  --source-index <dossier_public_1000G>/1kGP_high_coverage_Illumina.chr19.filtered.SNV_INDEL_SV_phased_panel.vcf.gz.tbi \
  --population-metadata <fichier_public_population_3202> \
  --genetic-map <carte_grch38_chr19.tsv> \
  --genetic-map-manifest <manifest_cache_cartes.json> \
  --partition-root calibration_workspace/<partition_d10> \
  --workspace-root calibration_workspace \
  --run-id <identifiant_nouveau> \
  --panel-role CALIBRATION \
  --chromosome 19 \
  --start-bp <début> \
  --end-bp <fin> \
  --random-seed <graine_haplotypes>
```

La D11b lie ensuite cette vérité dense au dessin exact de la puce. Elle lit
directement le paquet constructeur `CytoScan 750K Accel` `NA36.r1` conservé
hors du dépôt et hors du cache. Elle produit elle-même la table SNP normalisée
indépendante de tout échantillon dans le workspace privé :

```bash
.venv/bin/python -m effet_fondateur.calibration.materialize_chip_observation \
  --public-catalog calibration/references/public_reference_catalog.json \
  --dense-root calibration_workspace/<run_d11a> \
  --official-chip-package <paquet_officiel_accel_na36.zip> \
  --expected-official-package-sha256 <sha256_paquet> \
  --workspace-root calibration_workspace \
  --run-id <identifiant_nouveau_d11b>
```

Le format normalisé et les règles d'exclusion sont définis dans
`references/CYTOSCAN_750K_ACCEL_MARKERS.md`. D12 reste bloquée jusqu'à la
réussite d'un run D11b sur le paquet officiel. Le catalogue D09 doit d'abord
recevoir l'entrée `thermofisher_cytoscan750k_accel_na36_r1` avec le SHA-256
réel et une revue des conditions d'usage ; le matérialisateur refuse son
absence.

## Matérialisation D12

D12 réapparie la vérité phasée D11a selon une configuration explicite et
reconstruit séparément le VCF aveugle sur l'univers de marqueurs D11b :

```bash
.venv/bin/python -m effet_fondateur.calibration.materialize_population_structure \
  --d11a-root calibration_workspace/<run_d11a> \
  --d11b-root calibration_workspace/<run_d11b> \
  --population-config calibration/configurations/<configuration_d12>.json \
  --workspace-root calibration_workspace \
  --run-id <identifiant_nouveau_d12>
```

Les proportions inégales ne sont pas acceptées avant D14. Le diagnostic de
100 kb utilise la configuration `d12_diagnostic_100kb_calibration_v1.json` ;
il ne constitue ni un micro-panel, ni une qualification.

## Matérialisation D13

D13 ajoute des familles de cousins germains au panel D12 et transmet des
haplotypes parentaux entiers. La recombinaison reste interdite jusqu'à D14. La
parenté publiée est théorique, dérivée du pedigree ; l'ascendance de chaque
descendant est, elle, recalculée depuis les copies D12 réellement transmises.

```bash
.venv/bin/python -m effet_fondateur.calibration.materialize_pedigrees \
  --d12-root calibration_workspace/<run_d12> \
  --pedigree-config calibration/configurations/<configuration_d13>.json \
  --workspace-root calibration_workspace \
  --run-id <identifiant_nouveau_d13>
```

Le diagnostic de 100 kb utilise
`d13_diagnostic_100kb_calibration_v1.json`. Ses vérités et VCF restent dans le
workspace ignoré ; il ne vaut ni micro-panel, ni pilote, ni qualification.

## Matérialisation D14

D14 simule, sur la carte D11a, les premiers événements de recombinaison autour
d'un variant synthétique d'âge connu. Les trois racines sont fournies pour
vérifier explicitement toute la chaîne D11a→D12→D13 :

```bash
.venv/bin/python -m effet_fondateur.calibration.materialize_founder_recombination \
  --d13-root calibration_workspace/<run_d13> \
  --d12-root calibration_workspace/<run_d12> \
  --d11a-root calibration_workspace/<run_d11a> \
  --founder-config calibration/configurations/<configuration_d14>.json \
  --workspace-root calibration_workspace \
  --run-id <identifiant_nouveau_d14>
```

Les bras hors fenêtre restent explicitement censurés. Le diagnostic de 100 kb
à 150 générations vérifie ce comportement, mais ne qualifie pas la distribution
des longueurs ancestrales ni la datation.

## Matérialisation D14b

D14b dimensionne la fenêtre principale en `cM`, vérifie la cible publique, la
densité CytoScan et les seuils IBD avant D15. Une seconde fenêtre, volontairement
tronquée, reste réservée au contrôle de la censure.

```bash
.venv/bin/python -m effet_fondateur.calibration.materialize_simulation_windows \
  --window-config calibration/configurations/<configuration_d14b>.json \
  --source-vcf <vcf_public_1000G> \
  --source-index <index_public_1000G> \
  --genetic-map <carte_publique_grch38> \
  --normalized-chip-markers <table_snp_normalisee_d11b> \
  --workspace-root calibration_workspace \
  --run-id <identifiant_nouveau_d14b>
```

La chaîne D11a–D14 est ensuite rejouée sur la fenêtre principale. Un reçu neuf
lie les six manifests sans modifier rétroactivement un run finalisé :

```bash
.venv/bin/python -m effet_fondateur.calibration.materialize_window_replay_receipt \
  --d14b-root calibration_workspace/<run_d14b> \
  --d11a-root calibration_workspace/<run_d11a_large> \
  --d11b-root calibration_workspace/<run_d11b_large> \
  --d12-root calibration_workspace/<run_d12_large> \
  --d13-root calibration_workspace/<run_d13_large> \
  --d14-root calibration_workspace/<run_d14_large> \
  --workspace-root calibration_workspace \
  --run-id <identifiant_recu_nouveau>
```

Une fenêtre principale effectivement censurée reste dans la vérité ; elle
n'est jamais rejouée pour obtenir artificiellement un breakpoint interne.

## Matérialisation D16

D16 ajoute des ROH par copie H1→H2. La configuration distingue un témoin
sans injection, une ROH hors cible, une ROH sur la cible sans IBD injecté et
une ROH sur un sujet positif D15. Les H1, les négatifs/matrices IBD et les
transmissions sont préservés. Les vérités IBD et ROH restent séparées. Le
partage supplémentaire H2–H1 créé par le cas couplé est enregistré dans une
vérité IBD induite, à consommer en complément de la vérité D15 inchangée.

```bash
.venv/bin/python -m effet_fondateur.calibration.materialize_roh \
  --d15-root calibration_workspace/<run_d15> \
  --d14-root calibration_workspace/<run_d14> \
  --d11a-root calibration_workspace/<run_d11a> \
  --config-path calibration/configurations/d16_micro_chr5_calibration_v1.json \
  --workspace-root calibration_workspace \
  --run-id <identifiant_neuf_d16>
```

Les trois intervalles de cette configuration mesurent 2 Mb. Les contrôles de
longueur, nombre de sondes, densité et trous entre sondes reprennent le cadre
PLINK documenté ; ils ne mesurent pas encore la sensibilité de détection. Le
témoin négatif signifie absence d'injection, sans exclure une ROH de fond.

```bash
.venv/bin/python -m effet_fondateur.calibration.report_roh_injection \
  --run-root calibration_workspace/<run_d16> \
  --template-path <report_skill>/templates/recherche.html \
  --output-root calibration_workspace/<rapport_d16_neuf>
```

Le rapport `CAL-DEV-D16` reste une validation technique régionale. Les métriques
de détection seront mesurées par l'évaluateur D26 et le rapport `CAL-REP-16`.

## Matérialisation D17

D17 applique huit profils de mesure au même panel D16, avec graine
`measurement_error` explicite. La couche biologique D16 n'est jamais modifiée.
Le registre conserve chaque retrait de sonde, dosage erroné, inversion de
phase, frontière de commutation et appel manquant, avec valeurs avant/après.

```bash
.venv/bin/python -m effet_fondateur.calibration.materialize_measurement_errors \
  --d16-root calibration_workspace/<run_d16> \
  --config-path calibration/configurations/d17_micro_chr5_calibration_v1.json \
  --workspace-root calibration_workspace \
  --run-id <identifiant_neuf_d17>
```

Chaque profil part de la vérité sur sondes, pas du profil précédent. Le VCF
non phasé est distinct du candidat phasé privé, réservé aux évaluateurs de
résistance. Une phase injectée ne doit jamais être remise au pipeline comme
une sortie SHAPEIT. Les fichiers restent privés jusqu'à l'export D18.
Le manque et le dosage portent sur les cellules retenues ; la phase ponctuelle
sur les hétérozygotes ; les commutations sur leurs frontières successives par
sujet/chromosome. Les quotas sont arrondis avec demi supérieur ; les blocs
peuvent inverser davantage de GT qu'il n'existe de frontières. Les défauts
composés peuvent se masquer ou se compenser, sans effacer leur historique.

Le contrôle intact et les sept variantes sont des observations appariées,
pas huit réplications indépendantes. Les niveaux sont des stress techniques,
pas un modèle empirique de l'erreur de la puce. Le retrait de sondes ne modifie
ni les bornes vraies ni les seuils IBD/ROH. La cible est gardée si elle est
déjà observable ; son génotype n'est pas protégé des erreurs.

```bash
.venv/bin/python -m effet_fondateur.calibration.report_measurement_errors \
  --run-root calibration_workspace/<run_d17> \
  --template-path <report_skill>/templates/recherche.html \
  --output-root calibration_workspace/<rapport_d17_neuf>
```

`CAL-DEV-D17` ne mesure pas les performances du pipeline : appels réels et
évaluateurs restent requis. `calibrated` reste `false`. Aucun nouvel outil ni
aucune dépendance ne sont ajoutés.

## Export aveugle D18

D18 consomme uniquement les VCF observés D17 après contrôle de toutes les
empreintes. Les lots et sujets sont pseudonymisés, les ordres permutés et les
correspondances conservées exclusivement dans le contrôle privé.

```bash
.venv/bin/python -m effet_fondateur.calibration.export_blinded_inputs \
  --d17-root calibration_workspace/<run_d17> \
  --config-path calibration/configurations/d18_micro_chr5_calibration_v1.json \
  --workspace-root calibration_workspace \
  --run-id <identifiant_neuf_d18>
```

`blinded_inputs/` contient un inventaire opérateur et des lots opaques. Chaque
lot comporte `observed.vcf.gz`, `acpa_samples/`, `metadata/samples.tsv` et
`pipeline.yaml`. Le VCF est un gzip de contrôle, pas une entrée BGZF indexée ;
la configuration consomme les ACPA. Aucun registre d'erreurs, candidat phasé,
profil nommé ou chemin vers la vérité n'est livré au pipeline.

Le registre indique sexe/clinique UNKNOWN, groupe UNSPECIFIED, filiation et
génotype cible non renseignés. FID est une clé individuelle technique : il
ne prouve pas l'absence de parenté. Aucun diagnostic, annotation moléculaire
CONFIRMED, HGVS ou accord humain n'est fabriqué. Les étapes 03 et suivantes
restent désactivées : des informations d'observation autorisées et une revue
humaine restent nécessaires pour l'analyse familiale complète.

Après activation de `.venv` depuis la racine, l'opérateur peut lancer la
configuration depuis le répertoire du lot ; tous les chemins d'entrée sont
relatifs à ce répertoire :

```bash
cd calibration_workspace/<run_d18>/blinded_inputs/bundles/<lot_opaque>
effet-fondateur run --config pipeline.yaml \
  --runs-dir <chemin_absolu_de_sortie_hors_blinded_inputs>
```

Cette commande, non exécutée automatiquement, est limitée à 00–02 et n'accorde
pas l'approbation du registre. Les sorties doivent rester hors de l'arbre
aveugle immuable ; le répertoire par défaut n'est donc pas employé ici.
La portée chromosome unique est explicite ;
aucun succès genome-wide n'est revendiqué. L'opérateur ne doit recevoir que
l'arbre `blinded_inputs/`, pas le workspace entier ni les graines privées.
L'inventaire fermé refuse fichiers, dossiers et liens supplémentaires.

```bash
.venv/bin/python -m effet_fondateur.calibration.report_blinded_export \
  --run-root calibration_workspace/<run_d18> \
  --template-path <report_skill>/templates/recherche.html \
  --output-root calibration_workspace/<rapport_d18_neuf>
```

Le rapport `CAL-DEV-D18` reste agrégé, sans identifiant de sujet ni
correspondance de profil. Le scellement D06/D07 et les évaluateurs D19 sont
ultérieurs ; aucun outil externe ni aucune dépendance ne sont ajoutés.

## Évaluation technique D19

D19 sépare le rejeu aveugle et l'évaluation. Le runner diagnostique ci-dessous
refuse toute configuration autre que 00–02, tout outil génétique activé et
toute approbation automatique. Prévenir avant ce test d'intégration ; aucun
run d'étude n'est concerné. Les sorties restent hors des entrées D18 immuables.

```bash
.venv/bin/python -m effet_fondateur.calibration.run_technical_inputs \
  --d18-root calibration_workspace/<run_d18> \
  --workspace-root calibration_workspace \
  --run-id <rejeu_d19_neuf>

.venv/bin/python -m effet_fondateur.calibration.evaluate_technical \
  --d18-root calibration_workspace/<run_d18> \
  --run-index-path calibration_workspace/<rejeu_d19_neuf>/run_index.json \
  --workspace-root calibration_workspace \
  --run-id <evaluation_d19_neuve>

.venv/bin/python -m effet_fondateur.calibration.report_technical_evaluation \
  --run-root calibration_workspace/<evaluation_d19_neuve> \
  --template-path <report_skill>/templates/recherche.html \
  --output-root calibration_workspace/<rapports_d19_neufs>
```

Le reçu `technical_run_index` ancre le manifest D18 et, pour chaque lot,
le chemin relatif du run terminé, les empreintes de manifest et
d'environnement, et éventuellement une attente cible privée empreintée.
Le runner livré laisse cette attente absente. L'évaluateur refuse un état
RUNNING avant d'ouvrir la correspondance privée. Il n'exécute rien ; un reçu
absent produit `NOT_EVALUABLE`, une sortie déclarée réussie mais corrompue ou
manquante produit `FAIL`. Les résultats sans mesure portent `null`, pas zéro.

Les cinq fichiers privés `evaluation_00.json` à `evaluation_04.json` sont
scellés par un manifest. `load_checked_evaluation` contrôle leurs populations,
identités, empreintes, motifs et agrégats. Les six rapports HTML autonomes
conservent les statuts par étape et restent ignorés par Git ; leur export PDF
nécessite le contrôle de pagination habituel du skill report-formater.

Le rejeu des huit lots régionaux D18 donne PASS en 00–02 ; 03–04 restent non
exécutées et non évaluables. Les tests séparés vérifient le décodage BED,
l'absence de perte des génotypes/manques, l'ordre, les six champs FAM et la
préservation du fond lors de l'insertion d'une cible explicite. Ils ne
qualifient pas le pipeline génome entier. La cible D18 est déjà parmi les
marqueurs de base ; 04 exige une cible nouvelle et ne doit pas être activée
sans politique explicite, métadonnées adaptées et approbation humaine.
Les paramètres scientifiques et `calibrated=false` restent inchangés.
Les PASS portent sur la concordance D19 : ils ne remplacent pas les critères
transversaux D04/D08, notamment le journal système d'accès non instrumenté
ici, ni le scellement global D06/D07 requis avant qualification.

### Complément D19 : rendre 03 et 04 évaluables

`technical_scenarios.py` prépare un dérivé neuf sans modifier D18. La
configuration `configurations/d19_technical_target_v1.json` fixe une annotation
publique, une règle déterministe et l'approbation technique autorisée par
`CAL-DEC-003`. Ce n'est pas une revue clinique/familiale. Toute modification
de ces registres neutres ou du périmètre exige une nouvelle autorisation ;
le runner n'en fabrique jamais une.

La cible technique publique est APC rs1801155, GRCh38 chr5:112839514 T>A,
NM_000038.6:c.3920T>A, p.Ile1307Lys. Seule sa définition issue de
[ClinVar RCV001195214.19](https://www.ncbi.nlm.nih.gov/clinvar/RCV001195214.19/)
est utilisée : ni classification clinique ni génotypes de vrais patients ne
sont importés. Les confirmations du contrat concernent référence/annotation,
pas un examen de laboratoire des sujets simulés. Cette cible extérieure à la
fenêtre D14–D18 sert uniquement au test d'insertion 04 et ne remplace pas la
cible biologique du panel.

Les rangs du registre aveugle modulo trois reçoivent RR, RA, AA. Cette règle
est indépendante du groupe et de la clinique, sans nouvelle graine ni
prétention à une fréquence réelle. L'attente privée est scellée avant le run ;
seule sa sérialisation comme génotypes explicites entre dans le pipeline.
Les observations D18 préexistantes sont conservées.

Prévenir avant l'intégration réelle : seul PLINK est activé et les étapes
restent bornées à 00–04. Utiliser des identifiants neufs :

```bash
.venv/bin/python -m effet_fondateur.calibration.technical_scenarios \
  --d18-root calibration_workspace/<run_d18> \
  --configuration-path calibration/configurations/d19_technical_target_v1.json \
  --workspace-root calibration_workspace \
  --run-id <entrees_0304_neuves>

.venv/bin/python -m effet_fondateur.calibration.run_technical_scenarios \
  --panel-root calibration_workspace/<entrees_0304_neuves> \
  --workspace-root calibration_workspace \
  --run-id <rejeu_0304_neuf>

.venv/bin/python -m effet_fondateur.calibration.evaluate_technical \
  --d18-root calibration_workspace/<run_d18> \
  --technical-panel-root calibration_workspace/<entrees_0304_neuves> \
  --run-index-path calibration_workspace/<rejeu_0304_neuf>/run_index.json \
  --workspace-root calibration_workspace \
  --run-id <evaluation_0304_neuve>
```

Le scénario par défaut `TARGET_ABSENT_INSERTION` conserve les huit lots. Pour le négatif,
préparer un autre run avec `--scenario TARGET_PRESENT_REJECTION` : un lot
représentatif reçoit un ajout explicite de la cible dans une copie du VCF et
des ACPA. Exécuter et évaluer ce run avec les mêmes commandes et de nouveaux
identifiants. Le pipeline doit rester `FAILED/BLOCKED` en 04 avec code métier
V2=4 ; sa CLI retourne 2. Le contrôle n'est PASS que si la tentative empreintée
porte exactement `target_variant_already_present`, si la base est intacte et
si aucune sortie 04 n'est publiée. Un autre échec n'est jamais accepté.
Ce PASS du garde n'est pas une insertion réussie ; son taux de discordance
d'insertion reste `null`.

Le manifest des entrées ancre les fichiers, la configuration et la source
D18. L'évaluation vérifie aussi la dérivation indépendante des observations et
lie le reçu au bon panel. Générer les rapports avec la commande D19 précédente,
dans un dossier distinct pour chaque scénario.

Sans trio documenté, Mendel reste `NOT_APPLICABLE`. Ces diagnostics régionaux
ne remplacent pas les futurs panels génome entier, le pilote ou la
qualification scientifique ; `calibrated=false` reste inchangé.

## État actuel

D15 injecte une paire positive et une paire négative pour chacune des cinq
longueurs IBD figées dans le manifest D14b. Les 25 sujets et matrices sont
disjoints, hors pedigree et hors porteurs D14. La copie porte uniquement sur
H1 dans l'intervalle ; le fond négatif reste inchangé et son IBD historique
reste inconnu. La vérité enregistre des bornes nulles pour les négatifs.

```bash
.venv/bin/python -m effet_fondateur.calibration.materialize_ibd \
  --d14-root calibration_workspace/<run_d14> \
  --d11a-root calibration_workspace/<run_d11a> \
  --design-root calibration_workspace/<run_d14b> \
  --replay-root calibration_workspace/<recu_d14b> \
  --workspace-root calibration_workspace \
  --run-id <identifiant_neuf_d15>
```

Le manifest conserve les empreintes des quatre entrées, la dérivation de
graine, la vérité privée et les panels dense/observable. Un schéma est validé
avant publication atomique. Ce run régional ne qualifie pas 16C et conserve
`calibrated: false`.

Le rapport autonome de développement peut être produit avec :

```bash
.venv/bin/python -m effet_fondateur.calibration.report_ibd_injection \
  --run-root calibration_workspace/<run_d15> \
  --template-path <report_skill>/templates/recherche.html \
  --output-root calibration_workspace/<rapport_d15_neuf>
```

Les rapports `CAL-DEV-D15` et leurs PDF restent ignorés par Git. Leur statut
technique ne remplace pas le futur rapport de performance `CAL-REP-16C`.

Le périmètre initial D01 et l'arborescence isolée D02 sont établis. Le protocole
D03 est `APPROVED` et les critères D04 sont `APPROVED_PRE_PILOT`. Cette
approbation porte sur le cadre et les règles : les seuils de performance restent
non scellés et seront proposés sur le seul pilote de calibration selon une
règle écrite avant le run, puis figés avant le panel aveugle. Le schéma de
configuration D05, le manifest de vérité D06 et le manifest reproductible de
run D07 sont créés et testés. D06 impose neuf familles de vérité synthétique,
leurs champs, unités, empreintes, chemins privés et consommateurs autorisés.
D07 impose les sources, versions, commandes, graines, fichiers, empreintes et
instructions de rejeu nécessaires pour reproduire une exécution sans
connaissance orale. Les gardes D08 contrôlent les chemins physiques, liens
symboliques ou physiques, références croisées et collisions entre panels avant
toute lecture d'un artefact. D09 fige quatre sources publiques : panel phasé
1kGP, métadonnées de population, index non apparenté et cartes GRCh38
SHAPEIT4/HapMap. Chaque source possède version, origine, assemblage, conditions
d'usage et empreinte SHA-256 ; les charges utiles acquises ultérieurement
devront aussi recevoir un SHA-256 local avant usage. D10 ajoute le contrat et
le validateur des trois partitions : donneurs de calibration, donneurs de
validation et référence de phasage commune. Les groupes de pedigree ou de
parenté ne peuvent pas traverser deux partitions et les listes de membres
restent privées dans le workspace. D11a ajoute le contrat des haplotypes denses,
le réappariement stratifié des copies donneuses, la dérivation des génotypes et
la carte interpolée. D11b ajoute le contrat indépendant de liaison au dessin
CytoScan 750K Accel. Aucun panel complet ni rapport de calibration n'a encore
été produit. Une matérialisation D10 contrôlée existe
sous `calibration_workspace/d10_1kg_3202_20220422_v1/` ; elle contient 900
donneurs de calibration, 899 donneurs de validation et 1 403 références de
phasage, sans intersection d'identifiant ni de groupe. Un test d'intégration
D11a borné à 100 kb sur le chromosome 19 a produit 900 diploïdes synthétiques
et 2 484 variants avec les neuf contrôles denses réussis. Le run D11b borné
`d11b_integration_cal_chr19_11000000_11100000_v1` a ensuite validé le paquet
constructeur réel, normalisé 247 777 SNP autosomiques et lié sans ambiguïté
l'unique marqueur présent dans cette fenêtre. Les manifests D11a et D11b sont
valides et liés au même panel et au même domaine : le prérequis technique de
D12 est satisfait. Cette intégration de 100 kb reste diagnostique et ne
constitue pas le futur micro-panel D33 ni une validation de densité à l'échelle
du génome. D12 a ensuite conservé les 1 800 haplotypes par bijection et produit
20 mélanges CEU/YRI 50/50 parmi 900 diploïdes. Sa vérité dense compte 2 484
variants et son VCF aveugle conserve l'unique sonde D11b sans étiquette
d'ascendance. Les contrats, proportions, empreintes et permissions passent ;
ce diagnostic ne vaut toujours pas qualification. D13 ajoute ensuite deux
familles synthétiques, huit descendants et quatre degrés de relation au sens du
contrat, sans recombinaison. Le run de référence `v3` contrôle la transmission
et l'ascendance haplotypique ; le premier essai `v1` est invalidé et conservé
uniquement pour l'audit. D14 et D14b disposent maintenant d'une intégration
large indépendante de DOCK6 sur le chromosome 5 : 5,991 cM, 165 290 variants
denses, 487 marqueurs observables et cinq scénarios IBD évaluables. Le rejeu
D11a–D14 observe quatre breakpoints internes, sans censure ni erreur
mendélienne, et son reçu conserve toute la chaîne d'empreintes. Cette
intégration valide D14b techniquement mais ne remplace pas le micro-panel D33,
le pilote ni la qualification. Le développement doit suivre
`PLAN_DEVELOPPEMENT_PANEL_SIMULATION.md` et ses portes de passage.
