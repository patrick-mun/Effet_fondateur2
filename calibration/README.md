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
