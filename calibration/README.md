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

## État actuel

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
ce diagnostic ne vaut toujours pas qualification. Le développement doit suivre
`PLAN_DEVELOPPEMENT_PANEL_SIMULATION.md` et ses portes de passage.
