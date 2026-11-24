# Calibration synthétique unifiée

Les producteurs autorisés de `panel_master_generic_v1` sont définis par les
contrats de `calibration/configurations/contract_bundle_v1.yaml`.

Le noyau contractuel R09 comprend :

- `contracts.py` : schémas, index fermé des étapes, vérités et vues ;
- `release_guard.py` : chemins, inventaires aveugles, disjonctions,
  pseudonymes et paires à facteur unique ;
- `truth_router.py` : accès minimal à la vérité après clôture ;
- `evaluation.py` : ordre de décision R08 ;
- `readiness.py` : reçu statique ouvrant seulement la construction R10.

R10 ajoute :

- `micro_contract.py` : topologie fermée des 25 étapes et de leurs sorties ;
- `micro_panel.py` : construction et validation de la release synthétique ;
- `micro_stage_worker.py` : travailleur technique aveugle et isolé.

R11 ajoute `pilot_panel.py`, le harnais générique qui matérialise les tailles
pilotes fermées, résout les exécutables épinglés, mesure chaque appel réel et
produit les 25 Markdown. Il ne lit aucune donnée d'étude réelle et ne modifie
pas le pipeline de production.

R12 ajoute `r12_freeze.py`. Ce harnais dérive sans écrasement la release R11,
ajoute les tags `AC/AN` aux deux VCF régionaux avec `bcftools +fill-tags`,
reprend uniquement 12, 12A et 16C, puis scelle configuration, critères,
compatibilités, modules et outils avant toute validation aveugle.

R12b ajoute `r12b_validation.py`, `r13_runner.py` et `r13_evaluator.py`. La
nouvelle release matérialise une partition `VALIDATION_BLIND` indépendante,
scelle sa vérité privée et fige les deux composants R13 avant exécution. Son
reçu G5 autorise uniquement une validation d'intégrité technique ; R07 reste à
réviser avant toute qualification scientifique.

La release R10 se construit et se valide sans donnée réelle ni outil
scientifique externe :

```bash
.venv/bin/python -m effet_fondateur.calibration.micro_panel build \
  --repository-root . \
  --release-root calibration_workspace/releases/micro_panel_generic_v1_r10

.venv/bin/python -m effet_fondateur.calibration.micro_panel validate \
  --repository-root . \
  --release-root calibration_workspace/releases/micro_panel_generic_v1_r10
```

Le constructeur refuse l'écrasement. Les vérités privées restent séparées des
onze vues aveugles ; les travailleurs ne peuvent pas les lire et l'évaluateur
n'ouvre que les familles autorisées après clôture. Les décisions R10 portent
uniquement sur l'intégrité technique, jamais sur une performance scientifique.

Le pilote R11 se construit dans une release neuve et non écrasable :

```bash
.venv/bin/python -m effet_fondateur.calibration.pilot_panel build
.venv/bin/python -m effet_fondateur.calibration.pilot_panel validate
```

Le run R11 du 12 novembre 2026 a validé PLINK, KING et bcftools sur 500 sujets
et les densités 600 000/300 000/100 000. SHAPEIT5 a correctement refusé les
dix graines parce que l'entrée synthétique R11 omettait le champ INFO `AC`
requis. Cette défaillance du harnais ne constituait ni une qualification ni un
défaut démontré du pipeline de production.

La correction et le gel R12 s'exécutent et se revérifient ainsi :

```bash
.venv/bin/python -m effet_fondateur.calibration.r12_freeze build
.venv/bin/python -m effet_fondateur.calibration.r12_freeze validate
```

La release `panel_master_generic_v1_r12_frozen_candidate` a terminé les dix
graines SHAPEIT5 ainsi que Hap-IBD et Refined IBD, sans commande en échec.
Les grilles pilotes restent incomplètes : les étapes 12, 12A et 16C sont donc
`NOT_EVALUABLE`, et les trois drapeaux de qualification restent faux. La
partition `VALIDATION_BLIND` n'a pas été ouverte.

Le candidat R12 initial est conservé mais rejeté. La release R12b se construit
et se revérifie sans lancer R13 :

```bash
.venv/bin/python -m effet_fondateur.calibration.r12b_validation \
  --release-root calibration_workspace/releases/panel_master_generic_v1_r12b_validation_blind

.venv/bin/python -m effet_fondateur.calibration.r12b_validation \
  --release-root calibration_workspace/releases/panel_master_generic_v1_r12b_validation_blind \
  --validate-only --print-g5-receipt
```

Après décision explicite de lancer R13, le runner aveugle puis l'évaluateur
post-clôture utilisent des commandes distinctes :

```bash
.venv/bin/python -m effet_fondateur.calibration.r13_runner \
  --candidate-root calibration_workspace/releases/panel_master_generic_v1_r12b_validation_blind \
  --workspace calibration_workspace/r13_validation_blind

.venv/bin/python -m effet_fondateur.calibration.r13_evaluator \
  --candidate-root calibration_workspace/releases/panel_master_generic_v1_r12b_validation_blind \
  --workspace calibration_workspace/r13_validation_blind
```

Ces commandes ne qualifient ni une performance scientifique ni un effet
fondateur. Elles maintiennent `calibrated=false`,
`founder_effect_proven=false` et `production_promotion=false`.

Le run officiel `RUN-R13-VALIDATION-BLIND-001` est clos : 25/25 étapes et 25/25
évaluations ont le statut technique `PASS`. Le reçu
`calibration/readiness/R13_TECHNICAL_VALIDATION.json` conserve néanmoins
`scientific_performance_evaluated=false`, R07 à `REVISION_REQUIRED` et la porte
G6 scientifique fermée.

R14 produit ensuite le rapport de clôture technique, sans recalcul scientifique :

```bash
.venv/bin/python -m effet_fondateur.calibration.r14_report \
  --repository-root . \
  --candidate-root calibration_workspace/releases/panel_master_generic_v1_r12b_validation_blind \
  --r13-workspace calibration_workspace/r13_validation_blind \
  --output-root calibration_workspace/r14_technical_closure_attempt02 \
  --template /Users/utilisateur/.codex/skills/report_skill/templates/recherche.html \
  --chrome "/Applications/Google Chrome.app/Contents/MacOS/Google Chrome"
```

La version acceptée comporte 13 pages A4, quatre figures et huit tables. Le reçu
`calibration/readiness/R14_TECHNICAL_CLOSURE.json` fixe ses empreintes et sa
décision : intégrité technique `PASS`, qualification scientifique non accordée,
effet fondateur non démontré et promotion non autorisée.

La révision post-R14 de R07 est décrite par
`calibration/configurations/scientific_dimensioning_v2.yaml`. Elle conserve les
minima scientifiques, remplace la multiplication grossière `22 × 50` par un
calendrier explicite de 640 appels SHAPEIT5 et ferme 90 mondes IBD. Son reçu se
vérifie sans lancer de calcul scientifique :

```bash
.venv/bin/python -m effet_fondateur.calibration.r07_revision \
  --repository-root . --validate-receipt
```

Le statut `REVISED_READY_FOR_R11B_RESOURCE_CONFIRMATION` autorise seulement les
sondes synthétiques R11b. Il ne vaut ni G5 scientifique, ni qualification, ni
accès à des données réelles.

R11b utilise un workspace non écrasable et des sous-commandes reprenables :

```bash
.venv/bin/python -m effet_fondateur.calibration.r11b_resources \
  --repository-root . prepare
.venv/bin/python -m effet_fondateur.calibration.r11b_resources \
  --repository-root . run-ibd
.venv/bin/python -m effet_fondateur.calibration.r11b_resources \
  --repository-root . run-shapeit --probe dlow
.venv/bin/python -m effet_fondateur.calibration.r11b_resources \
  --repository-root . status
```

Les sondes `bridge` puis `target` sont longues et doivent rester séquentielles.
La commande `project` n'est valide qu'après leurs reçus et les 18 mesures IBD.
Après leur clôture, le reçu suivi se reconstruit et se vérifie avec :

```bash
.venv/bin/python -m effet_fondateur.calibration.r11b_closure \
  --repository-root .
```

Ce reçu lie R07 v2, les trois mesures SHAPEIT5, les 18 mesures IBD et la
décision finale par SHA-256. Son statut autorise la construction de R12c, mais
maintient G5 scientifique et les trois drapeaux de qualification à `false`.

Le contrat R12c et sa release gelée se contrôlent sans lancer R13 :

```bash
.venv/bin/python -m effet_fondateur.calibration.r12c_scientific \
  --repository-root .
```

Le préflight est désormais `READY_FOR_G5_AUDIT`. Le reçu séparé
`R12C_G5_READY.json` ouvre R13 scientifique sans déclarer de résultat.

Les critères exécutables sont figés dans
`calibration/configurations/scientific_metric_rules_v1.yaml`. Le manifeste lie
deux fragments par SHA-256 et le validateur contrôle 135 règles sur les 25
étapes, les unités indépendantes et la statistique décisionnelle (`LCB95`,
`UCB95`, `CI95` ou valeur déterministe). Cette validation de contrat ne produit
aucun résultat scientifique et ne rend pas l'évaluateur prêt.

Les modules `r13_scientific_runner.py` et `r13_scientific_evaluator.py` sont
gelés. Le runner conserve les 1 180 appels, les exécute par lots de huit avec
checkpoints et refuse toute vérité ; l'évaluateur est post-clôture uniquement.

Le fractionnement fin est géré par un registre séparé, à raison d'un
checkpoint par appel externe. Son initialisation ne lance ni SHAPEIT5, ni
Hap-IBD, ni Refined IBD :

```bash
R13_RUN_ROOT=calibration_workspace/r13_scientific/RUN-R13-SCIENTIFIC-001
.venv/bin/python -m effet_fondateur.calibration.r13_checkpoint_ledger init \
  --run-root "$R13_RUN_ROOT"
.venv/bin/python -m effet_fondateur.calibration.r13_checkpoint_ledger status \
  --run-root "$R13_RUN_ROOT"
```

Le JSON est le registre autoritatif. Les vues
`governance/external_call_checkpoints.tsv` et
`governance/external_call_checkpoints.md` donnent pour chacun des 1 180 appels
un état et une colonne explicite `done`/`Fait` (`OUI` ou `NON`). Un lot peut
contenir de un à huit appels.

La révision `R12C-EXECUTABLE-REFREEZE-002` ajoute le matérialiseur scientifique
et le worker unitaire qui manquaient au premier gel. Le matérialiseur construit
les entrées avant réservation du checkpoint et garde la vérité dans une racine
privée distincte. Le worker reçoit uniquement la vue aveugle :

```bash
R13_RELEASE_ROOT=calibration_workspace/releases/panel_master_generic_v1_r12c_scientific_validation
R13_TRUTH_ROOT=calibration_workspace/private_truth/RUN-R13-SCIENTIFIC-001
.venv/bin/python -m effet_fondateur.calibration.r13_scientific_materializer \
  materialize --repository-root . --release-root "$R13_RELEASE_ROOT" \
  --run-root "$R13_RUN_ROOT" --truth-root "$R13_TRUTH_ROOT" \
  --call-id SHAPEIT-DMAX-CHR01-SEED01
.venv/bin/python -m effet_fondateur.calibration.r13_scientific_call_worker \
  preflight --repository-root . --release-root "$R13_RELEASE_ROOT" \
  --run-root "$R13_RUN_ROOT" --call-id SHAPEIT-DMAX-CHR01-SEED01
```

Le fond privé contient 1 200 000 rangs autosomaux. Les vues strictement
emboîtées sélectionnent respectivement un rang sur 2, 4 et 12, soit exactement
600 000, 300 000 et 100 000 variants. Les haplotypes suivent un modèle de
Markov à trois populations avec
`p_refresh = 1 - exp(-delta_cM / 0,05)` ; 200 trios sont transmis selon
`p_recomb = 1 - exp(-delta_cM / 100)`. Les cartes sont des cartes synthétiques
linéaires figées dans le cadre de coordonnées GRCh38 : elles ne doivent pas
être décrites comme des cartes de recombinaison humaines empiriques.

```bash
.venv/bin/python -m effet_fondateur.calibration.r13_checkpoint_ledger claim \
  --run-root "$R13_RUN_ROOT" --batch-size 1 --worker-id LOCAL-8CPU
.venv/bin/python -m effet_fondateur.calibration.r13_scientific_call_worker \
  execute --repository-root . --release-root "$R13_RELEASE_ROOT" \
  --run-root "$R13_RUN_ROOT" --call-id SHAPEIT-DMAX-CHR01-SEED01
```

Le worker `execute` clôt lui-même le checkpoint ; les commandes manuelles
`complete` et `fail` restent réservées à une reprise ou un diagnostic contrôlé.

Avant d'éteindre la machine, l'action `interrupt` transforme les appels encore
`RUNNING` en échecs reprenables sans toucher aux appels `COMPLETED`. Après un
transport ou avant la clôture, `status --verify-output-hashes` vérifie les
empreintes de toutes les sorties terminées. Un appel déclaré terminé dont la
sortie manque est automatiquement repris au prochain `claim`.

```bash
.venv/bin/python -m effet_fondateur.calibration.r13_checkpoint_ledger interrupt \
  --run-root "$R13_RUN_ROOT"
.venv/bin/python -m effet_fondateur.calibration.r13_checkpoint_ledger status \
  --run-root "$R13_RUN_ROOT" --verify-output-hashes
```

Le runner développe le calendrier en 640 appels SHAPEIT5 et 540 appels IBD.
L'évaluateur applique les bornes non arrondies, les contextes, les effectifs
indépendants et l'ordre des statuts R08 à un document post-clôture.

`scientific_statistics.py` et `scientific_metric_producers.py` fournissent le
socle commun : Wilson unilatéral, Clopper–Pearson exact, bootstrap déterministe
de 10 000 grappes, Holm, Monte-Carlo plus-un, contrôles exacts et provenance
des graines. Le préflight vérifie que ces méthodes couvrent les 135 règles,
et confirme `stage_producers_complete=true`.

`scientific_stage_producers.py` raccorde les 25 étapes `00` à `19`. Les
proportions utilisent leurs unités indépendantes, les erreurs continues sont
bootstrapées par grappe et aucune valeur déclarative `PASS` ne remplace les
preuves élémentaires.

Les autres modules de ce répertoire proviennent des bancs exploratoires de la
campagne abandonnée. Ils restent consultables pour audit, mais ne sont ni
exportés par le paquet ni autorisés à générer, évaluer ou publier un artefact
de la nouvelle release sans preuve d'adéquation et tests distincts.
