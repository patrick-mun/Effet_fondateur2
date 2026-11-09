# Sous-système de calibration

Ce dossier contient le code de configuration, d'isolation, de simulation, de
vérité, d'évaluation et de production des rapports de calibration.

`isolation.py`, introduit à l'étape D08, contrôle les chemins physiques, liens
symboliques, références de manifeste et collisions entre panels avant toute
lecture d'un artefact. `partitioning.py`, introduit à l'étape D10, répartit de
façon déterministe les groupes indépendants entre donneurs de calibration,
donneurs de validation et référence de phasage commune. Il contrôle la
couverture de la source, les intersections, les empreintes et la confidentialité
du manifest sans écrire d'identifiant.

`materialize_partitions.py` vérifie les empreintes des métadonnées D09, refuse
l'écrasement d'un run et écrit atomiquement les trois listes privées ainsi que
leur manifest agrégé dans `calibration_workspace/`. Sa sortie standard ne
contient que des effectifs et statuts de contrôle.

`simulation/base_haplotypes.py`, introduit en D11, réapparie les copies
haplotypiques de donneurs distincts au sein de chaque strate publique. La
bijection conserve exactement fréquences alléliques et LD haplotypique ; les
génotypes non phasés sont ensuite dérivés des deux copies sans approximation.
`materialize_base_haplotypes.py` vérifie les MD5 fournisseur, les SHA-256
locaux et les liaisons aux catalogues D09, extrait les SNV bialléliques
complètement phasés avec `bcftools`, interpole la carte sans extrapolation et
écrit atomiquement les quatre artefacts privés D11.

`simulation/chip_markers.py` et `materialize_chip_observation.py` constituent
la couche D11b. Ils vérifient le paquet officiel CytoScan 750K Accel
`NA36.r1` et le manifest dense D11a. `simulation/chip_package.py` ouvre en
lecture seule sa base d'annotations SQLite, valide son identité et produit la
table SNP autosomique normalisée dans le workspace privé. D11b construit ensuite
l'intersection chromosome–position–allèles, un VCF observable et l'audit des
exclusions. Les dossiers connus d'exports d'échantillons sont refusés avant
lecture ; la vérité dense reste inchangée et le cache du pipeline n'est pas lu.

`simulation/population_structure.py` et
`materialize_population_structure.py` implémentent D12. Ils exigent des
manifests D11a/D11b finalisés, liés au même panel et au même domaine. Des paires
de diploïdes D11 provenant de deux strates donnent deux mélanges réciproques
`A1|B2` et `B1|A2`, ce qui conserve chaque haplotype exactement une fois. Sans
recombinaison, seules les proportions pures ou 50/50 sont exactes ; les
proportions inégales sont différées à D14. La vérité et la provenance restent
privées, tandis que le VCF observable ne contient aucune étiquette d'ascendance.

`simulation/pedigrees.py` et `materialize_pedigrees.py` implémentent D13. Ils
sélectionnent sans remplacement quatre fondateurs D12 par unité familiale,
construisent quatre trios jusqu'aux cousins germains et calculent la parenté
théorique depuis le pedigree. Chaque enfant reçoit une copie entière de chaque
parent ; l'ascendance est dérivée des copies D12 effectivement transmises et
contrôlée contre la provenance haplotypique. Les vérités restent privées et le
VCF observable ne contient aucun libellé familial ou d'ascendance.

`simulation/founder_recombination.py` et
`materialize_founder_recombination.py` implémentent D14. Ils vérifient la chaîne
empreintée D11a–D13, utilisent la carte génétique réelle du run D11a et tirent
indépendamment le premier breakpoint de chaque bras. La mosaïque fondatrice est
propagée par origine haplotypique afin de préserver Mendel. Les limites hors
domaine sont censurées, la mutation synthétique est injectée sur un SNV
existant et les vérités ne sont jamais écrites dans le VCF observable.

`simulation/windows.py` et `materialize_simulation_windows.py` implémentent
D14b. Ils calculent le quantile de censure par âge, inversent la carte en
arrondissant les bornes vers l'extérieur, refusent toute extrapolation et
harmonisent l'univers CytoScan avec le VCF public avant d'accepter les seuils
IBD. `materialize_window_replay_receipt.py` vérifie ensuite les liaisons
SHA-256 D14b→D11a→D11b→D12→D13→D14 et écrit un reçu privé atomique.

`simulation/ibd.py` et `materialize_ibd.py` implémentent D15 : dix paires
positives/négatives disjointes sur les cinq longueurs D14b, sélection hors
pedigree, copie H1 commune, invariance des négatifs, contrôle Mendel, empreintes
et schémas avant publication. Les matrices sont distinctes des sujets ; le
fond public peut conserver du partage inconnu. `report_ibd_injection.py`
produit un rapport technique français indépendant à partir d'un modèle
report-formater fourni explicitement, sans estimer de sensibilité IBD.

`simulation/roh.py` et `materialize_roh.py` implémentent D16. Ils copient H1
sur H2 pour réaliser trois origines autozygotes, dont une sur un sujet IBD
positif. H1 et les copies IBD sont contrôlées, les sujets de transmission et
les négatifs/matrices IBD sont protégés. Les seuils observables sont vérifiés
avant écriture. Les huit artefacts conservent séparément ROH, IBD original,
IBD induit H2–H1, pedigree, transmissions, audit de densité et VCF.
`report_roh_injection.py` publie le
rapport de développement à partir d'un modèle report-formater explicite.

`simulation/measurement_errors.py` et `materialize_measurement_errors.py`
implémentent D17 : vérité phasée projetée sur sondes, profils appariés,
retrait de marqueurs, dosage erroné, points et blocs de phase, puis manque.
Un registre chronologique privé permet un rejeu sans rappel de l'injecteur ;
les deux VCF sont relus et comparés intégralement avant publication atomique.
La densité peut être réduite sans changer les seuils scientifiques ; les
sources biologiques D16 sont vérifiées avant/après. Les candidats phasés ne
sont jamais présentés comme des sorties de phasage du pipeline.
`report_measurement_errors.py` produit le rapport autonome CAL-DEV-D17 depuis
les seuls agrégats du manifest et un modèle report-formater fourni en argument.

`export_blinded_inputs.py` implémente D18 à partir des seuls VCF observés D17.
Les colonnes sont pseudonymisées/permutées ; les champs autorisés sont
reconstruits dans des VCF et ACPA, un registre neutre et un profil pipeline
00–02. La correspondance, les profils et graines ne quittent pas le contrôle
privé. Relecture exhaustive, vrais lecteurs 01/02, inventaire fermé, schémas
et empreintes précèdent la publication atomique. Les entrées moléculaires et
familiales non disponibles restent explicitement manquantes.
`report_blinded_export.py` rend CAL-DEV-D18 sans exposer de sujet individuel.

`technical_evaluation.py` fournit les comparateurs purs et le lecteur BED
indépendant D19. `evaluate_technical.py` lie les sorties 00–04 aux références
D18, vérifie les empreintes/contrats, les six champs FAM, les axes et appels,
et distingue FAIL de NOT_EVALUABLE. Une attente cible 04 doit être explicite,
privée, liée au lot et à la définition configurée. `run_technical_inputs.py`
est un runner diagnostique séparé limité à 00–02, sans outils génétiques ni
approbation automatique ; l'évaluateur lui-même reste en lecture seule.
`report_technical_evaluation.py` produit CAL-DEV-D19 et cinq CAL-REP autonomes.
Le rejeu régional ne valide pas les étapes 03–04 non exécutées et ne change
aucun calcul scientifique.

`technical_scenarios.py` prépare des copies fermées pour 00–04 avec annotation
publique, génotypes synthétiques explicites et approbation technique fournie.
`run_technical_scenarios.py` exécute seulement 00–04 avec PLINK, distinctement
du runner historique 00–02. Le positif conserve D18 intact ; le négatif ajoute
une cible à une copie pour vérifier le refus de doublon.
`technical_rejection.py` exige code métier, motif exact, identité de tentative,
absence de publication et base inchangée. `evaluate_technical` accepte
`--technical-panel-root` uniquement si le reçu est lié par empreinte ; les
anciens formats D19 restent lisibles.

`qc_bench.py` prépare D20 : quatorze stress IID appariés autosomiques,
registres neutres et cohortes de fixture publiques, distincts des vérités
privées. `qc_oracle.py` calcule indépendamment manques, MAF, hétérozygotie,
F avec correction small-sample, HWE exact et décisions séquentielles.
`run_qc_bench.py` borne l'exécution à 00–05 puis au composant 10 isolé ; aucun
faux succès 09/10 n'est ajouté au manifest d'orchestration.
`evaluate_qc.py` contrôle clôture, contrats, motifs, métriques, matrices de
confusion et triplets retenus. `report_qc.py` publie les trois rapports
autonomes sans recalcul statistique. Les critères du pilote restent non
scellés : les PASS de conformité ne qualifient pas la performance biologique.

`simulation/kinship_bench.py` génère D21 par gene-drop mendélien IID, avec
coefficients de pedigree exacts et origines des copies. `kinship_bench.py`
sépare les observations neutres et la vérité privée ; `run_kinship_bench.py`
borne le vrai pipeline à 00–03/05–07, sans approbation des exclusions.
`kinship_oracle.py`, `kinship_panel_oracle.py` et `kinship_stage_oracle.py`
contrôlent indépendamment les génotypes, MAF/LD, équations KING, degrés,
proposition et vérité. `evaluate_kinship.py` exige clôture, empreintes et
rejeu exact. `report_kinship.py` publie les trois rapports indépendants,
sans calcul statistique. La conformité du code et l'exactitude de détection
restent deux verdicts distincts ; les profils appariés ne sont pas cumulés.

`structure_bench.py` prépare D22 et ses fixtures amont explicites, distinctes
de D21. `run_structure_bench.py` appelle 00–03 puis 08/09 isolées, avec revues
de test déclarées et gels répétés. `structure_oracle.py` reconstruit la PCA
par Gram/eigh et la queue chi-deux indépendante à deux axes ;
`cohort_oracle.py` reconstitue les sept cohortes sans helper de 09.
`evaluate_structure.py` exige clôture, rejeu, sources, signatures et entrées
cohérentes ; il vérifie les tableaux, invariances, keep et répétitions.
`report_structure.py` publie trois rapports agrégés, sans recalcul statistique.
Le FAIL 07 est conservé pour la revue globale, pas résolu par ce gel isolé.

Le rendu D22 v2 utilise `structure_report_evidence.py` pour extraire des
descriptions agrégées des sorties closes, liées à l'évaluation par empreintes.
`structure_report_graphics.py` les met en SVG ; `structure_report_math.py`
explique les équations réellement utilisées. Le moteur de rendu ne recalcule
aucune PCA ni décision. Les arguments `--panel-root` et `--pipeline-root`
désignent les sources closes des nouveaux graphiques, en plus du `--run-root`
d'évaluation. Les données visuelles ne contiennent aucun identifiant individuel.

`phasing_preparation.py` et `materialize_phasing_panel.py` préparent D23a
sur les 22 autosomes, sans lancer SHAPEIT5 : empreintes publiques, liaison
puce, copies fondatrices D10, transmissions Haldane, génotypes moléculaires
explicites, référence disjointe et relecture complète. La clôture v2 change
seulement les paramètres de phasage/chemins, jamais les génotypes du brouillon
v1 ; les graines de simulation restent privées. `report_phasing_preparation.py`
produit CAL-DEV-D23A avec deux graphiques et deux rappels mathématiques,
sans recalcul ni prétention de qualification du phasage.

`phasing_bench.py` exécute D23 : conversion VCF–PLINK et composants réels
10A/11/12/12A, avec fixtures QC/cohorte déclarées. `run_phasing_regional.py`
assure une reprise indépendante des anciens contrats mal enveloppés sans
changer les génotypes. `phasing_oracle.py` aligne la phase modulo H1/H2,
compte transitions, motifs isolés, erreurs rares et orientation hors cible ;
le bootstrap conserve les familles. `evaluate_phasing.py` contrôle clôture,
inventaires, signatures, dosages, cartes, univers et critères D04, puis ouvre
seul la vérité. `phasing_report_graphics.py` et `report_phasing.py` publient
quatre rapports indépendants depuis les seuls agrégats fermés. Un taux
mesurable reste distinct d'une performance qualifiée ; aucune calibration
scientifique ou validation externe n'est attribuée automatiquement.
