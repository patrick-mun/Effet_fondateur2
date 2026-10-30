# Schémas de calibration

Ce dossier contient les schémas versionnés des configurations, vérités
simulées, manifests de run, évaluations d'étapes et décisions de qualification.

`calibration_config.schema.json`, créé à l'étape D05, impose les graines, les
effectifs propres à chaque phase, les scénarios minimaux et les règles de
sécurité.

`truth_manifest.schema.json`, créé à l'étape D06, décrit les neuf familles de
vérité nécessaires aux évaluateurs : registre synthétique, génotypes et phase,
ascendance, parenté, variation cible, transmission fondatrice, segments IBD,
segments ROH et erreurs injectées. Chaque artefact doit rester sous
`private_truth/`, posséder une empreinte SHA-256, déclarer son contrat de
champs, ses unités et les rapports autorisés à le consommer. Les scénarios IBD
et ROH négatifs sont représentables avec des bornes nulles explicites.

`calibration_run_manifest.schema.json`, créé à l'étape D07, impose la
provenance nécessaire au rejeu d'un run : configuration et vérité empreintées,
registre privé des six graines, commit propre, environnement Python verrouillé,
sources et outils versionnés, commandes normalisées, limites de ressources,
artefacts, reprises, écarts et ordre de rejeu. Ce manifest est un document
privé de contrôle ; il n'est jamais remis au pipeline. Une commande de portée
`PIPELINE` ne peut pas référencer `private_truth/`.

`public_reference_catalog.schema.json`, créé à l'étape D09, impose une liste
fermée de sources publiques possédant chacune un identifiant, une version, une
origine HTTPS, un assemblage explicite, une décision d'usage local, une
empreinte SHA-256 et une liaison empreintée au catalogue technique du pipeline.
Il interdit les données d'étude réelles, le téléchargement pendant D09 et la
redistribution par le projet. Tout fichier acquis ultérieurement doit recevoir
un SHA-256 local enregistré dans le manifest D07. La cinquième source
`ARRAY_DESIGN_MANIFEST` Thermo Fisher a été ajoutée après acquisition,
empreinte et revue des conditions d'usage.

`reference_partition_manifest.schema.json`, créé à l'étape D10, décrit trois
partitions privées : donneurs de calibration, donneurs de validation et
référence de phasage commune. Il interdit les identifiants bruts dans le
manifest de contrôle, exige les listes sous `private_truth/partitions/`, leurs
SHA-256 et les comptes agrégés, puis impose une intersection nulle des
identifiants et des groupes de pedigree ou de parenté.

`base_haplotype_manifest.schema.json`, créé à l'étape D11, lie les haplotypes
phasés, les génotypes non phasés, la provenance donneur privée et la carte
génétique interpolée. Il impose une utilisation unique de chaque copie
haplotypique, interdit les paires issues du même donneur, exige une concordance
génotype–haplotype intégrale et refuse tout identifiant donneur ou chemin local
dans le manifest.

`chip_marker_harmonization_manifest.schema.json` complète D11 par la couche
D11b. Il lie le paquet constructeur CytoScan 750K Accel `NA36.r1`, sa table SNP
normalisée et la vérité dense D11a par SHA-256. Il impose GRCh38, la politique
d'intersection et d'exclusion des ambiguïtés, puis inventorie la table SNP
normalisée, la table de liaison, le VCF observable et l'audit agrégé. Toute
source dérivée d'un échantillon réel est déclarée interdite.

`population_structure_config.schema.json` et
`population_structure_manifest.schema.json` encadrent D12. Avant les mosaïques
de recombinaison D14, ils n'autorisent que les strates pures et les mélanges
diploïdes exacts à deux composantes 50/50. Le manifest lie par SHA-256 D11a,
D11b et la configuration, inventorie la vérité d'ascendance privée, la
provenance D12→D11, le panel phasé dense et sa copie observable sans étiquette.

`pedigree_config.schema.json` et `pedigree_manifest.schema.json` encadrent D13.
Ils figent le pedigree de cousins germains, la sélection sans remplacement dans
une même strate, la graine, l'interdiction des recombinaisons avant D14 et les
liaisons SHA-256 vers les quatre entrées D12. Le manifest inventorie séparément
les vérités de pedigree, parenté, transmission, provenance et ascendance, ainsi
que les panels dense et aveugle. Tous les compteurs d'incohérence mendélienne,
de copie fondatrice ou de chevauchement familial doivent être nuls.

`founder_simulation_config.schema.json` et
`founder_simulation_manifest.schema.json` encadrent D14. Ils figent l'âge, la
graine, le variant synthétique, le processus exponentiel des premiers
breakpoints et la censure aux bords. Le manifest lie D11a, D12 et D13 par
SHA-256, inventorie neuf artefacts privés et sépare segments, porteurs, vérité
du variant et mosaïque d'ascendance locale. Une recombinaison censurée ne peut
pas être comptée comme un événement observé.

`simulation_window_config.schema.json` et
`simulation_window_manifest.schema.json` encadrent D14b. Ils imposent le
dimensionnement primaire en `cM`, la séparation du défi de censure, une cible
publique présélectionnée, les seuils de marqueurs IBD et l'interdiction de
filtrer les tirages selon leur résultat. Le manifest conserve les empreintes
des trois références publiques, les bornes réalisées et l'audit de densité.

`simulation_window_replay.schema.json` scelle séparément le rejeu D11a–D14.
Il lie les six manifests par SHA-256, contrôle le domaine, la cible, le rôle de
panel, l'observation d'un chemin de recombinaison interne et l'absence d'erreur
mendélienne, sans modifier un manifest déjà finalisé.

Les références croisées entre identifiants et la résolution physique des
chemins seront contrôlées par les gardes D08. Aucun artefact généré ne doit être
stocké ici.

`ibd_simulation_manifest.schema.json` et `ibd_injection_truth.schema.json`
encadrent D15. Le premier lie D14b, le reçu, D11a et D14 par empreintes et
impose zéro erreur de copie, de modification négative ou de Mendel. Le second
décrit chaque réplication, paire, matrice, longueur demandée et bornes réalisées
en bp/cM ; les négatifs exigent des bornes et une matrice nulles.

`roh_simulation_config.schema.json`, `roh_injection_truth.schema.json`,
`roh_induced_ibd_truth.schema.json` et `roh_simulation_manifest.schema.json`
encadrent D16. Ils imposent quatre
scénarios, trois origines autozygotes et un témoin avec bornes nulles. Les
copies H1/H2 sont explicites et une ROH est couplée à un positif D15 nommé.
Le manifest lie configuration et sources, huit artefacts privés, la densité
observable et les compteurs d'erreur nuls. `gamma_inputs_modified` et
`calibrated` restent `false`.
La vérité IBD induite déclare explicitement le partage H2–H1 supplémentaire
créé sur l'intersection ROH/IBD de la paire couplée.

Les trois schémas `measurement_error_config`, `measurement_error_event` et
`measurement_error_manifest` encadrent D17 : profils stricts, graine de mesure,
registre chronologique des valeurs avant/après et manifest des 26 artefacts
pour huit observations appariées. Les frontières et GT inversés sont distincts.
Le validateur Python complète les schémas par les quotas, sous-graines,
dénominateurs, inventaires et le contrôle intact. Les vérités D16 restent
inchangées ; aucun candidat phasé n'est une entrée aveugle.

Les schémas `blinded_export_config`, `blinded_input_inventory`, `blinding_map`
et `blinded_export_manifest` encadrent D18. L'inventaire opérateur ne peut
référencer que les fichiers autorisés de lots opaques ; le manifest et la
correspondance privés conservent les sources, profils, graines et liaisons.
La validation Python impose aussi les bijections, dénominateurs, chemins par
rôle, artefacts exhaustifs et correspondances avec le registre exporté.

Les quatre schémas D19 `technical_run_index`, `technical_evaluation`,
`technical_evaluation_manifest` et `technical_target_expectation` encadrent
le reçu CLOSED, les cinq résultats par lot, leur scellement agrégé et l'attente
optionnelle de cible synthétique explicite. Aucune métrique sans mesure n'est
remplacée par zéro. La validation Python contrôle les identités des étapes,
les populations identiques, les empreintes et la cohérence verdict/motifs et
des compteurs ; elle refuse toute fuite de chemin ou source réelle.

`technical_scenario_config` et `technical_panel_manifest` ajoutent le dérivé
D19 autorisant 00–04. L'annotation passe aussi le contrat moléculaire existant ;
les génotypes restent synthétiques. Les fichiers, populations, copies D18 et
ajouts du négatif sont contrôlés en Python. Le reçu et le manifest d'évaluation
peuvent ancrer ce panel et préciser son scénario ; ces extensions restent
facultatives pour relire les anciens runs 00–02. Le contexte et la provenance
de la cible sont alors obligatoires ensemble.
