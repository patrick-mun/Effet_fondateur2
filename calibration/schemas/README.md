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

Les références croisées entre identifiants et la résolution physique des
chemins seront contrôlées par les gardes D08. Aucun artefact généré ne doit être
stocké ici.
