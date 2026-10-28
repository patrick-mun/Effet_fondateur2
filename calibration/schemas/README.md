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
un SHA-256 local enregistré dans le manifest D07.

`reference_partition_manifest.schema.json`, créé à l'étape D10, décrit trois
partitions privées : donneurs de calibration, donneurs de validation et
référence de phasage commune. Il interdit les identifiants bruts dans le
manifest de contrôle, exige les listes sous `private_truth/partitions/`, leurs
SHA-256 et les comptes agrégés, puis impose une intersection nulle des
identifiants et des groupes de pedigree ou de parenté.

Les références croisées entre identifiants et la résolution physique des
chemins seront contrôlées par les gardes D08. Aucun artefact généré ne doit être
stocké ici.
