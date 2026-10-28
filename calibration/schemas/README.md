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

Le manifest de run sera ajouté à l'étape D07. Aucun artefact généré ne doit
être stocké ici.
