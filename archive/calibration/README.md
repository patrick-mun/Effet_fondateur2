# Archive de calibration R00–R14b

Cette archive conserve le dispositif de calibration historique du pipeline V2.
Elle couvre la conception R00–R09, les panels et pilotes R10–R12c, les
validations R13, puis les clôtures R14 et R14b. Elle a été isolée du paquet
actif après la fusion de la PR #20 (`9dfc0bc`).

## Décision finale

R14b est clos avec le statut
`COMPLETED_WITH_TOOL_INTEGRATION_VERIFIED`. L'intégration de SHAPEIT5 5.1.1
est vérifiée ; Hap-IBD et Refined IBD sont soutenus par les contrôles locaux.
La performance scientifique exhaustive prévue par R07/R08 reste
`NOT_EVALUATED`. Aucun effet fondateur n'est démontré et aucune promotion en
production n'est autorisée : `calibrated=false`,
`founder_effect_proven=false` et `production_promotion=false`.

## Contenu suivi par Git

- `code/` : ancien paquet `effet_fondateur.calibration`, conservé sans être
  installé par le projet actif ;
- `configurations/` : contrats YAML et schémas JSON historiques ;
- `documentation/` : conventions et plan R00–R14b ;
- `receipts/` : reçus historiques, déplacés sans réécriture de leurs chemins
  ou empreintes ;
- `tests/` : tests historiques retirés de la suite active ;
- `manifests/` : provenance de la relocalisation.

Les anciens chemins inscrits dans les reçus restent intentionnellement
inchangés. Le manifeste de relocalisation permet de les résoudre vers leur
emplacement actuel sans modifier les preuves d'origine.

## Preuves locales ignorées

Sur la machine où l'archivage a été réalisé, les releases, vérités
synthétiques, BCF/CSI, journaux, mesures de ressources, toolchains et rapports
HTML/PDF sont sous :

`archive/calibration/data-calibration/calibration_workspace/`

Les caches Python relocalisés sont sous `data-calibration/local-caches/` et
les anciennes arborescences temporaires R12c, vides, sous
`data-calibration/tmp-backups/`.
Tout `data-calibration/` est ignoré explicitement par Git. Ces fichiers ne
sont donc pas présents dans un clone Git et ne doivent jamais être publiés.
Le rapport R14b accepté se trouve localement dans
`calibration_workspace/r14b_integration_closure/` à l'intérieur de ce dossier.

## Éléments restés actifs

Le pipeline V2 reste sous `src/effet_fondateur/`. Le schéma
`schemas/explicit_ibd_calibration.schema.json` et la table
`explicit_ibd_calibration.tsv` produite par l'étape 16C restent actifs : dans
ce contexte, « calibration » désigne le réglage interne de l'appel IBD et non
la campagne R00–R14b.

## Limites et interdiction d'usage

Cette archive sert uniquement à l'audit historique. Son code dépend d'anciens
chemins et contrats et ne doit pas être ajouté au `PYTHONPATH`, exécuté comme
pipeline, ni utilisé pour produire une conclusion scientifique ou clinique.
Les résultats locaux sont synthétiques ou techniques ; ils ne remplacent pas
une qualification de performance R07/R08.

Une nouvelle calibration exige au minimum un objectif approuvé, un plan et
des critères prospectifs, des données et partitions autorisées, des ressources
confirmées, des versions d'outils figées, des tests indépendants, un nouvel
espace de noms de release et une décision explicite. Elle doit être développée
hors de cette archive et ne peut modifier rétroactivement ses reçus.
