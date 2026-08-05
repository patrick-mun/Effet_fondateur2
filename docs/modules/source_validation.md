# Validation des sources ACPA V2

## Responsabilité

L'étape `01_validate_sources` inventorie les exports ACPA et contrôle leur
structure sans produire de génotype dérivé. Elle ne modifie jamais les fichiers
sources et ne conserve aucune valeur de génotype dans son audit.

## Entrées déclarées

L'orchestrateur résout `inputs.acpa_samples_dir`, refuse un répertoire absent ou
vide, puis déclare chaque fichier dans `stage_inputs.json` avec son empreinte
SHA-256. La signature de l'étape inclut ces descripteurs. Une source modifiée
après une publication bloque donc la reprise du run.

Les liens symboliques ne sont pas des sources acceptées. Les chemins relatifs
restent relatifs au dépôt, conformément aux profils de configuration actuels.

## Contrôles

Pour chaque export, l'étape vérifie :

- l'extension autorisée et l'encodage UTF-8 ;
- les métadonnées `Array Type Name` et `UCSC Genomic Version` ;
- les colonnes ChAS obligatoires ;
- la présence de lignes de données ;
- la couverture autosomique configurée et le chromosome cible ;
- l'empreinte calculée juste avant lecture.

Le jeu complet vérifie aussi l'effectif attendu lorsqu'il est configuré,
l'unicité des noms de fichiers et la cohérence de l'assemblage annoté.

## Sorties

- `source_inventory.tsv` : taille, empreinte, métadonnées et résumé de couverture ;
- `source_qc.tsv` : un code stable par contrôle, sans recopier de génotype ;
- `audit.json` : effectifs agrégés et provenance standard V2.

Une anomalie de source retourne le code `4`. Les sorties restent alors dans le
dossier de tentative échouée pour diagnostic et ne deviennent jamais des
artefacts scientifiques publiés.
