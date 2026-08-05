# Contrats de métadonnées V2 — maintenance et extension

## Objectif

Les tables `samples.master.tsv` et `cohorts.frozen.tsv` sont des frontières de
sécurité scientifique. Elles empêchent les scripts suivants d'interpréter
différemment les valeurs manquantes, les groupes ou les génotypes cibles.

La validation comporte trois niveaux :

1. structure TSV : encodage, en-tête, nombre et ordre des champs ;
2. schéma de ligne : type, vocabulaire contrôlé et champs conditionnels ;
3. règles métier : unicité, provenance, filiation, fichiers et références.

## Pourquoi les valeurs sont strictes

Les booléens acceptent uniquement `true` et `false`. Une cellule manquante est
vide. Les chaînes `NA`, `None`, `-` ou `0` ne sont pas converties implicitement,
car leur sens varie selon les outils et les colonnes.

Les statuts cliniques et les groupes restent descriptifs. Ils ne peuvent jamais
être utilisés comme source d'un génotype cible. Cette règle est contrôlée dans
le schéma et dans `contracts/samples.py`.

## Répartition des modules

| Module | Responsabilité |
|---|---|
| `contracts/tables.py` | parser et valider tout TSV versionné |
| `contracts/samples.py` | règles propres au registre maître |
| `contracts/cohorts.py` | règles propres aux cohortes figées |
| `samples_master.schema.json` | colonnes et types du registre |
| `cohorts_frozen.schema.json` | colonnes et types des cohortes |

Les erreurs indiquent la ligne et la contrainte, mais ne recopient pas la valeur
individuelle fautive dans les journaux ou la sortie CLI.

## Ajouter une colonne

Avant d'ajouter une colonne :

1. préciser sa définition scientifique et son vocabulaire ;
2. décider si elle est obligatoire ou nullable ;
3. définir une représentation unique de la valeur manquante ;
4. ajouter la colonne au tableau `x-tsv.columns` à sa position canonique ;
5. ajouter sa propriété au schéma de ligne ;
6. mettre à jour les producteurs et consommateurs ;
7. ajouter des tests de valeur valide, manquante et invalide ;
8. incrémenter la version si l'ancien fichier n'est plus accepté.

Une colonne ne doit jamais être ajoutée à la fin d'un fichier sans mise à jour du
schéma : le validateur refuse volontairement les champs supplémentaires.

## Ajouter un vocabulaire

Une nouvelle valeur de `ROLE`, `SEX` ou `CLINICAL_STATUS` constitue une décision
scientifique. Elle exige une justification et un test montrant son traitement
par les étapes consommatrices. Une valeur inconnue ne doit pas être convertie en
`OTHER` ou `UNKNOWN` sans règle explicite.

## Chemins sources

`SOURCE_FILE` est relatif à une racine fournie séparément. Le validateur refuse
les chemins absolus, les segments `..`, les séparateurs Windows et les liens
symboliques qui sortent de cette racine. Les fichiers sont seulement contrôlés
en lecture ; ils ne sont jamais modifiés.

## Extension à d'autres études

Les contrats ne contiennent ni nom de gène, ni chromosome, ni rsID imposé. Une
nouvelle étude fournit son profil de configuration et ses métadonnées, tout en
conservant les mêmes invariants de provenance et d'indépendance.

Les études X/Y, mitochondriales, structurales ou multi-variants nécessiteront un
nouveau contrat de génotype cible au lieu d'assouplir silencieusement le format
diploïde autosomique actuel.
