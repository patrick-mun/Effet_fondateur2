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

## Tables de résultats sans observation

Une table TSV vide est refusée par défaut. Un schéma peut déclarer
`x-tsv.allow_empty: true` uniquement lorsqu'une absence de ligne constitue un
résultat scientifique valide, par exemple lorsqu'aucune paire apparentée n'est
détectée. L'en-tête canonique reste obligatoire et le producteur doit consigner
le nombre de lignes dans son audit ; cette option ne doit jamais masquer une
entrée absente ou une étape incomplète.

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

## Étape 02 et approbation humaine

`02_build_sample_registry` reçoit la table utilisateur canonique par
`inputs.sample_metadata` et les artefacts `source_inventory` et `source_qc` de
l'étape `01`. Elle refuse une source absente de l'inventaire ou affectée à
plusieurs individus, puis publie une copie immuable `samples.master.tsv` et
`sample_registry_review.tsv`.

Sans approbation, l'étape peut réussir techniquement mais ajoute
`sample_registry_approval` aux décisions manuelles du run. L'étape `03` reste
alors bloquée. Une approbation doit être déclarée dans la configuration avant la
création d'un nouveau run :

```yaml
stages:
  build_sample_registry:
    enabled: true
    parameters:
      manual_approval:
        approved: true
        decision_id: registry_review_001
        reviewer_role: data_steward
        approved_at: "2026-08-05T10:00:00Z"
```

Le contrat enregistre un identifiant de décision et un rôle, pas le nom d'une
personne. Modifier l'approbation ou la table source exige un nouveau run ; la
configuration résolue d'un run existant reste immuable.
