# Étape 07 — Apparentement et concordance du pedigree

## Responsabilité

`07_infer_kinship` exécute KING sur le panel indépendant de l'étape `06`,
classe les relations observées, les compare aux liens parent-enfant déclarés et
publie une proposition révisable d'ensemble d'individus indépendants.

L'étape ne modifie aucun génotype et n'applique aucune exclusion. Toute
proposition d'exclusion ajoute la décision `kinship_exclusion_approval` au
manifest du run.

## Entrées

- `kinship_panel.bed/.bim/.fam` et son descripteur ;
- `samples.master.tsv` pour les identifiants et filiations déclarées ;
- `qc_individual_metrics.tsv` pour départager les représentants techniques.

Les individus du panel doivent être un sous-ensemble exact du registre et des
métriques QC. Un individu déjà exclu à l'étape `05` ne peut pas réapparaître.

## Paramètres et classification

La configuration versionne les seuils inférieurs suivants :

| Classe | Seuil de kinship par défaut |
|---|---:|
| duplicat ou jumeau monozygote | `0.354` |
| premier degré | `0.177` |
| deuxième degré | `0.0884` |
| troisième degré | `0.0442` |
| quatrième degré | `0.0221` |

Les seuils doivent être strictement décroissants. Pour une paire du premier
degré, `parent_child_ibs0_max` distingue une relation compatible parent-enfant
d'une autre relation du premier degré. Ce contrôle ne remplace pas une revue du
pedigree.

`king_degree` fixe la profondeur demandée à KING et
`relatedness_max_degree_for_independence` fixe les relations transformées en
arêtes du graphe de proposition. Le second degré maximal ne peut pas dépasser
le premier.

## Sorties

- `kinship_pairs.tsv` : paires publiées par KING, fichier génétique sensible ;
- `kinship_degree_summary.tsv` : comptes agrégés par degré ;
- `pedigree_concordance.tsv` : chaque parent déclaré, y compris s'il est absent ;
- `independent_set_proposal.tsv` : proposition par individu avec représentant
  et provenance de l'arête ;
- `kinship_report.json` ;
- `stage_outputs.json`, `audit.json` et `checksums.sha256`.

Les fichiers KING natifs `.kin` et `.kin0` sont validés puis supprimés du
dossier publié. Les résultats normalisés en conservent les métriques utiles et
leur provenance dans un contrat stable.

## Concordance

- `CONCORDANT` : premier degré et IBS0 compatible avec le parent-enfant déclaré ;
- `DISCORDANT_DEGREE` : degré observé incompatible ;
- `DISCORDANT_RELATION` : premier degré mais IBS0 non compatible ;
- `CRYPTIC_RELATEDNESS` : relation observée entre familles sans lien déclaré ;
- `NOT_EVALUATED_PARENT_ABSENT` : parent déclaré non génotypé ;
- `NOT_EVALUATED_KING_PAIR_MISSING` : relation déclarée sans ligne KING.

KING ne matérialise entre familles que les paires atteignant la profondeur
demandée. L'absence d'une paire inter-familles n'est donc pas publiée comme une
mesure individuelle de non-apparentement.

## Proposition indépendante

Le graphe contient les relations KING jusqu'au degré configuré et tous les
liens parent-enfant déclarés entre individus génotypés. Dans chaque composante,
un algorithme glouton déterministe privilégie successivement :

1. le statut QC `PASS` ;
2. le plus faible `F_MISS` ;
3. la plus faible valeur absolue du z-score d'hétérozygotie ;
4. l'identifiant stable pour départager les égalités.

Le résultat est un ensemble indépendant maximal orienté qualité. Il n'est pas
présenté comme l'ensemble de cardinalité maximum. `RELATION_BASIS` indique si
une exclusion proposée repose sur KING, le pedigree déclaré ou les deux.

## Codes de retour

- `0` : calcul et contrats valides, même si une revue manuelle est requise ;
- `2` : configuration, seuil ou entrée invalide ;
- `3` : échec de KING ou sortie native incohérente.
