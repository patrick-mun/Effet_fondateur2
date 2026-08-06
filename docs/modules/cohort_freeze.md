# Étape 09 — Gel des cohortes

## Responsabilité

`09_freeze_cohorts` transforme les décisions techniques et scientifiques des
étapes `02`, `04`, `05`, `07` et `08` en cohortes analytiques immuables. Elle
ne recalcule ni génotype, ni apparentement, ni structure populationnelle.

La table finale contient une ligne par individu et par cohorte. Une exclusion
reste donc visible avec son code et sa source au lieu de disparaître d'une
simple liste `--keep`.

## Cohortes publiées

- `controls_unrelated` : non-porteurs appartenant à un `GROUP_LABEL`
  explicitement déclaré comme témoin, retenus par l'ensemble indépendant et
  sans exclusion populationnelle approuvée ;
- `target_carriers_all` : tous les porteurs explicites inclus sur le chromosome
  cible, sans appliquer les exclusions KING/PCA ;
- `target_carriers_independent` : porteurs retenus par l'ensemble indépendant,
  sans exclusion populationnelle, avec un représentant technique par famille ;
- `family_noncarriers` : non-porteurs d'une famille contenant au moins un
  porteur explicite ;
- `affected_noncarriers` : individus `AFFECTED` dont le génotype explicite est
  non-porteur ;
- `target_chromosome_all_qc` : individus autorisés sur le jeu du chromosome
  cible et possédant un génotype accepté ;
- `genomewide_structure_reference` : référence exacte de l'étape `08`, moins
  les exclusions approuvées.

Le statut clinique n'est utilisé que conjointement à un génotype explicite pour
`affected_noncarriers`. Il ne produit jamais un génotype. Les groupes témoins
sont fournis par `control_group_labels` et restent une règle de rôle explicite.

## Unités indépendantes

Les contrôles et porteurs indépendants sont limités à un représentant par
`FID`. Le choix déterministe privilégie successivement le QC `PASS`, le plus
faible `F_MISS`, la plus faible valeur absolue du z-score d'hétérozygotie puis
`SAMPLE_ID`.

`INDEPENDENT_UNIT_ID` identifie l'unité de contrôle, la famille porteuse ou la
composante KING. `FAMILY_REPRESENTATIVE=true` n'est permis que sur une ligne
incluse possédant cette unité.

## Approbations

Lorsqu'une exclusion est proposée par `07` ou `08`, la configuration de
l'étape `09` doit fournir sous `decision_reviews` :

```yaml
kinship_exclusion_approval:
  review_id: kinship_review_001
  reviewer_role: scientific_reviewer
  reviewed_at: "2026-08-05T11:00:00Z"
  reviewed_artifact_sha256: <sha256 de independent_set_proposal.tsv>
  approved_exclusion_sample_ids: [sample_001]
```

La liste doit correspondre exactement aux exclusions proposées et l'empreinte
doit correspondre à l'artefact consommé dans le nouveau run. Une liste
incomplète, une empreinte différente ou une date sans fuseau bloque l'étape.
Changer un représentant ne relève pas de cette approbation et exige un nouveau
calcul amont.

Après succès et publication atomique, l'orchestrateur retire
`kinship_exclusion_approval` et/ou
`population_structure_exclusion_approval` du manifest. Une tentative échouée ne
résout aucune décision. Les identifiants approuvés sont conservés uniquement
dans `cohort_decisions.tsv`, classé `sensitive_genetic`; l'audit ne publie que
des comptes.

## Sorties

- `cohorts.frozen.tsv` : matrice complète des inclusions et exclusions ;
- `cohort_summary.tsv` : effectifs et unités agrégés pour les sept cohortes ;
- `cohort_decisions.tsv` : exclusions humaines appliquées et leur provenance ;
- `keep/<cohort>.keep` : couples PLINK `FID IID` inclus ;
- `cohort_freeze_report.json` : résumé sans identifiant individuel ;
- `stage_outputs.json`, `audit.json` et `checksums.sha256`.

`required_nonempty_cohorts` fixe les cohortes dont l'absence bloque le run. Les
cohortes biologiquement optionnelles, notamment `affected_noncarriers`, restent
publiées avec un effectif nul dans le résumé.

## Codes de retour

- `0` : cohortes gelées, revues validées et décisions résolues ;
- `2` : configuration, contrat ou cohérence inter-étapes invalide ;
- `4` : revue absente/incomplète, empreinte différente ou cohorte requise vide.
