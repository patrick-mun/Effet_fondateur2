# Étape 10 — QC final par cohorte

## Responsabilité

`10_qc_final` applique des politiques distinctes aux cohortes figées par
`09_freeze_cohorts`. Elle ne redéfinit aucune appartenance et vérifie que chaque
fichier PLINK `--keep` correspond exactement à `cohorts.frozen.tsv`.

Le calcul suit toujours le même ordre : extraction de la cohorte, missingness
individuelle, exclusion explicite des individus, recalcul des métriques variant
sur les individus restants, exclusion explicite des variants, puis publication
d'un nouveau jeu BED/BIM/FAM versionné.

## Politiques analytiques

- `controls_unrelated_qc`, issu du jeu genome-wide pré-QC : missingness
  individuelle et variant, MAF minimale et HWE chez les seuls témoins
  indépendants ;
- `target_carriers_independent_qc`, issu du jeu genome-wide pré-QC :
  missingness individuelle et variant uniquement ; ni MAF ni HWE ne peuvent
  exclure un signal chez les porteurs ;
- `target_chromosome_all_qc`, issu du jeu du chromosome cible avec génotypes
  explicites : missingness individuelle et variant, puis MAF locale.

Les différences de missingness entre lots produisent uniquement une alerte
descriptive. Elles ne retirent aucun individu ou variant automatiquement.

## Variant cible

Le variant déclaré dans `target_variant.yaml` doit être présent avant et après
le QC local. S'il dépasse un seuil standard, il reste dans le jeu final avec
`QC_STATUS=RETAINED_EXCEPTION`, les filtres concernés dans `FILTER_CODES` et
`FILTER_EXCEPTION=true`.

Cette conservation assure la traçabilité analytique ; elle ne valide pas un
génotype de mauvaise qualité. L'exception est comptée dans le rapport et dans
les avertissements de l'audit.

## Paramètres

Les seuils configurables sont des proportions ou probabilités comprises entre
`0` et `1` :

- `sample_missingness_max` ;
- `controls_variant_missingness_max`, `controls_maf_min` et
  `controls_hwe_p_min` ;
- `carriers_variant_missingness_max` ;
- `target_region_variant_missingness_max` et `target_region_maf_min` ;
- `batch_missingness_delta_alert`.

`plink_timeout_seconds` borne chaque appel externe. Modifier un seuil exige un
nouveau run, car la configuration résolue et la signature d'étape sont
immuables.

## Sorties

- trois jeux `*.bed`, `*.bim`, `*.fam` avec un descripteur `*.dataset.json` et
  de nouveaux `sample_set_id` et `variant_set_id` déterministes ;
- `sample_qc_final.tsv` : métriques, statut et motif individuel par jeu ;
- `variant_qc_final.tsv` : métriques et décision de chaque variant par jeu ;
- `batch_qc_final.tsv` : comparaison descriptive des lots ;
- `qc_final_report.json`, `stage_outputs.json`, `audit.json` et
  `checksums.sha256`.

## Codes de retour

- `0` : trois jeux finaux non vides publiés et variant cible conservé ;
- `2` : configuration, contrat, empreinte ou concordance de cohorte invalide ;
- `3` : PLINK absent, en échec ou rapport natif incohérent ;
- `4` : cohorte ou jeu de variants entièrement exclu, ou variant cible absent.
