# Étape 06 — Panel indépendant de marqueurs

## Responsabilité

`06_build_kinship_panel` construit, depuis `genomewide_pre_qc`, un jeu PLINK
autosomique dédié à KING et à la structure populationnelle. Ce jeu reste
strictement distinct du chromosome cible et ne prend aucune décision
d'exclusion d'individus.

La MAF est recalculée après les exclusions techniques de l'étape `05`. Les
alertes MAF préliminaires sont conservées dans l'audit des variants, mais ne
sont pas utilisées comme fréquence courante.

## Paramètres

- `maf_min` : fréquence allélique minimale recalculée sur le jeu pré-QC ;
- `ld_window_variants` : taille de fenêtre PLINK `--indep-pairwise` ;
- `ld_step_variants` : pas de déplacement de cette fenêtre ;
- `ld_r2_max` : seuil maximal de corrélation du pruning ;
- `required_autosomes` : chromosomes dont la couverture est obligatoire ;
- `min_markers_per_autosome` : minimum bloquant par chromosome requis ;
- `min_total_markers` : minimum bloquant dans le panel complet ;
- `max_autosome_marker_fraction` : concentration maximale sur un chromosome ;
- `excluded_regions` : régions complexes explicitement nommées en coordonnées bp ;
- `residual_ld_window_variants` et `residual_ld_window_kb` : fenêtre du contrôle
  de LD résiduel ;
- `plink_timeout_seconds` : délai maximal de chaque commande PLINK.

Les valeurs de l'exemple sont des valeurs de configuration, pas des seuils
scientifiques validés pour toutes les cohortes. Leur modification impose un
nouveau run.

## Traitements

1. validation du descripteur et des empreintes du jeu pré-QC ;
2. vérification que tous les variants présents ont passé les exclusions
   techniques de `05` ;
3. recalcul PLINK des fréquences sur les individus conservés ;
4. retrait des faibles MAF et des régions complexes configurées ;
5. pruning PLINK `--indep-pairwise` et validation de la partition
   `prune.in`/`prune.out` ;
6. création d'un nouveau BED/BIM/FAM sans changement de l'ordre des individus ;
7. contrôle du nombre, de la répartition et des distances entre marqueurs ;
8. calcul d'un LD résiduel puis publication de sa distribution agrégée.

## Sorties

- `kinship_panel.bed/.bim/.fam` et `kinship_panel.dataset.json` ;
- `pruned_variants.tsv`, une ligne par variant du jeu pré-QC ;
- `panel_coverage.tsv`, une ligne pour chacun des 22 autosomes ;
- `panel_residual_ld_bins.tsv`, distribution agrégée par classes de R² ;
- `kinship_panel_report.json` ;
- `stage_outputs.json`, `audit.json` et `checksums.sha256`.

Les paires individuelles utilisées pour le contrôle du LD résiduel ne sont pas
publiées. Cette agrégation fournit la provenance nécessaire aux figures tout en
limitant la diffusion d'informations génétiques détaillées.

## Critères bloquants

L'étape retourne le code V2 `4` si aucun variant n'est éligible, si le nombre
total de marqueurs est insuffisant, si un autosome requis est sous le minimum ou
si un chromosome concentre une fraction excessive du panel. Un échec PLINK
retourne le code `3` et une entrée ou un paramètre invalide le code `2`.

Le champ `effective_marker_count_operational` désigne uniquement le nombre de
marqueurs conservés après pruning. Il ne constitue pas une estimation spectrale
ou statistique de la dimension effective du panel.
