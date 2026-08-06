# Étape 11 — Préparation de la région cible

## Responsabilité

`11_prepare_target_region` extrait une fenêtre physique autour du variant cible
depuis `target_chromosome_all_qc`, associe chaque variant à une position
génétique issue d'une carte de référence et publie un jeu PLINK ordonné pour le
futur adaptateur de phasage.

Cette étape ne phase aucun génotype et ne choisit aucun logiciel de phasage.

## Carte d'entrée

`inputs.genetic_map` est un TSV conforme à `schemas/genetic_map.schema.json` :

```text
MAP_ID  ASSEMBLY  CHROMOSOME  POSITION_BP  POSITION_CM
```

La table doit utiliser un seul `MAP_ID`, déclarer exactement l'assemblage du
projet et fournir au moins deux ancres sur le chromosome cible. Sur ce
chromosome, les positions en bp doivent être strictement croissantes et les
positions cumulées en cM monotones.

L'origine, la population de référence et la méthode de construction associées
au `MAP_ID` doivent être validées avant une analyse scientifique réelle. Le
pipeline conserve l'identifiant et l'empreinte de la carte, mais ne peut pas en
évaluer la pertinence biologique.

## Fenêtre et interpolation

`window_left_bp` et `window_right_bp` définissent les distances physiques de
part et d'autre du variant cible. PLINK extrait cette fenêtre sans modifier la
liste d'individus du jeu final de l'étape 10.

Une position présente dans la carte reçoit `MAP_STATUS=EXACT`. Sinon, sa
position en cM est interpolée linéairement entre les deux ancres adjacentes.
L'intervalle d'ancres ne peut pas dépasser `max_interpolation_gap_bp`.

Les comportements suivants bloquent l'étape :

- variant hors des bornes de la carte ;
- carte génétique décroissante ;
- intervalle d'interpolation trop large ;
- doublon de position ou ordre physique non croissant ;
- allèle non A/C/G/T ou allèles identiques ;
- absence, mauvaise coordonnée ou mauvais couple REF/ALT du variant cible ;
- moins de `min_region_variants` après extraction.

Aucune extrapolation et aucune approximation `1 Mb = 1 cM` ne sont autorisées.

## Sorties

- `target_region.bed/.bim/.fam` : jeu régional, avec positions cM écrites dans
  la troisième colonne BIM ;
- `target_region.dataset.json` : descripteur et identifiants déterministes des
  ensembles d'individus et de variants ;
- `target_genetic_map.tsv` : ordre, positions bp/cM, ancres et taux local de
  chaque variant ;
- `phasing_input_manifest.json` : contrat générique
  `PLINK_BED_BIM_FAM_WITH_CM`, indépendant de l'outil de l'étape 12 ;
- `target_region_report.json`, `stage_outputs.json`, `audit.json` et
  `checksums.sha256`.

## Codes de retour

- `0` : région et carte publiées ;
- `2` : configuration, schéma, assemblage ou cohérence d'entrée invalide ;
- `3` : PLINK absent, en échec ou sortie incohérente ;
- `4` : région, carte ou variant cible scientifiquement inutilisable.
