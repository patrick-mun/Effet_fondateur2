# Étape 11 — Préparation de la région cible

## Responsabilité

`11_prepare_target_region` extrait une fenêtre physique autour du variant cible
depuis `target_chromosome_all_qc`, associe chaque variant à une position
génétique issue d'une carte de référence et publie un jeu PLINK ordonné pour le
futur adaptateur de phasage.

Cette étape ne phase aucun génotype et ne choisit aucun logiciel de phasage.

## Carte d'entrée et cache partagé

Deux modes exclusifs sont acceptés :

- `inputs.genetic_map` fournit directement un TSV validé ;
- `inputs.genetic_map_catalog` désigne un catalogue versionné, tandis que
  `genetic_map_id` choisit l'archive publique épinglée.

Le TSV normalisé est conforme à `schemas/genetic_map.schema.json` :

```text
MAP_ID  ASSEMBLY  CHROMOSOME  POSITION_BP  POSITION_CM
```

La table doit utiliser un seul `MAP_ID`, déclarer exactement l'assemblage du
projet et fournir au moins deux ancres sur le chromosome cible. Sur ce
chromosome, les positions en bp doivent être strictement croissantes et les
positions cumulées en cM monotones.

En mode catalogue, l'archive est téléchargée sous
`genetic_map_cache_dir`, contrôlée par SHA-256, extraite dans un dossier
temporaire puis publiée atomiquement. Les 22 cartes autosomiques sont
normalisées au premier téléchargement. L'archive d'origine, les cartes et un
manifest de leurs empreintes restent conservés : une nouvelle mutation du même
chromosome produit un cache `HIT`, et les autres autosomes sont déjà prêts.
`genetic_map_cache_offline=true` interdit tout accès réseau et bloque si le
cache manque. Un verrou empêche deux runs de peupler la même entrée en parallèle.

Le catalogue fait partie de la signature du run. Une incohérence du catalogue,
de l'assemblage, de l'archive, du manifest ou d'une carte normalisée bloque
l'étape. Le cache n'est jamais une source génétique du projet et ne contient
aucun échantillon. L'origine, la population de référence et la méthode restent
dans l'audit ; le pipeline ne peut pas en évaluer seul la pertinence biologique.

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
- archive ou fichier du cache dont l'empreinte diffère du catalogue/manifest ;
- archive absente lorsque le mode hors ligne est demandé.

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
