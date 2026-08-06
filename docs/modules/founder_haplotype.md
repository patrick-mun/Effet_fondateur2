# Haplotype fondateur candidat — étape 13

## Portée scientifique

L'étape 13 teste si les chromosomes porteurs fiables d'unités familiales
indépendantes partagent un haplotype autour du variant cible. La méthode primaire
`target_centered_exact_ibs_v1` mesure une identité allélique **IBS**. Elle ne
prouve pas une transmission **IBD** et ne conclut donc pas seule à un fondateur
unique.

Un appel IBD nécessiterait un adaptateur distinct, versionné et validé sur une
densité de marqueurs adaptée. La configuration `local_ibd_adapter: null` signifie
explicitement qu'aucun appel IBD n'est effectué.

## Sélection

- La cohorte primaire est `target_carriers_independent` gelée à l'étape 09.
- Seules les attributions `H1` ou `H2` marquées `PASS` à l'étape 12 sont incluses.
- Une phase non fiable reste visible avec `PHASING_UNRELIABLE`.
- Un porteur homozygote cible `BOTH` compte comme une seule unité familiale. Son
  segment est défini par le partage continu de marqueurs homozygotes identiques,
  conformément au cas récessif de la méthode Gamma ; les marqueurs hétérozygotes
  ne sont pas transformés artificiellement en allèle ancestral.
- La copie ALT du GT cible phasé est recontrôlée avant l'analyse.

## Détection du segment

Les variants sont ordonnés en bp et associés à la carte génétique de l'étape 11.
Depuis la cible, l'algorithme avance séparément à gauche et à droite tant que
tous les chromosomes porteurs sélectionnés possèdent le même allèle appelé. Le
premier manque ou la première discordance arrête le segment. Aucune tolérance de
mismatch ni interpolation allélique n'est appliquée.

Les limites individuelles sont les extensions maximales gauche et droite que le
chromosome partage continûment avec au moins un autre chromosome porteur
indépendant ; les deux bras peuvent avoir des comparateurs différents. Le
consensus primaire reste, lui, l'intersection exacte de tous les porteurs. Les
limites sont inclusives. `LEFT_LENGTH_CM` et `RIGHT_LENGTH_CM` sont des distances
positives ou nulles depuis la cible, jamais des positions absolues.
Le profil générique exige au moins trois unités indépendantes et deux marqueurs
informatifs de chaque côté. Une insuffisance publie `NO_FOUNDER_CONCLUSION`.

La matrice de partage contient la diagonale `SELF` et chaque paire indépendante.
Elle sert aux analyses de sensibilité ultérieures ; le résultat primaire reste
l'intersection exacte de tous les chromosomes retenus.

## Fréquence de fond

La signature consensus est comparée aux deux haplotypes appelés des témoins
indépendants et des non-porteurs QC du chromosome cible. Le variant cible est
exclu de cette signature : sa présence chez les porteurs et son absence chez les
non-porteurs ne doivent pas créer artificiellement une fréquence nulle.

## Sorties

- `founder_segments.tsv` : une ligne par unité porteuse indépendante, exclusions
  conservées ;
- `founder_consensus.tsv` : intervalle primaire, effectifs et fréquence de fond ;
- `founder_sharing_matrix.tsv` : segments IBS pairwise ;
- `founder_discordances.tsv` : statut allélique de chaque marqueur ;
- `founder_analysis_summary.json` : décision synthétique sans identifiant ;
- `audit.json` : paramètres, versions, limites et empreintes du run.

Toutes les tables individuelles sont `sensitive_genetic`. Le résumé est
`internal`. Une revue manuelle est toujours requise avant interprétation ou
datation.
