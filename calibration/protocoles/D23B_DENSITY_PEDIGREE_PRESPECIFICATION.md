# D23b — Sensibilité contrôlée densité × pedigree, version 1.0.0

Fixé avant la matérialisation du panel maximal et avant tout nouveau phasage,
le 9 novembre 2026, sur demande explicite de l'utilisateur. Cette investigation
complète D23 sans modifier ses runs, évaluations, rapports ou verdicts.

## Question et portée

D23b mesure séparément, sur les mêmes 40 sujets synthétiques, l'association
entre densité de marqueurs et erreurs de phase, puis l'association entre usage
du pedigree et erreurs de phase régionale. Elle ne prétend pas identifier une
cause biologique, qualifier SHAPEIT5, sceller un seuil `pilot_then_seal` ou
valider l'étape 13.

La cible reste l'annotation technique APC du chromosome 5 employée par D23a,
pas DOCK6. Les 32 fondateurs, huit trios, seize sujets isolés, haplotypes
donneurs, méioses, référence D10 de 1 403 sujets, cartes, graine de phasage
2303, paramètres SHAPEIT5 et vérité privée sont conservés. Aucun profil de
génotypes manquants n'est ajouté : la densité et le pedigree sont les seuls
facteurs modifiés.

## Détermination du maximum exploitable

Le maximum exploitable est défini avant phasage comme l'ensemble de tous les
SNP autosomiques CytoScan 750K Accel NA36.r1 qui satisfont simultanément les
gardes D23a : position unique, carte GRCh38 couvrante sans extrapolation,
orientation compatible non ambiguë, variant biallélique SNP présent dans la
source 1000G, génotypes complets chez les 64 donneurs D10 utilisés et
harmonisation REF/ALT non ambiguë. La cible technique supplémentaire est
conservée uniquement dans le parcours régional et n'est pas comptée comme une
sonde de puce.

La configuration de préparation ne sous-échantillonne aucun candidat : les
plafonds `candidate_markers_per_autosome` et
`maximum_markers_per_autosome` valent 1 000 000, au-dessus du catalogue
autosomal de la puce. Le nombre obtenu est un résultat d'harmonisation, pas un
objectif ajustable. Une source, carte ou empreinte différente invalide le
maximum et exige une nouvelle version.

## Profils de densité emboîtés

Dans chaque chromosome, les loci harmonisés sont ordonnés par position
physique, hors cible technique ajoutée. Pour un rang zéro-indexé `r`, les
ensembles sont :

- `DENSITY_25` : `r mod 4 = 0` ;
- `DENSITY_50` : `r mod 4 ∈ {0, 2}` ;
- `DENSITY_75` : `r mod 4 ∈ {0, 1, 2}` ;
- `DENSITY_100` : tous les rangs.

Ainsi, `DENSITY_25 ⊂ DENSITY_50 ⊂ DENSITY_75 ⊂ DENSITY_100`. La cible
technique est toujours conservée dans le chromosome 5. Aucun locus n'est
choisi selon une sortie de phasage. Les quatre profils sont des observations
appariées, pas quatre réplications indépendantes.

## Facteur pedigree et parcours

Le parcours régional 11/12 est exécuté pour chaque densité dans deux modes :

- `WITH_PEDIGREE` : PID/MID synthétiques D23a transmis à SHAPEIT5 ;
- `WITHOUT_PEDIGREE` : mêmes sujets, FID, génotypes et ordre, mais PID/MID
  remplacés par `0` dans une fixture dérivée et auditée.

Le parcours autosomal 10A/12A est exécuté une fois par densité. Son adaptateur
de production ne transmet pas de pedigree ; son mode est donc explicitement
`WITHOUT_PEDIGREE_BY_CONTRACT`. D23b ne modifie pas cet adaptateur pour créer
une comparaison artificielle.

Le primaire régional est `DENSITY_100 × WITH_PEDIGREE`. Le primaire autosomal
est `DENSITY_100 × WITHOUT_PEDIGREE_BY_CONTRACT`. Ce choix est indépendant des
résultats. Les autres cellules sont des sensibilités descriptives.

## Métriques, contrastes et incertitude

Les définitions D23 sont inchangées : taux d'erreur de commutation, Hamming de
phase, erreur rare MAF fondatrice ≤ 0,05, motifs isolés et orientation cible.
Les génotypes et loci sont alignés par chromosome, bp, REF/ALT et sujet, jamais
par simple rang.

Pour une métrique de taux dans l'unité familiale `u`, on conserve le numérateur
`e_u` et le dénominateur `n_u`. Le taux agrégé est :

`R = Σ_u e_u / Σ_u n_u`.

Les contrastes sont des différences de taux, en points de proportion :

- effet descriptif de densité :
  `Δ_density(d, p) = R(d, p) − R(DENSITY_100, p)` ;
- effet descriptif du retrait du pedigree :
  `Δ_pedigree(d) = R(d, WITHOUT_PEDIGREE) − R(d, WITH_PEDIGREE)`.

Un contraste positif indique davantage d'erreurs dans le premier terme. Les
IC descriptifs à 95 % utilisent 2 000 rééchantillonnages appariés des 24 unités
(huit familles entières et seize sujets isolés), graine 2311. Le même tirage
d'unités est appliqué aux deux termes du contraste. Un dénominateur nul dans
l'un des termes exclut le tirage et est compté. Aucun IC n'est causal : densité
et pedigree sont manipulés, mais ce petit panel et ce seul logiciel ne couvrent
pas tous les mécanismes possibles.

## Contrôles et décision

Les contrôles fixes exigent : sources et empreintes D10/D23a, mêmes sujets et
ordre, ensembles emboîtés, cible conservée, zéro discordance de dosage, cartes
concordantes, nombre de lignes de pedigree 8 en mode activé et 0 en mode retiré,
et absence de vérité dans les entrées des outils.

Une exécution d'outil échouée reste `FAIL`; une métrique sans dénominateur reste
`NOT_EVALUABLE`. Les métriques mesurables restent `REVIEW_REQUIRED`, car aucun
seuil D04 n'est scellé. Une densité plus élevée n'est pas déclarée meilleure
par principe et le parcours primaire n'est pas remplacé après observation.

Le rapport indépendant, produit avec `report-formater`, comporte les effectifs
par chromosome et densité, les taux avec IC, les contrastes densité et pedigree,
les équations ci-dessus, les limites d'indépendance et les empreintes. Le rendu
ne recalcule aucune statistique. `calibrated=false` et
`founder_effect_proven=false` sont obligatoires.
