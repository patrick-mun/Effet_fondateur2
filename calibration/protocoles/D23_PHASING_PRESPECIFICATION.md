# D23 — Évaluation du phasage, version 1.0.0

Fixé avant le premier phasage D23, le 30 octobre 2026, sur demande explicite
de l'utilisateur. Source exclusive : D23a v2, manifest SHA-256
`e91d686f5cc63d17a95f64dc52bb481d0836c898f23d03344a801b59b60d33cb`.
Le brouillon v1 est interdit. Le panel aveugle de validation reste fermé.

## Domaine et exécution

Banc de composants réels 10A, 11, 12 et 12A : PLINK, bcftools et SHAPEIT5
5.1.1, un thread de phasage, graine 2303, timeout explicite par commande.
La conversion VCF–PLINK et les contrats d'entrée QC/cohorte sont des adaptateurs
de banc identifiés `D23_UPSTREAM_CONTRACT_FIXTURE`, pas des exécutions 03–09.
Ils ne qualifient pas les étapes amont, notamment le point 07 laissé ouvert.
Les génotypes ne sont jamais déduits des statuts ou groupes.

Trois profils appariés, choisis avant les résultats : BASELINE ; DENSITY_50
(un marqueur sur deux en rang par chromosome, cible conservée) ; MISSING_1PCT
(masquage diploïde Bernoulli 0,01, graine 2304, cible conservée). Aucun profil
n'est une réplication indépendante. Les seuils de densité de production de
10A/12A ne sont pas abaissés pour obtenir un succès. Les échecs restent visibles.
12 utilise la fenêtre cible ±5 Mb ; 12A couvre les 22 autosomes.

Référence exacte D10 (1 403 sujets), vérifiée après extraction. Le cache
exclusif du banc est sous le workspace ; aucun cache interne d'étude n'est lu.
Le chromosome 5 public complet peut être copié depuis le Bureau dans ce cache
par l'API de production avec vérification MD5/SHA ; les autres acquisitions
12A utilisent directement les sources officielles du Bureau. Les cartes D23a
épinglées sont réutilisées sans modification. Ni les haplotypes vrais ni les
transmissions privées ne sont accessibles dans les commandes de phasage.

## Oracles indépendants et contrôles fixes

10A : l'ensemble attendu est calculé depuis les observations, par complétude
sur les 40 sujets et polymorphisme (MAF ≥ 0,000001), sans lire son résumé PASS.
La fréquence utilise les 32 fondateurs : PLINK n'inclut pas les non-fondateurs
sans `--nonfounders`, comme confirmé par l'aide locale de l'exécutable.
Comparer ordre des sujets, ensemble et ordre des loci, dosages et 22 autosomes.
11 : oracle de fenêtre physique et interpolation indépendante (NumPy) sur la
carte officielle ; contrôler bp/cM, cible unique, ordre et dosages. Tolérance
numérique absolue 1e-8 cM, hors performance scientifique.
12/12A : alignement par chromosome, position, REF/ALT et identifiant de sujet,
jamais par simple rang. Une permutation allélique est recodée explicitement.
Discordance de dosage, allèle incompatible, manque ou absence sont comptés,
jamais convertis en succès de phase. La cible ajoutée est attendue dans 12,
mais son absence dans 12A est une exclusion de référence déclarée, non un échec
de conservation de cible applicable à 12A.

## Métriques

Pour les hétérozygotes vrais, complets, phasés et de dosage concordant,
`o_j = 1[premier_allèle_prédit ≠ premier_allèle_vrai]`.
Une inversion globale H1/H2 est indétectable biologiquement et n'est pas une
erreur. Les homozygotes ne sont pas des observations de phase.

`SER = Σ 1[o_j ≠ o_(j-1)] / nombre_de_transitions_évaluables`.
Ne jamais relier deux chromosomes, blocs PS différents, ni deux hétérozygotes
séparés par un hétérozygote exclu. Les homozygotes intermédiaires sont ignorés.

`Hamming_phase = min(Σ o_j, n_het − Σ o_j) / n_het`, par sujet/chromosome.
L'alignement global est choisi sur tous les hétérozygotes sauf la cible,
séparément dans chaque bloc PS lorsqu'il existe (sinon par chromosome).
Une égalité rend l'orientation ponctuelle indéterminée, pas correcte.
Le taux rare applique cet alignement aux hétérozygotes MAF ≤ 0,05 ; MAF
calculée sur les 32 fondateurs vrais uniquement (64 copies), non sur les enfants.
Un point isolé est un état différent de ses deux voisins évaluables égaux ;
compter séparément ces motifs et les transitions. Cette description de motifs
ne prétend pas identifier la cause de l'erreur ni séparer sans ambiguïté un
bloc court d'une erreur ponctuelle.

Orientation cible : précision parmi hétérozygotes cibles dont l'ancrage
hors cible a une majorité stricte ; publier aussi éligibles, non évaluables
et couverture. Les homozygotes n'ont pas d'orientation porteuse unique.
Une mutation ajoutée à fréquence artificielle n'évalue pas la rareté biologique.

Strates : MAF fondatrice ≤0,05 / >0,05 ; densité de voisinage (distance
physique au marqueur précédent ≤100 kb / >100 kb) ; profil de manque 0 / 1 %.
Le premier marqueur d'un chromosome n'a pas d'écart précédent : il reste dans
ALL et dans sa strate MAF, mais n'est attribué à aucune strate d'écart.
Pour les transitions, la strate est celle du second locus. Toutes les cellules
ont numérateur, dénominateur et raison explicite si dénominateur nul.

IC descriptifs à 95 % : bootstrap de 24 unités (8 trios +16 sujets isolés),
2 000 rééchantillonnages, graine 2305. Rééchantillonner les familles entières,
pas les marqueurs, chromosomes ou les profils appariés. Publier les unités
contributrices ; un dénominateur bootstrap nul est exclu et compté. Aucun IC
ne constitue une garantie hors du sous-panel de puce ni une validation externe.

## Contrôles de l'évaluateur et décision

Tests synthétiques à vérité connue : phase parfaite, inversion globale,
commutation de bloc, erreur isolée, homozygotes, génotypes manquants,
discordance de dosage, rupture PS, cible seule/égalité d'ancrage, permutations
des sujets/loci, allèles inversés, doublons et absence de cible.
Contrôles d'intégrité : refus de symlink, run existant, source modifiée,
inventaire incomplet, liste de référence divergente et fuite de vérité.

Rapports indépendants 10A/11/12/12A, chacun avec résultats quantitatifs,
graphiques, équations et limites. Les critères fixes sont PASS/FAIL ; les
métriques mesurables restent évaluables même sans seuil scellé. La qualification
de performance reste REVIEW_REQUIRED (`pilot_then_seal`), jamais PASS par
défaut. Un exécutable échoué reste FAIL et ses métriques absentes restent
NOT_EVALUABLE. `calibrated=false`, `founder_effect_proven=false`.
