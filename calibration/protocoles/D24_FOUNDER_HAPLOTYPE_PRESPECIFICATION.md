# D24 — Évaluation de l’haplotype fondateur, version 1.0.0

Préparé le 10 novembre 2026 avant toute exécution D24, conformément à
`CAL-DEC-010`. Ce document fixe l’évaluateur de l’étape 13 ; il n’autorise ni
la matérialisation d’un nouveau panel, ni l’ouverture du panel de validation,
ni un changement de seuil scientifique.

## Question et portée

D24 mesure si `13_infer_founder_haplotype` récupère, autour d’une cible
synthétique, le segment IBS commun injecté par D14. L’évaluateur compare les
sorties closes de l’étape 13 à une vérité privée indépendante. Il ne relance
pas le pipeline et ne reçoit aucune donnée `dock6_real` ni aucun cache interne.

Le parcours primaire utilise la densité maximale exploitable déterminée en
D23b. `BASELINE`, `DENSITY_50` et `MISSING_1PCT` restent des sensibilités
appariées. Les mêmes scénarios, familles, cibles et graines doivent être
conservés entre profils. La qualité du phasage D23b reste
`REVIEW_REQUIRED` : D24 mesure l’effet résultant sans requalifier D23.

## Porte d’entrée obligatoire

Chaque scénario doit fournir avant ouverture de la vérité :

- une sortie close de l’étape 13 conforme aux quatre contrats de production ;
- au moins trois unités familiales indépendantes ;
- une cible SNV synthétique et une carte en bp/cM empreintées ;
- un rôle fixé `POSITIVE` ou `NEGATIVE` ;
- pour un positif, au moins trois segments vrais familiaux et une décision
  `detectable` fixée avant les sorties ;
- pour un négatif, aucun segment fondateur injecté et `detectable=false` ;
- une empreinte SHA-256 de la vérité et des sorties évaluées.

Le run D14 chromosome 5 existant ne satisfait pas cette porte : il contient
deux familles seulement. Les paires D15 ne le complètent pas, car leur portée
est une injection IBD technique par paire et leur IBD de fond est inconnu.
D23a fournit une vérité de phase, pas une vérité de segment fondateur.

Le banc de développement D24 exigera au minimum cinq scénarios positifs et
cinq négatifs. Ces effectifs servent à éprouver l’évaluateur ; ils ne
constituent pas le pilote D34 ni un seuil de qualification.

La matrice fermée comprend cinq positifs (50, 100, 150, 200 générations et
une seconde réplication à 100) et cinq négatifs à graines distinctes. Un
négatif porte la cible sur au moins trois lignées terminales indépendantes,
sans copier d’intervalle fondateur entre elles
(`INDEPENDENT_TARGET_CARRIERS_NO_SHARED_SEGMENT`). Il teste donc un faux
consensus autour de vrais porteurs synthétiques, et non le cas trivial d’une
cible absente. Chaque scénario biologique est rejoué sur les quatre profils,
soit 40 cas appariés. Les profils ne sont jamais comptés comme 40 réplications
biologiques indépendantes.

## Vérité et décision de détection

Pour un scénario positif comportant les intervalles familiaux vrais
`[L_f, R_f]`, le segment commun vrai est leur intersection :

`L* = max_f(L_f)` et `R* = min_f(R_f)`.

Le scénario n’est détectable que si cette intersection contient la cible et
si la décision privée `detectable` vaut vrai. Aucun tirage censuré ou difficile
n’est retiré après lecture des résultats.

Une détection prédite est définie exclusivement par
`status=SUPPORTED_IBS_CANDIDATE` et
`interpretation=IBS_SHARED_CANDIDATE`. Toute revendication IBD par l’étape 13
est un défaut de contrat, même si les bornes sont exactes.

## Métriques préspecifiées

Avec `TP`, `FN`, `TN` et `FP` comptés par scénario indépendant :

`sensibilité = TP / (TP + FN)`

`spécificité = TN / (TN + FP)`

Les critères principaux utilisent uniquement les dix scénarios indépendants
du profil `MAXIMUM`. Les trois autres profils décrivent la sensibilité de la
méthode sur exactement les mêmes vérités ; ils ne sont ni additionnés au
primaire ni traités comme des réplications. Les intervalles descriptifs à 95 %
utilisent Wilson, par scénario et non par marqueur. Un dénominateur nul produit
`NOT_EVALUABLE`, jamais zéro.

Pour chaque positif défini détectable avant analyse, la cible est contenue si
`L_prédit ≤ position_cible ≤ R_prédit`. Les erreurs signées et absolues des
bornes gauche et droite sont publiées séparément en bp et cM :

`e_L = L_prédit − L*` et `e_R = R_prédit − R*`.

Une absence de segment compte comme échec de containment. Une borne vraie
censurée reste visible mais n’est pas traitée comme un événement exact dans
l’erreur de borne. Avec `m` bornes vraies non censurées évaluables :

`MAE_bornes = Σ_j |e_j| / m`.

Les erreurs familiales sont également agrégées après appariement par famille,
sans publier les identifiants. Une famille absente, exclue ou mal phasée reste
comptée avec son motif et ne devient pas une borne correcte.

## Critères et décisions

Les cinq critères `CAL-FOUND-001` à `CAL-FOUND-005` sont repris sans
modification. Les quatre métriques `pilot_then_seal` restent
`REVIEW_REQUIRED` tant que D34–D35 n’ont pas fixé leurs seuils. Le critère fixe
`CAL-FOUND-005` exige zéro requalification IBS→IBD.

Un `PASS` D24 signifie seulement que l’évaluateur, les contrats, les calculs et
le rendu fonctionnent sur le banc de développement. Il ne qualifie pas la
performance de l’étape 13. `calibrated=false` et
`founder_effect_proven=false` sont obligatoires.

## Rapport indépendant

`CAL-REP-13` est produit en français avec `report-formater` sous forme HTML
autonome imprimable A4, puis PDF. Il comporte quatorze sections, au moins trois
graphiques quantitatifs, les équations ci-dessus, les numérateurs,
dénominateurs, intervalles, effectifs indépendants, anomalies, limites et
empreintes. Les graphiques ne recalculent aucune métrique et ne contiennent ni
identifiant individuel ni table de vérité.

Le rapport doit rappeler que l’étape 13 détecte un candidat IBS et ne prouve
pas une transmission IBD. Le FAIL 07 de D21 et les limites de phasage D23b
restent visibles et ne sont pas résolus par D24.
