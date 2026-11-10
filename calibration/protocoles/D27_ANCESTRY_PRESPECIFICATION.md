# D27 — Préspecification fermée de l’évaluateur d’ascendance 16A

Statut : `FINALIZED_BEFORE_RUN`
Autorisation : demande utilisateur du 10 novembre 2026
Rapport : `CAL-REP-16A`

## 1. Question et portée

D27 évalue le cœur réel de l’étape `16A_analyze_reference_ancestry` : PCA
ajustée exclusivement sur une référence, harmonisation allélique et projection
indépendante des sujets d’étude. Il ne réexécute pas toute la chaîne amont de
phasage. Le banc est régional sur le chromosome 5 ; son résultat diploïde est
donc un **proxy régional du positionnement global**, et non une validation de
la PCA génomique de production.

Les 1 403 références appartiennent à la partition D10
`PHASING_REFERENCE_SHARED`, disjointe des 900 donneurs D10 ayant produit D12.
Elles sont lues dans les actifs publics matérialisés de D23b, jamais dans le
cache interne. Aucune donnée `dock6_real` ou `data/input` n’est autorisée.

## 2. Appels fermés avant ouverture de la vérité

La phase d’appel lit le VCF aveugle D12 et, pour D17, les huit observations
appariées `control`, `missing_01`, `missing_05`, `genotype_01`,
`phase_point_02`, `phase_switch_02`, `density_50` et `combined`. Elle ajuste
les axes sur la seule référence indépendante puis scelle coordonnées,
valeurs propres, variants informatifs et empreintes. Les vérités D12 et D14
ne sont ouvertes qu’ensuite.

Les dosages diploïdes évaluent le positionnement régional. Les haplotypes D17
évaluent le positionnement local conditionnel à la phase candidate. Ce second
volet ne valide ni SHAPEIT5, ni l’étape 16A bout-en-bout.

## 3. Mélange simulé

La métrique primaire est limitée au contraste à deux strates effectivement
simulé en D12. Les codes sources sont remplacés dans le rapport par
`REFERENCE_STRATUM_A` et `REFERENCE_STRATUM_B`; ils désignent des strates du
simulateur, pas une identité ethnique ou géographique.

Les scores sont standardisés par la racine de la valeur propre :

\[
z_{ik}=s_{ik}/\sqrt{\lambda_k}.
\]

Avec les centroïdes de référence \(c_A\) et \(c_B\), la proportion estimée
sur le segment joignant les deux repères est :

\[
\hat\alpha_i=\operatorname{clip}_{[0,1]}
\frac{(z_i-c_B)^\top(c_A-c_B)}{\lVert c_A-c_B\rVert_2^2}.
\]

La métrique `CAL-ANC-001` est :

\[
\operatorname{MAE}_{mix}=n^{-1}\sum_i|\hat\alpha_i-\alpha_i|.
\]

Elle porte sur les sujets purs A/B et les mélanges réciproques 50/50. Un
intervalle bootstrap percentile à 95 %, 2 000 réplications et graine 2701 est
publié. Les sujets, et non les dix composantes, sont les unités de réplication.

## 4. Stabilité des positions répétées

Chaque profil est ajusté séparément comme en production. Pour rendre les axes
comparables, les scores standardisés des mêmes références sont alignés sur le
contrôle par transformation orthogonale de Procrustes. Aucun sujet simulé ne
participe à cet alignement.

Pour une entité \(i\), avec \(K=10\) composantes :

\[
d_i=\sqrt{K^{-1}\sum_{k=1}^{K}
(z^{aligné}_{ik,p}-z_{ik,0})^2}.
\]

Au niveau local, les deux copies d’un sujet sont appariées à celles du
contrôle par la permutation minimisant la somme des distances, afin de ne pas
confondre orientation arbitraire des copies et déplacement de position. La
métrique `CAL-ANC-002` est le maximum, parmi profils non témoins et portées
diploïde/locale, de la moyenne des \(d_i\). Les moyennes, médianes et 95e
percentiles par profil restent publiés.

## 5. Vérité locale et interprétation

Les six copies fondatrices D14 sont contrôlées descriptivement par rapport aux
deux centroïdes attendus dans la PCA locale du contrôle. Elles proviennent
d’un seul haplotype fondateur et ne constituent donc pas six réplications
indépendantes. Le positionnement local ne prouve ni ascendance locale, ni
ancêtre généalogique, ni IBD.

## 6. Verdicts autorisés

Le statut technique peut être `PASS_CORE_EVALUATOR` si toutes les sources sont
liées par SHA-256, les huit profils sont projetés, la séparation référence /
étude est vérifiée et les sorties sont scellées. Ce statut ne qualifie pas
l’intégration bout-en-bout de 16A.

Les seuils `CAL-ANC-001` et `CAL-ANC-002` sont `pilot_then_seal` et valent
toujours `null`. Leur statut reste `REVIEW_REQUIRED` quel que soit le résultat.
La qualification de performance est obligatoirement
`REVIEW_REQUIRED_UNSEALED`; `production_promotion=false`, `calibrated=false`
et `founder_effect_proven=false`.
