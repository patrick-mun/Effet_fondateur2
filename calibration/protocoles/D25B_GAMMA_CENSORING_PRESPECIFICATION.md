# D25b — Vraisemblance Gamma avec censure de type II

Version : `1.0.0`
État : `FINALIZED_BEFORE_RESULTS`
Date : 10 novembre 2026
Rapport : `CAL-REP-14B`

## 1. Motif et portée

D25a a établi trois faits avant D25b : la candidate Python V2 reproduit le
script R officiel, l'oracle indépendant sur points de recombinaison complets
est calibré, et l'extension de bord publiée ne calibre pas les longueurs
tronquées aux petits effectifs. D25b ne modifie donc ni V1 ni V2. Il évalue une
candidate distincte, `gamma_pairwise_type2_censored_v3`, dont la vraisemblance
représente explicitement le mécanisme d'observation paire-à-paire.

Le domaine qualifiable de D25b est une généalogie indépendante, des bras
exponentiels selon le modèle de Haldane, exactement un bras supérieur censuré
sur chaque côté et aucune correction de partage fortuit. Les généalogies
corrélées sont des tests de robustesse : elles ne peuvent pas être qualifiées
par cette vraisemblance indépendante.

## 2. Source biologique et contrat d'observation

Gandolfo, Bahlo et Speed (2014), DOI `10.1534/genetics.114.164616`, décrivent
chaque bras ancestral en Morgans comme une variable exponentielle de taux
`g`, l'âge en générations. Comme le partage n'est détectable qu'entre au moins
deux porteurs, le bras vrai le plus long dépasse généralement le deuxième plus
long, mais sa prolongation n'est pas observée. Sur chaque côté, les `n−1`
ruptures les plus courtes sont donc observées et le dernier bras est seulement
connu comme supérieur ou égal à la rupture `X_(n−1)`.

Le vecteur remis à la candidate conserve `n` expositions : les `n−1` valeurs
exactes et une copie de `X_(n−1)` pour le bras censuré. Il doit contenir au
moins deux familles, des longueurs strictement positives et un maximum répété
au moins deux fois sur chaque côté. D25b utilise 3, 5 et 10 familles.

## 3. Vraisemblance et estimateur préspecifiés

Pour un côté, avec `r=n−1` événements observés et une unité encore exposée au
temps `X_(r)`, le temps total sous risque est :

\[
T_s=\sum_{i=1}^{r}X_{(i)}+(n-r)X_{(r)}.
\]

Pour deux côtés indépendants, `T=T_L+T_R`, `k=2(n−1)` événements et :

\[
\ell(g)=k\log(g)-gT+C,\qquad \widehat g_{MLE}=\frac{k}{T}.
\]

Comme `gT ~ Gamma(k,1)`, l'estimateur sans biais et l'intervalle exact sont :

\[
\widehat g=\frac{k-1}{T},\qquad
IC_{1-\alpha}=\left[
\frac{G^{-1}_{k}(\alpha/2)}{T},
\frac{G^{-1}_{k}(1-\alpha/2)}{T}
\right].
\]

`G_k` est la fonction de répartition d'une Gamma de forme `k` et de taux 1.
Les longueurs d'entrée sont converties de cM en Morgans avant le calcul. La
candidate n'ajoute aucune longueur et refuse toute correction de partage
fortuit : ces deux mécanismes exigeraient une autre vraisemblance.

## 4. Panel indépendant de D25 et D25a

Les graines `25201` et `25202` sont nouvelles. La grille croise dix âges
(`5, 10, 15, 20, 30, 40, 50, 75, 100, 150` générations), trois effectifs
(`3, 5, 10` familles), trois profils généalogiques et 400 réplications par
cellule, soit 36 000 scénarios indépendants.

Le profil primaire `INDEPENDENT` simule les bras iid. Les profils
`CORRELATED_25` et `CORRELATED_50` reprennent le mécanisme de chemin partagé de
D25a avec de nouvelles graines ; ils servent uniquement à mesurer la rupture
d'hypothèse. Les côtés gauche et droit restent indépendants.

Chaque scénario produit quatre résultats avant ouverture de l'âge vrai :

- `TYPE2_CENSORED_INDEPENDENT_V3` sur l'observation tronquée ;
- `R_FAITHFUL_INDEPENDENT_V2` sur la même observation ;
- `R_FAITHFUL_CORRELATED_V2` sur la même observation ;
- `FULL_INDEPENDENT_ORACLE_V1` sur les ruptures complètes, comme contrôle.

Les copies de méthode sont appariées et ne constituent jamais des unités
indépendantes supplémentaires.

## 5. Mesures et critères fermés

Les résultats publient biais signé, IC95 bootstrap du biais, MAE, RMSE,
couverture exacte à 95 %, IC95 de Wilson, longueur moyenne d'intervalle,
résultats par âge et effectif, et stabilité leave-one-family-out. Les 1 000
bootstrap rééchantillonnent les scénarios, jamais les méthodes appariées.

Les critères sont fixés avant matérialisation :

- `CAL-AGE-B01` : erreur de la fixture analytique `≤ 1e-12` génération ;
- `CAL-AGE-B02` : biais absolu primaire agrégé `≤ 0,75` génération ;
- `CAL-AGE-B03` : couverture primaire agrégée entre `0,94` et `0,96` ;
- `CAL-AGE-B04` : biais absolu maximal entre effectifs `≤ 2,0` générations ;
- `CAL-AGE-B05` : couverture de chaque effectif entre `0,925` et `0,975` ;
- `CAL-AGE-B06` : aucun résultat primaire non estimé ;
- `CAL-AGE-B07` : les profils corrélés restent `REVIEW_REQUIRED` et ne peuvent
  provoquer une promotion automatique.

Un succès de B01–B06 qualifie uniquement le contrat indépendant censuré de
D25b. Il n'autorise pas à remplacer la méthode primaire corrélée du pipeline.
La promotion en production reste `NOT_PROMOTED` jusqu'à une décision humaine
et une solution dédiée à la dépendance généalogique.

## 6. Aveuglement, limites et décision

Les âges, réplications et graines par scénario restent sous `private_truth/`.
Le runner ne reçoit que les longueurs, l'effectif et le profil généalogique.
Tous les artefacts générés restent dans `calibration_workspace/`. Aucune donnée
`dock6_real`, aucun cache interne, aucun tirage D25/D25a et aucun pipeline
archivé ne sont utilisés.

D25b ne modélise pas la résolution des marqueurs, les erreurs de carte, le
phasage SHAPEIT5, le partage fortuit ni les faux positifs D24. Les égalités de
longueur créées par une grille de marqueurs réelle devront être distinguées de
la censure biologique avant toute utilisation. `calibrated=false` et
`founder_effect_proven=false` restent invariants.
