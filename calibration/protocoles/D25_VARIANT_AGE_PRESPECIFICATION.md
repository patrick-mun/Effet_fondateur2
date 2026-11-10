# D25 — Préspecification de l’évaluation de la datation Gamma

Version : `1.0.0`
État : `FINALIZED_BEFORE_RESULTS`
Date : 10 novembre 2026
Rapport : `CAL-REP-14`

## 1. Objet et séparation des questions

D25 évalue l’étape `14_estimate_variant_age` contre des âges simulés connus.
Le primaire mesure le moteur Gamma sur des longueurs exactes, sans erreur de
phasage, de carte ou de délimitation. Les sensibilités modifient ensuite une
seule composante à la fois. Les sorties D24 sont reprises séparément pour
mesurer la propagation réelle des bornes de l’étape 13 ; elles ne sont jamais
mélangées au primaire.

Un âge Gamma est conditionnel à l’hypothèse d’une origine fondatrice unique.
Il ne prouve ni identité par descendance, ni effet fondateur, ni origine
géographique. D24 restant non qualifié, D25 ne peut pas rendre la chaîne
scientifique qualifiée par héritage.

## 2. Grille, unités et répétitions

La grille ferme les âges `5, 10, 15, 20, 30, 40, 50, 75, 100, 150`
générations et `3, 5, 10` familles indépendantes. Chaque cellule reçoit 200
réplications, soit 6 000 scénarios indépendants et 600 réplications par âge.
La graine de simulation est fixée avant résultats.

Pour un âge vrai `g`, chaque longueur familiale gauche ou droite est tirée
indépendamment selon :

\[
X_{f,s}\sim\operatorname{Exp}(g/100),
\qquad E[X_{f,s}]=100/g\ \mathrm{cM},
\]

où `f` indexe la famille et `s` le côté. Le primaire `ORACLE_EXACT` transmet
ces longueurs exactes au même moteur Python que la production. La vérité reste
fermée jusqu’à la production des estimations.

## 3. Modèles et profils

Le modèle `CORRELATED` de `gamma_gandolfo_2014_v1` est primaire, conformément
au pipeline. Le modèle `INDEPENDENT` historique est publié comme référence non
modifiée. Les profils appariés sont :

- `ORACLE_EXACT` : longueurs simulées exactes ;
- `WINDOW_3CM_CAP` : chaque bras est plafonné à 3 cM pour exposer la censure
  d’une fenêtre régionale de 6 cM ;
- `MAP_SCALE_0P95` et `MAP_SCALE_1P05` : toutes les longueurs sont multipliées
  par 0,95 ou 1,05 pour mesurer une erreur d’échelle de carte.

Les sorties D24 positives sont datées séparément sur `MAXIMUM`, `BASELINE`,
`DENSITY_50` et `MISSING_1PCT`. Elles utilisent un phasage-vérité et ne
quantifient donc pas l’erreur de phase SHAPEIT5. Ce volet reste explicitement
`NOT_EVALUABLE`.

## 4. Mesures

Pour chaque scénario évaluable, avec âge estimé `ĝ` et âge vrai `g` :

\[
e=\hat g-g,\quad
\mathrm{biais}=n^{-1}\sum e_i,\quad
\mathrm{MAE}=n^{-1}\sum|e_i|,\quad
\mathrm{RMSE}=\sqrt{n^{-1}\sum e_i^2}.
\]

La couverture vaut `1` si l’intervalle Gamma à 95 % contient `g`. Son
incertitude est l’intervalle de Wilson à 95 %. Les intervalles du biais, de la
MAE et de la RMSE proviennent de 1 000 bootstrap iid fixés par une graine
distincte. La stabilité leave-one-family-out est le maximum, par scénario, de
`|ĝ_{−f}−ĝ|`, puis sa moyenne et son quantile 95 %. Elle n’est pas calculable
pour trois familles, car retirer une famille laisse moins de trois unités.

Les métriques sont publiées globalement, par âge et par effectif familial. Les
profils appariés et la propagation D24 restent descriptifs ; ils ne multiplient
pas les dénominateurs primaires.

## 5. Référence historique et décisions

Le cas à six individus publié dans les tests de référence doit reproduire
l’estimation corrélée 22,3 et l’estimation indépendante 27,5 générations avec
une erreur absolue maximale de 0,1 génération (`CAL-AGE-006`). Les critères
`CAL-AGE-001` à `005` restent `pilot_then_seal` et reçoivent
`REVIEW_REQUIRED`, quelle que soit la valeur observée.

La conformité technique exige les schémas, empreintes, effectifs, grilles,
séparation vérité/estimations et référence historique. Elle ne vaut pas
qualification de performance. `calibrated=false` et
`founder_effect_proven=false` sont invariants.

## 6. Confidentialité, provenance et rapport

Tous les artefacts sont confinés à `calibration_workspace/`. Aucune donnée
`dock6_real`, aucun cache interne et aucun pipeline archivé ne sont autorisés.
CAL-REP-14 consomme seulement l’évaluation agrégée close. Il publie au minimum
quatre graphiques explicatifs, les équations réellement utilisées, les
effectifs indépendants, l’incertitude, les non-évaluations et les empreintes.
