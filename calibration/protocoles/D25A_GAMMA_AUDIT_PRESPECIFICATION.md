# D25a — Audit de fidélité et de domaine de la datation Gamma

Version : `1.0.0`
État : `FINALIZED_BEFORE_NEW_RESULTS`
Date : 10 novembre 2026
Rapport : `CAL-REP-14A`

## 1. Motif et portée

D25a répond à l'anomalie révélée par D25 sans réinterpréter ni écraser ses
résultats. D25 transmettait au moteur Gamma des distances complètes jusqu'au
premier événement de recombinaison, puis appliquait la correction de bord
destinée aux deux bras les plus longs tronqués par le critère de partage. Il
appliquait aussi le modèle corrélé à des généalogies simulées indépendantes.

D25a sépare donc l'observation biologique, le mécanisme de mesure et le modèle
généalogique. Il vérifie en outre la fidélité de l'implémentation Python au
script R officiel de Gandolfo, Bahlo et Speed (2014). La méthode historique
`gamma_gandolfo_2014_v1` reste inchangée et visible. La candidate reçoit
l'identifiant distinct `gamma_gandolfo_2014_r_faithful_v2`.

## 2. Sources de référence et contrôle R

Le script officiel est
`https://bioinf.wehi.edu.au/software/ageofmutation/Mutation_Age_estimation.R`,
d'empreinte SHA-256
`9cdf7008aae60128005948299a87175ab23f54f93ddd1351d4d374359d20009e`.
Il est exécuté localement avec R sur cinq fixtures : l'exemple officiel à six
porteurs et quatre jeux synthétiques couvrant les intervalles de corrélation
`I1`, `I2`, `I3` et non négatif. La différence absolue maximale admise entre R
et la candidate Python est `1e-10` génération pour les estimations et bornes.

Le script R utilise un effectif `n*` borné pour la correction de biais, puis un
effectif éventuellement différent et strictement positif pour paramétrer
l'intervalle Gamma. Ces deux quantités doivent rester distinctes dans la
candidate. L'exemple historique seul ne suffit pas, car sa corrélation estimée
est positive.

## 3. Grille indépendante des résultats D25

Les nouvelles graines `25101` et `25102` sont distinctes de D25. La grille
croise dix âges (`5, 10, 15, 20, 30, 40, 50, 75, 100, 150` générations), trois
effectifs (`3, 5, 10` familles), trois profils généalogiques et 200 réplications
par cellule, soit 18 000 scénarios indépendants.

Pour un âge `g` et une fraction de chemin commune `q`, chaque côté reçoit une
distance commune `C` et des distances privées `P_f` :

\[
C\sim\operatorname{Exp}(gq/100),\qquad
P_f\sim\operatorname{Exp}(g(1-q)/100),\qquad
X_f=\min(C,P_f).
\]

Ainsi `X_f ~ Exp(g/100)` conserve la marginale attendue, tandis que `C` crée
une dépendance généalogique positive entre familles. `q=0` utilise directement
des tirages indépendants ; `q=0,25` et `q=0,50` définissent les deux profils
corrélés. Les côtés gauche et droit sont simulés indépendamment.

## 4. Deux contrats d'observation

`FULL_BREAKPOINTS` contient les distances biologiques complètes. Aucune
correction de bord n'y est autorisée. `PAIRWISE_SHARED_TRUNCATED` représente le
segment observable par partage entre au moins deux porteurs : sur chaque côté,
le bras le plus long est tronqué à la longueur du deuxième plus long. La
correction de bord Gamma est alors appliquée exactement une fois.

Les longueurs et le profil généalogique sont remis au runner sans âge vrai.
L'âge, la réplication et les paramètres de vérité restent dans
`private_truth/` jusqu'à fermeture des estimations.

## 5. Matrice des méthodes

Les estimations suivantes sont publiées en parallèle :

- `EXACT_INDEPENDENT_ORACLE_V1` : estimateur indépendant sans correction de
  bord, référence primaire uniquement pour `INDEPENDENT × FULL_BREAKPOINTS` ;
- `R_FAITHFUL_INDEPENDENT_V2` : version indépendante du script R, primaire pour
  `INDEPENDENT × PAIRWISE_SHARED_TRUNCATED` ;
- `R_FAITHFUL_CORRELATED_V2` : version corrélée fidèle au script R, primaire
  pour les profils corrélés, avec ou sans correction selon le contrat
  d'observation ;
- `LEGACY_INDEPENDENT_V1` et `LEGACY_CORRELATED_V1` : résultats historiques
  inchangés, toujours secondaires.

Aucune méthode ne peut recevoir une correction de bord sur `FULL_BREAKPOINTS`.
Les contrastes sont appariés par scénario, mais les copies de méthode ne sont
pas comptées comme réplications indépendantes supplémentaires.

## 6. Mesures et critères

Sont publiés par généalogie, observation, méthode, âge et nombre de familles :
biais signé, MAE, RMSE, couverture à 95 %, longueur moyenne de l'intervalle et
nombre non évaluable. Les intervalles du biais, de la MAE et de la RMSE
utilisent 1 000 bootstrap iid ; la couverture reçoit un intervalle de Wilson.

Les contrôles fixes sont :

- `CAL-AGE-A01` : parité R/Python maximale `≤ 1e-10` génération ;
- `CAL-AGE-A02` : biais absolu agrégé de l'oracle exact indépendant `≤ 0,75`
  génération ;
- `CAL-AGE-A03` : couverture agrégée de cet oracle entre `0,94` et `0,96` ;
- `CAL-AGE-A04` : biais absolu agrégé de la candidate indépendante sur les
  longueurs tronquées `≤ 0,75` génération ;
- `CAL-AGE-A05` : couverture agrégée correspondante entre `0,94` et `0,96`.

Les performances corrélées restent `pilot_then_seal` et donc
`REVIEW_REQUIRED`, quelle que soit leur valeur. Un succès technique de D25a ne
change pas automatiquement la méthode primaire du pipeline.

## 7. Confidentialité, décision et rapport

Tous les artefacts générés restent sous `calibration_workspace/`. Aucune donnée
`dock6_real`, aucun cache interne, aucun résultat D25 et aucun pipeline archivé
ne sert à générer le nouveau panel. CAL-REP-14A consomme uniquement l'évaluation
agrégée close et contient au moins quatre graphiques quantitatifs, cinq
équations, les effectifs, intervalles, anomalies et empreintes.

`calibrated=false` et `founder_effect_proven=false` restent invariants. D25a
n'évalue toujours ni l'erreur SHAPEIT5 ni la spécificité insuffisante de D24.
