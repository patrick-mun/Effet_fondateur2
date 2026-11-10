# D26 — Préspecification fermée des évaluateurs LD et ROH

Statut : `FINALIZED_BEFORE_RUN`
Autorisation : `CAL-DEC-014`
Rapports : `CAL-REP-15`, `CAL-REP-16`

## 1. Question et portée

D26 mesure la fidélité opérationnelle des étapes réelles `15_analyze_local_ld`
et `16_analyze_roh` sur le panel synthétique chr5 créé en D16 puis perturbé
en D17. D26 ne teste ni DOCK6, ni le chromosome 19, ni les données d'étude.
Une réussite technique ne vaut ni calibration scientifique, ni preuve d'effet
fondateur.

Les six observations non redondantes sont `control`, `missing_01`,
`missing_05`, `genotype_01`, `density_50` et `combined`. Les observations
non phasées de `phase_point_02` et `phase_switch_02` doivent être identiques au
contrôle ; elles sont vérifiées par SHA-256 et ne constituent pas deux
réplications supplémentaires.

## 2. Séparation des appels et de la vérité

La phase d'appel lit uniquement les VCF observés D17, la carte génétique D11a
et les paramètres finalisés ci-dessous. Elle exécute PLINK 1.9 par les mêmes
moteurs de publication que les étapes 15 et 16, puis scelle les tables et
leurs empreintes. Elle ne lit ni le VCF dense D16, ni `roh_truth.json`.

L'évaluation ouvre ensuite :

1. le VCF de contrôle D17 pour calculer directement les corrélations de dosage
   servant de référence LD sur les mêmes variants ;
2. le VCF phasé dense D16, converti en génotypes, pour produire une référence
   d'appels ROH dense avec les mêmes paramètres PLINK ;
3. la vérité d'injection D16, uniquement pour distinguer les trois contrôles
   positifs injectés du contrôle sans injection.

La référence dense ROH est une référence **opérationnelle** dépendante de
PLINK, pas une vérité biologique indépendante. Elle permet de mesurer la perte
due à l'observation sur puce, mais pas de valider la définition scientifique
d'un ROH. Cette limite interdit une qualification automatique.

## 3. LD : unités, équation et métriques

Le LD primaire est le carré de la corrélation de Pearson entre dosages
alléliques diploïdes, sur les individus appelés aux deux loci :

\[
r^2_{jk}=\left(\frac{\sum_i(g_{ij}-\bar g_j)(g_{ik}-\bar g_k)}
{\sqrt{\sum_i(g_{ij}-\bar g_j)^2\sum_i(g_{ik}-\bar g_k)^2}}\right)^2.
\]

La référence est calculée indépendamment du parseur PLINK sur le contrôle
intact. Les paires sont limitées à 1 Mb et 1 cM, avec MAF minimale 0,05 et au
moins cinq appels conjoints. Pour un profil `p`, l'erreur absolue moyenne est

\[
\operatorname{MAE}_{LD,p}=|E_p|^{-1}\sum_{(j,k)\in E_p}
|\hat r^2_{p,jk}-r^2_{0,jk}|.
\]

Sont aussi publiés le nombre de paires évaluables, les quantiles d'erreur et
la dégradation relative au contrôle. `CAL-LD-001` reste
`REVIEW_REQUIRED` : son seuil `pilot_then_seal` était nul avant D26 et ne sera
pas inventé après observation. `CAL-LD-002` exige zéro chromosome cible dans
un fond génomique ; le banc régional ne crée aucun fond génomique.

## 4. ROH : appariement et métriques

Les paramètres de production sont fixés à 1 500 kb, 50 SNP, 50 kb/SNP,
1 000 kb de trou maximal, fenêtres de 50 SNP, au plus un hétérozygote et cinq
manquants, seuil de fenêtre 0,05. Un appel observé et une référence dense du
même sujet sont appariés si leur intersection couvre au moins 50 % de la
référence. Chaque segment n'est apparié qu'une fois, par intersection maximale.

La sensibilité est `TP/(TP+FN)`. Le taux de faux positifs est défini au niveau
des sujets négatifs de la référence dense : nombre de sujets sans ROH dense
mais avec au moins un appel observé, divisé par le nombre de sujets sans ROH
dense. La proportion de segments non appariés est publiée séparément et ne
doit pas être appelée spécificité.

Pour les appariements, la longueur est convertie en cM par interpolation de la
carte D11a et l'erreur absolue moyenne vaut :

\[
\operatorname{MAE}_{L,p}=M_p^{-1}\sum_{m=1}^{M_p}
|L^{obs}_{m,p}-L^{dense}_{m}|.
\]

Les trois injections D16 forment un sous-ensemble descriptif. Les bornes
injectées ne sont pas traitées comme limites maximales d'appel. Le témoin
`NO_INJECTION` n'est pas supposé négatif sans contrôle de la référence dense.
`CAL-ROH-001` à `003` restent `REVIEW_REQUIRED`, leurs seuils n'étant pas
scellés. `CAL-ROH-004` exige qu'aucune longueur ROH ne soit utilisée comme
entrée directe de Gamma.

## 5. Verdicts autorisés

Le statut technique est `PASS` seulement si les six profils sont exécutés,
les sources et sorties sont liées par SHA-256, les duplicats de phase sont
confirmés, les tables sont lisibles et aucun usage Gamma ou donnée réelle
n'est observé. Le statut de performance reste `REVIEW_REQUIRED_UNSEALED` quel
que soit le résultat numérique. `calibrated` et `founder_effect_proven`
restent obligatoirement `false`.
