# Appel IBD explicite sur puce SNP — étape 16C

## Portée scientifique

`16C_call_explicit_ibd` recherche des segments contenant le variant cible avec
deux méthodes indépendantes, Hap-IBD et Refined IBD. Elle ne remplace ni l'IBS
strict de 13, ni la datation de 14, ni l'enrichissement empirique de 16B. Un
appel primaire concordant apporte un support IBD explicite ; il ne prouve pas à
lui seul un effet fondateur, un ancêtre unique ou une origine géographique.

Une absence d'appel n'est jamais une réfutation. Sur le run DOCK6 de référence,
l'intersection IBS mesure environ 1,436 cM et 18 marqueurs, sous le scénario
primaire à haute spécificité de 2 cM et 100 marqueurs.

## Préparation commune

Les deux outils consomment exactement le même VCF phasé et la même carte
génétique. Chaque variant est joint par identifiant, chromosome et position à
la carte validée. Un variant comportant un GT absent, partiellement absent ou
non phasé est retiré des deux entrées avec le motif `TRUE_MISSING_CALL`. Aucun
génotype n'est complété ou requalifié en observation.

Le premier jalon sait consommer le BCF produit par 12 et publie un rapport de
faisabilité. Pour l'analyse scientifique prévue, `input_scope` doit être
`TARGET_CHROMOSOME` et le producteur de phasage doit fournir le chromosome
cible complet. Un panel régional peut servir aux tests techniques, mais sa
densité insuffisante doit conduire à `NOT_EVALUABLE` ou à une sensibilité, pas
à `NO_PRIMARY_IBD_CALL` interprété biologiquement.

## Seuils et calibration

Le contrat refuse un scénario primaire inférieur à 2 cM ou 100 marqueurs. Les
sensibilités ne peuvent pas descendre sous 1 cM. Chaque scénario enregistre son
rôle, ses seuils, son statut de calibration, l'acceptabilité de sa spécificité
et le seuil maximal de fréquence de fond. Ces valeurs doivent être figées dans
une nouvelle configuration avant le run réel. Une calibration absente ou non
acceptable produit `NOT_EVALUABLE`.

La configuration d'exemple garde toutes les calibrations à `false`. Leur
passage à `true` exige une simulation séparée reproduisant la densité réelle,
les positions cM, les taux de données manquantes et les incertitudes de phase,
avec graines et résultats audités. Il est interdit de modifier ces seuils parce
qu'ils rendent DOCK6 positif.

## Attribution mutante et consensus

Les copies mutantes proviennent exclusivement de `carrier_haplotypes.tsv`,
lui-même confronté aux génotypes moléculaires explicites. `H1`, `H2` et `BOTH`
ne sont jamais déduits du phénotype ou du groupe. Pour les trois familles, 16C
énumère les assignations possibles et exige une assignation unique cohérente
sur toutes les paires et les deux outils. Une arête F1–F2 sur H1 et une arête
F1–F3 uniquement sur H2 ne suffisent donc pas.

Le support primaire exige les deux outils, toutes les paires familiales, la
cible dans chaque segment et dans l'intersection commune, des limites dans la
tolérance préspécifiée, une calibration acceptable et une fréquence de fond
sous le maximum configuré.

## Statuts

- `PRIMARY_CONCORDANT_IBD_SUPPORT` : tous les critères primaires sont remplis ;
- `SENSITIVITY_ONLY_IBD_SUPPORT` : seul un scénario secondaire calibré est
  concordant ;
- `METHOD_DISCORDANT` : appels présents mais limites non concordantes ;
- `NO_PRIMARY_IBD_CALL` : aucun appel primaire conforme dans une analyse
  évaluable ;
- `NOT_EVALUABLE` : densité, carte, phase, calibration, témoins ou outils ne
  permettent pas l'évaluation.

Le résumé conserve séparément `ibs_only_from_step13`, `ibd_proven`,
`founder_effect_proven: false`, `geographic_origin_inferred: false` et
`composite_score_calculated: false`.

## Outils et acquisition

Le bloc `tools.explicit_ibd_adapters` exige la commande Java, sa version majeure
attendue, les chemins des deux JAR, leurs versions déclarées et leurs SHA-256.
Les chemins sont refusés s'ils sont absents, symboliques ou si l'empreinte
diffère. Les appels sont bornés en mémoire, threads et temps ; stdout/stderr est
conservé par scénario. Aucun téléchargement ou installation n'est effectué par
le pipeline.

## Artefacts

16C publie le rapport de faisabilité, l'audit du même univers de marqueurs, les
segments normalisés, la matrice par paire, la concordance, la calibration, la
fréquence de fond et le résumé scientifique. Les tables contenant individus,
familles ou haplotypes sont `sensitive_genetic`. Les contrats sont versionnés
dans `schemas/explicit_ibd_*.schema.json`.

La table de segments et le résumé sont prêts pour une figure montrant
séparément Hap-IBD, Refined IBD, la cible et leur intersection. Le raccord à la
galerie consolidée 18 doit être activé en même temps que le panel chromosome
19 complet afin de ne pas rendre obligatoire un résultat 16C dans les runs
historiques.

## Tests et limites actuelles

Les tests ordinaires simulent les sorties des deux JAR et couvrent cible
incluse/absente, limites discordantes, paire ou outil manquant, sensibilité
seule, GT manquant/non phasé, chromosome mutant non assignable, fond trop
fréquent, calibration absente, seuil primaire post hoc, échec et timeout.

Java et les JAR ne sont pas des dépendances Python et ne figurent pas dans
`requirements.txt`. Un smoke test avec les vrais outils reste obligatoire
avant tout run scientifique. Le run de référence du 13 août 2026 est immuable
et ne doit jamais recevoir manuellement des artefacts 16C.
