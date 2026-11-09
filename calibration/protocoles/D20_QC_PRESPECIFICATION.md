# D20 — Banc de conformité QC 05 / 10

Version 1.0.1, 30 octobre 2026. Méthode `iid_autosomal_qc_component_bench_v1`.
À lire avant toute exécution ; aucun seuil ne sera réglé sur ses résultats.

La version initiale 1.0.0 a précédé la construction et les runs. La révision
1.0.1 consigne un renforcement du comparateur lors de la relecture : tolérance
absolue nulle pour HWE_P. Les premières comparaisons PASS ne suffisaient pas
à protéger les très petites p-values d'une assimilation à zéro. Un test de
régression et une nouvelle évaluation complète appliquent la règle plus
stricte ; aucun seuil scientifique, génotype ou défaut injecté n'est changé.

Ce contrôle renforcé a effectivement révélé un défaut de sérialisation en
10 : `_decimal` écrasait les valeurs de module inférieur à 5e-16, notamment
une p-value HWE brute PLINK de 2,373e-28. L'énumération entière indépendante
donne 2,373074178759543e-28. Le test de régression échoue avant correction ;
la correction normalise seulement le vrai zéro signé. Les décisions HWE,
leur seuil et les autres méthodes restent inchangés. Le premier run est
conservé en FAIL sur la métrique numérique ; un nouveau run vérifie le code
corrigé sans réécrire les résultats historiques.

## Préflight et portée

D18 est régional, mono-lot, sans cohorte figée issue de 09. Ses erreurs de
phase ou de dosage ne sont pas toutes identifiables par QC. Il ne permet pas
à lui seul de mesurer les règles de 05 et 10. D20 ajoute donc un **banc de
composants**, distinct du panel biologique 1000G : 192 sujets artificiels,
22 autosomes, 64 marqueurs indépendants par autosome, deux lots entrelacés,
graine explicite et génotypes diploïdes Bernoulli de fréquence 0,35.
Les coordonnées artificielles ne décrivent ni LD, ni démographie, ni IBD.
L'en-tête ACPA sert à tester le parseur : ce n'est pas une simulation de la
densité ou des propriétés de la puce réelle.

Les quatorze cas partagent le même fond : propre, manque individuel, manque
variant, hétérozygotie extrême, effet de lot, clone technique, allèle rare,
HWE témoins, exception cible rare, manque individuel régional et manque
variant limité aux témoins, manque variant limité aux porteurs, manque variant
régional et manque individuel devenu filtrant après 05. Ce sont des stress
appariés, pas quatorze réplications
indépendantes. Aucun intervalle populationnel n'est calculé sur ces cas.

## Injection préspecifiée

- Manque individuel : 30 % des marqueurs du premier sujet.
- Manque variant : premier variant absent chez le premier quart des sujets.
- Hétérozygotie : premier sujet hétérozygote à tous les marqueurs.
- Lot : cinq sujets du lot B absents sur quatre marqueurs sur cinq, rotation
  déterministe ; un marqueur sur cinq reste un contrôle négatif apparié ;
  manque variant global < 5 %, mais écart entre lots > 2 %.
- Clone : dernier sujet copie l'avant-dernier ; un seul appartient à la
  fixture de porteurs indépendants, les deux au jeu régional.
- Rare : un seul hétérozygote aux variants d'indices 12 et 268 (base zéro,
  respectivement chr1 et chr5).
- HWE : tous les témoins hétérozygotes au variant d'indice 13.
- Cible rare : un seul porteur hétérozygote de la cible technique publique
  de D19 ; les appels restent fictifs, sans inférence clinique.
- Manque individuel régional : premier sujet absent sur 30 % des marqueurs
  de chr5, sous le seuil global mais au-dessus du seuil local de 10.
- Manque variant témoins : variant d'indice 14 absent chez huit témoins,
  sous le seuil global de 05 et au-dessus du seuil témoins de 10.
- Manque variant porteurs : même principe, indice 15 chez les huit premiers
  porteurs ; contrôle positif du filtre de manque, mais pas de MAF/HWE.
- Manque variant régional : premier marqueur de chr5 absent chez le premier
  quart des sujets ; retiré en 05 mais encore présent à l'entrée locale de 10.
- Manque individuel après 05 : 200 premiers variants absents chez les dix
  premiers sujets ; 130 derniers variants absents chez les sujets de rangs
  20 et 100. Ceux-ci passent le manque global initial de 05 mais dépassent
  10 % après retrait des 200 variants, dans chacune des deux cohortes globales.

La cible ajoutée par 04 est normalement RR dans la première moitié et RA
dans la seconde. Les lots A/B alternent à l'intérieur de chaque moitié.
Pour le cas cible rare, seul le dernier sujet est RA ; les cohortes sont
définies à partir de ces appels explicites, jamais du statut clinique.
Sa cohorte de porteurs est donc un singleton : son résumé de lot teste le
comportement contractuel mono-lot, pas la détection d'un écart entre deux
lots. Les autres profils assurent les positifs et négatifs entre lots. Les
matrices publiées mesurent la conformité au contrat, pas une performance
populationnelle de détection sur cette cohorte singleton.

## Exécution et indépendance

Le vrai pipeline exécute 00–05, avec PLINK réel annoncé séparément. 10 est
appelée par son point d'entrée public avec des contrats de cohorte
**préspecifiés synthétiques**, producteur `D20_COHORT_CONTRACT_FIXTURE`.
Ce n'est pas une exécution de 06–09 : aucune validation de parenté,
structure ou sélection de cohortes n'est revendiquée. Les cohortes sont
préparées avant le run selon l'oracle d'entrée de 05, jamais ajustées aux
sorties observées. 09 reste à tester en D22. Le manifest d'orchestration ne
doit pas contenir de faux succès 09 ou 10.

Observations ACPA/VCF et métadonnées neutres : `blinded_inputs/`. Injection,
cohortes attendues et configuration source : `private_truth/`. Le processus
pipeline ne reçoit aucun chemin de vérité ; la fixture contractuelle de 10
est une entrée publique séparée et ne contient aucun motif de défaut.
Les empreintes des entrées et sorties, versions, codes de retour et reçus
CLOSED doivent être contrôlés avant l'ouverture par l'évaluateur.

## Oracle et réussite technique

L'oracle NumPy ne réutilise pas les fonctions QC de production. Il calcule
les manques entiers, fréquences, hétérozygotie, z-score de population et
écarts de lots. HWE est énuméré par distribution conditionnelle exacte,
sans mid-p. Références :
[Wigginton et al., 2005](https://pmc.ncbi.nlm.nih.gov/articles/PMC1199378/),
[PLINK statistiques](https://www.cog-genomics.org/plink/1.9/basic_stats),
[filtres](https://www.cog-genomics.org/plink/1.9/filter#hwe),
[PI_HAT](https://www.cog-genomics.org/plink/1.9/ibd).

05 exclut simultanément sujets et variants sur leurs manques initiaux.
MAF, hétérozygotie, lots et doublons sont des alertes, pas des exclusions.
10 retire d'abord les sujets, puis recalcule les métriques de variants.
HWE et MAF témoins ne s'appliquent pas aux porteurs. La cible régionale est
conservée par exception documentée si nécessaire. Les décisions et motifs
attendus doivent concorder exactement ; les mesures PLINK arrondies sont
comparées avec une tolérance déclarée de 5e-4 relative et 5e-6 absolue.
Pour HWE_P, la tolérance absolue est nulle : une très petite p-value ne doit
jamais être assimilée à zéro ou à une valeur franchissant le seuil HWE.
Les axes et tous les génotypes retenus sont comparés exactement.

Chaque famille de règle reçoit une matrice de confusion sur sa propre unité
(sujet, variant, lot ou paire). Le clone est une vérité de détection, pas un
oracle indépendant de l'estimateur PI_HAT : les autres alertes sont des faux
positifs descriptifs, non des erreurs d'implémentation présumées.
Les matrices de conformité à l'oracle et celles des défauts ciblés sont
distinguées. Dénominateur nul => taux null et motif explicite, jamais PASS
par défaut. La couverture agrégée exige au moins un positif et un négatif
pour chaque règle destinée à détecter un défaut.

Les contrôles sexuels ne sont pas évaluables avec un contrat autosomal :
statut explicite `NOT_EVALUATED_BY_AUTOSOMAL_CONTRACT`, sans invalider la
mesure des règles autosomiques. HWE en 05 et HWE porteurs en 10 sont
interdits par politique, pas des données manquantes.

Un PASS signifie conformité fonctionnelle dans ce domaine seulement.
Les seuils de performance du pilote ne sont pas scellés : décision
scientifique `REVIEW_REQUIRED`, `calibrated=false`,
`founder_effect_proven=false`. Rapports indépendants CAL-REP-05 et CAL-REP-10,
complétés par un rapport de développement et leurs preuves de rendu.
