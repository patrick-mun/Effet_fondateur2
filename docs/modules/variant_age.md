# Datation du variant — étape 14

## Méthode primaire

`gamma_gandolfo_2014_v1` réimplémente en Python l'estimateur décrit dans
Gandolfo, Bahlo et Speed, *Dating Rare Mutations from Small Samples with Dense
Marker Data*, Genetics 197:1315–1327 (2014), DOI
`10.1534/genetics.114.164616`.

Gamma est retenu parce qu'il utilise directement les longueurs génétiques des
segments ancestraux et a été conçu pour les petits échantillons typiques des
mutations rares. Il n'exige ni fréquence actuelle de la mutation, ni scénario
de croissance démographique. Le contexte insulaire motive le modèle corrélé en
primaire ; le modèle indépendant reste systématiquement visible.

L'implémentation reproduit l'exemple officiel à six individus : avec correction
de partage fortuit, elle retrouve 27,5 générations sous indépendance et 22,3
générations sous corrélation, ainsi que les intervalles publiés à l'arrondi.

L'audit D25a a toutefois montré que cet exemple positif ne couvre pas les
trois branches de corrélation négative du script R. La candidate séparée
`gamma_gandolfo_2014_r_faithful_v2` reproduit les estimations et intervalles R
sur cinq fixtures avec une erreur maximale de 7,11e-15 génération. Elle reste
expérimentale : le pipeline conserve V1 tant que le modèle d'observation des
extrémités tronquées n'est pas qualifié.

## Entrées

`variant_age_input.tsv` est dérivé de `founder_segments.tsv` et contient une
ligne par unité familiale indépendante. Seuls les segments
`target_centered_pairwise_max_ibs_v1` ayant deux bras strictement positifs sont
inclus. Les exclusions restent présentes et motivées.

Les positions absolues, fichiers MAP/RAW et effectifs d'individus non gelés ne
sont jamais transmis à Gamma. Les centiMorgans sont convertis en Morgans à la
frontière du calcul.

## Calculs

Pour chaque unité, la longueur totale est la somme des bras gauche et droit.
La correction de bord Gamma est réservée aux bras tronqués par le critère de
partage paire-à-paire : elle approxime la distance attendue jusqu'au premier
événement non observé. Elle ne doit pas être appliquée à des distances
biologiques complètes jusqu'au point de recombinaison. Le modèle indépendant
utilise l'estimateur sans biais et l'intervalle pivot exact de la loi Gamma. Le
modèle corrélé estime la corrélation et l'effectif effectif depuis la
distribution des longueurs.

D25a confirme que l'oracle indépendant sans correction est calibré lorsque les
points de recombinaison complets sont connus. En revanche, la correction
publiée appliquée aux longueurs tronquées surestime encore l'âge de 5,91
générations en moyenne et ne couvre que 91,88 % des âges vrais. Ces valeurs
interdisent d'interpréter la correction comme une modélisation validée de la
censure. Une future méthode devra intégrer explicitement le mécanisme
d'observation dans sa vraisemblance.

D25b implémente cette solution pour le seul domaine indépendant. Lorsque le
maximum vrai est remplacé par le deuxième maximum observable sur chaque côté,
la somme des longueurs est un temps total sous risque de censure de type II.
Avec `k=2(n−1)` ruptures, `gamma_pairwise_type2_censored_v3` utilise
`(k−1)/T` et l'intervalle Gamma exact. Le banc indépendant retrouve un biais de
−0,19 génération et une couverture de 95,11 %.

Cette candidate n'est pas une méthode corrélée. Les stress D25b montrent une
couverture de seulement 87,38 % et 76,83 % sous 25 % et 50 % de chemin
partagé. Elle reste donc non promue et ne peut être utilisée qu'en référence
indépendante explicitement étiquetée.

La correction de partage allélique fortuit est désactivée par défaut. Son
activation exige trois valeurs documentées pour le chromosome entier : fréquence
allélique médiane, nombre de marqueurs et longueur génétique en cM. Aucun défaut
silencieux n'est utilisé.

## Politique de conclusion

- moins de 3 unités incluses : `NOT_ESTIMATED` ;
- 3 ou 4 unités : `EXPLORATORY_ONLY` ;
- au moins 5 unités : `PRIMARY_ESTIMATE` ;
- le modèle corrélé V1 reste l'estimation principale historique, sans
  qualification scientifique ;
- le modèle indépendant et le leave-one-family-out sont des sensibilités ;
- 25, 28 et 30 ans par génération sont publiés séparément, jamais fusionnés en
  une date unique.

Une revue manuelle est toujours obligatoire. L'intervalle Gamma est conditionnel
aux limites observées : l'incertitude de phase, de carte et de résolution entre
le dernier marqueur partagé et le premier marqueur discordant reste une limite
scientifique explicite.

## Sorties

- `variant_age_input.tsv` : unités incluses et exclues ;
- `variant_age_estimates.tsv` : Gamma corrélé et indépendant ;
- `variant_age_scenarios.tsv` : modèles, leave-one-family-out et conversions en
  années ;
- `variant_age_summary.json` : résultat agrégé sans identifiant individuel ;
- `audit.json` : paramètres, versions, contrôles et limites.

EstiAge n'est pas exécuté à cette étape. Il pourra devenir une validation
secondaire lorsqu'une implémentation locale, versionnée et sans transfert de
données privées aura été qualifiée.
