# Enrichissement du partage haplotypique fondateur — étape 16B

## Question et unité primaire

L'étape `16B_evaluate_founder_haplotype_enrichment` teste si le partage IBS
exact centré sur la cible entre familles porteuses indépendantes est plus long
que le partage fortuit au même locus. L'unité primaire est la famille gelée à
l'étape 09. Les individus ou copies mutantes supplémentaires ne sont jamais
traités comme des réplications indépendantes.

La sélection des représentants est reprise sans modification de l'étape 13.
Il est interdit de rechercher après coup la copie familiale donnant le segment
le plus long. La méthode `target_centered_exact_ibs_v1` est reproduite avant le
test ; toute discordance avec les bras et comptes de marqueurs publiés par 13
bloque 16B.

## Statistique et fonds nuls

La statistique primaire préspécifiée est
`T_TOTAL_CM = LEFT_SHARED_CM + RIGHT_SHARED_CM`. La cible sert d'ancre mais est
exclue de la signature. Chaque bras s'arrête au premier manque ou désaccord.

Le fond interne énumère exhaustivement les combinaisons d'haplotypes non
porteurs fiables. Les unités d'un tirage appartiennent à des individus
distincts ; `H1` et `H2` d'un même individu ne peuvent jamais coexister.

Le fond externe tire sans remise des individus parmi les 2 504 références
1000 Genomes non apparentées confirmées par 16A, puis une seule copie par
individu. Il vise par défaut 100 000 tirages évaluables avec une graine fixe,
conserve les tentatives non évaluables et n'effectue aucun arrêt anticipé en
fonction du résultat. Les superpopulations sont des sensibilités descriptives,
jamais une attribution d'ascendance aux familles de l'étude.

La probabilité Monte-Carlo est `(1 + k) / (1 + N)`, avec `k` tirages tels que
`T_null >= T_observed`. Un intervalle binomial de Wilson documente l'incertitude
Monte-Carlo. L'énumération interne publie en plus la fraction exacte.

## Interprétation et séparation des domaines

Sans seuil préspécifié, le statut est `NOT_CLASSIFIED`. Les statuts autorisés
restent `NOT_EVALUATED`, `NOT_CLASSIFIED`,
`NO_UNUSUAL_SHARING_DETECTED`, `UNUSUAL_TARGET_CENTERED_SHARING` et
`MULTIPLE_CARRIER_BACKGROUNDS`.

L'étape 15 décrit séparément le LD populationnel de fond (`r²` et `D′`). Elle
n'alimente aucun calcul 16B. Le ROH, l'ascendance, la datation et l'IBS ne sont
ni additionnés ni transformés en score composite.

Même un partage inhabituel est seulement compatible avec un haplotype
ancestral commun. Il ne prouve ni IBD, ni effet fondateur, ni ancêtre unique,
ni origine géographique ou ethnique.

## Sorties et sécurité

Les unités, bornes et cohérences familiales restent `sensitive_genetic`. Les
tirages nuls gzip ne contiennent aucun identifiant. Le résumé JSON contient
uniquement des effectifs, longueurs, probabilités, intervalles, statuts et
empreintes de provenance. Le calcul utilise les BCF/VCF locaux et n'envoie
aucune donnée de l'étude sur le réseau.

La figure 18 présente les fonctions de survie interne et externe, la ligne de
`T_observed`, les deux bras et la mention « partage IBS centré cible, pas preuve
IBD ». Le rapport 19 consomme uniquement le résumé agrégé et la provenance de
la figure.
