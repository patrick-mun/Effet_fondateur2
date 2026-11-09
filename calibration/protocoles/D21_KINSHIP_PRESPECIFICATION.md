# D21 — Parenté et proposition indépendante : préspecification

Version 1.0.0, 30 octobre 2026. Méthode `iid_pedigree_kinship_bench_v1`.
Ce document précède la génération et les runs. Aucun seuil ne sera modifié
après lecture des résultats pour faire réussir le banc.

## Prérequis vérifiés et manques résolus

D13–D18 fournissent des pedigrees régionaux mais pas un panel autosomique
complet assez informatif pour tous les degrés. D20 couvre les autosomes mais
pas les familles, strates et degrés nécessaires. KING 2.3.2 et PLINK 1.9
sont disponibles ; leur usage réel est réservé à l'intégration annoncée.

D21 ajoute un banc mathématique indépendant : 22 autosomes, 512 loci par
autosome, deux strates artificielles de fréquences 0,25 et 0,65, huit unités
indépendantes par strate et type de relation. Les sept types sont : clone,
parent-enfant, germains, demi-germains, cousins germains, cousins germains
avec une génération de décalage et sujets non apparentés. Soit 112 unités
et 224 sujets observés. Les ancêtres non génotypés restent privés.

Les fondateurs sont tirés indépendamment, puis une copie de chaque parent
est transmise par locus. C'est un gene-drop mendélien IID, sans LD biologique
ni longueur de segment : il ne qualifie ni phasage, ni IBD, ni démographie,
ni la puce réelle. Les coefficients théoriques sont calculés par récurrence
du pedigree ; l'origine des copies permet un contrôle indépendant de l'IBD
réalisé dans ce modèle, sans l'assimiler à un appel IBD du pipeline.

## Couverture et entrées publiques

Trois profils appariés : complet, manque indépendant de 5 %, densité 25 %.
Ils partagent les mêmes unités ; ils ne sont pas trois réplications.
Les pseudonymes sont distincts par profil. Les strates, types et vérités de
parenté restent privés. Le registre public comporte seulement les liens
parent-enfant réellement observables et les identifiants techniques de groupe.
Un négatif partage volontairement un FID technique pour exiger une mesure
KING explicite de non-apparentement. Le premier clone de chaque strate est
réparti sur deux FID pour contrôler la voie inter-familles `.kin0`.

Une copie parfaite des deux premiers loci de chaque autosome éprouve le
pruning LD ; un locus monomorphe et une région complexe artificielle éprouvent
les autres motifs de retrait. Ces interventions sont propagées aux fondateurs
et transmissions, pas ajoutées comme erreurs mendéliennes aux descendants.
Les coordonnées, graines et seuils sont dans la configuration versionnée.

Le vrai pipeline est borné à 00, 01, 02, 03, 05, 06 et 07. Il reçoit seulement
`blinded_inputs/`. Il ne lit aucune vérité et n'active ni 04 ni 08–19. Le
runner laisse `kinship_exclusion_approval` en attente : 07 propose un ensemble,
elle n'applique pas les exclusions. Aucune approbation de cette proposition,
aucun gel de cohorte 09 ni aucune qualification n'est fabriqué.

## Oracle 06 : sélectionner et conserver les bons génotypes

Recalcul indépendant de la MAF sur les fondateurs PLINK déclarés, après 05 ;
vérification des motifs MAF/région, de la partition éligibles/retenus/LD,
du retrait d'au moins une copie des doublons LD, des 22 couvertures, des
distances et des comptes. Aucun tie-break particulier du pruning n'est imposé
si plusieurs panels satisfont la règle. Le jeu retenu doit conserver exactement
les individus, l'ordre et tous les génotypes. Le LD résiduel est recalculé
sur les mêmes fenêtres et fondateurs ; les classes et mesures sont confrontées
aux agrégats publiés, sans considérer une faible LD comme preuve d'IBD.

## Oracle 07 et métriques mesurables

Les métriques KING sont recalculées depuis les génotypes, indépendamment des
helpers de production : équation 9 en intra-groupe, équation 11 inter-groupes
de [Manichaikul et al., 2010](https://www.cog-genomics.org/static/pdf/Manichaikuletal2010.pdf).
Les comptes appariés excluent les manques dans l'un ou l'autre sujet. Le
précédent méthodologique n'est pas une garantie biologique pour ce banc.

Les degrés publiés doivent respecter les seuils configurés, sur les mesures
KING effectivement sérialisées. Leur concordance fonctionnelle est distincte
de l'exactitude par rapport au pedigree simulé. Les différences numériques
doivent respecter la précision native KING préspecifiée (5,1e-5 absolue).

La matrice de confusion comprend les six degrés et une colonne MISSING.
Toutes les paires cibles intra-FID, dont les négatifs, doivent avoir une mesure.
Une paire positive inter-FID omise compte comme manquée, pas UNRELATED.
Les négatifs inter-FID non publiés ne sont pas inventés comme mesures.
Les résultats sont stratifiés par strate, profil de manque, densité et degré.
Le bootstrap rééchantillonne les unités familiales, jamais les paires comme
si elles étaient indépendantes. Les graines, répétitions, quantiles et faibles
effectifs sont enregistrés ; une borne ponctuelle à 100 % ne remplace pas la
qualification sur au moins 200 paires indépendantes par classe.

L'oracle reconstruit le graphe des relations observées et déclarées, puis
vérifie la proposition gloutonne orientée qualité : choix, représentants,
motifs, indépendance et maximalité (pas cardinalité maximum). La proposition
est aussi comparée au graphe de vérité : toute paire théoriquement apparentée
conservée ensemble est comptée, même si KING l'a omise. Ce contrôle porte sur
la proposition 07, pas sur une cohorte gelée 09 encore absente.

## Décision et livrables

Rapports autonomes CAL-REP-06 et CAL-REP-07, plus CAL-DEV-D21, HTML/PDF selon
report-formater. Les entrées, sorties, versions et clôtures sont empreintées
avant ouverture de la vérité. Les tests unitaires simulent les exécutables.

Un défaut de contrat, un génotype modifié ou une paire interdite conservée
entraîne FAIL ; une erreur de classification reste visible dans les métriques.
Le seuil de performance du pilote n'est pas scellé : REVIEW_REQUIRED.
`calibrated=false`, `founder_effect_proven=false`. Les strates sont artificielles
et ne sont jamais présentées comme des ascendances humaines validées.

## Complément de lecture après le premier run — 30 octobre 2026

Sans changement des critères ou du panel : les rappels par degré reçoivent
un intervalle de Wilson descriptif ; le bootstrap est stratifié par degré et
strate. Ces intervalles sont conditionnels au panel commun sélectionné, dont
le pruning n'est pas réexécuté dans chaque tirage. Ils ne couvrent donc pas
l'incertitude de cette sélection. La BA avant pruning, recalculée par l'oracle
sur la sortie 05, et les moyennes φ théorique/réalisé/estimé sont ajoutées comme
diagnostics exploratoires, pas comme résultats préspecifiés ou secondes
exécutions KING. Les faux liens rapportés entre unités indépendantes sont
comptés sans inventer de taux sur les paires inter-FID censurées.
