# ROH secondaire — contrat scientifique V2

## Question scientifique

L'étape `16_analyze_roh` mesure l'autozygotie individuelle. Un run of
homozygosity est compatible avec la transmission de deux copies d'un segment
ancestral, mais ne démontre pas à lui seul que plusieurs familles partagent le
même haplotype fondateur. La méthode versionnée est
`plink19_array_roh_secondary_v1`.

## Deux périmètres séparés

`GENOMEWIDE_BURDEN` utilise uniquement `controls_unrelated_qc` et
`target_carriers_independent_qc`. L'univers de variants est leur intersection
exacte après vérification chromosome, position et allèles. Les scans restent
séparés par cohorte, avec le même univers et les mêmes paramètres.

`TARGET_CHROMOSOME` utilise `target_chromosome_all_qc`. Il conserve tous les
individus QC, leurs appartenances de cohortes et leur génotype cible accepté.
Pour chacun, l'appartenance du variant à un ROH est définie uniquement par
`POS1 <= TARGET_BP <= POS2` sur le chromosome cible.

Les non-porteurs familiaux ne participent pas au fardeau genome-wide, car aucun
jeu autosomique final comparable ne leur est actuellement publié. Ils restent
descriptibles sur le chromosome cible. Une comparaison de groupe ne peut
utiliser que des unités indépendantes.

## Profil PLINK exploratoire

Le profil par défaut fixe explicitement : 1 500 kb, 50 variants, 50 kb par
variant au maximum, un trou interne maximal de 1 000 kb, une fenêtre de 50
variants, un hétérozygote et cinq manquants au maximum par fenêtre, et un taux
minimal de fenêtres concordantes de 0,05. Ces paramètres sont adaptés à une
première exploration sur puce, mais leur stabilité doit être étudiée à l'étape
17 avec plusieurs longueurs et nombres de marqueurs.

La densité ou la couverture insuffisante ne déclenche jamais un assouplissement
automatique. Le périmètre devient `NOT_EVALUATED`. Les sorties natives PLINK et
les paramètres exacts restent liés à l'audit.

## Mesures

Chaque individu possède une ligne de fardeau même avec zéro segment : nombre de
ROH, longueur totale, moyenne, maximum et sommes par classes de longueur. Un
`F_ROH` n'est publié que si un dénominateur autosomique en kb et sa provenance
sont configurés ; aucune taille de génome implicite n'est utilisée.

`target_in_roh.tsv` contient également les individus sans ROH couvrant la cible.
Un génotype cible hétérozygote inclus dans un segment tolérant un hétérozygote
est signalé pour revue et n'est pas présenté comme une preuve d'autozygotie au
locus.

## Limites et articulation

La densité de la puce, les erreurs de génotypage, les trous entre marqueurs, la
consanguinité récente et la démographie modifient les ROH observés. Homozygotie,
IBS et IBD ne sont pas synonymes.

L'étape 16 ne redéfinit pas le partage IBS de l'étape 13, ne fournit aucune
longueur à Gamma en 14 et ne consomme pas le LD de 15. Une analyse conjointe ou
un changement de profil est un scénario de sensibilité de l'étape 17.
