# Configurations de calibration

Ce dossier recevra uniquement les configurations reproductibles et
versionnées des micro-panels, pilotes et panels de qualification.

Les configurations ne devront contenir ni chemin absolu local, ni donnée
génétique, ni identifiant réel. Elles ne deviendront exécutables qu'après
validation de leur schéma et des gardes d'isolation.

`d12_diagnostic_100kb_calibration_v1.json` fixe uniquement l'intégration
technique D12 liée au domaine D11 de 100 kb. Elle crée 20 diploïdes réciproques
CEU/YRI à 50/50 avec la graine haplotypique D05 `1002`. Cette configuration ne
vaut ni micro-panel D33, ni pilote, ni qualification scientifique.

`d13_diagnostic_100kb_calibration_v1.json` ajoute une famille CEU et une
famille `ADMIX_CEU_YRI` au panel D12 borné, avec la graine `pedigrees` `1003`.
Chaque unité contient quatre fondateurs, quatre descendants et quatre trios ;
deux paires indépendantes explicites complètent le contrôle négatif. Cette
configuration diagnostique ne vaut ni pilote, ni qualification scientifique.

`d14_diagnostic_100kb_age150_calibration_v1.json` utilise la graine de
recombinaison D05 `1004`, un âge vrai de 150 générations et le SNV synthétique
`chr19:11050057 C/T`. L'âge appartient à la grille scientifique minimale, mais
la fenêtre de 100 kb est trop courte pour attendre des breakpoints internes ;
elle contrôle principalement la censure et la propagation mendélienne.

`d14b_micro_chr5_age100_calibration_v1.json` fixe, hors chromosome 19, une
cible publique synthétique sur le chromosome 5, une fenêtre principale à 100
générations avec `p = 0,05` par bras et un défi de censure de `0,1 cM` par
côté. Il préspecifie les cinq longueurs IBD de 1 à 5 cM et leurs seuils de
marqueurs avant toute simulation.

`d14_micro_chr5_age100_calibration_v1.json` rejoue le moteur D14 avec la même
cible, l'âge de 100 générations et la graine `1004` sur la fenêtre principale
D14b. Cette configuration est un contrôle technique pré-micro-panel ; elle ne
vaut ni pilote ni qualification scientifique.

`d16_micro_chr5_calibration_v1.json` fixe quatre scénarios ROH : un témoin sans
injection, une ROH de 2 Mb hors cible, une ROH de 2 Mb sur la cible sans IBD
injecté et une ROH de 2 Mb couplée au positif D15 `primary_2cm_positive`.
La copie porte uniquement sur H2 ; H1, les négatifs IBD et les pedigrees sont
protégés. Le cadre observable exige 1 500 kb, 50 sondes, 50 kb/sonde au
maximum et un trou interne inférieur ou égal à 1 000 kb. Les marqueurs
flanquants sont exigés sans déplacer les bornes en cas de manque.

`d17_micro_chr5_calibration_v1.json` utilise la graine de mesure `1005` et
huit profils appariés : contrôle, manque 1 %/5 %, dosage 1 %, points de phase
2 %, frontières de phase 2 %, rétention de sondes 50 % et défauts combinés.
Les quotas et dénominateurs sont définis dans le plan D17. Les taux sont des
stress techniques présélectionnés, pas une mesure de l'erreur empirique de la
puce. Les graines de validation devront être distinctes ; aucune cible
absente de la puce n'est ajoutée.

`d18_micro_chr5_calibration_v1.json` fixe la graine de blinding `1006`, le rôle
de calibration, la cible d'analyse configurée et la politique de métadonnées
neutres. La portée est exclusivement régionale et la configuration pipeline
livrée est limitée à 00–02, sans approbation humaine automatique. Les profils
et sujets ne sont reconnaissables que via la correspondance privée.

`d19_technical_target_v1.json` décrit un scénario de concordance 03–04,
distinct des injections biologiques D14–D18. L'annotation publique APC
rs1801155 est épinglée à ClinVar RCV001195214.19 ; les génotypes ne proviennent
d'aucun patient et suivent une règle RR/RA/AA selon le rang aveugle modulo 3.
La décision `CAL-DEC-003` autorise seulement le registre technique neutre et
l'intégration régionale 00–04 avec PLINK. Aucun nouveau seuil ni graine ;
la distribution n'est pas une fréquence populationnelle. Les anciens fichiers
D18 et leurs configurations 00–02 restent immuables.
