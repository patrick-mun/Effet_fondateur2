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

`d20_qc_bench_v1.json` configure un banc mathématique IID de composants QC,
distinct du panel 1000G : graine 2026103020, 192 sujets, 64 loci par autosome,
fréquence initiale 0,35 et deux lots. Quatorze cas et leurs injections fixes
sont spécifiés dans `protocoles/D20_QC_PRESPECIFICATION.md` ; leurs constantes
appartiennent à la méthode versionnée, pas à un ajustement après résultats.
Les seuils 05/10 sont explicites ; la configuration cible D19 est épinglée par
SHA-256, mais la règle de génotypes et l'autorisation sont celles de D20.
CAL-DEC-004 autorise 00–05 et le test isolé de 10 avec fixtures contractuelles,
sans exécution de 06–09 ni qualification biologique.

`d21_kinship_bench_v1.json` configure le banc technique de pedigree IID :
graines distinctes de simulation/manque/bootstrap, sept relations, 112 unités,
224 sujets, 22 autosomes, deux strates artificielles, trois profils appariés.
Les contrôles LD/monomorphisme/région et les seuils 05/06/07 sont explicites.
CAL-DEC-005 autorise le registre synthétique seulement ; le runner garde
la validation des exclusions 07 en attente. Protocole :
`protocoles/D21_KINSHIP_PRESPECIFICATION.md`. Aucun réglage des seuils après
résultats, aucune référence aux données de patients ni au cache.

`d22_structure_bench_v1.json` fixe le banc de structure : 128 sujets dont
108 références et 20 projetés, 64 loci IID par autosome, contrastes 0,2/0,8,
graines distinctes, mélange et atypiques. L'annotation publique D19 est
épinglée par SHA-256, les génotypes sont explicitement synthétiques selon
le rang modulo 3. CAL-DEC-006 autorise seulement le registre et les revues
techniques de fixture ; le rôle `calibration_contract_fixture` ne peut pas
être remplacé par un prétendu rôle scientifique. Trois axes PCA et deux
axes d'alerte sont fixés ; modifier ce domaine exige une nouvelle version.

`d23a_phasing_panel_v2.json` est la configuration active de préparation :
8 trios, 16 sujets isolés, 22 autosomes liés à la puce, cible technique D19
épinglée, référence strictement D10. Les graines 2301/2302 restent privées,
2303 configure le futur phasage. La v1 est un brouillon conservé et refusé
par le chargement actif ; seule la clôture bornée `--from-preparation` peut
en reprendre les artefacts inchangés. Aucun outil de phasage n'est lancé.

`d23_phasing_bench_v1.json` fixe les trois observations appariées de D23 sur
la seule préparation v2, les graines distinctes de phasage/manque/bootstrap,
les strates MAF/écart physique, l'interpolation, les seuils de densité et les
timeouts. Il ne contient ni chemin local codé en dur, ni seuil scientifique
de commutation, ni approbation du point 07. Le protocole est sous
`protocoles/D23_PHASING_PRESPECIFICATION.md`.

`d24_founder_chr5_pedigree_v1.json` prépare trois unités familiales
indépendantes sur le domaine chr5 existant : deux CEU et une
`ADMIX_CEU_YRI`. Les cinq configurations positives fixent les âges 50, 100,
150 et 200 générations, dont deux graines indépendantes à 100. Elles ont été
matérialisées sous CAL-DEC-011 et ne doivent pas être remplacées par le panel
D14 historique à deux familles.

`d24_founder_haplotype_bench_v1.json` ferme le contrat d'évaluation de
l'étape 13 : profil primaire `MAXIMUM`, sensibilités `BASELINE`,
`DENSITY_50` et `MISSING_1PCT`, au moins trois familles par positif, cinq
positifs et cinq négatifs, IC de Wilson à 95 % et absence de reclassement
IBS→IBD. Les seuils de performance restent non scellés. CAL-DEC-011 autorise
l'exécution bornée du banc, sans qualification automatique.
La matrice ferme cinq positifs d'âges 50/100/150/200 générations (avec deux
réplications à 100) et cinq négatifs dont les porteurs de cible indépendants
ne partagent aucun segment injecté. Chaque scénario est apparié sur les quatre
profils ; seul `MAXIMUM` alimente les critères principaux.
La construction exécutée utilise 487 marqueurs régionaux pour `MAXIMUM`, le
sous-ensemble D23a plus la cible protégée pour `BASELINE`, une rétention
déterministe de 50 % et un manque de 1 % hors cible. Les haplotypes restent le
phasage-vérité D14 : les erreurs de phasage D23 ne sont pas injectées ici.

`d25_variant_age_bench_v1.json` ferme l'évaluation D25 : modèle primaire
`CORRELATED`, référence `INDEPENDENT`, dix âges, 3/5/10 familles, 200
réplications par cellule et graines distinctes pour simulation et bootstrap.
Il fixe le profil exact, la censure des bras à 3 cM, les échelles de carte 0,95
et 1,05, les cinq positifs D24 et la politique de phase
`NOT_EVALUABLE_WITHOUT_MATCHED_SHAPEIT5_OUTPUT`. Il interdit toute promotion
automatique de `calibrated` ou `founder_effect_proven`.

`d25a_gamma_audit_v1.json` ferme l'audit correctif D25a sur de nouvelles
graines. Il croise dix âges, 3/5/10 familles, trois fractions de chemin partagé
et les contrats `FULL_BREAKPOINTS`/`PAIRWISE_SHARED_TRUNCATED`. Il épingle le
script R officiel par URL et SHA-256, couvre ses quatre régions de corrélation
avec cinq fixtures et fixe les critères de parité, biais et couverture. La
méthode historique, la candidate fidèle à R et l'oracle complet ont des
identifiants distincts ; aucun résultat ne permet une promotion automatique.

`d25b_gamma_censoring_v1.json` ferme la nouvelle vraisemblance de censure sur
les graines 25201/25202 : dix âges, 3/5/10 familles, 400 réplications et trois
profils généalogiques. Il impose exactement un plus long bras censuré par côté,
interdit la correction de partage fortuit et fixe une fixture analytique ainsi
que B01–B07. La promotion de production et les deux verdicts globaux sont
forcés à `false`.

`d26_ld_roh_bench_v1.json` ferme les six observations non phasées distinctes,
les alias de phase, la cible chr5 et les paramètres réels PLINK des étapes 15
et 16. Les manifests D16/D17 sont épinglés par SHA-256. Le seuil d'appariement
ROH vaut 50 %, mais les seuils de performance LD/ROH restent explicitement
non scellés ; aucune promotion scientifique ou de production n'est possible.

`d27_ancestry_bench_v1.json` ferme le banc du cœur 16A : dix composantes,
référence D10 de 1 403 sujets disjointe des donneurs, contraste de mélange à
deux strates sous alias, huit profils appariés D17 et bootstrap de 2 000
réplications. Les six sources sont épinglées par SHA-256. Les seuils
CAL-ANC-001/002, la promotion de production et les verdicts globaux restent
explicitement non scellés ou faux.

`d28_enrichment_bench_v1.json` ferme le banc du cœur 16B sur les quarante cas
D24. Il fixe trois unités familiales, deux marqueurs par flanc, 10 000 tirages
nuls évaluables, la graine 2801 et alpha=0,05 pour les seules métriques du
banc. Le seuil de classification de production et les seuils d'acceptation de
performance restent nuls ; toute promotion globale est interdite.
