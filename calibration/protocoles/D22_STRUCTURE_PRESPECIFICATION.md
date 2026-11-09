# D22 — Structure et gel : préspecification

Version 1.0.0, 30 octobre 2026. Méthode `iid_structure_component_bench_v1`.
Protocole fixé avant génération et exécution. Aucun seuil ajusté pour PASS.

## Décision et domaine

L'utilisateur reporte explicitement l'étude de sensibilité de 07 après
l'ensemble de l'évaluation et autorise D22. Le FAIL D21 reste ouvert ;
CAL-DEC-006 ne l'approuve ni ne le corrige. Aucun fichier D21 n'est modifié.

D18 régional ne suffit pas aux contrôles de structure, et D21 ne contient
pas de mélanges et d'atypiques dédiés. D22 utilise un banc de composants
indépendant : 128 sujets, 22 autosomes, 64 loci IID par autosome ; 40 sujets
de fréquence 0,2, 40 de fréquence 0,8, 24 mélanges locus par locus 50/50,
quatre atypiques rares à contraste orthogonal et 20 sujets projetés dont
les génotypes initiaux copient des références. Les quatre atypiques de
référence et quatre copies projetées constituent les positifs connus.
Ce ne sont pas des populations humaines, segments d'ascendance ou données
1000G. La vérité de mélange est la probabilité de tirage, pas une proportion
génomique réalisée estimable sans erreur par PCA.

Quatre cas appariés : base, génotypes projetés inversés, métadonnées
descriptives permutées, manque 5 %. Le même identifiant technique synthétique
permet les comparaisons d'invariance ; aucune identité réelle n'est présente.
Le locus 0 est monomorphe dans la référence mais variable chez les projetés ;
le locus 1 manque dans toute la référence. Ils contrôlent les exclusions PCA.

## Exécution et fixtures honnêtes

00–03 sont exécutées réellement pour le registre et la conversion avec PLINK.
08 et 09 sont ensuite appelées via leurs points d'entrée publics, hors du
manifest 00–03. Le panel 06, la proposition 07, le QC 05 et l'audit cible 04
sont des **fixtures déclarées `D22_UPSTREAM_CONTRACT_FIXTURE`**, pas des étapes
exécutées. Le panel est la copie exacte du triplet 03, sans prétention de
pruning. La référence est préspecifiée (108 sujets), jamais sélectionnée
après lecture de la vérité ou des performances. Les 20 projetés sont exclus
de cette référence par construction, pas par un appel KING prétendu.

Les génotypes cibles T/T, T/A et A/A sont explicitement synthétiques selon
le rang modulo 3, indépendamment du groupe, de la clinique ou de l'ascendance.
L'annotation publique épinglée D19 est réutilisée, sans appel expérimental.
Les rôles témoins, clinique, lots et groupes de FID sont techniques et ne
révèlent pas les strates simulées. Les FID ne constituent pas ici des pedigrees.

Les revues de 09 sont des **fixtures de décision de test**, rôle explicite
`calibration_contract_fixture`, jamais des approbations scientifiques humaines.
L'autorisation préspecifie leur portée : accepter seulement les exclusions
techniques du nouveau banc pour vérifier le contrat. La liaison SHA aux
sorties 08 est renseignée après leur clôture, sans changer les règles de
revue ni utiliser les étiquettes vraies. Le registre des décisions permet
d'auditer ce mécanisme. Aucune décision réelle ou D21 n'est résolue.

Deux gels 09 identiques en données/paramètres sont exigés par cas, avec des
répertoires de publication distincts. Un contrôle sans revues, une mauvaise
empreinte, une liste incomplète et une date sans fuseau doivent être refusés,
sans artefact de cohorte publié. Les logs et codes d'erreur sont empreintés.

## Oracle 08

La PCA (analyse en composantes principales) est recalculée sans les helpers
de production : dosages physiques lus depuis BED, fréquences et imputation
sur la seule référence, matrice de Gram, décomposition symétrique `eigh`
au lieu de la SVD de production. Loadings et scores sont reconstruits ;
signe déterministe au loading maximal, tolérances numériques configurées.
Les p-values chi-deux à deux degrés de liberté sont vérifiées par la formule
indépendante `exp(-distance²/2)`. Les valeurs absentes ne valent jamais zéro.

Comparaisons : tous les axes, loadings, fréquences, statuts, effectifs,
distances et propositions ; invariance des axes et scores de référence quand
seuls les projetés changent ; invariance numérique aux rôles descriptifs.
La séparation A/B est une distance de centroïdes divisée par dispersion
intra-groupe ; la position des mélanges est exprimée sur l'axe des centroïdes.
Ce ne sont pas des classifications géographiques ni estimations ADMIXTURE.
Les comparaisons entre cas alignent les signes des axes par leurs scores de
référence et les loadings selon les allèles physiques : PLINK peut inverser
A1/A2 lorsque les génotypes projetés changent. Ce changement de codage ne
constitue pas un changement du sous-espace ; aucune rotation arbitraire
n'est autorisée pour masquer un écart d'axes non dégénérés.

CAL-POP-001/002 : sensibilité des huit atypiques injectés et taux d'alertes
des non-injectés, avec intervalles Wilson descriptifs. Les copies projetées
ne sont pas des réplications indépendantes ; publier aussi les résultats sur
les seuls fondateurs indépendants de référence, sans intervalle abusif sur
le total dépendant. Pas de seuil de qualification scellé. CAL-POP-003 porte
sur les assertions agrégées publiées, pas un journal exhaustif du système.

## Oracle 09 et réussite

Reconstruction indépendante des sept cohortes, des motifs ordonnés, sources,
unités et représentants ; contrôle exact des sept keep, comptes, décisions
et résumés. Les génotypes ne sont jamais inférés depuis les étiquettes.
CAL-POP-004 exige zéro différence entre les deux gels sur toutes les tables
de décisions/membres/résumés/keep et le rapport, hors audit horodaté et
enveloppe de publication. Une décision reproductible peut néanmoins être
scientifiquement erronée : ce PASS ne qualifie pas la sélection KING.

Zéro discordance fonctionnelle attendue ; défaut observé => FAIL, jamais
NOT_EVALUABLE maquillé. Les taux de détection restent REVIEW_REQUIRED ;
les faibles effectifs et limites sont publiés. Rapports indépendants
CAL-REP-08, CAL-REP-09 et CAL-DEV-D22, HTML/PDF selon report-formater.
`calibrated=false`, `founder_effect_proven=false`. Ni pipeline entier ni
cohorte d'étude validés ; 07 reste à reprendre après l'évaluation globale.

Sources : contrats `docs/modules/population_structure.md` et
`docs/modules/cohort_freeze.md` ; [PLINK 1.9, PCA et projection sur un
sous-ensemble](https://www.cog-genomics.org/plink/1.9/strat).
