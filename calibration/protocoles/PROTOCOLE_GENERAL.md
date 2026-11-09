# Protocole général de simulation et de calibration

**Identifiant :** `CAL-PROTOCOL-GENERAL`
**Version :** 1.0.0
**Date :** 28 octobre 2026
**Statut :** `APPROVED`
**Étape de développement :** D03

## 1. Objet du protocole

Ce protocole définit comment générer des données génétiques synthétiques à
vérité connue, les analyser à l'aveugle avec le pipeline Effet_fondateur2 et
mesurer séparément la performance de chaque étape activée.

Il poursuit trois objectifs :

1. vérifier l'intégrité technique du parcours de `00` à `19` ;
2. qualifier scientifiquement en priorité les étapes `12`, `13`, `14` et
   `16C` ;
3. produire les éléments nécessaires à une décision humaine de qualification,
   sans modifier rétroactivement un résultat réel.

La simulation valide un parcours analytique dans les scénarios étudiés. Elle ne
prouve pas un effet fondateur réel, ne remplace pas une validation biologique
externe et ne permet aucune conclusion clinique individuelle.

## 2. Sources d'autorité

L'ordre d'autorité est le suivant :

1. `AGENTS.md` pour la sécurité, le développement et la reproductibilité ;
2. `calibration/CONVENTIONS_DEVELOPPEMENT.md` pour leur application locale ;
3. le présent protocole pour le plan expérimental ;
4. `calibration/protocoles/CRITERES_ACCEPTATION.yaml` pour les décisions
   mesurables ;
5. `calibration/protocoles/REGISTRE_DECISIONS.md` pour les décisions humaines ;
6. `docs/modules/*.md` et les schémas versionnés pour les contrats de chaque
   étape ;
7. les configurations scellées et leurs empreintes pour chaque exécution.

Une source de rang inférieur ne peut pas assouplir une règle d'une source de
rang supérieur. Toute divergence bloque la porte de passage concernée jusqu'à
sa résolution documentée.

## 3. Périmètre scientifique

### 3.1 Qualification prioritaire

Conformément à `CAL-DEC-001`, la qualification approfondie porte d'abord sur :

- `12` : exactitude du phasage et orientation de la variation synthétique ;
- `13` : détection et délimitation du segment fondateur IBS ;
- `14` : biais, erreur et couverture de la datation Gamma ;
- `16C` : détection IBD explicite, faux supports et précision des bornes.

### 3.2 Validation de l'ensemble du pipeline

Toutes les autres étapes activées restent soumises à :

- un contrôle technique de leurs entrées, sorties, effectifs et empreintes ;
- au moins un scénario fonctionnel à vérité connue lorsque leur comportement
  est testable ;
- un rapport indépendant ;
- une conclusion explicite parmi `PASS`, `FAIL`, `NOT_EVALUABLE` et
  `REVIEW_REQUIRED`.

Aucune étape ne reçoit une qualification par héritage de la réussite d'une
autre étape.

### 3.3 Hors périmètre

Le protocole n'autorise pas :

- l'utilisation d'une donnée individuelle ou d'un résultat de `dock6_real` ;
- l'estimation d'un risque clinique ou d'une causalité ;
- l'attribution d'une origine géographique ou ethnique individuelle ;
- la transformation de IBS, LD, ROH ou proximité PCA en preuve IBD ;
- la recherche post hoc d'un seuil rendant un résultat réel positif ;
- l'écrasement ou la modification d'un run historique.

## 4. Responsabilités et séparation des rôles

Le processus distingue quatre responsabilités logiques, même si une même
personne peut les exercer à des moments différents :

1. **générateur** : crée les données synthétiques et la vérité privée ;
2. **opérateur aveugle** : exécute le pipeline sans accès à la vérité ;
3. **évaluateur** : compare les sorties immuables à la vérité après clôture ;
4. **responsable scientifique** : approuve le protocole, scelle les seuils et
   signe la décision finale.

Les actions de chaque rôle sont enregistrées. Le responsable scientifique ne
peut pas déclarer une méthode qualifiée si un rapport requis manque ou si une
violation d'aveuglement n'est pas résolue.

## 5. Architecture et isolation

Les documents, schémas et configurations sont versionnés sous `calibration/`.
Le code et les tests sont versionnés sous `src/effet_fondateur/calibration/` et
`tests/calibration/`. Tous les artefacts générés sont écrits sous :

```text
calibration_workspace/<calibration_run_id>/
├── blinded_inputs/
├── private_truth/
├── pipeline_runs/
├── evaluations/
├── rapports/
├── figures/
├── cache/
├── manifest.json
└── checksums.sha256
```

Le pipeline analysé reçoit uniquement `blinded_inputs/`. Il ne reçoit jamais
un chemin, une variable, un manifest ou un champ pointant vers
`private_truth/`. L'évaluateur reçoit la vérité uniquement après publication
atomique et empreinte du run aveugle.

L'accès à `data/input/dock6_real/`, direct ou par lien symbolique, est une
violation critique. Il arrête le run et interdit toute qualification issue de
ce run.

## 6. Dessin expérimental

### 6.1 Micro-panel technique

Le micro-panel comprend 24 à 60 individus synthétiques. Il contrôle formats,
contrats, erreurs, déterminisme, aveuglement et raccord entre étapes. Un
chromosome réduit peut être utilisé pour les tests ordinaires ; un test
d'intégration identifié couvre les 22 autosomes.

Le micro-panel ne produit aucun seuil de performance ni aucune conclusion de
qualification scientifique.

### 6.2 Pilote de calibration

Le pilote de calibration comprend 500 individus. Il sert à :

- vérifier la faisabilité des scénarios ;
- estimer temps, mémoire et stockage ;
- estimer les distributions des métriques ;
- proposer les seuils `pilot_then_seal` selon une règle écrite avant le run ;
- identifier les scénarios non évaluables ou insuffisamment représentés.

Les valeurs du pilote ne peuvent pas être utilisées pour évaluer ensuite ce
même panel comme preuve indépendante de performance.

### 6.3 Pilote de validation aveugle

Le pilote de validation comprend 500 individus nouveaux, avec graines,
haplotypes de base, pedigrees, transmissions et erreurs distincts. Le code, les
seuils et les critères sont scellés avant son lancement.

Après le début du run, une modification de code, seuil ou règle de décision
invalide le caractère confirmatoire. Le run peut rester diagnostique, mais un
nouveau panel aveugle est alors nécessaire.

### 6.4 Panels de qualification

La qualification comprend 2 500 individus de calibration et 2 500 individus de
validation. Elle n'est autorisée qu'après une décision `GO` documentée sur les
ressources et la qualité du pilote.

Chaque panel de qualification contient au minimum :

- cinq strates de référence et une composante mélangée ;
- 100 trios ;
- 200 paires indépendantes par classe de parenté évaluée ;
- 3 000 paires négatives indépendantes pour l'IBD ;
- 400 répétitions positives par longueur IBD ;
- 500 répétitions par âge simulé pour la datation.

Ces minima sont des conditions d'évaluabilité, pas des garanties de puissance.

## 7. Indépendance des panels

Les panels de calibration et de validation doivent être disjoints pour :

- les identifiants synthétiques ;
- les graines maîtresses et dérivées ;
- les haplotypes donneurs ;
- les pedigrees ;
- les événements de recombinaison ;
- les erreurs injectées ;
- les répétitions statistiques.

Les références publiques et cartes versionnées peuvent être communes, mais
leur partage est déclaré. Les donneurs simulés et la référence de phasage sont
partitionnés afin d'éviter une validation optimiste. Toute intersection
interdite produit `FAIL`.

## 8. Construction de la vérité simulée

### 8.1 Haplotype et structure de population

Les haplotypes de base doivent préserver, dans le domaine simulé, les
fréquences alléliques, la LD et la relation entre positions `bp` et `cM`. Les
sources publiques, assemblages, cartes, licences, versions et empreintes sont
figés avant génération.

Les strates et mélanges sont générés avec proportions connues. Ces étiquettes
servent uniquement à l'évaluation privée et ne sont jamais transmises comme
vérité au pipeline.

### 8.2 Pedigrees et transmissions

Les simulations couvrent les sujets non apparentés, trios et relations du
premier au troisième degré. Les transmissions mendéliennes sont exactes sauf
erreur explicitement injectée et inscrite dans la vérité.

L'unité familiale indépendante est définie avant l'analyse. Des paires qui
partagent un individu ou un même événement de transmission ne sont pas traitées
comme des répétitions statistiques indépendantes.

### 8.3 Fondateur et recombinaisons

Un haplotype fondateur synthétique est transmis pendant un âge connu. La vérité
enregistre, pour chaque copie pertinente :

- l'âge en générations ;
- les événements de recombinaison ;
- les bornes gauche et droite en `bp` et `cM` ;
- la longueur totale en `cM` ;
- l'unité familiale ;
- la copie portant la variation synthétique.

Les âges minimaux étudiés sont 5, 10, 15, 20, 30, 40, 50, 75, 100 et 150
générations, avec 3, 5 et 10 familles indépendantes.

### 8.4 IBD, scénarios négatifs et ROH

Les longueurs IBD minimales étudiées sont 1, 1,5, 2, 3 et 5 cM. Chaque
longueur possède des scénarios positifs et des scénarios négatifs sans segment
injecté à la cible.

Les ROH sont générées séparément de l'IBD interindividuel. Elles peuvent être
absentes, distinctes ou superposées au locus fondateur selon un scénario
déclaré. Une ROH ne devient jamais une longueur d'entrée Gamma.

### 8.5 Imperfections de mesure

Les profils d'erreur font varier au minimum :

- le taux de génotypes manquants ;
- le taux d'erreurs génotypiques ;
- les erreurs ponctuelles de phase ;
- les commutations de phase ;
- la densité des marqueurs ;
- l'incertitude de carte et de borne lorsqu'elle est étudiée.

Chaque défaut est versionné et enregistré dans `private_truth/`. Aucun défaut
non déclaré ne doit être introduit par la conversion ou l'export.

## 9. Aveuglement et déroulement d'un run

Le déroulement obligatoire est :

1. valider la configuration et toutes les graines ;
2. créer un identifiant de run immuable ;
3. générer séparément vérité et entrées aveugles ;
4. calculer leurs empreintes SHA-256 ;
5. vérifier l'absence de chemin ou champ de vérité dans les entrées ;
6. exécuter le pipeline sur les seules entrées aveugles ;
7. publier atomiquement les sorties, manifestes et empreintes ;
8. fermer le run aveugle ;
9. ouvrir la vérité à l'évaluateur ;
10. produire les évaluations et rapports indépendants ;
11. consolider sans recalcul dans `CAL-REP-E2E` ;
12. soumettre `CAL-REP-DEC` à une décision humaine.

L'ouverture de la vérité, les empreintes avant et après ouverture et toute
déviation sont consignées dans le manifest du run.

## 10. Populations d'analyse

Chaque métrique doit déclarer son dénominateur et sa population :

- **technique** : tous les fichiers, variants ou individus attendus ;
- **évaluable** : unités satisfaisant les conditions écrites avant analyse ;
- **positive** : unités avec événement vrai injecté ;
- **négative** : unités sans événement vrai injecté ;
- **primaire** : scénario, méthode et population préspécifiés ;
- **sensibilité** : variation déclarée du scénario primaire.

Une exclusion est autorisée uniquement si son motif était préspécifié ou si
elle est enregistrée comme déviation. Les exclus restent comptés et décrits.
Le dénominateur ne peut pas être modifié après observation du résultat.

## 11. Métriques et formules

### 11.1 Classification

À partir de `TP`, `TN`, `FP` et `FN` :

```text
sensibilité = TP / (TP + FN)
spécificité = TN / (TN + FP)
taux de faux positifs = FP / (FP + TN)
précision = TP / (TP + FP)
balanced accuracy = (sensibilité + spécificité) / 2
```

Un dénominateur nul produit `NOT_EVALUABLE`, jamais zéro ni un succès.

### 11.2 Mesures continues

Pour une vérité `theta` et une estimation `theta_hat` sur `n` répétitions
indépendantes :

```text
biais = moyenne(theta_hat - theta)
MAE = moyenne(abs(theta_hat - theta))
RMSE = sqrt(moyenne((theta_hat - theta)^2))
couverture à 95 % = proportion des intervalles à 95 % contenant theta
```

Le biais reste signé. Le MAE et la RMSE ne doivent pas être appelés « biais ».
Les erreurs de bornes sont publiées séparément à gauche et à droite, en `bp` et
en `cM`, puis éventuellement résumées par une moyenne absolue préspécifiée.

### 11.3 Phasage

Le taux de commutation est calculé sur les transitions hétérozygotes
évaluables contre la phase vraie. Les erreurs ponctuelles de variants rares et
l'orientation de la variation synthétique sont publiées séparément. Les
homozygotes ne sont pas transformés artificiellement en observations de phase.

### 11.4 Datation

La datation est évaluée par âge vrai et par effectif de familles indépendantes.
Sont obligatoires : biais signé, MAE, RMSE, couverture de l'intervalle à 95 %,
longueur d'intervalle et stabilité leave-one-family-out.

Le modèle Gamma historique reste inchangé comme référence. Toute correction
de calibration ou de biais reçoit un identifiant distinct et publie en parallèle
le résultat brut.

### 11.5 IBD explicite

L'unité primaire est la paire indépendante dans une répétition de simulation.
Le faux support concordant exige que les conditions préspécifiées des deux
outils soient remplies dans un scénario négatif. Sensibilité, spécificité,
erreur de borne et concordance sont publiées par scénario et longueur.

La fréquence de faux support est accompagnée de son intervalle binomial. Pour
`16C`, sa borne supérieure à 95 % doit rester inférieure ou égale à `0,05`
dans chaque scénario qualifié.

## 12. Incertitude et dépendance

- Les proportions reçoivent un intervalle de Wilson ou binomial exact à 95 %.
- Les métriques continues reçoivent un intervalle calculé sur l'unité
  indépendante de simulation, par méthode analytique justifiée ou bootstrap
  hiérarchique reproductible.
- Le bootstrap conserve ensemble les observations d'une même famille,
  transmission ou répétition.
- Les paires partageant un individu ne sont pas comptées comme observations
  indépendantes dans un intervalle naïf.
- Les résultats sont stratifiés par âge, longueur, densité, erreur, ascendance
  simulée et nombre de familles lorsque ces facteurs varient.
- Un résultat agrégé n'est publié que si sa pondération est préspécifiée ; les
  strates restent toujours visibles.
- Aucun score composite de preuve fondatrice n'est calculé.

## 13. Cycle des critères d'acceptation

### 13.1 Contrats fixes

Les critères de sécurité, intégrité, concordance exacte, absence de fuite,
nombre d'autosomes et seuil maximal de faux support `16C` sont fixés avant tout
pilote. Ils ne sont pas appris à partir des résultats.

### 13.2 Seuils `pilot_then_seal`

Avant le pilote, chaque seuil de performance doit posséder :

- une métrique et une unité ;
- un sens d'amélioration ;
- une population et des strates ;
- une taille minimale ;
- une méthode d'intervalle ;
- une règle écrite de dérivation du seuil ;
- une règle pour les valeurs manquantes.

Le résultat ponctuel observé au pilote ne peut pas être simplement recopié
comme seuil. Le choix doit tenir compte de l'incertitude, de la pertinence
scientifique et d'une marge explicite. Sans règle de dérivation approuvée, le
critère reste `REVIEW_REQUIRED`.

Les seuils candidats, leurs justifications et le SHA-256 du fichier de critères
sont inscrits dans le registre avant le panel aveugle. Après scellement, toute
modification invalide la validation confirmatoire et exige un nouveau panel.

### 13.3 Panel aveugle

Le panel aveugle est jugé uniquement contre les critères scellés. Un échec
n'autorise pas un assouplissement post hoc. Il conduit à `FAIL`,
`NOT_EVALUABLE` ou à une nouvelle version du système suivie d'une nouvelle
validation indépendante.

## 14. Valeurs manquantes et statuts

- `PASS` : tous les critères obligatoires évaluables respectent leurs règles ;
- `FAIL` : au moins un critère obligatoire évalué enfreint sa règle ;
- `NOT_EVALUABLE` : une entrée, vérité, taille ou mesure indispensable manque ;
- `REVIEW_REQUIRED` : une décision humaine ou un seuil non scellé est requis.

Une valeur manquante n'est jamais imputée pour produire un succès. La raison,
le dénominateur attendu, le dénominateur observé et l'action suivante sont
publiés.

La priorité d'agrégation est `FAIL`, puis `NOT_EVALUABLE`, puis
`REVIEW_REQUIRED`, puis `PASS`. Une violation de sécurité critique reste
`FAIL`, même si d'autres mesures sont non évaluables.

## 15. Arrêts, écarts et reprises

Le run s'arrête immédiatement en cas de :

- lecture ou copie d'une donnée réelle interdite ;
- fuite de vérité vers le pipeline ;
- collision entre calibration et validation ;
- empreinte incohérente ;
- tentative d'écrasement d'un run ;
- configuration ou schéma invalide ;
- assemblage, carte ou référence non concordants ;
- code ou seuil modifié pendant une validation aveugle.

Un échec d'outil, un timeout, une densité insuffisante ou un effectif trop faible
produit un artefact d'erreur audité et généralement `NOT_EVALUABLE`. Il ne doit
pas être masqué par une reprise silencieuse. Une reprise conserve l'identité du
run seulement si le contrat de reprise prouve que les entrées, paramètres et
signatures sont inchangés.

Tout écart au protocole reçoit un identifiant, une date, une justification, un
impact et une décision de conservation ou d'invalidation.

## 16. Rapports, provenance et conservation

Chaque étape activée produit son rapport autonome `CAL-REP-*`, avec :

- configuration, graines et versions ;
- vérité pertinente et population évaluée ;
- entrées aveugles consommées ;
- métriques, dénominateurs et intervalles ;
- comparaison aux critères préspécifiés ;
- exclusions, écarts, limites et décision ;
- empreintes des sources et produits.

Le rapport bout en bout consolide les décisions sans recalculer les statistiques
des rapports sources. Les rapports générés, panels et vérités restent hors Git.
Les documents versionnés conservent uniquement les contrats et décisions
agrégées sans identifiant individuel.

## 17. Décision et portes de passage

- **G0** autorise le développement après approbation du périmètre.
- **G1** autorise le micro-panel après validation des schémas, de l'isolation,
  des conventions et des tests unitaires.
- **G2** autorise le pilote de calibration après réussite du micro-panel.
- **G3** autorise le pilote aveugle après scellement des critères.
- **G4** autorise les 5 000 individus après revue des ressources et du pilote.
- **G5** autorise la décision scientifique lorsque tous les rapports existent.
- **G6** autorise l'intégration au pipeline après approbation de
  `CAL-REP-E2E` et `CAL-REP-DEC`.

Une porte échouée ne peut pas être contournée par modification manuelle d'un
statut. La décision finale précise pour chaque méthode : qualifiée, qualifiée
avec restrictions, non qualifiée ou non évaluable.

## 18. Limites générales

La validité est limitée aux assemblages, cartes, références, densités,
versions d'outils, profils d'erreur et scénarios simulés. Un changement de l'un
de ces éléments peut exiger une recalibration.

Une simulation ne couvre jamais toute la diversité biologique ou technique.
Une bonne performance simulée n'exclut pas un biais partagé entre générateur et
analyse. Les références publiques peuvent ne pas représenter certaines
histoires démographiques. Les petits segments IBD et les faibles effectifs
familiaux peuvent rester non évaluables malgré un protocole correct.

## 19. Références méthodologiques

1. Gandolfo LC, Bahlo M, Speed TP. *Dating Rare Mutations from Small Samples
   with Dense Marker Data*. Genetics. 2014;197:1315–1327.
   DOI : `10.1534/genetics.114.164616`.
2. Zhou Y, Browning SR, Browning BL. *A Fast and Simple Method for Detecting
   Identity-by-Descent Segments in Large-Scale Data*. American Journal of Human
   Genetics. 2020;106:426–437. DOI : `10.1016/j.ajhg.2020.02.010`.
3. Talts S, Betancourt M, Simpson D, Vehtari A, Gelman A. *Validating Bayesian
   Inference Algorithms with Simulation-Based Calibration*. 2018.
   `arXiv:1804.06788`.
4. Documentation interne versionnée sous `docs/modules/`.

## 20. Approbation

La revue humaine du 28 octobre 2026 a confirmé :

- le périmètre et les limites ;
- les unités indépendantes et dénominateurs ;
- l'aveuglement et les conditions d'arrêt ;
- la grille scientifique minimale ;
- le cycle de scellement des seuils ;
- les règles d'incertitude et de valeurs manquantes ;
- l'absence d'utilisation de données DOCK6 pour le réglage.

L'approbation est enregistrée par `CAL-DEC-002` dans
`REGISTRE_DECISIONS.md`. Elle autorise la poursuite du développement et des
tests unitaires. Elle n'autorise pas à elle seule un micro-panel, un pilote ou
une qualification : les portes G1 à G4 restent applicables.
