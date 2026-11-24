# Plan unifié du panel de simulation et de la calibration

**Projet :** Effet_fondateur2

**Version :** 1.0

**Date :** 12 novembre 2026

**Statut :** plan directeur approuvé ; R00 à R10 réalisés ; aucune nouvelle
analyse scientifique n'est autorisée avant validation du panel maître

## 1. Décision et périmètre

La campagne de calibration conduite jusqu'à D28 est abandonnée comme chaîne de
qualification. Elle a utilisé des bancs utiles séparément, mais construits avec
des individus, vérités, densités, régions et paramètres qui ne composent pas
une release scientifique commune.

La nouvelle campagne partira d'un panel maître conçu à partir des besoins
complets des étapes du pipeline `00` à `19`. Toutes les vues analysées seront
dérivées de cette même release et resteront reliées aux mêmes identités,
pedigrees, graines et vérités.

La remise à zéro concerne uniquement les contrats, configurations, protocoles
et résultats de calibration. Le code sous `src/`, les tests sous `tests/`, la
documentation fonctionnelle, les configurations du logiciel et les données du
projet sont conservés.

Après la dernière étape seulement, un rapport HTML/PDF unique rassemblera les
résultats, courbes, équations, limites et modifications logicielles proposées.
Pendant la campagne, chaque étape produira uniquement un fichier Markdown de
résultats et ses données sources.

## 2. Objectif et limites scientifiques

### 2.1 Objectif

Construire une validation reproductible capable de :

1. produire une vérité génétique synthétique cohérente ;
2. dériver des observations réalistes sans révéler cette vérité ;
3. exécuter le pipeline actif dans ses conditions réelles ;
4. comparer chaque sortie à une vérité indépendante ;
5. mesurer performances, biais, incertitudes et domaines non évaluables ;
6. qualifier ou refuser séparément chaque méthode ;
7. identifier les modifications nécessaires au logiciel sans les appliquer
   silencieusement pendant la validation.

### 2.2 Ce que la calibration ne démontrera pas

- Elle ne prouvera pas à elle seule un effet fondateur dans une étude réelle.
- Elle ne remplacera pas une validation biologique ou externe.
- Elle ne transformera pas IBS, ROH, LD ou proximité PCA en preuve IBD.
- Elle ne garantira pas des intervalles plus étroits.
- Elle ne permettra pas de régler un seuil à partir du résultat d'une étude
  réelle particulière.
- Elle ne permettra pas de déclarer un outil qualifié sur une réussite
  purement technique.

Jusqu'à une décision scientifique finale explicite :

- `calibrated=false` ;
- `founder_effect_proven=false` ;
- `production_promotion=false`.

## 3. Règles non négociables

### 3.1 Isolation

- Ne jamais lire ni copier de données réelles pour construire ou régler le
  panel ; `data/input/dock6_real/` est notamment un chemin réel protégé.
- Ne pas utiliser le cache interne comme source de simulation.
- Utiliser uniquement des identifiants synthétiques.
- Conserver les données générées, vérités privées et sorties volumineuses sous
  `calibration_workspace/`, ignoré par Git.
- Refuser les chemins interdits, les liens symboliques sortant du périmètre et
  les collisions d'identifiants.

### 3.2 Indépendance

- Les sujets de calibration et de validation auront des graines, transmissions
  et composantes de pedigree distinctes.
- Les références de phasage seront disjointes des sujets évalués.
- Les seuils seront réglés uniquement sur la calibration.
- Les critères seront scellés avant l'ouverture de la validation aveugle.
- Aucun seuil ni algorithme ne sera modifié après observation de cette
  validation.

### 3.3 Séparation vérité–analyse

- `private_truth/` contiendra la vérité inaccessible au pipeline.
- `blinded_inputs/` contiendra uniquement ce que le pipeline peut observer.
- Les manifests du pipeline ne référenceront jamais `private_truth/`.
- L'évaluateur ouvrira la vérité seulement après clôture et empreinte du run.
- La personne ou le processus qui exécute l'analyse ne recevra pas les classes
  positives, négatives ni les âges vrais.

### 3.4 Reproductibilité

Chaque release enregistrera :

- toutes les graines aléatoires ;
- les sources, licences, versions et assemblages ;
- les versions de Python, PLINK, KING, bcftools, SHAPEIT5, Hap-IBD,
  Refined IBD et Java lorsqu'ils sont utilisés ;
- les cartes génétiques et leurs unités ;
- les commandes et paramètres effectifs ;
- les SHA-256 des entrées, vérités, vues et sorties ;
- les relations producteur–artefact–consommateur ;
- les écarts au protocole et les reprises.

Un dossier finalisé ne pourra pas être écrasé ou complété silencieusement.

## 4. Sources d'autorité

En cas de divergence, l'ordre suivant s'applique :

1. `AGENTS.md` pour la sécurité, les modifications et la reproductibilité ;
2. le présent plan pour l'organisation de la calibration ;
3. `docs/modules/*.md` pour les contrats scientifiques du pipeline ;
4. les schémas et configurations recréés après conception du panel ;
5. le code réellement exécuté sous `src/effet_fondateur/` ;
6. les tests sous `tests/` ;
7. les résultats clos de la nouvelle campagne.

Le comportement actuel du logiciel est un objet à évaluer, pas une vérité
scientifique. Une divergence entre documentation, code et sortie sera notée
comme anomalie avant toute proposition de correction.

## 5. Produit central : une release de panel maître

### 5.1 Définition

Une release est un ensemble immuable qui relie :

- les sujets synthétiques et leurs rôles ;
- les haplotypes d'origine ;
- la structure de population et les mélanges ;
- les pedigrees et transmissions ;
- les mutations synthétiques ;
- les recombinaisons et segments fondateurs ;
- les segments IBD et ROH indépendants du fondateur ;
- les erreurs de mesure ;
- les vérités de génotype, phase, parenté, ascendance et segments ;
- toutes les observations dérivées ;
- les partitions de calibration et de validation ;
- les paramètres et empreintes.

Une étape ne pourra pas créer son propre panel parallèle. Si un besoin manque,
la release sera refusée avant exécution et son contrat sera corrigé de manière
prospective.

### 5.2 Couches obligatoires

La release contiendra quatre couches séparées :

1. **vérité biologique simulée** : haplotypes, transmissions, phases,
   recombinaisons, segments et âges ;
2. **vérité des événements de mesure** : marqueurs supprimés, génotypes
   perturbés, données manquantes, erreurs de phase et duplications ;
3. **observations aveugles** : fichiers réellement consommables par le
   pipeline ;
4. **sorties du pipeline** : artefacts produits sans accès aux vérités.

Les évaluations constitueront une cinquième couche créée seulement après
clôture du run.

### 5.3 Partitions

R06 distingue les partitions d'étude exclusives `CALIBRATION`,
`VALIDATION_BLIND` et `RESOURCE_ONLY`, ainsi que trois familles de références
`PHASING_REFERENCE`, `ANCESTRY_REFERENCE` et `NULL_REFERENCE`, chacune séparée
entre calibration et validation. `NEGATIVE_CONTROLS` est une strate privée
présente indépendamment dans les deux partitions évaluées, pas un pool de
sujets partagé.

L'indépendance est contrôlée par unité de séparation, qui réunit toute identité
biologique, composante de pedigree, événement fondateur, paire IBD et duplicat
technique reliés. Elle n'est jamais déduite des seuls identifiants :

`CALIBRATION ∩ VALIDATION_BLIND = ∅`

`(CALIBRATION ∪ VALIDATION_BLIND ∪ RESOURCE_ONLY) ∩ REFERENCES = ∅`

### 5.4 Vues dérivées

Les vues suivantes seront prévues avant la génération :

- génotypes denses et vérité haplotypique ;
- observation à la densité maximale réellement exploitable de la puce ;
- sous-densités contrôlées dérivées des mêmes marqueurs ;
- vue autosomale sur les 22 autosomes ;
- vues régionales autour de plusieurs cibles synthétiques ;
- vue avec pedigrees et vue sans pedigrees, toutes choses égales par ailleurs ;
- niveaux contrôlés de données manquantes et d'erreurs de génotypage ;
- phase vraie, génotypes non phasés et phase estimée ;
- scénarios fondateurs positifs et scénarios nuls ;
- segments IBD positifs et paires négatives ;
- ROH injectées et fonds sans ROH injectée ;
- vues d'ascendance globale et locale.

Chaque vue héritera des mêmes identifiants, rôles, vérités et empreintes de
release. Une comparaison de densité ou de pedigree utilisera exactement les
mêmes sujets, références et paramètres hors facteur étudié.

### 5.5 Décisions de conception R05

Le panel maître porte l'identifiant de conception
`panel_master_generic_v1`. Il couvre les 22 autosomes en `GRCh38` et contient
plusieurs variations cibles synthétiques. Un run du pipeline n'en sélectionne
qu'une ; aucune cible, aucun gène, chromosome,
allèle ou identifiant n'est codé en dur dans le moteur.

Tous les sujets évalués, ancêtres latents et sujets de référence de calibration
sont synthétiques. Aucun génotype individuel réel, public ou privé, ne sert de
donneur, de modèle génératif ou de fond nul. Une carte génétique publique et
des conventions d'assemblage peuvent être utilisées comme sources
méthodologiques versionnées, mais jamais des données individuelles réelles pour
construire ou régler le panel. Un profil de production utilisant une référence
publique reste extérieur à la release de calibration et ne peut pas fournir
une vérité.

R05 fige les entités, mécanismes et facteurs du panel. Il ne fixe pas les
effectifs, répétitions, taux d'erreur et budgets, qui relèvent de R07 ; les
partitions et matérialisations des vues relèvent de R06 ; les critères de
décision relèvent de R08 ; les chemins, schémas et configurations relèvent de
R09. Aucun élément laissé à ces étapes ne pourra ajouter un mécanisme
biologique après observation d'un résultat.

### 5.6 Composants immuables du panel maître

Les identifiants ci-dessous sont des objets de conception, pas des schémas déjà
créés.

| Composant | Contenu figé à R05 | Vérités R04 alimentées |
|---|---|---|
| `PM-SOURCE` | assemblage `GRCh38`, 22 autosomes, cartes, catalogues synthétiques valides et sentinelles de source | `T-SOURCE` |
| `PM-SUBJECT` | identités biologiques synthétiques, échantillons techniques, rôles, phénotypes orthogonaux et composantes | `T-SUBJECT` |
| `PM-HAPLOTYPE` | deux haplotypes denses ordonnés par sujet, variants bialléliques, LD et phase vrais | `T-GENOTYPE`, `T-LD`, `T-PHASE` |
| `PM-TARGET` | catalogue générique de cibles, allèles, copies porteuses et provenance moléculaire explicite | `T-TARGET` |
| `PM-PEDIGREE` | graphes d'ascendance observés et latents, transmissions et degrés 1 à 4 | `T-SUBJECT`, `T-KINSHIP`, `T-PHASE` |
| `PM-POPULATION` | trois populations latentes synthétiques, mélanges globaux, mosaïques locales et atypiques | `T-POPULATION`, `T-LD` |
| `PM-REFERENCE` | références synthétiques de phasage, positionnement et fond nul, mutuellement disjointes | `T-PHASE`, `T-POPULATION`, `T-FOUNDER` |
| `PM-FOUNDER` | événements ancestraux positifs, récurrents, nuls et censurés, reliés aux familles | `T-FOUNDER`, `T-TARGET`, `T-PHASE` |
| `PM-IBD` | segments IBD injectés indépendamment du fondateur, paires négatives et limites haplotypiques | `T-IBD` |
| `PM-ROH` | segments autozygotes, fardeaux, positions relatives à la cible et fonds sans injection | `T-ROH` |
| `PM-MEASURE` | erreurs de source, identité, génotype, lot, pedigree, carte et phase appliquées après la vérité | `T-MEASURE`, `T-SCENARIO` |
| `PM-SEED` | dérivation hiérarchique des graines par domaine, entité, partition future et répétition | `T-SCENARIO` et provenance de toutes les vérités |
| `PM-TECHNICAL` | cas d'intégrité, reprise, statuts et attentes de rendu sans donnée individuelle | `T-SOURCE`, `T-RENDER` |

Chaque objet biologique porte un identifiant stable indépendant de son futur
chemin. Les relations entre composants utilisent ces identifiants, jamais un
nom de fichier. Une release incomplète pour un composant est refusée : une
étape ne peut pas générer localement le morceau manquant.

### 5.7 Sujets et pedigrees

#### 5.7.1 Unités et identités

- une identité biologique possède exactement deux haplotypes par autosome et
  appartient à une seule composante d'ascendance ;
- un échantillon technique pointe vers une identité biologique ; les duplicats
  techniques ont des identifiants d'échantillon distincts mais le même donneur
  vrai ;
- une personne observée appartient à une seule composante connexe de pedigree ;
  les ancêtres latents ne sont jamais présentés comme échantillons ;
- les identifiants ne codent ni partition, ni population, ni statut, ni
  génotype cible ;
- le phénotype synthétique est généré orthogonalement au génotype cible et ne
  peut jamais servir à attribuer une copie ALT ;
- sexe, parents, génération, lot et rôle sont explicites ; une valeur inconnue
  est distincte d'une absence biologique.

#### 5.7.2 Bibliothèque de pedigrees

| Type | Structure figée | Besoin couvert |
|---|---|---|
| `PED-UNRELATED` | identité sans lien dans la profondeur simulée | contrôles indépendants et fonds nuls |
| `PED-TRIO` | deux parents et un enfant, avec variantes parent absent ou non génotypé | Mendel, phasage et concordance parent–enfant |
| `PED-SIBSHIP` | deux parents et au moins deux germains | premier degré hors parent–enfant et transmissions partagées |
| `PED-MULTIGEN` | au moins trois générations observées ou latentes | recombinaisons, phase et degrés éloignés |
| `PED-DEGREE` | paires connues de degrés 2, 3 et 4 | classification KING et indépendance |
| `PED-DUPLICATE` | duplicat technique d'une identité ou paire monozygote de deux identités | identité, QC et classe duplicat |
| `PED-CRYPTIC` | lien génétique entre familles déclarées distinctes | apparentement cryptique et revue humaine |
| `PED-FOUNDER-FAMILY` | famille observée récente descendant d'une branche ancestrale cachée | unités familiales indépendantes des étapes 13–16C |

Tous les graphes sont acycliques, les parents précèdent les enfants et chaque
méiose possède ses points de recombinaison vrais. Les familles d'un même
événement fondateur ne partagent aucun lien récent dans le pedigree observé :
elles sont reliées uniquement par l'ancêtre latent de l'événement. Les degrés
attendus proviennent du graphe vrai et restent séparés des classes estimées par
KING.

L'identifiant de composante récente et l'identifiant d'événement fondateur sont
distincts. Les familles d'un même événement sont les unités familiales des
analyses 13–16C, mais l'événement complet, avec toutes ses familles, formera une
unité indivisible lors de l'affectation aux partitions R06 afin d'empêcher une
fuite entre calibration et validation.

### 5.8 Haplotypes, cibles et références

#### 5.8.1 Fond haplotypique

Le fond dense contient des variants autosomiques bialléliques, des positions bp
et cM strictement ordonnées, des allèles REF/ALT non ambigus et une phase vraie
sur les 22 autosomes. Les haplotypes sont générés par un modèle paramétrique
synthétique de fréquences et de blocs de LD ; ses paramètres ne sont pas appris
sur des génotypes individuels réels. Les génotypes non phasés, densités de puce
et fichiers destinés au pipeline seront des dérivés R06, jamais la vérité
maître.

Les recombinaisons utilisent la carte déclarée, sans approximation
`1 cM = 1 Mb`, et conservent à la fois bp, cM, haplotype parental et génération.
Les vérités de LD sont calculées sur les haplotypes denses avant erreur de
mesure, séparément par population et ensemble de sujets défini.

#### 5.8.2 Catalogue des cibles

| Rôle de cible | Contrainte de sélection avant génération | Finalité |
|---|---|---|
| `TGT-STANDARD` | région interne au chromosome, carte dense et fenêtre faisable | parcours principal 04 et 11–16C |
| `TGT-LOW-RECOMB` | environnement à taux de recombinaison synthétique faible | segments longs et sensibilité à la carte |
| `TGT-HIGH-RECOMB` | environnement à taux élevé | segments courts et résolution des bornes |
| `TGT-LEFT-EDGE` | marge gauche insuffisante pour au moins un âge prévu | censure gauche et `INFEASIBLE_WINDOW` |
| `TGT-RIGHT-EDGE` | marge droite insuffisante pour au moins un âge prévu | censure droite et `INFEASIBLE_WINDOW` |

Les cinq rôles sont portés par des autosomes distincts. Les coordonnées exactes
seront résolues de façon déterministe à R09 depuis la carte et l'assemblage
alors empreintés, puis scellées avant toute génération. Chaque cible possède un
identifiant de projet synthétique, REF, ALT et une étiquette de gène optionnelle
qui peut rester vide. Dans une composante affectée à une cible, les quatre
autres cibles restent homozygotes REF afin d'éviter les porteurs multiples
implicites.

La résolution ne comporte aucun choix manuel : les ancres de la carte sont
triées par chromosome numérique puis bp, les rôles sont résolus dans l'ordre de
la table et la première ancre éligible d'un autosome encore inutilisé est
retenue. Le taux local est la pente cM/Mb dans une fenêtre physique centrée de
10 Mb ; `LOW` et `HIGH` utilisent respectivement les quartiles inférieur et
supérieur de toutes les ancres éligibles, tandis que `STANDARD` reste dans
l'intervalle interquartile et possède les deux marges nécessaires à l'âge de
5 générations. Les deux rôles `EDGE` prennent la première ancre qui échoue
uniquement la marge du côté indiqué. Une absence de candidat rend la release
`INFEASIBLE_TARGET_CATALOG` au lieu d'autoriser un remplacement discrétionnaire.

#### 5.8.3 Références de calibration

| Référence | Contenu | Séparation obligatoire |
|---|---|---|
| `REF-PHASE-SYN` | haplotypes synthétiques phasés et carte compatible SHAPEIT5 | aucun donneur commun avec les sujets évalués |
| `REF-ANCESTRY-SYN` | individus synthétiques étiquetés uniquement par population latente | disjointe de `REF-PHASE-SYN` et des sujets évalués |
| `REF-NULL-SYN` | individus et haplotypes non porteurs pour le fond externe 16B | disjointe des deux autres références et de tous les fonds internes |

Ces références sont générées par des flux de graines distincts et ne contiennent
aucune famille évaluée. La référence d'ascendance ajuste seule ses axes ; les
sujets d'étude sont uniquement projetés. Une référence publique éventuelle du
profil de production ne remplace aucune de ces trois références et ne contribue
ni aux seuils ni aux performances de calibration.

### 5.9 Structure de population et mélanges

Le modèle comporte trois populations latentes `POP-A`, `POP-B` et `POP-C`.
Ces codes désignent uniquement des processus synthétiques et ne correspondent
à aucune population humaine, origine géographique ou identité ethnique. Elles
diffèrent par des fréquences alléliques et une structure de LD connues, tout en
partageant le même assemblage et le même catalogue de variants.

Les états globaux comprennent les trois populations non mélangées et, pour
chaque paire de populations, les proportions synthétiques `25/75`, `50/50` et
`75/25`. Les mosaïques d'ascendance locale sont générées par segments avec
bornes vraies et peuvent différer du mélange global autour d'une cible. Une
classe atypique est toujours une population ou un mélange déjà défini, projeté
hors de l'ensemble de référence interne de la vue ; elle n'est jamais créée en
déplaçant artificiellement des coordonnées PCA.

Les sujets apparentés et leurs représentants indépendants partagent la même
vérité populationnelle. Le modèle PCA de l'étape 08 doit être ajusté uniquement
sur les représentants indépendants de la future vue. Pour 16A, les axes et
centroïdes sont ajustés uniquement sur `REF-ANCESTRY-SYN`. Les proportions et
segments vrais servent à évaluer la stabilité du positionnement ; ils
n'autorisent pas à exiger du logiciel une estimation de proportions que sa
méthode ne produit pas.

### 5.10 Événements fondateurs et variations cibles

Chaque événement possède un identifiant, une cible, un âge vrai en générations,
un haplotype ancestral, une copie mutante, des branches descendantes, des
recombinaisons, des familles observées et un mécanisme de censure. Les âges
figés sont `5`, `10`, `15`, `20`, `30`, `40`, `50`, `75`, `100` et `150`
générations ; les nombres de familles indépendantes étudiés sont `3`, `5` et
`10`. La section 7 fixe les répétitions de chaque cellule.

| Classe | Construction vraie | Interprétation attendue du scénario |
|---|---|---|
| `FND-POSITIVE` | une mutation ancestrale unique transmise sur un même haplotype à plusieurs familles récentes indépendantes | positif fondateur simulé |
| `FND-RECURRENT` | le même ALT apparaît indépendamment sur des haplotypes non apparentés | porteurs sans fondateur unique |
| `FND-BACKGROUND-IBS` | partage local fortuit sans mutation ancestrale commune | contrôle des faux partages IBS |
| `FND-SINGLE-FAMILY` | transmission dans une seule composante observée | insuffisance d'unités indépendantes, jamais faux positif multi-familles |
| `FND-NONE` | aucune mutation ni segment fondateur injecté | contrôle négatif complet |
| `FND-CENSORED` | événement positif dont un ou deux bras sortent de la fenêtre observable | censure gauche, droite ou double explicitement connue |

Les événements positifs sont simulés dans un graphe ancestral latent, puis
transmis jusqu'aux familles observées avec recombinaison. Les longueurs
biologiques complètes et les longueurs observables sont conservées séparément.
Les scénarios principaux possèdent deux bras observables lorsque la géométrie
le permet ; les défis censurés sont analysés séparément et ne contribuent pas
silencieusement aux métriques principales.

Chaque cellule inclut des copies porteuses `H1`, `H2` et, lorsque biologiquement
cohérent, `BOTH`. La copie ALT est déterminée exclusivement par le génotype et
la transmission synthétiques. Les familles supplémentaires, homozygotes ou
individus d'une même famille ne deviennent jamais des réplications
indépendantes. La signature de fond exclut toujours la cible.

### 5.11 IBD explicite et ROH

#### 5.11.1 IBD indépendant du scénario fondateur

Les événements `PM-IBD` sont générés après le fond haplotypique mais avant les
erreurs de mesure. Ils copient un segment d'un haplotype ancestral synthétique
vers les haplotypes désignés et enregistrent chromosome, bornes bp/cM, longueur,
copies, paire et présence de la cible. Leurs longueurs vraies sont `1`, `1,5`,
`2`, `3` et `5 cM`.

Les classes obligatoires sont : segment contenant la cible, segment sur le même
chromosome mais ne contenant pas la cible, segment sur un autre autosome,
segments multiples concurrents, limites légèrement différentes entre copies,
et absence de segment. Une paire négative ne possède aucun segment vrai d'au
moins `1 cM` dans la région testée ; ce statut est vérifié sur les haplotypes
denses, pas déduit de l'absence d'injection.

Les injections IBD et les événements fondateurs sont disjointes dans les
cellules principales. Leur coexistence n'est permise que dans une classe de
défi `COMPOUND`, explicitement exclue des estimations primaires. L'appel par
Hap-IBD ou Refined IBD reste aveugle à ces classes et à la variation cible.

#### 5.11.2 ROH

Un ROH vrai est créé en copiant un intervalle haplotypique entre les deux copies
d'un même sujet avant erreur de génotypage. Les classes structurelles figées
sont `NONE`, `SHORT`, `MEDIUM`, `LONG` et `MULTIPLE_BURDEN`; leurs longueurs
numériques et répétitions sont fixées en 7.4 en bp et cM compatibles avec la
densité. Chaque classe existe avec la cible à l'intérieur, à l'extérieur et sur
un autre chromosome lorsque cela est pertinent.

Le panel inclut aussi : fonds sans ROH injecté, ROH fragmenté par manque ou
erreur, deux ROH proches susceptibles d'être fusionnés, et cible hétérozygote
dans un appel tolérant un hétérozygote. Le fardeau vrai compte chaque segment,
sa longueur et le dénominateur autosomique explicite ; aucun dénominateur de
taille de génome n'est implicite. Les ROH ne redéfinissent ni l'IBS fondateur,
ni l'IBD interindividuel.

### 5.12 Événements de mesure et défis techniques

Les erreurs sont appliquées à une copie dérivée de la vérité biologique. Chaque
événement conserve valeur avant/après, sujet, variant ou fichier concerné,
lot, ordre d'application et graine. La vérité biologique n'est jamais modifiée
pour faire correspondre l'observation perturbée.

| Famille | Mécanismes figés | Étapes principalement exercées |
|---|---|---|
| `ERR-SOURCE` | colonne absente, encodage invalide, fichier tronqué, assemblage discordant, lien symbolique interdit | `00`, `01` |
| `ERR-IDENTITY` | doublon d'identifiant, source partagée, échange d'échantillons, duplicat technique non déclaré | `02`, `03`, `05`, `07` |
| `ERR-ALLELE` | REF/ALT inversés, allèle incompatible, position dupliquée, cible absente ou ambiguë | `03`, `04`, `11`, `12`, `16A` |
| `ERR-MISSING` | manque par cellule, sujet, variant, chromosome et lot | `03`, `05`, `06`, `10`, `12`, `15`, `16`, `16C` |
| `ERR-GENOTYPE` | dosage retourné, hétérozygotie excédentaire, erreur mendélienne et exception cible contrôlée | `04`, `05`, `10`, `12`, `16` |
| `ERR-BATCH` | taux d'appel et spectre d'erreurs différents entre lots sans changement biologique | `05`, `08`, `10`, `17` |
| `ERR-PEDIGREE` | parent absent, filiation erronée, famille déclarée séparée ou fusionnée | `02`, `07`, `09`, `12`, `17` |
| `ERR-PHASE` | commutations, copie porteuse inversée et région non fiable dans une vue de défi | `12`, `12A`, `13`, `16B`, `16C` |
| `ERR-MAP` | trou, ordre non monotone, fenêtre insuffisante et extrapolation demandée | `11`, `12`, `14`, `16C` |
| `ERR-INTEGRITY` | empreinte altérée, sortie absente, reprise incohérente et artefact dynamique non déclaré | `00`, clôture, `17`, `18`, `19` |

Chaque mécanisme possède un témoin `NONE`, les niveaux ordonnés fixés en 7.4
et au moins un cas de rejet déterministe. Les cellules
primaires ne modifient qu'une famille d'erreur à la fois. Les combinaisons sont
réservées aux défis `COMPOUND` et ne servent pas à attribuer causalement un
échec à un facteur unique. Les absences d'outil, timeouts et sorties externes
malformées relèvent du harnais simulé R09, pas d'une mutation des génotypes du
panel.

L'ordre immuable est : sources et fond biologique ; populations, haplotypes
ancestraux et pedigrees ; événements cibles, fondateurs, IBD et ROH à leur
niveau ancestral déclaré ; méioses et transmissions ; haplotypes et génotypes
denses des sujets observés ; échantillons techniques ; erreurs de mesure ;
dérivation des vues R06. Changer cet ordre crée une nouvelle version de
conception.

### 5.13 Graines et déterminisme

La règle publique de dérivation est :

`seed128 = first128(SHA-256("Effet_fondateur2|panel_master_generic_v1|" + partition_namespace + "|" + domain + "|" + entity_id + "|" + replicate_id))`.

La chaîne est encodée en UTF-8 avec des identifiants ASCII canoniques ;
`first128` désigne les 16 premiers octets du digest interprétés comme un entier
non signé big-endian. Ces conventions évitent deux graines différentes pour un
même identifiant logique.

Les générateurs Python utiliseront `PCG64DXSM` avec cette graine de 128 bits.
Un outil exigeant un entier positif sur 31 bits recevra
`1 + seed128 mod (2^31 - 2)`. La version de NumPy et l'algorithme effectifs
seront inscrits dans la release ; un changement d'algorithme invalide les
artefacts descendants.

Les domaines sont `SOURCE`, `HAPLOTYPE`, `POPULATION`, `REFERENCE`,
`PEDIGREE`, `MEIOSIS`, `TARGET`, `FOUNDER`, `IBD`, `ROH`, `MEASUREMENT`,
`VIEW`, `PHASING`, `NULL_DRAW` et `SENSITIVITY`. Aucun générateur séquentiel
n'est partagé entre domaines. Les namespaces des futures partitions de
calibration et validation sont distincts ; une même entité ne peut donc pas
recevoir la même graine dans les deux partitions.

L'identifiant de répétition est attribué avant génération selon la section 7.
Une répétition rejetée reste enregistrée avec son motif et sa graine ; elle
n'est pas remplacée sous le même identifiant. Les tirages nuls et sensibilités
n'utilisent jamais la graine d'un événement positif. Le pipeline reçoit
uniquement les graines nécessaires à ses propres outils, jamais les flux privés
ayant créé la vérité.

### 5.14 Couverture des vues et des étapes

#### 5.14.1 Support des vues R04

| Vue future | Composants maîtres nécessaires | Élément figé à R05 |
|---|---|---|
| `V-SOURCE` | `PM-SOURCE`, `PM-TECHNICAL` | exports synthétiques valides et sentinelles contrôlées |
| `V-TARGET` | `PM-TARGET`, `PM-PEDIGREE` | appels explicites, transmissions et provenances moléculaires |
| `V-DENSITY` | `PM-HAPLOTYPE`, `PM-SOURCE` | même vérité dense et mêmes sujets avant sous-échantillonnage |
| `V-MEASUREMENT` | `PM-MEASURE`, `PM-HAPLOTYPE` | événements ordonnés sur une copie dérivée |
| `V-PEDIGREE` | `PM-PEDIGREE`, `PM-SUBJECT` | même composante avec déclaration complète, absente ou erronée |
| `V-REGIONAL` | `PM-HAPLOTYPE`, `PM-TARGET`, `PM-FOUNDER` | fenêtres et vérités locales reliées à la cible sélectionnée |
| `V-AUTOSOMAL` | `PM-HAPLOTYPE`, `PM-PEDIGREE` | sujets et phase vrais identiques sur les 22 autosomes |
| `V-PHASING-REFERENCE` | `PM-REFERENCE`, `PM-HAPLOTYPE` | références synthétiques disjointes et haplotypes vrais |
| `V-FOUNDER-NULL` | `PM-FOUNDER`, `PM-REFERENCE` | fonds internes/externes, récurrents et sans événement |
| `V-ROH-IBD` | `PM-ROH`, `PM-IBD` | positifs, négatifs et défis composés étiquetés |
| `V-SENSITIVITY` | `PM-SEED`, `PM-MEASURE` | appariement déterministe ne changeant qu'un facteur |

R06 décidera quels sujets et artefacts apparaissent dans chaque vue ; elle ne
pourra pas inventer une vérité ou un mécanisme absent de cette table.

#### 5.14.2 Fermeture des besoins 00–19

| Étape | Composants R05 couvrant ses besoins |
|---|---|
| `00` | `PM-SOURCE`, `PM-SEED`, `PM-TECHNICAL` pour intégrité, altération et reprise |
| `01` | `PM-SOURCE`, `PM-TECHNICAL` pour sources valides et un rejet par règle |
| `02` | `PM-SUBJECT`, `PM-PEDIGREE`, `PM-MEASURE` pour registre, collisions et filiations |
| `03` | `PM-HAPLOTYPE`, `PM-SUBJECT`, `PM-MEASURE` pour conversion, ordre, allèles et manque |
| `04` | `PM-TARGET`, `PM-PEDIGREE`, `PM-MEASURE` pour appels REF/ALT et Mendel |
| `05` | `PM-HAPLOTYPE`, `PM-SUBJECT`, `PM-MEASURE` pour défauts individu, variant et lot |
| `06` | `PM-HAPLOTYPE`, `PM-POPULATION`, `PM-SOURCE` pour MAF, LD, chromosomes et densité |
| `07` | `PM-PEDIGREE`, `PM-SUBJECT` pour duplicats, degrés 1–4, cryptiques et indépendants |
| `08` | `PM-POPULATION`, `PM-PEDIGREE`, `PM-MEASURE` pour axes, mélanges, atypiques et projections |
| `09` | `PM-SUBJECT`, `PM-PEDIGREE`, `PM-TARGET` pour cohortes attendues et composantes |
| `10` | `PM-HAPLOTYPE`, `PM-TARGET`, `PM-MEASURE` pour QC par cohorte et exception cible |
| `10A` | `PM-HAPLOTYPE`, `PM-SUBJECT` pour univers complet sur 22 autosomes |
| `11` | `PM-TARGET`, `PM-SOURCE`, `PM-FOUNDER` pour carte, fenêtres faisables et censure |
| `12` | `PM-HAPLOTYPE`, `PM-PEDIGREE`, `PM-REFERENCE`, `PM-TARGET` pour phase régionale vraie |
| `12A` | `PM-HAPLOTYPE`, `PM-REFERENCE`, `PM-SEED` pour phase vraie sur 22 autosomes |
| `13` | `PM-FOUNDER`, `PM-TARGET`, `PM-HAPLOTYPE` pour IBS positif, nul, récurrent et censuré |
| `14` | `PM-FOUNDER`, `PM-PEDIGREE` pour âges, bras biologiques/observés et unités familiales |
| `15` | `PM-HAPLOTYPE`, `PM-POPULATION`, `PM-MEASURE` pour vérités `r²`, `D′`, MAF et manque |
| `16` | `PM-ROH`, `PM-TARGET`, `PM-MEASURE` pour segments, fardeau, fusion et fragmentation |
| `16A` | `PM-POPULATION`, `PM-REFERENCE`, `PM-HAPLOTYPE` pour positionnements global/local |
| `16B` | `PM-FOUNDER`, `PM-REFERENCE`, `PM-PEDIGREE`, `PM-SEED` pour fonds nuls et familles |
| `16C` | `PM-IBD`, `PM-HAPLOTYPE`, `PM-REFERENCE`, `PM-TARGET` pour appels aveugles genome-wide |
| `17` | `PM-MEASURE`, `PM-SEED`, tous composants affectés pour paires à facteur unique |
| `18` | `PM-TECHNICAL` et vérités agrégées pour statuts, unités, légendes et intégrité |
| `19` | `PM-TECHNICAL` et provenance de tous les composants pour faits et traçabilité |

Les familles d'erreurs `E-TOOL` sont couvertes par le futur harnais R09 ; elles
n'exigent pas de modifier le panel biologique. Toutes les autres familles
d'erreurs R03 possèdent au moins un composant et un scénario de conception.

### 5.15 Points transmis et porte R05

R05 n'autorise encore aucune génération. Les décisions suivantes sont
volontairement transmises aux étapes prévues :

- R06 : affectation disjointe des composantes aux partitions, filtres exacts,
  héritage et appariement des onze vues ;
- R07 : effectifs, répétitions, longueurs ROH numériques, taux d'erreur,
  densités, puissance et ressources ;
- R08 : métriques, seuils, intervalles, règles `NOT_EVALUABLE` et directions
  attendues ;
- R09 : coordonnées cibles résolues et scellées, schémas, configurations,
  générateur et harnais d'erreurs externes.

Deux écarts d'implémentation sont enregistrés sans modifier le logiciel :

1. les contrats actuels de phasage et de positionnement documentent une
   référence publique réelle ; le profil de calibration R09 devra accepter les
   trois références synthétiques de R05 sans modifier silencieusement le profil
   de production ni leur attribuer une identité publique fictive ;
2. les scripts exploratoires existants sous
   `simulation_genotype_famille/` et `src/effet_fondateur/calibration/` ne sont
   pas des générateurs autorisés de la nouvelle release ; leur éventuelle
   réutilisation à R09 exigera une adéquation démontrée à chaque composant R05.

La porte R05 est satisfaite au niveau de conception : sujets, références,
haplotypes, pedigrees, populations, cibles, événements fondateurs, IBD, ROH,
erreurs et dérivation des graines sont figés ; les quatorze vérités, onze vues
et 25 étapes R03–R04 sont couvertes avant génération. Toute extension de cette
grammaire crée `panel_master_generic_v2` et impose une nouvelle revue R05.

### 5.16 Espaces de partition R06

R06 fixe les espaces logiques exclusifs suivants. Ils sont attribués avant
toute génération ; aucun découpage aléatoire postérieur d'un panel déjà généré
n'est autorisé.

| Espace | Contenu autorisé | Usage |
|---|---|---|
| `STUDY/CALIBRATION` | sujets d'étude synthétiques et contrôles négatifs privés de calibration | réglage R11 uniquement |
| `STUDY/VALIDATION_BLIND` | sujets d'étude synthétiques et contrôles négatifs privés de validation | validation R12 uniquement |
| `STUDY/RESOURCE_ONLY` | micro-cas techniques sans estimation scientifique | ressources et faisabilité R07/R10 |
| `REFERENCE/PHASING/CALIBRATION` | `REF-PHASE-SYN` de calibration | phasage de calibration |
| `REFERENCE/PHASING/VALIDATION_BLIND` | `REF-PHASE-SYN` de validation | phasage de validation |
| `REFERENCE/ANCESTRY/CALIBRATION` | `REF-ANCESTRY-SYN` de calibration | positionnement de calibration |
| `REFERENCE/ANCESTRY/VALIDATION_BLIND` | `REF-ANCESTRY-SYN` de validation | positionnement de validation |
| `REFERENCE/NULL/CALIBRATION` | `REF-NULL-SYN` de calibration | fonds nuls de calibration |
| `REFERENCE/NULL/VALIDATION_BLIND` | `REF-NULL-SYN` de validation | fonds nuls de validation |
| `TECHNICAL/SENTINELS` | cas sans identité biologique ni métrique scientifique | intégrité des sources et rendus |

`NEGATIVE_CONTROLS` n'est pas une partition globale ni un synonyme du groupe
logiciel `controls_unrelated`. C'est une strate privée, générée indépendamment
dans `STUDY/CALIBRATION` et `STUDY/VALIDATION_BLIND`, qui couvre notamment les
classes `NEG-FND`, `NEG-IBD`, `NEG-ROH` et les témoins sans erreur injectée.
Un cache scientifique, un fichier de référence ou un tirage nul ne peut pas
être partagé entre ces deux partitions. Les installations et binaires des
outils peuvent être communs puisqu'ils ne contiennent aucune donnée
scientifique de la release.

### 5.17 Unité de séparation et affectation préalable

L'`isolation_unit_id` est la fermeture transitive, calculée comme des
composantes connexes, des relations suivantes : même identité biologique ou
duplicat technique ; même composante de pedigree récent ; même réplication
d'un événement fondateur, y compris ses familles descendantes ; même événement
IBD injecté et toutes ses paires ; même base biologique utilisée par des vues
appariées ; même donneur ou pool de référence. Une famille peut rester l'unité
d'analyse, mais toutes les familles d'une réplication fondatrice constituent
une seule unité de séparation.

La section 7 prescrit les nombres d'unités par
`(partition_id, scenario_cell_id, replicate_id)` avant R09. Le générateur R09
créera directement chaque unité dans l'espace de graines de sa partition, tel
que fixé en R05 ; il ne générera jamais un ensemble commun à scinder ensuite.
Une unité affectée est immuable et ne peut être déplacée pour équilibrer un
résultat observé. `RESOURCE_ONLY` est généré séparément et est exclu des seuils,
des estimations de performance et des décisions scientifiques.

Les invariants privés obligatoires sont :

`I_CAL ∩ I_VALID ∩ I_RESOURCE = ∅`, avec aussi chaque intersection deux à deux
vide pour les identités, composantes de pedigree, événements fondateurs,
événements IBD, duplicats et préfixes de graines ;

`R_TYPE,CAL ∩ R_TYPE,VALID = ∅` pour chacun des trois types de référence, et
tout donneur de référence est disjoint de tout sujet d'étude ;

`NEG_CAL ∩ NEG_VALID = ∅`, même si les deux strates suivent le même plan
factoriel. Des distributions identiques sont permises ; une identité, une
transmission, un événement, un segment vrai ou un tirage ne l'est jamais.

### 5.18 Pseudonymisation et séparation physique de la vérité

La release future devra séparer physiquement les espaces suivants, sans les
créer pendant R06 :

```text
release/
  governance/
  private_truth/<partition_id>/
  blinded_inputs/<partition_id>/<opaque_view_id>/
  runs/<closed_run_id>/
  evaluations/<evaluation_id>/
```

`governance/seed_ledger.tsv`, le registre des identités, les clés de
pseudonymisation, `master_truth_manifest.json` et les quatorze familles de
vérité sont privés. Le pipeline ne reçoit qu'un manifest de release assaini,
un inventaire exact, une configuration aveugle et une unique vue autorisée.
Les identifiants publics sont produits par une clé propre à la partition ; ils
restent stables entre vues appariées de cette partition mais ne permettent ni
le rapprochement entre partitions ni la récupération d'un identifiant maître.

Le processus du pipeline ne doit pouvoir ni monter, ni résoudre, ni lire
`private_truth/` ou `governance/`. R09 devra imposer une identité de processus
ou un espace de montage restreint, résoudre les chemins réels, refuser les
liens symboliques sortant des racines autorisées, assainir l'environnement et
scanner récursivement manifests et configurations pour les chemins ou champs
privés. Si cette isolation ne peut pas être démontrée, `READY_TO_RUN=false`.
L'inventaire d'entrée est une liste blanche fermée et le manifest du run prouve
les fichiers réellement consommés.

Sont interdits dans toute vue aveugle : âges fondateurs vrais, phase et
segments vrais, classes fondatrice/IBD/ROH, masque d'erreur injecté, verdict ou
direction attendus, population latente, graine maîtresse et chemin privé. Sont
autorisés seulement les observations que recevrait l'étape concernée : cible
moléculaire configurée, génotypes cibles mesurés et provenance permise,
pedigree déclaré éventuellement altéré, erreurs ou manque visibles, et phase
du seul panel de référence. Le statut clinique ou le groupe ne sert jamais à
inférer un génotype.

Un run est d'abord clos et empreinté. Ensuite seulement, un évaluateur isolé
reçoit les sorties closes et la famille de vérité strictement nécessaire. Il
ne reçoit jamais le magasin complet de vérité. Les identifiants de vues sont
opaques et leurs noms sémantiques restent dans la gouvernance privée.

### 5.19 Ordre canonique de dérivation

Toute vue est une transformation pure, déterministe, versionnée et empreintée
du panel maître. L'ordre immuable de la couche d'observation est : affecter la
partition, pseudonymiser, sélectionner la cible, appliquer les erreurs de
mesure au génotype dense, restreindre la densité, transformer la déclaration
du pedigree, encoder la source, puis dériver les projections régionale,
autosomale, de référence, de fond nul ou de sensibilité. Modifier cet ordre
crée une nouvelle version de release.

Les niveaux de densité sont emboîtés et héritent du même état biologique et du
même masque d'observation. Les génotypes d'étude restent non phasés ; seule la
référence synthétique peut exposer sa phase. Aucune transformation ne crée un
sujet, un allèle, une transmission ou un segment absent de la vérité maître.

### 5.20 Contrat exact des onze vues aveugles

| Vue | Filtres et transformations R06 | Interdits spécifiques | Étapes |
|---|---|---|---|
| `V-SOURCE` | partition sélectionnée, exports synthétiques ACPA/ChAS et métadonnées permises ; base valide ou défi `ERR-SOURCE` isolé | vérité, identité maître, classe attendue de l'erreur | `01`, `02`, `03` |
| `V-TARGET` | une cible configurée par vue, métadonnées moléculaires et appels observés issus de la source synthétique | génotype déduit du phénotype, copie porteuse vraie | `04`, `09`, `10` |
| `V-DENSITY` | `D_MAX ⊃ D_MID ⊃ D_LOW`, rang par hachage déterministe stratifié par chromosome et classe cM ; cible et ancres cartographiques protégées déclarées avant génération | sélection par résultat, changement de sujet ou de génotype | `03`, `05`, `06`, `08`, `10`, `15`, `16`, `17` |
| `V-MEASUREMENT` | base `NONE` puis une seule famille et un seul niveau d'erreur ; masques déterministes emboîtés | masque, taux vrai ou classe attendue exposés | `03`, `05`, `10`, `17` |
| `V-PEDIGREE` | `COMPLETE`, `PARENT_MISSING`, `NO_PEDIGREE` ou `MISDECLARED_CHALLENGE`, à génotypes constants | pedigree vrai des modes altérés | `02`, `07`, `09`, `12`, `17` |
| `V-REGIONAL` | une cible et une fenêtre déterministe dérivée de la carte ; impossibilité ou censure produisent un reçu explicite | extrapolation de carte, fenêtre fabriquée, phase vraie d'étude | `11`, `12`, `13`, `14`, `15`, `16`, `16A`, `16B`, `16C` |
| `V-AUTOSOMAL` | mêmes sujets, même ordre et observations non phasées sur les 22 autosomes ; garde de complétude des génotypes | sous-univers opportuniste ou dépendant de la cible | `05`, `06`, `07`, `08`, `10A`, `12A`, `16C` |
| `V-PHASING-REFERENCE` | sous-pool synthétique `REF-PHASE-SYN` ou `REF-ANCESTRY-SYN` propre à la partition, phasé, harmonisé et accompagné de sa carte | sujet d'étude comme donneur, référence partagée entre partitions, identité publique fictive | `12`, `12A`, `16A`, `16B` |
| `V-FOUNDER-NULL` | pool interne de non-porteurs et `REF-NULL-SYN` propre à la partition ; calendrier et graines de tirage préscrits | statistique observée, classe positive/nulle ou vérité fondatrice | `13`, `16B` |
| `V-ROH-IBD` | observations d'étude et sorties de phase estimée destinées aux appels aveugles genome-wide | segment vrai, paire positive, guidage par la cible | `16`, `16C` |
| `V-SENSITIVITY` | registre signé de deux runs clos ne différant que par un facteur autorisé | direction attendue, relance par l'étape `17`, paire inter-partitions | `17` |

Les comptes de variants, tailles de fenêtres, niveaux d'erreur et tailles de
pools ne sont pas fixés ici : R07 les dimensionne. Chaque vue contient son
parent, sa version de transformation, ses filtres, ses empreintes d'entrée et
de sortie, ainsi qu'un reçu explicite lorsque la dérivation est non évaluable.

### 5.21 Appariement à facteur unique

Deux membres d'une paire de sensibilité conservent partition, unités de
séparation, identifiants publics, cible, vérité biologique, graines de base,
versions logicielles et configuration, à l'exception exacte du facteur signé.
Les seules différences permises sont : sous-ensemble de variants pour
`DENSITY` ; masque ou intensité pour `MEASUREMENT` ; pedigree déclaré pour
`PEDIGREE` ; sous-pool de référence pour `REFERENCE` ; paramètre/graine de
phasage ou défi de phase aval pour `PHASE` ; état d'événement biologique prévu
par la cellule pour `FOUNDER`, `ROH` ou `IBD` ; paramètre ou graine d'outil pour
`TOOL/RNG`. Dans les trois facteurs biologiques, la paire est construite comme
un contrefactuel préscrit avant génération et toutes les autres vérités restent
identiques.

Un contrôle calcule le diff canonique des deux manifests : zéro ou plus d'un
facteur modifié invalide la paire. La direction attendue reste dans
`T-SCENARIO`, inaccessible à `17`. L'étape `17` ne relance rien et ne consomme
que des runs clos, y compris pour `16B` et `16C`. Une paire ne traverse jamais
les partitions.

### 5.22 Preuve de disjonction et porte R06

Avant tout `READY_TO_RUN`, R09 devra produire un reçu vérifiant au minimum :

1. unicité de la partition par `isolation_unit_id` et absence de partage des
   identités, pedigrees, événements, segments, duplicats et préfixes de graines ;
2. disjonction des six sous-pools de référence et des sujets d'étude ;
3. génération indépendante des contrôles négatifs de calibration et de
   validation ;
4. stabilité des pseudonymes dans une paire, mais impossibilité de les relier
   entre partitions ;
5. absence de champs de vérité, chemins privés et liens symboliques sortants
   dans l'inventaire aveugle ;
6. refus effectif de lecture de la vérité par le processus du pipeline ;
7. accès de l'évaluateur seulement après clôture, à la seule vérité autorisée ;
8. diff à facteur unique pour chaque paire et nouvelle vérification des
   sorties avant évaluation.

R06 est satisfaite au niveau de conception : les espaces sont exclusifs, les
unités liées sont indivisibles, les onze vues ont des transformations fermées,
les paires sont contrôlables et la vérité possède une frontière physique. La
réalisation et les tests de ces gardes restent exclusivement à R09. Trois
écarts sont transmis : implémenter l'isolation de processus et les gardes de
chemins (`GAP-R06-01`), les projections pseudonymisées et références
synthétiques propres aux partitions (`GAP-R06-02`), puis le registre signé et
le diff canonique à facteur unique couvrant `16B` et `16C` (`GAP-R06-03`).
R06 ne fixe aucun effectif, seuil ou verdict et ne génère aucun artefact.

## 6. Matrice détaillée des besoins des étapes 00 à 19

### 6.1 Règles communes de la matrice R03

La présente section ferme R03 au niveau des besoins. Elle ne fige ni les
effectifs finaux, qui relèvent de R07, ni les seuils d'acceptation, qui relèvent
de R08. Les identifiants couverts sont `00`, `01`, `02`, `03`, `04`, `05`,
`06`, `07`, `08`, `09`, `10`, `10A`, `11`, `12`, `12A`, `13`, `14`, `15`,
`16`, `16A`, `16B`, `16C`, `17`, `18` et `19`.

Toutes les étapes consomment des artefacts déclarés et empreintés. Toute étape
exécutée produit au minimum `stage_inputs.json`, `stage_outputs.json`,
`audit.json` et `checksums.sha256`. Après clôture du run, son évaluateur produit
le Markdown `calibration/resultats/ETAPE_<identifiant>.md`. Ces fichiers
communs ne sont pas répétés dans chaque ligne.

Les formats autorisés sont : YAML pour la configuration humaine, JSON
versionné pour les manifests et résumés structurés, TSV UTF-8 à en-tête
canonique pour les tables, BED/BIM/FAM pour les jeux PLINK, VCF bgzip ou BCF
indexé pour les génotypes phasés, cartes TSV ou PLINK MAP pour les positions
génétiques, SVG avec provenance JSON pour les figures, et texte brut ou gzip
uniquement pour les journaux et sorties natives déclarés. Aucun pickle, tableur
opaque ou format sans schéma ne constitue une source scientifique.

Les cardinalités suivantes seront enregistrées par vue et partition :

- `N_S` : sujets ; `N_F` : familles ou composantes indépendantes ;
- `N_V22` : variants autosomiques ; `N_VR` : variants régionaux ;
- `N_PAIR` : paires indépendantes ; `N_REF` : sujets de référence ;
- `N_REP` : répétitions indépendantes ; `N_DRAW` : tirages nuls évaluables.

Les mentions `micro`, `pilote` et `qualification` renvoient aux échelles de la
section 7. Les minima structurels indiqués ci-dessous sont obligatoires et les
valeurs numériques sont justifiées en 7.2–7.9. L'unité
indépendante est toujours déclarée : sujet, famille/composante, paire, variant,
segment, chromosome ou répétition. Un nombre de lignes n'est jamais assimilé
silencieusement à un effectif indépendant.

Les classes de ressources sont :

- `R0` : validation et agrégation Python, mémoire et disque linéaires dans la
  taille des tables ;
- `R1` : PLINK, KING ou bcftools sur puce, avec espace temporaire et temps
  proportionnels à `N_S × N_V` ou `N_PAIR` ;
- `R2` : références, phasage ou IBD sur 22 autosomes, parallélisme borné,
  cache contrôlé, mémoire/Java heap et espace temporaire mesurés ;
- `R3` : rendu déterministe depuis des tables closes, sans calcul scientifique.

La section 7.8 fixe les plafonds de CPU, RAM, stockage et temps. Dès le micro-panel,
chaque audit enregistrera durée murale, pic RSS, octets lus/écrits, volume de
sortie, nombre de threads, appels externes et statut du cache.

Les familles d'erreurs utilisées dans la matrice sont :

- `E-INT` : empreinte, schéma, chemin, ordre ou reprise incohérents ;
- `E-SRC` : source, version, licence, assemblage ou encodage invalide ;
- `E-ID` : doublon, collision, mauvais alignement ou fuite entre partitions ;
- `E-ALL` : coordonnée, allèle, orientation ou variant cible incohérent ;
- `E-MES` : manque, erreur de génotype, lot, hétérozygotie ou duplication ;
- `E-PED` : transmission, pedigree, degré ou indépendance incohérents ;
- `E-MAP` : carte absente, non monotone, extrapolée ou trop peu dense ;
- `E-PHA` : harmonisation, phase, commutation ou chromosome porteur erronés ;
- `E-SEG` : segment absent, fortuit, fusionné, fragmenté ou mal borné ;
- `E-CEN` : fenêtre ou observation censurée de façon non reconnue ;
- `E-LOWN` : effectif indépendant, densité ou fond nul insuffisant ;
- `E-TOOL` : outil absent, version incorrecte, timeout ou sortie malformée ;
- `E-LEAK` : vérité privée, sujet ou famille présent dans un domaine interdit ;
- `E-REN` : figure ou synthèse infidèle, incomplète ou recalculée.

### 6.2 Besoins détaillés 00 à 10

| Étape | Entrées aveugles et formats | Vérité privée | Unité et taille requises | Outils / ressources | Erreurs attendues | Données graphiques obligatoires |
|---|---|---|---|---|---|---|
| `00` initialisation | configuration YAML, inventaire JSON de release, environnement et versions | liste exacte des artefacts, signatures, SHA-256, états et résultats de reprise | un run par release et échelle ; au moins un arrêt/reprise et une altération contrôlée | Python/orchestrateur, `R0` | `E-INT`, `E-LEAK` | comptes d'artefacts par état et chronologie des transitions ; règles de décision, sans statistique inventée |
| `01` validation des sources | exports ACPA/ChAS synthétiques UTF-8 et catalogue JSON de sources | validité, assemblage, colonnes, couverture, version et empreinte de chaque source | sources valides et un cas de rejet par règle ; couverture des 22 autosomes et de la cible | Python, `R0` | `E-SRC`, `E-INT`, liens symboliques et fichiers tronqués | effectifs acceptés/rejetés par code, couverture par chromosome, aucune valeur génotypique |
| `02` registre | métadonnées synthétiques TSV, `source_inventory.tsv`, `source_qc.tsv` | identité, rôle, sexe, famille, filiation, lot et bijection sujet–source | `N_S` lignes attendues plus doublons, sources partagées, parents absents et discordances contrôlées | Python, `R0` | `E-ID`, `E-PED`, `E-LEAK`, génotype inféré du statut interdit | matrice attendu/observé des décisions, attrition par motif et composantes familiales agrégées |
| `03` conversion ACPA | exports validés, `samples.master.tsv` ; sorties attendues PLINK BED/BIM/FAM et descripteurs JSON | génotypes bruts, coordonnées, allèles, ordre des sujets et variants | matrice `N_S × N_V22`, mêmes sujets sur 22 autosomes ; vues intersection et union appariées | Python + PLINK 1.9, `R1` | `E-ID`, `E-ALL`, `E-MES`, `E-TOOL` | concordance génotypique, perte de variants, manque par sujet/chromosome et conflits d'allèles |
| `04` variant cible | métadonnées moléculaires YAML, génotypes explicites TSV, jeu PLINK du chromosome cible | génotype REF/ALT, source individuelle, transmission et cas non évaluables | tous les sujets exactement une fois ; porteurs, non-porteurs, manquants et trios dans chaque partition | Python + PLINK 1.9, `R1` | `E-ALL`, `E-ID`, `E-PED`, provenance clinique interdite, `E-TOOL` | matrice vérité/appel/exclusion, erreurs mendéliennes et effectifs par provenance explicite |
| `05` QC préliminaire | jeu PLINK autosomal, registre, vues appariées avec événements de mesure | événements injectés : manque, erreurs, duplications, hétérozygotie et lot | mêmes sujets/variants à plusieurs niveaux d'erreur ; positifs et négatifs pour chaque défaut | PLINK 1.9 + Python, `R1` | `E-MES`, `E-ID`, `E-LOWN`, `E-TOOL` | distributions des métriques, sensibilité/spécificité des décisions, faux rejets et non-évaluables avec dénominateurs |
| `06` panel de parenté | PLINK pré-QC, métriques variants TSV, régions exclues configurées | ensemble admissible, MAF, pruning et LD résiduel vrais | 22 autosomes, densités et minima par chromosome fixés en 7.4 | PLINK 1.9 + Python, `R1` | `E-MES`, `E-ALL`, `E-LOWN`, concentration chromosomique, `E-TOOL` | entonnoir de sélection, densité par chromosome, distribution du LD résiduel et fidélité de l'ensemble retenu |
| `07` parenté | panel PLINK, registre, métriques QC | pedigree, coefficients attendus, degrés et composantes indépendantes | trios, duplicats, degrés 1 à 4 et non-apparentés ; classes équilibrées selon R07 | KING + PLINK/Python, `R1` | `E-PED`, `E-ID`, apparentement cryptique, parent absent, `E-TOOL` | matrice de confusion des degrés, estimé contre vérité avec biais/IC, concordance pedigree et taille des composantes |
| `08` structure | panel indépendant, proposition d'indépendance, registre | strates, mélanges, proportions simulées et atypiques injectés | au moins deux strates, mélanges et apparentés projetés ; mêmes sujets entre densités | PLINK + PCA Python, `R1` | `E-PED`, `E-MES`, `E-LOWN`, fuite d'apparentés dans l'ajustement | valeurs propres, projections sur axes figés, distances d'atypiques et stabilité par densité/erreur ; pas d'étiquette ethnique |
| `09` gel des cohortes | registre, génotypes cibles, QC, propositions 07/08 et revues signées | appartenance attendue, composante familiale et unité indépendante | tous les sujets et toutes les composantes ; cohortes requises non vides, optionnelles explicitement nulles | Python, `R0` | `E-ID`, `E-PED`, `E-LEAK`, revue incomplète ou empreinte périmée | table attendu/observé des inclusions, attrition par règle et contrôle de fuite familiale ; règles logiques visibles |
| `10` QC final | jeux PLINK, cohortes gelées et variante cible | défauts par cohorte, décisions attendues et exception cible attendue | mêmes cohortes sous défauts stratifiés ; positifs/négatifs par politique de cohorte | PLINK 1.9 + Python, `R1` | `E-MES`, `E-ALL`, `E-ID`, cohorte vide, précision HWE, `E-TOOL` | attrition sujet/variant par cohorte, confusion des décisions, distributions MAF/HWE/manque et exceptions cible |

### 6.3 Besoins détaillés 10A à 19

| Étape | Entrées aveugles et formats | Vérité privée | Unité et taille requises | Outils / ressources | Erreurs attendues | Données graphiques obligatoires |
|---|---|---|---|---|---|---|
| `10A` univers autosomal | PLINK pré-QC et liste `target_chromosome_all_qc` | sujets, variants, ordre et complétude par autosome | mêmes sujets, génotypes complets et polymorphes sur exactement 22 autosomes | PLINK 1.9 + Python, `R1` | `E-ID`, `E-MES`, chromosome absent, doublon de variant, `E-LOWN` | nombre et densité de marqueurs par chromosome, complétude des sujets et motif des retraits |
| `11` région cible | PLINK final cible, métadonnées YAML, carte TSV/catalogue JSON | position/allèles cibles, carte vraie, fenêtre et variants attendus | plusieurs cibles autosomiques, fenêtres faisables et `INFEASIBLE_WINDOW`, densités appariées | PLINK + interpolation Python, `R1` | `E-MAP`, `E-ALL`, `E-CEN`, `E-LOWN`, `E-TOOL` | courbe bp–cM avec ancres, limites demandées/obtenues, densité locale et cas sans extrapolation possible |
| `12` phasage régional | BED/BIM/FAM régional, carte, pedigree, référence VCF/BCF et catalogues JSON | deux haplotypes vrais, recombinaisons, phase cible et copie porteuse | trios et familles sous densité, fréquence, manque et pedigree factoriels ; référence disjointe | PLINK, bcftools, SHAPEIT5, `R2` | `E-PHA`, `E-ALL`, `E-PED`, `E-MAP`, `E-LEAK`, `E-TOOL` | taux de commutation et erreur ponctuelle par MAF/densité/pedigree, confusion H1/H2/BOTH et régions non fiables |
| `12A` phasage autosomal | panel 10A, 22 cartes et 22 références indexées | phase vraie autosomale et ordre des haplotypes | mêmes sujets sur 22 autosomes, densité minimale par chromosome, références disjointes | PLINK, bcftools, SHAPEIT5, cache contrôlé, `R2` | `E-PHA`, `E-MAP`, `E-ID`, chromosome manquant, `E-LEAK`, `E-TOOL` | erreurs de phase et variants appariés/rejetés par chromosome, fréquence et densité ; coûts par chromosome |
| `13` haplotype fondateur | BCF régional, copie porteuse, cohortes, carte | segment fondateur, bornes, familles positives et partages fortuits | scénarios positifs et nuls, plusieurs cibles, âges et `N_F` = 3/5/10 | Python IBS centré cible, `R0` | `E-SEG`, `E-PHA`, `E-CEN`, `E-LOWN` | sensibilité/spécificité, erreurs gauche/droite en bp/cM, longueur estimée/vraie, exclusions et faux positifs |
| `14` datation | segments TSV et scénarios de censure | âge, recombinaisons, longueurs biologiques et mécanisme de censure | âges 5–150, `N_F` = 3/5/10, répétitions indépendantes et corrélées selon R07 | Python/SciPy Gamma, `R0` | `E-CEN`, `E-SEG`, dépendance familiale, `E-LOWN` | estimé contre vrai avec IC, biais, MAE, RMSE et couverture par âge/familles/censure ; leave-one-family-out |
| `15` LD local | PLINK régional, cohortes et carte commune | fréquences haplotypiques, `r²` génotypique/haplotypique et `D′` vrais | mêmes sujets sous densité, parenté et manque factoriels ; paires par classes bp/cM | PLINK 1.9 + Python d'audit, `R1` | `E-MES`, `E-PED`, faible MAF, paire non évaluable, `E-LOWN` | estimé contre vérité, erreur de `r²`/`D′`, décroissance par distance et effet séparé densité/parenté |
| `16` ROH | jeux PLINK par cohorte et chromosome cible | segments ROH, longueurs, fardeau et scénarios nuls | longueurs/niveaux d'autozygotie variés, fonds sans injection, densités appariées | PLINK 1.9 + Python, `R1` | `E-SEG`, fragmentation/fusion, `E-MES`, dénominateur absent, `E-LOWN` | sensibilité/faux positifs par longueur, erreurs de bornes/longueur, fardeau estimé/vrai et cible dans ROH |
| `16A` positionnement de référence | panel indépendant, BCF régional, références et métadonnées publiques | strates et mélanges globaux, haplotypes locaux et populations repères | références disjointes, plusieurs strates/mélanges et mêmes sujets sous erreurs appariées | PLINK, bcftools, PCA Python, cache contrôlé, `R2` | `E-ALL`, `E-PHA`, `E-LEAK`, référence apparentée, modèle influencé par l'étude, `E-TOOL` | axes/centroïdes ajustés sur référence, projections et distances/stabilité globales/locales ; aucune proportion inférée sans estimateur dédié |
| `16B` enrichissement | segments 13, BCF local, cohortes, ascendance et référence harmonisée | classe fondatrice/nulle, longueur partagée vraie et fond fortuit | familles indépendantes, scénarios positifs/nuls, fonds internes/externes, `N_DRAW` fixé à R07 | Python, tirages déterministes, `R1` | `E-SEG`, dépendance familiale, fond non évaluable, `E-LOWN`, `E-LEAK` | fonctions de survie nulles et valeur observée, puissance/FPR avec IC, taux d'évaluabilité et contributions gauche/droite |
| `16C` IBD explicite | 22 BCF phasés, cartes PLINK, phase régionale et cohortes | segments IBD, bornes, haplotypes et paires négatives | longueurs 1/1,5/2/3/5 cM, 500 paires positives par longueur et 3 000 paires négatives par partition | Hap-IBD + Refined IBD + Java/bcftools, `R2` | `E-SEG`, `E-PHA`, `E-MAP`, désaccord outils, calibration/fond insuffisants, `E-TOOL` | pistes de segments par outil, sensibilité/spécificité, erreur de bornes par longueur/densité et concordance des outils |
| `17` sensibilités | registre TSV de runs/vues appariés et manifests signés | facteur unique modifié, primaire et direction attendue | mêmes sujets, vérités et graines hors facteur ; répétitions inter-runs pour la stochasticité | Python, consolidation sans relance, `R0` | `E-INT`, scénario multifactoriel, appariement rompu, domaine manquant | deltas appariés avec IC, graphique en forêt par facteur/domaine, convergence et non-évaluables visibles |
| `18` visualisations | tables closes TSV/JSON et registre de valeurs autorisées | valeurs, unités, effectifs, statuts et légendes attendus | au moins un cas `PASS`, `FAIL`, `NOT_EVALUABLE` et intégrité bloquée par domaine | rendu SVG déterministe, `R3` | `E-REN`, `E-INT`, unité ou effectif discordant, donnée manquante transformée en zéro | figures explicatives par mesure, incertitudes/effectifs, provenance JSON et test pixel/structure ; aucun recalcul |
| `19` synthèse finale | Markdown clos des 24 identifiants amont (`00` à `18`, subdivisions incluses), index de figures et registre de décisions | faits autorisés, verdicts, limites, provenance et incompatibilités attendus | couverture de tous les résultats amont et de chaque statut ; aucune donnée individuelle ; le Markdown 19 est produit après l'assemblage | assemblage déterministe puis revue humaine, `R3` | `E-REN`, fait orphelin, verdict modifié, limite omise, donnée sensible ou source non autorisée | matrice de traçabilité faits–sources–figures–décisions et contrôle de complétude ; HTML/PDF seulement à R14 |

Une étape marquée `NOT_EVALUABLE` pour un besoin absent ne pourra pas être
requalifiée en `PASS`. Les intervalles et incertitudes utiliseront l'unité
indépendante appropriée. Les graphiques seront construits uniquement depuis les
tables sources empreintées décrites ci-dessus.

### 6.4 Fermeture producteur–artefact–consommateur

Cette table ferme la porte R03 à la granularité des étapes. R04 la détaille au
niveau des artefacts et vérifie les cycles.

| Étape | Producteur | Sorties scientifiques principales | Consommateurs ou justification terminale |
|---|---|---|---|
| `00` | orchestrateur/bootstrap | configuration résolue, environnement, manifest et journal | toutes les étapes, reprise et évaluateur |
| `01` | `validate_sources` | `source_inventory.tsv`, `source_qc.tsv` | `02`, audit de provenance |
| `02` | `build_sample_registry` | `samples.master.tsv`, `sample_registry_review.tsv` | `03–05`, `07–10`, `12–13`, `16`, `16A–16C` |
| `03` | `convert_acpa` | deux jeux PLINK, descripteurs, alignement et audit variants | `04`, `05` |
| `04` | `prepare_target_variant_dataset` | PLINK cible, audit génotypes, Mendel et rapport d'injection | `09`, `10`, contrôle cible aval |
| `05` | `qc_preliminary` | PLINK pré-QC, métriques individus/variants/lots et alertes | `06`, `07`, `10`, `10A` |
| `06` | `build_kinship_panel` | panel PLINK, variants prunés, couverture et LD résiduel | `07`, `08`, `16A` |
| `07` | `infer_kinship` | paires, degrés, concordance et proposition indépendante | `08`, `09`, `19` |
| `08` | `analyze_population_structure` | scores, valeurs propres, atypiques et loadings | `09`, `18` |
| `09` | `freeze_cohorts` | cohortes gelées, décisions, résumés et fichiers keep | `10`, `10A`, `13`, `15–17` |
| `10` | `qc_final` | trois jeux PLINK finaux et métriques finales | `10A`, `11`, `15`, `16` |
| `10A` | `prepare_autosomal_phasing_panel` | PLINK autosomal complet et résumé par chromosome | `12A` |
| `11` | `prepare_target_region` | PLINK régional, carte cible et manifest de phasage | `12`, `13`, `15`, `16`, `16B`, `16C` |
| `12` | `phase_target_region` | BCF/index, haplotypes porteurs, QC et harmonisation | `13`, `16A`, `16B`, `16C` |
| `12A` | `phase_autosomal_panel` | 22 BCF/index/cartes IBD et manifest autosomal | `16C` |
| `13` | `infer_founder_haplotype` | segments, consensus, matrice de partage et résumé | `14`, `16B`, `16C`, `18`, `19` |
| `14` | `estimate_variant_age` | entrées, estimations, scénarios et résumé Gamma | `18`, `19` |
| `15` | `analyze_local_ld` | paires, résumés et sorties PLINK natives | `18`, contexte terminal de `19` |
| `16` | `analyze_roh` | segments, fardeau, cible dans ROH et résumé | `16A`, `18`, `19` |
| `16A` | `analyze_reference_ancestry` | scores, modèle, audits variants, centroïdes et résumé | `16B`, `17`, `18`, `19` ; dépendance de contrôle pour `16C` |
| `16B` | `evaluate_founder_haplotype_enrichment` | unités, bornes, tirages nuls, statistique et résumé | `18`, `19` |
| `16C` | `call_explicit_ibd` | segments, paires, concordance, calibration, fond et résumé | `18`, `19` |
| `17` | `run_sensitivity_analyses` | comparaisons, stabilité, convergence et résumé | `18`, `19` |
| `18` | `build_visualizations` | SVG, provenance et index de figures | `19` ; galeries HTML/PDF interdites avant R14 |
| `19` | assemblage final contrôlé | faits, traçabilité, revue et rapport final | décision humaine R14 et archive de qualification |

### 6.5 Écarts de contrat à résoudre prospectivement

R03 identifie les écarts suivants sans modifier le logiciel :

1. `00` ne possède pas encore les contrats de release maître, vérité privée et
   vues dérivées définis à la section 9.
2. `16A` publie un positionnement PCA relatif. Il ne publie pas une proportion
   d'ascendance validée ; toute métrique de proportion restera
   `NOT_EVALUABLE` sauf création prospective d'une méthode distincte.
3. `17` doit couvrir les domaines `16B` et `16C` dans la nouvelle campagne,
   alors que le contrat documentaire actuel ne les consolide pas encore.
4. `18` produit actuellement une galerie HTML/PDF. La campagne unifiée
   n'autorise avant R14 que les SVG, provenances et tables sources ; le contrat
   devra être séparé à R09 sans réutiliser un rendu historique.
5. `19` produit actuellement un brouillon HTML. Pour la calibration, R10–R13
   évalueront uniquement les faits structurés, la traçabilité et la revue ; le
   rapport HTML/PDF de qualification ne sera matérialisé qu'à R14.
6. Les contrats et tests de l'ancienne calibration ne sont pas des producteurs
   autorisés. Les nouveaux schémas seront créés à R09 depuis cette matrice et le
   graphe R04, jamais restaurés depuis l'ancienne campagne.

Ces écarts ne bloquent pas R04–R08, qui sont des étapes de conception. Ils
bloquent toute création de contrat à R09 tant qu'une résolution explicite n'a
pas été adoptée.

### 6.6 Conventions du graphe R04

Le graphe ci-dessous est normatif pour la conception, mais ne crée encore
aucun fichier ni schéma. Les chemins logiques sont relatifs à la future racine
`calibration_workspace/releases/<release_id>/`. Une accolade désigne
exhaustivement une famille d'artefacts : chaque membre sera déclaré avec son
chemin, son schéma, sa cardinalité et son SHA-256 dans le manifest producteur.
Une sortie native ou dynamique n'est autorisée que si son registre fermé la
nomme ; un glob découvert après exécution n'est jamais une dépendance.

Chaque exécution d'étape `XX` produit aussi l'enveloppe commune
`runs/<run_id>/XX/{stage_inputs.json,command.json,stdout.log,stderr.log,stage_outputs.json,audit.json,checksums.sha256}`.
Cette enveloppe est consommée par la clôture du run, puis par l'évaluateur de
`XX`; elle n'est pas répétée dans les tables. Les artefacts racine du run
`config.resolved.yaml`, `environment.json`, `manifest.json` et `events.jsonl`
sont produits par `00`, consommés par toutes les étapes et clos avec le run.

Les producteurs sont de sept types : `GEN` génère la vérité maître privée ;
`DERIVE` construit une vue aveugle depuis cette vérité ; `PIPE-XX` est l'étape
logicielle aveugle ; `HUMAN` signe une décision de revue ; `CLOSE` fige un run
sans ouvrir la vérité ; `EVAL-XX` compare un run clos à la seule vérité
autorisée ; `RENDER` rend des SVG depuis des tables d'évaluation closes. Une
vérité privée ne peut jamais être consommée par `PIPE-*`, `CLOSE` ou `RENDER`.

### 6.7 Racines de gouvernance, vérités et vues

#### 6.7.1 Racines obligatoires

| Artefact logique | Producteur | Consommateurs autorisés |
|---|---|---|
| `governance/panel_release_manifest.json` | `GEN` puis signature humaine | `DERIVE`, `00`, `CLOSE`, tous les `EVAL-*`, R14 |
| `governance/source_catalog.json` | `GEN` | `DERIVE`, `01`, `11`, `12`, `12A`, `16A`, `CLOSE` |
| `governance/seed_ledger.tsv` | `GEN` | `DERIVE`, `EVAL-00`, `EVAL-17`, audit R14 |
| `governance/subject_registry.tsv` | `GEN` | `DERIVE` et évaluateurs explicitement autorisés ; jamais le pipeline |
| `private_truth/master_truth_manifest.json` | `GEN` | routeur des `EVAL-*` uniquement |
| `views/panel_view_manifest.json` | `DERIVE` | `00`, étapes déclarées, `CLOSE`, tous les `EVAL-*` |
| `blinded_inputs/blinded_input_inventory.json` | `DERIVE` | `00`, `01`, `CLOSE`, `EVAL-00`, `EVAL-01` |
| `governance/readiness_receipt.json` | contrôle humain préalable | `00`, `CLOSE`, R14 |

Le reçu de préparation ne peut être positif que si les sept artefacts qui le
précèdent existent, sont empreintés et compatibles. Il n'est donc pas une
source indépendante ni un raccourci vers la vérité.

#### 6.7.2 Familles de vérités privées

| Identifiant | Contenu privé exhaustif | Évaluateurs consommateurs |
|---|---|---|
| `T-SOURCE` | validité, assemblage, colonnes, couverture et empreintes attendues | `EVAL-01` |
| `T-SUBJECT` | identités, rôles, pedigree, familles, lots et partitions maîtres | `EVAL-02`, `EVAL-07`, `EVAL-09` |
| `T-GENOTYPE` | génotypes autosomiques denses, ordre et allèles | `EVAL-03`, `EVAL-05`, `EVAL-10`, `EVAL-10A` |
| `T-TARGET` | variants cibles, copies porteuses et transmissions | `EVAL-04`, `EVAL-09`, `EVAL-10`, `EVAL-12` |
| `T-MEASURE` | manque, erreurs, duplications, hétérozygotie et lots injectés | `EVAL-03`, `EVAL-05`, `EVAL-10`, `EVAL-17` |
| `T-KINSHIP` | paires, degrés, composantes et indépendance vrais | `EVAL-06`, `EVAL-07`, `EVAL-08`, `EVAL-09` |
| `T-POPULATION` | strates, mélanges, proportions et atypiques simulés | `EVAL-08`, `EVAL-16A` |
| `T-PHASE` | haplotypes, recombinaisons, phase et copie porteuse | `EVAL-11`, `EVAL-12`, `EVAL-12A`, `EVAL-13` |
| `T-FOUNDER` | segments fondateurs, familles, âges et censure | `EVAL-13`, `EVAL-14`, `EVAL-16B` |
| `T-LD` | fréquences haplotypiques, paires, `r²` et `D′` | `EVAL-15` |
| `T-ROH` | segments ROH, bornes, fardeaux et fonds nuls | `EVAL-16` |
| `T-IBD` | segments, haplotypes, bornes et paires négatives | `EVAL-16C` |
| `T-SCENARIO` | facteur modifié, appariement, métrique primaire et direction attendue | `EVAL-17` |
| `T-RENDER` | unités, légendes, statuts, faits et décisions attendus | `EVAL-18`, `EVAL-19` |

Les quatorze familles sont toutes produites par `GEN` et chacune possède au
moins un consommateur. Le `master_truth_manifest` accorde un chemin précis à
un évaluateur ; il ne livre jamais l'ensemble des vérités par défaut.

#### 6.7.3 Vues aveugles dérivées

| Vue | Dérivation aveugle | Étapes consommatrices |
|---|---|---|
| `V-SOURCE` | exports synthétiques et métadonnées publiques permises | `01`, `02`, `03` |
| `V-TARGET` | métadonnées moléculaires et appels cibles explicitement autorisés | `04`, `09`, `10` |
| `V-DENSITY` | intersections/unions et niveaux de densité appariés | `03`, `05`, `06`, `08`, `10`, `15`, `16`, `17` |
| `V-MEASUREMENT` | niveaux d'erreur et de manque sans étiquette de vérité | `03`, `05`, `10`, `17` |
| `V-PEDIGREE` | pedigree déclaré, incomplet ou altéré selon le scénario aveugle | `02`, `07`, `09`, `12`, `17` |
| `V-REGIONAL` | fenêtre cible, carte publique et génotypes régionaux | `11`, `12`, `13`, `14`, `15`, `16`, `16A`, `16B`, `16C` |
| `V-AUTOSOMAL` | mêmes sujets et variants sur les 22 autosomes | `05`, `06`, `07`, `08`, `10A`, `12A`, `16C` |
| `V-PHASING-REFERENCE` | références synthétiques disjointes, propres à la partition et harmonisées | `12`, `12A`, `16A`, `16B` |
| `V-FOUNDER-NULL` | fonds internes/externes et tirages aveugles préspecifiés | `13`, `16B` |
| `V-ROH-IBD` | scénarios positifs/nuls ROH et IBD sans labels privés | `16`, `16C` |
| `V-SENSITIVITY` | couples de runs ne différant que par un facteur | `17` |

Chaque vue est un descendant du `panel_release_manifest`, conserve les
identifiants pseudonymisés et l'empreinte de ses parents, et ne peut ni créer
un sujet ni changer une vérité. Les filtres, transformations et pertes
attendues sont fixés en 5.16–5.22.

### 6.8 Graphe détaillé des artefacts du pipeline

Dans les deux tables suivantes, `EVAL-XX` consomme toutes les sorties de
`PIPE-XX`, y compris lorsqu'un consommateur aval est aussi indiqué. Cette règle
donne une justification terminale aux audits, résumés et sorties natives sans
les transformer artificiellement en entrées scientifiques aval. Les noms sont
ceux du logiciel actif observé à R04 ; leur futur contrat sera recréé à R09.

| Producteur | Artefacts produits, groupés sans omission | Autres consommateurs déclarés |
|---|---|---|
| `PIPE-00` | `initialization_summary.json` et artefacts racine du run définis en 6.6 | toutes les étapes via le manifest ; `CLOSE` |
| `PIPE-01` | `source_inventory.tsv`, `source_qc.tsv` | `02`, `CLOSE` |
| `PIPE-02` | `samples.master.tsv`, `sample_registry_review.tsv` | `03`, `04`, `05`, `07`, `08`, `09`, `10`, `12`, `12A`, `13`, `16`, `16A`, `16B`, `16C` ; `HUMAN` pour la revue |
| `PIPE-03` | `genomewide_base.{bed,bim,fam}`, `genomewide_base.dataset.json`, `target_chromosome_base.{bed,bim,fam}`, `target_chromosome_base.dataset.json`, `sample_alignment.tsv`, `acpa_variant_audit.tsv`, `acpa_sample_chromosome_qc.tsv`, `acpa_conversion_report.json` | jeu génomique vers `05` ; jeu cible vers `04` ; audits vers `CLOSE` |
| `PIPE-04` | `target_variant.{bed,bim,fam}`, `target_variant.dataset.json`, `target_genotype_audit.tsv`, `target_variant_mendel.tsv`, rapport Mendel natif enregistré, `target_variant_injection_report.json` | `09`, `10`, `12`, `16` ; rapport natif via son registre |
| `PIPE-05` | `genomewide_pre_qc.{bed,bim,fam}`, `genomewide_pre_qc.dataset.json`, `qc_individual_metrics.tsv`, `qc_variant_metrics.tsv`, `qc_batch_summary.tsv`, `qc_alerts.tsv`, `qc_preliminary_report.json` | `06`, `07`, `10`, `10A` |
| `PIPE-06` | `kinship_panel.{bed,bim,fam}`, `kinship_panel.dataset.json`, `pruned_variants.tsv`, `panel_coverage.tsv`, `panel_residual_ld_bins.tsv`, `kinship_panel_report.json` | `07`, `08`, `16A` |
| `PIPE-07` | `kinship_pairs.tsv`, `kinship_degree_summary.tsv`, `pedigree_concordance.tsv`, `independent_set_proposal.tsv`, `kinship_report.json` | `08`, `09`, `17` |
| `PIPE-08` | `population_scores.tsv`, `population_eigenvalues.tsv`, `population_outliers.tsv`, `population_variant_loadings.tsv`, `population_structure_report.json` | `09`, `17` ; les figures passent par `EVAL-08` |
| `HUMAN-09` | décisions signées et empreintées de revue de parenté et de structure | `PIPE-09`, `CLOSE`, `EVAL-09` |
| `PIPE-09` | `cohorts_frozen.tsv`, `cohort_summary.tsv`, `cohort_decisions.tsv`, `cohort_freeze_report.json`, `{controls_unrelated,target_carriers_all,target_carriers_independent,family_noncarriers,affected_noncarriers,target_chromosome_all_qc,genomewide_structure_reference}.keep` | `10`, `10A`, `13`, `15`, `16`, `16A`, `16B`, `16C`, `17` selon cohorte |
| `PIPE-10` | `{controls_unrelated_qc,target_carriers_independent_qc,target_chromosome_all_qc}.{bed,bim,fam}`, les trois fichiers `.dataset.json`, `sample_qc_final.tsv`, `variant_qc_final.tsv`, `batch_qc_final.tsv`, `qc_final_report.json` | `10A`, `11`, `15`, `16`, `17` |

| Producteur | Artefacts produits, groupés sans omission | Autres consommateurs déclarés |
|---|---|---|
| `PIPE-10A` | `autosomal_phasing_panel.{bed,bim,fam}`, `autosomal_phasing_panel.dataset.json`, `autosomal_phasing_panel_summary.json` | `12A`, `17` |
| `PIPE-11` | `target_region.{bed,bim,fam}`, `target_region.dataset.json`, `target_genetic_map.tsv`, `phasing_input_manifest.json`, `target_region_report.json` | `12`, `13`, `15`, `16`, `16B`, `16C`, `17` |
| `PIPE-12` | artefacts logiques `{reference_window_vcf,reference_window_index,reference_window_manifest,harmonized_reference_vcf,harmonized_reference_index,variant_harmonization,reference_harmonization_manifest,shapeit5_study_vcf,shapeit5_study_index,shapeit5_genetic_map,shapeit5_pedigree,shapeit5_variant_selection,shapeit5_sample_mapping,shapeit5_mendel_exclusions,shapeit5_inputs_manifest,shapeit5_common_bcf,shapeit5_common_index,shapeit5_final_bcf,shapeit5_final_index,shapeit5_common_log,shapeit5_rare_log,carrier_haplotypes,phasing_transmissions,phasing_unreliable_regions,shapeit5_phasing_manifest,phasing_qc,carrier_haplotype_summary,phase_target_region_summary}` | `13`, `16A`, `16B`, `16C`, `17`; chaque chemin physique est nommé dans le manifest producteur |
| `PIPE-12A` | pour chacun des chromosomes `1..22` : BCF, index, carte IBD et journal ; `autosomal_phasing_manifest.json` | `16C`, `17` ; exactement 88 fichiers chromosomiques enregistrés |
| `PIPE-13` | `founder_segments.tsv`, `founder_consensus.tsv`, `founder_sharing_matrix.tsv`, `founder_discordances.tsv`, `founder_analysis_summary.json` | `14`, `16B`, `16C`, `17` ; les figures passent par `EVAL-13` |
| `PIPE-14` | `variant_age_input.tsv`, `variant_age_estimates.tsv`, `variant_age_scenarios.tsv`, `variant_age_summary.json` | `17` ; les figures passent par `EVAL-14` |
| `PIPE-15` | `local_ld_variants.tsv`, `local_ld_pairs.tsv`, `local_ld_summary.tsv`, `local_ld_analysis_summary.json`, sorties natives `local_ld_native_*` inscrites au registre | `17` ; les figures passent par `EVAL-15` |
| `PIPE-16` | `roh_segments.tsv`, `roh_burden.tsv`, `target_in_roh.tsv`, `roh_cohort_summary.tsv`, `roh_analysis_summary.json`, sorties natives `roh_native_*` inscrites au registre | `16A`, `17` ; les figures passent par `EVAL-16` |
| `PIPE-16A` | `ancestry_scores.tsv`, `ancestry_eigenvalues.tsv`, `ancestry_variant_loadings.tsv`, `ancestry_variant_audit.tsv`, `ancestry_population_centroids.tsv`, `reference_ancestry_summary.json` | `16B`, `16C` comme contrôle obligatoire, `17` ; les figures passent par `EVAL-16A` |
| `PIPE-16B` | `founder_haplotype_units.tsv`, `founder_haplotype_boundaries.tsv`, `founder_haplotype_family_consistency.tsv`, `founder_haplotype_null_draws.tsv.gz`, `founder_haplotype_variant_audit.tsv`, `founder_haplotype_enrichment_summary.tsv`, `founder_haplotype_enrichment_summary.json` | `17` ; les figures passent par `EVAL-16B` |
| `PIPE-16C` | `explicit_ibd_segments.tsv`, `ibd_marker_audit.tsv`, `target_ibd_pair_matrix.tsv`, `ibd_tool_concordance.tsv`, `ibd_calibration.tsv`, `ibd_control_frequency.tsv`, `target_chromosome_phase_alignment.tsv`, `explicit_ibd_feasibility.json`, `explicit_ibd_summary.json`, sorties et journaux bruts inscrits au registre | `17` ; les figures passent par `EVAL-16C` |
| `PIPE-17` | `sensitivity_comparisons.tsv`, `sensitivity_stability.tsv`, `sensitivity_analysis_summary.json`, `population_convergence.json` | `EVAL-17`, puis `RENDER-18` uniquement via les tables d'évaluation closes |
| `RENDER-18` | pour chacun des onze domaines `{population_structure,founder_ibs,variant_age,local_ld,roh,reference_ancestry_global,reference_ancestry_local,founder_haplotype_enrichment,explicit_ibd,population_convergence,sensitivity}` : SVG et JSON de provenance ; `figure_index.json`, `visualization_completeness.json` | `EVAL-18`, `19` ; aucune galerie HTML/PDF |
| `PIPE-19` | `interpretation_facts.json`, matrice de traçabilité faits–sources–figures–décisions et proposition de décision structurée | `EVAL-19`, revue humaine finale et R14 ; aucun brouillon HTML |

Les Markdown `calibration/resultats/ETAPE_<id>.md` sont produits par les
évaluateurs, non par les étapes aveugles. Pour `19`, l'assemblage déterministe
prépare les faits et la traçabilité, `EVAL-19` produit le Markdown, puis la
revue humaine le signe sans modifier les faits évalués.
Les sorties actuelles interdites avant R14 ne figurent donc pas comme arêtes
valides du nouveau graphe.

### 6.9 Clôture, évaluation et sensibilités

Pour chaque run, `CLOSE` vérifie les enveloppes communes, les sorties déclarées
et leurs empreintes, puis produit
`runs/<run_id>/pipeline_run_manifest.json`. Un manifest clos est immuable : une
correction crée un nouvel identifiant de run. Il contient les chemins et
empreintes, jamais la vérité privée.

Chaque `EVAL-XX` consomme exactement : le manifest du run clos, les sorties
`XX` qu'il référence, le sous-ensemble de vérité autorisé en 6.7.2 et les
paramètres de métriques qui seront figés à R08. Il produit :

- `evaluations/XX/ETAPE_XX.metrics.tsv`, source numérique des figures ;
- `evaluations/XX/ETAPE_XX.evaluation.json`, statut, effectifs indépendants,
  incertitudes, exclusions et empreintes ;
- `calibration/resultats/ETAPE_XX.md`, résultat lisible clos.

Pour `17`, `HUMAN-17` produit avant calcul
`governance/sensitivity_run_registry.tsv`, signé et empreinté. Chaque ligne
référence deux `pipeline_run_manifest.json` clos, la même release, les mêmes
sujets, vérités et graines, un unique facteur différent et un domaine parmi
`03`, `05`, `07`, `08`, `10`, `10A`, `11`, `12`, `12A`, `13`, `14`, `15`,
`16`, `16A`, `16B` ou `16C`. `PIPE-17` ne recherche aucun run par répertoire,
n'accepte aucun run courant et ne relance aucune étape. `EVAL-17` ouvre ensuite
`T-MEASURE` et `T-SCENARIO` pour vérifier l'appariement et les directions.

`RENDER-18` consomme uniquement les fichiers `*.metrics.tsv` et
`*.evaluation.json` clos de `EVAL-00` à `EVAL-17`, plus un registre explicite
des valeurs autorisées. Il ne lit ni vérité, ni table individuelle, ni sortie
native. `PIPE-19` consomme les 24 Markdown amont, l'index des figures, leurs
provenances et les décisions closes. R14 est seul autorisé à rendre le rapport
HTML/PDF depuis ces sources sans recalcul scientifique.

### 6.10 Écarts d'implémentation révélés par le graphe

| Écart | Constat R04 | Résolution prospective obligatoire |
|---|---|---|
| `GAP-R04-01` | `00` ne consomme pas encore les huit racines de gouvernance | créer les contrats et la porte de préparation à R09, après R05–R08 |
| `GAP-R04-02` | aucun routeur de vérités ni évaluateur séparé n'existe dans la nouvelle campagne | matérialiser `master_truth_manifest` et `EVAL-*` à R09–R10 |
| `GAP-R04-03` | les décisions manuelles de `09` sont partiellement portées par des paramètres | exiger des décisions signées, empreintées et liées au run |
| `GAP-R04-04` | la dépendance déclarée de `16C` envers `16A` ne nomme actuellement aucun artefact requis | rendre `reference_ancestry_summary.json` ou un reçu de contrôle `16A` obligatoire |
| `GAP-R04-05` | `17` résout des runs via un registre incomplet et ne couvre pas encore `16B`–`16C` | imposer le registre signé de 6.9 et couvrir les seize domaines listés |
| `GAP-R04-06` | certaines sorties PLINK, SHAPEIT5, ROH ou IBD ont une multiplicité dynamique | les fermer dans un manifest producteur exhaustif ; interdire la découverte par glob |
| `GAP-R04-07` | `18` matérialise actuellement une galerie HTML/PDF | séparer à R09 le rendu SVG d'évaluation du rendu final R14 |
| `GAP-R04-08` | `19` lit une sélection de sorties et produit un brouillon HTML | assembler les 24 Markdown et la traçabilité ; réserver HTML/PDF à R14 |

Ces écarts décrivent du travail futur et n'autorisent aucune modification du
logiciel pendant R04. Ils bloquent une exécution de calibration tant que R09
ne les a pas résolus, mais ils ne laissent aucune arête implicite dans le
graphe normatif : chaque résolution attendue, son producteur et son
consommateur sont désormais nommés.

### 6.11 Preuve de fermeture et d'acyclicité

Un ordre topologique valide est :

`GEN < DERIVE < HUMAN-readiness < 00 < 01 < 02 < 03 < {04,05} < 06 < 07 < 08 < HUMAN-09 < 09 < 10 < {10A,11} < {12,12A,15,16} < {13,16A} < {14,16B,16C} < CLOSE-runs < HUMAN-17 < 17 < CLOSE-17 < EVAL-00..17 < 18 < CLOSE-18 < EVAL-18 < 19 < CLOSE-19 < EVAL-19 < HUMAN-final < R14`.

Les éléments entre accolades sont des branches parallèles : leur ordre interne
n'est pas une arête. `16A` reçoit `12` et `16`; `16B` reçoit `13` et `16A`;
`16C` reçoit `12`, `12A`, `13` et le contrôle `16A`. `17` n'accepte que des
runs antérieurs clos. Les évaluateurs ne publient aucune vérité vers le
pipeline, `18` ne publie que des vues de résultats, `19` n'alimente aucune
étape scientifique et R14 est terminal. Aucune arête ne remonte donc cet ordre.

La porte R04 est satisfaite au niveau de conception :

- les 25 identifiants ont chacun un producteur, des sorties fermées et un
  évaluateur terminal ;
- les quatorze familles de vérité ont toutes un producteur et au moins un
  évaluateur consommateur, sans accès du pipeline ;
- les onze vues possèdent un parent de release et des consommateurs nommés ;
- les artefacts communs, décisions humaines, runs appariés, sorties dynamiques,
  figures et Markdown possèdent une règle de fermeture explicite ;
- l'ordre topologique ci-dessus exclut une dépendance circulaire.

La preuve devra être matérialisée en contrats testables à R09 ; elle n'autorise
ni génération de panel ni exécution scientifique avant les portes R05–R08.

## 7. Grille scientifique et dimensionnement R07

### 7.1 Unités de dimensionnement et mondes de scénario

R07 distingue trois quantités qui ne peuvent pas être additionnées ou
substituées silencieusement : `N_S_ACTIVE`, nombre d'identités biologiques
présentes dans un run ; `N_UNIT`, nombre d'unités indépendantes pour une
métrique ; `N_REP`, nombre de mondes de scénario générés par des graines
distinctes. Un variant, une paire de lignes ou plusieurs membres d'une famille
ne sont pas des répétitions indépendantes lorsque leur composante est commune.

Un `scenario_world_id` est une réalisation indépendante, contenue dans la
même release maître, avec ses propres `isolation_unit_id`, identités et graines.
Les mondes spécialisés de datation, ROH ou IBD ne sont donc ni des panels
parallèles ni des ajouts post hoc. Ils réutilisent les mécanismes R05, les
partitions R06 et les transformations R06, mais aucune identité biologique ne
passe d'un monde à l'autre. Les analyses appariées réutilisent une même unité
dans deux contrefactuels préscrits ; cette unité compte une seule fois.

Le run transversal mesure l'intégration des étapes `00`–`19`. Les mondes
spécialisés apportent la précision statistique d'un domaine sans imposer que
des milliers de répétitions légères soient simultanément phasées comme une
seule cohorte. Tout résultat agrégé conserve `scenario_world_id` et utilise une
variance groupée au niveau de l'unité indépendante pertinente.

### 7.2 Cibles de précision et justification des minima

Les calculs suivants servent uniquement au dimensionnement. Le niveau
bilatéral `0,05`, la puissance `0,90` et la correction de conception sur seize
domaines ne sont pas des seuils de qualification ; R08 fixera les métriques,
seuils et règles de multiplicité effectivement appliqués.

Pour une proportion, la demi-largeur de l'intervalle de Wilson à 95 % est
maximale autour de `p = 0,5`. Avec `N_UNIT = 500`, elle vaut `0,0437`, soit
environ 4,4 points de pourcentage. Les cellules binaires primaires utilisent
donc 500 unités. Une cellule secondaire à 250 unités a une précision moindre
et ne peut pas être promue implicitement en cellule primaire.

Pour un contrôle négatif sans événement observé, la borne supérieure exacte
unilatérale à 95 % est :

`u_0(n) = 1 - 0,05^(1/n)`.

Avec 3 000 unités négatives indépendantes, `u_0 = 0,000998`, ce qui permet
d'estimer un ordre de grandeur de faux positif rare sans préjuger du seuil R08.
Chaque domaine produisant une décision positive/négative conserve donc
3 000 unités nulles pour sa qualification primaire lorsque ce fond nul est
défini.

Pour une moyenne, 500 unités donnent sous approximation normale une
demi-largeur de `1,96 / sqrt(500) = 0,0877` écart-type. Les erreurs de bornes,
longueurs, âges et deltas appariés utilisent 500 unités primaires. Pour un
delta apparié standardisé `d = 0,20`, une puissance normale de 90 % avec une
allocation bilatérale `0,05 / 16` exige 449 paires ; 500 paires portent la
puissance de conception à environ 93,5 %. Les outils externes stochastiques
font exception : 50 graines indépendantes sont le minimum de qualification et
leurs intervalles seront présentés comme tels, sans prétendre à la précision
des 500 paires biologiques.

La variance continue inconnue sera estimée uniquement sur le pilote de
calibration. Avec `s` sur `n_p` unités et
`s_U = s × sqrt((n_p - 1) / chi2_{0,05;n_p-1})`, la règle prospective est :

`N_FINAL = max(N_MIN, ceil((z × s_U / h)^2))`,

où `h` est la précision absolue définie avant le pilote et `z` la valeur de
conception appropriée. La même règle et le même `N_FINAL` sont appliqués à la
validation encore scellée. Pour une proportion, l'inversion exacte de Wilson
remplace cette formule. Une règle de taille n'utilise jamais un résultat de la
validation.

### 7.3 Échelles, sujets et références

Une seule release de qualification contient dès l'origine des réserves
ordonnées par identifiant et graine. Le pilote ouvre uniquement le préfixe de
calibration prévu ; la validation reste inaccessible. Un `activation_rank`
est fixé avant génotype par rotation équilibrée des cellules de scénario, puis
par hachage des identifiants au sein de chaque cellule. L'activation d'une
réserve prend toujours les premiers rangs et ne dépend jamais d'un génotype ou
d'une performance observée.

| Échelle / vue | Étude active | Référence de phase | Référence d'ascendance | Référence nulle | Usage |
|---|---:|---:|---:|---:|---|
| micro transversal `RESOURCE_ONLY` | 60 | 90 | 90, soit 30 par population | 60 | intégration et erreurs déterministes, aucune performance |
| pilote de calibration | 500 | 500 | 300, soit 100 par population | 500 | variance, non-évaluables et ressources |
| qualification nominale, par partition | 4 000 | 3 000 | 1 800, soit 600 par population | 3 000 | estimation primaire et validation aveugle |
| réserve maximale, par partition | 8 000 | 6 000 | 3 600 | 6 000 | extension uniquement par la règle prospective |

Le micro-panel contient exactement quatre trios (`12` sujets), quatre fratries
de quatre (`16`), quatre composantes multigénérationnelles de cinq (`20`) et
douze sujets sans apparentement (`12`). Les duplicats sont des échantillons
techniques supplémentaires et ne gonflent pas les 60 identités. Les trois
populations synthétiques, les degrés 1–4, un positif et un négatif de chaque
mécanisme majeur sont représentés, mais aucun taux de performance n'est calculé.

Les 4 000 identités nominales de chaque partition de qualification comprennent
2 250 représentants indépendants, 250 pour chacun des trois états purs et des
six mélanges globaux. Deux cents de ces représentants ancrent des composantes
de huit sujets ; leurs 1 400 proches supplémentaires fournissent une paire
désignée de chacun des degrés 1–4 par composante et 200 trios. Les 350 identités
restantes sont des sujets non apparentés affectés avant génération aux défis
de cible, lot, QC et fonds négatifs. Quatre lots de 1 000 identités sont
équilibrés par état populationnel à une identité près. Deux cents duplicats
techniques sont ajoutés comme échantillons, jamais comme identités ou unités
indépendantes.

Les trois références restent disjointes entre elles, entre partitions et avec
l'étude. À la qualification nominale, la fréquence mineure de conception la
plus basse `0,01` représente environ 60 copies dans une référence de phase de
3 000 sujets ; toute strate réellement plus rare est rapportée séparément et
ne reçoit pas artificiellement des copies supplémentaires.

### 7.4 Niveaux numériques des facteurs

Les niveaux suivants sont fermés avant génération :

| Facteur | Niveaux R07 |
|---|---|
| cible | les cinq rôles R05 : `STANDARD`, `LOW-RECOMB`, `HIGH-RECOMB`, `LEFT-EDGE`, `RIGHT-EDGE` |
| âge fondateur | `5`, `10`, `15`, `20`, `30`, `40`, `50`, `75`, `100`, `150` générations |
| familles | sous-ensembles emboîtés des mêmes événements : `3`, `5`, `10` familles |
| fréquence ALT synthétique hors événement familial | `0,01`, `0,05`, `0,20` |
| densité autosomale observée | `D_MAX = 600 000`, `D_MID = 300 000`, `D_LOW = 100 000` variants |
| pedigree déclaré | `COMPLETE`, `PARENT_MISSING`, `NO_PEDIGREE`, `MISDECLARED_CHALLENGE` |
| populations | `POP-A`, `POP-B`, `POP-C` et les six mélanges par paires `25/75`, `50/50`, `75/25` |
| IBD vrai | `1`, `1,5`, `2`, `3`, `5 cM`, plus `NONE` |
| censure | `NONE`, `LEFT`, `RIGHT`, `DOUBLE` |

Le fond haplotypique privé comporte 1 200 000 variants autosomiques. Les trois
densités sont des sous-ensembles strictement emboîtés. Après protection de la
cible et des ancres déclarées, l'allocation suit la longueur génétique par
chromosome, avec au moins 5 000, 3 000 et 1 000 variants respectivement par
autosome pour `D_MAX`, `D_MID` et `D_LOW`. Dans une fenêtre principale
évaluable, les densités visées sont au moins 100, 50 et 20 variants par cM ;
une impossibilité est enregistrée avant génération et non compensée par un
autre chromosome choisi après observation.

Les événements de mesure numériques sont :

| Famille | Témoin, niveau bas, niveau haut |
|---|---|
| manque aléatoire par cellule | `0`, `0,01`, `0,05` |
| perte supplémentaire par sujet | `0`, `0,05`, `0,15` |
| perte supplémentaire par variant | `0`, `0,02`, `0,10` |
| retournement/dosage génotypique | `0`, `0,001`, `0,01` des appels observables |
| lot, appliqué à un des quatre lots | `0`, `+0,02`, `+0,10` de manque et `0`, `+0,001`, `+0,01` d'erreur au niveau correspondant |
| commutations de phase injectées en défi aval | `0`, `1`, `5` événements par 100 cM |
| trou de carte | `0`, `2`, `10 cM` |

`ERR-SOURCE`, `ERR-IDENTITY`, `ERR-ALLELE`, les incohérences catégorielles de
pedigree et `ERR-INTEGRITY` possèdent un témoin valide et un cas déterministe
par mécanisme de la section 5.12 ; répéter leur même fichier ne constitue pas
une puissance statistique. Dans les cellules primaires, une seule famille
d'erreur diffère du témoin. Les interactions autorisées sont limitées à celles
de la section suivante.

Les classes ROH numériques sont : `NONE` ; `SHORT`, 2 cM et 1–3 Mb ;
`MEDIUM`, 5 cM et 3–8 Mb ; `LONG`, 15 cM et 10–25 Mb ;
`MULTIPLE_BURDEN`, trois segments disjoints de 5 cM, chacun entre 3 et 8 Mb.
Les bornes bp sont obtenues par la carte et doivent satisfaire simultanément
la classe physique, la longueur cM à ±5 % et la densité requise. Deux segments
de 5 cM séparés de 0,25 cM forment le défi de fusion ; un intervalle perturbé
de 0,25 cM au centre d'un segment de 15 cM forme le défi de fragmentation. Une
cellule sans intervalle compatible est `INFEASIBLE_ROH`, jamais déplacée après
tirage.

### 7.5 Plans factoriels fermés par domaine

Un produit cartésien global serait ininterprétable et inutilement coûteux.
Les interactions biologiquement nécessaires sont complètes ci-dessous ; les
autres facteurs suivent un plan à un facteur variant autour du témoin
`STANDARD × D_MAX × COMPLETE × NONE × POP-A`. Les mêmes unités sont utilisées
pour les niveaux appariés.

| Domaine | Grille primaire | Répétitions de qualification par partition |
|---|---|---:|
| transversal `00`–`19` | une cible `STANDARD`, `D_MAX`, pedigree complet, aucune erreur ; puis chaque défi déterministe isolé | un run nominal et un run de reprise ; aucune estimation de performance sur le micro |
| QC `03`–`10` | trois densités appariées ; chaque famille de mesure aux trois niveaux ; interactions `densité × manque cellule` et `densité × lot` | 500 unités positives et 3 000 unités nulles par décision primaire |
| parenté `07` | `DUPLICATE`, degrés `1`–`4`, `UNRELATED`, sous les trois densités | 500 paires désignées indépendantes par classe ; 200 dans le run transversal et 300 dans des mondes spécialisés |
| structure `08` | neuf états populationnels sous les trois densités ; manque et lot testés séparément | 250 représentants indépendants par état dans le run nominal ; 500 paires pour chaque contraste primaire regroupé |
| phase régionale `12` | `densité 3 × fréquence ALT 3 × pedigree {COMPLETE, PARENT_MISSING, NO_PEDIGREE}` | mêmes 200 trios et 500 sujets indépendants par cellule appariée ; 50 graines d'outil |
| phase autosomale `12A` | 22 autosomes à `D_MAX` ; sensibilités `D_MID` et `D_LOW` sur trois sentinelles présélectionnées | 10 graines sur chacun des 22 autosomes à `D_MAX` ; 50 graines au total par sentinelle et par densité, soit 340 appels `D_MAX`, 150 `D_MID` et 150 `D_LOW` |
| fondateur et âge `13`–`14` | dix âges × sous-ensembles emboîtés de 3/5/10 familles | 500 événements indépendants par âge, chacun généré avec dix familles |
| LD `15` | trois densités × trois populations pures ; manque et parenté séparément | au moins 500 blocs haplotypiques indépendants par classe de distance, agrégation par bloc |
| ROH `16` | quatre classes positives × trois placements `{cible dedans, cible dehors même chromosome, autre chromosome}` × trois densités | 500 sujets-mondes par classe avec cible dedans ; 250 par autre placement ; 3 000 sujets-mondes `NONE` |
| ascendance `16A` | neuf états étude × vues globale/locale × trois densités | 250 représentants par état et référence de 600 sujets par population |
| enrichissement `16B` | âges, 3/5/10 familles, fonds interne/externe, classes positive et nulle | mêmes 500 événements par âge ; tirages nuls de 9 999 par statistique primaire |
| IBD `16C` | cinq longueurs × trois densités × deux outils ; localisation cible/hors cible/ autre autosome séparée | 500 paires positives indépendantes par longueur, 200 par longueur et classe de localisation secondaire, 3 000 paires négatives |
| sensibilités `17` | les seize domaines de 6.9, un facteur signé à la fois | 500 paires biologiques ou analytiques ; 50 paires de graines pour un outil externe stochastique |

Les trois chromosomes sentinelles de `12A` sont choisis avant génération : le
plus court, le médian inférieur et le plus long selon la longueur cM de la
carte empreintée, avec départage par numéro de chromosome. Les dix répétitions
sur 22 autosomes servent à vérifier l'hétérogénéité chromosomique ; elles ne
sont pas annoncées avec la précision d'une cellule à 500 unités.

### 7.6 Répétitions fondatrices, fonds nuls et IBD

Pour chaque âge, les 500 événements positifs possèdent dix familles dans un
ordre tiré avant transmission. Les analyses à 3 et 5 familles prennent les
préfixes emboîtés de cet ordre ; les comparaisons de taille familiale sont donc
appariées et leur unité reste l'événement. Elles ne produisent pas
`500 × 3` événements indépendants. Les copies `H1` et `H2` sont équilibrées à
une unité près ; `BOTH` constitue 10 % des événements et est stratifié dans les
résumés plutôt que compté comme une troisième répétition.

Par partition, le fond négatif fondateur comprend 1 000 événements
`FND-NONE`, 1 000 `FND-RECURRENT` et 1 000 `FND-BACKGROUND-IBS`. Les 500
`FND-SINGLE-FAMILY` sont une classe d'insuffisance distincte, non ajoutée au
fond nul primaire. Pour chaque âge et chacun des mécanismes `LEFT`, `RIGHT` et
`DOUBLE`, 100 événements de défi censuré sont générés ; ils servent aux règles
de censure et jamais aux métriques principales non censurées.

Les tirages de fond de `16B` sont `99` au micro-panel, `999` au pilote et
`9 999` à la qualification. L'estimation utilise la correction
`(b + 1) / (N_DRAW + 1)` ; le nombre de tirages est un niveau de résolution,
pas une unité biologique. Les fonds interne et externe sont exécutés
séparément et aucun tirage n'est partagé entre partitions.

Pour `16C`, les 500 paires positives de chacune des cinq longueurs sont
distinctes. Elles sont réparties en mondes d'au plus 100 paires afin de borner
les appels externes, sans partager de sujet entre mondes. Les 3 000 paires
négatives sont réparties en 30 mondes de 100 paires disjointes ; aucun sujet
n'apparaît dans deux paires négatives. Les scénarios `COMPOUND` comportent
100 paires par longueur et restent exclus des estimations primaires. Une paire
négative est confirmée sur la vérité dense avant dérivation, conformément à
R05.

La révision R07 v2 ferme aussi le nombre de mondes : 25 mondes positifs
primaires, 30 négatifs, 30 de localisation secondaire et 5 `COMPOUND`, soit
90 mondes indépendants. Les deux outils et trois densités produisent au plus
540 appels externes. Un monde peut être projeté vers plusieurs densités et
outils, mais son effectif biologique ne compte toujours qu'une fois.

### 7.7 Fenêtres génétiques et censure

Pour un âge vrai `g`, la distance `X` au premier événement sur un bras suit le
modèle `X ~ exponentielle(g / 100)`, avec `X` en cM. Pour une demi-fenêtre de
largeur `W`, `P(censure) = exp(-g × W / 100)`. La série principale fixe une
probabilité théorique maximale de 0,05 par bras, donc :

`W_min(g) = -100 × ln(0,05) / g = 299,573 / g cM`.

| Âge `g` | Demi-fenêtre minimale `W_min` |
|---:|---:|
| 5 | 59,915 cM |
| 10 | 29,957 cM |
| 15 | 19,972 cM |
| 20 | 14,979 cM |
| 30 | 9,986 cM |
| 40 | 7,489 cM |
| 50 | 5,991 cM |
| 75 | 3,994 cM |
| 100 | 2,996 cM |
| 150 | 1,997 cM |

Les limites bp sont obtenues par interpolation de la carte. L'extrapolation,
`1 cM = 1 Mb` et le rabotage aux bornes sont interdits. Un couple cible–âge
sans les deux marges est `INFEASIBLE_WINDOW` avant tirage.

La série naturelle utilise `W_min` sans conditionner les événements observés.
Les défis gauche et droit réduisent uniquement le bras concerné à
`0,5 × W_min`; le défi double réduit les deux bras à cette valeur. Les
sentinelles `LEFT-EDGE` et `RIGHT-EDGE` garantissent en plus au moins un cas
d'impossibilité de fenêtre, mais ne contribuent pas au biais, à la RMSE ou à
une courbe de correction. Aucun événement n'est conservé ou rejeté selon la
longueur estimée par le logiciel.

### 7.8 Budget de ressources et modèle d'extrapolation

Les plafonds sont des portes, pas des objectifs à consommer. Ils incluent les
sorties, temporaires et caches propres à la partition.

| Échelle | Parallélisme | RSS / Java heap | Stockage alloué | Durée murale / appels | Budget agrégé |
|---|---|---|---|---|---|
| micro | 8 threads maximum | 16 Gio RSS, heap ≤ 12 Gio | 100 Gio | 12 h pour la chaîne | 100 heures-cœur |
| pilote calibration | 16 threads maximum | 64 Gio RSS, heap ≤ 48 Gio | 1 Tio | 7 jours pour la campagne pilote | 12 000 heures-cœur |
| qualification complète R07 v2 | profil SHAPEIT5 mesuré à 1 thread, au plus 24 workers R2 simultanés | 192 Gio RSS cumulée, heap d'un worker R2 ≤ 144 Gio | 8 Tio pour la release complète | 120 jours calendaires | 120 000 heures-cœur |

Un appel `R0` ou `R3` est limité à 2 h, un appel `R1` à 12 h, un appel `R2`
régional à 24 h et un appel `R2` autosomal/IBD à 168 h en qualification R07
v2. Cette dernière limite dépasse la projection sûre SHAPEIT5 de 126,70 h sans
créditer un gain de parallélisation non mesuré. Un seul nouvel essai est
autorisé pour une panne d'infrastructure documentée ; il reçoit un nouvel
identifiant de run et ne remplace pas la répétition rejetée. Aucun cache
scientifique ne traverse une partition.

Le budget de qualification reste 120 000 heures-cœur, allouées avant nouveau
résultat : 75 000 au phasage, 30 000 à l'IBD explicite et 15 000 aux autres
étapes. Le calendrier SHAPEIT5 exact comporte 640 appels. La projection R11,
pondérée par les densités, vaut 53 399,59 heures-cœur ; avec une réserve
supplémentaire de 25 %, elle vaut 66 749,48 heures-cœur et tient dans
l'allocation de phasage. Les 24 workers consommeraient au plus 160,06 Gio selon
la borne sûre de 6,669 Gio par worker, sous le plafond cumulé de 192 Gio.

Le pilote mesure séparément chaque étape et chromosome. Pour `R0/R3`, le coût
est projeté par nombre de lignes ; pour PLINK/bcftools par `N_S × N_V`; pour
KING et les tables de paires par le nombre de paires effectivement évaluées ;
pour le phasage et l'IBD par chromosome, outil, densité et taille de référence.
La borne supérieure à 95 % de chaque projection est multipliée par un facteur
de sécurité `2`. Une hypothèse linéaire non vérifiée ne peut pas ouvrir la
porte de qualification.

Si le besoin prévu excède un plafond, l'ordre obligatoire est : réduire le
parallélisme sans changer les données ; supprimer un cache reconstructible ;
planifier des lots de mondes indépendants ; sinon déclarer la cellule
`RESOURCE_BLOCKED`. Réduire `N_UNIT`, une densité ou une longueur pour faire
tenir le calcul est interdit sans nouvelle version de R07.

### 7.9 Adaptation prospective et porte R07

Le pilote utilise les 500 premières identités de calibration, les 100 premières
répétitions de chaque cellule légère et les dix premières graines des outils.
Après mesure, `N_FINAL` peut augmenter jusqu'aux réserves suivantes : 8 000
identités actives, 1 000 répétitions par cellule primaire, 6 000 unités nulles
et 100 graines d'outil. Il ne peut jamais descendre sous les minima nominaux.
Le même préfixe de taille finale est alors activé dans les deux partitions ;
la validation reste scellée jusqu'à R13.

Si la formule exige davantage que la réserve ou si la projection avec facteur
de sécurité dépasse les ressources, le domaine est bloqué avant validation.
Il faut une nouvelle version du dimensionnement et, si les identités ou
mécanismes changent, une nouvelle release ; aucun seuil R08 ne peut masquer
l'insuffisance.

R07 v2 est satisfaite au niveau du dimensionnement prospectif : unités
indépendantes, tailles nominales et réserves, niveaux de facteurs, plans
complets ou fractionnaires, répétitions, fenêtres, calendrier de 640 appels de
phasage, 90 mondes IBD et plafonds de ressources sont fixés sans réduire les
effectifs scientifiques. La configuration normative est
`calibration/configurations/scientific_dimensioning_v2.yaml` et le reçu
`calibration/readiness/R07_SCIENTIFIC_DIMENSIONING_V2.json`.

Cette clôture de R07 ne vaut pas qualification. R11b a depuis mesuré les deux
ponts d'échelle SHAPEIT5 dont la cible 4 000 + 3 000, la densité `D_LOW`, les
lots IBD de 100 paires pour les deux outils et trois densités, le stockage
complet et la faisabilité du calendrier à 24 workers. Les cinq contrôles de la
décision R11b passent, y compris la projection des 90 mondes avec borne à 95 %
et facteur de sécurité 2. Le reçu suivi
`calibration/readiness/R11B_RESOURCE_CONFIRMATION.json` autorise uniquement la
construction de R12c ; `scientific_g5_passed=false` jusqu'au gel de la nouvelle
release, de sa partition aveugle et de ses composants R13 scientifiques.

## 8. Préspecification des critères R08

### 8.1 Statuts et ordre de décision

R08 fixe des critères de qualification du logiciel sur données synthétiques ;
ils ne deviennent pas automatiquement des seuils de production. Chaque étape
reçoit un statut unique selon l'ordre suivant :

1. `FAIL` si une fuite de vérité, une violation d'intégrité, une sortie requise
   absente, une erreur silencieuse ou une valeur non finie non autorisée est
   observée ;
2. `NOT_EVALUABLE` si une impossibilité structurelle préspecifiée empêche de
   calculer la métrique malgré une exécution conforme ;
3. `FAIL` si l'estimation ponctuelle est du mauvais côté d'un seuil primaire ;
4. `NOT_EVALUABLE` si l'estimation ponctuelle est acceptable mais que sa borne
   de confiance ne franchit pas le seuil après activation de toute la réserve
   R07 ;
5. `REVIEW_REQUIRED` si toutes les métriques calculables passent mais qu'une
   adjudication humaine explicitement prévue reste ouverte ;
6. `PASS` seulement si toutes les métriques primaires obligatoires passent et
   qu'aucune condition précédente n'est vraie.

Une métrique secondaire ou exploratoire ne peut ni sauver une métrique
primaire en échec, ni transformer `NOT_EVALUABLE` en `PASS`. Une étape amont
requise en `FAIL` interdit le `PASS` aval. Une étape aval peut néanmoins être
évaluée sur une branche indépendante si son manifest démontre qu'elle ne
consomme pas l'artefact défaillant.

Le micro-panel ne produit que des verdicts d'intégrité technique. Les verdicts
scientifiques de calibration sont provisoires jusqu'au gel R12 ; ceux de
validation sont définitifs pour la release évaluée.

### 8.2 Intervalles, unités et multiplicité

Les notations sont `LCB95`, borne inférieure unilatérale à 95 %, `UCB95`, borne
supérieure unilatérale à 95 %, et `CI95`, intervalle bilatéral à 95 %. Toutes
les bornes utilisent l'unité indépendante R07, jamais le nombre brut de lignes.

- une proportion utilise Wilson ; une proportion nulle rare utilise la borne
  exacte binomiale unilatérale ;
- une moyenne, médiane, erreur absolue, RMSE, corrélation ou différence
  appariée utilise 10 000 rééchantillonnages par grappes, déterministes, avec
  `isolation_unit_id` ou `scenario_world_id` comme grappe ;
- une cellule dégénérée déterministe publie sa valeur exacte et n'invente pas
  un intervalle ;
- les intervalles Gamma produits par `14` sont évalués comme sorties de la
  méthode ; ils ne remplacent pas l'intervalle de calibration de leur taux de
  couverture ;
- un tirage nul utilise `p_MC = (b + 1) / (N_DRAW + 1)` et enregistre aussi
  l'incertitude binomiale de Monte-Carlo ;
- les graines de bootstrap et d'évaluation appartiennent au domaine privé
  `EVALUATION`, dérivé comme en 5.13, et ne sont jamais choisies après lecture
  du résultat.

Les métriques co-primaires d'une étape sont toutes exigées : leur décision est
un test d'intersection–union et n'est pas assouplie par une moyenne entre
critères. Lorsqu'une conclusion peut être déclenchée par l'une de plusieurs
comparaisons, les valeurs `p` sont corrigées par Holm à `alpha = 0,05` dans la
famille de l'étape. Les 25 étapes correspondent à des qualifications séparées ;
aucune conclusion globale ne fusionne leurs valeurs `p`. Les analyses de
strates secondaires publient leurs intervalles sans créer de nouvelle
qualification.

Pour une exigence minimale `T`, `PASS` requiert `LCB95 >= T`. Pour une erreur
maximale, il requiert `UCB95 <= T`. Pour une équivalence de marge `M`, tout le
`CI95` doit être inclus dans `[-M, +M]`. Les égalités aux seuils passent. Les
valeurs affichées sont arrondies seulement après décision sur leur précision
complète.

### 8.3 Valeurs manquantes et non-évaluables

Une absence est encodée `NA` avec un code, jamais `0`, chaîne vide, infini ou
exclusion silencieuse. Chaque métrique porte `N_EXPECTED`, `N_ELIGIBLE`,
`N_EVALUABLE`, `N_EXCLUDED` et les motifs. Sauf règle plus stricte ci-dessous,
`LCB95(N_EVALUABLE / N_ELIGIBLE) >= 0,95` est co-primaire.

Les seuls motifs structurels autorisant `NOT_EVALUABLE` sont :

- `INFEASIBLE_TARGET_CATALOG`, `INFEASIBLE_WINDOW` ou `INFEASIBLE_ROH`
  déterminé avant tirage ;
- effectif indépendant, densité, nombre de copies mineures ou transitions
  informatives sous le minimum R07 malgré activation de la réserve ;
- cellule optionnelle sans positif ou négatif par conception signée ;
- métrique de proportion d'ascendance en `16A`, car la méthode actuelle ne
  contient aucun estimateur de proportion ;
- scénario censuré de `14` si la méthode identifiée ne possède pas de modèle
  de censure, à condition que la censure soit correctement signalée ;
- domaine bloqué avant exécution par la porte de ressources R07.

Une dépendance absente, un outil manquant après `READY_TO_RUN=true`, un timeout,
un parseur défaillant, une sortie vide, une exclusion excessive produite par
la méthode ou un dénominateur perdu est un `FAIL`, pas un non-évaluable. Un
motif inconnu donne `REVIEW_REQUIRED` et bloque la qualification.

### 8.4 Critères des étapes 00 à 10

| Étape | Métriques primaires et unité | Critères `PASS` |
|---|---|---|
| `00` | complétude des huit racines, empreintes, isolation, reprise ; unité = artefact/run | 100 % des artefacts attendus présents et valides ; 100 % des altérations et chemins privés refusés ; reprise donnant les mêmes SHA-256 scientifiques ; aucune écriture dans un run clos |
| `01` | classification des sources et couverture ; unité = source/sentinelle | sensibilité et spécificité exactes de 100 % pour les cas déterministes ; 22 autosomes et cible couverts ; assemblage, licence, version et empreinte tous renseignés |
| `02` | bijection source–échantillon–identité, pedigree déclaré et défauts d'identité ; unité = identité/composante | 100 % des lignes valides conservées ; 100 % des collisions, échanges, sources partagées et filiations impossibles détectés ; zéro génotype inféré d'un groupe ou phénotype |
| `03` | concordance sujet/variant/génotype observé et manque ; unité = sujet et variant, grappée par monde | concordance exacte de 100 % avec l'observation synthétique, compte de manque exact, aucune perte non déclarée, aucun allèle réorienté silencieusement ; les sentinelles incompatibles sont toutes rejetées |
| `04` | matrice vérité–appel cible, provenance et Mendel ; unité = identité/trio | concordance de 100 % sur les appels non manquants ; manque préservé à 100 % ; provenance moléculaire présente pour chaque appel ; toutes les erreurs mendéliennes injectées signalées et aucun génotype déduit du statut |
| `05` | détection des défauts QC et fidélité des métriques ; unité = sujet, variant ou lot selon le défaut | niveau haut : sensibilité `LCB95 >= 0,90` ; niveau bas : `LCB95 >= 0,80` ; spécificité `LCB95 >= 0,99` ; erreur absolue des taux de manque/appel `UCB95 <= 0,005` ; évaluabilité `LCB95 >= 0,95` |
| `06` | fidélité de l'ensemble de parenté, couverture et LD résiduel ; unité = variant/bloc | rappel et précision de l'éligibilité `LCB95 >= 0,99` ; 22 chromosomes et minima 7.4 respectés ; fraction de paires dépassant le seuil de pruning configuré `UCB95 <= 0,01` ; aucun variant interdit ou cible protégée utilisé contrairement au contrat |
| `07` | classification duplicat/degrés 1–4/non-apparenté et erreur de coefficient ; unité = paire désignée, grappée par composante | par classe liée, sensibilité et VPP `LCB95 >= 0,90` ; spécificité non-apparentée `LCB95 >= 0,99` ; MAE du coefficient de parenté `UCB95 <= 0,03` ; détection cryptique `LCB95 >= 0,90` |
| `08` | stabilité PCA, séparation des états purs et détection d'atypiques ; unité = représentant indépendant | zéro apparenté dans l'ajustement ; exactitude du centroïde le plus proche pour états purs `LCB95 >= 0,90` ; sensibilité atypique `LCB95 >= 0,90` et spécificité `LCB95 >= 0,95` ; déplacement Procrustes médian `UCB95 <= 0,10` écart-type et percentile 95 `<= 0,25` |
| `09` | concordance des cohortes et indépendance familiale ; unité = identité/composante | 100 % des inclusions/exclusions conformes à la politique signée ; zéro composante divisée entre cohortes incompatibles ; cohortes obligatoires non vides ; empreintes de revue valides. Une revue prévue mais non signée donne `REVIEW_REQUIRED` |
| `10` | décisions QC finales par cohorte et préservation de la cible ; unité = sujet/variant/cohorte | défauts hauts : sensibilité `LCB95 >= 0,90` ; défauts bas : `LCB95 >= 0,80` ; spécificité `LCB95 >= 0,99` ; 100 % des exceptions cible conformes ; aucune modification de cohorte non auditée |

Les étapes `00`–`04` et `09` sont des contrats déterministes : une seule erreur
est un échec, sans intervalle destiné à diluer le défaut. Pour `05` et `10`,
les seuils évaluent la politique configurée sur les défauts injectés ; ils ne
créent pas une nouvelle valeur de filtre de production.

### 8.5 Critères des étapes 10A à 16

| Étape | Métriques primaires et unité | Critères `PASS` |
|---|---|---|
| `10A` | complétude autosomale et ordre ; unité = chromosome/sujet/variant | exactement 22 autosomes ; 100 % des sujets dans le même ordre ; zéro génotype structurellement absent, variant dupliqué ou chromosome vide ; tous les retraits justifiés |
| `11` | exactitude cible, carte, fenêtre et faisabilité ; unité = cible–âge | cible/allèles exacts ; erreur cM `<= 10^-6 cM` et erreur de borne `<= 1 bp` par rapport à l'interpolation attendue ; aucune extrapolation ; densités 7.4 respectées ; classification des cas infaisables exacte à 100 % |
| `12` | SER, erreur ponctuelle de phase, copie porteuse et concordance GT ; unité = sujet/famille, transitions en sous-unités | concordance GT = 100 % ; à `D_MAX/D_MID` avec pedigree complet : SER commun `UCB95 <= 0,01`, SER MAF 0,01–0,05 `<= 0,03`, erreur ponctuelle `<= 0,02` ; à `D_LOW` ou sans pedigree : limites `0,03`, `0,05`, `0,05` ; copie H1/H2/BOTH `LCB95 >= 0,95` au témoin et `>= 0,90` en stress ; évaluabilité `LCB95 >= 0,95` |
| `12A` | phase par chromosome, harmonisation et complétude ; unité = sujet/chromosome | 22/22 chromosomes et 100 % des sujets produits ; zéro discordance allèlique acceptée ; SER commun `UCB95 <= 0,02` à `D_MAX/D_MID` et `<= 0,05` à `D_LOW` ; erreur ponctuelle `<= 0,03` et `<= 0,07` respectivement |
| `13` | détection du segment fondateur, faux positifs, recouvrement et bornes ; unité = événement fondateur | pour âge `<= 50` et `N_F >= 5`, sensibilité `LCB95 >= 0,90` et Jaccard médian `LCB95 >= 0,80` ; pour âge `>= 75`, seuils `0,80` et `0,65` ; avec `N_F = 3`, retrancher `0,10` avec planchers `0,70/0,55` ; pour chaque classe nulle, FPR `UCB95 <= 0,01` ; erreur normalisée médiane de borne `UCB95 <= 1`, avec chaque erreur absolue divisée par `max(0,5 cM, 0,25 × longueur vraie)` ; cible contenue dans 100 % des segments positifs appelés |
| `14` | biais multiplicatif, MdAPE, couverture Gamma et stabilité familiale ; unité = événement | `N_F=10` : biais relatif médian `<= 0,10`, MdAPE `<= 0,25`, ratio médian borne haute/basse `<= 4` ; `N_F=5` : `0,15/0,35/6` ; `N_F=3` : `0,25/0,50/10` ; couverture nominale 95 % avec `LCB95 >= 0,90` ; variation leave-one-family-out médiane `UCB95 <= 0,20` ; censure reconnue `LCB95 >= 0,90` |
| `15` | erreur de `r²`, `D′`, biais et stabilité par densité ; unité = bloc haplotypique | pour MAF `>= 0,05`, MAE `r² <= 0,02`, MAE `D′ <= 0,05`, biais absolus `<= 0,01` et `<= 0,02` ; delta apparié médian entre densités `<= 0,03` ; évaluabilité `LCB95 >= 0,95`. Une strate avec moins de 20 copies mineures est `NOT_EVALUABLE` |
| `16` | détection ROH, FPR, bornes, longueur, fardeau et cible ; unité = sujet-monde | sensibilité `LCB95 >= 0,80` pour `SHORT` et `>= 0,95` pour `MEDIUM/LONG/MULTIPLE_BURDEN` ; FPR `UCB95 <= 0,01` ; MAE de borne `<= 0,5 cM` ; erreur relative médiane de longueur `<= 0,20` courte et `<= 0,10` sinon ; exactitude cible-dans-ROH `LCB95 >= 0,95` ; erreur relative du fardeau `<= 0,10` ; fusion/fragmentation correctement classée avec `LCB95 >= 0,90` |

Pour `12` et `12A`, une cellule exige au moins 20 copies mineures et 100
transitions informatives indépendantes ; sinon elle suit 8.3. Une région
marquée non fiable par la méthode reste dans le dénominateur d'évaluabilité.
Pour `14`, le biais relatif médian vaut
`abs(exp(median(log(g_estime / g_vrai))) - 1)` ; une estimation non positive
est un `FAIL`. Les cellules censurées sans modèle compatible ne contribuent pas
aux métriques non censurées.

### 8.6 Critères des étapes 16A à 19

| Étape | Métriques primaires et unité | Critères `PASS` |
|---|---|---|
| `16A` | indépendance du modèle, positionnement global/local et stabilité ; unité = représentant indépendant | ajustement exclusivement sur `REF-ANCESTRY-SYN` ; exactitude du centroïde pour états purs `LCB95 >= 0,95`, exactitude locale `>= 0,90`, corrélation de rang des mélanges `LCB95 >= 0,90` ; déplacement Procrustes médian `<= 0,10` écart-type et percentile 95 `<= 0,25` ; toute proportion d'ascendance reste `NOT_EVALUABLE` |
| `16B` | FPR nul, puissance, évaluabilité et concordance des fonds ; unité = événement/famille indépendante | `N_DRAW = 9 999`, formule plus-un exacte et décision positive si `p_MC <= 0,05` ; FPR `UCB95 <= 0,05` séparément pour `NONE`, `RECURRENT` et `BACKGROUND-IBS` ; puissance `LCB95 >= 0,80` pour âge `<= 50`, `N_F >= 5`, et `>= 0,70` pour autres cellules faisables ; concordance de décision fonds interne/externe `LCB95 >= 0,90` ; évaluabilité `LCB95 >= 0,95` |
| `16C` | sensibilité par longueur, FPR cible, bornes et concordance outils ; unité = paire indépendante | sensibilité `LCB95 >= 0,70/0,80/0,90/0,95/0,95` pour `1/1,5/2/3/5 cM` ; FPR sur paires négatives `UCB95 <= 0,001` ; erreur médiane de chaque borne `<= 0,50 cM` pour 1–1,5 cM et `<= 0,25 cM` sinon ; Jaccard Hap-IBD/Refined IBD `LCB95 >= 0,60` pour 1–1,5 cM et `>= 0,80` sinon ; alignement de phase cible exact à 100 % ; contrôle `16A` disponible |
| `17` | validité des paires, delta signé, équivalence et convergence ; unité = paire de runs/unité indépendante | les seize domaines présents, manifests clos et diff à facteur unique exact ; `DEGRADATION_DETECTABLE` : delta standardisé signé `LCB95 >= 0,20` et concordance de signe `LCB95 >= 0,80` ; `MONOTONIC` : `LCB95 >= -0,20` ; `EQUIVALENT` : `CI95` inclus dans `[-0,20,+0,20]` ; graines d'outil : écart 50–100 graines `<= 0,10` écart-type |
| `18` | fidélité des SVG, unités, statuts et provenance ; unité = figure/valeur autorisée | 100 % des valeurs, unités, effectifs, intervalles, légendes et statuts identiques aux tables closes ; provenance et SHA-256 complets ; cas `PASS/FAIL/NOT_EVALUABLE` rendus ; zéro recalcul, donnée individuelle, HTML ou PDF |
| `19` | couverture, traçabilité et fidélité des faits ; unité = fait/Markdown/figure | 24/24 Markdown amont présents ; 100 % des faits reliés à une source et une figure ou justification ; verdicts et limites reproduits à 100 % ; zéro fait orphelin, donnée sensible, conclusion d'effet fondateur ou promotion automatique ; revue humaine signée |

Dans `16A`, les centroïdes et corrélations sont des métriques d'évaluation
calculées après clôture ; ils n'ajoutent pas un estimateur de proportion au
pipeline. Dans `16C`, les outils sont qualifiés séparément en plus de leur
concordance : un accord entre deux appels faux ne permet pas un `PASS`.

### 8.7 Directions des sensibilités

Chaque paire `17` reçoit avant génération l'une des classes privées suivantes :

| Facteur | Direction préspecifiée | Classe |
|---|---|---|
| densité `D_MAX -> D_LOW` | une densité moindre ne doit pas améliorer une erreur ou une détection | `MONOTONIC` |
| manque, erreur génotypique ou lot `NONE -> HIGH` | les métriques de fidélité doivent se dégrader | `DEGRADATION_DETECTABLE` |
| pedigree `COMPLETE -> NO_PEDIGREE` | la phase et la résolution familiale ne doivent pas s'améliorer | `MONOTONIC` |
| référence nominale -> réduite | la phase et le positionnement ne doivent pas s'améliorer | `MONOTONIC` |
| erreur de phase `NONE -> HIGH` | les métriques dépendant de la phase doivent se dégrader | `DEGRADATION_DETECTABLE` |
| nombre de familles `3 -> 10` | détection, précision et stabilité ne doivent pas se dégrader | `MONOTONIC` |
| longueur IBD ou ROH croissante | sensibilité et erreur relative ne doivent pas se dégrader | `MONOTONIC` |
| graine d'outil ou ordre de lots sans changement scientifique | aucun changement matériel | `EQUIVALENT` |

Une métrique d'erreur est orientée avant calcul pour que le delta positif soit
toujours favorable. Une direction ne peut pas être inversée après ouverture de
`T-SCENARIO`. Les contrastes supplémentaires d'une même famille utilisent
Holm ; les huit contrastes principaux ci-dessus sont tous obligatoires dans
leurs domaines applicables.

### 8.8 Gel, amendements et porte R08

Les critères de cette section sont la version initiale préspecifiée. Après le
pilote de calibration, R12 peut proposer un amendement seulement si : la
validation demeure inaccessible ; la justification scientifique est écrite ;
une nouvelle version de critères est signée ; tous les résultats de calibration
affectés sont invalidés puis recalculés ; et le changement n'est pas une simple
réécriture destinée à faire passer une sortie observée. Aucun amendement n'est
possible après ouverture de la validation.

R08 est satisfaite au niveau de conception : les 25 identifiants possèdent une
unité, des métriques primaires, des seuils, une méthode d'intervalle, une règle
de multiplicité, un traitement des absences et des domaines non évaluables.
R09 devra les encoder sans les modifier dans les contrats d'évaluation ; toute
divergence donnera `READY_TO_RUN=false`.

## 9. Contrats à créer après la conception

Les schémas ne seront recréés qu'après validation des besoins. Le socle minimal
comprendra :

1. `panel_release_manifest` : identité, version, sources, graines et couches ;
2. `subject_registry` : rôles, partitions et composantes familiales ;
3. `master_truth_manifest` : inventaire typé de toutes les vérités ;
4. `panel_view_manifest` : filtre, transformation et héritage de chaque vue ;
5. `blinded_input_inventory` : seules entrées accessibles au pipeline ;
6. `pipeline_run_manifest` : commandes, versions, sorties et clôture ;
7. `stage_evaluation_manifest` : vérité autorisée et métriques par étape ;
8. `readiness_receipt` : preuve que tous les prérequis sont satisfaits ;
9. `calibration_decision` : qualification, restrictions et compatibilité.

La porte d'exécution sera :

`READY_TO_RUN = contrat ∧ données ∧ effectifs ∧ indépendance ∧ outils ∧ aveuglement ∧ ressources ∧ résultats_attendus`.

Chaque terme devra être démontré par un artefact empreinté. Une valeur
implicite, un chemin absent ou un effectif insuffisant donnera `false`.

### 9.1 Réalisation R09

Le socle est matérialisé dans `calibration/schemas/` par exactement les neuf
schémas ci-dessus. Les configurations humaines versionnées sont :

- `contract_bundle_v1.yaml`, fermeture des 25 étapes, huit racines, quatorze
  vérités, onze vues, sorties dynamiques et onze écarts R04/R06 ;
- `stage_criteria_v1.yaml`, transcription des critères R08 liée à l'empreinte
  exacte de la section 8 ;
- `micro_panel_v1.yaml`, dimensionnement `RESOURCE_ONLY` R07 sans donnée
  générée ;
- `isolation_policy_v1.yaml`, séparation physique, pseudonymisation par
  partition, accès minimal à la vérité et registre à facteur unique ;
- `tool_harness_v1.yaml`, contrats simulés des outils et échecs attendus, sans
  autorisation d'exécution externe à R09.

Le noyau autorisé sous `src/effet_fondateur/calibration/` se compose de
`contracts.py`, `release_guard.py`, `truth_router.py`, `evaluation.py` et
`readiness.py`. Il refuse les chemins réels ou privés, les liens sortants, les
recouvrements de partition, l'accès à une vérité avant clôture et les demandes
de vérité hors périmètre. Les sorties dynamiques doivent être nommées et
empreintées ; `16C` exige le contrôle `16A`, `17` les seize domaines clos,
`18` reste limité aux SVG et `19` aux faits/Markdown avant R14.

Les anciens bancs exploratoires restent présents pour audit, mais sont exclus
du paquet exporté et de la liste des producteurs autorisés. Leurs 62 modules de
test, leur fixture et leur README ont été remplacés par la suite R09 ; aucun
ancien schéma, configuration ou protocole D n'a été restauré.

Le reçu `readiness/R09_READY_TO_RUN.json` est positif uniquement pour le
périmètre `R10_MICRO_PANEL_CONSTRUCTION`. Il prouve la complétude statique des
contrats, spécifications, effectifs, gardes, ressources et sorties attendues ;
il ne prétend pas qu'un panel existe et n'autorise ni verdict scientifique ni
outil externe. La vérification physique des données et empreintes devient
obligatoire à R10 avant toute exécution aveugle du pipeline.

### 9.2 Réalisation R10

La release `panel-master-generic-v1-r10-resource` matérialise le micro-panel
transversal entièrement synthétique sous
`calibration_workspace/releases/micro_panel_generic_v1_r10/`. Sa cible
`STANDARD` et quatre cibles de bord ou de recombinaison sont résolues depuis
`micro_panel_v1.yaml` ; aucune cible n'est codée en dur dans le constructeur.

La release comprend 60 sujets d'étude répartis en quatre trios, quatre
fratries de quatre, quatre composantes multigénérationnelles de cinq et douze
non-apparentés. Six sous-références disjointes totalisent 240 identités. Les
169 variants couvrent les 22 autosomes et une région cible dense. Quatorze
familles de vérité privée alimentent onze vues opaques, sans recouvrement de
référence ni exposition des identités maîtres.

Les 25 étapes contractuelles ont traversé quatre runs successivement clos :
`00–16C`, puis `17` après registre signé des seize domaines, `18` après les
évaluations closes, et `19` après les 24 Markdown amont. Chaque étape possède
son enveloppe, son inventaire d'entrée fermé, ses sorties et ses SHA-256 ; les
25 évaluations et Markdown concluent seulement sur
`TECHNICAL_INTEGRITY_ONLY`. L'accès direct du travailleur à la vérité privée
est refusé par l'isolation de processus et l'évaluateur ne l'ouvre qu'après
clôture selon la liste blanche de l'étape.

La validation indépendante compte 484 fichiers, 14 vérités, 11 vues, 25
étapes, 25 évaluations, 25 Markdown et aucun HTML/PDF. Un rejeu intégral dans
un répertoire temporaire a produit exactement les mêmes chemins et empreintes,
avec l'agrégat SHA-256
`5ff90ef966cbc7bdd29bc46414b6375a95c946b1667ee7e7cb53fae7fdecd01f`.
La release occupe 2 740 224 octets alloués après clôture, très sous la porte
R07 de 100 Gio. Aucun outil scientifique externe ni aucune donnée réelle n'a
été lu ou exécuté : R10 démontre l'intégrité technique et la faisabilité du
micro-panel, pas une performance scientifique.

## 10. Résultat obligatoire de chaque étape

Chaque étape exécutée produira :

`calibration/resultats/ETAPE_<identifiant>.md`

Ce fichier contiendra :

1. identifiant, date et version de méthode ;
2. question et hypothèse testées ;
3. identifiant et empreinte de la release ;
4. vues, partitions et effectifs réellement consommés ;
5. commandes, paramètres, graines et versions ;
6. vérité ouverte par l'évaluateur ;
7. métriques, unités et intervalles d'incertitude ;
8. comparaison aux critères préspecifiés ;
9. décision `PASS`, `FAIL`, `NOT_EVALUABLE` ou `REVIEW_REQUIRED` ;
10. exclusions, anomalies et écarts au protocole ;
11. limites d'interprétation ;
12. modifications du logiciel à envisager, sans mise en œuvre implicite ;
13. chemins et SHA-256 des tables nécessaires aux courbes finales.

Le Markdown ne contiendra aucune donnée génétique individuelle ni vérité
sensible. Les tables volumineuses, données graphiques et sorties resteront
dans le workspace de la release.

Aucun rapport HTML/PDF intermédiaire ne sera créé.

## 11. Déroulement définitif de la nouvelle calibration

### R00 — Fixer le périmètre

**État :** réalisé.

Inventorier les anciens espaces, distinguer calibration et logiciel, et
valider la remise à zéro sans archive.

### R01 — Nettoyer l'ancienne campagne

**État :** réalisé.

Supprimer définitivement tout le contenu de `calibration/`, sauf le présent
plan, et tout le contenu de `calibration_workspace/`. Ne rien supprimer sous
`src/`, `tests/`, `docs/`, `config/`, `data/` ou
`simulation_genotype_famille/`.

**Porte :** seul ce plan subsiste sous `calibration/`, le workspace est vide
et aucun fichier logiciel n'a été touché.

### R02 — Vérifier l'intégrité du logiciel

**État :** réalisé.

Exécuter les imports, les tests ciblés puis la suite de référence adaptée. Les
échecs préexistants seront distingués des régressions.

**Porte :** aucune régression imputable au nettoyage.

### R03 — Compléter la matrice des besoins

**État :** réalisé.

Pour chaque étape de la section 6, préciser formats, unités, vérités, tailles,
outils, ressources, erreurs attendues et données graphiques.

**Porte :** chaque entrée possède un producteur et chaque sortie un
consommateur ou une justification.

### R04 — Établir le graphe des dépendances

**État :** réalisé.

Décrire le chemin producteur–artefact–consommateur de la génération jusqu'au
résultat final.

**Porte :** aucun artefact implicite, aucune vérité orpheline et aucune
dépendance circulaire.

### R05 — Concevoir le panel maître

**État :** réalisé.

Figer sujets, références, haplotypes, pedigrees, événements fondateurs, IBD,
ROH, structure, erreurs et graines.

**Porte :** tous les besoins R03–R04 sont couverts avant génération.

### R06 — Définir les partitions et les vues

**État :** réalisé.

Figer les partitions indépendantes et toutes les transformations dérivées du
panel maître.

**Porte :** disjonctions démontrées, vues appariées et vérité inaccessible au
pipeline.

### R07 — Dimensionner les scénarios

**État :** révisé en v2 après R14 ; dimensionnement scientifique fermé et
confirmation de ressources R11b terminée. Une nouvelle release indépendante
reste requise avant G5 scientifique.

Justifier les effectifs, répétitions, précisions cibles, ressources et plans
factoriels.

**Porte :** dimensionnement signé avant observation d'une performance.

### R08 — Préspecifier les critères par étape

**État :** réalisé.

Fixer métriques, seuils, intervalles, règles de valeurs manquantes et domaines
non évaluables.

**Porte :** aucune décision ne dépend d'un résultat encore inconnu.

### R09 — Recréer les contrats et configurations

**État :** réalisé.

Créer seulement les schémas, configurations et gardes nécessaires au panel
unifié. Réutiliser le code existant uniquement après vérification de son
adéquation.

**Porte :** validation statique, tests ciblés et reçu `READY_TO_RUN` du
micro-panel.

### R10 — Construire le micro-panel transversal

**État :** réalisé.

Générer une petite release capable de traverser réellement les étapes 00–19,
y compris les 22 autosomes lorsque l'étape l'exige.

**Porte :** chaîne techniquement complète, vérités séparées, artefacts
traçables ; aucune conclusion de performance.

### R11 — Exécuter le pilote de calibration

**État :** rouvert le 13 novembre 2026 ; porte non franchie et ressources
`RESOURCE_BLOCKED` au profil mesuré.

**R11b :** confirmation post-R07 v2 terminée le 14 novembre 2026 dans la portée
`RESOURCE_CONFIRMATION_ONLY`. Les trois sondes SHAPEIT5 et les 18 sondes IBD
sont closes. La décision `RESOURCE_READY_FOR_NEW_SCIENTIFIC_RELEASE` autorise
la construction d'une nouvelle release, sans qualifier rétroactivement R12b,
R13 ou R14.

Exécuter les outils réels sur la partition de calibration et produire les
Markdown par étape.

**Porte :** résultats exploitables, ressources mesurées, causes des
non-évaluables comprises.

Le pilote initial a matérialisé 500 identités et les trois densités emboîtées. PLINK,
KING et bcftools ont réussi et leurs ressources sont mesurées. Après un premier
essai rejeté pour incompatibilité du wrapper `/usr/bin/time`, l'unique nouvel
essai autorisé a lancé dix graines SHAPEIT5 5.1.1. Elles ont toutes échoué de
façon déterministe avant phasage, car le VCF synthétique du harnais R11 ne
portait pas le champ INFO `AC` requis. Les outils IBD dépendants n'ont donc pas
été exécutés scientifiquement. Les 25 Markdown ont été clos avec 3 `PASS`
techniques pilotes, 20 `NOT_EVALUABLE` et 2 `FAIL` ; aucun de ces statuts ne
vaut qualification.

Après l'audit pré-R13, R11 a été rouvert pour mesurer les ressources manquantes
sur un chromosome synthétique médian de 26 000 variants communs. Cinq points
égaux `125+125`, `187+187`, `250+250`, `375+375` et `500+500`, un contrôle
d'allocation `286+214` et un contrôle apparié à demi-densité ont terminé avec
15 itérations chacun. Le modèle puissance log–log préspecifié projette
4 000 sujets et 3 000 références avec une borne prédictive unilatérale à 95 %,
la correction observée du ratio `4:3` et le facteur de sécurité 2. Il donne
126,70 h murales et 121,82 heures-cœur par appel médian, avec 6,67 Gio de RSS.
Le temps dépasse le plafond R2 de 72 h ; `22 autosomes × 50 graines` atteignent
134 000 heures-cœur pour le seul phasage commun, au-delà du budget global de
120 000 heures-cœur avant les autres outils. R11 est donc
`RESOURCE_BLOCKED`. Aucun crédit de parallélisation non mesuré n'est appliqué,
aucun effectif ou niveau de densité R07 n'est réduit. La qualification
scientifique reste fermée ; une validation technique indépendante est traitée
séparément par R12b/R13.

### R12 — Corriger puis figer

**État :** candidat initial rejeté après audit pré-R13 ; remplacé pour la seule
validation technique par R12b.

Les défauts démontrés peuvent conduire à des propositions de correction du
logiciel. Chaque correction autorisée sera testée et obligera à reprendre les
étapes affectées. Ensuite, méthodes et seuils seront empreintés et scellés.

R12 a ajouté `AC/AN` aux deux VCF régionaux synthétiques par
`bcftools +fill-tags`, sans modifier R11, les effectifs, les densités ni les
critères R08. Les dix graines SHAPEIT5 5.1.1, Hap-IBD et Refined IBD ont
terminé sans échec. Les trois étapes reprises (12, 12A et 16C) restent
`NOT_EVALUABLE` parce que le pilote ne contient pas leurs grilles complètes ;
ce statut ne vaut pas qualification. Le candidat, ses méthodes, critères,
compatibilités, outils et modules sont scellés par SHA-256. La validation
aveugle n'a pas été ouverte et aucun rapport HTML/PDF n'a été produit.

L'audit préalable à R13 a ensuite invalidé le reçu G5 : ce candidat ne contient
ni partition `VALIDATION_BLIND`, ni vérité de validation scellée, ni exécuteur
et évaluateur R13 inclus dans les empreintes. De plus, R11 est resté
explicitement sous sa porte avec vingt étapes `NOT_EVALUABLE` et n'a pas produit
les projections R07 requises pour 4 000 sujets, 3 000 références de phase,
22 autosomes, 50 graines et les mondes IBD de qualification. Le plan interdit
une extrapolation linéaire non vérifiée. Le candidat demeure conservé pour
audit mais est rejeté pour R13 ; aucune donnée de validation n'a été ouverte.

**Porte :** code candidat, critères et compatibilités figés avant validation.

#### R12b — Gel indépendant pour la validation technique

**État :** réalisé le 13 novembre 2026 ; G5 franchie dans une portée technique
explicite.

R12b matérialise une nouvelle release non écrasable avec la partition
`VALIDATION_BLIND`, onze vues aveugles fermées et quatorze familles de vérité
privée scellées. Son espace de noms, ses identités synthétiques et ses
pseudonymes sont indépendants des entrées R10–R12. Le runner R13,
l'évaluateur post-clôture, le worker, la topologie et leurs configurations sont
figés par SHA-256 avant toute exécution.

Le reçu `calibration/readiness/R12B_G5_READY.json` est positif uniquement pour
`R13_TECHNICAL_BLIND_VALIDATION`. R12b ne contient encore aucun run ni aucune
évaluation, et la vérité n'a pas été ouverte. Les drapeaux de qualification
restent faux.

R07 est marqué `REVISION_REQUIRED` pour les effectifs de qualification,
l'enveloppe SHAPEIT5, le budget des 22 autosomes × graines et les mondes IBD.
Ces lacunes interdisent une qualification scientifique complète, mais ne
bloquent pas le contrôle d'intégrité technique R13 prescrit par R12b.

#### R12c — Nouvelle release de validation scientifique

**État :** release scientifique gelée et G5 obtenu le 14 novembre 2026 ; un
checkpoint R13 grandeur réelle a été terminé, puis la campagne exhaustive a
été remplacée par l'amendement d'intégration décrit en R13.

Le contrat `scientific_validation_v1.yaml` reprend sans réduction les minima
R07 v2, les 640 appels SHAPEIT5, les 90 mondes IBD et les 25 critères R08. Il
utilise un nouvel espace de noms `VALIDATION_BLIND_R12C_V1`, ne réutilise aucun
sujet ou artefact aveugle R10–R12b et conserve une cible autosomale résolue par
profil plutôt que codée en dur.

Le préflight `R12C_SCIENTIFIC_PREFLIGHT.json` lie par SHA-256 R07, R08 et le
reçu de ressources R11b. La projection est réordonnancée sans réduction sur
huit workers, soit 90,802 jours sous la limite de 120 jours.

Les contrats initiaux `r13_scientific_runner.py` et
`r13_scientific_evaluator.py` valident respectivement la fermeture du plan de
640 appels SHAPEIT5 et 540 appels IBD, puis la couverture de chaque métrique
primaire R08. Le runner développe désormais les 640 appels en
`340 D_MAX + 150 D_MID + 150 D_LOW` et les 90 mondes IBD en 540 appels
outil–densité ; toute dérive de cette liste est bloquante. L'évaluateur
applique aussi l'ordre `FAIL_INTEGRITY` → `NOT_EVALUABLE_STRUCTURAL` →
`FAIL_POINT_ESTIMATE` → `NOT_EVALUABLE_PRECISION` → `REVIEW_REQUIRED` →
`PASS` aux métriques post-clôture. Ils sont gelés en états
`FULL_EXECUTION_READY` et `FULL_EVALUATION_READY`.

Le fichier `scientific_metric_rules_v1.yaml` fige désormais 135 règles, dont
21 contextuelles, couvrant les 25 étapes R08 par deux fragments empreintés. Les
statistiques de décision suivent strictement 8.2 : `LCB95` pour un minimum,
`UCB95` pour un maximum, `CI95` pour une équivalence et valeur ponctuelle
uniquement pour un contrôle déterministe. Cette fermeture retire le composant
`scientific_metric_rules` de la liste des manques R12c, sans rendre
l'évaluateur productif et sans ouvrir G5.

Le socle producteur fixe aussi la méthode de chaque règle avant résultat : 65
contrôles déterministes sans intervalle, 39 proportions Wilson unilatérales,
4 FPR rares avec borne exacte de Clopper–Pearson et 27 statistiques continues
avec 10 000 bootstraps par `isolation_unit_id` ou `scenario_world_id`. Les
graines privées sont dérivées par HMAC-SHA-256 vers PCG64DXSM ; la correction
de Holm et `p_MC=(b+1)/(N_DRAW+1)` sont implémentées. Cet audit vaut seulement
`METHODS_AND_STAGE_BINDINGS_FROZEN`.

Le raccordement spécifique couvre les 25 étapes `00` à `19`. Une preuve
manquante bloque la production et une discordance devient une observation en
échec. L'état est `FULL_STAGE_PRODUCERS_FROZEN` (`25/25`). La release contient
4 000 identités pseudonymisées et 14 familles de vérité scellées. Le reçu
`R12C_G5_READY.json` a autorisé le plan exhaustif R13 ; il reste une preuve
historique et n'est pas réécrit par l'amendement ultérieur.

Le fractionnement R13 est fermé avant exécution par un registre de 1 180
checkpoints, un par appel externe. Les réservations contiennent de un à huit
appels sans réduire le volume scientifique. Le JSON autoritatif et ses vues
TSV/Markdown indiquent explicitement `PENDING`, `RUNNING`, `COMPLETED` ou
`FAILED`, ainsi que `Fait=OUI/NON`. Les appels interrompus restent reprenables
et les sorties terminées sont liées par taille et SHA-256.

Le premier gel a été déclaré exécutable à tort alors qu'il ne contenait ni
matérialiseur de shard ni worker unitaire. La révision corrective
`R12C-EXECUTABLE-REFREEZE-002` remplace explicitement ce G5 incomplet. Le
matérialiseur couvre les 640 appels SHAPEIT5 et 540 appels IBD, sépare la
vérité privée de la vue aveugle et produit les trois densités comme
sous-ensembles stricts d'un fond de 1 200 000 rangs. Le worker exécute un seul
appel sans shell, mesure ses ressources et ferme atomiquement son checkpoint.

Le modèle haplotypique synthétique est un processus de Markov à trois
populations, avec rafraîchissement
`1 - exp(-delta_cM / 0,05)` et transmissions mendéliennes de 200 trios avec
recombinaison `1 - exp(-delta_cM / 100)`. Les longueurs et cartes utilisées
sont synthétiques et figées ; elles ne constituent pas une estimation de la
carte humaine. Cette approximation doit rester visible dans l'interprétation
scientifique finale.

### R13 — Exécuter la validation aveugle

**État :** validation technique terminée le 13 novembre 2026 ; 25/25 étapes
`PASS`. Qualification scientifique non ouverte.

Exécuter la partition indépendante sans accès à la vérité, clore le run, puis
ouvrir la vérité pour l'évaluation.

Cette exécution était limitée à la portée
`TECHNICAL_BLIND_VALIDATION_ONLY`. Elle peut démontrer la fermeture des 25
étapes, les empreintes, l'isolement des entrées et le routage post-clôture de
la vérité. Elle ne peut ni mesurer les performances scientifiques prévues par
R07/R08, ni conclure à un effet fondateur, ni autoriser une promotion.

Le run `RUN-R13-VALIDATION-BLIND-001` a été clos avant toute ouverture de la
vérité. Le runner a déclaré 26 sorties pour les 25 étapes, sans outil
scientifique externe ni donnée réelle. L'évaluateur figé a ensuite produit 25
évaluations `PASS` de portée `TECHNICAL_INTEGRITY`. Le reçu
`calibration/readiness/R13_TECHNICAL_VALIDATION.json` conserve
`scientific_performance_evaluated=false`, G6 scientifique à `false` et les
trois drapeaux de qualification à `false`.

#### Amendement R13 allégé — intégration des outils scientifiques

**État :** amendement `R13-LIGHTWEIGHT-INTEGRATION-001` clos le 15 novembre
2026 dans la seule portée `R13_TOOL_INTEGRATION_VERIFICATION`.

La validité intrinsèque de SHAPEIT5, Hap-IBD et Refined IBD relève de leurs
publications méthodologiques. La vérification locale porte donc sur leur bon
raccordement au pipeline : identité et version, paramètres et entrées,
contrôles synthétiques appropriés, sorties non vides et intègres, séparation
de la vérité et provenance reproductible. Cette décision ne remplace pas une
mesure de performance R07/R08 et ne qualifie aucune conclusion biologique.

La preuve principale SHAPEIT5 est le checkpoint
`SHAPEIT-DMAX-CHR01-SEED01`, tentative 3 : SHAPEIT5 5.1.1 a terminé avec code
0 sur 4 000 sujets, 3 000 références, 46 524 variants et 200 trios, avec zéro
erreur mendélienne. Le BCF final contient 4 000 échantillons et 46 524
variants ; le BCF, son index CSI, le journal et le reçu sont liés par SHA-256.
L'exécution a duré 21 548,10 secondes et n'a donné au worker aucun accès à la
vérité privée.

L'intégration locale de Hap-IBD et Refined IBD repose sur les contrôles
positifs synthétiques R12 déjà clos : 23/25 duplicats désignés retrouvés et
178 segments pour Hap-IBD ; 22/25 et 122 segments pour Refined IBD. Ces
résultats démontrent que les deux adaptateurs s'exécutent et produisent les
sorties attendues ; ils ne satisfont pas la grille de sensibilité,
spécificité, précision des bornes et concordance prévue pour une
qualification complète.

Le plan historique de 1 180 appels reste intact pour l'audit. Sa disposition
amendée est de 1 appel terminé et retenu, 1 179 appels non exécutés retirés et
zéro nouvel appel autorisé. Un overlay de gouvernance bloque toute nouvelle
réservation ; les lignes `PENDING` ne sont ni effacées ni faussement marquées
`COMPLETED`. Le reçu suivi est
`calibration/readiness/R13_INTEGRATION_AMENDMENT.json`.

Les références externes liées au contrat sont Hofmeister et al. (SHAPEIT5,
DOI `10.1038/s41588-023-01415-w`), Zhou, Browning et Browning (Hap-IBD, DOI
`10.1016/j.ajhg.2020.02.010`) et Browning et Browning (Refined IBD, DOI
`10.1534/genetics.113.150029`). Elles justifient de ne pas répéter une
validation intrinsèque des méthodes ; la preuve d'intégration reste locale.

La vérité `VALIDATION_BLIND` n'est pas ouverte. La performance scientifique
complète reste `NOT_EVALUATED` et les drapeaux `calibrated`,
`founder_effect_proven` et `production_promotion` restent faux.

**Porte de la validation technique historique :** aucun changement de code,
méthode ou seuil entre sa clôture et son évaluation. L'amendement ultérieur
est une décision de portée distincte et traçable.

### R14 — Décider la qualification et produire le rapport final

**État :** clôture technique terminée le 13 novembre 2026 ; qualification
scientifique non accordée et promotion non autorisée.

Décider séparément la qualification de chaque méthode, ses restrictions et sa
durée de validité. Produire ensuite un unique rapport HTML/PDF avec les courbes,
équations, commentaires d'étapes, limites et modifications logicielles.

Le rapport `R14-TECHNICAL-CLOSURE-001` matérialise uniquement les preuves
disponibles : intégrité technique R13 `PASS`, performances scientifiques non
évaluées, blocage de ressources R11 et révision R07 obligatoire. La version
acceptée contient 13 pages A4, quatre figures et huit tables. Le reçu
`calibration/readiness/R14_TECHNICAL_CLOSURE.json` en fixe les artefacts et les
empreintes. Il maintient G6 scientifique et G7 promotion à `false` ainsi que
les trois drapeaux de qualification à `false`.

**Porte :** décision humaine explicite reliée aux versions et empreintes.

### R14b — Clore l'intégration scientifique amendée

**État :** terminé le 15 novembre 2026 après la clôture de l'amendement
`R13-LIGHTWEIGHT-INTEGRATION-001`. Le rendu accepté contient 16 pages A4,
quatre figures, onze tables et aucun espace réservé.

R14b produit un rapport final autonome, détaillé et strictement limité à
la bonne intégration locale des outils. Il consolide comme preuve principale
le checkpoint SHAPEIT5 grandeur réelle terminé, puis comme preuves de soutien
les contrôles positifs Hap-IBD et Refined IBD et les mesures de ressources
R11b. Les graphiques montrent les mesures observées, leurs effectifs et leurs
limites ; les équations employées sont rappelées avec leurs hypothèses et
leurs unités.

Les intervalles binomiaux IBD sont descriptifs : l'indépendance des 25
duplicats désignés n'est pas établie. Le temps SHAPEIT5 est une observation
unique et ne reçoit donc pas d'intervalle de performance inventé. R14b ne
rouvre pas la vérité privée, ne relance aucun outil scientifique et ne
reclasse pas les 1 179 appels retirés comme exécutés.

La décision autorisée est `INTEGRATION_VERIFIED_PERFORMANCE_NOT_EVALUATED`.
Elle maintient G6 scientifique, G7 promotion, `calibrated`,
`founder_effect_proven` et `production_promotion` à `false`.

**Porte :** HTML autonome, PDF A4 paginé, manifeste et contrôles qualité
validés ; reçu suivi liant les empreintes des artefacts finaux.

## 12. Portes globales

| Porte | Autorise | Conditions |
|---|---|---|
| G0 | conception | plan unifié approuvé |
| G1 | création des contrats | besoins, dépendances, panel, vues, dimensionnement et critères validés |
| G2 | micro-panel | schémas, isolation et dimensionnement validés |
| G3 | pilote | micro-panel transversal techniquement complet |
| G4 | gel | pilote interprétable et corrections closes |
| G5 | validation aveugle | code, méthodes, seuils et composants d'exécution scellés ; portée nommée dans le reçu |
| G6 | qualification | validation close sans violation |
| G7 | rapport et promotion éventuelle | décision humaine et traçabilité complète |

Une porte échouée bloque la suivante. Elle ne peut pas être contournée par la
modification manuelle d'un statut.

Le reçu R12b satisfait G5 uniquement pour la portée technique explicitement
nommée `R13_TECHNICAL_BLIND_VALIDATION`. Il ne transforme pas l'échec de la
porte scientifique R11/G4 en succès et n'ouvre pas G6 pour une qualification
scientifique. Un futur reçu de qualification complète exigera d'abord la
révision R07 et la résolution du blocage de ressources.

## 13. Tests et contrôles

### 13.1 Tests unitaires

- déterminisme des graines ;
- transmissions mendéliennes ;
- conversions bp/kb/cM ;
- injection et absence des événements ;
- séparation des partitions et de la vérité ;
- refus des chemins réels ;
- validation des schémas ;
- calcul des métriques et intervalles ;
- censure et cas non évaluables ;
- fidélité du résultat Markdown.

PLINK, KING, bcftools, SHAPEIT5, Hap-IBD et Refined IBD seront simulés dans les
tests unitaires.

### 13.2 Intégration

- petit jeu entièrement synthétique ;
- outils externes réels uniquement dans les tests identifiés ;
- dossiers temporaires ou workspace de calibration ;
- reprise sans écrasement ;
- contrôle des producteurs et signatures ;
- exécution sur 22 autosomes lorsque le contrat l'exige.

### 13.3 Contrôles avant une porte majeure

1. `git diff --check` ;
2. imports du package et de l'orchestrateur ;
3. tests ciblés ;
4. suite complète avant pilote, gel et validation ;
5. vérification des SHA-256 ;
6. audit de fuite de vérité ;
7. contrôle des effectifs indépendants ;
8. revue humaine des critères.

## 14. Analyse des risques

| Risque | Conséquence | Prévention |
|---|---|---|
| panel différent par étape | résultats non composables | release maître et vues dérivées |
| fuite de vérité | performances artificiellement élevées | séparation physique et clôture avant ouverture |
| suradaptation | échec en validation | partitions indépendantes et seuils scellés |
| densité confondue avec pedigree | causalité erronée | plan factoriel sur mêmes sujets |
| simulation trop simple | calibration optimiste | erreurs, structure, nulls et défis prospectifs |
| effectifs dépendants | intervalles trop étroits | unité indépendante déclarée par métrique |
| fenêtre trop courte | censure prise pour recombinaison | dimensionnement en cM et série de défi séparée |
| phase vraie donnée à l'outil | évaluation conditionnelle trompeuse | entrée non phasée et phase vraie privée |
| correction post hoc | biais de validation | version nouvelle et reprise des étapes affectées |
| ressources insuffisantes | qualification inachevée | micro-panel, pilote et porte de ressources |
| rapport final infidèle | conclusion non traçable | valeurs issues uniquement des Markdown clos |

## 15. Modifications du logiciel

Les résultats d'étape peuvent proposer :

- correction d'un bug démontré ;
- clarification d'un contrat ;
- nouveau seuil issu uniquement de la calibration ;
- nouvel identifiant de méthode pour une modification statistique ;
- élargissement d'un intervalle sous-couvrant ;
- meilleure gestion des non-évaluables ;
- amélioration de provenance ou de rendu ;
- test de non-régression supplémentaire.

Une proposition sera inscrite dans le Markdown de l'étape. Elle ne sera
appliquée que dans un lot séparé, avec test ciblé et reprise explicite des
résultats qu'elle invalide. Une amélioration de rendu ne modifiera jamais une
statistique close.

## 16. Rapport final unique

Après R14, le rapport final contiendra :

- architecture et provenance du panel ;
- séparation calibration–validation ;
- résultats de chaque étape ;
- effectifs indépendants et incertitudes ;
- courbes et tables dérivées des sources empreintées ;
- équations réellement utilisées, variables, unités et hypothèses ;
- critères, décisions et domaines non évaluables ;
- analyse des échecs et sensibilités ;
- liste des modifications logicielles réalisées ou encore proposées ;
- restrictions de qualification et conditions de recalibration.

Le rapport sera produit avec `report-formater`, contrôlé visuellement et
exporté en HTML autoportant et PDF. Il ne recalculera aucune statistique.

## 17. Clôture et invalidation

La calibration sera close uniquement si :

1. toutes les étapes requises possèdent un Markdown clos ;
2. la validation aveugle a été exécutée sans modification post hoc ;
3. les décisions sont reliées aux versions et empreintes ;
4. les limites et non-évaluables sont visibles ;
5. le rapport final est fidèle aux sources ;
6. aucune donnée réelle n'a été utilisée pour régler le système ;
7. une décision humaine autorise ou refuse explicitement la promotion.

Un changement de carte, panel, densité, outil, seuil, algorithme, format,
contrat scientifique ou domaine d'usage invalidera tout ou partie de la
qualification selon une matrice de compatibilité à produire en clôture.

## 18. Prochaine action

R12c et G5 sont terminés. L'amendement R13 allégé conserve comme preuve
principale le checkpoint SHAPEIT5 grandeur réelle terminé et les contrôles IBD
R12 ; les 1 179 appels restants sont retirés sans être déclarés exécutés. Le
verrou d'exécution interdit toute nouvelle réservation. Une qualification de
performance R07/R08 complète nécessiterait une nouvelle décision et un nouveau
plan explicites ; elle ne doit pas être inférée de cet amendement.

R14b est terminé et complète l'ancien R14 sans le réinterpréter : il qualifie
l'intégration locale, pas les performances R07/R08, un effet fondateur ou la
production. La prochaine décision éventuelle porte donc soit sur l'archivage de
ce cycle, soit sur un nouveau plan explicite de qualification de performance.
Aucune donnée réelle n'est utilisée et les trois drapeaux de qualification
restent faux.
