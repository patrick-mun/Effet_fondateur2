# Registre des décisions de calibration

Ce registre conserve les décisions humaines qui modifient le périmètre, les
critères ou l'interprétation du projet de calibration. Une décision ne vaut ni
résultat scientifique ni qualification automatique d'une étape du pipeline.

## CAL-DEC-016 — Exécution bornée de D28

Le 10 novembre 2026, l'utilisateur demande explicitement de passer à D28.
Statut : `DEMANDEE_EXPLICITEMENT`. Cette décision autorise la
préspecification, les appels aveugles du cœur d'enrichissement 16B sur le banc
synthétique D24, l'évaluation après ouverture de la vérité et le rapport
CAL-REP-16B.

L'autorisation reste bornée à `calibration_workspace/`, hors `dock6_real`,
cache interne et `data/input`. Les probabilités doivent être scellées avant la
vérité. Le seuil alpha du banc sert uniquement à mesurer faux positifs et
puissance ; le seuil de classification de production reste nul.

Les seuils d'acceptation CAL-ENRICH-002/003 restent non scellés. Le banc à
trois familles et fond synthétique n'est pas une validation bout-en-bout de la
référence 1000 Genomes. `calibrated=false`, `founder_effect_proven=false` et
`production_promotion=false` restent imposés.

## CAL-DEC-015 — Exécution bornée de D27

Le 10 novembre 2026, l'utilisateur demande explicitement de passer à D27.
Statut : `DEMANDEE_EXPLICITEMENT`. Cette décision autorise la
préspecification, l'appel aveugle du cœur PCA/projection de 16A, l'évaluation
sur les vérités D12/D14 et le rapport CAL-REP-16A.

L'autorisation reste bornée à `calibration_workspace/`, hors `dock6_real`,
cache interne et `data/input`. La référence D10 doit être disjointe des
donneurs simulés et les coordonnées doivent être scellées avant ouverture de
la vérité. L'alignement des profils ne peut utiliser que les références.

Le banc régional ne vaut pas validation genome-wide et la phase candidate D17
ne vaut pas exécution bout-en-bout de SHAPEIT5 puis 16A. Aucun code de strate
ne peut être interprété comme une identité géographique ou ethnique. Les
seuils CAL-ANC-001/002 restent non scellés ; `calibrated=false` et
`founder_effect_proven=false` restent imposés.

## CAL-DEC-014 — D25b puis D26

Le 10 novembre 2026, l'utilisateur demande explicitement « faisons le D25B,
puis le D26 ». Statut : `DEMANDEE_EXPLICITEMENT`. Cette décision autorise
d'abord la préspecification, le développement, la simulation indépendante,
l'évaluation close et le rapport de D25b. D26 ne commence qu'après fermeture
et interprétation de D25b.

D25b doit créer une nouvelle méthode de vraisemblance sans modifier V1/V2 et
ne peut promouvoir automatiquement une méthode en production. L'autorisation
reste limitée aux données synthétiques sous `calibration_workspace/`, hors
`dock6_real`, cache interne et pipeline archivé. Le passage ultérieur à D26
fera l'objet de son propre protocole fermé avant exécution.

D26 est ensuite exécutée sous
`D26_LD_ROH_PRESPECIFICATION.md`. Les appels observés sont scellés avant la
référence dense. Le statut technique est `PASS`, mais les critères LD/ROH
`pilot_then_seal` restent `REVIEW_REQUIRED` faute de seuil antérieur ; aucun
seuil n'est ajusté sur les résultats. Les deux rapports indépendants sont
`CAL-REP-15` et `CAL-REP-16`.

`calibrated=false` et `founder_effect_proven=false` restent imposés pendant
les deux étapes.

## CAL-DEC-013 — Exécution bornée de l'audit D25a

Le 10 novembre 2026, l'utilisateur demande explicitement « fais la D25a ».
Statut : `DEMANDEE_EXPLICITEMENT`. Cette décision autorise la
préspecification, la simulation indépendante, la vérification contre le script
R officiel, l'évaluation close et le rapport CAL-REP-14A.

L'autorisation reste bornée à `calibration_workspace/` et exclut les données
`dock6_real`, le cache interne, le pipeline archivé et toute modification
rétrospective des résultats D25. La méthode V1 doit rester disponible ; une
candidate distincte ne peut la remplacer qu'après succès des critères
préspecifiés de fidélité et de performance.

Un `PASS` de parité R/Python ne vaut pas qualification de l'observation
tronquée ni du modèle corrélé. `calibrated=false`,
`founder_effect_proven=false`, les faux positifs D24 et la non-évaluation des
erreurs SHAPEIT5 sont maintenus.

## CAL-DEC-012 — Exécution bornée de D25

Le 10 novembre 2026, l'utilisateur demande explicitement « passons à D25 ».
Statut : `DEMANDEE_EXPLICITEMENT`. Cette décision autorise la préspecification,
la matérialisation puis l'évaluation close de la datation Gamma sur la grille
minimale complète, la reprise descriptive des sorties positives D24 et le
rapport CAL-REP-14.

L'autorisation reste bornée à `calibration_workspace/` et exclut
`dock6_real`, le cache interne et le pipeline archivé. La vérité d'âge doit
rester fermée jusqu'à la production des estimations. Les erreurs SHAPEIT5 ne
peuvent pas être inventées : D24 utilisant le phasage-vérité, leur effet reste
`NOT_EVALUABLE` sans sortie appariée dédiée.

La conformité technique ne scelle aucun seuil `pilot_then_seal`. D25 évalue
un estimateur conditionnel à une origine unique ; elle ne prouve ni partage
IBD ni effet fondateur et ne corrige pas les faux positifs D24.
`calibrated=false`, `founder_effect_proven=false` et le problème 07 différé
sont maintenus.

## CAL-DEC-011 — Exécution bornée de D24

Le 10 novembre 2026, l'utilisateur demande explicitement « lance le D24 ».
Statut : `DEMANDEE_EXPLICITEMENT`. Cette décision autorise la matérialisation
du panel synthétique D13–D14 à trois familles, des cinq scénarios positifs et
cinq négatifs préspecifiés, leur projection sur les quatre profils, l'appel du
moteur réel de l'étape 13, l'évaluation close et le rapport CAL-REP-13.

L'autorisation reste bornée à `calibration_workspace/`. Elle exclut toute
donnée `dock6_real`, tout cache interne, tout panel de validation et toute
modification rétrospective des critères. Le phasage-vérité D14 est utilisé
pour isoler l'étape 13 ; l'exécution ne prétend donc pas mesurer la propagation
des erreurs de phasage D23/D23b.

Un `PASS` des contrats techniques ne vaut pas qualification scientifique. Les
critères `pilot_then_seal` restent non scellés, les faux positifs sont conservés
et aucune optimisation sur ces mêmes vérités ne pourra servir de validation.
`calibrated=false`, `founder_effect_proven=false` et le problème 07 différé
sont maintenus.

## CAL-DEC-010 — Préparation de D24 sans exécution

Le 10 novembre 2026, l'utilisateur demande explicitement de préparer D24.
Statut : `DEMANDEE_EXPLICITEMENT`. Cette décision autorise la
préspecification de l'évaluation de l'étape 13, ses configurations, schémas,
oracles, évaluateur agrégé, contrat de rapport et tests. Elle n'autorise ni la
matérialisation du nouveau panel, ni l'exécution du pipeline, ni l'ouverture
d'un panel de validation.

L'audit d'entrée conclut que D14 ne fournit que deux unités familiales
indépendantes, que D15 décrit un partage IBD paire-à-paire dont le fond reste
inconnu et que D23a fournit une vérité de phase, non une vérité de consensus
fondateur. D24 requiert donc un panel propre comportant au moins trois unités
familiales indépendantes par scénario positif, de vrais scénarios négatifs et
des bornes fondatrices privées indépendantes des sorties de l'étape 13.

Les critères `CAL-FOUND-001` à `CAL-FOUND-004` restent
`pilot_then_seal`. Seul `CAL-FOUND-005` conserve son seuil fixe : aucune
reclassification IBS→IBD sans preuve IBD explicite. La préparation maintient
`calibrated=false`, `founder_effect_proven=false` et le problème scientifique
07 différé. Aucun résultat D23/D23b ne qualifie D24, et aucune donnée
`dock6_real` ni aucun cache interne n'est autorisé.

## CAL-DEC-009 — Investigation D23b densité × pedigree

Le 9 novembre 2026, après revue des résultats D23, l'utilisateur demande
explicitement d'exécuter les parties 1 à 3 de la reprise proposée :
préspecification d'une comparaison contrôlée, détermination du maximum de
marqueurs exploitable et utilisation de ce maximum comme parcours primaire,
avant de reprendre séparément les parties 4 et 5 relatives à D24. Statut :
`DEMANDEE_EXPLICITEMENT`.

L'autorisation porte sur le protocole
`D23B_DENSITY_PEDIGREE_PRESPECIFICATION.md`, une nouvelle préparation issue
des mêmes donneurs, sujets, transmissions, références et paramètres D23a, les
phasages bornés 10A/11/12/12A nécessaires, leur évaluateur et leur rapport
indépendant. Les profils de densité sont emboîtés et le facteur pedigree est
manipulé uniquement dans le parcours régional qui le prend réellement en
charge. La densité maximale reste primaire quel que soit le résultat.

Cette décision n'autorise aucune donnée `dock6_real`, aucun cache interne de
l'étude, aucune préparation ou exécution D24, aucun réglage sur le panel de
validation, aucun seuil scellé et aucun passage à `calibrated=true` ou
`founder_effect_proven=true`. Le FAIL 07 reste différé. Les nouveaux artefacts
générés restent sous `calibration_workspace/`; aucun commit, push ou merge
n'est autorisé par cette demande.

## CAL-DEC-008 — Exécution bornée D23

Le 30 octobre 2026, l'utilisateur demande explicitement : « lance l'étape
D23 ». Cette demande autorise le banc de calibration séparé 10A/11/12/12A,
ses évaluateurs et ses quatre rapports indépendants, conformément au plan.
Statut : `DEMANDEE_EXPLICITEMENT`.

Le protocole D23 est fixé avant les performances. D23a v2 et la référence D10
restent les seules sources autorisées ; ni le panel aveugle synthétique de
validation, ni l'étude réelle ne sont ouverts. Le défaut scientifique 07 reste
reporté. Cette décision n'approuve aucun résultat, aucune exclusion d'étude,
aucun seuil `pilot_then_seal` et aucun passage à `calibrated=true`.

## CAL-DEC-007 — Préparation D23a du phasage autosomique

Le 30 octobre 2026, l'utilisateur approuve explicitement la préparation D23a
proposée avant D23 : panel sur 22 autosomes lié à la puce, raccord de la liste
de référence D10 aux parcours 12/12A, métadonnées synthétiques observables et
contrôles d'absence de fuite. Statut : `APPROUVEE`.

Le nouveau banc borné et ses règles sont décrits dans
`D23A_PHASING_PREPARATION.md` et `d23a_phasing_panel_v2.json` (v1 conservée
comme brouillon interdit d'usage, graine de phasage séparée en v2). L'annotation
technique D19 est réutilisée ; les génotypes sont explicitement simulés,
sans inférence clinique. Les sources restent publiques, hors cache interne,
les données générées dans le workspace ignoré. Le panel de validation n'est
pas ouvert et les panels existants ne sont pas modifiés.

Cette décision n'approuve ni une cohorte d'étude, ni des exclusions 07, ni des
seuils de performance. Aucun run complet, passage à `calibrated=true` ou
qualification des étapes 10A/11/12/12A n'en découle. Un rapport de préparation
indépendant doit distinguer ressources disponibles et performances non mesurées.

## CAL-DEC-001 — Périmètre de la première tranche

| Champ | Valeur |
|---|---|
| Date | 28 octobre 2026 |
| Étape de développement | `D01` |
| Statut | `APPROUVEE` |
| Effet | immédiat |
| Porte concernée | `G0` |

### Décision

La première tranche de qualification scientifique approfondie porte sur les
étapes suivantes :

- `12` — qualité du phasage autour de la variation cible ;
- `13` — délimitation de l'haplotype fondateur ;
- `14` — estimation de l'âge de l'haplotype ;
- `16C` — détection IBD explicite et calibration de ses scénarios.

Cette priorité définit l'ordre de développement des scénarios scientifiques et
des évaluateurs. Elle ne change pas l'ordre d'exécution du pipeline et ne rend
aucune étape amont facultative.

### Traitement des autres étapes

Toutes les étapes activées de `00` à `19` restent dans le périmètre de
validation du pipeline :

1. chaque étape reçoit un contrôle technique de ses entrées, sorties, effectifs,
   formats, statuts et empreintes ;
2. chaque comportement scientifique testable reçoit au minimum un scénario
   fonctionnel à vérité connue ;
3. chaque étape fait l'objet du rapport indépendant prévu par le plan ;
4. une étape non évaluée est déclarée `NOT_EVALUABLE` ou non qualifiée ;
5. aucune étape n'est considérée comme validée du seul fait que `12`, `13`,
   `14` ou `16C` réussit.

Les étapes `15`, `16`, `16A` et `16B` seront qualifiées séparément avant que
leurs résultats puissent renforcer une conclusion scientifique. Les étapes
`17` à `19` devront démontrer qu'elles restituent les résultats sources sans
inventer, masquer ou recalculer une conclusion.

### Justification

Les quatre étapes prioritaires forment la chaîne la plus directement liée à
l'hypothèse d'effet fondateur : qualité de la phase, récupération du segment
partagé, datation sur une longueur estimée et confirmation indépendante par
IBD explicite. Une erreur à l'une de ces étapes peut modifier directement
l'interprétation scientifique.

Les autres étapes restent nécessaires pour prévenir les erreurs de format, de
qualité, de parenté, de structure de population, de LD, de ROH, d'ascendance ou
de restitution. Leur validation est donc obligatoire, avec une profondeur
adaptée à leur rôle.

### Hors périmètre de cette décision

Cette décision :

- n'autorise encore aucun développement de simulateur ;
- n'autorise aucune exécution de panel ;
- ne fixe aucun seuil scientifique ;
- ne qualifie aucune méthode existante ;
- ne permet pas de passer `16C` à `calibrated: true` ;
- ne permet pas d'utiliser des données ou résultats individuels de
  `dock6_real`.

### Conditions de révision

Une nouvelle décision est requise si une étape prioritaire est ajoutée ou
retirée, si une étape activée est dispensée de rapport, ou si le niveau de
validation prévu est réduit. La nouvelle entrée devra référencer
`CAL-DEC-001` sans modifier rétroactivement son contenu.

### Validation humaine

Décision approuvée explicitement par le responsable du projet au cours de la
session du 28 octobre 2026, après examen de la différence entre priorité
scientifique et validation complète du pipeline.

## CAL-DEC-002 — Approbation de D03 et D04 avant pilote

| Champ | Valeur |
|---|---|
| Date | 28 octobre 2026 |
| Étapes de développement | `D03`, `D04` |
| Statut | `APPROUVEE` |
| Effet | immédiat |
| Porte concernée | préparation de `G1` |

### Décision

Le protocole général D03 version `1.0.0` est approuvé avec le statut
`APPROVED`. Les critères d'acceptation D04 version `0.2.0` sont approuvés avec
le statut `APPROVED_PRE_PILOT`.

L'approbation de D04 porte sur les métriques, unités, sens d'amélioration,
règles de valeurs manquantes, méthodes d'incertitude et cycle de décision. Les
critères `pilot_then_seal` restent volontairement sans valeur numérique.

### Empreintes au moment de l'approbation

- `PROTOCOLE_GENERAL.md` :
  `108e8b877091f44e26ad0a8637f911388af559da86680ff685e43569ecfb7062` ;
- `CRITERES_ACCEPTATION.yaml` :
  `90f65e9a685fa388bcab813603c6b037daeb04badf6cf2024a2bbd6187c1cec2`.

### Portée et restrictions

Cette décision autorise la poursuite du développement documentaire, des
schémas, des gardes d'isolation et des tests unitaires. Elle :

- ne scelle aucun seuil de performance issu du pilote ;
- n'autorise pas encore le micro-panel, le pilote ou la qualification ;
- ne qualifie aucune méthode scientifique existante ;
- ne permet pas de passer `16C` à `calibrated: true` ;
- ne permet aucune utilisation de donnée ou résultat individuel de
  `dock6_real`.

Les autorisations de run restent soumises aux portes G1 à G4 du plan.

### Validation humaine

Décision approuvée explicitement par le responsable du projet dans la session
du 28 octobre 2026 après lecture du protocole D03 et de la révision D04.

## CAL-DEC-003 — Autorisation du diagnostic technique 03–04 de D19

Date : 30 octobre 2026. Portée : intégration régionale synthétique, hors pilote
et qualification. Autorité : demande explicite de l'utilisateur « traite c'est
de point pour qu'il devienne evaluable », après présentation de l'approbation
du registre technique, de PLINK et des scénarios cible absente/présente.

Cette autorisation est enregistrée pour les nouveaux registres synthétiques
strictement identiques aux registres neutres D18, après contrôle de leurs
identités, sources, effectifs et ordres. Elle ne prétend ni à une revue clinique
ni à une validation des familles. Le rôle enregistré est `calibration_owner`,
sans nom ou identifiant personnel. Toute modification de ces registres ou du
périmètre exige une nouvelle décision ; les anciens runs restent immuables.

Les nouveaux runs sont limités aux étapes 00–04 avec PLINK. Le scénario
principal conserve toutes les observations D18 et utilise une cible technique
publique documentée absente de leurs sondes. Ses génotypes sont synthétiques,
équilibrés par une règle déterministe explicite sans statut clinique. La
confirmation de référence et HGVS concerne seulement l'annotation publique
ClinVar, jamais une mesure sur les sujets simulés. Un scénario négatif dérivé
ajoute cette cible à une copie des observations et doit être refusé par 04.

Aucune donnée d'étude, aucun nouveau seuil, aucune interprétation clinique,
aucune qualification ni passage à `calibrated=true` n'est autorisé ici.

## CAL-DEC-004 — Diagnostic de composants QC 05 / 10 de D20

Date : 30 octobre 2026. Autorité : demande explicite de l'utilisateur de
passer à D20 en vérifiant les manques pour rendre les étapes évaluables.
Portée : nouveau registre neutre entièrement synthétique du banc D20,
intégration 00–05 et appel isolé de 10 avec fixtures contractuelles de
cohortes préspecifiées. Ces cohortes ne sont pas des résultats de 09.

La construction mathématique IID, les défauts, la graine et les attentes
sont spécifiés dans `D20_QC_PRESPECIFICATION.md` et la configuration D20 avant
le run. Cette décision n'étend pas CAL-DEC-003 au panel D18 modifié ; elle
concerne ce nouveau banc de composants et sa revue technique de registre,
sans revue clinique ou validation de familles prétendue. Aucun pilote,
qualification, changement des seuils de production ou passage à
`calibrated=true` n'est autorisé. Les étapes 06–09 restent hors exécution.

## CAL-DEC-005 — Diagnostic D21 de parenté 06 / 07

Date : 30 octobre 2026. Autorité : demande explicite de passer à D21 en
vérifiant les prérequis comme pour D20. Portée : nouveau banc entièrement
synthétique, son registre technique et ses pedigrees générés connus ;
exécution bornée 00–03 puis 05–07 avec PLINK et KING. Le protocole D21 et
sa configuration fixent les cas avant les runs.

Cette décision n'approuve aucune exclusion issue de 07. Le besoin de revue
`kinship_exclusion_approval` reste en attente et aucune cohorte 09 n'est gelée.
Aucune donnée d'étude, aucun seuil réglé après résultats, aucun pilote ou
qualification ni passage à `calibrated=true` n'est autorisé.

## CAL-DEC-006 — Diagnostic D22 de structure et gel 08 / 09

Date : 30 octobre 2026. Autorité : l'utilisateur reporte l'étude du FAIL 07
après l'ensemble de l'évaluation et demande l'étape suivante, D22.
Portée : nouveau registre synthétique, 00–03 puis composants 08 et 09 isolés.
Les contrats amont non exécutés et les revues d'exclusion sont des fixtures
techniques déclarées. Leurs règles sont fixées avant le run ; leur rôle est
`calibration_contract_fixture`, jamais une revue scientifique humaine simulée.
Les empreintes des sorties peuvent être liées après clôture pour tester 09,
sans ouvrir les strates vraies ou choisir les exclusions selon la réussite.

Cette décision n'approuve aucune exclusion de D21 ou de l'étude, ne corrige
pas 07 et ne qualifie aucune cohorte réelle. Le FAIL 07 reste ouvert pour la
revue globale. Aucun seuil changé après résultats, aucune donnée réelle,
aucun pilote, aucun `calibrated=true`.
