# Registre des décisions de calibration

Ce registre conserve les décisions humaines qui modifient le périmètre, les
critères ou l'interprétation du projet de calibration. Une décision ne vaut ni
résultat scientifique ni qualification automatique d'une étape du pipeline.

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
