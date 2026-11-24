# Conventions de développement — calibration unifiée

**Version :** 1.0.0 — 12 novembre 2026
**Périmètre :** `panel_master_generic_v1`, étapes R00–R14

Ce document complète le plan unifié sans créer de protocole parallèle. En cas
d'écart, `PLAN_DEVELOPPEMENT_PANEL_SIMULATION.md` prévaut. Les anciens
protocoles D01–D42 ne sont ni une source ni une référence compatible.

## Artefacts et clôture

- Tous les chemins sont relatifs à une racine de release et résolus avant
  lecture. Les chemins absolus, liens symboliques sortants et liens physiques
  sont refusés.
- Chaque artefact consommé ou produit possède un identifiant, un producteur,
  un format, une taille et un SHA-256. Une sortie dynamique est exhaustive dans
  le manifest de son producteur ; aucune découverte par glob n'est admise.
- Une correction ne modifie jamais un run clos. Elle crée un nouvel identifiant
  et rend explicitement caduques les évaluations dépendantes.
- Les données réelles, notamment tout chemin sous `data/input/`, sont interdites
  dans une release de calibration synthétique.

## Aveuglement

- Le processus du pipeline ne peut lire ni `private_truth/` ni `governance/`.
  Il reçoit seulement une vue opaque et son inventaire fermé.
- La vérité n'est ouverte qu'après clôture et seulement par le routeur de
  vérité, pour les familles autorisées à l'évaluateur de l'étape.
- Aucune identité, graine, référence, événement, segment ou contrôle négatif
  ne traverse les partitions `CALIBRATION`, `VALIDATION_BLIND` et
  `RESOURCE_ONLY`.
- Un génotype cible n'est jamais déduit d'un groupe, d'un phénotype ou d'un
  statut clinique.

## Résultats et graphiques

- Chaque résultat d'étape est un Markdown clos ; aucun HTML ou PDF
  intermédiaire n'est créé. Le rendu HTML/PDF unique est réservé à R14.
- Les SVG proviennent uniquement de tables d'évaluation closes et publient
  effectifs indépendants, unités, intervalles, seuils, statut et provenance.
- Les équations réellement utilisées, leurs variables, unités, hypothèses et
  limites sont rappelées. Une règle logique est décrite comme telle et ne reçoit
  pas de formule statistique inventée.
- Une valeur absente est `NA` avec motif et dénominateurs. Elle n'est jamais
  remplacée par zéro, une chaîne vide ou une exclusion silencieuse.
- Aucune donnée génétique individuelle, vérité privée ou précision non mesurée
  n'apparaît dans un Markdown ou une figure.

## Décisions

Les seuls statuts sont `PASS`, `FAIL`, `NOT_EVALUABLE` et
`REVIEW_REQUIRED`, dans l'ordre de décision R08. Une métrique secondaire ne
compense jamais une métrique primaire. Les décisions sont prises avant arrondi
et les seuils de `stage_criteria_v1.yaml` sont liés à l'empreinte de la section
R08 du plan.

Le micro-panel R10 ne peut produire qu'un verdict d'intégrité technique. Il ne
peut établir ni calibration, ni effet fondateur, ni promotion en production.
