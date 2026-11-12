# Continuité de session — Effet_fondateur2

**Date :** 12 novembre 2026

**Branche :** `agent/calibration-d24-d28`

**Décision active :** R01 et R02 sont terminées. La campagne D01–D42 est
abandonnée et ne doit pas être restaurée. La prochaine action autorisée est
R03, la matrice détaillée des besoins des étapes `00` à `19`.

## État décisionnel

- Le plan directeur est
  `calibration/PLAN_DEVELOPPEMENT_PANEL_SIMULATION.md`, version unifiée 1.0 du
  12 novembre 2026, organisée en R00–R14.
- Aucun panel, outil scientifique, pipeline, rapport HTML ou PDF n'a été lancé.
- Le logiciel sous `src/` et ses tests n'ont pas été modifiés.
- Maintenir `calibrated=false`, `founder_effect_proven=false` et
  `production_promotion=false`.

## R01 — remise à zéro terminée

- `calibration/` ne contient plus que le plan unifié.
- 119 anciens fichiers suivis ont été supprimés : 2 fichiers à la racine,
  31 configurations, 16 protocoles, 4 références, 65 schémas et 1 modèle.
- Le fichier non suivi `calibration/AUDIT_PREVISIONNEL_D29_D42.md` a également
  été supprimé. Au total, R01 a supprimé 120 fichiers et 5 sous-répertoires de
  `calibration/`.
- Tout le contenu de `calibration_workspace/` a été supprimé définitivement :
  29 623 fichiers, 1 809 sous-répertoires et 146 entrées de premier niveau.
  Le volume alloué avant suppression était de 9 199 776 Kio, soit environ
  8,8 Gio ; le workspace occupe maintenant 0 Kio et est vide.
- Les reliquats en lecture seule de deux caches D23 ont été rendus inscriptibles
  uniquement dans `calibration_workspace/`, puis supprimés.
- Aucune archive, copie, branche, stash, tag ou réécriture d'historique n'a été
  créé. Aucun chemin sous `src/`, `tests/`, `docs/`, `config/`, `data/` ou
  `simulation_genotype_famille/` n'a été supprimé ou modifié.

## R02 — intégrité du logiciel vérifiée

- `git diff --check` : réussi.
- Imports `effet_fondateur` et `effet_fondateur.orchestrator.pipeline` :
  réussis.
- Tests ciblés de configuration, démarrage et orchestration : 20 réussis en
  6,21 s.
- Suite logicielle de référence, exécutée avec
  `--ignore=tests/calibration` : 334 tests réussis en 723,04 s.
- Les 62 modules de `tests/calibration/` ont été exclus de la référence car ils
  appartiennent à l'ancien contrat de calibration ; 41 référencent directement
  les configurations, protocoles, références, schémas ou modèles supprimés.
  Ils sont classés obsolètes par dépendance, sans restauration ni modification,
  et seront réévalués/adaptés à R09 après conception du nouveau panel.
- Aucun échec logiciel ou régression imputable au nettoyage n'a été observé.

## État Git à préserver

- `SESSION.md` et le plan unifié restent modifiés.
- Les 119 suppressions suivies sous `calibration/` sont attendues.
- Le fichier d'audit non suivi a disparu et `calibration_workspace/` reste
  vide et ignoré par Git.
- Aucun commit, push ou merge n'est autorisé sans nouvelle demande.

## Prochaine action

Exécuter R03 uniquement après nouvelle demande : compléter, pour chaque étape
`00` à `19`, les formats, unités, vérités, tailles, outils, ressources, erreurs
attendues et données graphiques. Ne pas commencer R04 ni générer de panel avant
validation de la matrice R03.
