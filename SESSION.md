# Continuité de session — Effet_fondateur2

**Date :** 24 novembre 2026

**Branche :** `agent/archive-calibration`

**Base :** `origin/main` à `9dfc0bc`

## État

La campagne de calibration R00–R14b est close. Sa décision finale reste
`COMPLETED_WITH_TOOL_INTEGRATION_VERIFIED` : SHAPEIT5 est vérifié, Hap-IBD et
Refined IBD sont soutenus, mais la performance scientifique exhaustive R07/R08
n'a pas été évaluée. `calibrated`, `founder_effect_proven` et
`production_promotion` restent à `false`.

Le dispositif de calibration a été retiré du pipeline actif et déplacé sous
`archive/calibration/` :

- code historique sous `code/` ;
- contrats et schémas sous `configurations/` ;
- plan et conventions sous `documentation/` ;
- reçus inchangés sous `receipts/` ;
- tests historiques sous `tests/` ;
- manifeste de provenance sous `manifests/`.

Les 1,7 Gio de releases, résultats, vérités synthétiques, journaux, toolchains,
BCF/CSI et rapports locaux ont été déplacés sans suppression sous
`archive/calibration/data-calibration/calibration_workspace/`. Les caches
Python propres à la calibration sont sous `data-calibration/local-caches/`.
Les cinq arborescences temporaires R12c vides ont également été regroupées
sous `data-calibration/tmp-backups/`. Tout ce dossier est ignoré explicitement
par Git.

Le paquet actif ne contient plus `effet_fondateur.calibration` et la CLI V2
n'expose aucun point d'entrée de calibration. Le schéma actif
`schemas/explicit_ibd_calibration.schema.json` reste en place car il appartient
à l'étape scientifique 16C, pas à la campagne R00–R14b.

## Décisions

- Ne pas modifier rétroactivement les reçus historiques ; leurs anciens chemins
  sont expliqués par le manifeste de relocalisation.
- Ne pas exécuter le code archivé ni l'utiliser comme pipeline de production.
- Toute nouvelle calibration exige un nouveau plan et une autorisation
  explicite.
- `websockets`, utilisé uniquement par le rendu R14 archivé, est retiré des
  dépendances actives.
- Aucun commit ni push n'est autorisé sans demande explicite.

## Validation

- `git diff --check` réussi ;
- imports `effet_fondateur` et `effet_fondateur.orchestrator.pipeline`
  réussis ;
- `config/pipeline.example.yaml` valide et aide CLI active vérifiée ;
- aucune importation ni ancien chemin R00–R14b dans le code actif ;
- 20 tests ciblés de configuration et d'orchestration réussis en 7,24 s ;
- suite active complète : 335 tests réussis en 734,89 s ;
- 21 080 fichiers, 1 225 liens relatifs et 237 caches présents sous la
  destination ignorée, sans fichier local suivi par Git ;
- empreintes du rapport R14b et du BCF SHAPEIT5 principal identiques avant et
  après déplacement.

## Prochaine action

Présenter le compte rendu final. Ne pas committer ni pousser sans demande
explicite. Ne pas lancer le pipeline complet, SHAPEIT5, Hap-IBD ou Refined IBD.
