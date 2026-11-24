# Manifeste de relocalisation de la calibration

Date de relocalisation : 24 novembre 2026 (Indian/Reunion)

Base Git : `origin/main` à `9dfc0bc`

Les tailles ci-dessous sont les tailles logiques relevées avant déplacement.
Les déplacements suivis ont été réalisés avec `git mv`. Les répertoires ne
possèdent pas d'empreinte SHA-256 portable ; Git conserve l'identité des
fichiers suivis et les artefacts importants sont empreintés séparément.

| Ancien chemin | Nouveau chemin | Taille | Statut | Justification |
|---|---|---:|---|---|
| `calibration/` | `archive/calibration/{configurations,documentation,receipts}/` | 273 420 octets, 37 fichiers | suivi | contrats, plan et reçus R00–R14b |
| `src/effet_fondateur/calibration/` | `archive/calibration/code/src/effet_fondateur/calibration/` | 2 237 293 octets, 121 fichiers | suivi | code exclusivement historique |
| `tests/calibration/` | `archive/calibration/tests/` | 142 644 octets, 28 fichiers | suivi | tests retirés de la suite active |
| `calibration_workspace/` | `archive/calibration/data-calibration/calibration_workspace/` | 1 807 334 572 octets, 21 080 fichiers et 1 225 liens | ignoré | releases, preuves, rapports et outils locaux |
| caches Python des sources et tests de calibration | `archive/calibration/data-calibration/local-caches/` | 5 085 152 octets, 237 fichiers | ignoré | conserver sans laisser de reliquat actif |
| cinq sauvegardes R12c sous `/private/tmp/` | `archive/calibration/data-calibration/tmp-backups/` | 0 octet, arborescences vides | ignoré | regrouper tous les reliquats sans suppression |

## Empreintes d'artefacts importants

| Ancien artefact | Nouveau chemin relatif à `data-calibration/` | SHA-256 |
|---|---|---|
| `calibration_workspace/r14b_integration_closure/R14B_DECISION.json` | `calibration_workspace/r14b_integration_closure/R14B_DECISION.json` | `8040f33f74d0e6aeeaca13e2b805c3344393c6256fd76c4fee5467423cd5f922` |
| `calibration_workspace/r14b_integration_closure/R14B_REPORT_MANIFEST.json` | `calibration_workspace/r14b_integration_closure/R14B_REPORT_MANIFEST.json` | `6a490048a97eb45f1f1cbd81bbb01cb6b24267c4aa1a693a230cc346921e044b` |
| `calibration_workspace/r14b_integration_closure/rapport_final_r14b.html` | `calibration_workspace/r14b_integration_closure/rapport_final_r14b.html` | `a5d29e94fd0db71c14a7429ed9214144d2a26ec69fccb549c558169fbab34025` |
| `calibration_workspace/r14b_integration_closure/rapport_final_r14b.pdf` | `calibration_workspace/r14b_integration_closure/rapport_final_r14b.pdf` | `c1c092c04ae9e6241ea59cf530b8e4eac658c58ab798ad66122863f5805ad3bb` |
| `calibration_workspace/r13_scientific/RUN-R13-SCIENTIFIC-001/outputs/external_calls/SHAPEIT-DMAX-CHR01-SEED01/attempt-03/phased.bcf` | `calibration_workspace/r13_scientific/RUN-R13-SCIENTIFIC-001/outputs/external_calls/SHAPEIT-DMAX-CHR01-SEED01/attempt-03/phased.bcf` | `a046db99337a561b9e9598da42cf1d4277296115edf396858c302ea3e2967cc8` |

Les fichiers de `data-calibration/` sont locaux et absents de Git. Ce
manifeste ne contient aucun génotype, identifiant individuel ou secret.
