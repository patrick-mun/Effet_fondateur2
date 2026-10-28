# Sous-système de calibration

Ce dossier contient le code de configuration, d'isolation, de simulation, de
vérité, d'évaluation et de production des rapports de calibration.

`isolation.py`, introduit à l'étape D08, contrôle les chemins physiques, liens
symboliques, références de manifeste et collisions entre panels avant toute
lecture d'un artefact. `partitioning.py`, introduit à l'étape D10, répartit de
façon déterministe les groupes indépendants entre donneurs de calibration,
donneurs de validation et référence de phasage commune. Il contrôle la
couverture de la source, les intersections, les empreintes et la confidentialité
du manifest sans écrire d'identifiant.

`materialize_partitions.py` vérifie les empreintes des métadonnées D09, refuse
l'écrasement d'un run et écrit atomiquement les trois listes privées ainsi que
leur manifest agrégé dans `calibration_workspace/`. Sa sortie standard ne
contient que des effectifs et statuts de contrôle.
