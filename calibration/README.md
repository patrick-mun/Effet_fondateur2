# Calibration par simulation

Ce dossier contient uniquement les éléments reproductibles qui doivent être
versionnés avec le pipeline : plan, protocoles, configurations, schémas et
modèles de rapports.

Les développements suivent obligatoirement
`CONVENTIONS_DEVELOPPEMENT.md`, déclinaison locale de `AGENTS.md` pour la
simulation, l'aveuglement, les tests et le calcul scientifique.

## Séparation des espaces

- `calibration/` : documents et contrats suivis par Git ;
- `src/effet_fondateur/calibration/` : futur code de simulation et
  d'évaluation suivi par Git ;
- `tests/calibration/` : futurs tests et minuscules fixtures synthétiques ;
- `calibration_workspace/` : panels, vérités, runs, figures et rapports
  générés, entièrement ignorés par Git.

Les chemins historiques éventuels `calibration/runs/` et
`calibration/local/` sont également ignorés par précaution.

## Règles de sécurité

- Ne jamais copier de donnée de `data/input/dock6_real/` dans ces espaces.
- Ne jamais utiliser un résultat individuel DOCK6 pour régler un seuil.
- Ne jamais committer un VCF, BCF, fichier de vérité ou rapport généré.
- Fournir le chemin du workspace au lancement ; ne pas enregistrer de chemin
  absolu local dans une configuration versionnée.
- Vérifier les graines, versions et SHA-256 avant toute qualification.

## État actuel

Le périmètre initial D01 et l'arborescence isolée D02 sont établis. Le protocole
D03 est `APPROVED` et les critères D04 sont `APPROVED_PRE_PILOT`. Cette
approbation porte sur le cadre et les règles : les seuils de performance restent
non scellés et seront proposés sur le seul pilote de calibration selon une
règle écrite avant le run, puis figés avant le panel aveugle. Le schéma de
configuration D05 et le manifest de vérité D06 sont créés et testés. Le
manifest D06 impose neuf familles de vérité synthétique, leurs champs, unités,
empreintes, chemins privés et consommateurs autorisés. Aucun simulateur, panel
ou rapport de calibration n'a encore été produit. Le développement doit suivre
`PLAN_DEVELOPPEMENT_PANEL_SIMULATION.md` et ses portes de passage.
