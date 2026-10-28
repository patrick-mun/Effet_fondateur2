# Calibration par simulation

Ce dossier contient uniquement les éléments reproductibles qui doivent être
versionnés avec le pipeline : plan, protocoles, configurations, schémas et
modèles de rapports.

Les développements suivent obligatoirement
`CONVENTIONS_DEVELOPPEMENT.md`, déclinaison locale de `AGENTS.md` pour la
simulation, l'aveuglement, les tests et le calcul scientifique.

## Séparation des espaces

- `calibration/` : documents et contrats suivis par Git ;
- `calibration/references/` : liste fermée des références publiques D09 ;
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

## Matérialisation D10

La commande suivante construit un nouveau jeu de partitions à partir des deux
métadonnées publiques vérifiées par D09 :

```bash
.venv/bin/python -m effet_fondateur.calibration.materialize_partitions \
  --public-catalog calibration/references/public_reference_catalog.json \
  --population-metadata <fichier_public_population_3202> \
  --unrelated-index <index_public_2504> \
  --workspace-root calibration_workspace \
  --partition-run-id <identifiant_nouveau>
```

Le dossier cible doit être neuf. La commande n'affiche aucun identifiant et
écrit les listes avec des permissions privées sous
`private_truth/partitions/`.

## État actuel

Le périmètre initial D01 et l'arborescence isolée D02 sont établis. Le protocole
D03 est `APPROVED` et les critères D04 sont `APPROVED_PRE_PILOT`. Cette
approbation porte sur le cadre et les règles : les seuils de performance restent
non scellés et seront proposés sur le seul pilote de calibration selon une
règle écrite avant le run, puis figés avant le panel aveugle. Le schéma de
configuration D05, le manifest de vérité D06 et le manifest reproductible de
run D07 sont créés et testés. D06 impose neuf familles de vérité synthétique,
leurs champs, unités, empreintes, chemins privés et consommateurs autorisés.
D07 impose les sources, versions, commandes, graines, fichiers, empreintes et
instructions de rejeu nécessaires pour reproduire une exécution sans
connaissance orale. Les gardes D08 contrôlent les chemins physiques, liens
symboliques ou physiques, références croisées et collisions entre panels avant
toute lecture d'un artefact. D09 fige quatre sources publiques : panel phasé
1kGP, métadonnées de population, index non apparenté et cartes GRCh38
SHAPEIT4/HapMap. Chaque source possède version, origine, assemblage, conditions
d'usage et empreinte SHA-256 ; les charges utiles acquises ultérieurement
devront aussi recevoir un SHA-256 local avant usage. D10 ajoute le contrat et
le validateur des trois partitions : donneurs de calibration, donneurs de
validation et référence de phasage commune. Les groupes de pedigree ou de
parenté ne peuvent pas traverser deux partitions et les listes de membres
restent privées dans le workspace. Aucun simulateur, panel ou rapport de
calibration n'a encore été produit. Une matérialisation D10 contrôlée existe
sous `calibration_workspace/d10_1kg_3202_20220422_v1/` ; elle contient 900
donneurs de calibration, 899 donneurs de validation et 1 403 références de
phasage, sans intersection d'identifiant ni de groupe. Le développement doit
suivre `PLAN_DEVELOPPEMENT_PANEL_SIMULATION.md` et ses portes de passage.
