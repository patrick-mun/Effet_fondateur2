# Suivi de session

Dernière mise à jour : 4 août 2026

## État actuel

- Dépôt GitHub : `patrick-mun/Effet_fondateur2`
- Branche de travail : `agent/restore-project-environment`
- Pull request : `#1` — `Restore project environment and documentation`
- Branche cible : `main`
- Environnement local : Python 3.12.13 dans `.venv`
- PLINK installé : 1.9.0-b.7.11
- KING, Gamma, Rscript et les packages R requis sont disponibles.

## Travail terminé

- Inventaire complet du projet et du rôle de ses scripts.
- Suppression des caches Python et fichiers `.DS_Store`.
- Mise à jour de `.gitignore`.
- Création et installation d'un environnement Python reproductible.
- Ajout des dépendances réelles dans `requirements.txt`.
- Compatibilité macOS Intel stabilisée avec NumPy 1.26.4 et PyArrow 14.0.2.
- Réécriture du README avec installation, activation et commandes actuelles.
- Correction de l'interface pour relancer le pipeline avec `sys.executable`.
- Correction de `run_all_tests.sh` pour utiliser `.venv`.
- Installation et test réel de PLINK sur les données PED/MAP existantes.

## Validations effectuées

- `pip check` : aucune dépendance cassée.
- Import de `run_pipeline` : réussi.
- Imports Python principaux : réussis.
- Packages R `adegenet`, `ggplot2`, `ape`, `ade4` et `poppr` : disponibles.
- Test PLINK : 2 726 variants et 127 individus lus correctement.
- Tests Python : 12 réussis, 4 défaillants.

## Problèmes connus

1. `data/input/complex_simulation/cas.txt` et `temoins.txt` sont absents.
2. Le test Gamma appelle une confirmation interactive sans `auto_confirm=True`.
3. Le test d'âge Gamma n'envoie pas l'argument obligatoire `n`.
4. Deux tests de prétraitement simulent PLINK sans créer `filtered_data.ped`.
5. Les fichiers uploadés par Streamlit sont enregistrés sous `user_input.*`,
   mais le pipeline lit encore `genotype_data.*`.
6. La datation Gamma reçoit actuellement des positions de SNP séparées en deux
   moitiés, alors que le modèle attend des longueurs gauche/droite pertinentes
   autour de la mutation.
7. Le README et les modules Wiki peuvent encore contenir quelques descriptions
   scientifiques historiques à confronter au comportement réel.

## Prochaines étapes recommandées

1. Générer automatiquement `cas.txt` et `temoins.txt` depuis le PED et les
   groupes, avec validation des identifiants PLINK.
2. Actualiser les quatre tests défaillants et obtenir une suite entièrement
   verte.
3. Raccorder les uploads Streamlit aux fichiers réellement consommés.
4. Revoir les entrées mathématiques de l'estimation Gamma.
5. Exécuter le pipeline complet sur une copie contrôlée des données.
6. Examiner le rapport produit avant de fusionner la PR.

## Décisions à conserver

- Utiliser le nom standard `AGENTS.md` pour les instructions Codex.
- Conserver les deux convertisseurs : ACPA vers VCF, puis VCF vers MAP/allèles.
- Ne pas supprimer les rapports historiques, le binaire KING embarqué ou les
  anciens scripts avant une revue dédiée.
- Garder la PR en brouillon tant que les tests et le pipeline complet ne sont
  pas validés.

## Format des prochaines mises à jour

Lors d'une prochaine session, mettre à jour en priorité :

- `Dernière mise à jour` ;
- `État actuel` si la branche, la PR ou l'environnement changent ;
- `Travail terminé` avec seulement les nouveaux résultats importants ;
- `Problèmes connus` en supprimant ceux qui sont résolus ;
- `Prochaines étapes recommandées` avec l'action suivante la plus concrète.

