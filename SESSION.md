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
- Création du convertisseur multi-échantillons `acpa_to_plink.py` avec modèle
  de métadonnées, contrôles QC et génération PED/MAP/groupes/cas/témoins.
- Ajout de l'annotation dbSNP hg38 après sélection du chromosome, avec
  propagation depuis un export ACPA annoté et résolution des SNV par allèles.
- Documentation de la conversion ACPA et ajout de tests unitaires dédiés.
- Ajout dans le README d'un protocole ACPA numéroté, incluant la préparation
  détaillée de `samples.tsv`, l'annotation, le contrôle des rapports et PLINK.
- Création d'un simulateur dédié de témoins indépendants sous HWE, reproductible
  et non destructif, avec jeu combiné, listes cas/témoins et rapport de limites.
- Création d'un injecteur de mutation non destructif : insertion PED/MAP triée,
  affectations explicites par groupe ou individu, métadonnées de rapport,
  empreintes SHA-256 et audit complet des génotypes.
- Suppression des artefacts versionnés et lanceurs historiques devenus
  obsolètes : binaire KING dupliqué, ancien rapport de tests, ancien script
  shell et plan structurel désynchronisé.

## Validations effectuées

- `pip check` : aucune dépendance cassée.
- Import de `run_pipeline` : réussi.
- Imports Python principaux : réussis.
- Packages R `adegenet`, `ggplot2`, `ape`, `ade4` et `poppr` : disponibles.
- Test PLINK : 2 726 variants et 127 individus lus correctement.
- Conversion temporaire des 14 exports ACPA réels sur le chromosome 19 :
  2 725 marqueurs, 2 725 rsID résolus, aucun rsID non résolu et validation PLINK
  réussie avec un taux de génotypage de 98,04 %.
- Tests du convertisseur ACPA : 8 réussis.
- Simulation temporaire de 50 témoins : 59 individus et 2 725 marqueurs lus par
  PLINK, puis analyse KING réussie sans modification des données source.
- Tests du simulateur de témoins : 3 réussis.
- Tests de l'injecteur de mutation : 5 réussis.
- Tests Python : 28 réussis, 4 anciens tests défaillants.
- Injection technique temporaire sur le jeu combiné : 59 individus, passage de
  2 725 à 2 726 marqueurs et lecture PLINK réussie à 98,25 % de génotypage.
- Contrôle mendélien temporaire : 10 erreurs au total, dont une au marqueur
  injecté pour `F2/E_82303707` avec les règles de groupe d'exemple.

## Problèmes connus

1. `data/input/complex_simulation/cas.txt` et `temoins.txt` restent absents tant
   que les métadonnées réelles ne sont pas renseignées et la conversion validée.
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
8. Les allèles génomiques REF/ALT de la mutation DOCK6 doivent être confirmés
   sur GRCh38 avant une injection réelle ; la notation HGVS du transcrit ne
   suffit pas à déterminer le brin écrit dans le PED.
9. Les règles d'exemple produisent une incompatibilité mendélienne au marqueur
   injecté pour `F2/E_82303707` (`*/* x T/T -> G/G`) ; utiliser le résultat
   moléculaire individuel plutôt que le seul groupe clinique.

## Prochaines étapes recommandées

1. Confirmer la position GRCh38, le transcrit et les allèles génomiques REF/ALT,
   puis compléter `mutation_info.json` et les génotypes individuels réels.
2. Utiliser les témoins synthétiques uniquement pour des validations techniques
   ou exploratoires, jamais pour conclure biologiquement sur le LD ou les ROH.
3. Actualiser les quatre tests défaillants et obtenir une suite entièrement
   verte.
4. Raccorder les uploads Streamlit aux fichiers réellement consommés.
5. Revoir les entrées mathématiques de l'estimation Gamma.
6. Exécuter le pipeline complet sur une copie contrôlée des données.
7. Examiner le rapport produit avant de fusionner la PR.

## Décisions à conserver

- Utiliser le nom standard `AGENTS.md` pour les instructions Codex.
- Conserver les deux convertisseurs : ACPA vers VCF, puis VCF vers MAP/allèles.
- Conserver les anciens convertisseurs et le simulateur familial jusqu'à
  extraction de leur éventuelle logique de mutation encore utile.
- Garder la PR en brouillon tant que les tests et le pipeline complet ne sont
  pas validés.

## Format des prochaines mises à jour

Lors d'une prochaine session, mettre à jour en priorité :

- `Dernière mise à jour` ;
- `État actuel` si la branche, la PR ou l'environnement changent ;
- `Travail terminé` avec seulement les nouveaux résultats importants ;
- `Problèmes connus` en supprimant ceux qui sont résolus ;
- `Prochaines étapes recommandées` avec l'action suivante la plus concrète.
