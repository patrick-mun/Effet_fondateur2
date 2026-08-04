# Instructions Codex — Effet_fondateur2

## Objectif du projet

Ce dépôt contient un pipeline de recherche pour étudier un possible effet
fondateur autour d'une mutation de `DOCK6` sur le chromosome 19. Le moteur
combine PLINK, KING, Python, R/Adegenet et une estimation Gamma.

## Environnement de travail

- Travailler depuis la racine du dépôt.
- Utiliser Python 3.12 dans `.venv`.
- Activer l'environnement avec `source .venv/bin/activate` ou appeler
  explicitement `.venv/bin/python`.
- Installer les dépendances avec `.venv/bin/python -m pip install -r requirements.txt`.
- Les exécutables `plink`, `king`, `bcftools` et `Rscript` doivent être dans le
  `PATH`.
- `Gamma` est optionnel pour le pipeline principal actuel.

## Points d'entrée

- Pipeline : `.venv/bin/python run_pipeline.py`
- Interface : `.venv/bin/streamlit run interface_effet_fondateur.py`
- Tests : `.venv/bin/python -m pytest test/`
- Conversion multi-échantillons ACPA vers PLINK :
  `.venv/bin/python -m simulation_genotype_famille.acpa_to_plink`
- Conversion historique ACPA vers VCF sans génotypes individuels :
  `.venv/bin/python simulation_genotype_famille/acpa_to_vcf2.py`
- Simulation exploratoire de témoins indépendants :
  `.venv/bin/python -m simulation_genotype_famille.simulate_unrelated_controls`
- Simulation familiale :
  `.venv/bin/python simulation_genotype_famille/ped_generator_precod.py`

## Données et sécurité

- Ne pas supprimer ou remplacer les données de `data/input/` sans demande
  explicite de l'utilisateur.
- Ne pas considérer les fichiers de `data/output/` comme des sources : ce sont
  des résultats générés, parfois historiques.
- Ne pas lancer le pipeline complet sans prévenir l'utilisateur : il réécrit
  des sorties et ouvre le rapport HTML dans le navigateur.
- Ne jamais inclure de données génétiques, jetons ou secrets dans les commits,
  les rapports de session ou les messages GitHub.
- Préserver les modifications locales existantes et inspecter `git status`
  avant toute édition ou opération Git.

## Règles de modification

- Corriger la cause d'un problème avec des changements ciblés.
- Conserver les chemins actuels tant qu'une refactorisation n'est pas demandée.
- Utiliser `apply_patch` pour modifier les fichiers texte.
- Ne pas committer les caches, `.venv`, données d'entrée, résultats ou rapports
  générés.
- Mettre à jour `README.md`, `requirements.txt` et `SESSION.md` lorsque les
  commandes, dépendances ou priorités changent.

## Règles de codage et maintenabilité

### Nommage et structure

- Utiliser des noms explicites pour les variables, fonctions, classes et
  fichiers ; éviter les noms d'une lettre hors indices locaux évidents.
- Conserver une langue et une terminologie cohérentes dans un même module.
- Limiter chaque fonction à une responsabilité clairement identifiable.
- Préférer plusieurs petites fonctions testables à une fonction monolithique.
- Éviter la duplication et extraire les traitements réellement partagés.
- Remplacer progressivement les chemins et paramètres codés en dur par une
  configuration explicite, sans modifier le comportement sans validation.
- Ajouter des annotations de types aux nouvelles fonctions et lors des
  refactorisations ciblées, en priorité pour les chemins, tableaux et résultats.

### Commentaires et documentation

- Commenter le pourquoi d'une décision scientifique ou technique, pas une
  simple traduction de ce que fait le code.
- Supprimer les commentaires devenus faux en même temps que le code concerné.
- Ajouter une docstring aux fonctions publiques ou complexes en précisant les
  entrées, sorties, unités, hypothèses et erreurs possibles.
- Documenter toute approximation biologique, statistique ou génétique à
  proximité de son implémentation et dans la documentation adaptée.

### Gestion des erreurs et débogage

- Ne pas ignorer silencieusement une exception ou un code de retour externe.
- Produire des erreurs indiquant l'étape, la commande et le fichier concernés,
  sans exposer de donnée sensible.
- Utiliser `logging` pour le diagnostic durable et réserver `print` aux sorties
  utilisateur intentionnelles des scripts interactifs.
- Retirer les traces temporaires, fichiers de diagnostic et points d'arrêt avant
  un commit.
- Reproduire un bug avec le plus petit jeu de données possible avant de le
  corriger, puis vérifier la cause racine.

### Tests unitaires

- Ajouter ou actualiser un test ciblé pour chaque correction de bug et chaque
  nouveau comportement testable.
- Tester le comportement observable, les cas limites et les erreurs attendues,
  sans dépendre des détails internes inutiles.
- Simuler PLINK, KING, Gamma ou R dans les tests unitaires ; réserver leur
  exécution réelle à des tests d'intégration identifiés.
- Utiliser des dossiers temporaires et ne jamais écrire dans les données ou
  résultats réels pendant les tests.
- Ne pas affaiblir une assertion uniquement pour faire passer un test.

### Refactoring

- Séparer autant que possible le refactoring d'un changement fonctionnel.
- Préserver le comportement existant avec des tests avant de déplacer ou
  simplifier du code.
- Refactoriser par petites étapes vérifiables et relancer les tests ciblés après
  chaque étape significative.
- Ne pas élargir une refactorisation à des modules sans rapport avec la tâche.

### Calcul scientifique et reproductibilité

- Nommer ou documenter les unités utilisées (`bp`, `kb`, `cM`, générations,
  fréquences) et vérifier les conversions aux frontières des fonctions.
- Documenter les formules, hypothèses, seuils et références des méthodes
  statistiques avant de les modifier.
- Fixer et enregistrer les graines aléatoires des simulations reproductibles.
- Enregistrer les paramètres et versions d'outils nécessaires pour reproduire
  une analyse.
- Ne pas modifier une méthode scientifique ou son interprétation sans test,
  comparaison avec un résultat de référence et signalement explicite.

## Validation

- Commencer par `git diff --check`.
- Vérifier les imports avec `.venv/bin/python -c "import run_pipeline"`.
- Lancer les tests ciblés avant la suite complète.
- Pour une modification liée à PLINK, utiliser une sortie temporaire et ne pas
  écraser les résultats existants.
- Signaler clairement les tests déjà défaillants sans corriger des problèmes
  sans rapport avec la tâche.

## Continuité des sessions

- Lire `SESSION.md` au début d'une nouvelle session.
- Mettre à jour uniquement les sections concernées avant de terminer une tâche
  importante ou lorsque le contexte risque d'être compacté.
- Garder `SESSION.md` synthétique : état, décisions, validations, blocages et
  prochaine action concrète.
- Remplacer les informations devenues fausses au lieu d'accumuler un journal
  exhaustif de toutes les conversations.
