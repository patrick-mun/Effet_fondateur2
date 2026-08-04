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
- Les exécutables `plink`, `king` et `Rscript` doivent être dans le `PATH`.
- `Gamma` est optionnel pour le pipeline principal actuel.

## Points d'entrée

- Pipeline : `.venv/bin/python run_pipeline.py`
- Interface : `.venv/bin/streamlit run interface_effet_fondateur.py`
- Tests : `.venv/bin/python -m pytest test/`
- Conversion ACPA vers VCF :
  `.venv/bin/python simulation_genotype_famille/acpa_to_vcf2.py`
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

