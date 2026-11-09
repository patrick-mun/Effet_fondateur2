# Conventions de développement de la calibration

**Version :** 1.1
**Date :** 30 octobre 2026
**Périmètre :** `calibration/`, `src/effet_fondateur/calibration/` et
`tests/calibration/`

## 1. Autorité et objectif

`AGENTS.md`, à la racine du dépôt, reste la source d'autorité. Le présent
document reformule ses exigences pour le sous-système de calibration et ajoute
les précautions propres à la simulation génétique. En cas de divergence,
`AGENTS.md` prévaut et ce document doit être corrigé dans le même changement.

Ces conventions s'appliquent aux documents exécutables, schémas, fichiers de
configuration, modules Python, tests, commandes externes, évaluateurs et
générateurs de rapports de calibration.

## 2. Environnement et emplacements

- Travailler depuis la racine du dépôt.
- Utiliser Python 3.12 via `.venv/bin/python` et les outils installés dans
  `.venv`.
- Versionner les protocoles, configurations, schémas, modèles, sources et
  tests dans leurs dossiers prévus.
- Écrire les panels, vérités, sorties, figures, caches et rapports générés
  uniquement sous `calibration_workspace/`.
- Ne jamais écrire un test dans `data/input/`, `data/output/`, `data/runs/` ou
  dans un résultat scientifique existant.
- Ne pas coder de chemin absolu local. Les chemins autorisés doivent venir
  d'une configuration validée ou d'un argument explicite.
- Ne pas lancer le pipeline complet ou un outil externe réel sans prévenir
  l'utilisateur et identifier l'exécution comme test d'intégration.

## 3. Données, confidentialité et aveuglement

- Ne jamais lire, copier ou dériver une donnée depuis `data/input/dock6_real/`
  pour construire ou régler un panel de calibration.
- N'introduire dans Git aucun génotype, segment individuel, identifiant
  personnel, secret, jeton, panel généré ou vérité simulée volumineuse.
- Employer uniquement des identifiants synthétiques sans correspondance avec
  une personne réelle.
- Séparer physiquement `private_truth/` de `blinded_inputs/` et interdire toute
  référence à la vérité dans les entrées remises au pipeline.
- Ne jamais déduire un génotype à partir d'un statut clinique, d'un groupe ou
  d'une étiquette simulée cachée.
- Interdire le partage d'individus, de graines et de transmissions entre panel
  de calibration et panel de validation.
- Régler les seuils uniquement sur le panel de calibration. Une fois scellés,
  ils ne peuvent pas être modifiés après ouverture du panel aveugle.
- Publier uniquement des résultats agrégés dans les rapports versionnables ou
  partageables.

## 4. Nommage, structure et types

- Utiliser des noms explicites pour les modules, classes, fonctions,
  paramètres, variables et fichiers.
- Éviter les noms d'une lettre, sauf pour un indice local évident et limité.
- Conserver une terminologie cohérente dans chaque module. Les identifiants de
  code peuvent rester en anglais ; les explications et rapports sont rédigés
  en français.
- Donner une responsabilité identifiable à chaque fonction et préférer de
  petites fonctions testables à une fonction monolithique.
- Extraire un traitement partagé seulement lorsqu'il est réellement commun ;
  ne pas créer d'abstraction anticipée sans usage démontré.
- Ajouter des annotations de types aux nouvelles fonctions, notamment pour les
  `Path`, configurations, tableaux, graines, unités et résultats.
- Utiliser des structures immuables ou explicitement copiées lorsque plusieurs
  scénarios partagent une configuration de base.
- Garder les paramètres scientifiques et chemins hors du code : ils doivent
  être déclarés dans une configuration validée.
- Refuser les clés inconnues dans les configurations et schémas scientifiques.

## 5. Documentation du code scientifique

- Documenter le pourquoi d'une décision, d'un seuil ou d'une approximation,
  sans paraphraser chaque instruction Python.
- Ajouter une docstring aux fonctions publiques et aux fonctions scientifiques
  complexes.
- Préciser dans les docstrings les entrées, sorties, unités, hypothèses,
  valeurs manquantes et erreurs possibles.
- Nommer les unités dans les variables ou les documenter sans ambiguïté :
  `bp`, `kb`, `cM`, Morgans, générations, proportions et fréquences.
- Vérifier les conversions d'unités aux frontières des fonctions.
- Documenter les formules, références, hypothèses et domaines de validité avant
  de modifier une méthode scientifique.
- Supprimer ou corriger un commentaire devenu faux en même temps que le code
  concerné.
- Toute modification de méthode reçoit un identifiant de méthode nouveau ou
  une justification de compatibilité explicitement testée.

## 6. Reproductibilité

- Exiger et enregistrer toutes les graines aléatoires ; aucun générateur ne
  doit employer silencieusement une graine issue de l'heure ou du système.
- Enregistrer les paramètres, versions d'outils, commandes, références,
  assemblages, cartes et empreintes SHA-256 nécessaires à la reproduction.
- Produire un identifiant immuable par run et refuser l'écrasement d'un run
  existant.
- Garantir qu'une même configuration, les mêmes entrées et les mêmes graines
  produisent les mêmes artefacts déterministes.
- Séparer les partitions de donneurs, références de phasage, calibration et
  validation ; vérifier explicitement que les intersections interdites sont
  vides.
- Une correction de biais ne doit jamais remplacer silencieusement le résultat
  brut : publier les deux avec des identifiants de méthode distincts.

## 7. Erreurs, journaux et commandes externes

- Ne jamais ignorer silencieusement une exception ou un code de retour.
- Les erreurs doivent nommer l'étape, l'opération et le fichier concernés sans
  exposer d'identifiant ni de donnée sensible.
- Utiliser `logging` pour le diagnostic durable. Réserver `print` aux sorties
  utilisateur intentionnelles d'une commande interactive.
- Borner les commandes externes en temps, mémoire et nombre de threads lorsque
  ces paramètres sont disponibles.
- Enregistrer la commande normalisée, la version de l'outil, le code de retour
  et les chemins d'artefacts, sans consigner de secret.
- PLINK, KING, bcftools, SHAPEIT5, Hap-IBD et Refined IBD sont simulés dans les
  tests unitaires. Leur exécution réelle est réservée aux tests d'intégration
  annoncés et écrit uniquement dans un dossier temporaire ou un run neuf.
- Retirer traces de débogage, fichiers temporaires et points d'arrêt avant tout
  commit.

## 8. Tests

- Ajouter ou mettre à jour un test ciblé pour chaque correction, contrat ou
  comportement nouveau.
- Tester le comportement observable, les cas limites et les erreurs attendues,
  sans figer inutilement l'implémentation interne.
- Utiliser `tmp_path` ou une autre zone temporaire lorsqu'un test écrit des
  fichiers.
- Les fixtures sont petites, synthétiques et dépourvues de correspondance avec
  une donnée réelle.
- Tester explicitement les valeurs manquantes, graines absentes, effectifs
  incohérents, unités invalides, chemins interdits, collisions et tentatives de
  fuite de vérité.
- Vérifier les scénarios positifs et négatifs. Une absence de détection ne doit
  jamais être transformée en preuve d'absence.
- Ne pas affaiblir une assertion ni modifier un seuil scientifique uniquement
  pour faire réussir un test.
- Reproduire un bug avec le plus petit cas synthétique possible avant de
  corriger sa cause racine.
- Identifier clairement les tests d'intégration réels afin qu'ils ne soient pas
  déclenchés comme des tests unitaires ordinaires.

## 9. Refactorisation et portée des changements

- Préserver les modifications locales existantes et vérifier `git status`
  avant toute édition ou opération Git.
- Corriger la cause avec un changement ciblé et conserver les chemins existants
  tant qu'une refactorisation n'est pas demandée.
- Séparer autant que possible une refactorisation d'un changement fonctionnel
  ou scientifique.
- Établir des tests de comportement avant de déplacer ou simplifier du code.
- Refactoriser par petites étapes vérifiables et relancer les tests ciblés
  après chaque étape significative.
- Ne pas étendre une refactorisation à des modules sans rapport avec la tâche.
- Utiliser `apply_patch` pour les modifications textuelles ciblées.

## 10. Rapports et interprétation

Exigence permanente de l'utilisateur, 30 octobre 2026 : chaque rapport doit
faciliter l'explication scientifique par des graphiques de résultats et une
représentation mathématique adaptée, pas seulement documenter l'audit.

- Prévoir deux à quatre visuels utiles, avec au moins deux graphiques
  quantitatifs quand les données s'y prêtent : distributions, projections,
  matrices de confusion, effectifs, comparaisons ou intervalles.
- Les cartes de chiffres et schémas de parcours ne remplacent pas les
  graphiques scientifiques. Une étape non quantitative expose sa logique et
  justifie explicitement l'absence de graphique quantitatif.
- Présenter les équations réellement implémentées avec symboles définis,
  unités, hypothèses, seuils, domaine de validité et interprétation. Pour une
  sélection de cohortes, privilégier ensembles et règles de décision.
- Expliquer chaque graphique et chaque formule dans le texte. Légendes,
  axes, unités, effectifs, type d'incertitude et sources sont obligatoires.
- Séparer dispersion et intervalle de confiance, données indépendantes et
  copies/profils appariés. Aucun identifiant individuel dans un graphique.
- Les agrégations destinées aux figures sont préparées séparément du rendu,
  empreintées et rattachées aux sorties closes. Le rendu ne refait ni PCA,
  estimation, test ou qualification.
- Relire les figures et formules dans le PDF réel, en plus du contrôle HTML.

- Les évaluateurs calculent les métriques ; le moteur de rendu ne recalcule
  aucune statistique scientifique.
- Distinguer systématiquement vérité simulée, observation, erreur,
  incertitude, limite et décision.
- Conserver les statuts `PASS`, `FAIL`, `NOT_EVALUABLE` et
  `REVIEW_REQUIRED` sans les requalifier pendant le rendu.
- Ne jamais présenter IBS, LD, ROH ou proximité PCA comme une preuve IBD.
- Conserver `founder_effect_proven: false` tant qu'une décision scientifique
  humaine ne s'appuie pas sur les rapports qualifiés prévus.
- Une absence d'appel IBD ne réfute pas automatiquement un effet fondateur.
- Toute conclusion issue du panel porte la mention `SIMULATION_CALIBRATED` et
  ne constitue pas une validation biologique externe.

## 11. Validation avant commit

La validation est proportionnée au changement et suit cet ordre :

1. `git diff --check` ;
2. vérification des imports :
   `.venv/bin/python -c "import effet_fondateur; import effet_fondateur.orchestrator.pipeline"` ;
3. validation syntaxique des JSON, YAML et schémas modifiés ;
4. tests ciblés du comportement ajouté ou corrigé ;
5. tests des modules dépendants ;
6. suite complète lorsque le changement peut affecter le pipeline au-delà du
   sous-système de calibration.

Pour une modification liée à PLINK ou à un autre outil externe, les sorties de
validation sont écrites dans un dossier temporaire et ne remplacent jamais un
résultat existant. Toute défaillance antérieure sans rapport avec la tâche est
signalée séparément et n'est pas corrigée sans autorisation.

## 12. Définition de terminé

Un changement de calibration est terminé uniquement si :

1. son périmètre et sa responsabilité sont clairs ;
2. les types, unités, hypothèses et erreurs sont documentés ;
3. les tests positifs, négatifs et de sécurité pertinents réussissent ;
4. aucune donnée réelle, vérité générée ou sortie de run n'est suivie par Git ;
5. les graines, versions, paramètres et empreintes nécessaires sont publiés ;
6. les résultats manquants ne peuvent pas devenir un succès ;
7. la séparation calibration-validation et vérité-analyse est préservée ;
8. les documents concernés et `SESSION.md` sont cohérents avec le changement ;
9. la validation exécutée et les validations non exécutées sont signalées ;
10. aucun seuil ou résultat scientifique n'a été modifié pour favoriser le cas
    DOCK6.
