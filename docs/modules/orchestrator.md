# Orchestrateur V2 — maintenance et extension

## Responsabilité

L'orchestrateur contrôle l'ordre et la traçabilité des étapes. Il ne contient
aucune formule scientifique et ne transforme aucun génotype. Une méthode
scientifique appartient toujours à un module d'étape séparé et testable.

Les invariants suivants doivent rester vrais :

- la configuration résolue est immuable pendant un run ;
- une étape est un module explicitement autorisé par le catalogue ;
- la commande utilise une liste d'arguments et jamais `shell=True` ;
- une tentative écrit uniquement dans son dossier temporaire ;
- les sorties ne sont publiées qu'après validation des schémas et empreintes ;
- une tentative échouée reste disponible pour le diagnostic ;
- une étape publiée n'est réutilisée qu'après contrôle d'intégrité ;
- le manifest ne contient aucune donnée génétique individuelle.

## Répartition des modules

| Module | Responsabilité |
|---|---|
| `bootstrap.py` | créer atomiquement le run et l'étape `00` |
| `catalog.py` | autoriser et valider les définitions d'étapes |
| `environment.py` | inventorier Python, la plateforme et les outils |
| `pipeline.py` | parcourir les étapes activées dans la configuration |
| `runner.py` | gérer le cycle d'une tentative et le sous-processus |
| `signatures.py` | calculer la signature déterministe d'une étape |
| `integrity.py` | contrôler chemins, documents et empreintes |
| `inputs.py` | déclarer les sources externes et artefacts dépendants |
| `state.py` | lire et écrire le manifest et le journal |
| `models.py` | définir les objets stables partagés |

Cette séparation évite que l'ajout d'un format scientifique modifie la logique
de reprise ou que l'ajout d'une règle d'intégrité modifie le lancement des
sous-processus.

## Cycle d'une étape

1. Le pipeline vérifie l'empreinte de `config.resolved.yaml`.
2. Le catalogue refuse les noms, identifiants ou dépendances incohérents.
3. Le runner vérifie l'état des dépendances et résout leurs artefacts requis.
4. La signature lie code, paramètres, entrées, configuration et contrats.
5. Le runner crée un dossier temporaire et `stage_inputs.json`.
6. Le script d'étape produit ses artefacts, `stage_outputs.json` et `audit.json`.
7. `integrity.py` vérifie l'identité, les chemins et toutes les empreintes.
8. Le dossier temporaire est renommé atomiquement vers son nom définitif.
9. Le manifest puis le journal sont actualisés.

Une sortie présente mais non déclarée n'est pas un artefact. Une sortie déclarée
mais absente, déplacée ou modifiée provoque un échec.

## Ajouter une étape scientifique

L'ajout d'une étape ne doit pas nécessiter de modifier `runner.py`.

1. Créer `src/effet_fondateur/stages/<nom>.py`.
2. Fournir une fonction testable qui reçoit `stage_inputs.json` et
   `--output-dir`.
3. Valider `stage_inputs.json` avant toute lecture scientifique.
4. Déclarer chaque artefact avec `build_file_artifact()`.
5. Produire `stage_outputs.json`, `audit.json` et `checksums.sha256`.
6. Retourner un code V2 documenté ; ne jamais ignorer un outil externe en échec.
7. Ajouter une `StageDefinition` au catalogue de production avec un ID unique,
   ses dépendances, entrées configurées, artefacts requis et sa criticité.
8. Ajouter la section de configuration et ses paramètres versionnés.
9. Ajouter des tests unitaires avec outils externes simulés.
10. Ajouter au moins un test d'intégration synthétique du contrat complet.

Le nom du module n'est jamais lu directement depuis le YAML. Cette liste blanche
est volontaire : une configuration ne doit pas pouvoir importer arbitrairement
du code Python.

Une étape peut déclarer `blocking_manual_decision_ids`. Le runner refuse alors
son lancement tant que le manifest contient l'une de ces décisions. Les étapes
`03_convert_acpa` et `04_prepare_target_variant_dataset` utilisent ce mécanisme
pour exiger l'approbation du registre maître avant toute transformation
génétique.

## Modifier un contrat

Une modification incompatible exige :

- une nouvelle version de schéma ;
- une mise à jour de la version incluse dans la signature ;
- un test de rejet de l'ancien format ou une migration explicite ;
- une mise à jour du producteur et de tous les consommateurs ;
- une note dans `SESSION.md` si le comportement du pipeline change.

Ne jamais assouplir un schéma uniquement pour accepter un fichier défectueux.

## Reprise et diagnostic

Les dossiers `stages/attempts/*.failed` conservent les entrées déclarées, la
commande structurée et les journaux de la tentative. Ils ne sont jamais lus
comme résultats scientifiques.

Une reprise conserve la même configuration. Si un paramètre scientifique doit
changer, il faut créer un nouveau run : modifier `config.resolved.yaml` bloque
volontairement le run existant.

## Limites actuelles connues

Le socle minimal ne fournit pas encore :

- verrou interprocessus empêchant deux reprises simultanées ;
- délai maximal ou protocole d'annulation des outils externes longs ;
- cache partagé entre deux runs distincts ;
- registre de production des étapes scientifiques `06` à `19` ;
- artefacts composés BED/BIM/FAM ;
- migration de version automatique des manifests.

Ces limites ne doivent pas être contournées silencieusement. Elles seront
traitées avant l'utilisation en production des étapes concernées.
