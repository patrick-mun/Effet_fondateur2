# Contribuer à Effet_fondateur2

Ce projet est un pipeline de recherche. Une modification techniquement correcte
peut changer une population analysée, une unité statistique ou l'interprétation
d'un résultat. La revue doit donc distinguer explicitement le code
d'orchestration, les contrats de données et les méthodes scientifiques.

## Périmètre actif

Le pipeline maintenu est le paquet `src/effet_fondateur/`. Le contenu de
`archive/` est une preuve historique : il ne doit pas être importé, exécuté ou
recopié dans la V2. Les fichiers de `data/output/` et `data/runs/` sont des
résultats générés, jamais des sources de code ou de vérité scientifique.

La cible génétique appartient au profil d'étude. Aucun gène, chromosome,
variant, allèle ou position propre à une étude ne doit être codé en dur dans le
moteur.

## Lire le projet

- `src/effet_fondateur/stages/` contient les points d'entrée des étapes. Une
  fonction `execute()` orchestre les entrées, le calcul, la validation et la
  publication ; elle ne devrait pas cacher une nouvelle méthode scientifique.
- Les paquets `ancestry/`, `dating/`, `explicit_ibd/`, `founder/`, `ld/`,
  `phasing/`, `references/`, `roh/` et `sensitivity/` portent la logique métier
  testable indépendamment de l'orchestrateur.
- `contracts/` et `resources/schemas/` définissent les frontières de données.
  Un changement de colonne, de vocabulaire ou de sensibilité est un changement
  de contrat, pas une simple modification de sérialisation.
- `external_tools.py` est la frontière commune vers les exécutables externes.
  Les modules consommateurs traduisent ses erreurs techniques en erreurs métier.
- `docs/modules/` décrit les méthodes, hypothèses, entrées, sorties et limites de
  chaque étape. Cette documentation fait partie du contrat de maintenance.

## Environnement et contrôles rapides

Utiliser Python 3.12 dans l'environnement du dépôt :

```bash
.venv/bin/python -m pip install -r requirements.txt
.venv/bin/python -m pytest tests/
.venv/bin/ruff check src/effet_fondateur tests
.venv/bin/mypy
```

PLINK, KING, bcftools, SHAPEIT5, Java, Hap-IBD et Refined IBD doivent être
simulés dans les tests unitaires. Leur exécution réelle appartient à des tests
d'intégration explicitement identifiés et ne doit jamais écraser un run ou une
donnée d'entrée existante.

## Documenter le code

Une docstring publique ou complexe doit préciser, selon le cas :

- la responsabilité de la fonction et ce qu'elle ne conclut pas ;
- la forme des entrées et sorties, notamment l'ordre des échantillons ou des
  variants ;
- les unités (`bp`, `cM`, Morgans, générations, années, fréquences) ;
- l'unité statistique indépendante (individu, famille ou haplotype) ;
- les hypothèses, seuils et règles d'exclusion ;
- les erreurs qui bloquent le calcul ;
- le document de `docs/modules/` qui porte la méthode complète.

Les commentaires internes expliquent une décision non évidente : prévention
d'un biais, contrainte d'un outil, sécurité, provenance ou raison scientifique.
Ils ne doivent pas paraphraser une affectation. Un long orchestrateur peut
utiliser quelques commentaires de phase pour séparer validation, calcul,
contrôles et publication.

Toute formule nouvelle ou modifiée doit être accompagnée de sa référence, de
ses unités, de ses hypothèses et d'un test sur une valeur de référence. Une
ligne non couverte et un commentaire absent ne constituent jamais, seuls, une
preuve de code mort.

## Préserver le sens scientifique

- Ne jamais déduire le génotype cible depuis le statut clinique ou le groupe.
- Ne pas confondre IBS, IBD, LD, ROH, ascendance et effet fondateur.
- Ne pas compter plusieurs personnes d'une même famille comme unités
  indépendantes lorsque le contrat retient la famille.
- Ajuster les PCA de référence uniquement sur les références prévues ; les
  échantillons d'étude sont projetés dans des axes figés.
- Conserver les graines, paramètres, versions d'outils, empreintes et motifs
  d'exclusion nécessaires à la reproduction.
- Une approximation ou une limite doit rester visible dans le code, le document
  de méthode et l'audit publié.

## Préparer une modification

1. Identifier le contrat et le document de méthode concernés.
2. Ajouter ou sélectionner un test qui observe le comportement actuel.
3. Modifier un seul module ou une seule responsabilité à la fois.
4. Relancer les tests ciblés, puis les contrôles statiques.
5. Mettre à jour la documentation dans le même changement si une commande, une
   hypothèse, une unité, une sortie ou une interprétation évolue.
6. Vérifier que le diff ne contient ni donnée génétique, ni résultat de run, ni
   secret, ni chemin absolu propre à une machine.

Un refactoring doit conserver les fichiers publiés, codes d'erreur, audits et
résultats numériques. Une modification scientifique demande en plus une
comparaison avec une référence et doit être annoncée explicitement dans la
revue.
