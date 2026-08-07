# Étape 19 — Rapport révisable et validation humaine

L'étape 19 assemble les résultats validés sans calcul scientifique. Elle produit
d'abord un rapport HTML de travail. Les commentaires proposés restent séparés
des faits contrôlés et doivent être corrigés puis approuvés par une personne.

## Sources autorisées

- artefacts signés de parenté KING produits par l'étape 07 ;
- index, complétude, figures et manifest de rendu produits par l'étape 18 ;
- manifest, configuration résolue et audits du run courant, après contrôle de
  leurs empreintes.

Les exports ACPA, les sorties historiques de `data/output/`, les génotypes
individuels et les identifiants sources sont interdits. Les paires KING sont
pseudonymisées spécifiquement pour le rapport et restent `sensitive_genetic`.

## Assistance rédactionnelle

`interpretation_facts.json` est la seule charge utile autorisée pour une IA. Le
prompt impose une rédaction descriptive, la séparation observation / prudence /
limites et l'absence de conclusion pour `NOT_EVALUATED`. Une réponse doit être
un JSON conforme au schéma et passe des contrôles déterministes avant affichage.
Le fournisseur externe est désactivé par défaut ; une proposition déterministe
prudente permet de tester tout le parcours sans réseau.

## Relecture et finalisation

Le HTML de brouillon contient une zone modifiable par section, les faits sources
non modifiables et une checklist finale. Le bouton de validation télécharge un
`report_review.json`; il ne modifie pas le run depuis le navigateur. Le
finaliseur vérifie la décision, les empreintes, les mots interdits et la
complétude, puis produit le HTML verrouillé avant le PDF. La provenance conserve
le relecteur, la date, la révision et la méthode de génération, sans secret.

`READY_FOR_SCIENTIFIC_REVIEW` signifie que le rapport est techniquement complet
et approuvé pour relecture scientifique. Il ne signifie pas que l'hypothèse
d'effet fondateur est démontrée.
