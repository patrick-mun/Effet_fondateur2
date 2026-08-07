# Analyses de sensibilité — étape 17

L'étape 17 consolide des runs V2 immuables. Elle ne relance pas les étapes
scientifiques et ne modifie jamais leur résultat primaire.

## Unité d'analyse

Un scénario correspond à un run distinct avec une configuration signée. Le
registre désigne exactement un run primaire et zéro ou plusieurs runs de
sensibilité. Il enregistre le facteur modifié, les domaines attendus et les
empreintes nécessaires pour vérifier les manifestes et résumés.

Le modèle `config/sensitivity/scenarios.example.tsv` contient volontairement
des chemins et empreintes factices. Il doit être copié puis renseigné avant
d'activer `run_sensitivity_analyses`. Une empreinte de manifeste obsolète bloque
la consolidation au lieu de réutiliser silencieusement un run modifié.

Un scénario à un facteur permet d'attribuer une variation à cet axe. Un scénario
multifactoriel reste descriptif et doit être étiqueté comme tel. Les retraits se
font par unité indépendante ou par porteur explicitement identifié depuis les
génotypes moléculaires acceptés, jamais depuis le phénotype.

## Domaines comparés

- `FOUNDER_IBS` : statut technique et limites du candidat IBS de l'étape 13 ;
- `VARIANT_AGE` : statut, estimation et intervalle de datation de l'étape 14 ;
- `LOCAL_LD` : statuts descriptifs par cohorte de l'étape 15 ;
- `ROH` : statuts par périmètre et indicateurs secondaires de l'étape 16.

Les domaines restent séparés. Leur concordance peut renforcer la cohérence du
faisceau d'arguments, mais l'étape 17 ne calcule ni score global ni probabilité
d'effet fondateur.

## Statuts

La stabilité catégorielle compare le statut du scénario à celui du primaire :

- `STABLE` si tous les scénarios attendus et évaluables conservent le statut ;
- `VARIABLE` si au moins un scénario évaluable change le statut ;
- `INCONCLUSIVE` si le primaire n'est pas concluant ou si les scénarios requis
  sont systématiquement non évaluables ;
- `NOT_EVALUATED` si aucun résultat comparable n'est disponible.

Les métriques numériques sont publiées avec leur minimum, maximum et changement
relatif lorsque cela a un sens. Sans tolérance préspécifiée, leur stabilité est
`NOT_CLASSIFIED`. Une tolérance ajoutée après lecture des résultats invaliderait
la nature préspécifiée de l'analyse et exige un nouveau run.

## Limites

La robustesse observée ne traite que les scénarios enregistrés. Elle ne corrige
pas un biais commun à tous les runs, une erreur moléculaire, une faible densité
de marqueurs ou une mauvaise carte. Elle ne constitue pas une validation dans
une cohorte externe et ne transforme pas un partage IBS en IBD.
