# Étape 12.6 — Phasage et chromosome porteur

## Exécution

La sous-étape `12.6` consomme exclusivement le manifeste validé de `12.5`.
Elle sonde les deux exécutables et exige SHAPEIT5 `5.1.1`, exécute
`phase_common` avec le panel de référence, puis `phase_rare` avec le scaffold
commun. Les BCF et leurs index CSI sont publiés avec les journaux des deux
passes. Le journal de `phase_rare` est la capture stdout/stderr, car cette
version ne prend pas d'option `--log`.

La graine, les threads, la taille effective `Ne`, le délai et le seuil de
confiance sont explicites dans la configuration et le manifeste. Un pedigree
vide omet entièrement l'argument `--pedigree`.

## Contrôles scientifiques

Avant le phasage, les empreintes des entrées, l'ordre des individus et les
compatibilités mendéliennes régionales sont vérifiés. Après les deux passes, le
pipeline exige :

- le même ordre d'individus dans l'étude, le scaffold et le BCF final ;
- un scaffold commun phasé, sans variant étranger à l'étude ;
- la conservation exacte de tous les variants et génotypes de l'étude ;
- des GT finaux diploïdes et phasés ;
- aucune erreur mendélienne introduite ;
- une cible unique concordant avec le génotype moléculaire explicite audité.

## Attribution et confiance

Pour chaque individu, `0|1` attribue l'allèle ALT à `H2`, `1|0` à `H1`,
`1|1` aux deux haplotypes et `0|0` à aucun. Cette attribution provient du GT
phasé et n'utilise jamais le statut clinique ou le groupe.

L'option expérimentale `--score-singletons` ajoute le FORMAT `PP`, compris
entre `0,5` et `1`, pour les singletons. Un hétérozygote est fiable si son PP
atteint `minimum_phase_confidence`, fixé à `0,9` par défaut. Un score inférieur
ou absent conserve l'attribution mais la marque `UNRELIABLE`, crée une zone
non fiable au variant cible et exige une revue manuelle. Les homozygotes ne
nécessitent pas de score de phase.

Les trios distinguent une orientation `DIRECT`, `SWAPPED` ou `AMBIGUOUS`; les
duos sont signalés `DUO_COMPATIBLE` sans inventer l'origine du second
haplotype.

Les contrôles mendéliens distinguent les appels diploïdes complets des appels
contenant un allèle manquant (`.`). Une transmission avec un génotype manquant
est comptée `NOT_EVALUATED` et n'est jamais transformée en erreur mendélienne.
Un allèle autre que `0`, `1` ou `.` reste bloquant. Les nombres de transmissions
évaluables et non évaluables sont publiés avant et après phasage.

SHAPEIT5 peut compléter temporairement un appel manquant pour construire la
phase. Cette valeur n'est pas une observation : le contrôle de conservation
compare séparément tous les GT initialement observés, à l'allèle près et sans
tenir compte de l'ordre phasé. Après la passe rare, le masque de l'entrée est
réappliqué au seul champ `GT` des sorties commune et finale ; les autres champs
`FORMAT` sont conservés. Le BCF commun n'est remasqué qu'après avoir servi de
scaffold. Les deux BCF sont ensuite réindexés et relus afin de vérifier la
restauration exacte. Un variant, un échantillon ou un `GT` absent, un doublon,
un GT observé modifié ou un remasquage incomplet bloque la publication.

Le manifeste et le QC publient le nombre de GT manquants en entrée, le nombre
de complétions internes SHAPEIT5, les nombres remasqués dans les BCF commun et
final, ainsi que l'indicateur explicite
`completed_genotypes_published_as_observed: false`. Le variant cible reste
soumis à une exigence plus stricte : tous ses GT doivent être complets et
évaluables avant comme après le phasage.

## Sorties

- `common.phased.bcf` et son index ;
- `target.phased.bcf` et son index ;
- `common.phase.log` et `rare.phase.log` ;
- `carrier_haplotypes.tsv` ;
- `phasing_transmissions.tsv` ;
- `phasing_unreliable_regions.tsv` ;
- `shapeit5_phasing_manifest.json`.

Avec plusieurs threads, SHAPEIT5 avertit que la graine seule ne garantit pas
une reproduction bit à bit. Le défaut reste donc un thread.
