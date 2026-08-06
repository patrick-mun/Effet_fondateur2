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
