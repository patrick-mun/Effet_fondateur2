# Étape 12.4 — Harmonisation référence/étude

## Responsabilité

`harmonize_reference_window` prépare la fenêtre publique de `12.3` et compare
chaque marqueur de la région ACPA de `11` à la référence avec la clé scientifique
`assemblage + chromosome + position + REF + ALT`. Les rsID et identifiants de
sonde restent des annotations et ne participent jamais à la jointure.

Cette sous-étape ne construit pas encore le VCF des individus de l'étude et ne
lance pas SHAPEIT5. Ces responsabilités appartiennent respectivement à `12.5`
et `12.6`.

## Normalisation de la référence

La fenêtre est décomposée avec `bcftools norm --multiallelics -any`, puis limitée
aux SNV et indels bialléliques dont REF et ALT sont des séquences `ACGT`. Les
représentations non minimales, les allèles symboliques et les sorties vides sont
bloquants. Tous les échantillons de la référence et leur ordre doivent être
strictement conservés, l'index tabix est interrogé et les génotypes appelés
doivent rester phasés.

Aucun gauche-alignement dépendant d'une séquence FASTA n'est effectué. La
décision de projet pour ce contrat est d'exiger des coordonnées et allèles déjà
minimaux sur GRCh38. Une normalisation avec FASTA constituerait une évolution de
méthode séparée, avec catalogue, cache, comparaison de référence et nouveaux
tests scientifiques ; elle ne bloque pas `12.4`.

## Statuts d'harmonisation

Une ligne de `variant_harmonization.tsv` est produite dans l'ordre du BIM :

- `MATCHED_DIRECT` : A1/A2 correspondent à REF/ALT ;
- `MATCHED_SWAPPED` : A1/A2 correspondent à ALT/REF ;
- `MATCHED_REFERENCE_COMPLETED` : un allèle PLINK vaut `0` et l'allèle manquant
  est défini par une seule correspondance de référence possible ;
- `STUDY_ONLY` : aucun variant de référence à cette position ;
- `ALLELE_MISMATCH` : position commune, allèles différents ;
- `AMBIGUOUS_REFERENCE_MATCH` : plusieurs enregistrements de référence possibles ;
- `DUPLICATE_STUDY_MATCH` : plusieurs marqueurs ACPA visent la même clé canonique.

Les trois statuts `MATCHED_*` sont éligibles au scaffold commun. La complétion
concerne uniquement la définition REF/ALT du marqueur, jamais un génotype
individuel.
Une complémentation de brin n'est jamais appliquée automatiquement : elle est
signalée par `STRAND_COMPLEMENT_NOT_APPLIED`. Cette règle repose sur le contrat
amont `forward_strand_calls` et empêche une résolution silencieuse des SNP
palindromiques.

Le variant cible peut être `STUDY_ONLY`, car `phase_rare` doit le conserver même
s'il est absent du panel. Une collision allélique ou une correspondance ambiguë
du variant cible est bloquante.

## Sorties et reproductibilité

La publication atomique contient :

- `reference.harmonized.chrN.START-END.vcf.gz` et son index tabix ;
- `variant_harmonization.tsv`, validé par un schéma versionné ;
- `reference_harmonization_manifest.json`, avec les effectifs, versions,
  paramètres scientifiques et SHA-256 des deux chaînes d'entrée.

Une sortie existante n'est jamais remplacée. Toute modification du VCF, de son
index ou du manifeste `12.3` bloque avant la publication. L'étape orchestrée
`12_phase_target_region` compare aussi le BIM, le manifeste de phasage, les
métadonnées cible et le catalogue à leurs descripteurs d'artefacts amont.

Elle publie séparément les trois artefacts de fenêtre `12.3` et les quatre
artefacts harmonisés `12.4`, puis les référence dans `stage_outputs.json`,
`audit.json` et `checksums.sha256`. Une reprise ne les réutilise qu'après
recalcul des signatures et empreintes.

Paramètre scientifique :

```yaml
phase_target_region:
  parameters:
    minimum_common_variants: 1
```

Une valeur de production plus exigeante devra être justifiée selon la densité
de la puce et les besoins de SHAPEIT5 ; `1` sert uniquement aux fixtures
minimales et ne constitue pas un seuil biologique validé.
