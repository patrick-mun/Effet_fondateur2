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

Aucun gauche-alignement dépendant d'une séquence FASTA n'est inventé. Le contrat
exige donc déjà des coordonnées et allèles minimaux dans le panel GRCh38. Une
future acquisition de FASTA devra passer par un catalogue et un cache immuables
avant d'autoriser une normalisation dépendante de cette séquence.

## Statuts d'harmonisation

Une ligne de `variant_harmonization.tsv` est produite dans l'ordre du BIM :

- `MATCHED_DIRECT` : A1/A2 correspondent à REF/ALT ;
- `MATCHED_SWAPPED` : A1/A2 correspondent à ALT/REF ;
- `STUDY_ONLY` : aucun variant de référence à cette position ;
- `ALLELE_MISMATCH` : position commune, allèles différents ;
- `AMBIGUOUS_REFERENCE_MATCH` : plusieurs enregistrements de référence possibles ;
- `DUPLICATE_STUDY_MATCH` : plusieurs marqueurs ACPA visent la même clé canonique.

Seuls `MATCHED_DIRECT` et `MATCHED_SWAPPED` sont éligibles au scaffold commun.
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
index ou du manifeste `12.3` bloque avant la publication. Le BIM et le manifeste
de phasage sont validés ensemble et leurs SHA-256 sont enregistrés ; leur
comparaison aux descripteurs d'artefacts amont sera ajoutée lors de l'intégration
de cette fonction à l'étape `12`.
