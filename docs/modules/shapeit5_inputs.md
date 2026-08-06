# Étape 12.5 — Entrées validées pour SHAPEIT5

## Responsabilité

La sous-étape `12.5`, intégrée à `12_phase_target_region`, transforme le jeu
PLINK régional et les artefacts harmonisés de `12.4` en entrées prêtes pour
l'adaptateur `shapeit5_phase_common_rare_v1`. Elle ne lance pas SHAPEIT5.

## Sélection des variants

Les variants communs dont l'identité canonique est confirmée par la référence
sont conservés pour `phase_common`. Le variant cible est toujours conservé : il
est marqué `COMMON_TARGET` s'il appartient à la référence ou `RARE_TARGET` s'il
est propre à l'étude. Les autres variants propres à l'étude sont exclus, car
leur orientation REF/ALT ne peut pas être confirmée par le panel public.

Le VCF d'étude est produit par PLINK avec les identifiants maître `SAMPLE_ID`.
L'allèle REF canonique est imposé avec `--a2-allele`; PLINK interdit de combiner
`--a1-allele` et `--a2-allele`. Le pipeline recontrôle ensuite avec bcftools
l'ordre exact des individus et l'identité `CHROM/POS/ID/REF/ALT` de chaque
variant. Cette vérification confirme aussi ALT sans l'imposer séparément.
BCFtools calcule ensuite les tags INFO `AC` et `AN`, exigés par
`SHAPEIT5_phase_common`, puis le pipeline vérifie qu'ils sont présents et
cohérents pour chaque variant.

## Carte et pedigree

La carte génétique SHAPEIT est un fichier gzip déterministe à trois colonnes :
`pos`, `chr`, `cM`. Les positions physiques doivent être strictement
croissantes, les cM monotones et la carte doit couvrir tous les variants
conservés. Aucune approximation `1 cM/Mb` n'est utilisée.

Le pedigree est sans en-tête et contient `enfant`, `père`, `mère`. Les trois
valeurs utilisent les mêmes `SAMPLE_ID` que le VCF et `NA` représente un parent
absent du jeu d'étude. Un enfant n'est écrit que si au moins un parent est
présent. Le fichier peut donc être vide ; `12.6` devra alors omettre l'argument
pedigree.

## Sorties

Le dossier `shapeit5_inputs/` publie atomiquement :

- `study.shapeit5.vcf.gz` et son index tabix ;
- `shapeit5.genetic_map.tsv.gz` ;
- `shapeit5.pedigree.tsv` ;
- `shapeit5_variant_selection.tsv` ;
- `shapeit5_sample_mapping.tsv` ;
- `shapeit5_inputs_manifest.json`.

Le manifeste lie ces fichiers aux empreintes de la référence harmonisée,
enregistre les versions PLINK/bcftools, les régions, les effectifs et le rôle du
variant cible. Toute collision d'identifiant, discordance d'ordre, absence du
variant cible, absence de variant commun ou divergence allélique bloque avant
le phasage.

## Limites

- Les entrées GRCh38 doivent déjà utiliser une représentation minimale ; aucune
  FASTA ni correction automatique de gauche-alignement n'est ajoutée.
- Aucun génotype cible n'est inféré depuis le statut clinique ou le groupe.
- La séparation `phase_common`/`phase_rare`, l'attribution du chromosome
  porteur et le QC du phasage relèvent de `12.6–12.7`.

Références officielles :

- [`phase_common`](https://odelaneau.github.io/shapeit/docs/documentation/phase_common/)
- [`phase_rare`](https://odelaneau.github.io/shapeit/docs/documentation/phase_rare/)
