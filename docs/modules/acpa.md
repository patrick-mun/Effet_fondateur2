# Conversion ACPA V2

## Responsabilité

`03_convert_acpa` consomme le registre maître approuvé et uniquement les exports
ACPA qui y sont référencés. Chaque export source est lu une fois. Les appels
autosomes sont écrits dans un fichier temporaire propre à l'échantillon, puis
réalignés dans l'ordre canonique du registre sans relire l'export.

L'étape ne modifie jamais les sources et ne renseigne aucun génotype cible. Elle
convertit seulement les appels présents dans les fichiers ChAS.

## Dépendances et blocages

- `01_validate_sources` doit avoir validé les exports ;
- `02_build_sample_registry` doit avoir publié `samples.master.tsv` ;
- `sample_registry_approval` doit être résolue ;
- PLINK 1.9 ou compatible doit être configuré et disponible dans le `PATH`.

PLINK est lancé avec une liste d'arguments et sans shell. Un timeout, un code de
retour non nul ou un triplet BED/BIM/FAM absent ou incohérent retourne le code
V2 `3`. Aucun appel à bcftools, KING, R ou Gamma n'est requis à cette étape.

## Politique des marqueurs

Les chromosomes autosomiques `1` à `22` sont conservés. Le paramètre
`marker_mode` vaut :

- `intersection` : une sonde doit être présente chez tous les échantillons ;
- `union` : une sonde absente reçoit un génotype PLINK manquant `0 0`.

Une sonde est exclue si ses coordonnées diffèrent entre exports ou si plus de
deux allèles non manquants sont observés. Un appel invalide devient explicitement
manquant et est compté dans l'audit. L'absence ou le conflit d'un rsID est annoté
mais ne provoque pas l'exclusion d'une sonde cohérente.

Les appels `Forward Strand Base Calls` restent identifiés comme allèles du brin
forward. Ils ne sont jamais présentés comme REF/ALT de GRCh38 ; cette validation
moléculaire appartient à l'étape `04`.

## Sorties

L'étape publie deux triplets PLINK avec le même ordre d'individus :

- `genomewide_base.bed/.bim/.fam` pour les 22 autosomes ;
- `target_chromosome_base.bed/.bim/.fam` pour le chromosome cible.

Chaque jeu possède un descripteur `*.dataset.json` avec assemblage, effectifs,
`sample_set_id`, `variant_set_id` et empreintes. Les tables
`sample_alignment.tsv` et `acpa_variant_audit.tsv` documentent respectivement
l'ordre des individus et la décision appliquée à chaque sonde. La table
`acpa_sample_chromosome_qc.tsv` conserve les effectifs appelés et le taux de
génotypage par individu et chromosome. Le rapport JSON agrège les lectures de
sources, génotypes invalides et lignes ignorées sans contenir les appels
individuels.

## Limites actuelles

Les graphiques de QC de conversion seront produits par la couche de
visualisation auditée prévue plus loin dans la migration V2. Les fichiers
PED/MAP temporaires et les spools sont supprimés uniquement après validation des
deux triplets binaires.
