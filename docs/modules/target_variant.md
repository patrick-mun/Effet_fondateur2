# Étape 04 — jeu du variant cible

## Responsabilité

`04_prepare_target_variant_dataset` injecte un variant moléculaire confirmé dans
une copie du jeu PLINK du chromosome cible produit par l'étape `03`. Le jeu de
base n'est jamais modifié.

L'étape dépend directement de `02_build_sample_registry` pour le registre
approuvé et de `03_convert_acpa` pour le triplet BED/BIM/FAM cible et son
descripteur. Elle reste bloquée tant que `sample_registry_approval` est en
attente.

## Entrées

`inputs.target_variant_metadata` désigne un document YAML conforme à
`schemas/target_variant_metadata.schema.json`. Il contient assemblage,
coordonnée, REF/ALT, transcrit, HGVS c./p., identifiant projet et deux
confirmations auditables : référence génomique et annotation transcript/HGVS.

`inputs.target_variant_genotypes` désigne une table TSV conforme à
`schemas/target_variant_genotypes.schema.json`. Chaque échantillon du registre
doit apparaître exactement une fois avec :

- un génotype explicite composé uniquement de REF et ALT ;
- une source individuelle documentée ;
- une méthode de mesure ;
- aucune provenance issue du statut clinique ou du groupe.

Si un génotype est déjà présent dans `samples.master.tsv`, sa valeur et sa
source doivent correspondre à la table dédiée.

## Traitement et contrôles

L'étape vérifie les empreintes du registre, du descripteur et des trois fichiers
PLINK de base. Elle refuse un variant déjà présent, une coordonnée déjà occupée,
un ordre FID/IID différent ou une définition moléculaire incompatible avec la
configuration résolue.

Un jeu PLINK monomarque temporaire est construit depuis les génotypes
individuels, puis fusionné par PLINK avec le jeu du chromosome cible. Les
effectifs, l'ordre des individus, le nombre de variants, l'en-tête BED et les
allèles du variant injecté sont contrôlés après fusion.

PLINK `--mendel` est limité au variant cible. Une incompatibilité mendélienne
bloque la publication avec le code V2 `4`. En l'absence de trio complet, le
contrôle est enregistré `NOT_APPLICABLE`, et non assimilé à un succès testé.
Un échec, timeout ou résultat PLINK malformé utilise le code V2 `3`.

## Sorties

- `target_variant.bed/.bim/.fam` et leur descripteur versionné ;
- `target_genotype_audit.tsv`, sensible et individuel ;
- `target_variant_mendel.tsv` et le rapport PLINK brut ;
- `target_variant_injection_report.json`, avec empreintes avant/après ;
- `stage_outputs.json`, `audit.json` et `checksums.sha256`.

Les fichiers intermédiaires PLINK sont supprimés uniquement après validation.
Une tentative bloquée reste dans `stages/attempts/` pour diagnostic et n'est
jamais publiée comme résultat scientifique.

## Dépendances

Cette étape n'ajoute aucun package. Elle utilise Python 3.12, PyYAML et
jsonschema déjà déclarés par le projet, ainsi que PLINK 1.9. KING, bcftools, R,
Gamma, les adaptateurs de phasage et d'IBD local ne sont pas utilisés.
