# Étape 12.1 — Catalogue du panel phasé

## Responsabilité

La sous-étape `12.1` sélectionne de manière déterministe un panel public phasé
compatible avec l'assemblage et le chromosome de l'étude. Elle ne télécharge
aucun VCF et ne lit aucune donnée génétique privée.

Le catalogue versionné est
`config/references/1000g_high_coverage_grch38_phased.json`. Son schéma est
`schemas/reference_panel_catalog.schema.json`.

## Panel retenu

Le panel `1kg_3202_high_coverage_20220422` correspond aux appels phasés
SNV/INDEL/SV de haute couverture produits pour `3 202` échantillons 1000
Genomes sur `GRCh38`. Le profil principal conserve tous les échantillons de la
release et toutes les populations disponibles ; aucune population n'est
sélectionnée à cette étape.

Le catalogue contient pour chacun des 22 autosomes :

- le nom exact du VCF bgzip et de son index tabix ;
- leurs URL HTTPS sous le domaine officiel EBI/IGSR ;
- les MD5 publiés dans `20220804_manifest.txt` ;
- l'URL et le SHA-256 observé du manifeste officiel ;
- l'URL et le SHA-256 observé du README de la release.

Le SHA-256 du catalogue lui-même accompagne chaque résolution. Une modification
du catalogue change donc l'identité de la référence résolue.

## Résolution

`resolve_reference_panel` exige :

- l'identifiant exact du panel ;
- l'assemblage exact `GRCh38` ;
- un autosome entre `1` et `22`.

Le chargement bloque un doublon ou un autosome manquant, un MD5 mal formé, un
panel ambigu et toute URL hors de `ftp.1000genomes.ebi.ac.uk`. La résolution est
entièrement hors ligne : elle ne remplace jamais silencieusement une release ou
un miroir indisponible.

## Frontière avec 12.2

Les MD5 du catalogue sont les empreintes publiées par le fournisseur. `12.2`
devra télécharger atomiquement le VCF, l'index et les métadonnées, vérifier les
MD5 officiels, calculer leurs SHA-256 locaux et les placer dans un cache
immuable. `12.1` ne garantit pas encore qu'un fichier distant est disponible ou
que son contenu local est intact.

Sources officielles :

- [collection haute couverture 1000 Genomes](https://www.internationalgenome.org/data-portal/data-collection/1000genomes_30x)
- [FAQ IGSR sur les appels phasés](https://www.internationalgenome.org/faq/are-the-variant-calls-in-igsr-phased/)
- [répertoire FTP de la release](https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/working/20220422_3202_phased_SNV_INDEL_SV/)
- [manifeste MD5 officiel](https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/working/20220422_3202_phased_SNV_INDEL_SV/20220804_manifest.txt)
