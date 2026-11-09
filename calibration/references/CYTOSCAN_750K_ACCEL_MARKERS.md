# Référence des marqueurs CytoScan 750K Accel — D11b

## Décision de source

L'univers observable du panel simulé doit provenir du paquet constructeur
`CytoScan 750K Accel Array Analysis Files-hg38-NA36.r1`. Il ne doit jamais être
reconstruit depuis les exports ACPA d'un patient, d'un témoin ou d'une famille.

La page produit Thermo Fisher du modèle Accel est la source fournisseur :

- <https://www.thermofisher.com/order/catalog/product/952539>
- paquet attendu : `CytoScan 750K Accel Array Analysis Files-hg38-NA36.r1` ;
- fichier fournisseur :
  `CytoScan750K Accel_Array_Analysis_Files_hg38_NA36.zip`, version `r1`,
  environ 509 MB ;
- conception : `CytoScan750K_Accel` ;
- assemblage : `GRCh38` ;
- version d'annotation : `NA36.r1`.
- SHA-256 du paquet acquis et vérifié le 29 octobre 2026 :
  `6a696ef3d2e61e9db74ead91c6410e15efc90da6a8bfc814f5da034e364607a3`.

L'ancienne `CytoScan 750K Array` utilise un autre dessin et le paquet
`NA36.r7`. Elle ne peut donc pas être substituée à la version Accel.

Le fichier est déclaré `secureAsset` par l'API Thermo Fisher. Son URL
authentifiée est :

<https://sec-assets.thermofisher.com/TFS-Assets/LSG/Support-Files/CytoScan750K%20Accel_Array_Analysis_Files_hg38_NA36.zip>

Le téléchargement nécessite une connexion Thermo Fisher ; aucune identité,
session ou information de compte ne doit être enregistrée dans le dépôt ou
dans un manifest.

## Conservation locale et provenance

Le paquet fournisseur et toute table normalisée restent hors Git, dans un
dossier public de référence ou sous `calibration_workspace/`. Le paquet de
référence est conservé hors du dépôt et hors du cache du pipeline. D11b le lit
directement et écrit sa table normalisée dans son propre dossier privé.

La base constructeur retenue est
`CytoScan750k_Accel_Array.na36.annot.db`, table `Annotations`. Son enregistrement
`Information` identifie `CytoScan750k_Accel`, `hg38` et le build NetAffx `36`.
Elle contient 960 877 annotations, dont 254 822 SNP et 247 777 SNP
autosomiques normalisables. Le fichier
`CytoScan750K_Accel_Array.snplist.r3.txt` ne contient que 10 000 SNP destinés
au calcul du SNPQC ; il ne doit jamais servir d'univers observable.

Avant usage, le matérialisateur :

1. calcule le SHA-256 du paquet officiel téléchargé ;
2. exige l'entrée D09
   `thermofisher_cytoscan750k_accel_na36_r1` avec la même empreinte ;
3. vérifie l'intégrité SQLite, le schéma, les métadonnées et les effectifs ;
4. extrait uniquement les SNP autosomiques génotypables et leurs annotations
   de brin ;
5. normalise la table au format ci-dessous, sans colonne d'échantillon ;
6. calcule son SHA-256 et la conserve comme artefact privé D11b.

La redistribution du paquet fournisseur ou de sa table dérivée dans le dépôt
n'est pas autorisée. Le paquet renvoie aux conditions Affymetrix/Thermo Fisher ;
les conditions officielles revues le 29 octobre 2026 permettent le
téléchargement pour un usage personnel non commercial et interdisent la
distribution publique des matériaux. La décision du projet est donc limitée à
la calibration locale, sans redistribution ni revendication de licence plus
large. Le code, le schéma et les petites fixtures entièrement synthétiques sont
seuls versionnés.

Conditions revues :
<https://www.thermofisher.com/us/en/home/global/terms-of-use.html>.

## Contrat du TSV normalisé

Le séparateur est une tabulation et l'en-tête est strictement :

```text
MARKER_ID	ARRAY_DESIGN	DESIGN_VERSION	ASSEMBLY	CHROMOSOME	POSITION_BP	ALLELE_A	ALLELE_B	STRAND
```

Chaque ligne doit respecter les règles suivantes :

- `MARKER_ID` : identifiant de sonde SNP du dessin, unique sans tenir compte de
  la casse ;
- `ARRAY_DESIGN` : `CytoScan750K_Accel` ;
- `DESIGN_VERSION` : `NA36.r1` ;
- `ASSEMBLY` : `GRCh38` ;
- `CHROMOSOME` : autosome `1` à `22` ;
- `POSITION_BP` : coordonnée GRCh38, base 1, fermée ;
- `ALLELE_A` et `ALLELE_B` : deux bases distinctes parmi `A`, `C`, `G`, `T` ;
- `STRAND` : `+`, `-` ou `UNKNOWN`.

Le fichier ne doit contenir ni génotype, ni appel, ni statut clinique, ni nom
de fichier d'échantillon, ni identifiant pouvant être relié à un individu.

## Règles d'harmonisation

La vérité D11a reste dense et orientée comme le VCF 1000G. D11b sélectionne
les variants observables selon l'intersection exacte chromosome–position
GRCh38 puis vérifie les allèles directs, inversés ou complémentaires selon le
brin. Les cas suivants sont exclus et comptés :

- position absente de 1000G ;
- plusieurs variants 1000G à la même position ;
- plusieurs marqueurs de puce à la même position ;
- allèles incompatibles ;
- paire palindromique `A/T` ou `C/G` avec brin inconnu.

Le VCF observable conserve l'orientation REF/ALT de 1000G. La relation avec
les allèles de la sonde est enregistrée dans
`chip_reference_marker_map.tsv`; elle devra être appliquée explicitement par
l'export ACPA ultérieur.

## Porte avant D12

D12 ne peut commencer qu'après un run D11b réel dont :

- le SHA-256 du paquet correspond et son contenu SQLite est validé ;
- l'entrée Thermo Fisher du catalogue public est présente et liée à l'empreinte
  du paquet ;
- la table normalisée est produite directement par D11b et liée par SHA-256 ;
- le manifest passe son schéma ;
- le nombre de variants du VCF observable égale le nombre de marqueurs liés ;
- les exclusions et relations alléliques sont auditées ;
- aucune source dérivée d'un échantillon réel n'a été lue.
