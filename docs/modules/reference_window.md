# Étape 12.3 — Fenêtre phasée de référence

## Responsabilité

`extract_reference_window` extrait avec `bcftools` la région cible depuis le VCF
chromosomique immuable produit par `12.2`. Cette sous-étape conserve tous les
échantillons et toutes les populations du panel ; elle ne sélectionne ni
population, ni individu, ni variant commun avec la puce ACPA.

L'harmonisation allélique avec les marqueurs de l'étude appartient à `12.4`.

## Entrées et sortie

L'entrée associe obligatoirement la résolution du catalogue `12.1` et l'entrée
de cache vérifiée `12.2`. Avant d'appeler `bcftools`, les identifiants, chemins,
tailles et SHA-256 du cache sont recalculés et comparés au manifeste immuable.

La sortie contient :

- `reference.chrN.START-END.vcf.gz`, compressé bgzip ;
- son index tabix `.tbi` ;
- `reference_window_manifest.json`, validé par
  `schemas/reference_window_manifest.schema.json`.

Les trois fichiers sont préparés dans un dossier temporaire du même système de
fichiers puis publiés ensemble par renommage atomique. Une sortie existante
n'est jamais réutilisée ni remplacée silencieusement.

## Contrôles bloquants

L'extraction exige :

- un autosome GRCh38 et des bornes incluses dans sa longueur canonique ;
- le contig `chrN` avec la longueur GRCh38 attendue dans les en-têtes source et
  extrait ;
- un champ `FORMAT/GT` ;
- exactement l'effectif déclaré par le catalogue, soit `3 202` pour le panel
  principal ;
- les mêmes identifiants d'échantillons dans le même ordre avant et après
  extraction ;
- au moins un variant dans la fenêtre ;
- aucun génotype appelé utilisant le séparateur non phasé `/` ;
- un index tabix produit et vérifiable par `bcftools`.

Les génotypes entièrement manquants ne sont pas interprétés comme une preuve de
phase. Ils sont conservés ; leur traitement scientifique interviendra dans les
contrôles et l'harmonisation suivants.

## Reproductibilité

Le manifeste enregistre la release, le catalogue, la clé et les empreintes du
cache, la région exacte en bp, l'effectif, le nombre de variants, les SHA-256
des sorties et la version de `bcftools`. L'extraction n'altère jamais le cache.

Paramètres de configuration :

```yaml
phase_target_region:
  parameters:
    reference_extract_timeout_seconds: 7200
    reference_extract_threads: 1
```

Les tests unitaires simulent `bcftools` avec des sorties minimales. Un smoke
local avec `bcftools 1.21` valide l'expression de détection des génotypes non
phasés. Aucun VCF 1000 Genomes complet n'est téléchargé pendant `12.3`.
