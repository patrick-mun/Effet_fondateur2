# Étape 12.0 — Contrat de phasage SHAPEIT5

## Choix figé

L'adaptateur `shapeit5_phase_common_rare_v1` utilise SHAPEIT5 `5.1.1`, sous
licence MIT. Il cible les autosomes en `GRCh38` et produit des BCF. La version
est exacte : un autre numéro bloque l'adaptateur jusqu'à une validation et une
mise à jour explicites du contrat.

Deux exécutables composent l'adaptateur :

- `SHAPEIT5_phase_common` construit le scaffold des marqueurs communs partagés
  entre l'étude et le panel phasé de référence ;
- `SHAPEIT5_phase_rare` phase ensuite les variants rares de l'étude, dont le
  variant cible lorsqu'il est absent du panel, sur ce scaffold.

La documentation officielle précise que `phase_common --reference` ne considère
pas les sites de l'étude absents du panel. Le pipeline ne doit donc jamais
prendre la sortie commune comme preuve que le variant cible a été conservé.

## Entrées et paramètres

Le contrat exige :

- des génotypes d'étude en VCF, VCF bgzip ou BCF ;
- une référence phasée dans l'un de ces formats ;
- une carte génétique SHAPEIT correspondant à `GRCh38` ;
- une région autosomique explicite ;
- une sortie BCF et un journal distincts ;
- une graine positive explicite, `15052011` par défaut ;
- un nombre de threads positif ;
- `Ne=15000` par défaut pour `phase_rare` ;
- un pedigree optionnel au format `enfant père mère`, avec `NA` pour un parent
  inconnu.

La carte est obligatoire même si SHAPEIT5 sait utiliser une approximation de
`1 cM/Mb` lorsqu'elle manque. Cette approximation reste interdite par le
pipeline. Les étapes `12.1–12.5` vérifient l'assemblage, les allèles, les
index, les individus et les régions avant de construire ces commandes.

## Configuration

```yaml
tools:
  phasing_adapter:
    adapter_id: shapeit5_phase_common_rare_v1
    phase_common_command: SHAPEIT5_phase_common
    phase_rare_command: SHAPEIT5_phase_rare
    expected_version: 5.1.1
```

Les commandes peuvent être des chemins absolus. Sur la machine de développement,
elles sont disponibles après `conda activate effet-fondateur-shapeit5`.

## Limites du contrat 12.0

Cette sous-étape ne télécharge aucune référence, ne convertit pas PLINK en BCF et
ne lance aucun phasage. Le catalogue 1000 Genomes, le cache, l'harmonisation et
les fichiers de pedigree sont implémentés dans `12.1–12.5`. L'exécution et
l'attribution du chromosome porteur sont implémentées dans `12.6`; la
publication consolidée du QC relève de `12.7`. Les contrats concrets sont
décrits dans `docs/modules/shapeit5_inputs.md` et
`docs/modules/shapeit5_execution.md`.

Références officielles :

- [SHAPEIT5 et licence](https://odelaneau.github.io/shapeit/)
- [`phase_common`](https://odelaneau.github.io/shapeit/docs/documentation/phase_common/)
- [`phase_rare`](https://odelaneau.github.io/shapeit/docs/documentation/phase_rare/)
