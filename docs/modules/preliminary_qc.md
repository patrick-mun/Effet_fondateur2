# Étape 05 — QC préliminaire genome-wide

## Responsabilité

`05_qc_preliminary` produit un jeu PLINK autosomique techniquement fiable pour
la construction du panel indépendant de l'étape `06`. Il consomme le registre
approuvé de l'étape `02` et le jeu genome-wide de l'étape `03` sans modifier
leurs fichiers.

Cette étape intervient avant la résolution de la parenté, de la structure de
population et des cohortes. Elle n'applique donc aucun filtre HWE et ne prend
aucune décision clinique ou populationnelle.

## Paramètres

Les paramètres sont enregistrés dans la configuration résolue et l'audit :

- `sample_missingness_max` : exclusion technique individuelle ;
- `variant_missingness_max` : exclusion technique d'un variant ;
- `kinship_maf_alert_threshold` : alerte préparant l'étape `06`, sans exclusion ;
- `heterozygosity_z_alert` : alerte individuelle, sans exclusion ;
- `batch_missingness_delta_alert` : alerte de différence maximale entre lots ;
- `duplicate_pi_hat_alert` : seuil du scan préliminaire de paires ;
- `duplicate_scan_min_polymorphic_variants` : information minimale requise ;
- `plink_timeout_seconds` : délai maximal de chaque commande PLINK.

Les valeurs doivent être explicites. Un changement de seuil impose un nouveau
run, conformément à l'immutabilité de `config.resolved.yaml`.

## Métriques et décisions

PLINK produit les rapports de données manquantes par individu et variant, les
fréquences alléliques et l'hétérozygotie autosomique. Un second rapport
`--missing --within` stratifie chaque variant par `ARRAY_BATCH`; une valeur
absente dans le registre devient le lot explicite `UNREPORTED`.

Les exclusions automatiques sont limitées aux seuils de données manquantes.
Elles sont figées dans des listes `--remove` et `--exclude` calculées depuis les
rapports initiaux, puis appliquées à une copie PLINK. Cette procédure évite que
l'ordre interne de recalcul de plusieurs filtres modifie la décision auditée.

Les faibles MAF, valeurs extrêmes d'hétérozygotie et différences de missingness
entre lots restent des alertes. Elles ne retirent aucun échantillon ou variant.

Le scan `--genome` des duplicats potentiels n'est lancé que si le nombre minimal
de variants polymorphes est atteint. Une paire détectée reste une alerte : la
décision d'apparentement appartient à l'étape `07` sur un panel indépendant.

## Contrôles non applicables

- Le HWE est volontairement absent avant le gel des témoins indépendants.
- La discordance de sexe est `NOT_EVALUATED`, car le jeu canonique de l'étape
  `03` contient uniquement les autosomes.
- La comparaison entre lots est `NOT_EVALUATED` lorsqu'un seul lot est présent.
- Le scan de duplicats est `NOT_EVALUATED` si le jeu est trop peu polymorphe.

Ces absences sont enregistrées dans `qc_alerts.tsv` et `audit.json`; elles ne
sont jamais présentées comme des contrôles réussis.

## Sorties

- `genomewide_pre_qc.bed/.bim/.fam` et son descripteur versionné ;
- `qc_individual_metrics.tsv` ;
- `qc_variant_metrics.tsv` ;
- `qc_batch_summary.tsv` ;
- `qc_alerts.tsv` ;
- `qc_preliminary_report.json` ;
- `stage_outputs.json`, `audit.json` et `checksums.sha256`.

Les visualisations attendues sont déclarées dans l'audit mais restent séparées
du calcul ; elles seront produites par l'étape consolidée `18` depuis les tables
du run courant.

## Dépendances

Cette étape n'ajoute aucun package. Elle utilise Python 3.12, PyYAML et
jsonschema déjà présents, ainsi que PLINK 1.9. KING, bcftools, R, Gamma et les
adaptateurs de phasage ou d'IBD local ne sont pas utilisés.
