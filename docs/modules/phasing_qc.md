# Étape 12.7 — Publication du QC de phasage

## Responsabilité

La sous-étape `12.7` clôt l'étape 12 sans recalculer le phasage. Elle vérifie
les empreintes des sorties `12.6`, consolide les contrôles scientifiques et
publie un résumé agrégé des haplotypes. Les données individuelles restent dans
`carrier_haplotypes.tsv`; le résumé final ne contient aucun identifiant.

## Contrôles consolidés

`phasing_qc.tsv` contient douze contrôles ordonnés :

- intégrité des entrées et présence des tags `AC/AN` ;
- nombres de variants du scaffold commun et du BCF final ;
- ordre des individus et conservation des génotypes ;
- conservation unique du variant cible et concordance moléculaire explicite ;
- absence d'erreur mendélienne avant et après phasage ;
- transmissions pedigree, ou `NOT_APPLICABLE` sans pedigree ;
- nombre de porteurs dont la confiance atteint le seuil configuré.

Une confiance insuffisante est un `WARN`, jamais un succès silencieux ni une
suppression automatique. Elle est répercutée dans le résumé et impose une revue
manuelle. Toute incohérence de compte, d'empreinte ou de warning bloque la
publication.

## Sorties

- `phasing_qc.tsv` : contrôles observés, attendus et codes de détail ;
- `carrier_haplotype_summary.tsv` : effectifs agrégés par `NONE/H1/H2/BOTH` et
  statut de fiabilité ;
- `phase_target_region_summary.json` : synthèse finale, empreintes des deux
  manifestes sources, effectifs et décision de revue.

La publication est atomique. La reprise de l'étape réutilise le dossier publié
uniquement si la signature et toutes les empreintes restent valides.

## Visualisations attendues

Les figures individuelles ou familiales ne sont pas générées ici afin de ne
pas anticiper l'étape `18`. L'audit déclare les produits attendus : haplotypes
phasés, transmissions et confiance. Ils devront employer les pseudonymes et les
seules sorties du run courant.
