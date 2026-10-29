# Sous-système de calibration

Ce dossier contient le code de configuration, d'isolation, de simulation, de
vérité, d'évaluation et de production des rapports de calibration.

`isolation.py`, introduit à l'étape D08, contrôle les chemins physiques, liens
symboliques, références de manifeste et collisions entre panels avant toute
lecture d'un artefact. `partitioning.py`, introduit à l'étape D10, répartit de
façon déterministe les groupes indépendants entre donneurs de calibration,
donneurs de validation et référence de phasage commune. Il contrôle la
couverture de la source, les intersections, les empreintes et la confidentialité
du manifest sans écrire d'identifiant.

`materialize_partitions.py` vérifie les empreintes des métadonnées D09, refuse
l'écrasement d'un run et écrit atomiquement les trois listes privées ainsi que
leur manifest agrégé dans `calibration_workspace/`. Sa sortie standard ne
contient que des effectifs et statuts de contrôle.

`simulation/base_haplotypes.py`, introduit en D11, réapparie les copies
haplotypiques de donneurs distincts au sein de chaque strate publique. La
bijection conserve exactement fréquences alléliques et LD haplotypique ; les
génotypes non phasés sont ensuite dérivés des deux copies sans approximation.
`materialize_base_haplotypes.py` vérifie les MD5 fournisseur, les SHA-256
locaux et les liaisons aux catalogues D09, extrait les SNV bialléliques
complètement phasés avec `bcftools`, interpole la carte sans extrapolation et
écrit atomiquement les quatre artefacts privés D11.

`simulation/chip_markers.py` et `materialize_chip_observation.py` constituent
la couche D11b. Ils vérifient le paquet officiel CytoScan 750K Accel
`NA36.r1` et le manifest dense D11a. `simulation/chip_package.py` ouvre en
lecture seule sa base d'annotations SQLite, valide son identité et produit la
table SNP autosomique normalisée dans le workspace privé. D11b construit ensuite
l'intersection chromosome–position–allèles, un VCF observable et l'audit des
exclusions. Les dossiers connus d'exports d'échantillons sont refusés avant
lecture ; la vérité dense reste inchangée et le cache du pipeline n'est pas lu.

`simulation/population_structure.py` et
`materialize_population_structure.py` implémentent D12. Ils exigent des
manifests D11a/D11b finalisés, liés au même panel et au même domaine. Des paires
de diploïdes D11 provenant de deux strates donnent deux mélanges réciproques
`A1|B2` et `B1|A2`, ce qui conserve chaque haplotype exactement une fois. Sans
recombinaison, seules les proportions pures ou 50/50 sont exactes ; les
proportions inégales sont différées à D14. La vérité et la provenance restent
privées, tandis que le VCF observable ne contient aucune étiquette d'ascendance.
