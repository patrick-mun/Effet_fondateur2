# Positionnement par rapport aux références d'ascendance

## Portée générale

L'étape `16A_analyze_reference_ancestry`, méthode
`reference_only_global_local_pca_v1`, est générique pour toute
variation cible autosomique décrite par la configuration du run. `DOCK6` et la
coordonnée du chromosome 19 constituent un cas d'étude, pas une dépendance de
la méthode.

L'analyse ne doit contenir aucun chromosome, gène, identifiant de variant,
position ou allèle codé en dur. Elle résout systématiquement la cible depuis
les métadonnées moléculaires et les artefacts versionnés du run. La fenêtre
locale est celle préparée autour de cette cible par les étapes régionales du
pipeline ; ses bornes physiques et génétiques doivent être publiées dans
l'audit.

## Deux analyses séparées

La branche globale ajuste une PCA sur les seuls individus 1000 Genomes non
apparentés retenus comme références. Les individus de l'étude sont projetés
dans ces axes et ne peuvent ni les orienter ni les recalculer.

La branche locale travaille sur les haplotypes phasés de la région entourant
la variation cible configurée. Elle harmonise les variants de cette région avec
la référence, ajuste le modèle sur les haplotypes de référence, puis projette
séparément les haplotypes de l'étude. Les chromosomes porteurs et non porteurs
de l'allèle cible doivent rester identifiables dans les sorties pseudonymisées.
Changer de variation cible ou de région doit donc suffire à réutiliser la même
implémentation sans modification du code.

Cette branche mesure un positionnement haplotypique relatif dans une fenêtre.
Elle ne constitue ni une méthode générale d'inférence d'ascendance locale par
segments, ni une attribution ethnique, ni une preuve d'ancêtre généalogique ou
d'identité par descendance.

## Référence et cache

Les métadonnées des 3 202 individus 1000 Genomes permettent de documenter les
populations de référence. L'ajustement principal est limité aux 2 504 individus
non apparentés de la liste officielle. Le premier passage télécharge et publie
atomiquement uniquement les extraits nécessaires aux variants analysés. Chaque
entrée de cache est liée à l'identité et au MD5 officiels du VCF source ainsi
qu'à un SHA-256 local ; un passage ultérieur doit la vérifier et la réutiliser
sans réseau. Une entrée absente en mode hors ligne ou une empreinte discordante
bloque l'étape.

Pour éviter les lectures HTTP indexées, les VCF complets et leurs index peuvent
être exposés localement dans
`data/cache/references/source_panels/<panel_id>/`, en conservant exactement les
noms du catalogue. Le pipeline recalcule les MD5 du VCF et du TBI et bloque si
l'un d'eux ne correspond pas à la publication officielle. L'URL officielle
reste l'identité auditée de la source et la copie locale ne change donc ni la
clé ni le contenu scientifique du cache d'extraits.

Le cache contient uniquement les variants demandés et les 2 504 références
publiques. Aucun génotype, identifiant ou fichier de l'étude n'est transmis au
serveur 1000 Genomes. Sa clé lie le panel, l'assemblage, le chromosome, l'URL et
les MD5 officiels du VCF et de son index, ainsi que les SHA-256 des positions et
de la liste d'échantillons. Les fichiers publiés sont en lecture seule et leur
SHA-256 est revalidé à chaque accès.

## Harmonisation et modèle

La PCA globale utilise les variants autosomiques du panel indépendant de
l'étape 06. PLINK exporte seulement les dosages de l'étude ; bcftools extrait
les mêmes positions chez les références publiques, chromosome par chromosome.
La PCA locale utilise le BCF final phasé de l'étape 12 et décompose chaque
individu en `H1` et `H2`. La variation cible doit être présente et complètement
évaluable dans le BCF d'étude, mais elle peut être absente de la référence et
ne doit pas être inventée.

L'harmonisation exige une égalité de coordonnée GRCh38 et une concordance
directe ou inversée REF/ALT. Une inversion corrige le dosage ; les compléments
de brin et les correspondances fondées uniquement sur un identifiant sont
refusés. Lorsque la référence publie plusieurs enregistrements bialléliques à
une même position, la paire REF/ALT de l'étude doit en désigner exactement un ;
une absence ou plusieurs correspondances compatibles restent bloquantes. Les
variants monomorphes, trop manquants, absents ou incompatibles sont
audités séparément. Les valeurs manquantes restantes sont imputées uniquement
à la fréquence de la référence pour le calcul des coordonnées PCA ; elles ne
deviennent jamais des observations publiées.

Les fréquences, échelles et axes sont ajustés exclusivement sur les références.
Les individus de l'étude, puis leurs haplotypes locaux, sont projetés dans ce
modèle figé. Les centroïdes de population et superpopulation sont descriptifs.

## Sorties versionnées

- `ancestry_scores.tsv` : références utilisées et entités d'étude projetées,
  globales ou locales, avec copie porteuse explicitement distinguée ;
- `ancestry_eigenvalues.tsv` et `ancestry_variant_loadings.tsv` : modèle de
  référence reproductible ;
- `ancestry_variant_audit.tsv` : chaque variant candidat et sa décision ;
- `ancestry_population_centroids.tsv` : repères agrégés 1000 Genomes ;
- `reference_ancestry_summary.json` : effectifs, région, cache, contrôles et
  interdictions d'interprétation ;
- `audit.json`, `stage_outputs.json` et `checksums.sha256` : provenance et
  intégrité orchestrées.

Les coordonnées individuelles restent `sensitive_genetic`. Les figures 18
séparent obligatoirement `REFERENCE_ANCESTRY_GLOBAL` et
`REFERENCE_ANCESTRY_LOCAL`, et l'étape 17 compare le domaine technique
`REFERENCE_ANCESTRY` sans score composite.

## Contrôles obligatoires

L'étape doit bloquer en cas d'assemblage discordant, de cible absente, de
variant ou d'échantillon dupliqué, d'allèles incompatibles, d'ordre
d'haplotypes ambigu, de référence apparentée dans l'ajustement principal ou de
projection calculée avec un modèle différent de celui publié.

L'audit doit enregistrer au minimum la variation cible effectivement utilisée,
les bornes de la région, les variants candidats, harmonisés et rejetés, les
effectifs de référence et d'étude, les composantes retenues, les empreintes du
catalogue et du cache, et la confirmation que les échantillons de l'étude n'ont
pas influencé les axes.

## Interprétation

Les populations 1000 Genomes sont des repères externes larges. Elles ne
représentent pas exhaustivement les histoires démographiques réunionnaise,
malgache ou comorienne. Les proximités avec leurs centroïdes sont donc
exploratoires et doivent être rapportées comme un positionnement relatif, sans
étiquette identitaire individuelle. Les résultats globaux et locaux doivent
rester séparés : une proximité globale ne détermine pas l'origine d'un
haplotype local, et une proximité locale ne prouve pas un effet fondateur.
