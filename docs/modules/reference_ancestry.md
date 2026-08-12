# Positionnement par rapport aux références d'ascendance

## Portée générale

L'étape planifiée `16A_analyze_reference_ancestry` est générique pour toute
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
