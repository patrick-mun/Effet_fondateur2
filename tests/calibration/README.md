# Tests de calibration

Ce dossier est réservé aux futurs tests de configuration, simulation,
aveuglement, évaluation, rapports et parcours bout en bout.

Les fixtures devront être petites, entièrement synthétiques et sans
correspondance avec une personne ou une donnée de `dock6_real`.

Les premiers contrôles couvrent le schéma de configuration D05, le manifest de
vérité D06, le manifest reproductible de run D07 ainsi que la cohérence
documentaire de D03 et D04. Ils vérifient notamment le refus des clés
dupliquées, des seuils pilote prématurément renseignés, des chemins de vérité
non privés, des vérités non typées, des consommateurs manquants, des
provenances incomplètes, des commandes non portables et de toute violation
déclarée de l'aveuglement. Les tests D08 utilisent uniquement des fichiers
synthétiques temporaires et démontrent que les chemins réels, liens
symboliques ou physiques et collisions calibration–validation sont bloqués
avant l'ouverture d'un artefact. Les tests D09 valident le catalogue fermé des
références publiques, l'unicité des sources, les versions, origines HTTPS,
assemblages, droits d'usage local et empreintes, puis recalculent les SHA-256
des trois catalogues techniques liés sans télécharger de donnée génétique. Les
tests D10 démontrent sur des identifiants entièrement synthétiques que les
groupes mixtes quittent la référence de phasage, qu'aucun groupe ne traverse
les partitions, que toutes les sources sont affectées exactement une fois et
que toute incohérence de chemin, empreinte, effectif ou catalogue est refusée.
Ils contrôlent aussi la lecture des deux formats publics, la production
atomique sous permissions privées, l'absence d'identifiant dans le manifest,
le refus d'écrasement et le rejet d'une entrée corrompue.

Les tests D11b construisent une fausse archive fournisseur et une base SQLite
entièrement synthétiques. Ils vérifient la version `NA36.r1`, GRCh38, la
normalisation depuis `Annotations` plutôt que depuis la liste SNPQC, les quatre
relations alléliques, les exclusions ambiguës, les empreintes, le refus des
dossiers d'exports individuels, le filtrage du VCF et le schéma du manifest.
Aucun fichier ACPA réel n'est lu.

Les tests D11 utilisent exclusivement quatre donneurs synthétiques. Ils
contrôlent le schéma du manifest, le déterminisme de la graine, la bijection
des copies, l'absence de paire issue du même donneur, la conservation exacte du
LD haplotypique, la concordance des génotypes, l'interpolation de carte et la
matérialisation atomique avec un exécuteur `bcftools` simulé.

Les tests D12 utilisent uniquement des identifiants et VCF synthétiques. Ils
contrôlent le schéma de configuration, la liaison empreintée D11a–D11b, le
déterminisme, la bijection des haplotypes, les proportions 50/50 et leur somme,
la conservation des allèles, l'univers de marqueurs D11b, les permissions
privées et l'absence de strate ou de donneur dans le panel aveugle.

Les tests D13 utilisent huit fondateurs entièrement synthétiques. Ils vérifient
les deux schémas, le déterminisme de la graine, les quatre trios du modèle de
cousins germains, les coefficients de parenté dérivés, la transmission exacte
de chaque copie et l'absence de recombinaison. Un scénario métissé contrôle que
l'ascendance du descendant suit les haplotypes réellement transmis, et non la
moyenne déclarée de ses parents. La matérialisation temporaire vérifie aussi les
empreintes D12, les permissions et l'absence de vérité dans le VCF aveugle.

Les tests D14 couvrent la carte avec positions dupliquées cohérentes,
l'interpolation réversible, les breakpoints internes, la censure des fenêtres
courtes, le déterminisme, la disjonction familiale et l'injection du variant.
La matérialisation temporaire lie trois manifests synthétiques, contrôle Mendel,
valide le schéma final et vérifie l'aveuglement ainsi que les permissions.

Les tests D14b couvrent la formule aux dix âges, l'arrondi extérieur des
bornes, les fenêtres infaisables, les cartes non monotones, la séparation de
la série principale et du défi, la densité insuffisante, le déterminisme de la
matérialisation et les deux schémas. Le reçu de rejeu est testé avec une chaîne
de manifests entièrement synthétique et refuse une empreinte D12 altérée.

Les tests D15 contrôlent les cinq longueurs, le déterminisme, les 25 sujets
disjoints, les copies positives, H2 et les variants extérieurs conservés, les
négatifs inchangés, les bornes nulles, les GT invalides et les transmissions.
Le matérialisateur est testé dans un workspace temporaire : schémas,
reproductibilité des empreintes, refus d'écrasement, altération d'entrée,
chemin réel interdit et rendu du rapport sans vérité individuelle.

Les tests D16 couvrent la copie H1→H2, H1 conservée, l'autozygotie des trois
intervalles, la coexistence ROH/IBD, le témoin sans injection, les copies
extérieures inchangées, les limites de carte, les défauts de densité et les GT
invalides. Ils refusent les modifications des négatifs/matrices IBD ou d'un
sujet de transmission, ainsi qu'un partage IBD source devenu incohérent.
Une régression dédiée vérifie les bornes et les copies du partage H2–H1
supplémentaire induit par le cas couplé.
Les matérialisations temporaires vérifient schémas, empreintes,
déterminisme, permissions, absence d'écrasement et rapport agrégé sans nom
d'échantillon.
