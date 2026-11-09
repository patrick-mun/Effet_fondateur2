# D23a — Préparation du phasage, version 1.0.0

Décision CAL-DEC-007 : préparation approuvée explicitement par l'utilisateur
le 30 octobre 2026. Ce document est fixé avant la génération. Aucun seuil
scientifique de performance n'est fixé ou ajusté ici.

Révision de clôture avant tout phasage : la configuration active est
`d23a_phasing_panel_v2.json`, avec des graines distinctes 2301 (fondateurs),
2302 (méioses privées) et 2303 (futur phasage). Le brouillon v1 reste conservé,
mais n'est pas utilisable : il proposait la graine des méioses dans les
paramètres aveugles. La révision autorise uniquement l'identifiant de
configuration, la graine de phasage et les chemins du nouveau run ; tous les
artefacts génétiques doivent conserver leurs SHA-256 et être relus. Aucun
paramètre de simulation ou seuil n'est adapté après observation.

## Domaine et sources

Nouveau banc de calibration uniquement, sans modifier D10–D22 ni ouvrir les
génotypes de validation. Les 22 autosomes utilisent les VCF 1000G du Bureau,
leurs MD5 officiels et leurs SHA-256 locaux. L'annotation vient uniquement du
paquet officiel CytoScan 750K Accel NA36.r1, jamais d'un échantillon, du cache
interne ou de données d'étude. Cartes HapMap/SHAPEIT GRCh38 issues de l'archive
épinglée D09, avec interpolation sans extrapolation.

Le panel borné contient huit trios et seize fondateurs supplémentaires, soit
40 sujets et 32 fondateurs issus de 64 donneurs D10 calibration. Un seul
donneur est pris par groupe d'indépendance D10, et chaque donneur contribue
une seule copie haplotypique. Les deux copies d'un fondateur proviennent de
donneurs de même population source, distincts. Les familles synthétiques ne
prétendent pas représenter la généalogie des donneurs ou l'indépendance KING.
Les unités prévues pour l'incertitude sont 8 familles et 16 sujets isolés,
pas 40 réplications indépendantes, ni les chromosomes comme sujets.

Sélection avant inspection des performances : jusqu'à 4 096 candidats répartis
en rangs physiques par autosome, puis jusqu'à 1 024 SNP harmonisés. S'ajoutent
les marqueurs de puce harmonisés dans la fenêtre cible ±5 Mb sur chr5.
Doublons de position, orientations ambiguës, discordances alléliques, GT de
donneur incomplets et positions hors carte sont exclus et comptés. Minimum
256 marqueurs par autosome et 10 000 au total. C'est un sous-échantillon
autosomal de la puce, pas une reproduction de sa densité complète.

Les transmissions suivent un modèle de Haldane sans interférence :
N ~ Poisson((cM_fin − cM_début)/100) ; les points de recombinaison conditionnels
sont uniformes en distance génétique, la copie initiale est tirée à 0/1.
Les sous-graines SHA-256 dépendent de la graine, du chromosome et du rôle.
La phase du descendant est père|mère, conservée uniquement dans la vérité.
Ce modèle simple ne qualifie ni cartes sexuées ni interférence des crossing-over.

## Cible et métadonnées

Réutiliser l'annotation technique publique D19 épinglée (APC chr5), pas DOCK6,
et jamais un statut clinique. La règle explicite injecte 0|0, 0|1, 1|1 selon
le rang du fondateur modulo trois ; l'enfant reçoit les allèles des copies
parentales réellement transmises. Le locus est conservé même s'il est absent
de la puce : c'est alors une observation synthétique supplémentaire déclarée,
pas une sonde de puce prétendue. Le génotype moléculaire non phasé est exporté
avec une provenance de simulation ; clinique UNKNOWN, groupe UNSPECIFIED.
FID/PID/MID décrivent seulement le pedigree synthétique construit.

La vérité phasée et les correspondances donneurs restent privées. Les entrées
du pipeline contiennent les observations non phasées, la filiation déclarée,
les génotypes moléculaires explicites et la seule liste de références D10.
Aucun haplotype vrai, point de recombinaison ou correspondance donneur n'y figure.

## Raccord de référence

Nouvelle entrée optionnelle `inputs.phasing_reference_selection`, intégrée
aux signatures et audits de 12/12A. Elle contient exactement la référence D10
(1 403 sujets), son ordre canonique, ses empreintes et l'engagement du manifest
D10. La préparation vérifie les trois listes, les groupes et leur disjonction.
L'option `require_reference_sample_selection: true` empêche tout retour
silencieux au panel complet. Une sélection inconnue, modifiée, dupliquée,
incomplète dans le VCF ou non liée à la configuration est bloquante.

12 filtre lors de l'extraction, avant harmonisation. 12A filtre avant le cache
d'extraction dont la clé inclut la liste de sujets. Les sorties sont relues
pour exiger les sujets et l'ordre exacts. Les références publiques disponibles
ne deviennent jamais implicitement la référence autorisée de calibration.
Sans cette nouvelle entrée, le comportement des profils historiques reste
inchangé. Aucun profil d'étude n'est modifié.

## Contrôles et limites

Contrôler copies fondatrices, transmissions, dosages exportés, carte, présence
unique de la cible, pedigrees, génotypes moléculaires et absence d'intersection
D10. Inventorier et empreinter toutes les sorties avant publication atomique.
Refuser tout écrasement ; ne pas écrire dans les données d'étude.

D23a ne lance ni SHAPEIT5 ni pipeline complet, ne résout pas 07 et ne qualifie
aucune performance. Les quatre évaluateurs et rapports 10A/11/12/12A restent
à développer dans D23. Leur protocole séparera commutation, erreur ponctuelle
et orientation de cible modulo permutation globale H1/H2, avec ancrage
indépendant de la cible et strates fréquence/densité/manque. Le rapport de
préparation comporte graphiques de comptes et formules du modèle, pas de
sensibilité, taux de phase ou intervalle de performance inventés.

`calibrated=false`, `founder_effect_proven=false`. Le panel de validation
reste fermé ; ce banc ne remplace pas le pilote D33.
