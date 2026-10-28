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
