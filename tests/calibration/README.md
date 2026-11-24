# Tests de calibration unifiée R07 v2–R14

Cette suite couvre exclusivement les contrats de `panel_master_generic_v1` :
schémas, garde d'aveuglement, disjonctions, routeur de vérité, ordre de décision
R08 et reçu statique de préparation. Les tests R10 couvrent aussi la topologie
fermée, les dimensions synthétiques, la séparation des graines et références,
les formats PLINK/VCF, la reproductibilité binaire, l'inventaire exhaustif et
les règles d'isolation. Ils n'appellent aucun outil scientifique externe.

La révision post-R14 de R07 possède en plus un validateur déterministe. Ses
tests vérifient la conservation des minima scientifiques, le calendrier exact
des appels SHAPEIT5, la fermeture des mondes IBD, les allocations de ressources
et le refus d'une réduction silencieuse de l'effectif primaire. Ils ne lancent
aucune sonde R11b ni campagne scientifique.

Les tests R11b couvrent le dessin fermé des trois sondes SHAPEIT5 et des 18
sondes IBD, les bornes prédictives et le statut en lecture seule. Les appels
externes réels restent exclusivement dans le workspace de ressources.

Les tests de la campagne D01–D42 ont été supprimés à R09 : leurs schémas,
configurations et protocoles ont été définitivement retirés à R01 et ne sont
pas des dépendances compatibles avec la release unifiée.
