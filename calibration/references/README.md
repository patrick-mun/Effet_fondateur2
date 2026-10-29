# Références publiques de calibration

`public_reference_catalog.json` est la liste fermée D09 des sources publiques
autorisées pour préparer les futurs panels synthétiques. Ce catalogue ne
contient ni génotype, ni identifiant d'étude, ni donnée de `dock6_real`.

Chaque entrée fixe un identifiant, une version, une origine HTTPS, la portée de
l'assemblage, les conditions d'usage revues et une empreinte SHA-256. Elle est
également liée par SHA-256 au catalogue technique existant sous
`config/references/`, afin qu'une modification de ce dernier invalide le
contrôle plutôt que de changer silencieusement la référence.

## Politique d'acquisition

- D09 ne télécharge aucun VCF, BCF ni panel génétique.
- Une ressource téléchargée ultérieurement doit être contrôlée contre
  l'empreinte publiée lorsqu'elle existe, puis recevoir un SHA-256 local avant
  toute utilisation.
- Pour le panel phasé 1kGP, IGSR publie des MD5 par chromosome. Le SHA-256 du
  manifeste officiel fige donc la version distante ; chaque VCF et index acquis
  devra en plus être vérifié par MD5 puis enregistré par SHA-256 dans le
  manifest D07 du run.
- Les fichiers acquis et leurs empreintes d'exécution restent sous
  `calibration_workspace/` et ne sont jamais ajoutés à Git.

## Conditions d'usage

Les données 1000 Genomes sont utilisables sans embargo sous les conditions
IGSR et EMBL-EBI, avec vérification des déclarations de réutilisation et des
citations associées. Le catalogue interdit leur redistribution par le projet.

Les cartes GRCh38 sont épinglées sur un commit SHAPEIT4 et documentées par les
auteurs comme provenant de HapMap. La licence MIT du dépôt SHAPEIT4 est
vérifiée, mais elle ne suffit pas à établir sans réserve tous les droits sur les
données dérivées. Leur usage est donc limité localement à la calibration et
leur redistribution est interdite sans revue complémentaire.

D10 définit désormais la séparation déterministe par groupe d'indépendance
entre donneurs de calibration, donneurs de validation et référence de phasage
commune. Une première matérialisation vérifiée a été produite à partir du cache
public sous `calibration_workspace/d10_1kg_3202_20220422_v1/` ; ses listes
restent ignorées par Git.
L'autorisation D09 d'une source ne constitue pas à elle seule cette séparation.

`CYTOSCAN_750K_ACCEL_MARKERS.md` fixe la source et le format D11b de l'univers
observable. Il distingue explicitement la conception Accel `NA36.r1` de
l'ancienne CytoScan 750K `NA36.r7`, interdit de dériver la liste de marqueurs
depuis un échantillon et exige les SHA-256 du paquet constructeur et du TSV
normalisé. L'entrée Thermo Fisher est inscrite au catalogue D09 depuis
l'acquisition et la revue du 29 octobre 2026. D11b lit le ZIP conservé hors du
dépôt, produit directement le TSV sous `calibration_workspace/` et refuse une
empreinte, une base SQLite ou des métadonnées incompatibles. Ces deux charges
utiles restent locales et hors Git.
