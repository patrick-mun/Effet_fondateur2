# Demandes d'amélioration — interface SaaS du pipeline V2

Ce fichier consigne les demandes exprimées lors de l'échange du 20 août 2026
sur l'évolution du pipeline V2 vers un usage SaaS. Il complète
`PIPELINE_V2_PRECODE.md` et `SESSION.md` sans les remplacer : il capture des
intentions produit/architecture encore à trancher, pas des contrats d'étape
validés.

## 1. Interface utilisateur (maquettes)

Une première direction visuelle a été validée comme point de départ :
- Page d'accueil publique (vitrine, présentation du logiciel).
- Tableau de bord (études, runs récents, statuts).
- Préparation des données (import des sources, table maître des échantillons,
  validation avant run).
- Assistant de lancement d'un run (sélection étude/profil, vérification,
  lancement).
- Suivi d'un run étape par étape (groupes d'audit `00–04`, `05–09`, `10–14`,
  `15–19`).
- Rapport de résultats (PCA d'ascendance, ROH, parenté, datation).

Maquette publiée (Artifact Claude) : à retrouver via `/artifacts` dans la
session ayant servi à la générer. Statique, direction éditoriale/clinique
(Newsreader + Public Sans + IBM Plex Mono, accent sarcelle) — non figée,
modifiable.

**Reste à faire** : valider la navigation entre écrans, décider si un
prototype cliquable est nécessaire avant développement.

## 2. Support de plusieurs types de puces en entrée

**Constat actuel** : le seul convertisseur d'entrée existant
(`03_convert_acpa`) est spécifique aux exports ACPA/ChAS (Affymetrix
ChromosomeAnalysisSuite, `Forward Strand Base Calls`). Rien d'autre n'est
branché aujourd'hui.

**Ce qui joue en faveur du changement** : à partir de l'étape `04`, le
pipeline identifie chaque variant par assemblage/chromosome/position/REF/ALT
et non par identifiant de sonde propriétaire — l'architecture aval est déjà
agnostique de la puce.

**Demande** : pouvoir constituer les fichiers d'entrée du pipeline à partir
de formats d'export différents (Illumina GenomeStudio, autres panels
Affymetrix, VCF de séquençage direct, etc.), pas seulement ACPA/ChAS.

**Ce que ça implique concrètement** :
- Un convertisseur dédié par format source, produisant la même sortie
  normalisée que `03_convert_acpa` (mêmes triplets PLINK, mêmes tables
  d'audit `sample_alignment.tsv` / `*_variant_audit.tsv`).
- Chaque convertisseur doit gérer les pièges propres à son format :
  orientation de brin, build génomique natif de la puce (liftover éventuel
  vers GRCh38), mapping sonde → coordonnée physique.
- Le contrat de sortie commun (`schemas/plink_dataset.schema.json` et
  assimilés) ne change pas — seul l'amont est à multiplier.

**Priorité et formats à couvrir** : à préciser (quelle puce/format est le
plus urgent après ACPA ?).

## 3. Architecture cloud / SaaS — confidentialité des données

### Proposition initiale de l'utilisateur

- Aucune donnée d'entrée ou de sortie stockée en base de données.
- La base de données ne conserve que le nom du run (et un pointeur vers le
  dossier local de l'utilisateur ayant les droits de connexion).
- Le logiciel s'exécute dans le cloud.
- Les fichiers produits à chaque étape sont rapatriés en local dès leur
  création, et placés en RAM serveur si le run en cours en a besoin pour
  l'étape suivante.

### Lecture de cette proposition

**L'intention est la bonne et cohérente avec les décisions déjà actées dans
le dépôt** (`SESSION.md` : *« ne jamais envoyer les données privées de
l'étude vers un service externe »*, données individuelles classées
sensibles). Une base de données qui ne contient que des métadonnées de run
(nom, statut, horodatage, empreintes — pas les génotypes) est un bon
principe, et correspond déjà à ce que le manifest de run V2
(`schemas/run_manifest.schema.json`) fait conceptuellement.

**Deux points de la mise en œuvre littérale posent un problème pratique,
indépendamment de l'intention** :

1. **« Tout en RAM »** est irréaliste au sens strict. Le pipeline manipule
   des fichiers volumineux (VCF/BCF phasés, extraits 1000G, panels de
   parenté) sur ~20 étapes séquentielles, et les outils externes (PLINK,
   KING, SHAPEIT5, bcftools, R) lisent et écrivent des fichiers réels sur
   disque — ils ne savent pas travailler sur un objet mémoire abstrait. Un
   run réel peut nécessiter des dizaines de Go cumulés ; les garder tous en
   RAM pour toute la durée d'un run de plusieurs heures est coûteux et
   fragile (perte totale en cas de redémarrage du worker).
   → **Raffinement réaliste** : monter le répertoire de travail du run sur
   un `tmpfs` (système de fichiers en RAM) le temps du run. Les outils
   externes voient de vrais fichiers, mais rien ne touche jamais un disque
   physique, et tout disparaît à l'arrêt du conteneur. Ça conserve l'esprit
   de la demande sans casser la compatibilité avec les outils externes.

2. **« Rapatrié en local dès la création de chaque fichier »** (à chaque
   étape, pas seulement à la fin) crée une dépendance forte à la connexion
   de la machine locale de l'utilisateur pendant tout le run : si elle est
   déconnectée, ou si une étape ultérieure a besoin de relire un artefact
   d'une étape bien antérieure, il faut soit le garder aussi côté serveur,
   soit le retélécharger depuis le poste local — ce qui double les
   transferts et rend le pipeline dépendant d'un aller-retour réseau à
   chaque étape (20 étapes × fichiers parfois volumineux).
   → **Raffinement réaliste** : stockage éphémère isolé par run côté cloud
   (bucket ou volume dédié, chiffré, accessible uniquement au run et à son
   propriétaire), qui vit uniquement pendant l'exécution et est purgé
   immédiatement à la fin (succès ou échec). Le rapatriement vers
   l'utilisateur se fait une fois, à la fin du run (ou à la demande), pas à
   chaque étape.

**Point non couvert par la proposition, à trancher** : des données
génétiques individuelles hébergées dans le cloud, même de façon éphémère,
peuvent relever de régimes réglementaires spécifiques (données de santé/
génétiques au sens RGPD, hébergement de données de santé en France si le
contexte médical l'exige). À vérifier avant tout choix d'hébergeur, même
pour du transitoire.

### Recommandation de synthèse

- Base de données : métadonnées de run uniquement (nom, statut par étape,
  empreintes, chemin de sortie). Jamais de génotype, jamais de contenu.
- Répertoire de travail par run : `tmpfs` ou volume chiffré isolé par
  tenant/run, jamais partagé entre études.
- Purge automatique et systématique du répertoire de travail à la fin du
  run (succès ou échec), sans action manuelle requise.
- Rapatriement des résultats vers l'utilisateur en une fois à la fin du run
  (ou export à la demande), pas fichier par fichier à chaque étape.
- Vérifier le cadre réglementaire applicable (RGPD données sensibles,
  hébergement de données de santé) avant choix d'infrastructure.

**Reste à faire** : choisir l'infrastructure cible (fournisseur, isolation
par tenant, dimensionnement du `tmpfs`/volume selon la taille des runs
réels), et confirmer le cadre réglementaire applicable.

## 4. Estimation de coût — hébergement, installation, fonctionnement

Ordre de grandeur, pas un devis : dépend surtout de la fréquence réelle des
runs et de l'usage visé (labo interne vs SaaS multi-clients). Découpé en
postes de nature différente.

### Partie toujours active (app web + API + base de métadonnées)

Légère par construction : aucune donnée génétique n'y transite, seulement
l'interface et les statuts de run. Un petit serveur géré (Scaleway, OVH,
Hetzner) suffit.

**~15–50 €/mois**

### Calcul à la demande (le run lui-même) — poste dominant et variable

Aucune machine allumée en permanence : un worker éphémère se lève pour la
durée du run puis disparaît. Un run avec phasage SHAPEIT5 + KING + PCA sur
une instance correcte (8–16 vCPU, 32–64 Go RAM pour le `tmpfs`) coûte de
l'ordre de 1 à 6 € par run (quelques heures de calcul facturées à l'usage).
Pour un usage recherche ponctuel (quelques runs/semaine) le total mensuel
suit linéairement le nombre de runs, pas le nombre d'utilisateurs inactifs.

**~20–150 €/mois** pour un usage recherche modéré (quelques runs/semaine)

### Stockage

Quasi nul si la purge automatique fonctionne comme prévu (le stockage
éphémère par run n'existe que le temps du run). Prévoir un petit espace
tampon pour que l'utilisateur récupère son résultat final avant purge.

**Quelques €/mois**

### Outils externes (PLINK, KING, SHAPEIT5, bcftools, R)

Open source, aucun coût de licence.

**0 €**

### Facteur pouvant multiplier significativement la facture : conformité HDS

Si l'hébergement doit être certifié HDS (hébergeur de données de santé,
potentiellement requis en France pour des données génétiques individuelles
selon le contexte), les offres certifiées coûtent nettement plus cher que
l'offre standard. Si le projet devient lui-même responsable de traitement
plutôt que de s'appuyer sur un hébergeur déjà certifié, l'audit de
conformité peut représenter plusieurs milliers d'euros. **À vérifier avant
tout chiffrage définitif** — c'est le facteur qui peut multiplier le budget
par un facteur important, indépendamment du volume d'usage.

### Non inclus dans cette estimation

Le temps de développement pour construire l'API et brancher l'orchestrateur
V2 derrière (`effet-fondateur run/resume`) n'est pas un coût d'hébergement :
c'est un travail d'ingénierie one-shot, à chiffrer séparément.

### Total, ordre de grandeur

Usage labo/recherche modéré, hors HDS et hors développement :

**~50–250 €/mois**

**Reste à faire** : préciser le scénario réel (nombre de runs/mois, usage
interne ou multi-clients) pour affiner la fourchette, et confirmer si un
hébergement certifié HDS est requis.

## 5. Monétisation — internationalisation, facturation au run publiable

**Piste évoquée** : facturer ~1500 $ par run finalisé/publiable (résultat
utilisable dans une publication scientifique), plutôt qu'à l'exécution brute
(cf. §3 sur l'itération normale d'un run).

**Ordre de grandeur du prix** : pas déraisonnable en soi pour une analyse de
cohorte complète (QC, parenté, phasage, datation, ascendance, rapport
publiable) — probablement inférieur au coût d'un bio-informaticien faisant
l'équivalent manuellement pendant plusieurs semaines, et cohérent avec les
tarifs de services d'analyse génomique en labo/CRO. À moduler selon la
taille de cohorte plutôt qu'un tarif plat unique (voir §3).

### Risque de contournement identifié

Facturer sur une action déclenchée côté client (« bouton valider le rapport
final ») est contournable : les artefacts intermédiaires (`*.parquet`,
`*.json` — segments founder, discordances, panels de parenté, etc.) sont
déjà accessibles en téléchargement à chaque étape dans le suivi de run,
avant même la génération du rapport final. Ce sont ces fichiers, pas le
HTML, qui portent la valeur scientifique. Un utilisateur peut donc les
récupérer sans jamais déclencher l'action de paiement si le verrou est
uniquement dans l'interface.

**Principe à respecter** : le contrôle de paiement doit être **côté
serveur, sur l'accès aux données**, jamais sur une action côté client.

- Le suivi de progression (statuts, étapes, logs) peut rester visible
  gratuitement — n'expose aucun résultat scientifique.
- Chaque endpoint de téléchargement/export d'artefact vérifie le paiement à
  chaque appel, indépendamment de l'état de l'interface (pas de drapeau
  « payé » stocké côté client).
- Aperçu dégradé possible avant paiement (valeurs tronquées/floutées,
  watermark), export haute précision débloqué uniquement après paiement
  confirmé côté serveur.

### Risque structurel additionnel — pipeline ouvert

Le pipeline V2 (`src/effet_fondateur/`) reste un logiciel installable et
exécutable localement dans ce dépôt. À l'échelle internationale, un
utilisateur technique peut l'installer chez lui sans jamais passer par la
version hébergée payante. La monétisation ne peut donc pas reposer sur « le
code est payant », mais sur la valeur de l'infrastructure gérée : cache de
références déjà peuplé, calcul à la demande sans installation d'outils
externes (PLINK/KING/SHAPEIT5/bcftools/R), support, interface. Modèle
viable (courant en open source commercial), mais à assumer dès la
conception plutôt qu'à découvrir après coup.

**Reste à faire** : concevoir le modèle d'autorisation serveur (paiement →
déblocage export) avant tout développement de la facturation, et décider du
niveau d'aperçu gratuit acceptable.
