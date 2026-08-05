# Suivi de session

Dernière mise à jour : 5 août 2026

## État du dépôt

- Dépôt GitHub : `patrick-mun/Effet_fondateur2`.
- Branche active : `main`, synchronisée avec `origin/main` avant cette mise à
  jour.
- PR `#1` fusionnée dans `main` le 4 août 2026.
- Commit de fusion : `f2fec8f`.
- Environnement : Python 3.12.13 dans `.venv`.
- Outils disponibles : PLINK 1.9, KING, bcftools, Rscript, Gamma et packages R
  du projet.

## Bilan des avancées

- Inventaire du dépôt réalisé et artefacts manifestement obsolètes supprimés.
- Environnement Python reproductible installé et documenté dans `README.md`.
- Instructions de codage, tests, maintenance et sécurité ajoutées à
  `AGENTS.md`.
- Convertisseur multi-échantillons `acpa_to_plink.py` créé :
  - lecture des exports ACPA ;
  - sélection du chromosome ;
  - annotation dbSNP hg38 ;
  - génération de PED, MAP, groupes, cas, témoins et rapports QC ;
  - refus d'écrasement sans `--force`.
- Protocole ACPA complet et préparation de `samples.tsv` documentés.
- Simulateur `simulate_unrelated_controls.py` créé pour générer des témoins
  techniques reproductibles sous HWE, avec graine, audit et limites explicites.
- Injecteur `inject_mutation.py` créé :
  - insertion triée d'un marqueur `rsMUT...` dans une copie PED/MAP ;
  - métadonnées du gène, variant, HGVS c./p., transcrit, position et allèles ;
  - génotypes obligatoirement affectés par règle explicite de groupe ou par
    individu ;
  - audit TSV, rapport JSON, empreintes SHA-256 et données pour le rapport final.
- README mis à jour avec les commandes de conversion, simulation, injection et
  validation PLINK/Mendel.

## Jeux de données actuellement disponibles

### Conversion ACPA réelle

- Dossier conseillé pour la dernière relance : `data/output/acpa_chr19_relance/`.
- 9 individus réels.
- 2 725 marqueurs du chromosome 19.
- Aucun rsID non résolu dans le rapport de conversion.
- PED/MAP et fichiers `groupes.txt`, `cas.txt`, `temoins.txt` présents.

### Jeu de test avec témoins simulés

- Dossier : `data/output/acpa_chr19_avec_temoins_simules/`.
- 9 individus ACPA et 50 témoins synthétiques, soit 59 individus.
- 2 725 marqueurs avant injection de la mutation.
- Groupes : 5 `ATTEINT`, 3 `HTZ`, 1 `SAINS`, 50 `TEMOIN_SIMULE`.
- Graine de simulation : `20260804`.
- Réservé aux tests techniques et exploratoires, sans interprétation biologique
  du LD, des haplotypes ou des ROH.

### Entrée actuellement lue par le pipeline

- `data/input/complex_simulation/genotype_data.ped` contient 127 individus.
- `data/input/complex_simulation/genotype_data.map` contient 2 726 marqueurs,
  dont `rsMUT_11237780`.
- Ce jeu est historique et ne correspond pas au nouveau jeu ACPA de 9 individus
  avec 50 témoins simulés.
- `cas.txt` et `temoins.txt` sont absents de ce dossier, donc le pipeline actuel
  ne doit pas encore être lancé sur cette entrée.
- Aucune injection définitive n'a été enregistrée dans
  `data/output/acpa_chr19_mutation/` ; la validation de l'injecteur a été faite
  uniquement dans un dossier temporaire.

## Validations effectuées

- `pip check` : aucune dépendance cassée.
- Imports Python principaux et `run_pipeline` réussis.
- Packages R `adegenet`, `ggplot2`, `ape`, `ade4` et `poppr` disponibles.
- Conversion ACPA réelle validée par PLINK : 9 individus, 2 725 marqueurs et
  taux de génotypage de 98,04 %.
- Jeu avec témoins simulés validé par PLINK et KING : 59 individus et 2 725
  marqueurs.
- Injection technique temporaire validée par PLINK : 59 individus, 2 726
  marqueurs et taux de génotypage de 98,25 %.
- Contrôle mendélien temporaire : 10 erreurs, dont une au marqueur injecté pour
  `F2/E_82303707` avec les règles de groupe d'exemple.
- Tests ciblés de préparation des données : 16 réussis.
- Suite Python globale : 28 réussis et 4 échecs historiques.

## Problèmes connus

1. Les allèles génomiques REF/ALT de la mutation DOCK6 doivent être confirmés
   sur le brin de référence GRCh38 ; la notation HGVS du transcrit ne suffit pas.
2. Les génotypes individuels réels de la mutation restent à renseigner. Les
   règles de groupe d'exemple créent une incompatibilité mendélienne pour
   `F2/E_82303707` (`*/* x T/T -> G/G`).
3. `run_pipeline.py` utilise encore des chemins fixes dans
   `data/input/complex_simulation/` et ouvre automatiquement le rapport HTML.
4. L'interface Streamlit écrit `user_input.ped/map`, tandis que le pipeline lit
   `genotype_data.ped/map`.
5. La datation Gamma sépare actuellement les positions en deux moitiés au lieu
   d'utiliser des longueurs gauche/droite pertinentes autour de la mutation.
6. Quatre tests historiques échouent encore :
   - confirmation interactive non neutralisée dans un test Gamma ;
   - argument `n` absent d'un test d'estimation d'âge ;
   - deux mocks de prétraitement ne créent pas `filtered_data.ped`.
7. Le simulateur familial historique et les anciens convertisseurs sont encore
   conservés tant que leur dernière logique utile n'a pas été vérifiée.

## Priorités de la prochaine session

1. Intégrer les vrais témoins ACPA lorsqu'ils sont disponibles, mettre à jour
   `samples.tsv`, puis reconstruire un jeu réel homogène du chromosome 19.
2. Confirmer pour la mutation DOCK6 : assemblage, position, transcrit, allèles
   génomiques REF/ALT, HGVS c., HGVS p. si connu et génotype de chaque individu.
3. Lancer l'injecteur sur le jeu définitif, examiner
   `mutation_genotype_audit.tsv`, puis résoudre les erreurs mendéliennes.
4. Valider le PED/MAP injecté avec PLINK avant de remplacer de façon contrôlée
   l'entrée historique de `data/input/complex_simulation/`.
5. Raccorder les chemins du pipeline et de Streamlit au jeu réellement choisi,
   sans copie manuelle ambiguë.
6. Corriger les quatre tests historiques, puis revoir les entrées scientifiques
   de Gamma.
7. Exécuter le pipeline complet sur une copie contrôlée et examiner le rapport
   avant toute interprétation biologique.

## Décisions à conserver

- Ne jamais déduire silencieusement un génotype de mutation du statut clinique.
- Utiliser les témoins synthétiques uniquement pour tester le fonctionnement.
- Ne pas écraser `data/input/` ou `data/output/` sans demande explicite et
  `--force` lorsque l'outil le prévoit.
- Considérer `data/output/` comme reproductible, mais pas comme source primaire.
- Vérifier PLINK, l'audit de mutation et les erreurs mendéliennes avant le
  pipeline complet.
- Ne pas lancer le pipeline complet sans prévenir : il écrit de nombreuses
  sorties et ouvre le navigateur.

## Reprise rapide

Au début de la prochaine session :

```bash
cd /Users/utilisateur/Documents/python_programme/Effet_fondateur2
source .venv/bin/activate
git status -sb
```

Lire ensuite `AGENTS.md`, ce fichier et la section de protocole correspondante
dans `README.md`. La première décision attendue concerne l'utilisation des vrais
témoins ACPA et la confirmation moléculaire de la mutation.
