# Visualisations consolidées — contrat de l'étape 18

## Portée

L'étape `18_build_visualizations` applique la méthode
`validated_current_run_figures_v1`. Elle rend des vues de résultats déjà
calculés ; elle ne recalcule aucune statistique, ne modifie aucun artefact
scientifique et ne combine jamais IBS, datation, LD, ROH ou sensibilités en un
score d'effet fondateur.

La sortie minimale comprend six domaines séparés : `POPULATION_STRUCTURE`,
`FOUNDER_IBS`, `VARIANT_AGE`, `LOCAL_LD`, `ROH` et `SENSITIVITY`. Une
ressemblance entre leurs figures ne constitue ni une preuve causale, ni une
preuve IBD, ni une validation croisée implicite.

## Entrées autorisées

| Domaine | Étape | Artefacts validés |
|---|---:|---|
| Structure populationnelle | `08` | `population_scores`, `population_eigenvalues`, `population_outliers` |
| IBS | `13` | `founder_segments`, `founder_analysis_summary` |
| Datation | `14` | `variant_age_estimates`, `variant_age_scenarios` |
| LD | `15` | `local_ld_summary` |
| ROH | `16` | `roh_cohort_summary` |
| Sensibilités | `17` | `sensitivity_comparisons`, `sensitivity_stability` |

Tous les artefacts sont résolus depuis `stage_inputs.json`. Leur chemin doit
rester dans le run courant, leur SHA-256 doit correspondre au descripteur et
leur producteur, signature et schéma doivent être ceux attendus. Aucun script
de visualisation ne parcourt `data/input/`, les exports ACPA, `data/output/` ou
un autre run. Les données natives PLINK/SHAPEIT ne sont pas des entrées de
l'étape 18.

## Contrôles et blocage local

Avant de rendre un domaine, le module valide :

1. l'identité du producteur et sa signature ;
2. l'empreinte du fichier et son schéma versionné ;
3. le confinement au run courant ;
4. les ensembles d'échantillons/variants et l'assemblage lorsqu'ils sont
   applicables ;
5. les effectifs et statuts répétés entre tables, résumés et audits ;
6. les unités : bp pour les positions, kb pour le fardeau ROH, cM pour les bras
   IBS et générations pour la datation ;
7. la pseudonymisation de tout identifiant affiché.

Une incohérence d'effectif, de signature, d'empreinte, de schéma ou d'unité
bloque uniquement la figure concernée. Aucun SVG partiel ou ancien n'est alors
réutilisé. L'index conserve un statut `BLOCKED` et un code de raison sans valeur
individuelle. Un domaine bloqué rend l'étape incomplète pour le rapport final.

## Non-évaluation et petits effectifs

`NOT_EVALUATED`, `NOT_ESTIMATED`, les cellules manquantes et les exclusions ne
sont jamais transformés en zéro. Un domaine scientifiquement non évaluable mais
intègre produit un panneau neuf qui affiche le statut, les effectifs, les
raisons contrôlées et les limites. Ce panneau est distinct d'une figure bloquée
par intégrité et d'une étape désactivée.

Les seuils et statuts producteurs sont repris tels quels. L'étape 18 ne relâche
aucun seuil et ne requalifie pas un résultat exploratoire en résultat primaire.

## Contenu minimal par domaine

- `population_structure.svg` : scree plot, projection PC1–PC2 pseudonymisée,
  référence indépendante, apparentés projetés et outliers exploratoires
  distingués ; effectifs et coordonnées manquantes visibles, sans attribution
  automatique d'ascendance ni exclusion automatique.
- `founder_ibs.svg` : bras gauche/droit en cM, unités pseudonymisées, statut des
  unités, exclusions, confiance absente, effectif et mention « IBS, pas preuve
  IBD ».
- `variant_age.svg` : modèle corrélé primaire, modèle indépendant et scénarios
  séparés, intervalles en générations, petits effectifs et hypothèses Gamma.
- `local_ld.svg` : `r²` génotypique et `D′` haplotypique séparés, cohortes et
  classes de distance distinctes, effectifs et paires évaluables.
- `roh.svg` : fardeau genome-wide distinct du chromosome cible, effectifs
  observés/évalués et rappel que l'autozygotie individuelle n'est pas l'IBS/IBD
  inter-individus.
- `sensitivity.svg` : primaire explicitement étiqueté et scénarios
  exploratoires distincts pour chaque domaine, non-évaluations visibles, aucun
  score composite.

## Provenance et classification

Chaque SVG possède un fichier `.figure.json` qui enregistre le run, le domaine,
le script et sa version, les artefacts sources et leurs SHA-256, les effectifs,
les statuts, les valeurs manquantes, les exclusions, les limites, la légende et
la pseudonymisation. Les figures et leurs provenances sont classées
`sensitive_genetic` par défaut.

`figure_index.json` référence les cinq domaines sans dupliquer de données
individuelles. `visualization_completeness.json` compte les figures rendues,
non évaluées et bloquées et indique si l'étape 19 peut assembler un rapport
complet. Ces documents n'autorisent aucune lecture d'un autre run.
