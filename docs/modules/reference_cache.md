# Étape 12.2 — Cache immuable de référence

## Responsabilité

`cache_reference_panel` télécharge ou réutilise les quatre fichiers publics
résolus en `12.1` : VCF phasé, index tabix, README et manifeste fournisseur. Le
cache est partagé entre les runs mais chaque entrée est liée à un catalogue, une
release, un assemblage, un chromosome et des empreintes exactes.

Le chemin configuré par défaut est `data/cache/references`. Ce dossier est exclu
de Git et ne doit jamais contenir de données privées.

## Identité d'une entrée

La clé SHA-256 du cache dépend de :

- l'empreinte du catalogue ;
- l'identifiant du panel, du fournisseur et de la release ;
- l'assemblage et le chromosome ;
- chaque nom, URL et empreinte attendue.

Deux releases ou deux catalogues différents ne partagent donc jamais
silencieusement la même entrée.

## Téléchargement et publication

Un verrou `flock` exclusif sérialise les demandes portant sur la même entrée.
Les fichiers sont écrits dans un dossier temporaire du même système de fichiers,
puis contrôlés avant une unique publication par renommage atomique.

Les contrôles sont :

- MD5 officiel pour le VCF et l'index ;
- SHA-256 épinglé pour le README et le manifeste fournisseur ;
- SHA-256 local et taille pour les quatre fichiers ;
- schéma et identité du manifeste de cache ;
- absence de liens symboliques et noms de chemin sûrs.

Après publication, les fichiers passent en `0444` et le dossier en `0555`.
Chaque réutilisation recalcule les empreintes et vérifie les permissions. Une
entrée publiée corrompue ou redevenue inscriptible bloque : elle n'est jamais
remplacée automatiquement.

## Reprise et mode hors ligne

Après acquisition du verrou, une reprise supprime uniquement les dossiers
temporaires portant la clé attendue et laissés par une interruption. Elle ne
supprime jamais une entrée publiée.

Avec `reference_cache_offline: true`, une entrée valide peut être réutilisée,
mais une absence produit `reference_cache_offline_miss` sans appel réseau.

## Paramètres

```yaml
phase_target_region:
  parameters:
    reference_cache_dir: data/cache/references
    reference_cache_offline: false
    reference_lock_timeout_seconds: 60
    reference_download_timeout_seconds: 7200
    reference_download_chunk_size: 1048576
```

## Limite de validation

Les tests utilisent des fichiers synthétiques minuscules et simulent le
téléchargement externe. Aucun VCF 1000 Genomes complet n'est téléchargé pendant
`12.2`. Le premier téléchargement réel sera volontaire, volumineux et devra être
annoncé avant exécution.
