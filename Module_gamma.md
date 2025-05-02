# 🧬 Wiki DOCK6 – Module : Gamma

## 📌 Contexte d'utilisation

Le module **Gamma** permet d’estimer la date d’apparition d’une mutation génétique à partir des fréquences alléliques observées autour d’un locus. Dans le pipeline DOCK6, il est utilisé après génération des fichiers `.raw` et `.map` via PLINK pour :

- Calculer les fréquences alléliques de SNPs autour d’une mutation
- Produire un graphique illustrant le déséquilibre de liaison
- Estimer le nombre de générations depuis l’effet fondateur

---

## 📖 Récapitulatif du document

- **🔍 Présentation de Gamma**
  - Objectifs et principe
- **📁 Formats de fichiers supportés**
- **🔄 Fichiers utilisés dans le pipeline DOCK6**
- **⚙️ Commande utilisée (DOCK6)**
- **📊 Exemples de fichiers de sortie**
- **🔍 Interprétation des résultats**
- **🧪 Intégration dans DOCK6**
- **✅ Avantages et ⚠️ Limites**
- **💡 Bonnes pratiques**

---

## 🔍 Présentation de Gamma

**Objectif** : Estimer le temps écoulé depuis une mutation fondatrice à partir de la forme du pic de fréquence allélique.

**Principe** : Plus le pic est large, plus l'événement est récent. Un pic étroit suggère un événement ancien. Gamma modélise cette dynamique à l’aide d’un algorithme de régression sur les fréquences.

---

## 📁 Formats de fichiers supportés

Gamma utilise un fichier tabulé de type `.txt` contenant :
- Identifiants SNP
- Position génétique (en cM)
- Fréquence de l’allèle 1

---

## 🔄 Fichiers utilisés dans le pipeline DOCK6

### 📥 Fichiers d’entrée :
```
Nom : genotype_data.gamma_input.txt
Format : tabulé (.txt) généré automatiquement par le script Python
Contenu :
SNP	Position_cM	Freq_Allèle_1
rs123	2.1	0.25
rs124	2.3	0.35
rs125	2.5	0.40
```
Généré à partir de `genotype_data.raw` et `genotype_data.map` avec la fonction `generate_gamma_input()`.

### 📤 Fichiers de sortie :
- `gamma_output.txt` : estimation du nombre de générations (et intervalle de confiance)
- `frequence_allelique.png` : graphique du profil de fréquence
- `rapport_gamma.pdf` : synthèse avec le graphique + statistiques

---

## ⚙️ Commande utilisée (DOCK6)

```bash
Gamma -i genotype_data.gamma_input.txt -o gamma_output.txt
```

🔗 [Documentation officielle](https://www.statgen.org/gamma/)

---

## 📊 Exemples de fichiers de sortie

### `gamma_output.txt`
```
Mean generations: 14.3
95% Confidence interval: [11.8, 16.7]
```
- **Mean generations** : estimation centrale de la date d'apparition
- **Intervalle de confiance** : borne inférieure et supérieure (générations)

### `frequence_allelique.png`
Graphique en nuage de points :
- Axe X : position génétique (cM)
- Axe Y : fréquence de l’allèle 1
- Forme attendue : courbe symétrique centrée sur la mutation

*(insérer ici le visuel `frequence_allelique.png`)*

---

## 🔍 Interprétation des résultats

- Une **courbe étroite** suggère un événement ancien
- Une **courbe large** suggère un effet fondateur plus récent
- Un pic central autour de 0.5 est typique d’un déséquilibre de liaison

Ces éléments confirment ou complètent les analyses ROH et KING.

---

## 🧪 Intégration dans DOCK6

L’étape Gamma est appelée automatiquement dans `pipeline_dock6.py` par la fonction `run_gamma()`. Elle repose sur :
- Le script Python `generate_gamma_input()`
- Le binaire `Gamma`

Les résultats sont sauvegardés dans `output/gamma_results/` et utilisés dans :
- Le rapport PDF final
- L’interface Streamlit (onglet Gamma)

---

## ✅ Avantages
- Méthode robuste pour dater un événement fondateur
- Graphique facile à interpréter
- Complémentaire aux outils ROH/KING

## ⚠️ Limites
- Dépend fortement de la qualité du fichier `.map`
- Moins fiable si peu de SNPs ou si frq. extrêmes (<0.01 ou >0.99)

---

## 💡 Bonnes pratiques
- Filtrer les SNPs invariants avant génération du `.gamma_input.txt`
- Vérifier que la carte génétique contient bien les positions en cM
- Utiliser au moins 20 SNPs centrés autour de la mutation cible
- Inspecter visuellement la forme du pic pour détecter d’éventuelles anomalies

