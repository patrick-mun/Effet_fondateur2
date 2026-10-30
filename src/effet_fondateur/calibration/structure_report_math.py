"""Rappels mathématiques des méthodes D22, sans calcul ni dépendance réseau."""

from __future__ import annotations

from html import escape


def _fraction(numerator: str, denominator: str) -> str:
    return f'<span class="fraction"><span>{numerator}</span><span>{denominator}</span></span>'


def _equation(number: int, title: str, formula: str, explanation: str) -> str:
    """Formule typographiée et explication indivisibles dans le PDF."""
    return f'<aside class="equation" data-equation="{number}"><p class="equation-title">Équation {number} — {escape(title)}</p><div class="math" role="math" aria-label="{escape(title, quote=True)}">{formula}</div><p>{explanation}</p></aside>'


def pca_equations() -> str:
    formula = '<div>p<sub>j</sub> = ' + _fraction('Σ<sub>i ∈ R appelé</sub> g<sub>ij</sub>', '2 n<sub>j</sub>')
    formula += ' ; z<sub>ij</sub> = ' + _fraction('g<sub>ij</sub> − 2p<sub>j</sub>', '√[2p<sub>j</sub>(1 − p<sub>j</sub>)]') + '</div>'
    formula += '<div>X = Z / √m ; X<sub>R</sub> = UΣV<sup>T</sup> ; S = XV</div>'
    formula += '<div>λ<sub>k</sub> = σ<sub>k</sub><sup>2</sup> / (n<sub>R</sub> − 1)</div>'
    first = _equation(1, "Standardiser, ajuster et projeter la PCA", formula,
                       'g<sub>ij</sub> : dosage A1 (0, 1 ou 2) ; R : référence ; n<sub>j</sub> : sujets appelés dans R au locus j ; '
                       'm : loci informatifs ; n<sub>R</sub> : taille de R. Un manque est remplacé par 2p<sub>j</sub>, donc z=0. '
                       'V porte les axes, S les scores sans unité, σ les valeurs singulières et λ les valeurs propres. '
                       'Fréquences et axes sont ajustés sur R seulement : les projetés ne doivent pas déplacer R. '
                       'Source : contrat 08 et module analyze_population_structure.')
    formula = '<div>D<sub>i</sub><sup>2</sup> = Σ<sub>k=1</sub><sup>2</sup> ' + _fraction('S<sub>ik</sub><sup>2</sup>', 'λ<sub>k</sub>')
    formula += ' ; p<sub>i</sub> = exp(−D<sub>i</sub><sup>2</sup>/2)</div>'
    formula += '<div>Alerte ⇔ p<sub>i</sub> &lt; α ⇔ D<sub>i</sub><sup>2</sup> &gt; −2 ln(α)</div>'
    return first + _equation(2, "Distance et alerte exploratoire", formula,
                             'Deux axes sont utilisés ici ; α=0,01. La queue exp(−D²/2) est celle du χ² à deux degrés de liberté. '
                             'Son utilisation pour ces scores est une approximation exploratoire, pas une p-value biologiquement calibrée. '
                             'La figure des distances montre les sujets par classes et le seuil ; elle ne prouve aucune ascendance. Source : contrat 08.')


def wilson_equation() -> str:
    formula = '<div>Sensibilité = VP / (VP + FN) ; taux de fausses alertes = FP / (FP + VN)</div>'
    numerator = 'p̂ + z²/(2n) ± z √[p̂(1 − p̂)/n + z²/(4n²)]'
    formula += '<div>IC<sub>Wilson</sub> = ' + _fraction(numerator, '1 + z²/n') + '</div>'
    return _equation(3, "Performance et incertitude descriptive", formula,
                     'VP/FN : atypiques détectés/manqués ; FP/VN : négatifs alertés/non alertés. p̂=x/n ; z≈1,96 pour 95 %. '
                     'n=4 pour la sensibilité, n=104 pour les fausses alertes dans la référence. '
                     'Les copies projetées et profils appariés ne sont pas de nouvelles unités indépendantes. '
                     'Un taux observé de 100 % ou de 0 % ne signifie pas certitude ; Wilson est descriptif sur ce banc. '
                     'Source : kinship_oracle.wilson_interval et CAL-REP-08.json.')


def cohort_equations() -> str:
    formula = '<div>C<sub>témoins</sub> = repr<sub>FID</sub>(N ∩ G ∩ W ∩ I<sub>07</sub> ∖ (E<sub>07</sub> ∪ E<sub>08</sub>))</div>'
    formula += '<div>repr<sub>FID</sub> = arg min<sub>i dans le FID</sub> (rang QC, F<sub>MISS</sub>, |Z<sub>HET</sub>|, ID technique)</div>'
    first = _equation(1, "Règle des témoins et représentant de famille technique", formula,
                      'N : non-porteurs selon génotype explicite ; G : rôle témoin autorisé ; W : inclusion genome-wide ; '
                      'I<sub>07</sub> : retenus dans la proposition d’entrée ; E<sub>07</sub>/E<sub>08</sub> : exclusions acceptées dans le test. '
                      'Le minimum est lexicographique : PASS avant ALERT, puis moindre manque, moindre |Z| et identifiant comme dernier départage. '
                      'Les porteurs indépendants suivent une règle analogue avec inclusion sur le chromosome cible. '
                      'Les cohortes « tous porteurs » et « chromosome cible » ne subissent pas ces exclusions. '
                      'Ici, les revues et FID sont des fixtures techniques ; cette règle ne qualifie pas la parenté de D21. Source : contrat 09.')
    formula = '<div>Δ<sub>répétition</sub> = Σ<sub>f ∈ F</sub> 𝟙[octets<sub>a</sub>(f) ≠ octets<sub>b</sub>(f)] = 0</div>'
    formula += '<div>Pour chaque cohorte : inclus + non inclus = 128</div>'
    return first + _equation(2, "Reproductibilité et conservation des décisions", formula,
                             'F contient les 11 artefacts scientifiques d’un gel : membres, résumé, décisions, rapport et sept fichiers keep. '
                             'Les audits horodatés et enveloppes de publication sont exclus de cette égalité. '
                             'Le graphique des effectifs utilise les 128 sujets, pas seulement les candidats ; chaque non-inclusion a un motif prioritaire unique. '
                             'Reproduire une décision n’en démontre pas la justesse biologique. Source : evaluate_structure.evaluate_freeze et contrats 09.')
