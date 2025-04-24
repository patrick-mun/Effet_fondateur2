#!/bin/bash

# Lancer tous les tests et générer un rapport HTML
PYTHONPATH=$(pwd) pytest test/ --html=rapport_tests.html --self-contained-html

# Ouvrir le rapport HTML automatiquement (MacOS)
open rapport_tests.html
