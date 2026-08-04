#!/bin/bash
set -euo pipefail

# Lancer tous les tests et générer un rapport HTML
PROJECT_DIR=$(cd "$(dirname "$0")" && pwd)
cd "$PROJECT_DIR"
TEST_PATHS=(tests)
if [[ -d test ]]; then
    TEST_PATHS+=(test)
fi
PYTHONPATH="$PROJECT_DIR" .venv/bin/python -m pytest "${TEST_PATHS[@]}" --html=rapport_tests.html --self-contained-html

# Ouvrir le rapport HTML automatiquement (MacOS)
open rapport_tests.html
