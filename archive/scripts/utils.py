# utils.py
import os
import subprocess
import logging

def safe_mkdir(path):
    """Crée un dossier s'il n'existe pas."""
    os.makedirs(path, exist_ok=True)

def run_cmd(command, description=None):
    """Exécute une commande shell avec logs."""
    try:
        if description:
            logging.info(description)
        logging.info(f"Commande : {command}")
        subprocess.run(command, shell=True, check=True)
    except subprocess.CalledProcessError as e:
        logging.error(f"Erreur : {e}")
        raise
