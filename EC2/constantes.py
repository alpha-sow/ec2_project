# -*- coding: Utf-8 -*-
# Nom de fichier: constantes.py
"""
 Module des constantes
"""
__version__ = '0.1'
# béton
CMAX = 90.  # MPa, compression limite pour le béton
GAMMAC = 1.5  # coefficient de sécurité sur le béton
ALPHAC = 1.  # coefficient effet à long terme en compression
ALPHACT = 1.  # coefficient effet à long terme en traction acier
GAMMAS = 1.15  # coefficient de sécurité sur l'acier
RHOS = 7850.  # kg/m3, masse volumique de l'acier
ES = 200000.  # MPa module d'élasticité de l'acier #numérique
ZERO = 1e-9  # utile pour les comparaisons numériques
GAMMACE = 1.2  # coefficient pour les instabilités
