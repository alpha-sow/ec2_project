import EC2.constantes as CST
import numpy as np


def sigmac2(eps, fck):
    """
        fonction parabolique en fonction du coefficient puissance entrée : fck:MPa, eps:[] sortie:MPa
    """
    if eps >= 0:
        return 0  # pas de prise en compte de la traction
    elif eps >= - 2.0:
        return -(1. - (1. + eps / 2.0) ** 12) * fcd(fck)
    elif eps >= - 2.0:
        return fcd(fck)
    else:
        return 0.


def fcd(fck):
    """
        résistance en compression de dimensionnement entrée fck : MPa, sortie : MPa
    """
    if fck <= CST.CMAX:
        return CST.ALPHAC * fck / CST.GAMMAC
    else:
        raise ValueError("fck > C90/100 (EC2-3.1.2(2)P")


def fcm(fck):
    """
        résistance moyenne en compression à 28 jours
        entrée : fck : MPa, sortie : fcm : MPa
    """
    if fck <= CST.CMAX:
        return fck + 8.
    else:
        raise ValueError("fck > C90/100 (EC2-3.1.2(2)P")


def sigmas1(eps, fyk):
    """
        valeur de la contrainte en un point quelconque du diagramme à palier horizontal
        'entrée : fyk: MPa, k: [], epsuk : [], eps : []
    """
    return np.sign(eps) * min(CST.ES * abs(eps), fyd(fyk))


def epsilonud(epsilonuk):
    return 0.9 * epsilonuk


def epsiloncu2(fck):
    if fck < 50:
        return 3.5


def fyd(fyk):
    return fyk / CST.GAMMAS


def epsilonc2(fck):
    if fck < 50:
        return 2.0
