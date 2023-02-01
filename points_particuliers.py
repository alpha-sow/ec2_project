import numpy as np
import EC2.materiaux as MAT


def PointsParticuliersAbaquesNM(v, vp, d, dp, fck, epsilonuk, bNonSymetrique=False):
    """
        calcul des paramètres particuliers aux changements de zone pour les 3 pivots
        entrée : section_height_gravity, section_width_gravity, d, dp [m], resistance_beton [MPa],
            deformation_epsilonuk [],
        sortie : eps0 [], lambda [m-1]
    """
    epsud = MAT.epsilonud(epsilonuk)
    epscu2 = - MAT.epsiloncu2(fck)
    epsc2 = - MAT.epsilonc2(fck)
    h = v + vp

    # partie supérieure
    # zone 1, pivot A, pivot B
    eps_p = [epsud, epscu2]
    ki = [(eps_p1 - epsud) / d for eps_p1 in eps_p]
    eps0 = [(epsud * v - epsP1 * (v - d)) / d for epsP1 in eps_p]

    # zone 2, pivot B, flexion simple, flexion composée
    eps_m2 = (epsud * v + epscu2 * (v - d)) / d + (epscu2 + epsud) / d * vp
    eps_m = [eps_m2, 0]
    ki = np.append(ki, [(epscu2 - epsM21) / h for epsM21 in eps_m])
    eps0 = np.append(eps0, [(epscu2 * vp + epsM21 * v) / h for epsM21 in eps_m])

    # zone 3, pivot C, flexion composée, compression simple
    dc = h * (epscu2 - epsc2) / epscu2
    eps_p = [epscu2, epsc2]
    ki = np.append(ki, [(epsc2 - epsP1) / dc for epsP1 in eps_p])
    eps0 = np.append(eps0, [(epsc2 * v - epsP1 * (v - dc)) / dc for epsP1 in eps_p])

    if bNonSymetrique:
        # partie inférieure si besoin
        # zone 3
        dcp = h * epsc2 / epscu2
        eps_p = [epsc2, 0]
        ki = np.append(ki, [(epsP1 - epsc2) / dcp for epsP1 in eps_p])
        eps0 = np.append(eps0, [(epsc2 * v - epsP1 * (v - dcp)) / dcp for epsP1 in eps_p])

        # zone 2
        eps_p2 = (epsud * vp + epscu2 * (v - dp)) / (h - dp) + (epsud - epscu2) / (h - dp) * v
        eps_p = [0, eps_p2]
        ki = np.append(ki, [(epsP21 - epscu2) / h for epsP21 in eps_p])
        eps0 = np.append(eps0, [(epscu2 * v + epsP21 * vp) / h for epsP21 in eps_p])

        # zone 1
        eps_m = [epscu2, epsud]
        ki = np.append(ki, [(epsud - epsM1) / (h - dp) for epsM1 in eps_m])
        eps0 = np.append(eps0, [(epsud * vp + epsM1 * (v - dp)) / (h - dp) for epsM1 in eps_m])

    return [eps0, ki]
