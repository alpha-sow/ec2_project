import numpy as np
import EC2.materiaux as MAT


def VecteurAbaquesNM(v, vp, d, dp, fck, epsilonuk, bNonSymetrique=False, nbpts=[100, 100, 100]):
    """
        calcul des paramètres eps0 et ki de la droite de déformation (3 pivots)
        entrée : section_height_gravity, section_width_gravity, d, dp [m], resistance_beton [MPa],
            deformation_epsilonuk [],
        sortie : eps0 [], lambda [m-1]
    """
    epsud = MAT.epsilonud(epsilonuk)
    epscu2 = - MAT.epsiloncu2(fck)
    epsc2 = - MAT.epsilonc2(fck)
    h = v + vp

    # partie supérieure
    nbpts_z1, nbpts_z2, nbpts_z3 = nbpts

    # zone 1, pivot A, traction simple, traction excentrée, flexion simple
    eps_p = np.linspace(epsud, epscu2, nbpts_z1)
    ki = [(epsP1 - epsud) / d for epsP1 in eps_p]
    eps0 = [(epsud * v + epsP1 * (v - d)) / d for epsP1 in eps_p]

    # zone 2, pivot B, flexion simple, flexion composée
    epsM2 = (epsud * v + epscu2 * (v - d)) / d + (epscu2 + epsud) / d * vp
    epsM = np.linspace(epsM2, 0, nbpts_z2)
    ki = np.append(ki, [(epscu2 - epsM21) / h for epsM21 in epsM])
    eps0 = np.append(eps0, [(epscu2 * vp + epsM21 * v) / h for epsM21 in epsM])

    # zone 3, pivot C, flexion composée, compression simple
    dc = h * (epscu2 - epsc2) / epscu2
    epsP = np.linspace(epscu2, epsc2, nbpts_z3)
    ki = np.append(ki, [(epsP1 - epsc2) / dc for epsP1 in epsP])
    eps0 = np.append(eps0, [(epsc2 * v - epsP1 * (v - dc)) / dc for epsP1 in epsP])

    if bNonSymetrique:
        # partie inférieure si besoin
        # zone 3
        dcp = h * epsc2 / epscu2
        epsP = np.linspace(epsc2, 0, nbpts_z3)
        ki = np.append(ki, [(epsP1 - epsc2) / dcp for epsP1 in epsP])
        eps0 = np.append(eps0, [(epsc2 * v - epsP1 * (v - dcp)) / dcp for epsP1 in epsP])

        # zone 2
        epsP2 = (epsud * vp + epscu2 * (v - dp)) / (h - dp) + (epsud - epscu2) / (h - dp) * v
        epsP = np.linspace(0, epsP2, nbpts_z2)
        ki = np.append(ki, [(epsP21 - epscu2) / h for epsP21 in epsP])
        eps0 = np.append(eps0, [(epscu2 * v + epsP21 * vp) / h for epsP21 in epsP])

        # zone 1
        epsM = np.linspace(epscu2, epsud, nbpts_z1)
        ki = np.append(ki, [(epsud - epsM1) / (h - dp) for epsM1 in epsM])
        eps0 = np.append(eps0, [(epsud * vp + epsM1 * (v - dp)) / (h - dp) for epsM1 in epsM])

    return [eps0, ki]
