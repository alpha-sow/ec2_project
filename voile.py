import matplotlib.pylab as plt
import numpy as np
import EC2.materiaux as MAT
from scipy.integrate import quad

from points_particuliers import PointsParticuliersAbaquesNM
from vecteur_abaques import VecteurAbaquesNM


def TraceAbaques(v, vp, geom, fck, fyk, epsilonuk, Aciers, soll=[[], []], bNonSymetrique=False, nbpts=[100, 100, 100]):
    """
        tracé de l'abaque d'interaction pour la section
        entrée : section_height_gravity, section_width_gravity, geometrie_section [m], resistance_beton,
            resistance_acier [MPa], deformation_epsilonuk [],
        sollicitation = [[MN], [MN.m]]
        sortie : NMr, NMrpt [MN, MM.m, MN, MN.m]
        geometrie_section est un dictionnaire pour
        geometrie_section = {"RECT": [section_height_gravity, section_width_gravity, section_width]}
    """
    d, dp = v - min(Aciers[:, 0]), vp - max(Aciers[:, 0])
    h = v + vp
    [eps0, ki] = VecteurAbaquesNM(v, vp, d, dp, fck, epsilonuk, bNonSymetrique, nbpts)
    [eps0p, kip] = PointsParticuliersAbaquesNM(v, vp, d, dp, fck, epsilonuk, bNonSymetrique)
    NMr = np.array([NM_C_S(eps0C, kiC, v, vp, geom, fck, fyk, Aciers) for eps0C, kiC in zip(eps0, ki)])
    NMrPt = np.array([NM_C_S(eps0C, kiC, v, vp, geom, fck, fyk, Aciers) for eps0C, kiC in zip(eps0p, kip)])
    plt.figure()
    plt.plot(NMr[:, 0], NMr[:, 1], 'b')
    plt.plot(NMrPt[:, 0], NMrPt[:, 1], 'g*')
    plt.plot(NMr[:, 0] + NMr[:, 2], NMr[:, 1] + NMr[:, 3], 'r-')
    plt.plot(NMrPt[:, 0] + NMrPt[:, 2], NMrPt[:, 1] + NMrPt[:, 3], 'g*')
    # placement des sollicitations
    NEdu, MEdu = soll
    plt.plot(NEdu, MEdu, 'r*')
    plt.grid('on')
    plt.xlabel(r'$N_{Rd}$ [MN]')
    plt.ylabel(r'$M_{Rd}$ [MN.m]')
    if "RECT" in geom:
        v, vp, bw = geom["RECT"]
        plt.title(u"Abaque d'interaction d'une section rectangulaire - {: .0f}x{: .0f} mm2".format(bw * 1e3, h * 1e3))
        plt.savefig('./images/AbaqueInteractionRect-{:.0f}x{:.0f}.pdf'.format(bw * 1e3, h * 1e3))
    return [NMr, NMrPt]


def NM_C_S(eps0, ki, v, vp, geom, fck, fyk, Aciers):
    """
        efforts résistants béton et acier
        entrée : eps0 [], lambda [m-1], section_height_gravity,
            section_width_gravity, d, dp, [m], resistance_beton, resistance_acier,[MPa], Aciers [m2, m],
        sortie : Nrc [MN], Mrc [MN.m] NRs [MN] MRs [MN.m]
    """
    nrc = Nc(eps0, ki, v, vp, geom, fck)
    mrc = Mc(eps0, ki, v, vp, geom, fck)
    nrs = sum([MAT.sigmas1(eps0 + ki * ysi, fyk) * Asi for Asi, ysi in Aciers])
    mrs = sum([- ysi * MAT.sigmas1(eps0 + ki * ysi, fyk) * Asi for Asi, ysi in Aciers])
    return [nrc, mrc, nrs, mrs]


def b(y, geom):
    """
        largeur de la section à traiter en fonction de geometrie_section
        entrée : y [m], geometrie_section
        sortie : [m]
    """
    if "RECT" in geom:
        v, vp, bw = geom["RECT"]
        return bw
    if "TE" in geom:
        v, vp, hf, bw, beff = geom["TE"]
        if y > v - hf:
            return beff
        else:
            return bw
    if "CIRC" in geom:
        v, vp, R = geom["CIRC"]
        return 2. * (R * R - y * y) ** .5
    raise ValueError("Cle non reconnue dans b")


def Nc(eps0, ki, v, vp, geom, fck):
    """
        Effort normal résistant béton pour une section rectangulaire
        entrée : eps0 [], ki [m-1], section_height_gravity, section_width_gravity, geometrie_section [m],
            resistance_beton [MPa]
        sortie : [MN]
    """

    def integrand(y, eps0, ki, fck, geom):
        return b(y, geom) * MAT.sigmac2(eps0 + ki * y, fck)

    return quad(integrand, -vp, v, args=(eps0, ki, fck, geom))[0]


def Mc(eps0, ki, v, vp, geom, fck):
    """
        Moment fléchissant résistant béton pour une section rectangulaire
        entrée : eps0 [], ki [m-1], section_height_gravity, section_width_gravity [m], resistance_beton [MPa]
        sortie : [MN]
    """

    def integrand(y, eps0, ki, fck, geom):
        return -(b(y, geom) * MAT.sigmac2(eps0 + ki * y, fck) * y)

    return quad(integrand, -vp, v, args=(eps0, ki, fck, geom))[0]

# def TraceAbaquesNormalise(v, vp, geom, fck, fyk, epsilonuk, Aciers, soll=[[], []], bNonSymetrique=False,
#                           nbpts=[100, 100, 100]):
#     fcd1, fyd1 = MAT.fcd(fck), MAT.fyd(fyk)
#     d, dp = v - min(Aciers[:, 0]), vp - max(Aciers[:, 0])
#     h = v + vp
#     [eps0, ki] = VecteurAbaquesNM(v, vp, d, dp, fck, epsilonuk, bNonSymetrique, nbpts)
#     [eps0p, kip] = PointsParticuliersAbaquesNM(v, vp, d, dp, fck, epsilonuk, bNonSymetrique)
#     NMr = np.array([NM_C_S(eps0C, kiC, v, vp, geom, fck, fyk, Aciers) for eps0C, kiC in zip(eps0, ki)])
#     NMrPt = np.array([NM_C_S(eps0C, kiC, v, vp, geom, fck, fyk, Aciers) for eps0C, kiC in zip(eps0p, kip)])
#     if ("RECT" in geom) or ("CIRC" in geom):
#         Astot = sum(Aciers[:, 0])
#         if "RECT" in geom:
#             v, vp, bw = geom["RECT"]
#         if "CIRC" in geom:
#             v, vp, R = geom["CIRC"]
#             bw = math.pi * R / 2.
#         [nurc, murc, nurs, murs] = [NMr[:, 0] / bw / h / fcd1,
#                                     NMr[:, 1] / bw / h / h / fcd1,
#                                     NMr[:, 2] / Astot / fyd1,
#                                     NMr[:, 3] / Astot / h / fyd1]
#         [nurcPt, murcPt, nursPt, mursPt] = [NMrPt[:, 0] / bw / h / fcd1, NMrPt[:, 1] / bw / h / h / fcd1,
#                                             NMrPt[:, 2] / Astot / fyd1,
#                                             NMrPt[:, 3] / h / Astot / fyd1]
#         n1, n2 = [nurc, murc, nurs, murs], [nurcPt, murcPt, nursPt, mursPt]
#         # courbes normalisées
#         plt.figure()
#         omega1 = Astot * fyd1 / (bw * h * fcd1)
#         nbCourbes = 20.
#         nu1 = nurc + omega1 * nurs
#         mu1 = murc + omega1 * murs
#         plt.plot(nu1, mu1, 'r-')
#         for i in range(int(nbCourbes) + 1):
#             omega = i / nbCourbes
#         # courbes
#         nu1 = nurc + omega * nurs
#         mu1 = murc + omega * murs
#         plt.plot(nu1, mu1, 'b-')  # points particuliers
#         nu1Pt = nurcPt + omega * nursPt
#         mu1Pt = murcPt + omega * mursPt
#         plt.plot(nu1Pt, mu1Pt, 'r.')
