import matplotlib.pylab as plt
import numpy as np
import EC2.materiaux as MAT
from scipy.integrate import quad

from points_particuliers import PointsParticuliersAbaquesNM
from vecteur_abaques import VecteurAbaquesNM


def TraceAbaques(v, vp, geom, fck, fyk, epsilonuk, Aciers, soll=None, bNonSymetrique=False, nbpts=None):
    """
        tracé de l'abaque d'interaction pour la section
        entrée : section_height_gravity, section_width_gravity, geometrie_section [m], resistance_beton,
            resistance_acier [MPa], deformation_epsilonuk [],
        sollicitation = [[MN], [MN.m]]
        sortie : NMr, NMrpt [MN, MM.m, MN, MN.m]
        geometrie_section est un dictionnaire pour
        geometrie_section = {"RECT": [section_height_gravity, section_width_gravity, section_width]}
    """
    if nbpts is None:
        nbpts = [100, 100, 100]
    if soll is None:
        soll = [[], []]
    d, dp = v - min(Aciers[:, 0]), vp - max(Aciers[:, 0])
    h = v + vp
    [eps0, ki] = VecteurAbaquesNM(v, vp, d, dp, fck, epsilonuk, bNonSymetrique, nbpts)
    [eps0p, kip] = PointsParticuliersAbaquesNM(v, vp, d, dp, fck, epsilonuk, bNonSymetrique)
    n_mr = np.array([NM_C_S(eps0C, kiC, v, vp, geom, fck, fyk, Aciers) for eps0C, kiC in zip(eps0, ki)])
    n_mr_pt = np.array([NM_C_S(eps0C, kiC, v, vp, geom, fck, fyk, Aciers) for eps0C, kiC in zip(eps0p, kip)])
    plt.figure()
    plt.plot(n_mr[:, 0], n_mr[:, 1], 'b')
    plt.plot(n_mr_pt[:, 0], n_mr_pt[:, 1], 'g*')
    plt.plot(n_mr[:, 0] + n_mr[:, 2], n_mr[:, 1] + n_mr[:, 3], 'r-')
    plt.plot(n_mr_pt[:, 0] + n_mr_pt[:, 2], n_mr_pt[:, 1] + n_mr_pt[:, 3], 'g*')
    # placement des sollicitations
    n_edu, m_edu = soll
    plt.plot(n_edu, m_edu, 'r*')
    plt.grid('on')
    plt.xlabel(r'$N_{Rd}$ [MN]')
    plt.ylabel(r'$M_{Rd}$ [MN.m]')
    if "RECT" in geom:
        v, vp, bw = geom["RECT"]
        plt.title(u"Abaque d'interaction d'une section rectangulaire - {: .0f}x{: .0f} mm2".format(bw * 1e3, h * 1e3))
        plt.savefig('./images/AbaqueInteractionRect-{:.0f}x{:.0f}.pdf'.format(bw * 1e3, h * 1e3))
    return [n_mr, n_mr_pt]


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
        v, vp, r = geom["CIRC"]
        return 2. * (r * r - y * y) ** .5
    raise ValueError("Cle non reconnue dans b")


def Nc(eps0, ki, v, vp, geom, fck):
    """
        Effort normal résistant béton pour une section rectangulaire
        entrée : eps0 [], ki [m-1], section_height_gravity, section_width_gravity, geometrie_section [m],
            resistance_beton [MPa]
        sortie : [MN]
    """

    def integrand(y, eps0_i, ki_i, fck_i, geom_i):
        return b(y, geom_i) * MAT.sigmac2(eps0_i + ki_i * y, fck_i)

    return quad(integrand, -vp, v, args=(eps0, ki, fck, geom))[0]


def Mc(eps0, ki, v, vp, geom, fck):
    """
        Moment fléchissant résistant béton pour une section rectangulaire
        entrée : eps0 [], ki [m-1], section_height_gravity, section_width_gravity [m], resistance_beton [MPa]
        sortie : [MN]
    """

    def integrand(y, eps0_i, ki_i, fck_i, geom_i):
        return -(b(y, geom_i) * MAT.sigmac2(eps0_i + ki_i * y, fck_i) * y)

    return quad(integrand, -vp, v, args=(eps0, ki, fck, geom))[0]
