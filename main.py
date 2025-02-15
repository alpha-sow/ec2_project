#!/usr/bin/python
#-*- Encoding: Utf-8 -*-
from scipy.optimize import fsolve
from scipy.integrate import quad
import numpy as np
import matplotlib.pyplot as plt
from math import pi

import EC2.BA as EC2
from voile import TraceAbaques

#matériaux
fck = 25. #MPa, résistance caractéristique du béton
fyk = 500. #MPa, résistance caractéristique de l'acier

#géométrie
bw, h = 0.25, 0.60 #m
d, dp = 0.510, 0.05 #m

#chargement Mg, Mq en MN.m
sol = [[0.127, 0.059],
[0.127, 0.250],]

psi2 = 0.3
print("ELU : ")
for [Mg, Mq] in sol:
    print("-----------------------------------------")
    print("Mg = {0:.3f} MN.m - Mq = {1:.3f} MN.m".format(Mg, Mq))
    MEdu = 1.35 * Mg + 1.5 * Mq
    EC2.AsRectELU(MEdu, bw, d, dp, fck, fyk, True)
print("-----------------------------------------")
EC2.AsminRect(bw, d, fck, fyk, True)
EC2.AsmaxRect(bw, h, True)

print("Voile")
bw = .25
h = 2.500
v, vp = h / 2., h / 2.
geom = {"RECT": [v, vp, bw]}
fck = 30.0
fyk = 500.0
epsilonuk = 5e-2
NEdu = np.array([-2, -3.871])
MEdu = np.array([-1.15, 0.])  # 0.05558 #0.043 #MN.m
soll = [NEdu, MEdu]
Aciers = np.array([[1608e-6, -1.186],
                   [1608e-6, -1.111],
                   [1608e-6, -1.036],
                   [1608e-6, -0.961],
                   [1608e-6, -0.886],
                   [1608e-6, -0.811],
                   [157e-6, -0.651],
                   [157e-6, -0.491],
                   [157e-6, -0.331],
                   [157e-6, -0.171],
                   [157e-6, -0.0011],
                   [157e-6, 0.149],
                   [157e-6, 0.309],
                   [157e-6, 0.469],
                   [157e-6, 0.629],
                   [157e-6, 0.789],
                   [1608e-6, 0.864],
                   [1608e-6, 0.939],
                   [1608e-6, 1.014],
                   [1608e-6, 1.089],
                   [1608e-6, 1.164],
                   [1608e-6, 1.239]])
[NMr, NMrpt] = TraceAbaques(v, vp, geom, fck, fyk, epsilonuk, Aciers, soll, True, [100, 500, 100])
plt.show()

plt.close("all")
### chapitre 7 : flexion à l'état limmite de service
print("\nExemple p.216")
bw, h = .2, .4
Ac, u = bw * h, 2. * (bw + h)
fck = 25.
t0 = 28.
RH = 80.
print(EC2.Ecm(fck))
typeCiment = "N"
phi = EC2.phiinfto(RH, Ac, u, fck, t0, typeCiment, True)

print("\nExemple de section rectangulaire p.221")
bw, h = .3, .5
d, dp = 0.9 * h, 0.05
Mg, Mq = .07, .04
psi2 = 0.3
fck, fyk = 25. , 500.
As, Asp = EC2.HA(6, 14), 0.
MEdc = Mg + Mq
MEdqp = Mg + psi2 * Mq
phic = MEdqp / MEdc * phi
EC2.VerifContRect(MEdc, As, Asp, bw, h, d, dp, fck, fyk, phic, 0.6, True)
#quelques petites différences dues à phi/1.05 et le module calculé et non lu au Tab3.1

### chapitre 8 : Flexion à l'état limite ultime
print("\nsection rectangulaire p.249")
bw, h = 0.3, 0.5
d, dp = 0.45, 0.05
fck, fyk = 25. , 500.
Mg, Mq = 0.075, 0.037
MEdu = 1.35 * Mg + 1.5 * Mq
EC2.AsminRect(bw, d, fck, fyk, True)
[As1, Asp1] = EC2.AsRectELU(MEdu, bw, d, dp, fck, fyk, True)
d = .46
[As1, Asp1] = EC2.AsRectELU(MEdu, bw, d, dp, fck, fyk, True)

print("\nsuite exemple p.250")
k, epsilonuk = 1.05, 0.025
d = 0.46
[As2, Asp2, sigc, sigs] = EC2.AsRectELU2(MEdu, bw, d, dp, fck, fyk, k, epsilonuk, True)
print(As1 * 0.46 / 0.46)

print("\nexemple 1, p.265")
bw, d = 0.3, 0.45
dp = 0.0315
fck, fyk = 25., 500.
MEdu = 0.278
EC2.AsRectELU(MEdu, bw, d, dp, fck, fyk, True)
print("\nExemple 2, p.266")
MEdu = 0.422
EC2.AsRectELU(MEdu, bw, d, dp, fck, fyk, True)

#### chapitre 9 : effort tranchant torsion poiçonnement
print("\nExemple p.274")
bw, d = .25, .4
fck, fyk = 25., 500.
VEdu = 0.2031
theta1 = np.arctan(1. / 2.5)
EC2.AswsV(VEdu, d, fck, fyk, True, theta = theta1)
EC2.Aswsmin(bw, fck, fyk, True)

print("\nExemple p.297")
vEdu = .725
hf = .15
thetaf1 = np.arctan(1. / 2)
fyk = 500.
EC2.Asfsf(vEdu, hf, fck, fyk, thetaf1, True)
EC2.Aswsmin(hf, fck, fyk, True)

print("\nExemple de poutre calculée  la torsion p.307")
bw, h = 0.3, .5
d, dp = 0.45, 0.05
fck, fyk = 25., 500.
MEdu = .160
EC2.AsRectELU(MEdu, bw, d, dp, fck, fyk, True)
VEdu = .160
VRdmax1 = EC2.VRdmaxV(bw, d, fck, True)
EC2.AswsV(VEdu, d, fck, fyk, True)
TEdu = .0464
tefi, Ak, uk = 0.094, 0.0836, 1.22
TRdmax1 = EC2.TRdmax(Ak, tefi, fck, True)
print(VEdu / VRdmax1 + TEdu / TRdmax1)
EC2.AswsT(TEdu, Ak, fyk, True)
EC2.AslT(TEdu, Ak, uk, fyk, True)

print("\nExemple de poinçonnement d'un plancher-dalle")
bw, d = 0.3, 0.198
fck, fyk = 25., 500.
vEdu = 0.645
s = .15
u1 = 3.69
Aslx = EC2.HA(1, 12) / s * bw
Asly = Aslx
typeBA1 = "dallepoutre"
EC2.vRdcP(bw, d, Aslx, Asly, typeBA1, fck, True)
EC2.vRdmaxP(fck, True) #petite évolution EC2
EC2.AswsuP(vEdu, bw, d, Aslx, Asly, typeBA1, fck, fyk, True)
EC2.AswminuP(fck, fyk, True)
EC2.DispoPoinc(d)

### chapitre 11
print("\nPoteau (Exemple)")
bw, h = 0.4, 0.2
Geom = {"RECT": [h / 2, h / 2., bw]}
dp = 0.031
d = h - dp
Ac, u = bw * h, 2. * (bw + h)
fck, fyk = 25., 500.
As = EC2.HA(3, 10)
Asp = As
Aciers = np.array([[As, (h / 2 - d)],
                    [Asp, h / 2 - dp]])
Ng, Nq = -0.36, -0.16
NEdu = 1.35 * Ng + 1.5 * Nq
NEdc, NEdqp = Ng + Nq, Ng + psi2 * Nq
L0 = 2.6
psi2 = 0.3
Ecm1 = EC2.Ecm(fck)
RH = 50.
t0 = 28.
typeCiment = "N"
phi = EC2.phiinfto(RH, Ac, u, fck, t0, typeCiment, True)
phieff = phi * NEdqp / NEdu
print("1+phie = {:.3f}".format(1+phieff))

print("\nMéthode générale")
e0 = EC2.e0i(L0)
lambda1 = L0 / EC2.RayonGirationGeom(Geom)
MEdu0 = NEdu * e0
v, vp = h / 2.,  h / 2.
MEduMG = EC2.MethodeGenerale(e0, L0, v, vp, Geom, phi, fck, Aciers, fyk, NEdu,
                    kimax = 80e-3, kimin = 0., nbpts = 50, aff = True)
print("\nTest avec valeur de NRdu pour la tangence des 2 courbes")
MEduMG = EC2.MethodeGenerale(e0, L0, v, vp, Geom, phi, fck, Aciers, fyk, -0.9215,
                    kimax = 80e-3, kimin = 0., nbpts = 50, aff = True)
# print("")
# EC2.RigiditeNominale(lambda1, L0, Geom, phieff, fck, Aciers, fyk, NEdu, MEdu0, 1, True)
# EC2.CourbureNominale(lambda1, L0, Geom, phieff, fck, Aciers, fyk, NEdu, MEdu0, 8, True)
# rho = (As + Asp) / bw / h
# print("NRdu = {:.3f} MN".format(EC2.NRdRect(rho, bw, h, dp, fck, fyk, lambda1, 1)))