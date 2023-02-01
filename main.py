import matplotlib.pylab as plt
import numpy as np
from voile import TraceAbaques

if __name__ == "__main__":
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
