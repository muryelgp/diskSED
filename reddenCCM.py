import math
import numpy as np
import astropy.units as u

np.seterr(over='ignore')


def reddenCCM_info():
    reddenCCM_info = ("ebmv             \"\"        0.001       0.0  0.0       1.0    1.0   1e-3",)
    return reddenCCM_info


# Constants
rv = 3.1  # assumed R_V value


def cardelli(E, ebmv):
    """
    Calculates the extinction for a given E(B-V) and wavelength according to
    Cardelli et al. (1989, ApJ, 345, 245). Units of wavelength are angstroms.
    """

    rlambda = (E * u.keV).to('AA', equivalencies=u.spectral()).value
    av = ebmv * rv
    x = 1.e4 / rlambda
    y = x - 1.82

    if 0.3 <= x < 1.1:
        ax = 0.574 * x ** 1.61
        bx = -0.527 * x ** 1.61
    elif 1.1 <= x <= 3.3:
        ax = 1. + 0.17699 * y - 0.50447 * y ** 2 - 0.02427 * y ** 3 + 0.72085 * y ** 4 + 0.01979 * y ** 5 - 0.77530 * y ** 6 + 0.32999 * y ** 7
        bx = 1.41338 * y + 2.28305 * y ** 2 + 1.07233 * y ** 3 - 5.38434 * y ** 4 -0.62251 * y ** 5 + 5.30260 * y ** 6 - 2.09002 * y ** 7
    elif 3.3 < x <= 8:
        if x >= 5.9:
            ax = 1.752 - 0.316 * x - 0.104 / ((x - 4.67) ** 2 + 0.341) - 0.04473 * (x - 5.9) ** 2 - 0.009779 * (x - 5.9) ** 3
            bx = -3.090 + 1.825 * x + 1.206 / ((x - 4.62) ** 2 + 0.263) + 0.2130 * (x - 5.9) ** 2 + 0.1207 * (x - 5.9) ** 3
        else:
            ax = 1.752 - 0.316 * x - 0.104 / ((x - 4.67) ** 2 + 0.341)
            bx = -3.090 + 1.825 * x + 1.206 / ((x - 4.62) ** 2 + 0.263)
    elif 8 < x <= 11.2:
        ax = -1.073 - 0.628 * (x - 8) + 0.137 * (x - 8) ** 2 - 0.07 * (x - 8) ** 3
        bx = 13.670 + 4.257 * (x - 8) - 0.42 * (x - 8) ** 2 + 0.374 * (x - 8) ** 3
    else:
        ax, bx = 0, 0

    al = av * (ax + bx / rv)
    return 10. ** (-al / 2.512)


def reddenCCM(engs, params, flux):
    ebv = params[0]

    for i in range(len(engs) - 1):
        E = 0.5 * (engs[i + 1] + engs[i])
        flux[i] = cardelli(E, ebv)
