import math
import numpy as np
import astropy.units as u

np.seterr(over='ignore')


def reddenSF_info():
    reddenSF_info = ("ebmv             \"\"        0.001       0.0  0.0       1.0    1.0   1e-3",)
    return reddenSF_info


# Constants
rv = 4.05 # assumed R_V value


def calzetti(E, ebmv):
    """
    Calculates the attenuation/extinction for a given E(B-V) and wavelength according to Calzetti et al.
    (2000, ApJ, 533, 682) from 2.20 um to 0.15um, and down to 0.09 um using Reddy et al. (2016, Apj, 828, 107).
    """

    rlambda = (E * u.keV).to('um', equivalencies=u.spectral()).value
    x = 1 / rlambda

    if 0.64 <= rlambda <= 2.5:
        k_wl = 2.659 * (-1.857 + 1.040 * x) + rv
    elif 0.15 <= rlambda < 0.64:
        k_wl = 2.659*(-2.156 + 1.509*x - 0.198*x**2 + 0.011*x**3) + rv
    elif 0.085 <= rlambda < 0.15:
        k_wl = 4.126 + 0.931*x
    else:
        k_wl = 0

    al = k_wl*ebmv
    return 10. ** (-al / 2.512)


def reddenSF(engs, params, flux):
    ebv = params[0]

    for i in range(len(engs) - 1):
        E = 0.5 * (engs[i + 1] + engs[i])
        flux[i] = calzetti(E, ebv)
