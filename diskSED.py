import math
import numpy as np
import astropy.units as u
from astropy.constants import h, c, k_B

# Ignore overflows in numpy operations
np.seterr(over='ignore')


def diskSED_info():
    """
    Function to provide information about the parameters pyXSPEC needs this.

    Returns:
        diskRout1_info (tuple): Contains strings with parameter names, units, and ranges.
    """
    diskRout1_info = ("R_in*             km       1e6       1   1        1    1e10     1e10",
                      "T_p             Kelvin       5e5        1e4     1e4         1e7     1e7    1e7",
                      "R_ratio             \"\"         100        0.1    5         5     1e5     1e5",
                      "D_Mpc             Mpc         10        10    10       1e5     1e5     1"
                      )
    return diskRout1_info


# Precompute constants outside the function for efficiency

# Constants for color correction function
f_inf = 2.3
nu_b = delt_nu = 5e15
h_kev_Hz = h.to('keV/Hz').value  # Planck constant in keV/Hz
K_to_keV = (1*u.K).to('keV', equivalencies=u.temperature_energy()).value

# Disk emission constant, converting to desired units
consts = (u.km ** 2 / (h ** 3 * c ** 2 * u.Mpc ** 2)).cgs.to('1/(s cm**2 keV**3)').value
k_B = k_B.to('keV/Kelvin')


def f_col_T(T):
    """
    Color correction function.

    Analytical expression given by Chiang, J. 2002 (ApJ, 572, 79), based on Hubeny, I.,
    et  al.,2001, (ApJ, 559, 680) numerical results.

    Args:
        T (float): Temperature in keV.
    Returns:
        f_col (float): Color correction factor.
    """
    f_inf = 2.3
    nu_b = delt_nu = 5e15
    h_kev_Hz = h.to('keV/Hz').value  # Planck constant in keV/Hz
    K_to_keV = (1 * u.K).to('keV', equivalencies=u.temperature_energy()).value

    T = T*K_to_keV
    nu_p = 2.82 * T / h_kev_Hz
    f_col = f_inf - (((f_inf - 1) * (1 + np.exp(-nu_b / delt_nu))) / (1 + np.exp((nu_p - nu_b) / delt_nu)))
    return f_col


def diskSED_sub(E, R_in, T_p, R_ratio, D_mpc):
    """
    Subroutine to compute the disk spectrum for a given energy.

    Expression given in section 5.5 of Frank et al. 2002, but including temperature dependent color correction.

    Args:
        E (float): Energy in keV.
        R_in* (float): Inner radius in km times sqrt(cos i), where i is the inclination.
        T_p (float): Peak temperature in Kelvin.
        R_ratio (float): Ratio of outer to inner disk radius.
        D_mpc (float): Distance to the source in Mpc.
    Returns:
        dbb_sub_E (float): Computed disk blackbody spectrum at energy E.
    """

    log_R_ratio = np.log10(R_ratio)
    rolog = np.log(10 ** log_R_ratio)
    nstep = math.ceil(rolog / 0.05)
    dr = rolog / float(nstep)
    i = np.arange(1, nstep + 1)
    r2log = i * dr
    r1log = (i - 1) * dr
    rlog = 0.5 * (r1log + r2log)
    x = np.exp(rlog)
    tau = ((1. - 1. / np.sqrt(x)) / x ** 3) ** (1. / 4.)
    T_x = 2.05 * T_p * tau
    amp = (4 * np.pi * consts / f_col_T(T_x) ** 4) * (R_in ** 2 / D_mpc ** 2)
    new_T = f_col_T(T_x) * T_x
    y = E / (k_B*new_T).value
    factor = 1.0 / (np.exp(y) - 1.0)
    fac = factor * x ** 2
    s = np.sum(fac * amp) * dr
    dbb_sub_E = s * E ** 2
    return dbb_sub_E


def diskSED(engs, params, flux):
    """
    Function to compute the disk spectral energy distribution (SED).

    This is the format that pyXSPEC needs, need to have these three array/tuple and cannot return anything,
    just modifies the input flux array.


    Args:
        engs (array): Energy grid in keV.
        params (array): Array containing the parameters [R_in, T_p, R_ratio, D_mpc].
        flux (array): Array to store the computed flux values.
    """
    R_in = params[0]  # Inner radius in km
    T_p = params[1]  # Peak temperature in Kelvin
    R_ratio = params[2]  # Ratio of outer to inner disk radius
    D_mpc = params[3]  # Distance to the source in Mpc
    norm = params[4]  # I don't know how to get read of the normalization, should be set to 1 in the fitting.
    if norm == 0:
        flux = np.zeros(len(flux))
    else:
        # Compute the flux for each energy bin. This still has a for loop, but could be vectorized for efficiency.
        for i in range(len(engs) - 1):
            E = 0.5 * (engs[i + 1] + engs[i])  # Midpoint energy of the bin
            dE = engs[i + 1] - engs[i]  # Width of the energy bin
            flux[i] = diskSED_sub(E, R_in, T_p, R_ratio, D_mpc) * dE  # Compute and store flux
