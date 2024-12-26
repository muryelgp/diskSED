import matplotlib.pyplot as plt
import numpy as np
import astropy.units as u
from astropy.constants import c, G
from .diskSED import diskSED, diskSED_sub
# Convert constants to CGS units
c = c.cgs
G = G.cgs

def get_bolometric(R_in_star, T_p, R_ratio):
    """
    Calculate the bolometric luminosity.

    Parameters:
    R_in_star (float, int, np.float64, or np.ndarray): Inner radius times sqrt(cos i) in km, where i is the inclination.
    T_p (float, int, np.float64, or np.ndarray): Peak temperature in keV.
    R_ratio (float, int, np.float64, or np.ndarray): Ratio of outer to inner radius.

    Returns:
    np.ndarray: Bolometric luminosity in erg/s.
    """
    # Ensure inputs are numpy arrays
    if isinstance(R_in_star, (float, int, np.float64)):
        R_in_star, T_p, R_ratio = np.array([R_in_star]), np.array([T_p]), np.array([R_ratio])
    elif isinstance(R_in_star, list):
        R_in_star, T_p, R_ratio = np.array(R_in_star), np.array(T_p), np.array(R_ratio)

    size = len(R_in_star)
    Lbol = np.zeros(size)
    E_range = np.logspace(-4, 1, 2000)  # Energy range in keV
    A_range = (E_range * u.keV).to('AA', equivalencies=u.spectral()).value  # Convert energy to Angstrom

    for i in range(size):
        # Calculate spectrum for each energy
        ct_range_kev = [diskSED_sub(E_range[j], R_in_star[i], T_p[i], R_ratio[i], 100) for j in range(len(E_range))]
        EF_E = ((ct_range_kev * E_range ** 2) * (u.keV ** 2 / (u.s * u.cm ** 2 * u.keV))).to('erg/(cm**2 s)')
        scale = 4 * np.pi * (100 * u.Mpc.to('cm')) ** 2  # Scaling factor for luminosity distance
        F_wl = (EF_E / A_range).value

        Lbol[i] = np.trapz(np.flip(F_wl), np.flip(A_range)) * scale  # Integrate over wavelength to get luminosity

    return Lbol


def get_bolometricGR(R_in, T_p, R_ratio, a, incl):
    """
    Calculate the bolometric luminosity.

    Parameters:
    R_in_star (float, int, np.float64, or np.ndarray): Inner radius times sqrt(cos i) in km, where i is the inclination.
    T_p (float, int, np.float64, or np.ndarray): Peak temperature in keV.
    R_ratio (float, int, np.float64, or np.ndarray): Ratio of outer to inner radius.

    Returns:
    np.ndarray: Bolometric luminosity in erg/s.
    """
    # Ensure inputs are numpy arrays
    if isinstance(R_in, (float, int, np.float64)):
        R_in_star, T_p, R_ratio = np.array([R_in]), np.array([T_p]), np.array([R_ratio])
    elif isinstance(R_in, list):
        R_in_star, T_p, R_ratio = np.array(R_in), np.array(T_p), np.array(R_ratio)

    size = len(R_in)
    Lbol = np.zeros(size)
    E_range = np.logspace(-3.3, 0.7, 100)  # Energy range in keV
    A_range = (E_range * u.keV).to('AA', equivalencies=u.spectral()).value  # Convert energy to Angstrom

    for i in range(size):
        # Calculate spectrum for each energy
        ct_range_kev = kerrSED_sub(E_range, R_in[i], T_p[i], R_ratio[i], a[i], incl[i], 100)
        EF_E = ((ct_range_kev * E_range ** 2) * (u.keV ** 2 / (u.s * u.cm ** 2 * u.keV))).to('erg/(cm**2 s)')
        scale = 4 * np.pi * (100 * u.Mpc.to('cm')) ** 2  # Scaling factor for luminosity distance
        F_wl = (EF_E / A_range).value

        Lbol[i] = np.trapz(np.flip(F_wl), np.flip(A_range)) * scale  # Integrate over wavelength to get luminosity
    return Lbol


def get_Edd_ratio(L_bol, R_in,  i_min=0, i_max=80, a_min=0, a_max=0.998, i_dist=None, a_dist=None):

    Edd_ratio = []
    for j in range(len(L_bol)):
        if a_dist is None:
            M_bh_j = get_Mbh(R_in, i_min=i_min, i_max=i_max, a_min=a_min, a_max=a_max)
            L_bol_j = np.full(len(M_bh_j), L_bol[j])
            Edd_ratio.append( L_bol_j / (1.26e38*M_bh_j))
        else:
            M_bh_j = get_Mbh(R_in, i_dist=i_dist, a_dist=a_dist)
            L_bol_j = np.full(len(M_bh_j), L_bol[j])
            Edd_ratio.append(L_bol_j / (1.26e38 * M_bh_j))

    Edd_ratio = np.array(Edd_ratio).flatten()
    return Edd_ratio

'''
def get_bolometricGR(R_in_star, T_p, R_ratio, a, incl):
    """
    Calculate the bolometric luminosity.

    Parameters:
    R_in_star (float, int, np.float64, or np.ndarray): Inner radius times sqrt(cos i) in km, where i is the inclination.
    T_p (float, int, np.float64, or np.ndarray): Peak temperature in keV.
    R_ratio (float, int, np.float64, or np.ndarray): Ratio of outer to inner radius.


    Returns:
    np.ndarray: Bolometric luminosity in erg/s.
    """
    # Ensure inputs are numpy arrays
    if isinstance(R_in_star, (float, int, np.float64)):
        R_in_star, T_p, R_ratio, a, incl = np.array([R_in_star]), np.array([T_p]), np.array([R_ratio]), np.array([a]), np.array([incl])
    elif isinstance(R_in_star, list):
        R_in_star, T_p, R_ratio, a, incl = np.array(R_in_star), np.array(T_p), np.array(R_ratio), np.array(a), np.array(i)

    size = len(R_in_star)
    Lbol = np.zeros(size)
    E_range = np.logspace(-3.3, 0.7, 200)  # Energy range in keV
    A_range = (E_range * u.keV).to('AA', equivalencies=u.spectral()).value  # Convert energy to Angstrom

    for i in range(size):
        # Calculate spectrum for each energy
        ct_range_kev = kerrSED_sub(E_range, R_in_star[i], T_p[i], R_ratio[i], a[i], incl[i], 100)
        EF_E = ((ct_range_kev * E_range ** 2) * (u.keV ** 2 / (u.s * u.cm ** 2 * u.keV))).to('erg/(cm**2 s)')
        scale = 4 * np.pi * (100 * u.Mpc.to('cm')) ** 2  # Scaling factor for luminosity distance
        F_wl = (EF_E / A_range).value

        Lbol[i] = np.trapz(np.flip(F_wl), np.flip(A_range)) * scale  # Integrate over wavelength to get luminosity

    return Lbol
'''
def get_Lx(R_in_star, T_p, R_ratio, Emin=0.3, Emax=10.0):
    """
    Calculate the X-ray luminosity between Emin and Emax keV.

    Parameters:
    R_in_star (float, int, np.float64, or np.ndarray): Inner radius in km.
    T_p (float, int, np.float64, or np.ndarray): Peak temperature in keV.
    R_ratio (float, int, np.float64, or np.ndarray): Ratio of outer to inner radius.
    Emin (float): Minimum energy in keV.
    Emax (float): Maximum energy in keV.

    Returns:
    np.ndarray: X-ray luminosity in erg/s.
    """
    if isinstance(R_in_star, (float, int, np.float64)):
        R_in_star, T_p, R_ratio = np.array([R_in_star]), np.array([T_p]), np.array([R_ratio])
    elif isinstance(R_in_star, list):
        R_in_star, T_p, R_ratio = np.array(R_in_star), np.array(T_p), np.array(R_ratio)

    size = len(R_in_star)
    Lx = np.zeros(size)
    E_range = np.logspace(np.log10(Emin), np.log10(Emax), 50)  # Energy range in keV
    A_range = (E_range * u.keV).to('AA', equivalencies=u.spectral()).value  # Convert energy to Angstrom

    for i in range(size):
        # Calculate spectrum for each energy
        ct_range_kev = [diskSED_sub(E_range[j], R_in_star[i], T_p[i], R_ratio[i], 100) for j in range(len(E_range))]
        EF_E = ((ct_range_kev * E_range ** 2) * (u.keV ** 2 / (u.s * u.cm ** 2 * u.keV))).to('erg/(cm**2 s)')
        scale = 4 * np.pi * (100 * u.Mpc.to('cm')) ** 2  # Scaling factor for luminosity distance
        F_wl = (EF_E / A_range).value
        Lx[i] = np.trapz(np.flip(F_wl), np.flip(A_range)) * scale  # Integrate over wavelength to get luminosity

    return Lx

def get_Mbh(R_in, i_min=0, i_max=80, a_min=0, a_max=0.998, a_dist=None, i_dist=None):
    """
    Calculate the black hole mass given the inner radius and inclination.

    Parameters:
    R_in (float, int, np.float64, or np.ndarray): Inner radius times sqrt(cos i) in km, where i is the inclination.
    i_min (float): Minimum inclination in degrees.
    i_max (float): Maximum inclination in degrees.
    a_min (float): Minimum dimensionless spin parameter.
    a_max (float): Maximum dimensionless spin parameter.
    a_dist (np.ndarray, optional): Custom distribution for spin parameter 'a'.
    i_dist (np.ndarray, optional): Custom distribution for inclination 'i'.

    Returns:
    np.ndarray: Black hole mass in solar masses.
    """
    if isinstance(R_in, (float, int, np.float64)):
        R_in = np.array([R_in])
    elif isinstance(R_in, list):
        R_in = np.array(R_in)

    # Case 1: Both a_dist and i_dist are None (default distributions)
    if (a_dist is None) and (i_dist is None):
        cos_i_list = np.linspace(np.cos(np.radians(i_min)), np.cos(np.radians(i_max)), 500)
        i_factor = np.sqrt(cos_i_list)
        a_list = np.linspace(a_min, a_max, 500)
        f_a_list = f_a(a_list)
    
    # Case 2: Only a_dist is provided
    elif (a_dist is not None) and (i_dist is None):
        i_factor = np.sqrt(np.random.uniform(np.cos(np.radians(i_min)), np.cos(np.radians(i_max)), len(a_dist)))
        a_list = a_dist
        f_a_list = f_a(a_list)
    
    # Case 3: Only i_dist is provided
    elif (a_dist is None) and (i_dist is not None):
        i_factor = np.sqrt(np.cos(np.radians(i_dist)))
        a_list = np.random.uniform(a_min, a_max, len(i_dist))
        f_a_list = f_a(a_list)
    
    # Case 4: Both a_dist and i_dist are provided
    else:
        if len(a_dist) != len(i_dist):
            raise ValueError("a_dist and i_dist must have the same length when both are provided.")
        i_factor = np.sqrt(np.cos(np.radians(i_dist)))
        a_list = a_dist
        f_a_list = f_a(a_list)
    
    # Calculate M_BH
    Mbh = (R_in[:, np.newaxis] * u.km * ((c ** 2) / (i_factor * G * f_a_list))).to('solMass').value

    return Mbh.flatten()


def get_Mbh_stacking(R_in_list, i_min=0, i_max=80, a_min=0, a_max=0.998, a_dist=None, i_dist=None):
    """
    Calculate the black hole mass by stacking multiple inner radius posteriors.

    Parameters:
    R_in_list (list or np.ndarray): List or array of inner radius values times sqrt(cos i) in km, where i is the inclination.
    i_min (float): Minimum inclination in degrees.
    i_max (float): Maximum inclination in degrees.
    a_min (float): Minimum dimensionless spin parameter.
    a_max (float): Maximum dimensionless spin parameter.
    a_dist (np.ndarray, optional): Custom distribution for spin parameter 'a'.
    i_dist (np.ndarray, optional): Custom distribution for inclination 'i'.

    Returns:
    np.ndarray: Black hole mass in solar masses.
    """
    # Ensure R_in_list is an array of arrays
    if isinstance(R_in_list, (float, int, np.float64)):
        raise ValueError("R_in_list must be a list or array of inner radius arrays.")
    elif isinstance(R_in_list, list):
        R_in_list = np.array(R_in_list, dtype=object)

    # Find the minimum size among all input lists
    sizes = np.array([len(r) for r in R_in_list])
    min_size = np.min(sizes)

    # Truncate all inner radius arrays to match the smallest size
    new_posterior = np.concatenate([r[:min_size] for r in R_in_list])

    # Handle inclination and spin distributions
    if (a_dist is None) and (i_dist is None):
        # Case 1: Default distributions for both
        cos_i_list = np.linspace(np.cos(np.radians(i_min)), np.cos(np.radians(i_max)), 50)
        i_factor = np.sqrt(cos_i_list)
        a_list = np.linspace(a_min, a_max, 50)
        f_a_list = f_a(a_list)

    elif (a_dist is not None) and (i_dist is None):
        # Case 2: a_dist provided, i_dist default
        i_factor = np.sqrt(np.random.uniform(np.cos(np.radians(i_min)), np.cos(np.radians(i_max)), len(a_dist)))
        a_list = a_dist
        f_a_list = f_a(a_list)

    elif (a_dist is None) and (i_dist is not None):
        # Case 3: i_dist provided, a_dist default
        i_factor = np.sqrt(np.cos(np.radians(i_dist)))
        a_list = np.random.uniform(a_min, a_max, len(i_dist))
        f_a_list = f_a(a_list)

    else:
        # Case 4: Both a_dist and i_dist are provided
        if len(a_dist) != len(i_dist):
            raise ValueError("a_dist and i_dist must have the same length when both are provided.")
        i_factor = np.sqrt(np.cos(np.radians(i_dist)))
        a_list = a_dist
        f_a_list = f_a(a_list)

    # Calculate black hole mass
    Mbh = (new_posterior[:, np.newaxis] * u.km * ((c ** 2) / (i_factor * G * f_a_list))).to('solMass').value

    return Mbh.flatten()



def f_a(a):
    """
    Calculate the function f(a), which is the defined as the ration betwen the ISCO (R_isco) and a gravitational
    radii, and is a function of the spin parameter: f(0) = 6, f(1) = 1.

    Parameters:
    a (float or np.ndarray): Dimensionless spin parameter.

    Returns:
    float or np.ndarray: Value of the function f(a).
    """
    Z_1 = 1 + (1 - a ** 2) ** (1 / 3) * ((1 + a) ** (1 / 3) + (1 - a) ** (1 / 3))
    Z_2 = np.sqrt(3 * a ** 2 + Z_1 ** 2)
    return (3 + Z_2 - np.sign(a) * np.sqrt((3 - Z_1) * (3 + Z_1 + 2 * Z_2)))

def get_T_profile(R_in, T_p, R_ratio):
    """
    Calculate the temperature profile of the disk.

    Parameters:
    R_in (float, int, np.float64, or np.ndarray): Inner radius in km.
    T_p (float, int, np.float64, or np.ndarray): Peak temperature in keV.
    R_ratio (float, int, np.float64, or np.ndarray): Ratio of outer to inner radius.

    Returns:
    tuple: Radii (divided by sqrt(cos i)) and corresponding temperatures.
    """
    if isinstance(R_in, (float, int, np.float64)):
        R_in, T_p, R_ratio = np.array([R_in]), np.array([T_p]), np.array([R_ratio])
    elif isinstance(R_in, list):
        R_in, T_p, R_ratio = np.array(R_in), np.array(T_p), np.array(R_ratio)

    r = np.logspace(np.log10(np.median(R_in)), np.log10(np.median(R_in * R_ratio)), 1000)
    x = r / np.median(R_in)
    size = len(T_p)
    T_r = np.zeros((size, 1000))
    for i in range(size):
        T_r[i, :] = ((2.05 * T_p[i]) ** 4 * ((1. - 1. / np.sqrt(x)) / x ** 3)) ** 0.25

    return r, T_r

def get_Rout_Rg(R_ratio, a_min=0, a_max=0.998, a_dist=None):
    """
    Calculate the outer radius in gravitational radii.

    Parameters:
    R_ratio (float, int, np.float64, or np.ndarray): Ratio of outer to inner radius.
    a_min (float): Minimum spin parameter.
    a_max (float): Maximum spin parameter.

    Returns:
    np.ndarray: Outer radius in gravitational radii.
    """
    if isinstance(R_ratio, (float, int, np.float64)):
        R_ratio = np.array([R_ratio])
    elif isinstance(R_ratio, list):
        R_ratio = np.array(R_ratio)

    size = len(R_ratio)
    if a_dist is None:
        a_list = np.linspace(a_min, a_max, 100)
        Rout_Rg = np.zeros((size, 100))
        f_a_list = f_a(a_list)
        for i in range(size):
            Rout_Rg[i, :] = R_ratio[i] * f_a_list

    else:
        a_list = a_dist
        f_a_list = f_a(a_list)
        Rout_Rg = R_ratio * f_a_list


    return Rout_Rg.flatten()

def get_R_circ_Rg(Mbh, beta=1, M_star=1):
    from astropy.constants import c, G, R_sun, M_sun
    if isinstance(Mbh, (float, int, np.float64)):
        Mbh = np.array([Mbh])

    if isinstance(beta, (float, int, np.float64)) and isinstance(M_star, (float, int, np.float64)):
        R_circ = ((2*R_sun*c**2)/(beta*G*M_sun))*(Mbh)**(-2./3.)*M_star**(7./15.)

    '''
    elif (beta == 'beta_pdf') and isinstance(M_star, (float, int)):

        xis = np.random.random(500)
        betas = _icdf_beta(xis)
        R_circ = np.multiply.outer(((2*R_sun*c**2)/(betas*G*M_sun))*M_star**(7./15.) , Mbh**(-2./3.))
        R_circ = R_circ.flatten()
    '''
    return R_circ



def _icdf_beta(xi, beta_max=4):
    """
    Generates penetration factor (beta) distribution as PDF(beta) \propto beta^-2

    Parameters:
    xii (float): random number between 0.0 and 1.0

    Returns:
    beta (float): penetration factor
    """
    C = 1 - (1/beta_max)
    betas = 1 / (1 - (C*xi))

    return betas