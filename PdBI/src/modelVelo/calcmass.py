'''

calculate dynamical mass


Last modified: 21 June 16

History:
21 June 2016
    - add func. to calc dyn mass using disk model (Neri+03)
17 May 2016
    - get mass using quanties reported in paper
13 May 2016
    - extracted from modelGradient.py

'''
import numpy as np

kpc_to_m = 3.086e+19
b_arcsec = 1.8       # from C06
a_arcsec = 3.25
i_rad = np.arccos(b_arcsec/a_arcsec)
i_deg = i_rad * 180./np.pi
Msun = 1.989e30      # kg


def M_encl(R_kpc, v):
    """ calc mass enclose in spherical model, as in Solomon+97

    Parameters
    ----------
    v: float
        in km/s

    Returns
    -------
    M10: float
        same as M = 2.32e5 * V_km/s^2 * R_kpc
        in unit of 1e10 Msun
        """
    from astropy.constants import G

    R = kpc_to_m * R_kpc
    v2 = (v * 1e3)**2
    v2oG = v2 / G.value
    v2oGoMsun = v2oG / Msun
    M = R * v2oGoMsun           # per Msun
    M10 = M / 1e10              # in units of x 1e10 Msun
    return M10

# calc. mass within some physical dist from line center
R_in_kpc = 6.22
V_rot = 345
# no correction for incl. angle
print("Mass enclosed within {:.2f} kpc: {:.2f} x 1e+10 Msun".format(R_in_kpc, M_encl(R_in_kpc, V_rot)))
# corrected
print M_encl(R_in_kpc, V_rot)/(np.sin(i_rad)**2)


# use line peak separation
R_in_kpc = 6.22
V_rot = 400/2.
# no correction for incl. angle
print("Mass enclosed within {:.2f} kpc: {:.2f} x 1e+10 Msun".format(R_in_kpc, M_encl(R_in_kpc, V_rot)))

# with incl. correction
print M_encl(R_in_kpc, V_rot)/(np.sin(i_rad)**2)



def Mdyn_virial(R_kpc, v_kms):
    """
        assuming virialized "isotropic virial estimator" (Spitzer 1987),
            M = 5 * sigma^2 * R/G,
            --> the difference in prefactor as in spherical

        the scaling factor appropriate for a rotating disk at an average inclination is a factor of ∼1.5 smaller (Bothwell et al. 2010), i.e. if the above formula were applied to a disk galaxy, it would somewhat overestimate the dynamical mass.

        v_kms: FWHM in km/s

    Returns
    -------
    M10: float
        in unit of 1e10 Msun
    """
    M = 2.8e5 * (v_kms)**2 * R
    M10 = M / 1e10
    return M10


def Mdyn_disk(R_kpc, FWHM_kms):
    """ Calc. dyn. mass assuming dsk model: M sin^2 = 4e4 * FWHM^2 * R
        (Neri+03) "global rotating disk estimator"
        In a merger model, Mdyn would be a factor of 2 bigger (Genzel+03)

        the numerical constant incoporates a factor of 2.4 between FWHM and the product of rot. velo and sin i (estimated from model disks taking into acc. local line broadening, beam & spectral smearing)
            i.e. for M = \gamma * R * FWHM^2 / G; \gamma = 1 for spherical & 0.3 here

    Returns
    -------
    M10: float
        in unit of 1e10 Msun
    """
    M = 4e4 * (FWHM_kms)**2 * R_kpc
    M10 = M / 1e10
    return M10

