'''

calculate dynamical mass


Last modified: 17 May 16

History:
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

def M_encl(R_kpc, v):
    """ calc mass enclosed """     # v in km/s
    from astropy.constants import G

    R = kpc_to_m * R_kpc
    Msun = 1.989e30      # kg
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

