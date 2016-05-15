'''

calculate dynamical mass


Last modified: 13 May 16

History:
13 May 2016
    - extracted from modelGradient.py

'''
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
R_in_kpc = xdata.max()
# print("Mass enclosed within {:.2f} kpc: {:.2f} x 1e+10 Msun".format(R_in_kpc, M_encl(R_in_kpc, pfit[0]/np.sin(i_rad))))

