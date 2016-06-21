
'''
get unc. for decomposed Mgas for RXJ1131 and companion based on flux in Chan 126-130

Note:
    i.e. Mgas here is not the final total Mgas for RXJ1131

'''
import numpy as np

# In chan 126-130, lensed
Mtot = 3.034e10      # Msun

# RXJ
mu1 = 7.2
mu1_e = 5.6
S1 = 4.818        # mJy
S1_e = 3.456

# Companion
mu2 = 6.7
mu2_e = 2.5
S2 = 5.955
S2_e = 3.565

def M_gas(mu1, mu2, s1, s2, Mtot):
    '''
    separate Mgas for companion and RXJ1131

    Parameters
    ----------
    mu1, mu2:
        mag. factor for RXJ1131 or companion
    s1, s2:
        intrisnic flux for RXJ1131 or companion from model

    Returns
    -------
    Mgas in unit of Msun for RXJ1131 or companion, dep. on order when input parameters

    '''
    totalFlux = mu1 * s1 + mu1 * s2        # order does not matter
    s1Flux = mu1 * s1                      # order matters
    return Mtot * s1Flux / totalFlux


# get distribution on Mgas based on unc. in mu and intrinsic flux
niters = 1e3
mu1_iters = abs(np.random.normal(loc=mu1, scale=mu1_e, size=niters))
mu2_iters = abs(np.random.normal(loc=mu2, scale=mu2_e, size=niters))
S1_iters = abs(np.random.normal(loc=S1, scale=S1_e, size=niters))
S2_iters = abs(np.random.normal(loc=S2, scale=S2_e, size=niters))

# M_1
mass_gas_iters = M_gas(mu1_iters, mu2_iters, S1_iters, S2_iters, Mtot)
M_gas_mean = np.mean(mass_gas_iters)
M_gas_median = np.median(mass_gas_iters)

CI = 1.        # 1 <=> 68.3 %
Nsigma = CI
e_M_gas = Nsigma * np.std(mass_gas_iters)
print " RXJ"
print "corrected for lensing"
print M_gas_mean/mu1/1e10, M_gas_median/mu1/1e10, e_M_gas/mu1/1e10

# M_2
mass_gas_iters = M_gas(mu2_iters, mu1_iters, S2_iters, S1_iters, Mtot)
M_gas_mean = np.mean(mass_gas_iters)
M_gas_median = np.median(mass_gas_iters)
e_M_gas = Nsigma * np.std(mass_gas_iters)
print " Companion "
print M_gas_mean/mu2/1e10, M_gas_median/mu2/1e10, e_M_gas/mu2/1e10


