'''

If we want to plot on the mbb_fit?... but the data points on that was not lensing corrected...

Maybe want to make a bigger plot to also show the delensed?

'''

import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from astropy.cosmology import WMAP9

z = 0.6537
# (1+z) / (4 * pi * DLmpc_to_cm^2) * L = S
dl = WMAP9.luminosity_distance(z).value  # in Mpc
Mpc_to_cm = 3.086e+24
prefac = 1./(4*np.pi*Mpc_to_cm**2 * dl**2) # * ( 1+z)

mpl.rcdefaults()

x = np.loadtxt('daisy.sed', skiprows=10)
log_lamb_AA, atten_lambLum, unatten_lambLum = x[:, 0], x[:, 1], x[:, 2]

x = 10**log_lamb_AA     # wavelength in AA
xx = x / 1.e+4  # log_lamb_AA in microns

Lsun_to_cgs = 3.839e33

# ergs/s
y_at = x * 10**(atten_lambLum) * Lsun_to_cgs      # total attenuated SED
y_un = x * 10**(unatten_lambLum) * Lsun_to_cgs    # unattenuated SED


# need to divide by Area^2 and Hz^-1?
# ergs/s/cm
y_at_ = y_at / (x*1e-8)
y_un_ = y_un / (x*1e-8)

# from 1/cm to /hz
y_at_ *= (x*1e-8)**2 /3.e10
y_un_ *= (x*1e-8)**2 /3.e10


flux_at = y_at_ * prefac/ 1.e-26    # to milliJy
flux_unatt = y_un_ * prefac/ 1.e-26

plt.loglog(xx, flux_at)
plt.ylim(1e-3, 100)
plt.show()