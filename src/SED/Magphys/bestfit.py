'''

convert the best-fit (.sed) values into flux [mJy] versus microns.

plot unlensed data + MAGPHYS fit, and residual, make it look like MAGPHYS output plot

Last Modfied: Aug 1 2016

History:
Aug 1 2016:
    - created


'''

import numpy as np
import matplotlib
from astropy.cosmology import WMAP9
import astropy.units as u


# Boiler-plate settings for producing pub-quality figures
# 1 point = 1/72 inch
matplotlib.rcParams.update({'figure.figsize': (8, 5),    # inches
                            'font.size': 12,      # points
                            'legend.fontsize': 16,      # points
                            'lines.linewidth': 1.5,       # points
                            'axes.linewidth': 1.5,       # points
                            'font.family': "simplex"
                            })

import matplotlib.pyplot as plt      # must come after rcParams.update()


sedfilename = 'daisy.sed'
fitfilename = 'daisy.fit'
filterfilename = 'filters.dat'
inputfile = 'de-lensed-stellarMass.txt'
foutname = 'de-lensed-magphys.pdf'

# read in original
wavelg_um, flux_mJy, flux_err_mJy = np.loadtxt(
    inputfile, skiprows=2, comments="#", unpack=True)


# read in fluxes from model output
fitfile = open(fitfilename)
fitinfo = fitfile.readlines()
fitfile.close()

# strip out new lines (\r\n)
for i in range(len(fitinfo)):
    fitinfo[i] = fitinfo[i].strip()

filternames = fitinfo[1].strip("#")
obs_filters = filternames.split()
# observed flux L_nu [Lsun / Hz] through each filter
obs_flux = np.array(fitinfo[2].split(), dtype=float)
obs_flux_err = np.array(fitinfo[3].split(), dtype=float)
obs_predict = np.array(fitinfo[12].split(), dtype=float)

# read in best-fit Chi2 and redshift
bestfitmodel = fitinfo[8].split()
bestfit_chi2 = float(bestfitmodel[2])
z = float(bestfitmodel[3])


# read in effective wavelength at obs. flux
lambda_eff = np.loadtxt(filterfilename, skiprows=1,
                        comments="#", usecols=[1], unpack=True)


# read in best-fit SED (the one that minimizes chi-2)
x = np.loadtxt(sedfilename, skiprows=10)
log_lamb_AA, atten_lambLum, unatten_lambLum = x[:, 0], x[:, 1], x[:, 2]

x = 10**log_lamb_AA     # wavelength in AA
xx = x / 1.e+4  # log_lamb_AA in microns

Lsun_to_cgs = 3.839e33

# ergs/s
y_at = x * 10**(atten_lambLum) * Lsun_to_cgs      # total attenuated SED
y_un = x * 10**(unatten_lambLum) * Lsun_to_cgs    # unattenuated SED

# need to divide by Hz^-1 or cm^-1
# ergs/s/cm
y_at_ = y_at / (x * 1e-8)
y_un_ = y_un / (x * 1e-8)

# dlam / dnu = c/nu^2 = lam^2/c
# L_nu = L_lam * dlam / dnu = L_lam * lam^2 /c
# from 1/cm to /hz
y_at_ *= (x * 1e-8)**2 / 3.e10
y_un_ *= (x * 1e-8)**2 / 3.e10


dl = WMAP9.luminosity_distance(z).value  # in Mpc
Mpc_to_cm = 3.086e+24
prefac = 1. / (4 * np.pi * Mpc_to_cm**2 * dl**2)  # * (1+z)

flux_at = y_at_ * prefac / 1.e-26    # to milliJy
flux_unatt = y_un_ * prefac / 1.e-26


# convert obs flux L_nu to mJy
obs_flux_mJy = obs_flux * prefac / 1.e-26 * Lsun_to_cgs
obs_flux_err_mJy = obs_flux_err * prefac / 1.e-26 * Lsun_to_cgs
obs_predict_mJy = obs_predict * prefac / 1.e-26 * Lsun_to_cgs


# plot
fig1 = plt.figure(1)
frame1 = fig1.add_axes((.1, .3, .8, .6))

plt.loglog(xx, flux_at, 'k', xx, flux_unatt, 'b')
plt.ylim(1e-5, 300)
plt.xlim(0.05, 3e5)
frame1.set_xticklabels([])

# These have been convolved with filter.. don't think we wanna show
## plt.errorbar(lambda_eff, obs_flux_mJy, obs_flux_err_mJy, fmt='r+', label='convolved')
# plt.errorbar(lambda_eff, obs_predict_mJy, obs_flux_err_mJy, fmt='r+', label='Model prediction')

plt.errorbar(wavelg_um, flux_mJy, flux_err_mJy, fmt='r+', markersize='10', label='original data')
plt.ylabel('Flux [mJy]')
frame1.annotate('RXJ 1131-1231 at $z$ = %.4f\n $\chi^2$ = %.2f'
                % (z, bestfit_chi2),
                xy=(.05, .95),       # fraction of the plot
                xycoords='axes fraction',
                ha="left", va="top")
#                bbox=dict(boxstyle="round", fc='1'))

frame2 = fig1.add_axes((.1, .1, .8, .2))
# [lower_error, upper_error]
plt.errorbar(x=lambda_eff, y=(obs_predict_mJy - obs_flux_mJy),
             yerr=obs_flux_err_mJy, fmt="k+", label="Residuals")
plt.hlines(0., 0.05, 3e5, linestyles='dotted')
plt.xscale("log")
plt.ylabel('Residuals')
plt.xlabel('Observed-frame wavelength [$\mu m$]')
frame2.set_xlim(frame1.get_xlim())

from matplotlib.ticker import MaxNLocator
frame2.yaxis.set_major_locator(MaxNLocator(prune='upper'))

plt.show(block=False)   # so it doesn't hold off the next command

response = raw_input('Save fig?: (y/n)')

if response.lower() in ['y', 'yes']:
    plt.savefig(foutname, bbox_inches='tight')
else:
    print('...Exiting...')
    plt.close()

