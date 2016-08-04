'''

Plot SED with MIPS extrapolation to 24 um and MBB_EMCEE fit
+ magphys and de-lensed photometry

Last Modified: Aug 4 2016

History
-------
Aug 4 2016:
    - copied from SEDplot.py
    - added MAGPHYS best-fit

Note
----
- Twin axes work perfectly, even on loglog
- hard coded indices to separate out photometry for different sources, so that we can plot them with different styles

'''

import numpy as np
import matplotlib
from astropy.cosmology import WMAP9

matplotlib.rcParams.update({#'figure.figsize': (8, 5),    # inches
                            'font.size': 18,      # points
                            'legend.fontsize': 16,      # points
                            'text.usetex': True,
                            'font.family': "Simplex",
                            'font.weight': 'bold',
                            'lines.linewidth': 1.5,       # points
                            'axes.linewidth': 1.5,       # points
                            'font.family': "simplex"
                            })

import matplotlib.pyplot as plt
from os.path import join

# Input observed photometry
photoFile = 'RXJ1131photometry.dat'

# output save plot file
plotName = 'FullSED_withMagphys'

# MAGPHYs output
MagphysPath = "/Users/admin/Research/RXJ1131/src/SED/Magphys/"
sedfilename = 'daisy.sed'
fitfilename = 'daisy.fit'
filterfilename = 'filters.dat'
inputfile = 'de-lensed-stellarMass.txt'

# ---------------- MAGPHYS -----------------
# read in original
wavelg_um, flux_mJy, flux_err_mJy = np.loadtxt(
    join(MagphysPath, inputfile), skiprows=2, comments="#", unpack=True)


# read in fluxes from model output
fitfile = open(join(MagphysPath, fitfilename))
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
lambda_eff = np.loadtxt(join(MagphysPath, filterfilename), skiprows=1,
                        comments="#", usecols=[1], unpack=True)


# read in best-fit SED (the one that minimizes chi-2)
x = np.loadtxt(join(MagphysPath, sedfilename), skiprows=10)
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
# ------------------- MAGPHYS -------------------


def read_file(filename):
    """
    Read in rows that are not commented out

    Paramters
    ---------
    filename: str

    Return
    ------
    dataTable: astropy.table.table.Table

    """

    import astropy.io.ascii as ascii
    dataTable = ascii.read(filename, comment='^#')
    # data = ascii.read(filename, format='no_header', data_start=0)

    #dataTable = Table.read(file, format='ascii.commented_header')
    if len(dataTable) == 0:
        errstr = "No data read from %s" % filename
        raise IOError(errstr)
    return dataTable

tbl = read_file(photoFile)

# separate foreground and background
# total = [6:17], [25:33], [35], [38]
# host = [3:6], [17:21],[33], [37]
# fg = [0:3], [21:25],[34], [36]
#

idx_total = range(6, 17)
for i in range(25, 33):
    idx_total.append(i)
idx_total.append(35)
idx_total.append(38)

idx_host = range(3, 6)
for i in range(17, 21):
    idx_host.append(i)
idx_host.append(33)
idx_host.append(37)

idx_fg = range(0, 3)
for i in range(21, 26):
    idx_fg.append(i)
idx_fg.append(34)
idx_fg.append(36)

# separate out photometry and those as upper limits
# print tbl.colnames
# total
t = tbl[idx_total]
idx, = np.where(np.isnan(t['Flux Err [mJy]']) == True)
tblUpLimit_total = t[idx]
idx, = np.where(np.isnan(t['Flux Err [mJy]']) == False)
tblWithErr_total = t[idx]

# separate out those that are upper limits
waveUpLimit_total = np.asarray(tblUpLimit_total['Wavelength [micron]'])
fluxUpLimit_total = np.asarray(tblUpLimit_total['Flux Density [mJy]'])

# non-limits
wave_total = np.asarray(tblWithErr_total['Wavelength [micron]'])
flux_total = np.asarray(tblWithErr_total['Flux Density [mJy]'])
fluxErr_total = np.asarray(tblWithErr_total['Flux Err [mJy]'])


# Host:
h = tbl[idx_host]
idx, = np.where(np.isnan(h['Flux Err [mJy]']) == True)
tblUpLimit_host = h[idx]
idx, = np.where(np.isnan(h['Flux Err [mJy]']) == False)
tblWithErr_host = h[idx]

if len(tblUpLimit_host):
    # separate out those that are upper limits
    waveUpLimit_host = np.asarray(tblUpLimit_host['Wavelength [micron]'])
    fluxUpLimit_host = np.asarray(tblUpLimit_host['Flux Density [mJy]'])
else:
    print(" No upper limits ")

# non-limits
wave_host = np.asarray(tblWithErr_host['Wavelength [micron]'])
flux_host = np.asarray(tblWithErr_host['Flux Density [mJy]'])
fluxErr_host = np.asarray(tblWithErr_host['Flux Err [mJy]'])


# Foreground
fg = tbl[idx_fg]
idx, = np.where(np.isnan(fg['Flux Err [mJy]']) == True)
tblUpLimit_fg = fg[idx]
idx, = np.where(np.isnan(fg['Flux Err [mJy]']) == False)
tblWithErr_fg = fg[idx]

if len(tblUpLimit_fg):
    # separate out those that are upper limits
    waveUpLimit_fg = np.asarray(tblUpLimit_fg['Wavelength [micron]'])
    fluxUpLimit_fg = np.asarray(tblUpLimit_fg['Flux Density [mJy]'])
else:
    print(" No upper limits ")

# non-limits
wave_fg = np.asarray(tblWithErr_fg['Wavelength [micron]'])
flux_fg = np.asarray(tblWithErr_fg['Flux Density [mJy]'])
fluxErr_fg = np.asarray(tblWithErr_fg['Flux Err [mJy]'])


# --------------- Fit slope -------------------------
def powerf(x, prefact, index):
    '''
    parameters
    ----------
    x: np.array
        wavelength in micron
    prefact: float
    index: float
    '''
    # convert wavelength to freq, s.t. we can compare index with literature
    from astropy.constants import c
    x_hz = c.value / (x * 1e-6)
    x_ghz = x_hz / 1.e9
    return prefact * x_ghz**index

def MC_simulation(ydata, yerr, func, p0, xdata, nsam=500, CI=1., debug=False):
    '''
    monte-carlo: generates random data points starting from the given data plus a random variation based on the systematic error, take into acct systematic uncertainties

    Parameters
    ----------
    CI: float
        1 <=> 68.3%

    '''

    ps = []
    for i in range(nsam):
        randomDelta = np.random.normal(0., yerr, size=len(ydata))
        randomdataY = np.array(ydata) + randomDelta

        if debug:
            # verify it's drawing from distribution we expect
            randomDelta = np.random.normal(0., yerr, size=1000)
            count, bins, ignored = plt.hist(randomDelta, 30, normed=True)
            plt.plot(bins, 1/(yerr * np.sqrt(2 * np.pi)) *
                     np.exp(-(bins)**2 / (2 * yerr**2)), linewidth=2, color='r')
            plt.show()
        randomfit, randomcov = optimize.curve_fit(func,
                                                  xdata, randomdataY, p0=p0)
        ps.append(randomfit)

    ps = np.array(ps)
    mean_pfit = np.median(ps, 0)
    Nsigma = CI
    err_pfit_MC = Nsigma * np.std(ps, 0)
    return ps, mean_pfit, err_pfit_MC


def extrapolate_24_MC(ps):
    '''
    get a distribution of fluxes at 24 based on the MC fit parameters
    '''

    out_24 = []
    for i in ps:
        out_24.append(powerf(24., *i))
    return np.array(out_24)


ppower = [1.0, -1.]
import scipy.optimize as optimize
# get initial parameters for MC using curve_fit
p_fg, pcov_fg = optimize.curve_fit(powerf, wave_fg[3:7],
                                   flux_fg[3:7], p0=[0.1, -1.9],
                                   sigma=fluxErr_fg[3:7], absolute_sigma=True)
# perr_fg = np.sqrt(np.diag(pcov_fg))
p_host, pcov_host = optimize.curve_fit(powerf, wave_host[3:7], flux_host[3:7],
                                       sigma=fluxErr_host[3:7], p0=ppower, absolute_sigma=True)
# perr_host = np.sqrt(np.diag(pcov_host))
fg_MC, fg_MC_fit, fg_err_MC = MC_simulation(flux_fg[3:7], fluxErr_fg[3:7], powerf, p_fg, wave_fg[3:7])
print "\n MC method :"
print "pfit = ", fg_MC_fit
print "perr = ", fg_err_MC
host_MC, host_MC_fit, host_err_MC = MC_simulation(flux_host[3:7], fluxErr_host[3:7], powerf, p_host, wave_fg[3:7])
print "\n MC method :"
print "pfit = ", host_MC_fit
print "perr = ", host_err_MC
# ----> err similar from curve_fit, but best-fit differ because curve_fit is very sensitive to intial choice of parameter.. therefore use MC best-fit

# spectral index from MC
print("Spectral index alpha for foreground: {:.2f}+/-{:.2f}").format(fg_MC_fit[1], fg_err_MC[1])
print("Spectral index alpha for RXJ1131: {:.2f}+/-{:.2f}").format(host_MC_fit[1], host_err_MC[1])


# decompose 24 um point:
flux_fg_24 = powerf(24, *fg_MC_fit)
flux_host_24 = powerf(24, *host_MC_fit)
# MC uncertainty
fluxerr_fg_24 = np.std(extrapolate_24_MC(fg_MC))
fluxerr_host_24 = np.std(extrapolate_24_MC(host_MC))


print(" Total flux at MIPS {:.1f} um: {:.3f} mJy").format(wave_total[11], flux_total[11])
print(" Host flux \t\t: {:.3f} +/- {:.3f} mJy").format(flux_host_24, fluxerr_host_24)
print(" Fg flux \t\t: {:.3f} +/- {:.3f} mJy").format(flux_fg_24, fluxerr_fg_24)
print(" Sum from extrapolated {:.3f} mJy").format((flux_fg_24 + flux_host_24))

# extrapolate
space = 0.2
xspan = wave_fg[3:7].max() - wave_fg[3:7].min()
xspace = space * xspan
xarray = np.linspace(wave_fg[3:7].min()-space, 24+xspace, 500)

# --------------- Plot Data points ----------------------
def savefigure(figure, f, verbose=True, dirStr='../../Figures/'):
    """
    Save figure in designatied directory and make directory if it doesn't alraedy exist

    Paramters
    ---------
    figure: matplotlib.figure.Figure
        since plt.savefig() will be empty is calling after plt.show(). using figure Object will solve this problem

    f: str
        filename without extension
    dirStr: str
        default to the ../../Figures/
    """
    import os
    if verbose is True:
        if not os.path.isdir(dirStr):
            print "... Making directory {:s}...".format(dirStr)
            os.mkdir(dirStr)
        else:
            print "... Directory {:s} exists...".format(dirStr)
    figure.savefig(dirStr + f + '.pdf', dvi=600, bbox_inches='tight')
    print('... Saved Figure: {:s}').format(dirStr + f + '.eps')

f, ax = plt.subplots(figsize=(10, 8))
f.subplots_adjust(left=0.10, bottom=0.09, top=0.85, right=0.98)

# plot photometry
ax.errorbar(wave_total, flux_total, fluxErr_total, fmt=".r", ms=8,
            ecolor='gray', label='Total', capsize=8, elinewidth=4,
            capthick=1.5)
ax.errorbar(wave_host, flux_host, fluxErr_host, fmt=".b", ms=8, ecolor='gray',
            label='Host', capsize=8, elinewidth=4, capthick=1.5)
ax.errorbar(wave_fg, flux_fg, fluxErr_fg, fmt=".g", ms=8, ecolor='gray',
            label='Foreground', capsize=8, elinewidth=4, capthick=1.5)
ax.errorbar(waveUpLimit_total, fluxUpLimit_total, uplims=True, fmt="vk", ms=8,
            ecolor='gray', capsize=8, elinewidth=4, capthick=1.5)  # label='upper limits'

# Plot the extrapolate flux at 24 um with error bars
plt.errorbar(24., flux_host_24, fluxerr_host_24, fmt=".b", ms=8, ecolor='gray',
             capsize=8,
             elinewidth=4, capthick=1.5)    # label='Extrapolated from fit to decomposed IRAC; Host',

plt.errorbar(24., flux_fg_24, fluxerr_fg_24, fmt='.g', ms=8, ecolor='gray',
             capsize=8, elinewidth=4, capthick=1.5)
             # label='Extrapolated from fit to decomposed IRAC; Foreground',

# plot MBB results
import mbb_emcee
filename_withMIPS = '8pts/thickPower_priorBeta_limT.h5'
filename_noMIPS = '7pts/thickPower_priorBeta_limT.h5'

res = mbb_emcee.mbb_results(h5file=filename_withMIPS)
res_noMIPS = mbb_emcee.mbb_results(h5file=filename_noMIPS)

wave_SED, flux_SED, flux_unc_SED = res.data
p_wave = np.linspace(wave_SED.min() * 0.5, wave_SED.max() * 1.5, 200)

ax.plot(p_wave, res.best_fit_sed(p_wave), 'k--',
        label='MBB fit with 24 um', alpha=1.0, zorder=0.65)
ax.plot(p_wave, res_noMIPS.best_fit_sed(p_wave), '-.',
        color='#152863', label='MBB fit without 24 um',
        alpha=0.9, zorder=0.5)


# plot MAGPHYS results
ax.plot(xx, flux_at, 'k', label='attenuated MAGPHYS fit')
ax.plot(xx, flux_unatt, color='#439EAA', linewidth=1.0, alpha=0.8)

# de-lensed
plt.errorbar(wavelg_um, flux_mJy, flux_err_mJy, color='#335DE5', mec='#335DE5', ecolor='gray', fmt='o', ms=5, capsize=8, elinewidth=3, capthick=1.5, label='De-lensed')

# ------- pretty plot configuration below -------
# W/o de-lenesed;
# ax.set_ylim(0.001, 2500.)
# ax.set_xlim(0.1, 3e5)

# W/ de-lensed & magphys
ax.set_ylim(1e-4, 2500.)
ax.set_xlim(0.05, 3e5)

ax.set_yscale("log")
ax.set_xscale("log")

ax.set_ylabel(r'S$_{\nu}$ [mJy]', fontsize=16)
ax.set_xlabel(r'$\lambda_{\rm obs}$ [$\mu$m]', fontsize=16, fontweight='bold')

led = plt.legend(loc='upper right', fontsize=14, numpoints=1,
                 fancybox=True, borderpad=0.8,
                 handlelength=1.25, labelspacing=0.3)
ax.tick_params(length=14, pad=5)

# get the existing axes limits and convert to an array
_axes1_range = np.array(ax.axis())
axes2_range = _axes1_range.copy()

# set second x-axis, convert limits
import astropy.units as u
axes2_range[0:2] = (_axes1_range[0:2]*u.micron).to(u.GHz, equivalencies=u.spectral()).value

# instantiate, Mpl weird thing: x and y are inverse
# ax_topy = ax.twinx()
ax_topx = ax.twiny()

# set limits for twin axes
ax_topx.set_xlim(axes2_range[:2])
ax_topx.set_xscale('log')
# ax_topy.set_yscale('log')
ax_topx.tick_params(length=14)
ax_topx.set_xlabel(r'$\nu_{\rm obs}$ [GHz]', size=16)

# plt.tight_layout()
plt.show(block=False)   # so it doesn't hold off the next command
response = raw_input('Save fig?: (y/n)')

if response.lower() in ['y', 'yes']:
    savefigure(f, plotName)
else:
    print('...Exiting...')
    plt.close()

