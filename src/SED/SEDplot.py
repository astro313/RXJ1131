'''

Plot SED

Last Modified: May 13 2016

History
-------
May 15 2016:
    - curve_fit to power_law working
    - but extrapolated flux > MIPS point
May 13 2016:
    - update hardcode part because added HST
    - replace least sq with curve_fit, but curve_fit to powerlaw is not working properly
May 10 2016:
    - plot foreground vs. background with different style, markers
    fit decomposed IRAC points and extracted to "decompose" MIPS
    - put mbb_emcee fit
    - added power law fit

Jan 09 2016:
    - created

Note
----
- Twin axes work perfectly, even on loglog
- hard coded indices to separate out photometry for different sources, so that we can plot them with different styles

'''

import numpy as np
import matplotlib.pyplot as plt
import matplotlib

photoFile = 'RXJ1131photometry.dat'
plotName = 'FullSED'

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
    return prefact * x**index
ppower = [1.0, 1.]
import scipy.optimize as optimize
p_fg, pcov_fg = optimize.curve_fit(powerf, wave_fg[3:7],
                                   flux_fg[3:7], p0=[0.9, 1],
                                   sigma=fluxErr_fg[3:7], absolute_sigma=True, maxfev=10000)
perr_fg = np.sqrt(np.diag(pcov_fg))

p_host, pcov_host = optimize.curve_fit(powerf, wave_host[3:7], flux_host[3:7],
                                       sigma=fluxErr_host[3:7], p0=ppower, absolute_sigma=True)
perr_host = np.sqrt(np.diag(pcov_host))

# decompose 24 um point:
flux_fg_24 = powerf(24, *p_fg)
flux_host_24 = powerf(24, *p_host)
# propagate uncertainty
# sqrt((24^index * delta_amp)**2 + (24^(index-1) * delta_index * index * amp)**2)
fluxerr_host_24 = np.sqrt((24.**p_host[1] * perr_host[0])**2 + (p_host[1] * p_host[0] * 24.**(p_host[1]-1) * perr_host[1])**2)
fluxerr_fg_24 = np.sqrt((24.**p_fg[1] * perr_fg[0])**2 + (p_fg[1] * p_fg[0] * 24.**(p_fg[1]-1) * perr_fg[1])**2)

print(" Total flux at MIPS {:.1f} um: {:.3f} mJy").format(wave_total[11], flux_total[11])
print(" Host flux \t\t: {:.3f} +/- {:.3f} mJy").format(flux_host_24, fluxerr_host_24)
print(" Fg flux \t\t: {:.3f} +/- {:.3f} mJy").format(flux_fg_24, fluxerr_fg_24)
print(" Sum from extrapolated {:.3f} mJy").format((flux_fg_24 + flux_host_24))


# extrapolate
space = 0.2
xspan = wave_fg[3:7].max() - wave_fg[3:7].min()
xspace = space * xspan
xarray = np.linspace(wave_fg[3:7].min()-space, 24+xspace, 500)

plt.figure()
plt.errorbar(wave_fg[3:7], flux_fg[3:7], fluxErr_fg[3:7], fmt='.g', ms=8, ecolor='gray', capsize=8, elinewidth=4, capthick=1.5)
plt.errorbar(wave_host[3:7], flux_host[3:7], fluxErr_host[3:7], fmt='.b', ms=8,
             ecolor='gray',
             capsize=8,
             elinewidth=4, capthick=1.5)

plt.plot(xarray, powerf(xarray, *p_fg), 'r--', label='fit slope to FG IRAC')
plt.plot(xarray, powerf(xarray, *p_host), 'b--', label='fit slope to host IRAC')

# Plot the extrapolate flux at 24 um with error bars
plt.errorbar(24., flux_host_24, fluxerr_host_24, fmt=".b", ms=8, ecolor='gray',
             label='Extrapolated from fit to decomposed IRAC; Host', capsize=8,
             elinewidth=4, capthick=1.5)
plt.errorbar(24., flux_fg_24, fluxerr_fg_24, fmt='.g', ms=8, ecolor='gray',
             label='Extrapolated from fit to decomposed IRAC; Foreground',
             capsize=8, elinewidth=4, capthick=1.5)
plt.xscale('log')
plt.yscale('log')
plt.show(block=False)

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
    figure.savefig(dirStr + f + '.eps', dvi=600, box_inches='tight')
    print('... Saved Figure: {:s}').format(dirStr + f + '.eps')


font = {'family': 'Arial Narrow',
        'weight': 'bold',
        'size': 15}
matplotlib.rc('font', **font)


f, ax = plt.subplots()
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
            ecolor='gray', label='upperlimits', capsize=8, elinewidth=4,
            capthick=1.5)

# plot least sq result, fit to decomposed IRAC points
# extrapolate
space = 0.2
xspan = wave_fg[3:7].max() - wave_fg[3:7].min()
xspace = space * xspan
xarray = np.linspace(wave_fg.min()-xspace, 24+xspace, 500)

plt.plot(xarray, powerf(xarray, p_fg), 'r--', label='fit slope to FG IRAC')
plt.plot(xarray, powerf(xarray, p_host), 'b--', label='fit slope to host IRAC')

# Plot the extrapolate flux at 24 um with error bars
plt.errorbar(24., flux_host_24, fluxerr_host_24, fmt=".b", ms=8, ecolor='gray',
             label='Extrapolated from fit to decomposed IRAC; Host', capsize=8,
             elinewidth=4, capthick=1.5)
plt.errorbar(24., flux_fg_24, fluxerr_fg_24, fmt='.g', ms=8, ecolor='gray',
             label='Extrapolated from fit to decomposed IRAC; Foreground',
             capsize=8, elinewidth=4, capthick=1.5)

# plot MBB results
import mbb_emcee
filename_withMIPS = '8pts/thickPower_priorBeta_limT.h5'
filename_noMIPS = '7pts/thickPower_priorBeta_limT.h5'

res = mbb_emcee.mbb_results(h5file=filename_withMIPS)
res_noMIPS = mbb_emcee.mbb_results(h5file=filename_noMIPS)

wave_SED, flux_SED, flux_unc_SED = res.data
p_wave = np.linspace(wave_SED.min() * 0.5, wave_SED.max() * 1.5, 200)

ax.plot(p_wave, res.best_fit_sed(p_wave), 'b--', lw=2.5,
        label='MBB fit with 24 um', alpha=0.45)
ax.plot(p_wave, res_noMIPS.best_fit_sed(p_wave), '-',
        color='c', lw=2.5, label='MBB fit w/o 24 um', alpha=0.45)


# ------- pretty plot configuration below -------
ax.set_ylim(0.1, 1700.)
ax.set_yscale("log")
ax.set_xscale("log")

ax.set_ylabel(r'$\rm S_{\nu}$ [mJy]', fontsize=16)
ax.set_xlabel(r'$\lambda_{\rm obs}\ [\mu$m]', fontsize=16, fontweight='bold')

led = plt.legend(loc='best', fontsize=15, numpoints=1,
                 fancybox=True, borderpad=0.5,
                 handlelength=2, labelspacing=0.1)
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

plt.tight_layout()
plt.show(block=False)   # so it doesn't hold off the next command
response = raw_input('Save fig?: (y/n)')

if response.lower() in ['y', 'yes']:
    savefigure(f, plotName)
else:
    print('...Exiting...')
