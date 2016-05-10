'''

Plot SED

Last Modified: May 10 2016


TODO
- put mbb_emcee fit

History
-------
May 10 2016:
    plot foreground vs. background with different style, markers
    fit decomposed IRAC points and extracted to "decompose" MIPS
Jan 09 2016:
    created

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
# total = [:11], [19:27], [29], [32]
# host = [11:15],[27], [31]
# fg = [15:19],[28], [30]
#

idx_total = range(0, 11)
for i in range(19, 27):
    idx_total.append(i)
idx_total.append(29)
idx_total.append(32)

idx_host = range(11, 15)
idx_host.append(27)
idx_host.append(31)

idx_fg = range(15, 19)
idx_fg.append(28)
idx_fg.append(30)

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
def linef(x, p):
    p0, p1 = p
    return p0 * x + p1


def fit_function(p0, datax, datay, yerr, function, **kwargs):

    import scipy.optimize as optimize
    errfunc = lambda p, x, y, err: (function(x, p) - y) / err
    pfit, pcov, infodict, errmsg, success = optimize.leastsq(errfunc, p0, args=(datax, datay, yerr), full_output=1)

    if (len(datay) > len(p0)) and pcov is not None:
        reducechi2_sq = (errfunc(pfit, datax, datay, yerr)**2).sum() / (len(datay) - len(p0))
        pcov = pcov * reducechi2_sq
    else:
        pcov = inf

    error = []
    for i in range(len(pfit)):
        try:
            error.append(np.absolute(pcov[i][i])**0.5)
        except:
            error.append(0.00)
    pfit_leastsq = pfit
    perr_leastsq = np.array(error)
    return pfit_leastsq, perr_leastsq

# IRAC points
p0 = [0.5, 5.]
p_fg, perr_fg = fit_function(p0, wave_fg[:4],
                             flux_fg[:4], fluxErr_fg[:4], linef)
p_host, perr_host = fit_function(p0, wave_host[:4], flux_host[:4],
                                 fluxErr_host[:4], linef)

# decompose 24 um point:
flux_fg_24 = linef(24, p_fg)
flux_host_24 = linef(24, p_host)
# propagate uncertainty
fluxerr_host_24 = np.sqrt((24. * perr_host[0])**2 + (perr_host[1])**2)
fluxerr_fg_24 = np.sqrt((24. * perr_fg[0])**2 + (perr_fg[1])**2)

print(" Total flux at MIPS {:.1f} um: {:.3f} mJy").format(wave_total[11], flux_total[11])
print(" Host flux \t\t: {:.3f} +/- {:.3f} mJy").format(flux_host_24, fluxerr_host_24)
print(" Fg flux \t\t: {:.3f} +/- {:.3f} mJy").format(flux_fg_24, fluxerr_fg_24)
print(" Sum from extrapolated {:.3f} mJy").format((flux_fg_24 + flux_host_24))

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
xspan = wave_fg[:4].max() - wave_fg[:4].min()
xspace = space * xspan
xarray = np.linspace(wave_fg.min()-xspace, 24+xspace, 500)

plt.plot(xarray, linef(xarray, p_fg), 'r--', label='fit slope to FG IRAC')
plt.plot(xarray, linef(xarray, p_host), 'b--', label='fit slope to host IRAC')

# Plot the extrapolate flux at 24 um with error bars
plt.errorbar(24., flux_host_24, fluxerr_host_24, fmt=".b", ms=8, ecolor='gray',
             label='Extrapolated from fit to decomposed IRAC; Host', capsize=8,
             elinewidth=4, capthick=1.5)
plt.errorbar(24., flux_fg_24, fluxerr_fg_24, fmt='.g', ms=8, ecolor='gray',
             label='Extrapolated from fit to decomposed IRAC; Foreground',
             capsize=8, elinewidth=4, capthick=1.5)

ax.set_ylim(0.1, 1700.)
ax.set_yscale("log")
ax.set_xscale("log")

ax.set_ylabel(r'$\rm S_{\nu}$ [mJy]', fontsize=16)
ax.set_xlabel(r'$\rm \lambda [\mu$m]', fontsize=16)

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
