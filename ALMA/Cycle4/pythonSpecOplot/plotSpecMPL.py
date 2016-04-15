'''
Overplot spectra at different spatial location for ALMA C4 proposal


History:
--------
Appril 14 2016:
    - edit labels and text annotations
April 11 2016:
    - copied from Zsearch project: ~/Research/zsearch_Conley/science_goal.uidA001_X121_X1fc/combine

Note
----

'''

import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import scipy


def gauss(x, a, b, c):
    '''
    x: array
    a,b,c: float
    '''
    return a * np.exp(-(x - b) ** 2 / (2 * c ** 2))


def gauss_withcont(x, a, b, c, d):
    '''
    x: array
    a,b,c,d: float
    '''
    return a * np.exp(-(x - b) ** 2 / (2 * c ** 2)) + d


def double_gauss_withcont(x, a, b, c, d, e, f, g):
    '''
    x: array
    a,b,c,d,e,f,g: float
    '''
    return a * np.exp(-(x - b)**2 / (2 * c**2)) + d + e * np.exp(-(x - f)**2 / (2 * g**2))


def read_specfile(specFile):
    """read spectrum file

    Parameters
    ----------
    specFile: str
        file name containing 2 or 3 columns
        velo [km/s]/freq [GHz], flux [??], rms [optional]
    Returns
    -------
    x, flux, width, rms [optional]

    """

    a = np.loadtxt(specFile)
    x, flux = a[:, 0], a[:, 1]
    width = abs(x[1] - x[0])

    if a.shape[1] == 2:
        return x, flux, width
    elif a.shape[1] == 3:
        rms = a[:, 2]
        return x, flux, width, rms


def setup_spec_plot(ax, xax, xlabel=False, ylabel=False):
    '''

    ax: axis object
    xax: str
        'velocity' or 'frequency'
    xlabel: boolean
        True to turn xaxis label on
    ylabel: Boolean
        True to turn yaxis label on
    '''

    from matplotlib.ticker import AutoMinorLocator

    # specify a fixed number of minor intervals per major interval, e.g.:
    # minorLocator = AutoMinorLocator()
    minorLocator = AutoMinorLocator()
    ax.xaxis.set_minor_locator(minorLocator)
    ax.yaxis.set_minor_locator(minorLocator)

    if ylabel:
        ax.set_ylabel('Flux density [mJy/Beam]')
    if xax == 'velocity':
        if xlabel:
            ax.set_xlabel('Velocity [km/s]')
    elif xax == 'frequency':
        if xlabel:
            ax.set_xlabel(r'$\nu_{\rm obs}$ [GHz]')

    ax.minorticks_on()
    ax.tick_params(which='both', width=1.2)
    ax.tick_params(which='major', length=10)
    ax.tick_params(which='minor', length=5)


def plot_spec_base(ax, xarray, flux, width, bottom=None, rms=None, ec='yellow', c='yellow', color='black', fill=True, ls='-'):
    """

    Parameters
    ----------
    ax: axis object
    xarray: array
    flux: array
        same size as xarray
    width: float
    bottom: None or float
        if there's continuum, this would be the start y of bar
    rms: array [Optional]
        same size as xarray, same unit as flux
    color: string [Optional]
        line color
    fill: Bool
        fill bar if True
    ls: str
        mpl linestyle
    Returns
    -------

    """

    if fill:
        if rms is None:
            l = ax.bar(xarray, flux, width, bottom=bottom, color=c, ec=ec, align='center')
        else:
            l = ax.bar(xarray, flux, width, yerr=rms, bottom=bottom, color=c, ec=ec, align='center')

    ax.plot(xarray, flux,
            drawstyle='steps-mid',
            linestyle='solid',
            color=color,
            ls=ls,
            lw=1.25)

def ylim(ax, fluxMin, fluxMax, ylowerpad=0.5, yupperpad=1):
    """

    Parameters
    ----------
    ax: axis object
    fluxMin: float
    fluxMax: float
    ylowerpad: float
        to set lower yaxis limit
    yupperpad: float
        to set upper yaxis limit
    Returns
    -------

    """
    ax.set_ylim(fluxMin - ylowerpad, fluxMax + yupperpad)


def put_baseline(ax, ybase=0):
    ax.axhline(y=ybase, linewidth=1, color='k', ls='dashed')


# ==================================================================
#
#

specF1 = 'sup127_155_2ndcln_noCont_sym.spec'
specF2 = 'sup127_155_2ndcln_noCont_blue2.spec'
specF3 = 'sup127_155_2ndcln_noCont_red.spec'
xaxis = 'velocity'
x1, flux1, width1 = read_specfile(specF1)
x2, flux2, width2 = read_specfile(specF2)
x3, flux3, width3 = read_specfile(specF3)

# ===========================================================================
# Plot setup
fig = plt.figure()
fig.subplots_adjust(top=0.8)

ax1 = fig.add_subplot(111)
setup_spec_plot(ax1, xaxis, xlabel=True, ylabel=True)
plot_spec_base(ax1, x1[70:200], flux1[70:200], width1, color='green', fill=False)
plot_spec_base(ax1, x2[70:200], flux2[70:200], width2, color='blue', c=None, fill=False) # , ls='dashed')
plot_spec_base(ax1, x3[70:200], flux3[70:200], width3, color='red', c=None, fill=False) #, ls='dashed')
put_baseline(ax1)

fluxMin = min(flux1.min(), flux2.min(), flux3.min())
fluxMax = max(flux1.max(), flux2.max(), flux3.max())
ylim(ax1, fluxMin, fluxMax)
ax1.set_xlim(x1[200], x1[70])

plt.figtext(0.7, 0.7, "CO(2-1)", fontsize=20, fontweight='bold')
plt.figtext(0.15, 0.7, "Position vs color", fontsize=20, fontweight='bold')

plt.savefig('olayCO21Spectra.eps', bbox_inches='tight')
plt.show()

