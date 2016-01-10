'''

Plot points on SED, put MIR and radio data on plot for background source but don't include in fit, maybe use IRAS 60, 100 micron to constrain SED peak
- mark continuum for foreground and background differently

Last Modified: Jan 09 2016


TODO
- put mbb_emcee fit
-


History
-------
Jan 09 2016:
    created


Note
----
- Twin axes work perfectly, even on loglog

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

# separate out photometry and those as upper limits
print tbl.colnames
idx, = np.where(np.isnan(tbl['Flux Err [mJy]']) == True)
tblUpLimit = tbl[idx]
idx, = np.where(np.isnan(tbl['Flux Err [mJy]']) == False)
tblWithErr = tbl[idx]

waveUpLimit = np.asarray(tblUpLimit['Wavelength [micron]'])
fluxUpLimit = np.asarray(tblUpLimit['Flux Density [mJy]'])

wave = np.asarray(tblWithErr['Wavelength [micron]'])
flux = np.asarray(tblWithErr['Flux Density [mJy]'])
fluxErr = np.asarray(tblWithErr['Flux Err [mJy]'])


# Plot
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
#    os.system('rm -rf ' + f + '.png' + ' ' + f + '.eps')
#    figure.savefig(dirStr + f + '.png')
    figure.savefig(dirStr + f + '.eps', dvi=600, box_inches='tight')
    print('... Saved Figure: {:s}').format(dirStr + f + '.eps')


font = {'family': 'Arial Narrow',
        'weight': 'bold',
        'size': 15}
matplotlib.rc('font', **font)

f, ax = plt.subplots()
f.subplots_adjust(left=0.10, bottom=0.09, top=0.85, right=0.98)
ax.errorbar(wave, flux, fluxErr, fmt=".r", ms=8, ecolor='gray',
            label='photo', capsize=8, elinewidth=4, capthick=1.5)
ax.errorbar(waveUpLimit, fluxUpLimit, uplims=True, fmt="vk", ms=8, ecolor='gray', label='upperlimit', capsize=8, elinewidth=4, capthick=1.5)

ax.set_ylim(min(flux)-0.5, 3*max(flux))
ax.set_yscale("log")
ax.set_xscale("log")


ax.set_ylabel(r'$\rm S_{\nu}$ [mJy]', fontsize=16)
ax.set_xlabel(r'$\rm \lambda [\mu$m]', fontsize=16)
ax.set_title("SED ", x=0.5, y=1.125)

led = plt.legend(loc='best', fontsize=15, numpoints=1,
                 fancybox=True, borderpad=0.5,
                 handlelength=2, labelspacing=0.1)
led.get_frame().set_alpha(0)
led.get_frame().set_edgecolor('white')
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

plt.show(block=False)   # so it doesn't hold off the next command
response = raw_input('Save fig?: (y/n)')

if response.lower() in ['y', 'yes']:
    savefigure(f, plotName)
else:
    print('...Exiting...')
