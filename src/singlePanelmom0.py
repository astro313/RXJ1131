"""

Plot mom0 of all channels overlay on HST (single panel)


Author: Daisy Leung


Last edited: 29 Aug 2016


History:
--------
29 Aug 2016
    - changing save filename to exclude . between name and ext, else would mess up the prepare to submit codes
20 July 2016
    - change zoom & HST contrast to match RedBlue.py
28 May 2016:
    - increase font size
    - longer ticks
    - zoom in a little to avoid tick label overlaps
17 May 2016:
    - pythonic
16 May 2016:
    - update font
    - replace mom0 wiith centralizedCube4GILDAS-CASA_ch126-160_mom0.fits from 15May16/
    - update sigma
01 Jan 2016:
    - update marker
31 Dec 2015:
    - use the linear shifted HST F555W image

Note:
-----
- option to plot inverted HST


"""

from astropy import log
log.setLevel('ERROR')
import glob
import matplotlib.pyplot as plt
import matplotlib as mpl
from APLpySetup import *

path = '../HST/'
Plotpath = '../Figures/'

mpl.rcdefaults()
font = {'family': 'Arial Narrow',
        'weight': 'normal',
        'size': 13.5}
xtick = {'major.size': 10,
         'minor.size': 6.5}
ytick = {'major.size': 10,
         'minor.size': 6.5}
mpl.rc('font', **font)
mpl.rc('xtick', **xtick)

label = dict.fromkeys(['F555W_drzcroplinearShift'])
for k in label.iterkeys():
    files = glob.glob(path + '*' + k + '*.fits')
    label[k] = files
# print label

label_pdbi = dict.fromkeys(['CASA_ch126-160'])
#
path_pdbi = '../PdBI/data/15May16/'

for k in label_pdbi.iterkeys():
    file_pdBI = glob.glob(path_pdbi + '*' + k + '*.fits')
    label_pdbi[k] = file_pdBI

figC = plt.figure(1)
figC.clf()
########################################
# user define area
########################################
ra_center = 172.96434563888
dec_center = -12.5328629
sizep = 0.002

ra_cross, dec_cross = ra_center, dec_center
row_a = 0.1
width = 0.35
full_width = 0.65
x_gap = 0.05
x0 = 0.2
dy = 0.9
sigma = 0.305       # 0.43    # 0.32e-3 Jy/B

line_min = 0.035668
line_max = 7.56445

inverted_HSTleft = True
tickcolor = 'k' if inverted_HSTleft else 'white'
########################################
# intialize base figure
########################################
fig1 = aplpy.FITSFigure(label['F555W_drzcroplinearShift'][0],
                        figure=figC, subplot=[x0, row_a, full_width, dy])
fig1.show_grayscale(stretch='log', vmin=-0.00256, vmax=1.922, vmid=-0.02)
if inverted_HSTleft: fig1.set_theme('publication')
########################################
# Contours
########################################
fig1.show_contour(label_pdbi['CASA_ch126-160'][0], colors=tickcolor, levels=sigma_contour_array(sigma), linewidths=1.42)  # , layer='fg')

########################################
# beam
########################################
fig1.show_beam(major=4.44/3600, minor=1.95/3600, angle=13., edgecolor='grey', facecolor='grey', linestyle='solid', linewidth=3, frame=False, alpha=0.65)

########################################
# scale bar
########################################
lg_1arcsec = 1. / 3600
# lg_20kpc_fg = lg_1arcsec * 20./scale_radio
# lg_20kpc_bg = lg_1arcsec * 20./scale_SMG
# setup_scalebar(fvla, lg_20kpc_fg, str('20kpc'))

########################################
# axes
########################################
standard_plot_setup(fig1, ra_center, dec_center, sizep, tickc=tickcolor)

########################################
# markers
########################################
markers_cross(fig1, ra_cross, dec_cross, layer='marker_set_1', ec='red', lw=1.2)

########################################
# Labels
########################################
# if '_' in sym[:-1]: symf = sym.replace('_', ' ')
put_label(fig1, 0.82, 0.90, 'HST F555W', 'titleBand1', c=tickcolor)
put_label(fig1, 0.82, 0.95, 'CO (2-1)', 'titleBand2', c=tickcolor)

########################################
########################################
if __name__ == '__main__':
    """
    run script.py True
    sys.argv[1] determines whether to save all the plots or not
    """
    import sys
    import os
    if len(sys.argv) < 2:
        errmsg = "Invalid number of arguments: {0:d}\n  run script.py Save_True"
        raise IndexError(errmsg.format(len(sys.argv)))
    saveFig = True if sys.argv[1].lower() == 'true' else False
    if saveFig:
        outname = 'F555WCO21_mom0_single.eps'
        if inverted_HSTleft: outname = outname.replace('.eps', '_invertedgray.eps')
        figC.savefig(Plotpath + outname, dpi=600, bbox_inches='tight')
    else:
        #        figC.canvas.draw()
        plt.show()
