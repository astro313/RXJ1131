"""

Plot mom0 of all channels overlay on HST (single panel)


Author: Daisy Leung


Last edited: Dec 31 2015


History:
--------
31 Dec 2015: use the linear shifted HST F555W image

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

label = dict.fromkeys(['F555W_drzcroplinearShift'])
for k in label.iterkeys():
    files = glob.glob(path + '*' + k + '*.fits')
    label[k] = files
# print label

label_pdbi = dict.fromkeys(['mom0'])
path_pdbi = '../PdBI/data/14Oct15/'

for k in label_pdbi.iterkeys():
    file_pdBI = glob.glob(path_pdbi + '*' + k + '.fits')
    label_pdbi[k] = file_pdBI

figC = plt.figure(1, figsize=(7, 7))
figC.clf()
########################################
# user define area
########################################
ra_center = 172.96418
dec_center = -12.5330955
sizep = 0.0025

ra_cross, dec_cross = ra_center, dec_center
row_a = 0.10
width = 0.35
full_width = 0.8
x_gap = 0.05
x0 = 0.15
dy = 0.90
sigma = 0.43    # 0.32e-3 Jy/B

line_min = 0.035668
line_max = 7.56445

inverted_HSTleft = True
tickcolor = 'k' if inverted_HSTleft else 'white'
########################################
# intialize base figure
########################################
fig1 = aplpy.FITSFigure(label['F555W_drzcroplinearShift'][0],
                        figure=figC, subplot=[x0, row_a, full_width, dy])
fig1.show_grayscale(stretch='log', vmin=-0.01869, vmax=151, vmid=-0.025)
if inverted_HSTleft: fig1.set_theme('publication')
########################################
# Contours
########################################
fig1.show_contour(label_pdbi['mom0'][0], colors=tickcolor, levels=sigma_contour_array(sigma), linewidths=2)  # , layer='fg')

########################################
# beam
########################################
fig1.show_beam(major=4.44/3600, minor=1.95/3600, angle=13., edgecolor='grey', facecolor='orange', linestyle='solid', linewidth=3, frame=True, alpha=0.65)


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
markers_cross(fig1, ra_cross, dec_cross, layer='marker_set_1', ec=tickcolor)

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
    if saveFig == True:
        outname = 'F555WCO21_mom0_single.eps'
        if inverted_HSTleft: outname = outname.replace('.eps', '.invertedgray.eps')
        figC.savefig(Plotpath + outname, dpi=600)
    else:
        #        figC.canvas.draw()
        plt.show()
