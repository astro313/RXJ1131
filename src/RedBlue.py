"""
Author: Daisy Leung

Last edited: 03 May 2016

Purpose:
    Plot Red wing and Blue wing different color overlay on HST
    channels: blue: [150, 142]; velocity [-423.79, -251.39]
              red: [142, 134]; velocity [-251.39, -90.22]


History:
03 May 2016: change label color from yellow to black, since we inverted HST image, else can't read them
30 April 2016: Change contour colors, inverted color for HST, change HST contrast based on ds9
01 Jan 2016: update marker
31 Dec 2015: Use linear shifted F555W image

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

label_pdbi = dict.fromkeys(['red', 'blue'])
path_pdbi = '../PdBI/data/5Nov15/'    # copied from 14Oct15

for k in label_pdbi.iterkeys():
    file_pdBI = glob.glob(path_pdbi + '*' + k + '.fits')
    label_pdbi[k] = file_pdBI
# print label_pdbi

figC = plt.figure(1, figsize=(7, 7))
########################################
# user define area
########################################
ra_center = 172.96434563888
dec_center = -12.5328629
sizep = 0.002                # 0.00158066

ra_cross, dec_cross = ra_center, dec_center
row_a = 0.10
width = 0.35
full_width = 0.8
x_gap = 0.05
x0 = 0.15
dy = 0.90
sigma_blue = 0.1
sigma_red = 0.15

red_line_min = -0.488125
red_line_max = 2.19431
blue_line_min = -0.343676
blue_line_max = 1.72505

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
fig1.show_contour(label_pdbi['blue'][0], colors="blue", levels=sigma_contour_array(
    sigma_blue), linewidths=2)  # , layer='fg')

fig1.show_contour(label_pdbi['red'][0], colors="red",
                  levels=sigma_contour_array(sigma_red), linewidths=2)

########################################
# beam
########################################
fig1.show_beam(major=4.44/3600, minor=1.95/3600, angle=13., edgecolor='grey', facecolor='orange', linestyle='solid', linewidth=3, frame=False, alpha=0.8)


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
standard_plot_setup(fig1, ra_center, dec_center, sizep, tickc='white')

########################################
# markers
########################################
markers_cross(fig1, ra_cross, dec_cross, layer='marker_set_1')

########################################
# Labels
########################################
# if '_' in sym[:-1]: symf = sym.replace('_', ' ')
put_label(fig1, 0.82, 0.90, 'HST F555W', 'titleBand1', c='black')
put_label(fig1, 0.82, 0.95, 'CO (2-1)', 'titleBand2', c='black')

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
        figC.savefig(Plotpath + 'F555W_REDBLUE.png', dpi=600)
    else:
        #        figC.canvas.draw()
        plt.show()
