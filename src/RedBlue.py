"""
Author: Daisy Leung

Last edited: November 5 2015

Purpose:
- Plot Red wing and Blue wing different color overlay on HST

History:

"""

from astropy import log
log.setLevel('ERROR')
import glob
import matplotlib.pyplot as plt
import matplotlib as mpl
from APLpySetup import *

path = '../HST/'
Plotpath = '../Figures/'

label = dict.fromkeys(['F555W'])
for k in label.iterkeys():
    files = glob.glob(path + '*' + k + '*.fits')
    label[k] = files
# print label

label_pdbi = dict.fromkeys(['red', 'blue'])
path_pdbi = '../PdBI/data/5Nov15/'

for k in label_pdbi.iterkeys():
    file_pdBI = glob.glob(path_pdbi + '*' + k + '.fits')
    label_pdbi[k] = file_pdBI
# print label_pdbi

figC = plt.figure(1, figsize=(7, 7))
########################################
# user define area
########################################
ra_center = 172.96418        # 172.96426
dec_center = -12.5330955    # -12.533213
sizep = 0.00158066

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

########################################
# intialize base figure
########################################

fig1 = aplpy.FITSFigure(label['F555W'][0],
                        figure=figC, subplot=[x0, row_a, full_width, dy])
fig1.show_grayscale(stretch='log', vmin=-0.01869, vmax=151, vmid=-0.025)
fig1.set_theme('pulication')
########################################
# Contours
########################################
fig1.show_contour(label_pdbi['blue'][0], colors="lime", levels=sigma_contour_array(
    sigma_blue), linewidths=2)  # , layer='fg')

fig1.show_contour(label_pdbi['red'][0], colors="blue",
                  levels=sigma_contour_array(sigma_red), linewidths=2)

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

put_label(fig1, 0.20, 0.95, 'HST F555W', 'titleBand')
put_label(fig1, 0.20, 0.85, 'CO (2-1)', 'titleBand')
########################################
# Colorbar
########################################
# axisflin = fig.add_axes([0.92,0.19,0.02,0.68])
# normflin = mpl.colors.Normalize(vmin=min_line, vmax=max_line)
# cbflin = mpl.colorbar.ColorbarBase(axisflin, cmap=mpl.cm.jet, norm=normflin, orientation='vertical')
# cbflin.set_label('mJy')
# fig.canvas.draw()
# fig_line.canvas.draw()
# plt.show()

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
        #        os.system('rm -rf ' + C[:-1] + '.png' + ' ' + C[:-1] + '.eps')
        #        os.system('rm -rf ' + D[:-1] + '.png' + ' ' + D[:-1] + '.eps')
        #        os.system('rm -rf ' + CD[:-1] + '.png' + ' ' + CD[:-1] + '.eps')
        #        figC.savefig(Plotpath + C[:-1] + '.eps', dpi=600)
        figC.savefig(Plotpath + 'F555W_pdbi.eps', dpi=600)
#         figRed.savefig(Plotpath + 'F555W_pdbiRed.eps', dpi=600)
    else:
        #        figC.canvas.draw()
        plt.show()
