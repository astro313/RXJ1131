"""
Author: Daisy Leung

Last edited: Sept 19 2015

Purpose:
- Plot HST image F555W

History:

"""

from astropy import log
log.setLevel('ERROR')
import glob
import matplotlib.pyplot as plt
from APLpySetup import *

path = '../HST/'
Plotpath = '../Figures/'

label = dict.fromkeys(['F555W'])
for k in label.iterkeys():
    files = glob.glob(path + '*' + k + '*.fits')
    label[k] = files
# print label

label_pdbi = dict.fromkeys(['red', 'blue'])
path_pdbi = '../PdBI/data/10Sep15/'   ## << This flux is wrong

for k in label_pdbi.iterkeys():
    file_pdBI = glob.glob(path_pdbi + k + '_average.fits')
    label_pdbi[k] = file_pdBI
# print label_pdbi

figC = plt.figure(1, figsize=(12, 7))
figRed = plt.figure(2, figsize=(12, 7))
########################################
# user define area
########################################
ra_center = 172.96418 #172.96426
dec_center = -12.5330955 #-12.533213
sizep = 0.00258066

ra_cross, dec_cross = ra_center, dec_center
row_a = 0.10
width = 0.35
full_width = 0.8
x_gap = 0.05
x0 = 0.15
dy = 0.90
sigma_blue = 0.328e-3    # 0.32e-3 Jy/B
sigma_red = 0.31205e-3


########################################
# intialize base figure
########################################
# fig1 = aplpy.FITSFigure(label['F555W'][0],
#                         figure=figC, subplot=[x0, row_a, full_width, dy])#, north=True)
# # fig1.show_grayscale(stretch='log', vmin=-0.45869, vmax=250, vmid=-0.85)
# fig1.show_grayscale(stretch='log', vmin=-0.01869, vmax=151, vmid=-0.025)
# #fig1.set_theme('publication')   # inverted grayscale

fig1 = aplpy.FITSFigure(label['F555W'][0], \
        figure=figC, subplot=[x0,row_a,width,dy])
fig1.show_grayscale(stretch='log', vmin=-0.01869, vmax=151, vmid=-0.025)
figO = aplpy.FITSFigure(label['F555W'][0], \
        figure=figC, subplot=[x0+width+2*x_gap, row_a, width, dy])
figO.show_grayscale(stretch='log', vmin=-0.01869, vmax=151, vmid=-0.025)


figHST = aplpy.FITSFigure(label['F555W'][0], \
        figure=figRed, subplot=[x0,row_a,width,dy])
figHST.show_grayscale(stretch='log', vmin=-0.01869, vmax=151, vmid=-0.025)
figred = aplpy.FITSFigure(label['F555W'][0], \
        figure=figRed, subplot=[x0+width+2*x_gap, row_a, width, dy])
figred.show_grayscale(stretch='log', vmin=-0.01869, vmax=151, vmid=-0.025)

########################################
# Contours
########################################
figO.show_contour(label_pdbi['blue'][0], colors="lime", levels=sigma_contour_array(sigma_blue), linewidths=2)#, layer='fg')

figred.show_contour(label_pdbi['red'][0], colors="lime", levels=sigma_contour_array(sigma_red), linewidths=2)

########################################
# Compass, bug in APLpy
########################################
# fig1_com = aplpy.overlays.Compass()
# # import pdb;pdb.set_trace()
# fig1_com.show_compass(color='white', corner=1, length=0.05) # corner=1 Top-right
# fig1_com.compass._compass[0].set_arrowstyle('-')       # Remove head from East arrow

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
standard_plot_setup(figO, ra_center, dec_center, sizep, tickc='white')
standard_plot_setup(figHST, ra_center, dec_center, sizep, tickc='white')
standard_plot_setup(figred, ra_center, dec_center, sizep, tickc='white')
figO.tick_labels.hide()
figO.axis_labels.hide()
figred.tick_labels.hide()
figred.axis_labels.hide()

########################################
# markers
########################################
markers_cross(fig1, ra_cross, dec_cross, layer='marker_set_1')
markers_cross(figO, ra_cross, dec_cross, layer='marker_set_1')
markers_cross(figHST, ra_cross, dec_cross, layer='marker_set_1')
markers_cross(figred, ra_cross, dec_cross, layer='marker_set_1')

########################################
# Labels
########################################
# if '_' in sym[:-1]: symf = sym.replace('_', ' ')

put_label(fig1, 0.20, 0.95, 'HST F555W', 'titleBand')
put_label(fig1, 0.1825, 0.9, 'RXJ1131', 'titleObj')
put_label(figO, 0.40, 0.95, 'HST F555W, CO Blue wing', 'titleBand')
put_label(figO, 0.4025, 0.9, 'RXJ1131', 'titleObj')
put_label(figHST, 0.20, 0.95, 'HST F555W', 'titleBand')
put_label(figHST, 0.1825, 0.9, 'RXJ1131', 'titleObj')
put_label(figred, 0.40, 0.95, 'HST F555W, CO Red wing', 'titleBand')
put_label(figred, 0.4025, 0.9, 'RXJ1131', 'titleObj')

labsize = 'xx-large'
labc = 'white'
#put_label(fvla, 0.80, 0.925, '(a)', 'ref', c=labc, s=labsize)
#put_label(fcont, 0.80, 0.925, '(b)', 'ref', c=labc, s=labsize)
#put_label(flin, 0.80, 0.925, '(c)', 'ref', c=labc, s=labsize)
#put_label(fSMA, 0.80, 0.925, '(d)', 'ref', c=labc, s=labsize)

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
        figC.savefig(Plotpath + 'F555W_pdbiBlue.eps', dpi=600)
        figRed.savefig(Plotpath + 'F555W_pdbiRed.eps', dpi=600)
    else:
#        figC.canvas.draw()
        plt.show()
