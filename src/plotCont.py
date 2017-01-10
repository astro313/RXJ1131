"""
Author: Daisy Leung

Last edited: 9 Jan 2016

Purpose:
- Overlay 2mm Cont contour on HST

History:
Jan 9 2016:
    - zoom in
    - use np.sort() w/ contour level, to satisfy mpl version > 1.5.3
Nov 21 2016:
    - add inverted HST
    - no marker, zoom in
July 20 2016:
    - change zoom to match VLA image (HSTVLA.py)
May 27 2016:
    - update aspect ratio
May 17 2016:
    - add beam
Jan 01 2016:
    - update marker
Dec 31 2015:
    - use linear shifted F555W image
Dec 30 2015:
    - make pretty layout
Sept 19 2015:
    - created script

"""

from astropy import log
log.setLevel('ERROR')
import glob
import matplotlib.pyplot as plt
import matplotlib as mpl
import numpy as np
from APLpySetup import *

path = '../HST/'
Plotpath = '../Figures/'

mpl.rcdefaults()
font = {'family': 'Arial Narrow',
        'weight': 'normal',
        'size': 11.5}
xtick = {'major.size': 10,
         'minor.size': 6}
ytick = {'major.size': 10,
         'minor.size': 6}
mpl.rc('font', **font)
mpl.rcParams.update({'font.size': 11.5})
mpl.rc('xtick', **xtick)

label = dict.fromkeys(['F555W_drzcroplinearShift'])
for k in label.iterkeys():
    files = glob.glob(path + '*' + k + '*.fits')
    label[k] = files
# print label

label_pdbi = dict.fromkeys(['cont'])
path_pdbi = '../PdBI/data/10Nov15/'

for k in label_pdbi.iterkeys():
    file_pdBI = glob.glob(path_pdbi + '*' + k + '.fits')
    label_pdbi[k] = file_pdBI
# print label_pdbi

figC = plt.figure(1, figsize=(12, 8))
# figRed = plt.figure(2, figsize=(12, 7))
########################################
# user define area
########################################
ra_center = 172.96434563888
dec_center = -12.5328629
sizep = 0.0015 # 0.002058066

ra_cross, dec_cross = ra_center, dec_center
row_a = 0.1
width = 0.35
full_width = 0.65
x_gap = 0.05
x0 = 0.2
dy = 0.9
sigma = 0.085e-3     # theoretical: 0.0816E-03; off-region in map: 0.0888e-3

inverted_HSTleft = True
tickcolor = 'k' if inverted_HSTleft else 'white'
########################################
# intialize base figure
########################################
fig1 = aplpy.FITSFigure(label['F555W_drzcroplinearShift'][0],
                        figure=figC, subplot=[x0, row_a, full_width, dy])
fig1.show_grayscale(stretch='log', vmin=-0.00256, vmax=1.922, vmid=-0.02)
if inverted_HSTleft: fig1.set_theme('publication')

# fig1 = aplpy.FITSFigure(label['F555W_drzcroplinearShift'][0], \
#         figure=figC, subplot=[x0,row_a,full_width,dy])
# fig1.show_grayscale(stretch='log', vmin=-0.01869, vmax=125, vmid=-0.025)
#figO = aplpy.FITSFigure(label['F555W_drzcroplinearShift'][0], \
#        figure=figC, subplot=[x0+width+2*x_gap, row_a, width, dy])
#figO.show_grayscale(stretch='log', vmin=-0.01869, vmax=151, vmid=-0.025)


# figHST = aplpy.FITSFigure(label['F555W_drzcroplinearShift'][0], \
#         figure=figRed, subplot=[x0,row_a,width,dy])
# figHST.show_grayscale(stretch='log', vmin=-0.01869, vmax=151, vmid=-0.025)
# figred = aplpy.FITSFigure(label['F555W_drzcroplinearShift'][0], \
#         figure=figRed, subplot=[x0+width+2*x_gap, row_a, width, dy])
# figred.show_grayscale(stretch='log', vmin=-0.01869, vmax=151, vmid=-0.025)

########################################
# Contours
########################################
fig1.show_contour(label_pdbi['cont'][0], colors="black", levels=np.sort(sigma_contour_array(sigma)), linewidths=1.62)#, layer='fg')

# figred.show_contour(label_pdbi['red'][0], colors="lime", levels=sigma_contour_array(sigma_red), linewidths=2)

########################################
# beam
########################################
fig1.show_beam(major=4.44/3600, minor=1.95/3600, angle=13., edgecolor='gray', facecolor='gray', linestyle='solid', linewidth=3, frame=False, alpha=0.8)

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
#standard_plot_setup(figO, ra_center, dec_center, sizep, tickc='white')
# standard_plot_setup(figHST, ra_center, dec_center, sizep, tickc='white')
# standard_plot_setup(figred, ra_center, dec_center, sizep, tickc='white')
#figO.tick_labels.hide()
# fig0.axis_labels.hide()
# figred.tick_labels.hide()
# figred.axis_labels.hide()

########################################
# markers
########################################
# markers_cross(fig1, ra_cross, dec_cross, layer='marker_set_1')
#markers_cross(figO, ra_cross, dec_cross, layer='marker_set_1')
# markers_cross(figHST, ra_cross, dec_cross, layer='marker_set_1')
# markers_cross(figred, ra_cross, dec_cross, layer='marker_set_1')

########################################
# Labels
########################################
# if '_' in sym[:-1]: symf = sym.replace('_', ' ')

put_label(fig1, 0.15725, 0.95, 'RXJ1131', 'titleObj', c=tickcolor)
put_label(fig1, 0.38, 0.9, 'HST F555W, 2 mm Continuum', 'titleBand', c=tickcolor)
# put_label(figHST, 0.20, 0.95, 'HST F555W', 'titleBand')
# put_label(figHST, 0.1825, 0.9, 'RXJ1131', 'titleObj')
# put_label(figred, 0.40, 0.95, 'HST F555W, CO Red wing', 'titleBand')
# put_label(figred, 0.4025, 0.9, 'RXJ1131', 'titleObj')

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
    if saveFig:
        figC.savefig(Plotpath + 'F555W_ContPdBI.eps', dpi=600, bbox_inches='tight')
#         figRed.savefig(Plotpath + 'F555W_pdbiRed.eps', dpi=600)
    else:
#        figC.canvas.draw()
        plt.show()
