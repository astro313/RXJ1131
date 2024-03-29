"""
Author: Daisy Leung

Last edited: 01 Jan 2016

Purpose:
- Plot continuum as grayscale overlay CO(2-1) mom0

History:
01 Jan 2016: update marker postion
30 Dec 2015: initial file

"""

from astropy import log
log.setLevel('ERROR')
import glob
import matplotlib.pyplot as plt
import matplotlib as mpl
from APLpySetup import *

Plotpath = '../Figures/'

path = '../PdBI/data/10Nov15/'
label = dict.fromkeys(['cont'])
for k in label.iterkeys():
    files = glob.glob(path + '*' + k + '*.fits')
    label[k] = files
# print label

label_pdbi = dict.fromkeys(['mom0'])
path_pdbi = '../PdBI/data/14Oct15/'

for k in label_pdbi.iterkeys():
    file_pdBI = glob.glob(path_pdbi + '*' + k + '.fits')
    label_pdbi[k] = file_pdBI
# print label_pdbi

figC = plt.figure(1, figsize=(7, 7))
# figRed = plt.figure(2, figsize=(12, 7))
########################################
# user define area
########################################
ra_center = 172.96434563888
dec_center = -12.5328629
sizep = 0.00408066

ra_cross, dec_cross = ra_center, dec_center
row_a = 0.10
width = 0.35
full_width = 0.8
x_gap = 0.05
x0 = 0.15
dy = 0.90
sigma = 0.085e-3     # theoretical: 0.0816E-03; off-region in map: 0.0888e-3
sigma_mom0 = 0.43    # 0.32e-3 Jy/B

line_min = 0.035668
line_max = 7.56445

########################################
# intialize base figure
########################################
fig1 = aplpy.FITSFigure(label['cont'][0], \
        figure=figC, subplot=[x0,row_a,full_width,dy])
# fig1.show_grayscale(stretch='arcsinh', vmin=-0.0005, vmax=0.001, vmid=0)
#figO = aplpy.FITSFigure(label['F555W'][0], \
#        figure=figC, subplot=[x0+width+2*x_gap, row_a, width, dy])
#figO.show_grayscale(stretch='log', vmin=-0.01869, vmax=151, vmid=-0.025)


########################################
# Contours
########################################
fig1.show_contour(label['cont'][0], colors="lime", levels=sigma_contour_array(sigma), linewidths=2)
fig1.show_contour(label_pdbi['mom0'][0], colors="orange", levels=sigma_contour_array(sigma_mom0), linewidths=2)

# figred.show_contour(label_pdbi['red'][0], colors="lime", levels=sigma_contour_array(sigma_red), linewidths=2)

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
standard_plot_setup(fig1, ra_center, dec_center, sizep, tickc='black')
#standard_plot_setup(figO, ra_center, dec_center, sizep, tickc='white')
# standard_plot_setup(figHST, ra_center, dec_center, sizep, tickc='white')
# standard_plot_setup(figred, ra_center, dec_center, sizep, tickc='white')
#figO.tick_labels.hide()
#figO.axis_labels.hide()
# figred.tick_labels.hide()
# figred.axis_labels.hide()

########################################
# markers
########################################
markers_cross(fig1, ra_cross, dec_cross, layer='marker_set_1', ec='k')
#markers_cross(figO, ra_cross, dec_cross, layer='marker_set_1')
# markers_cross(figHST, ra_cross, dec_cross, layer='marker_set_1')
# markers_cross(figred, ra_cross, dec_cross, layer='marker_set_1')

########################################
# Labels
########################################
# if '_' in sym[:-1]: symf = sym.replace('_', ' ')

# put_label(fig1, 0.20, 0.95, 'HST F555W', 'titleBand')
put_label(fig1, 0.2750, 0.95, '2 mm Continuum, CO(2-1)', 'titleBand',c='k')
put_label(fig1, 0.15725, 0.9, 'RXJ1131', 'titleObj',c='k')
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
    if saveFig == True:
        figC.savefig(Plotpath + 'ContCO21.eps', dpi=600)
#         figRed.savefig(Plotpath + 'F555W_pdbiRed.eps', dpi=600)
    else:
#        figC.canvas.draw()
        plt.show()
