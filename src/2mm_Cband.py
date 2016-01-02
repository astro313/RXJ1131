"""
Author: Daisy Leung

Last edited: 01 Jan 2016

Purpose:
- Plot VLA C band continuum as grayscale, overlay 2 mm PdBI continuum contours

History:
01 Jan 2016: update marker
30 Dec 2015: initial
"""

from astropy import log
log.setLevel('ERROR')
import glob
import matplotlib.pyplot as plt
import matplotlib as mpl
from APLpySetup import *

Plotpath = '../Figures/'

path = '../Radio/Cband/'
label = dict.fromkeys(['C_R0'])
for k in label.iterkeys():
    files = glob.glob(path + '*' + k + '*')
    label[k] = files
# print label

path_pdbi = '../PdBI/data/10Nov15/'
label_pdbi = dict.fromkeys(['cont'])

for k in label_pdbi.iterkeys():
    file_pdBI = glob.glob(path_pdbi + '*' + k + '.fits')
    label_pdbi[k] = file_pdBI
# print label_pdbi

figC = plt.figure(1, figsize=(7, 7))
fig2 = plt.figure(2, figsize=(12, 5))
########################################
# user define area
########################################
ra_center = 172.96434563888
dec_center = -12.5328629
sizeVLA = 0.00135
sizep = 0.00235066

ra_cross, dec_cross = ra_center, dec_center
row_a = 0.10
width = 0.35
full_width = 0.8
x_gap = 0.05
x0 = 0.15
dy = 0.875
sigma_Cband = 1.3e-5
sigma_2mm = 0.085e-3     # theoretical: 0.0816E-03; off-region in map: 0.0888e-3

########################################
# intialize base figure
########################################
# Single panel: 2mm cont contour on 5 GHz Cont
fig1 = aplpy.FITSFigure(label['C_R0'][0], \
        figure=figC, subplot=[x0,row_a,full_width,dy])
fig1.show_grayscale(stretch='arcsinh', vmin=-4e-5, vmax=7e-4, vmid=-1e-4)
# [-6.10554e-05, 0.000739422]

# Double panel: left: VLA only ; right: same as above
figVLA = aplpy.FITSFigure(label['C_R0'][0], \
        figure=fig2, subplot=[x0,row_a,width,dy])
figVLA.show_grayscale(stretch='arcsinh', vmin=-4e-5, vmax=7e-4, vmid=-1e-4)
figO = aplpy.FITSFigure(label['C_R0'][0], \
        figure=fig2, subplot=[x0+width+2*x_gap, row_a, width, dy])
figO.show_grayscale(stretch='arcsinh', vmin=-4e-5, vmax=7e-4, vmid=-1e-4)

########################################
# Contours
########################################

fig1.show_contour(label_pdbi['cont'][0], colors="lime", levels=sigma_contour_array(sigma_2mm), linewidths=2)

figVLA.show_contour(label['C_R0'][0], colors="blue", levels=sigma_contour_VLA(sigma_Cband), linewidths=2)
figO.show_contour(label_pdbi['cont'][0], colors="lime", levels=sigma_contour_array(sigma_2mm), linewidths=2)

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
standard_plot_setup(figVLA, ra_center, dec_center, sizeVLA, tickc='w')
standard_plot_setup(figO, ra_center, dec_center, sizep, tickc='white')
# figO.tick_labels.hide()
figO.axis_labels.hide()

########################################
# markers
########################################
markers_cross(fig1, ra_cross, dec_cross, layer='marker_set_1', ec='red')
markers_cross(figVLA, ra_cross, dec_cross, layer='marker_set_1', ec='red')
markers_cross(figO, ra_cross, dec_cross, layer='marker_set_1', ec='red')

########################################
# Labels
########################################
# if '_' in sym[:-1]: symf = sym.replace('_', ' ')

put_label(fig1, 0.250, 0.95, 'VLA 5 GHz, PdBI 2 mm', 'titleBand',c='y')
put_label(fig1, 0.15725, 0.9, 'RXJ1131', 'titleObj',c='y')
put_label(figVLA, 0.170, 0.95, 'VLA 5 GHz', 'titleBand')
put_label(figO, 0.4025, 0.95, 'VLA 5 GHz, PdBI 2 mm', 'titleObj')

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
        figC.savefig(Plotpath + 'Cont2mm_5GHz_single.eps', dpi=600)
        fig2.savefig(Plotpath + 'Cont2mm_5GHz_double.eps', dpi=600)
    else:
#        figC.canvas.draw()
        plt.show()
