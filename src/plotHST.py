"""
Author: Daisy Leung

Last edited: 01 Jan 2016

Purpose:
- Plot HST image F555W (left), plot mom0 of all channels overlay on HST (right)

History:
Jan 01 2016:
    update marker
Dec 31 2015:
    use the linear shifted F555W image
    added option to plot left panel HST image in inverted grayscale

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
print label

label_pdbi = dict.fromkeys(['mom0'])
path_pdbi = '../PdBI/data/14Oct15/'

for k in label_pdbi.iterkeys():
    file_pdBI = glob.glob(path_pdbi + k + '.fits')
    label_pdbi[k] = file_pdBI
# print label_pdbi

figC = plt.figure(1, figsize=(12, 7))
# figRed = plt.figure(2, figsize=(12, 7))
########################################
# user define area
########################################
ra_center = 172.96434563888
dec_center = -12.5328629
sizep = 0.00178066

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
########################################
# intialize base figure
########################################
# fig1 = aplpy.FITSFigure(label['F555W_drzcroplinearShift'][0],
#                         figure=figC, subplot=[x0, row_a, full_width, dy])
# fig1.show_grayscale(stretch='log', vmin=-0.01869, vmax=151, vmid=-0.025)
# #fig1.set_theme('publication')   # inverted grayscale

fig1 = aplpy.FITSFigure(label['F555W_drzcroplinearShift'][0], \
        figure=figC, subplot=[x0,row_a,width,dy])
fig1.show_grayscale(stretch='log', vmin=-0.01869, vmax=151, vmid=-0.025)
if inverted_HSTleft: fig1.set_theme('publication')


figO = aplpy.FITSFigure(label['F555W_drzcroplinearShift'][0], \
        figure=figC, subplot=[x0+width, row_a, width, dy])
figO.show_grayscale(stretch='log', vmin=-0.01869, vmax=151, vmid=-0.025)


########################################
# Contours
########################################
figO.show_contour(label_pdbi['mom0'][0], colors="lime", levels=sigma_contour_array(sigma), linewidths=2)#, layer='fg')


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
tickcolor = 'k' if inverted_HSTleft else 'white'
standard_plot_setup(fig1, ra_center, dec_center, sizep, tickc=tickcolor)
standard_plot_setup(figO, ra_center, dec_center, sizep, tickc='white')

figO.tick_labels.hide()
figO.axis_labels.hide()
# figred.tick_labels.hide()
# figred.axis_labels.hide()

########################################
# markers
########################################
markers_cross(fig1, ra_cross, dec_cross, layer='marker_set_1')
markers_cross(figO, ra_cross, dec_cross, layer='marker_set_1')
# markers_cross(figHST, ra_cross, dec_cross, layer='marker_set_1')
# markers_cross(figred, ra_cross, dec_cross, layer='marker_set_1')

########################################
# Labels
########################################
put_label(fig1, 0.20, 0.95, 'HST F555W', 'titleBand', c=tickcolor)
put_label(figO, 0.40, 0.95, 'HST F555W, CO (2-1)', 'titleBand')
put_label(figO, 0.4025, 0.9, 'RXJ1131', 'titleObj')
# put_label(figHST, 0.20, 0.95, 'HST F555W', 'titleBand')
# put_label(figHST, 0.1825, 0.9, 'RXJ1131', 'titleObj')

labsize = 'xx-large'
labc = 'white'
#put_label(fvla, 0.80, 0.925, '(a)', 'ref', c=labc, s=labsize)
#put_label(fcont, 0.80, 0.925, '(b)', 'ref', c=labc, s=labsize)


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
        outname = 'F555WCO21_mom0.eps'
        if inverted_HSTleft: outname = outname.replace('.eps', '.invertedgray.eps')
        figC.savefig(Plotpath + outname, dpi=600)
    else:
#        figC.canvas.draw()
        plt.show()
