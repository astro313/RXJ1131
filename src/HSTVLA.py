"""
Author: Daisy Leung

Last edited: 11 Nov 2016

Purpose:
- Plot single panel 5GHz VLA continuum on HST F555W, investigate continuum offset

History:
Nov 11 2016:
    - update contour intervals to +/-10%, 20%,... peak
Nov 2 2016:
    - change lw from 2 to 1.2
May 27 2016:
    - aspect ratio
May 17 2016:
    - update cosmetics
    - add beam
Jan 01 2016:
    - update marker
Dec 31 2015:
    - don't plot HST comparison, only use the astrometric correct F555W image
    - plot as well the astrometric corrected image - F555W_ContVLA_astrometry.eps
Dec 30 2015:
    create script

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
        'size': 11.5}
xtick = {'major.size': 10,
         'minor.size': 6}
ytick = {'major.size': 10,
         'minor.size': 6}
mpl.rc('font', **font)
mpl.rc('xtick', **xtick)

label = dict.fromkeys(['F555W_drzcroplinearShift'])
for k in label.iterkeys():
    files = glob.glob(path + '*' + k + '*fits')
    label[k] = files

path_vla = '../Radio/Cband/'
label_vla  = dict.fromkeys(['C_R0'])

for k in label_vla.iterkeys():
    file_vla = glob.glob(path_vla + '*' + k + '*')
    label_vla[k] = file_vla
# print label_vla

figC = plt.figure(1)
# figRed = plt.figure(2, figsize=(12, 7))
########################################
# user define area
########################################
ra_center = 172.96434563888
dec_center = -12.5328629
sizep = 0.00215

ra_cross, dec_cross = ra_center, dec_center
row_a = 0.1
width = 0.35
full_width = 0.65
x_gap = 0.05
x0 = 0.2
dy = 0.9
sigma = 1.3e-5
Smax = 0.74e-3 # VLA

########################################
# intialize base figure
########################################
# fig1 = aplpy.FITSFigure(label['F555W_drzcroplinearShift'][0],
#                         figure=figC, subplot=[x0, row_a, full_width, dy])#, north=True)
# # fig1.show_grayscale(stretch='log', vmin=-0.45869, vmax=250, vmid=-0.85)
# fig1.show_grayscale(stretch='log', vmin=-0.01869, vmax=151, vmid=-0.025)
# #fig1.set_theme('publication')   # inverted grayscale

fig1 = aplpy.FITSFigure(label['F555W_drzcroplinearShift'][0], \
        figure=figC, subplot=[x0,row_a,full_width,dy])
fig1.show_grayscale(stretch='log', vmin=-0.01869, vmax=121, vmid=-0.025)
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
fig1.show_contour(label_vla['C_R0'][0], colors="blue", levels=sigma_contour_tenpercent(Smax), linewidths=1.2)#, layer='fg')

# figred.show_contour(label_vla['red'][0], colors="lime", levels=sigma_contour_array(sigma_red), linewidths=2)

########################################
# beam
########################################
fig1.show_beam(major=0.49/3600, minor=0.35/3600, angle=0.18, edgecolor='white', facecolor='white', linestyle='solid', linewidth=3, frame=False, alpha=0.8)

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
markers_cross(fig1, ra_cross, dec_cross, layer='marker_set_1')
#markers_cross(figO, ra_cross, dec_cross, layer='marker_set_1')
# markers_cross(figHST, ra_cross, dec_cross, layer='marker_set_1')
# markers_cross(figred, ra_cross, dec_cross, layer='marker_set_1')

########################################
# Labels
########################################
# if '_' in sym[:-1]: symf = sym.replace('_', ' ')

# put_label(fig1, 0.20, 0.95, 'HST F555W', 'titleBand')
put_label(fig1, 0.38, 0.9, 'HST F555W, 5GHz Continuum', 'titleBand')
put_label(fig1, 0.15725, 0.95, 'RXJ1131', 'titleObj')
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
        figC.savefig(Plotpath + 'F555W_ContVLA.eps', dpi=600, bbox_inches='tight')
    else:
#        figC.canvas.draw()
        plt.show()
