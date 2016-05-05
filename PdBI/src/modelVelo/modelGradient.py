'''

plot source locations from lens model of various channel as markers on observed 1st moment map

Last Modified: 05 May 16

History:
05 May 16:
  - added error bars to positions
04 May 16:
  - created code


'''
import matplotlib
import matplotlib.pyplot as plt
import pyfits, pywcs
import wcsaxes
import numpy as np

font = {'family': 'Arial Narrow',
        'weight': 'bold',
        'size': 18}
matplotlib.rc('font', **font)

Plotpath = '/Users/admin/Research/RXJ1131/Figures/'
firstMom_obs = '/Users/admin/Research/RXJ1131/PdBI/data/30Apr16/centralizedCube4GILDAS-momCLIP5sigma.velo.fits'
interval = 50.        # km/s between velocity contour
DecCentroid = -12.5328629
RACentroid = 172.96434563
CellSize = 0.5        # arcsec per pixel in 1D

# source positions from models of different channel, offset in arcsec
p_list = [(-0.20, 0.12), (-0.02, 0.11), (0.20, -0.13), (0.29, -0.01), (0.57, -0.76), (0.87, -0.57), (1.04, -0.56), (1.32, -0.77)]
# corresponding error in arcsec
p_list_err = [(0.39, 0.3), (0.25, 0.2), (0.06, 0.21), (0.14, 0.22), (0.12, 0.29), (0.09, 0.15), (0.06, 0.12), (0.34, 0.45)]

# -----------------------------------------------
#  Read fits file and Set up WCS
# -----------------------------------------------
a = pyfits.open(firstMom_obs)
hdr = a[0].header
wcs = pywcs.WCS(hdr).sub([1, 2])
plotwcs = wcsaxes.WCS(hdr, naxis=2)
LensRA_PX, LensDEC_PX = wcs.wcs_sky2pix(RACentroid, DecCentroid, 1)
dat = a[0].data.reshape(a[0].data.shape[2], a[0].data.shape[3])

# -----------------------------------------------
#  Set up Figure
# -----------------------------------------------
fig = plt.figure(figsize=(12, 15), dpi=100)
plt.subplots_adjust(top=0.9)
ax = fig.add_subplot(1, 1, 1, projection=plotwcs)
ra, dec = ax.coords
ra.set_major_formatter('hh:mm:ss.s')
ra.display_minor_ticks(True)
dec.display_minor_ticks(True)
ra.set_minor_frequency(5)
ra.set_axislabel(" RA (J2000)", minpad=0.5)
dec.set_axislabel(" Dec (degrees)", minpad=-0.4)
ra.set_separator((':', ':'))
ra.set_ticks(size=10, width=1.5)
dec.set_ticks(size=10, width=1.5)
ax.set_xlim(120, 145)
ax.set_ylim(110, 140)

cmap = plt.cm.jet
import matplotlib.colors as mcolors
norm = mcolors.Normalize()

# -----------------------------------------------
# 1st moment map
# -----------------------------------------------
im = ax.imshow(dat, origin="lower", norm=norm, cmap=cmap)
cont_list = [interval * i for i in range(-10, 10) if i != 0]
ax.contour(dat,
           cont_list,
           colors='black')
ax.plot([LensRA_PX], [LensDEC_PX], marker='x', zorder=1, mec='k', ms=8)

from matplotlib.pyplot import cm
color = iter(cm.rainbow_r(np.linspace(0, 1, len(p_list))))
for i, p1 in enumerate(p_list):
    c = next(color)
    ra_offset_arcsec, dec_offset_arcsec = p1

    # convert into coord then pixel on image
    ra_px, dec_px = wcs.wcs_sky2pix(ra_offset_arcsec/3600.+RACentroid, dec_offset_arcsec/3600.+DecCentroid, 1)

    # convert arcsec to pixel offset on image
    ra_err = p_list_err[i][0]/CellSize
    dec_err = p_list_err[i][1]/CellSize
    ax.errorbar(x=[ra_px], y=[dec_px], xerr=[ra_err], yerr=[dec_err],
                    marker='+', ms=8, mec=c, mew=1.25, zorder=3.1, ecolor=c,
                    capsize=5, capthick=2, lw=2, label=str(i))

# matplotlib.rcParams['legend.handlelength'] = 0
# matplotlib.rcParams['legend.markerscale'] = 0
# ax.legend(loc=2, bbox_to_anchor=(0.9, 0.9), borderaxespad=0., fontsize=15, numpoints=1)

plt.show()
a.close()

User_input = raw_input('Save figure? (Y/N): ')
if User_input == 'Y':
    filename = "veloGradient_markers.eps"
    fig.savefig(Plotpath + filename, dpi=100,
                bbox_inches="tight", pad_inches=0.1)
    print "-- Saved figure as : %s --" % (Plotpath + filename)