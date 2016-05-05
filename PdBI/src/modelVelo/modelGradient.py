'''

plot source locations from lens model of various channel as markers on observed 1st moment map

Last Modified: 04 May 16

History:
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
interval = 50.    # km/s between velocity contour
DecCentroid = -12.5328629
RACentroid = 172.96434563
# source positions from models of different channel, offset in arcsec
p_list = [(-0.20, 0.12), (-0.02, 0.11), (0.20, -0.13), (0.29, -0.01), (0.57, -0.76), (0.87, -0.57), (1.04, -0.56), (1.32, -0.77)]


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
for p1 in p_list:
    c = next(color)
    ra_offset_arcsec, dec_offset_arcsec = p1

    # convert into coord then pixel on image
    ra_px, dec_px = wcs.wcs_sky2pix(ra_offset_arcsec/3600.+RACentroid, dec_offset_arcsec/3600.+DecCentroid, 1)

    l1, = ax.plot([ra_px], [dec_px],
                  marker='+',
                  mec=c, mew=1.25, ms=8,
                  zorder=3.1)    # lower zorder are drawn first
#    plt.pause(0.05)

# plt.legend(['c{}'.format(i) for i in range(len(p_list))], loc=2, bbox_to_anchor=(1.05, 1), borderaxespad=0., fontsize=15)

plt.show()
a.close()

User_input = raw_input('Save figure? (Y/N): ')
if User_input == 'Y':
    filename = "veloGradient_markers.eps"
    fig.savefig(Plotpath + filename, dpi=100,
                bbox_inches="tight", pad_inches=0.1)
    print "-- Saved figure as : %s --" % (Plotpath + filename)