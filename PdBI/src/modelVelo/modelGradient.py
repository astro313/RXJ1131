'''


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

firstMom_obs = '/Users/admin/Research/RXJ1131/PdBI/data/30Apr16/centralizedCube4GILDAS-momCLIP5sigma.velo.fits'
interval = 50.    # km/s between velocity contour

a = pyfits.open(firstMom_obs)
hdr = a[0].header
wcs = pywcs.WCS(hdr).sub([1, 2])
plotwcs = wcsaxes.WCS(hdr, naxis=2)

fig = plt.figure(figsize=(12, 15), dpi=100)
plt.subplots_adjust(top=0.9)
ax = fig.add_subplot(1, 1, 1, projection=plotwcs)
ra, dec = ax.coords
ra.set_major_formatter('hh:mm:ss.s')
ra.display_minor_ticks(True)
dec.display_minor_ticks(True)
ra.set_minor_frequency(5)
ra.set_axislabel(" RA (J2000)", minpad=0.5)
dec.set_axislabel("DEC(degrees)", minpad=-0.4)
ra.set_separator((':', ':'))
ra.set_ticks(size=10, width=1.5)
dec.set_ticks(size=10, width=1.5)

cmap = plt.cm.jet
import matplotlib.colors as mcolors
norm = mcolors.Normalize()

dat = a[0].data.reshape(a[0].data.shape[2], a[0].data.shape[3])

# first moment map
im = ax.imshow(dat, origin="lower", norm=norm, cmap=cmap)
cont_list = [interval * i for i in range(-10, 10) if i != 0]
ax.contour(dat,
           cont_list,
           colors='black')

DecCentroid = -12.5328629
RACentroid = 172.96434563
LensRA_PX, LensDEC_PX = wcs.wcs_sky2pix(RACentroid, DecCentroid, 1)

p_list = [(-0.20, 0.12), (-0.02, 0.11), (0.20, -0.13), (0.29, -0.01), (0.57, -0.76), (0.87, -0.57), (1.04, -0.56), (1.32, -0.77)]
for p1 in p_list:
    ra_px, dec_px = p1
    l1, = ax.plot([ra_px+LensRA_PX], [dec_px+LensDEC_PX],
                  color='white',
                  marker='+',
                  mec="black", mew=1.25, ms=6,
                  zorder=3.1)    # lower zorder are drawn first

ax.set_xlim(110, 150)
ax.set_ylim(110, 140)
plt.show()
a.close()
