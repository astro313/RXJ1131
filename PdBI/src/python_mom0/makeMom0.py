#!/usr/bin/env python
'''

look through histogram again for each channel

add the channels manually to create moment-0 map, and look at historgram after

since the noise properties look ok, make a FITS file using this with header from the GILDAS moment0 map (which can potentially cause trouble in the future, but run with it for now) for APLpy plot.


Last modified: May 16 2016

History:
--------
16 May 2016:
    - created script, working fine
'''

import matplotlib
import matplotlib.pyplot as plt
import matplotlib.mlab as mlab
import pywcs
import pyfits
import numpy as np
from scipy.stats import norm


fits_cube = pyfits.open("/Users/admin/Research/RXJ1131/PdBI/data/15May16/centralizedCube4GILDAS.fits")

header = fits_cube[0].header
cdelt1 = np.abs(header['CDELT1'] * 3600)       # arcsec/pix
cdelt2 = np.abs(header['CDELT2'] * 3600)
major_px = header['BMAJ'] * 3600 / cdelt1                 # arcsrc/pix
minor_px = header['BMIN'] * 3600 / cdelt2
BPA_deg = header['BPA']         # should be 13 deg

# check histogram of each slice
start_channel = 124 - 1        # python indexing from 0
nchan = 157-124+1
for i in range(nchan):
    channel_number = start_channel + i
    channel = fits_cube[0].data[0][channel_number]
    n, bins, patches = plt.hist(channel[15:220, 15:100].flatten(), bins=100, facecolor='green', alpha=0.75, normed=True)
    # best fit of data
    (mu, sigma) = norm.fit(channel[15:220, 15:100])
    # add a 'best fit' line
    y = mlab.normpdf(bins, mu, sigma)
    l = plt.plot(bins, y, 'r--', linewidth=2)
    #plot
    plt.xlabel('flux Jy per channel')
    plt.ylabel('counts')
    plt.title(r'$\mathrm{Histogram\ of\ offregion for chan:%i :}\ \mu=%.3f,\ \sigma=%.1f$ mJy' %(channel_number, mu, sigma*1e3))
    plt.grid(True)
    plt.show()
    raw_input('')
    plt.close()


# add them and check histo
mom0_holder = np.zeros_like(fits_cube[0].data[0][0])
for i in range(nchan):
    channel_number = start_channel + i
    mom0_holder += fits_cube[0].data[0][channel_number]
    print mom0_holder.mean()
    raw_input(" ")

n, bins, patches = plt.hist(mom0_holder.flatten(), bins=100, facecolor='green', alpha=0.75, normed=True)
(mu, sigma) = norm.fit(mom0_holder[15:220, 15:100])
y = mlab.normpdf(bins, mu, sigma)
l = plt.plot(bins, y, 'r--', linewidth=2)
plt.show()
# --> looks expected


cmap = plt.cm.jet
import matplotlib.colors as mcolors
norm = mcolors.Normalize()
images = []

fig, ax = plt.subplots()
fig.subplots_adjust(top=0.9)
im = ax.imshow(mom0_holder, origin="lower", norm=norm, cmap=cmap)
cont_list = [sigma * i for i in range(-6, 24, 3) if i != 0]
ax.contour(mom0_holder, cont_list, colors='black')
plt.show()
# --> save this to a new fits file for APLpy plot

