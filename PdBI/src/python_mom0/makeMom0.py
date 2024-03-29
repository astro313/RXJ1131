#!/usr/bin/env python
'''

look through histogram again for each channel

add the channels manually to create moment-0 map, and look at historgram after

since the noise properties look ok, make a FITS file using this with header from the GILDAS moment0 map for APLpy plot. --> unfortunately can't cheat APLpy with the GILDAS header and I don't know how to edit it appropriately.. Use instead CASA to make the moment-0 map, see compare_mom0_CASA.py

Last modified: May 16 2016

History:
--------
16 May 2016:
    - updated mom0 value to multiplied by chan width
    - fix bug, write out new FITS using the highest SNR array
    - look for ranges of channel that gives highest SNR
    - created script, working fine


Note:
-----
    - useful to investigate the noise distribution per channel and if also see it in the moment-0 map, like in GILDAS
    - useful to find the chan range with highest SNR
    - but can't use to plot with APLpy, because header is not consistent with data
    - altho. can use to plot with Matplotlib, but I can't figure out how to match coord. sys. with the HST image
'''

import matplotlib
import matplotlib.pyplot as plt
import matplotlib.mlab as mlab
from os.path import join
import pyfits
import numpy as np
from scipy.stats import norm

font = {'family': 'Arial Narrow',
        'weight': 'bold',
        'size': 18.75}
matplotlib.rc('font', **font)

specRes = 21.528154         # km/s
path = '/Users/admin/Research/RXJ1131/PdBI/data/15May16/'
fits_cube = pyfits.open(join(path, 'centralizedCube4GILDAS.fits'))

header = fits_cube[0].header
cdelt1 = np.abs(header['CDELT1'] * 3600)       # arcsec/pix
cdelt2 = np.abs(header['CDELT2'] * 3600)
major_px = header['BMAJ'] * 3600 / cdelt1                 # arcsrc/pix
minor_px = header['BMIN'] * 3600 / cdelt2
BPA_deg = header['BPA']         # should be 13 deg

# check histogram of each slice, offsource region
start_channel = 127 - 1        # python indexing from 0
nchan = 160 - start_channel + 1
_plot = False

for i in range(nchan):
    channel_number = start_channel + i
    channel = fits_cube[0].data[0][channel_number]
    # best fit of data
    (mu, sigma) = norm.fit(channel[15:220, 15:100])
#     print mu, sigma
    if _plot:
        n, bins, patches = plt.hist(channel[15:220, 15:100].flatten(), bins=100, facecolor='green', alpha=0.75, normed=True)
        # add a 'best fit' line
        y = mlab.normpdf(bins, mu, sigma)
        l = plt.plot(bins, y, 'r--', linewidth=2)

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
    mom0_holder += fits_cube[0].data[0][channel_number] * specRes
(mu, sigma) = norm.fit(mom0_holder[15:220, 15:100].flatten())
print mu, sigma

n, bins, patches = plt.hist(mom0_holder[15:220, 15:100].flatten(), bins=100, facecolor='green', alpha=0.75, normed=True)
y = mlab.normpdf(bins, mu, sigma)
l = plt.plot(bins, y, 'r--', linewidth=2)
plt.show()


# check by summing different number of channels
sigma_ch = np.std(fits_cube[0].data[0][10:120, :, :].flatten())      # in Jy
# from GILDAS 124-156 == Python 123-156
start_ = 123
end_ = 156
plt.hist(fits_cube[0].data[0][start_:end_, 15:220, 15:100].sum(axis=0).flatten()*specRes, bins=100, normed=True)
(mu, sigma) = norm.fit(fits_cube[0].data[0][start_:end_, 15:220, 15:100].sum(axis=0).flatten()*specRes)
print (" from array {}").format(fits_cube[0].data[0][start_:end_, 15:220, 15:100].sum(axis=0).flatten().std()*specRes)
print (" from bestfit {}").format(sigma)
print (" expected: {}").format(sigma_ch * np.sqrt(end_-start_)*specRes)
# --> some channels have higher RMS

cmap = plt.cm.jet
import matplotlib.colors as mcolors
normC = mcolors.Normalize()

fig, ax = plt.subplots()
fig.subplots_adjust(top=0.9)
im = ax.imshow(mom0_holder, origin="lower", norm=normC, cmap=cmap)
cont_list = [sigma * i for i in range(-6, 30, 3) if i != 0]
ax.contour(mom0_holder, cont_list, colors='black')
plt.show()

# ------------------------------------------------------------------
# look for channel range that gives the highest SNR
# weird indexing because python start from 0, and also index slicing is up to n-1
# case 1: 124-156 <=> python 123 - 156; SNR = 26.33
# case 2: 124-157 <=> python 123 - 157; SNR = 26.82
# case 3: 127-155 <=> python 126 - 155; SNR = 24.22
# case 4: 126-160 <=> python 125 - 160; SNR = 27.48 with sigma = 0.013731 Jy/B

start_chan = 125
end_chan = 160

# only need to test SNR in a small region
mom0_SNR = np.zeros_like(fits_cube[0].data[0][0][110:145, 110:145])
for i in range((end_chan-start_chan+1)):
    channel_number = start_chan + i
    mom0_SNR += fits_cube[0].data[0][channel_number][110:145, 110:145] * specRes
print ("max SNR: {:.2f}").format(mom0_SNR.max()/(fits_cube[0].data[0][start_chan:end_chan, 15:220, 15:100].sum(axis=0).flatten().std()*specRes))
print (" map STD {} ").format(fits_cube[0].data[0][start_chan:end_chan, 15:220, 15:100].sum(axis=0).flatten().std()*specRes)

fig, ax = plt.subplots()
fig.subplots_adjust(top=0.9)
im = ax.imshow(mom0_SNR, origin="lower", norm=normC, cmap=cmap)
cont_list = [fits_cube[0].data[0][start_chan:end_chan, 15:220, 15:100].sum(axis=0).flatten().std()*specRes * i for i in range(-6, 30, 3) if i != 0]
ax.contour(mom0_SNR, cont_list, colors='black')
plt.show()

# the full image of the highest SNR range
final = fits_cube[0].data[0][start_chan:end_chan, :, :].sum(axis=0) * specRes


# ------------------------------------------------------------------
# --> save this to a new fits file for APLpy plot
# get header from GILDAS moment0 file

fits_gildas_mom0 = pyfits.open(join(path, "centralizedCube4GILDAS-mom0.fits"))
hdr = fits_gildas_mom0[0].header
# patch up header, somehow info is missing in the GILDAS creater header
hdr['CTYPE3'] = header['CTYPE3']
hdr.set('RESTFREQ', header['RESTFREQ'], 'patched from original cube')
hdr.set('VELO-LSR', header['VELO-LSR'], 'patched from original cube')
hdr.set('VELREF', header['VELREF'], 'patched from original cube')
hdr.set('SPECSYS', header['SPECSYS'], 'patched from original cube')
hdr.add_comment("This Moment-0 map is created with python, but I stole the header from GILDAS the moment0 map.")
new_mom0 = 'centralizedCube4GILDAS-python_ch' + str(start_chan+1) + '-' + str(end_chan) + '_mom0.fits'
pyfits.writeto(join(path, new_mom0), final, hdr, clobber=True)
fits_cube.close()


