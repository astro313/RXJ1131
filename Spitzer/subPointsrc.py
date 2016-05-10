#!/usr/bin/env python

'''

Subtract foreground profile to extract fluxes


'''

import numpy as np
from astropy.io import fits
from astropy.modeling import models, fitting
import matplotlib.pyplot as plt
import pywcs, pyfits

# 3.6, 4.5, 5.8, and 8.0 um
IRACim1 = 'IRAC/ch1/pbcd/maic.fits'       # use Gauss
IRACim2 = 'IRAC/ch2/pbcd/maic.fits'     # use Gauss
IRACim3 = 'IRAC/ch3/pbcd/maic.fits'     # use Gauss
IRACim4 = 'IRAC/ch4/pbcd/maic.fits'     # use gauss

srcRadPix = 5.    # radius in px, within which we sum up fluxes; archive flux used aperture ~5.8"diameter --> ~5 px radius
# to follow aperture correction; http://irsa.ipac.caltech.edu/data/SPITZER/docs/irac/iracinstrumenthandbook/27/
# 12-20 in 1.2" px <-> 24 - 40 in 0.6" px
innerRadPix = 24
outerRadPix = 40
aperture_correction_4px = np.array([1.208, 1.220, 1.349, 1.554])   # based on R = 4px of .6"
# but I am using R = 5 px of .6", correction should be lower
# scale it down using aperture_coorection of 6px of .6"
aperture_correction_6px = np.array([1.112, 1.112, 1.118, 1.213])
# take average between 4px and 6px to get correction factor for R=5px
aperture_correction = (aperture_correction_4px + aperture_correction_6px)/2.

# to subtract PSF
starFWHM_ch1 = [3.51, 4.15, 3.21, 2.84, 3.61]
starFWHM_ch2 = [3.42, 3.1, 3.56, 3.19, 2.77]
starFWHM_ch3 = [3.37, 3.87, 3.87, 3.27, 3.67]
starFWHM_ch4 = [4.40, 4.01, 3.57, 3.79, 3.7]

# ---------------------------------------
# change this as necessary
FWHM = np.mean(starFWHM_ch1)
IRACim = IRACim1
aper_cor = aperture_correction[0]       # [channel # -1]
# ---------------------------------------

xstd = FWHM/2.355
ystd = xstd
# centroid position on Foreground, guided by HST image by matching crosshair
xcen = 1242         # 1241.24
ycen = 563          # 560.46
nsersic = 4
q = 0.95     # from C06
PA = -60*3600./206265     # radian
Reff = 1.136          # 8.33 kpc / (4.4 kpc/1") * (0.6"/1px)

# ------------------------------------------------------------------------
a = fits.open(IRACim)
b = a[0].data
header = pyfits.getheader(IRACim)
wcs = pywcs.WCS(header)
yy, xx = np.indices(b.shape)

g_init = models.Gaussian2D(b[ycen, xcen], xcen, ycen, xstd, ystd)
sersic_init = models.Sersic2D(b[ycen, xcen], Reff, nsersic, xcen, ycen, 1-q, PA)

fitter = fitting.LevMarLSQFitter()

g = fitter(g_init, xx, yy, b)
# g = fitter(sersic_init, xx, yy, b)

res_zoom = b[ycen-50:ycen+50, xcen-50:xcen+50] # (b - g(xx, yy))[ycen-50:ycen+50, xcen-50:xcen+50]
pos_res_zoom = np.ma.masked_where(res_zoom < 0., res_zoom)
# show over subtracting
# plt.figure(1)
# plt.imshow(b[ycen-50:ycen+50, xcen-50:xcen+50])
plt.figure(2)
plt.imshow(pos_res_zoom)
plt.show()

plt.figure(figsize=(10, 10))
plt.subplot(1, 3, 1)
plt.imshow(b[ycen-50:ycen+50, xcen-50:xcen+50], origin='lower', interpolation='nearest')
plt.title("Data")
plt.subplot(1, 3, 2)
plt.imshow(g(xx, yy)[ycen-50:ycen+50, xcen-50:xcen+50], origin='lower', interpolation='nearest')
plt.title("Model")
plt.subplot(1, 3, 3)
plt.imshow(res_zoom, origin='lower', interpolation='nearest')
plt.title("Residual")
plt.show()


data = pos_res_zoom   # after sub
xcen, ycen = int(len(data[1])/2), int(len(data[0])/2)  # since we subtracted fg & do aperturn pho, doesn't matter if we center on fg or bg

import math
sumSrc = 0.0
npixSrc = 0
srcDataLst = []

sumSky = 0.0
npixSky = 0
skyDataLst = []

# range to iterate in for loop
minXPix = xcen - outerRadPix
minYPix = ycen - outerRadPix
maxXpix = xcen + outerRadPix
maxYpix = ycen + outerRadPix

minXPix = int(math.floor(minXPix))
minYPix = int(math.floor(minYPix))
maxXPix = int(math.ceil(maxXpix))
maxYPix = int(math.ceil(maxYpix))

for i in range(minXPix, maxXPix + 1):
    for j in range(minYPix, maxYPix + 1):
        dist = np.sqrt(math.pow((i-xcen), 2)+math.pow((j-ycen), 2))

        # within circular aperture
        if dist <= srcRadPix:
            if not np.isnan(data[j-1, i-1]) and not (data[j-1, i-1] is np.ma.masked):
                sumSrc += data[j-1, i-1]
                npixSrc += 1
                srcDataLst.append(data[j-1, i-1])

        # annulus, sky
        if dist > innerRadPix and dist <= outerRadPix:
            if not np.isnan(data[j-1, i-1]):
                sumSky += data[j-1, i-1]
                npixSky += 1
                skyDataLst.append(data[j-1, i-1])

medianSky = np.median(skyDataLst)
meanSky = np.mean(skyDataLst)
stdevSky = np.std(skyDataLst)
maxSky = np.max(skyDataLst)
minSky = np.min(skyDataLst)

medianSrc = np.median(srcDataLst)
meanSrc = np.mean(srcDataLst)
stdevSrc = np.std(srcDataLst)
maxSrc = np.max(srcDataLst)
minSrc = np.min(srcDataLst)


def get_pixscale(header):

    xdelt = abs(header['PXSCAL1'])
    ydelt = abs(header['PXSCAL2'])
    return xdelt, ydelt


#-----------------------------------------------------------------------------#
# Calculate the conversion from MJy sr^-1 to flux density
#-----------------------------------------------------------------------------#
xpixScale, ypixScale = get_pixscale(header)        # in arcsec/px
pxScale = xpixScale*ypixScale                      # arcsec^2/px
conversion = pxScale * 1./206265**2 * 1e6 * 1e6    # Mjy/sr/px to uJy/px
print("Aperture correction Factor: {:8.5f}").format(aper_cor)

# integrated flux
# for IRAC images, unit = MJy sr^-1
integFlux = (sumSrc - (medianSky * npixSrc)) * conversion / 1e3
print '\n sum inside aperture: \t', sumSrc, ' Myr/sr \n'
print '\n sum inside aperture - annulus: \t', (sumSrc - (medianSky * npixSrc)), ' Myr/sr \n'
print '\n sum inside aperture - annulus: \t', integFlux, 'mJy \n'
integFlux_corrected = integFlux * aper_cor
print '\n aperture corrected flux inside aperture - annulus: \t', integFlux_corrected, 'mJy \n'

# uncertainty in the integrated flux
# F.Masci, IPAC: 'Flux-Uncertainty from Aperture Photometry'
dintegFlux = math.sqrt(npixSrc * math.pow(stdevSky, 2.0) *
                       (1.0 + float(npixSrc) / npixSky)) * conversion / 1e3 * aper_cor

a.close()

print("on host gal: {:.4f} +/- {:.4f}").format(integFlux_corrected, dintegFlux)

