'''
Last Modified: 25 Dec 2015


History:
--------
25 Dec 2015: fix minor bugs
24 Dec 2015: cropped image has correct wcs


Note:
-----
using astropy wcs is simplier than pywcs to write new header to fits
    -> can just parse a cropped header to PrimaryHDU()
- but this set of HST fits file is partciularly difficult to dealt with because of the multiple extensions, and the header contains characters that astropy.wcs doesn't like
'''

#!/usr/bin/env python
import pyfits as pf
import pywcs
import numpy as np
from os.path import join
from sys import argv

# for the HST images of this object
ext = 1

PATH = '/Users/admin/Research/RXJ1131/HST/'
filename = 'HST_9744_75_ACS_WFC_F555W_drz.fits'
file2 = 'HST_9744_74_ACS_WFC_F814W_drz.fits'
file3 = 'HST_9744_41_NIC_NIC2_F160W_drz.fits'
ra_center = 172.96418  # 172.96426
dec_center = -12.5330955  # -12.533213
# ra_center_px = 1942
# dec_center_px = 2406
#
size = 200  # px, one side


def load_fits_image(filename, xDeg, yDeg, ext=0):
    header = pf.getheader(filename, ext=ext)
    wcs = pywcs.WCS(header)
    [[xPix, yPix]] = wcs.wcs_sky2pix([(xDeg, yDeg)], 1)
    print xPix, yPix
    return xPix, yPix

def cropfits(dirname, xctr, yctr, xbdr, ybdr, infits=None, extension=0):
    import os
    files = [infits]

    for i, fname in enumerate(files):
        name, ext = os.path.splitext(fname)
        shortfiles = name + 'crop' + ext
        hdulist = pf.open(fname)
#        print hdulist
        header = hdulist[0].header.copy()
        if extension != 0:
            header = hdulist[extension].header.copy()

        if header['NAXIS'] > 2:
            raise DimensionError("Too many (%i) dimensions!" % header['NAXIS'])
        cd1 = header.get('CDELT1') if header.get(
            'CDELT1') else header.get('CD1_1')
        cd2 = header.get('CDELT2') if header.get(
            'CDELT2') else header.get('CD2_2')
        if cd1 is None or cd2 is None:
            raise Exception("Missing CD or CDELT keywords in header")
        wcs = pywcs.WCS(header)
        im = hdulist[extension].data

        print(" Note the image dimensions .. ")
        print im.shape
        print ("should be same as header.. ")
        print str(header.get('NAXIS1')) + ' times ' + str(header.get('NAXIS2'))
        print("")

        # crop region
        ymin, ymax = yctr - ybdr, yctr + ybdr
        xmin, xmax = xctr - xbdr, xctr + xbdr
        New_xctr = int((xmax - xmax)/2)
        New_yctr = int((y_max - ymin)/2)

        # shift the ref. pix s.t. the ref. wcs corresponds to where this pixel
        # is w.r.t to the cropped image, where xmin --> 0, ymin --> 0
        header['CRPIX1'] -= xmin
        header['CRPIX2'] -= ymin

        # Another way to do this might be ....
        # shift the ref. point to new x, y ctr that we are interested in.
        # header['CRPIX1'] = New_xctr
        # header['CRPIX2'] = New_yctr
        # header['CRVAL1'] = ra_center
        # header['CRVAL2'] = dec_center

        header['NAXIS1'] = int(xmax - xmin)
        header['NAXIS2'] = int(ymax - ymin)

        # always row, column
        im = im[ymin:ymax, xmin:xmax]
        print im.shape
        cutout = pf.PrimaryHDU(data=im, header=header)
        outname = join(dirname, shortfiles)
        cutout.writeto(outname, clobber=True)
        print "Cut image %s with dims %s to %s.  xrange: %f:%f, yrange: %f:%f \n" % (fname, hdulist[extension].data.shape, im.shape, xmin, xmax, ymin, ymax)
        print outname, 'successfully cropped!'

#F555W
xctr, yctr = load_fits_image(
    join(PATH, filename), ra_center, dec_center, ext=ext)
cropfits(PATH, int(round(xctr)), int(round(yctr)), size,
         size, infits=join(PATH, filename), extension=ext)

#F814W
xctr, yctr = load_fits_image(
    join(PATH, file2), ra_center, dec_center, ext=ext)
cropfits(PATH, int(round(xctr)), int(round(yctr)), size,
         size, infits=join(PATH, file2), extension=ext)
#F160W
xctr, yctr = load_fits_image(
    join(PATH, file3), ra_center, dec_center, ext=ext)
cropfits(PATH, int(round(xctr)), int(round(yctr)), size,
         size, infits=join(PATH, file3), extension=ext)
