'''
Shift PdBI image to desire WCS

Last Modified: 25 Dec 2015


History:
--------
25 Dec 2015: copied from RXJ1131/HST/src/cropfits.py, modify with astropy (cropfits.py handles wcs in HST images. astropy.wcs.WCS deals better with radio line cube)


Note:
-----
- for GILDAS go spectrum to have (0,0) offset at the desire wcs -> need to change header as well, not just shifting the image to center on that wcs, which only cause the central pixel of the image to land on the desire WCS. But the (0,0) in go spectrum corresponds to HEADER

'''

#!/usr/bin/env python
import numpy as np
from os.path import join
from astropy.io import fits
import astropy
from astropy.wcs import WCS

ext = 0

PATH = '/Users/admin/Research/RXJ1131/PdBI/data/25Dec15'
filename = 'sup127_155_2ndcln_noCont.fits'

ra_center = 172.96418
dec_center = -12.5330955
size = 50  # px, one side


def load_fits_image(filename, xDeg, yDeg, ext=0):
    f = fits.open(filename)
    header = f[ext].header
    wcs = WCS(header, naxis=2)
    [[xPix, yPix]] = wcs.wcs_world2pix([(xDeg, yDeg)], 1)
    print xPix, yPix
    return xPix, yPix

def cropfits(dirname, xctr, yctr, xbdr, ybdr, infits=None, extension=0):
    import os
    files = [infits]

    for i, fname in enumerate(files):
        name, ext = os.path.splitext(fname)
        shortfiles = name + 'crop' + ext
        f = fits.open(fname)
#        print f
        header = f[0].header.copy()
        if extension != 0:
            header = f[extension].header.copy()

        cd1 = header.get('CDELT1') if header.get(
            'CDELT1') else header.get('CD1_1')
        cd2 = header.get('CDELT2') if header.get(
            'CDELT2') else header.get('CD2_2')
        if cd1 is None or cd2 is None:
            raise Exception("Missing CD or CDELT keywords in header")
        wcs = WCS(header)
        im = f[extension].data

        print(" Note the image dimensions .. ")
        print im.shape
        print ("should be same as header.. ")
        print str(header.get('NAXIS1')) + ' times ' + str(header.get('NAXIS2'))+ '; nchan: ' + str(header.get('NAXIS3'))
        print("")

        # crop region
        ymin, ymax = yctr - ybdr, yctr + ybdr
        xmin, xmax = xctr - xbdr, xctr + xbdr
        New_xctr = int((xmax-xmin)/2)
        New_yctr = int((ymax-ymin)/2)

        # shift the ref. pix s.t. the ref. wcs corresponds to where this pixel
        # is w.r.t to the cropped image, where xmin --> 0, ymin --> 0
        header['CRPIX1'] -= xmin
        header['CRPIX2'] -= ymin

        # Also need to do the following for GILDAS go spectrum to have (0,0) offset on the desire WCS
        # shift the ref. point to new x, y ctr that we are interested in.
        header['CRPIX1'] = New_xctr
        header['CRPIX2'] = New_yctr
        header['CRVAL1'] = ra_center
        header['CRVAL2'] = dec_center

        header['NAXIS1'] = int(xmax - xmin)
        header['NAXIS2'] = int(ymax - ymin)

        # always row, column
        im = im[0][:, ymin:ymax, xmin:xmax]
        print im.shape
        cutout = fits.PrimaryHDU(data=im, header=header)
        outname = join(dirname, shortfiles)
        cutout.writeto(outname, clobber=True)
        print "Cut image %s with dims %s to %s.  xrange: %f:%f, yrange: %f:%f \n" % (fname, f[extension].data.shape, im.shape, xmin, xmax, ymin, ymax)
        print outname, 'successfully cropped!'

xctr, yctr = load_fits_image(
    join(PATH, filename), ra_center, dec_center, ext=ext)
cropfits(PATH, int(round(xctr)), int(round(yctr)), size,
         size, infits=join(PATH, filename), extension=ext)

