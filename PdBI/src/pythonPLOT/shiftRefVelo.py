#!/usr/bin/env python
'''

Change HEADER CRPIX3 (channel) and write a new fits file, so that channel map labels can be symmetric across systemic velocity

Last Modified: 30 Apr 16

Author: Daisy Leung

History:
30 Apr 2016
- created, works with COchan.py

Note:
-

'''


import pywcs
import pyfits

final_image = '/Users/admin/Research/RXJ1131/PdBI/data/30Apr16/sup127_155_2ndcln_noCont.fits'
fits_cube = pyfits.open(final)
header = fits_cube[0].header
# Use channel 145 as 0 km/s, chosen using `go view` and also corresponds to z~0.654
header['CRPIX3'] = 145.0

# sanity check
wcs = pywcs.WCS(header)
x = wcs.sub([3])
x.wcs_pix2sky([[130]], 0)/1e3
x.wcs_pix2sky([[145]], 0)/1e3

outfile = '/Users/admin/Research/RXJ1131/PdBI/data/30Apr16/centralizedCube.fits'
fits_cube.writeto(outfile, clobber=True)