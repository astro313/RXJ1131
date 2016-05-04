'''

Check by how many channels is shifted when we center 0km/s as the systemic velocity. It was chan 130, now chan 145 (systemic velocity) correspond to 0 km/s.
Also to find out what are the velocity ranges to use for making moment maps in GILDAS.

Last modified: 04 May 2016

History:
04 May 2016
    - created file

Note:
04 May 2016: was missing emission if just use chan 126-155 to make moment maps --> will use 124-156


'''
import pywcs
import pyfits

before = '/Users/admin/Research/RXJ1131/PdBI/data/30Apr16/sup127_155_2ndcln_noCont.fits'
after = '/Users/admin/Research/RXJ1131/PdBI/data/30Apr16/centralizedCube.fits'

fits_cube = pyfits.open(before)
header = fits_cube[0].header
print("reference channel was: {} for {} km/s").format(header['CRPIX3'], header['CRVAL3'])

wcs = pywcs.WCS(header)
x = wcs.sub([3])
print x.wcs_pix2sky([[130]], 0)/1e3        # corresponding velocity
print x.wcs_pix2sky([[145]], 0)/1e3
print x.wcs_sky2pix([[0]], 0)              # what channel corresponds to 0 kms
print

cube2 = pyfits.open(after)
header2 = cube2[0].header
print("reference channel is: {} for {} km/s").format(header2['CRPIX3'], header2['CRVAL3'])

# channel ranges for high-O maps:
first, last = 124, 156
# their velocity ranges:
wcs = pywcs.WCS(header2)
x = wcs.sub([3])
x.wcs_pix2sky([[first]], 0)/1e3        # corresponding velocity
x.wcs_pix2sky([[last]], 0)/1e3

