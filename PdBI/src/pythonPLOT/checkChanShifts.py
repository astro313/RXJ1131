'''

Check by how many channels is shifted when we center 0km/s as the systemic velocity. It was chan 130, now chan 145 (systemic velocity) correspond to 0 km/s.

Last modified: 04 May 2016

History:
04 May 2016
    - created file


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
x.wcs_pix2sky([[130]], 0)/1e3        # corresponding velocity
x.wcs_pix2sky([[145]], 0)/1e3
x.wcs_sky2pix([[0]], 0)              # what channel corresponds to 0 kms

cube2 = pyfits.open(after)
header2 = cube2[0].header
print("reference channel is: {} for {} km/s").format(header2['CRPIX3'], header2['CRVAL3'])
