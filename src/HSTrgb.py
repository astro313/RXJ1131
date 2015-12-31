'''
!!!!
31 Dec 2015: The astrometry of these HST images are not corrected
!!!

Use APLpy to generate RGB cube, and crop image if needed
Then generate a publication quality image


Last Modified: 24 Dec 2015

History:
--------


Note:
-----
- works for HST images for RXJ1131, need pre-process image to have just one extension using script: ../HST/src/cropimage.py
- another way to make RGB image is use numpy and append the array then imshow(), and the axes are generated using wcsaxes.WCS()
'''

import numpy as np
import os
import aplpy
import astropy
from astropy.io import fits
import pylab as pl
from APLpySetup import *

# since the HST files have multiple extensions, and are large --> crop and make it to just one extension before running this script

path = '../HST/'
Plotpath = '../Figures/'
rgbfile = 'hst_rgb_cube.fits'
ra_center = 172.96418
dec_center = -12.5330955

# if we need to crop image, for this set, croped
crop = False

if not os.path.isfile(os.path.join(path, rgbfile)):
    # ordered R->G->B
    aplpy.make_rgb_cube([
                         path+'HST_9744_41_NIC_NIC2_F160W_drzcrop.fits',
                         path+'HST_9744_75_ACS_WFC_F555W_drzcrop.fits',
                         path+'HST_9744_74_ACS_WFC_F814W_drzcrop.fits',
                        ], os.path.join(path, rgbfile))

rgb_cube = os.path.join(path, rgbfile)
f = fits.open(rgb_cube)
data = f[0].data
header = f[0].header
wcs = astropy.wcs.WCS(header)

if crop:
    ra_center_px, dec_center_px = 1942, 2406

    cx, cy = ra_center_px, dec_center_px
    # how large of a crop
    wx, wy = 70, 70

    cutout = fits.PrimaryHDU(data=data[:, cy-wy/2:cy+wy/2, cx-wx/2:cx+wx/2],
                             header=wcs[:, cy-wy/2:cy+wy/2, cx-wx/2:cx+wx/2].to_header())
    outfile = 'hst_rgb_cutout.fits'
    cutout.writeto(os.path.join(path, outfile), clobber=True)
    rgb_cube = os.path.join(path, outfile)


# save RGB image
aplpy.make_rgb_image(rgb_cube, rgb_cube.replace('.fits', '_nowcs.png'),
                     stretch_b='log', # vmin_b=-1.060e-02, vmax_b=2.320, vmid_b=0.001,
                     stretch_g='log', # vmin_g=-1.067e-02, vmax_g=6.028e-1, vmid_g=0.001,
                     stretch_r='log', # vmin_r=-1.730e-02, vmax_r=7.717e-01, vmid_r=0.001,
                    )

# Display RGB
pl.figure(1).clf()
# first get the axes setup with wcs via the original HST image
F = aplpy.FITSFigure(path+'HST_9744_75_ACS_WFC_F555W_drzcrop.fits', figure=pl.figure(1, figsize=(12, 8)))
F.show_rgb(rgb_cube.replace('.fits', '_nowcs.png'))
F.recenter(ra_center, dec_center, radius=0.0025)  # small radius = zoom in
markers_cross(F, ra_center, dec_center, ec='red')
put_label(F, 0.5, 0.95, 'HST F160W, F555W, F814W', 'titleBand')


if __name__ == '__main__':
    """
    run script.py True
    sys.argv[1] determines whether to save all the plots or not
    """
    import sys
    import os
    if len(sys.argv) < 2:
        errmsg = "Invalid number of arguments: {0:d}\n  run script.py Save_True"
        raise IndexError(errmsg.format(len(sys.argv)))
    saveFig = True if sys.argv[1].lower() == 'true' else False
    if saveFig:
        print "saved at: ", (Plotpath + rgbfile.replace('.fits', '.eps'))
        F.save(Plotpath + rgbfile.replace('.fits', '.eps'), dpi=600)
    else:
        pl.show()
