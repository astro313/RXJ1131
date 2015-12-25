'''
Last Modified: 25 Dec 15


Author: Daisy Leung


History:
24 Dec 15
- added compass
- added super title
- show markers
- add beam

09 Dec 15
- works with CO(2-1) cube
- added contours overlay on channel map
- will zoom on image with correct axes
- only show axes label in one panel


Note:
- depends on package - astrolib.coords for putting up markers based on Coords
-

'''

import matplotlib.pyplot as plt
import pywcsgrid2
import pywcs
import pyfits
import numpy as np

import mpl_toolkits.axes_grid1.axes_grid as axes_grid
#from mpl_toolkits.axes_grid.colorbar import colorbar
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
from mpl_toolkits.axes_grid1.anchored_artists import AnchoredText

class Velo(object):
    def __init__(self, header):
        wcs = pywcs.WCS(header)
        self.wcs_vel = wcs.sub([3])        # velocity is 3rd axis

    def to_vel(self, p):
        v = self.wcs_vel.wcs_pix2sky([[p]], 0)
        return v[0][0]


def setup_axes(fig, header, nx, ny):

    gh = pywcsgrid2.GridHelper(wcs=header)
    gh.locator_params(nbins=4)

    g = axes_grid.ImageGrid(fig, 111,
                            nrows_ncols=(ny, nx),
                            ngrids=None,
                            direction='row',
                            # axes_pad=0.02,
                            axes_pad = 0.0,    # between panels
                            add_all=True,
                            share_all=True,
                            # if share_all = False, need the following
                            # share_x = False,
                            # share_y = False,
                            aspect=True,
                            label_mode='1',  # 'L'
                            cbar_mode=None,
                            axes_class=(pywcsgrid2.Axes, dict(grid_helper=gh)))

    # only show xylabel at that panel, depending on nrows_ncols in axes_grid.ImageGrid()
    # if comment out the following block, will shorten to Greek Letters
    ax = g[-1*nx]     # display in the bottom left corner only
    ax.set_xlabel("Right Ascension (J2000)")
    ax.set_ylabel("Declination (J2000)")

    # make colorbar
    ax = g[-1]
    cax = inset_axes(ax,
                     width="8%", # # width = 8% width = 10% of parent_bbox width
                     height="100%", # height : 50%
                     loc=3,    # 4
                     bbox_to_anchor=(1.01, 0, 1, 1),
                     bbox_transform=ax.transAxes,
                     borderpad=0.
                     )

    # super title
    at = AnchoredText(
                     "RXJ 1131-1231 CO(2-1) Channel Map",
                     loc=3,
                     prop=dict(size=24),
                     frameon=False,
                     bbox_to_anchor=(0.05, 0.925, 0.75, 1),
                     bbox_transform=ax.transAxes
                     # bbox_transform=gcf().transFigure
                     )
    g[-4].add_artist(at)

    return g, cax


fits_cube = pyfits.open("/Users/admin/Research/RXJ1131/PdBI/data/04Sep15/sup127_155_2ndcln_noCont.fits")
dxy = 125                   # starting channel
nx = 6
ny = 5
nxy = nx * ny               # total number of channels
sigma = 1.451e-3            # sigma in map for contour plots
# mark HST knots location on channel map
HST = True


header = fits_cube[0].header
cdelt1 = np.abs(header['CDELT1'] * 3600)       # arcsec/pix
cdelt2 = np.abs(header['CDELT2'] * 3600)
major_px = header['BMAJ'] * 3600 / cdelt1                 # arcsrc/pix
minor_px = header['BMIN'] * 3600 / cdelt2
BPA_deg = header['BPA']         # should be 13 deg

wcs = pywcs.WCS(header).sub([1, 2])       # .sub([1,2]) needed for wcs.wcs_pix2sky(px, py, 1) later in the code
vel = Velo(header)
fig = plt.figure(1, figsize=(12, 12), dpi=100)      # don't change dpi, disastrous
g, cax = setup_axes(fig, header, nx, ny)

# draw images
i = 0
# cmap = plt.cm.gray_r
cmap = plt.cm.jet
import matplotlib.colors as mcolors
norm = mcolors.Normalize()
images = []
start_channel = i*nxy+dxy

for i, ax in enumerate(g):
    channel_number = start_channel + i
    channel = fits_cube[0].data[0][channel_number]
    im = ax.imshow(channel, origin="lower", norm=norm, cmap=cmap)
    ax.contour(channel, [4*sigma, 8*sigma, 12*sigma, 16*sigma, 20*sigma],
               colors='black')

    # zoom-in on pixels
    ax.set_xlim(100, 150)
    ax.set_ylim(100, 150)


    if HST:   # mark lensing knots
        import astrolib.coords as coords
        p_list = [coords.Position("11:31:51.3  -12:31:59").j2000(),
              coords.Position("11:31:51.6  -12:31:59.6").j2000()]
        for p1 in p_list:
            ra, dec = p1
            ra_px, dec_px = wcs.wcs_sky2pix(ra, dec, 1)
            l1, = ax.plot([ra_px], [dec_px],
                          "k+", mec="black", mew=1.5, ms=8,
                          zorder=3.1)    # lower zorder are drawn first
    images.append(im)


if HST:
   # put up legend on last panel --> HST markers
    g[0].legend([l1],
                ["HST QSO lensing knots"],
                loc=0,
                numpoints=1,
                handlelength=0.75,
                frameon=False,
                prop=dict(size=8),
                fancybox=True)


# label with velocities
use_path_effect = False         # Fancy text
try:
    from matplotlib.patheffects import withStroke
except ImportError:
    use_path_effect = False

for i, ax in enumerate(g):
    channel_number = start_channel + i
    v = vel.to_vel(channel_number) / 1.e3
    t = ax.add_inner_title(r"$v=%4.1f$ km s$^{-1}$" % (v),
                           loc=2,
                           frameon=False
                           )
    t.patch.set_alpha(1.0)
    if use_path_effect:
        t.txt._text.set_path_effects([withStroke(foreground="w",
                                                 linewidth=3)])
    else:
        # show velocity text in some color
        t.txt._text.set_color('k')

# Beam
g[-1].add_beam_size(minor_px, major_px, BPA_deg,      # y, x, angle
                    loc=4, frameon=True,
                    patch_props=dict(facecolor='gray', alpha=0.8, edgecolor='k', fill=True))
# g[-1].add_compass(loc=3)

# make colorbar
cb = plt.colorbar(im, cax=cax)
cb.set_label("Flux Density [mJy B"+r"$^{-1}$]")
cb.set_ticks([0, 0.01, 0.02, 0.03])

# adjust norm
# norm.vmin = -0.007      # min in map
norm.vmin = -0.015        # show more contrast
norm.vmax = 0.03          # max in map
for im in images:
    im.changed()


plt.show()
fits_cube.close()

if 0:
    plt.savefig("co_channel_maps.eps", dpi=70, bbox_inches="tight")