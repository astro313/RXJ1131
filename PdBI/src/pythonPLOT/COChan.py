'''
Last Modified: 27 May 16


Author: Daisy Leung


History:
27 May 2016
- change flux label unit
- update font inside panel

30 Apr 2016
- change velocity labels, reference systemic velocity
- use data cube with velo. axis 0 km/s == line center (see shiftRefVelo.py)

27 Mar 2016
- change symbol for foreground galaxy
- increment contour by 3sigma instead of 4 sigma
- add channels before first and after last
- don't put legend for markers (don't fit in there now with more panels)
- change colorbar ticks to mJy B^-1
- change xaxis tick frequency to avoid overlap labeling

31 Dec 2015
- update HST coordinates, see HSTmarkers.txt

25 Dec 2015
- Markers corresponds to HST features
- removed super title code.. doesn't work with saved figure

24 Dec 2015
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

'''
import matplotlib
import matplotlib.pyplot as plt
import pywcsgrid2
import pywcs
import pyfits
import numpy as np

import mpl_toolkits.axes_grid1.axes_grid as axes_grid
#from mpl_toolkits.axes_grid.colorbar import colorbar
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
from mpl_toolkits.axes_grid1.anchored_artists import AnchoredText

font = {'family': 'Arial Narrow',
        'weight': 'bold',
        'size': 8.75}
matplotlib.rc('font', **font)


class Velo(object):

    def __init__(self, header):
        wcs = pywcs.WCS(header)
        self.wcs_vel = wcs.sub([3])        # velocity is 3rd axis

    def to_vel(self, p):
        v = self.wcs_vel.wcs_pix2sky([[p]], 0)
        return v[0][0]


def setup_axes(fig, header, nx, ny):

    gh = pywcsgrid2.GridHelper(wcs=header)
    gh.locator_params('x', nbins=2)

    g = axes_grid.ImageGrid(fig, 111,
                            nrows_ncols=(ny, nx),
                            ngrids=None,
                            direction='row',
                            # axes_pad=0.02,
                            axes_pad=0.0,    # between panels
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
    ax = g[-1 * nx]     # display in the bottom left corner only
    ax.set_xlabel("Right Ascension (J2000)")
    ax.set_ylabel("Declination (J2000)")
    # other attributes:
    # http://matplotlib.org/mpl_toolkits/axes_grid/api/axis_artist_api.html

    # make colorbar
    ax = g[-1]
    cax = inset_axes(ax,
                     width="8%",  # width = 8% width = 10% of parent_bbox width
                     height="100%",  # height : 50%
                     loc=3,    # 4
                     bbox_to_anchor=(1.01, 0, 1, 1),
                     bbox_transform=ax.transAxes,
                     borderpad=0.
                     )

    return g, cax

Plotpath = '/Users/admin/Research/RXJ1131/Figures/'
# fits_cube = pyfits.open("/Users/admin/Research/RXJ1131/PdBI/data/04Sep15/sup127_155_2ndcln_noCont.fits")
fits_cube = pyfits.open("/Users/admin/Research/RXJ1131/PdBI/data/30Apr16/centralizedCube.fits")

dxy = 124                   # starting channel
nx = 7
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

# .sub([1,2]) needed for wcs.wcs_pix2sky(px, py, 1) later in the code, flattens cube axes
wcs = pywcs.WCS(header).sub([1, 2])
vel = Velo(header)
fig = plt.figure(figsize=(12, 15), dpi=100)      # don't change dpi, disastrous
plt.subplots_adjust(top=0.9)
g, cax = setup_axes(fig, header, nx, ny)

# draw images
i = 0
# cmap = plt.cm.gray_r
cmap = plt.cm.jet
import matplotlib.colors as mcolors
norm = mcolors.Normalize()
images = []
start_channel = i * nxy + dxy

for i, ax in enumerate(g):
    channel_number = start_channel + i
    channel = fits_cube[0].data[0][channel_number]
    im = ax.imshow(channel, origin="lower", norm=norm, cmap=cmap)
    cont_list = [sigma * i for i in range(-6, 24, 3) if i != 0]
    ax.contour(channel,
               cont_list,
               colors='black')

    # zoom-in on pixels
    ax.set_xlim(110, 150)
    ax.set_ylim(107, 147)

    if HST:   # mark lensing knots
        import astrolib.coords as coords
        # to convert sexidesimal coord --> deg in tuple
        p_list = [
            coords.Position("11:31:51.5769 -12:31:58.9708").j2000(),  # A
            coords.Position("11:31:51.5678 -12:31:57.8143").j2000(),  # B
            coords.Position("11:31:51.5366 -12:32:00.1192").j2000(),  # C
            coords.Position("11:31:51.3658 -12:31:58.1183").j2000()    # D
        ]
        for p1 in p_list:
            ra, dec = p1
            ra_px, dec_px = wcs.wcs_sky2pix(ra, dec, 1)
            l1, = ax.plot([ra_px], [dec_px],
                          color='white',
                          marker='+',
                          mec="black", mew=1.25, ms=6,
                          zorder=3.1)    # lower zorder are drawn first

        # fg Galaxy
        fg_list = [coords.Position("11:31:51.44295 -12:31:58.30659").j2000()]
        for p in fg_list:
            ra, dec = p
            ra_px, dec_px = wcs.wcs_sky2pix(ra, dec, 1)
            l2, = ax.plot([ra_px], [dec_px],
                          color='white',
                          marker="*", mec="black", mew=1.5, ms=6,
                          zorder=3.1)    # lower zorder are drawn first

    images.append(im)


# if HST:
#    # put up legend on last panel --> HST markers
#     g[0].legend([l1, l2],
#                 ["QSO lensing knots", "Lensing Galaxy"],
#                 loc=3,
#                 numpoints=1,
#                 handlelength=1,
#                 frameon=True,
#                 framealpha=0.8,
#                 prop=dict(size=8),
#                 fancybox=True)


# label with velocities
use_path_effect = False         # Fancy text
try:
    from matplotlib.patheffects import withStroke
except ImportError:
    use_path_effect = False

for i, ax in enumerate(g):
    channel_number = start_channel + i
    v = vel.to_vel(channel_number) / 1.e3
    t = ax.add_inner_title(r"$v=%4.1f\ {\rm km s}^{-1}$" % (v),
                           # (u'$\\mathrm{km/s}$')
                           loc=2,
                           frameon=False
                           )
    t.patch.set_alpha(1.0)
    if use_path_effect:
        t.txt._text.set_path_effects([withStroke(foreground="w",
                                                 linewidth=2.5)])
    else:
        # show velocity text in some color
        t.txt._text.set_color("#243106")
        t.txt._text.set_fontweight('black')

# Beam
g[-1].add_beam_size(minor_px, major_px, BPA_deg,      # y, x, angle
                    loc=4, frameon=True,
                    patch_props=dict(facecolor='gray', alpha=0.8, edgecolor='k', fill=True))
# g[-1].add_compass(loc=3)

# make colorbar
cb = plt.colorbar(im, cax=cax)
cb.set_label("Flux Density [mJy beam" + r"$^{-1}$]")
cb.set_ticks([-0.01, 0, 0.01, 0.02, 0.03])   # place tick locations

# customize to mJy B^-1 instead of Jy B^-1, uncomment to use same as cb.set_ticks
cb.set_ticklabels(['-10', '0', '10', '20', '30'])

# adjust norm
# norm.vmin = -0.007      # min in map
norm.vmin = -0.015        # show more contrast
norm.vmax = 0.03          # max in map
for im in images:
    im.changed()

# fig.suptitle("CO(2-1) Channel Maps", fontsize=22, x=0.5, y=0.92)
fig1 = plt.gcf()
plt.show()
fits_cube.close()

User_input = raw_input('Save figure? (Y/N): ')
if User_input == 'Y':
    filename = "co_channel_maps.eps"
    fig1.savefig(Plotpath + filename, dpi=100,
                bbox_inches="tight", pad_inches=0.1)
    print "-- Saved figure as : %s --" % (Plotpath + filename)
