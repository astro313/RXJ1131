#!/usr/bin/env python
'''

plot source locations from lens model of various channel as markers on observed 1st moment map.
Kinematics.

Last Modified: 28 May 16

History:
28 May 2016
  - add beam to 1st moment + markers plot
  - add scale bar
  - add color bar
  - use patheffects to make markers stand out from background
17 May 16:
  - move var kpc_to_m to calcmass.py
15 May 16:
  - adjust residual fig outlook
13 May 16:
  - moved Mdyn calc to another script
  - limit r_t to 15kpc
12 May 16:
  - replace r_s with r_t since it is the transition radius and not necessarily the scale length
11 May 16:
  - use curve fit, problem is non-linear in parameters
  - note, we are also interested in reporting the correlation between the degenerate parameters --> need pearson R coefficients
  - removed least sq fitting to arctan, becuase non-linear in param.
  - arctangent form v is asymptotic v, careful when interpreting results
10 May 16:
  - remove ODR fit to rot. curve, since xerr drives poor fit
  - fix R function, major axis is not affected by inclination angle in circular
  - update error bars on velocity to using velocity range of models
  - add error on central position to allow fit for position error on central chan
  - remove code to fit model where R dep. on i
  - remove deriving Mdyn using old models
  - fix bug RA_major_err, Dec_major_err
  - add xerr bar on v_0
  - add ODR to fit arctangent model with xerr
  - fix bug in PA calc.
09 May 16:
  - fix a bug in run_odr_for_model()
  - added func to calc offset taking into acct inclination effect on y
  - added code to fit inclination altogether, updating R on-the-fly, due to dependence on inclination on Y-axis
  - update code to calc. distance using arcsec, not pixel !!
  - plot Mdyn as a func of R, of all models, all bug fixed, basically showing identical trend
  - fix bug in calc dist for model where R dep. on i, [3] instead of [4]
  - fix bug in calc error on physical dist.
  - removed ODR for fitting to velocity models, becuase isn't working as expected.
  - added fit to tangent model, which is highly sensitive to the characteristic radius...
08 May 16:
  - plot PV, and rotation curve along major axis
  - fit to rotation curve, just solid body
  - calc. enclosed dyn mass using velocity from fit
07 May 16:
  - fit through data to get `major axis` using fmin, leastsq, ODR
05 May 16:
  - added error bars to positions
04 May 16:
  - created code

Note
----
Can only save 1st plot (1st moment + source marker) as pdf, not eps

Fitting data w/ unc. on both x, y using ordinary least squares can lead to bias in the solution. Use Orthogonal Distance Regression (ODR), rather than minimizing the sum of squared errors in the dependent variable, ODR minimizes the orthogonal distance from the data to the fitted curve by adjusting both the model coefficients and an adjustment to the values of the independent variable.

 The results of the covariance matrix, as implemented by optimize.curvefit and optimize.leastsq rely on assumptions regarding the probability distribution of the errors and the interactions between parameters; interactions which may exist, depending on the fit function f(x).

The best way to deal with a complicated f(x) may be to use bootstrap.

'''
import matplotlib
import matplotlib.pyplot as plt
import pyfits
import pywcs
import wcsaxes
import numpy as np
from astropy.cosmology import WMAP9


font = {'family': 'Arial Narrow',
        'weight': 'bold',
        'size': 18}
matplotlib.rc('font', **font)
matplotlib.rc('text', usetex=True)

Plotpath = '/Users/admin/Research/RXJ1131/Figures/'
firstMom_obs = '/Users/admin/Research/RXJ1131/PdBI/data/30Apr16/centralizedCube4GILDAS-momCLIP5sigma.velo.fits'

interval = 50.        # km/s between velocity contour
DecCentroid = -12.5328629
RACentroid = 172.96434563
CellSize = 0.5        # arcsec per pixel in 1D
redshift = 0.6357
b_arcsec = 1.8       # from C06
a_arcsec = 3.25
i_rad = np.arccos(b_arcsec/a_arcsec)
i_deg = i_rad * 180./np.pi

# source positions from models of different channel, offset in arcsec
p_list = [(-0.20, 0.12), (-0.02, 0.11), (0.20, -0.13), (0.29, -0.01),
          (0.57, -0.76), (0.87, -0.57), (1.04, -0.56), (1.32, -0.77)]
# corresponding error in arcsec
p_list_err = [(0.39, 0.3), (0.25, 0.2), (0.06, 0.21), (0.14, 0.22),
              (0.12, 0.29), (0.09, 0.15), (0.06, 0.12), (0.34, 0.45)]
# corresponding velocity in km/s for each points in p_list
z = [344.46, 344.46, 236.82, 129.16, 21.54, -86.08, -193.76, -301.38]

plotMajorFit = False
plotMajor_FirstMom = False
plotPV = False
plotRot = False
plotODR_arctan = False
plotMdyn = False
# -----------------------------------------------
#  Read fits file and Set up WCS
# -----------------------------------------------
a = pyfits.open(firstMom_obs)
hdr = a[0].header
cdelt1 = np.abs(hdr['CDELT1'] * 3600)       # arcsec/pix
cdelt2 = np.abs(hdr['CDELT2'] * 3600)
major_px = hdr['BMAJ'] * 3600 / cdelt1                 # arcsrc/pix
minor_px = hdr['BMIN'] * 3600 / cdelt2
BPA_deg = hdr['BPA']         # should be 13 deg
wcs = pywcs.WCS(hdr).sub([1, 2])
plotwcs = wcsaxes.WCS(hdr, naxis=2)
LensRA_PX, LensDEC_PX = wcs.wcs_sky2pix(RACentroid, DecCentroid, 1)
dat = a[0].data.reshape(a[0].data.shape[2], a[0].data.shape[3])

# -----------------------------------------------
#  Set up Figure
# -----------------------------------------------
fig = plt.figure(1)
ax = fig.add_subplot(1, 1, 1, projection=plotwcs)
ra, dec = ax.coords
ra.set_major_formatter('hh:mm:ss.s')
ra.display_minor_ticks(True)
dec.display_minor_ticks(True)
ra.set_minor_frequency(5)
ra.set_axislabel(" R.A. (J2000)", minpad=0.5)
dec.set_axislabel(" Dec", minpad=-0.4)
# ra.set_separator((':', ':'))
ra.set_ticks(size=10, width=1.5)
dec.set_ticks(size=10, width=1.5)
ax.set_xlim(120, 143)
ax.set_ylim(110, 140)


# beam
def draw_ellipse(ax):
    """ Draw beam in data coordinate
    """
    from mpl_toolkits.axes_grid1.anchored_artists import AnchoredEllipse
    ae = AnchoredEllipse(ax.transData, width=minor_px, height=major_px,
                         angle=BPA_deg,
                         loc=3, pad=0.5, borderpad=0.4, frameon=False)
    ax.set_alpha(0.6)
    ax.add_artist(ae)


def draw_sizebar(ax, z, pxScale, arcsec=1.0):
    """ Draw a scale bar with length of 0.1 in Data coordinate
    (ax.transData) with a label underneath.

    z: float
        redshift
    pxScale: float
        arcsec per 1D-pixel
    arcsec: float
        how many arcsec to draw
    """
    from mpl_toolkits.axes_grid1.anchored_artists import AnchoredSizeBar

    kpc = WMAP9.kpc_proper_per_arcmin(z).value / 60. * arcsec

    # num px
    length = arcsec / pxScale
    asb = AnchoredSizeBar(ax.transData,
                          length,
                          str(round(kpc))+" kpc at " + r"$z = $" + (("{0:.3f}").format(z)),
                          loc=8,
                          pad=0.1, borderpad=0.5, sep=5,
                          frameon=False)
    ax.add_artist(asb)

# setup color
cmap = plt.cm.jet
import matplotlib.colors as mcolors
norm = mcolors.Normalize()

# -----------------------------------------------
# 1st moment map
# -----------------------------------------------
im = ax.imshow(dat, origin="lower", norm=norm, cmap=cmap, alpha=0.65)
cont_list = [interval * i for i in range(-10, 10) if i != 0]
ax.contour(dat,
           cont_list,
           colors='black', alpha=0.5)
ax.plot([LensRA_PX], [LensDEC_PX], marker='x', zorder=5,    # lens galaxy
         mec='black', ms=8, mew=1.75)

from matplotlib.pyplot import cm
color = iter(cm.rainbow_r(np.linspace(0, 1, len(p_list))))
import matplotlib.patheffects as path_effects
ra_px_ls = []
dec_px_ls = []
ra_coord_ls = []
dec_coord_ls = []
ra_err_ls = []
dec_err_ls = []
for i, p1 in enumerate(p_list):
    c = next(color)
    ra_offset_arcsec, dec_offset_arcsec = p1

    # convert into coord in deg
    ra_coord = ra_offset_arcsec / 3600. + RACentroid
    dec_coord = dec_offset_arcsec / 3600. + DecCentroid
    ra_coord_ls.append(ra_coord)
    dec_coord_ls.append(dec_coord)

    # convert to pixel on image
    ra_px, dec_px = wcs.wcs_sky2pix(ra_coord, dec_coord, 1)
    ra_px_ls.append(ra_px[0])
    dec_px_ls.append(dec_px[0])

    # convert arcsec to pixel offset on image
    ra_err = p_list_err[i][0] / CellSize
    dec_err = p_list_err[i][1] / CellSize
    ra_err_ls.append(ra_err)
    dec_err_ls.append(dec_err)
    m = ax.errorbar(x=[ra_px], y=[dec_px], xerr=[ra_err], yerr=[dec_err],
                marker='+', ms=8, mec=c, mew=2, zorder=3.1, ecolor=c,
                capsize=1.5, capthick=1.5, lw=2, label=str(i))

    # to make marker stand out from background
    m[0].set_path_effects([path_effects.Stroke(linewidth=2.5, foreground='black'), path_effects.Normal()])
    m[1][0].set_path_effects([path_effects.Stroke(linewidth=2.5, foreground='black'), path_effects.Normal()])
    m[1][1].set_path_effects([path_effects.Stroke(linewidth=2.5, foreground='black'), path_effects.Normal()])
    m[1][2].set_path_effects([path_effects.Stroke(linewidth=2.5, foreground='black'), path_effects.Normal()])
    m[1][3].set_path_effects([path_effects.Stroke(linewidth=2.5, foreground='black'), path_effects.Normal()])
    m[2][0].set_path_effects([path_effects.Stroke(linewidth=2.5, foreground='black'), path_effects.Normal()])
    m[2][1].set_path_effects([path_effects.Stroke(linewidth=2.5, foreground='black'), path_effects.Normal()])
#     print ra_px, dec_px

# ------ additional setup -----
# matplotlib.rcParams['legend.handlelength'] = 0
# matplotlib.rcParams['legend.markerscale'] = 0
# ax.legend(loc=2, bbox_to_anchor=(0.9, 0.9), borderaxespad=0., fontsize=15, numpoints=1)
draw_ellipse(ax)
draw_sizebar(ax, redshift, cdelt1)

# make colorbar
cb = fig.colorbar(im, fraction=0.1, pad=0.0)
cb.set_label("Velocity [km s" + r"$^{-1}$]")

plt.tight_layout(h_pad=0.4, w_pad=0.1)
plt.show()
a.close()
User_input = raw_input('Save figure? (Y/N): ')
if User_input == 'Y':
    filename = "veloGradient_markers.pdf"
    fig.savefig(Plotpath + filename, dpi=100,
                bbox_inches="tight", pad_inches=0.1)
    print "-- Saved figure as : %s --" % (Plotpath + filename)

# -------------------------------
theta_rad = abs(np.arctan((dec_coord_ls[0] - dec_coord_ls[-1]) / ((ra_coord_ls[0] - ra_coord_ls[-1]) * np.cos(np.mean(dec_coord_ls) * np.pi/180.))))
theta_deg = theta_rad * 180. / np.pi
PA_deg = theta_deg + 90.
print(" \n PA along the `major axis`: {} deg. \n ").format(PA_deg)

# --------------------------------------------
# fit a major axis through the px coordinate
# --------------------------------------------
xdata = np.array(ra_px_ls)
ydata = np.array(dec_px_ls)
xerr = np.array(ra_err_ls)
yerr = np.array(dec_err_ls)
import scipy.optimize as optimize
#
# using yerr only
#
fun = lambda x, m, b: m*x + b
def residuals(theta, x, y, err):
    '''
    minimize the residual, which, for fmin() is a scalar

    Parameters
    ----------
    theta: list
        [m, b]

    Returns
    -------
    chisq: float
        scalar
    '''
    model = fun(x, *theta)
    chisq = np.sum(((y - model) / err)**2)
    return chisq

x0 = [10., 90.]
best = optimize.fmin(residuals, x0, args=(
    xdata, ydata, yerr))

if plotMajorFit:
    plt.errorbar(ra_px_ls, dec_px_ls, xerr=ra_err_ls, yerr=dec_err_ls,
             marker='+', ms=8, mew=1.25, capsize=5, capthick=2, linestyle='None')
    plt.plot(xdata, fun(xdata, best[0], best[1]))
    plt.tight_layout()
    plt.show()

#
# using leastsq, just to practice writing fitting codes
#
def fun2(p, x):
    " Use instead a unpacked parameter p for input"
    m, b = p
    return m * x + b
# for least-square minimization, errorfunc should return an array not scalar
errfunc = lambda p, x, y, err: (y - fun2(p, x)) / err

pfit, pcov, infodict, errmsg, success = optimize.leastsq(errfunc, x0, args=(xdata, ydata, yerr), full_output=1)

if (len(ydata) > len(x0)) and pcov is not None:
    # reduced-chi
    s_sq = (errfunc(pfit, xdata, ydata, yerr)**2).sum() / (len(ydata) - len(x0))
    # covar = fractional covar * reduced-chi
    pcov = pcov * s_sq
else:
    pcov = inf
error = []
for i in range(len(pfit)):
    try:
        error.append(np.absolute(pcov[i][i])**0.5)
    except:
        error.append(0.00)
perr = np.array(error)


# ------------------------------------------
# take into account err in both x, y directions
from scipy.odr import *
_model = Model(fun2)
data = RealData(xdata, ydata,
                sx=xerr, sy=yerr)
odr = ODR(data, _model, beta0=x0)
out = odr.run()
param = out.beta
cov = out.cov_beta
# print("\nStandard Covariance Matrix : \n", cov, "\n")
uncertainty = out.sd_beta
res_variance = out.res_var
# if res_variance < 1.0 :
#     uncertainty = uncertainty/np.sqrt(res_variance)
# if False: # debugging print statements
#     print("sd_beta",out.sd_beta)
#     print("cov_beta",out.cov_beta)
#     print("delta",out.delta)
#     print("eps",out.eps)
#     print("res_var",out.res_var)s
#     print("rel_error",out.rel_error)
# out.pprint()

# Calculate initial residuals and the 'adjusted error' for each data point
delta = out.delta  # estimated x-component of the residuals
eps   = out.eps    # estimated y-component of the residuals
# (xstar,ystar) is the point where the 'residual line' (in black)
#   intersects the 'ellipse' created by xerr & yerr.
xstar = (xerr*np.sqrt(((yerr*delta)**2) / ((yerr*delta)**2 + (xerr*eps)**2)))
ystar = (yerr*np.sqrt(((xerr*eps)**2) / ((yerr*delta)**2 + (xerr*eps)**2)))
adjusted_err = np.sqrt(xstar**2 + ystar**2)
residual = np.sign(ydata-fun2(param, xdata))*np.sqrt(delta**2 + eps**2)

## Plot
if plotMajorFit:
    fig = plt.figure(facecolor="0.98")   # with light gray background
    fig.subplots_adjust(hspace=0)
    fit = fig.add_axes((.1, .3, .8, .6))
    fit.set_xticklabels(())
    plt.ylabel("Dec")

    fit.plot(xdata, ydata, 'ko', xdata, fun2(param, xdata))
    fit.errorbar(xdata, ydata, xerr=xerr, yerr=yerr, fmt='k+')
    fit.set_yscale('linear')
    a = np.array([out.xplus, xdata])   # out.xplus = x + delta
    b = np.array([out.y, ydata])
    # to show where the fit is relative to the points
    fit.plot(np.array([a[0][0], a[1][0]]), np.array([b[0][0], b[1][0]]),
             'k-', label='Residuals')
    for i in range(1, len(ydata)):
        fit.plot(np.array([a[0][i], a[1][i]]), np.array([b[0][i], b[1][i]]),
             'k-')
    fit.legend(loc='upper left')
    fit.grid()
    plt.minorticks_on()
    # separate plot to show residuals
    residuals = fig.add_axes((.1, .1, .8, .2))
    residuals.errorbar(x=xdata, y=residual, yerr=adjusted_err,
                                fmt="k+", label="Residuals")
    # make sure residual plot has same x axis as fit plot
    residuals.set_xlim(fit.get_xlim())
    plt.axhline(y=0, color='b')
    plt.xlabel("RA ")
    plt.ylabel("Residuals")
    residuals.grid()
    plt.minorticks_on()
    plt.show()

youtfit = fun2(out.beta, xdata)

# -------------------------------------------
# Plot velocity field again, with major axis
# -------------------------------------------
if plotMajor_FirstMom:
    fig = plt.figure(figsize=(8, 12), dpi=100)
    plt.subplots_adjust(top=0.9)
    ax = fig.add_subplot(1, 1, 1, projection=plotwcs)
    ra, dec = ax.coords
    ra.set_major_formatter('hh:mm:ss.s')
    ra.display_minor_ticks(True)
    dec.display_minor_ticks(True)
    ra.set_minor_frequency(5)
    ra.set_axislabel(" RA (J2000)", minpad=0.5)
    dec.set_axislabel(" Dec (degrees)", minpad=-0.4)
    ra.set_separator((':', ':'))
    ra.set_ticks(size=10, width=1.5)
    dec.set_ticks(size=10, width=1.5)
    ax.set_xlim(120, 145)
    ax.set_ylim(110, 140)

    cmap = plt.cm.jet
    import matplotlib.colors as mcolors
    norm = mcolors.Normalize()

    # -----------------------------------------------
    # 1st moment map
    # -----------------------------------------------
    im = ax.imshow(dat, origin="lower", norm=norm, cmap=cmap)
    cont_list = [interval * i for i in range(-10, 10) if i != 0]
    ax.contour(dat,
               cont_list,
               colors='black')
    ax.plot([LensRA_PX], [LensDEC_PX], marker='x', zorder=1, mec='k', ms=8)

    ax.errorbar(x=ra_px_ls, y=dec_px_ls, xerr=ra_err_ls, yerr=dec_err_ls,
                marker='+', ms=8, linestyle='None', mew=1.25, zorder=3.1, capsize=5 , capthick=2, lw=2)
    ax.plot(xdata, youtfit, 'k--', lw=1.2)
    plt.tight_layout()
    plt.show()


# -------------------------------------------------------
# PV and rot. curve exclude 2nd component in redmost channel
# -------------------------------------------------------
# convert the px positions along the fitted major axis to RA, Dec to calc. separation
deg_to_arcsec = 3600.
RA_major, Dec_major = wcs.wcs_pix2sky(ra_px_ls, youtfit, 1)       # degree

# arcsec
RA_major_err = np.array([p_list_err[i][0] for i in range(len(p_list_err))])
Dec_major_err = np.array([p_list_err[i][1] for i in range(len(p_list_err))])

def fitted_err_on_y(deltaP, x):
    """ calc y_error of point on fit

    """
    deltam, deltab = deltaP
    dy = np.sqrt((x * deltam)**2 + (deltab)**2)
    return dy

def calc_dist(xloc, yloc, delta_x, delta_y): # , xcenter_err, ycenter_err):
    """ calc distance in arcsec offset along major axis

    Parameters
    ----------
    xloc: array
        RA in deg of source position
    yloc: array
        dec in deg of source position along the fitted major axis
    delta_x: array
        RA in arcsec
    delta_y: array
        Dec in arcsec
    Return
    ------
    offset: array
        in arcsec along the major axis, w.r.t the central emission position

    """
    delta_x0 = delta_x[4]
    delta_y0 = delta_y[4]

    diff_x = np.array([xloc[i] - xloc[4] for i in range(len(xloc))])
    diff_y = np.array([yloc[i] - yloc[4] for i in range(len(yloc))])
    offset_sq = (diff_x * np.cos(np.mean(yloc/3600. * np.pi/180.)))**2 + diff_y**2
    R = np.sqrt(offset_sq)
    offset_err = np.sqrt((diff_x*delta_x)**2 + (diff_y*delta_y)**2 + (diff_x*delta_x0)**2 + (diff_y*delta_y0)**2) / R

    # put error bar on center point as well, in arcsec
    xcenter_err, ycenter_err = p_list_err[4]
    offset_err[np.where(np.isnan(offset_err))[0][0]] = np.sqrt(delta_x0**2 +delta_y0**2)

    return R, offset_err

# youfitErr = fitted_err_on_y(uncertainty, xdata)     # very big..., because of x
# use uncertainty on ydata instead
# reasonble because people slice along their map, they don't fit for the major axis

offset, offset_err = calc_dist(RA_major*deg_to_arcsec,
                               Dec_major*deg_to_arcsec,
                               RA_major_err, Dec_major_err)
# to show blue, red on left and right of central
offset[5:] = -offset[5:]
off = list(offset)
offset_err = list(offset_err)
off.pop(1)
z.pop(1)
offset_err.pop(1)
z_err = 21.5*5/2.     # range of channels combined to make models/2. --> error bar
if plotPV:
    f, ax = plt.subplots()
    ax.errorbar(off, z, yerr=z_err, fmt='ko', xerr=offset_err)
    plt.ylabel('v_r = v sin i')
    plt.xlabel(' offset from line center position ["] ')
    plt.title('PV along major axis at PA {:.2f} deg'.format(PA_deg))
    plt.tight_layout()
    plt.minorticks_on()
    plt.show(block=False)
    User_input = raw_input('Save figure? (Y/N): ')
    if User_input == 'Y':
        filename = "PV_major.eps"
        f.savefig(Plotpath + filename, dpi=100,
                    bbox_inches="tight", pad_inches=0.1)
        print "-- Saved figure as : %s --" % (Plotpath + filename)

offset, offset_err = calc_dist(RA_major*deg_to_arcsec, Dec_major*deg_to_arcsec, RA_major_err, Dec_major_err)
off = list(offset)
offset_err = list(offset_err)
off.pop(1)
offset_err.pop(1)
if plotRot:
    f, ax = plt.subplots()
    ax.errorbar(off[3], z[3], xerr=offset_err[3], yerr=z_err, fmt='go', label='line center')
    ax.errorbar(off[:3], z[:3], yerr=z_err, fmt='ro', label='red', xerr=offset_err[:3])

    plt.ylabel('|v_r = v sin i|')
    ax.errorbar(off[4:], abs(np.array(z[4:])), yerr=z_err, fmt='bo', label='blue', xerr=offset_err[4:])
    plt.xlabel(' radial offset from line center position ["] ')
    plt.title('"Rotation Curve" taken along major axis at PA {:.2f} deg'.format(PA_deg))
    ax.set_ylim(-50, 450)
    ax.set_xlim(-0.2, 1.2)
    plt.legend()
    plt.tight_layout()
    plt.minorticks_on()
    plt.show(block=False)
    User_input = raw_input('Save figure? (Y/N): ')
    if User_input == 'Y':
        filename = "RotCurve_major.eps"
        f.savefig(Plotpath + filename, dpi=100,
                    bbox_inches="tight", pad_inches=0.1)
        print "-- Saved figure as : %s --" % (Plotpath + filename)

# --------------------------------------------------------------
#  disk model
# ---------------------------------------------------------------
def velo_disk(V_obs, theta=0., i=30.):
    """ simple axisymmetric thin disk model

    Parameters
    ----------
    V_obs: float
        km/s
    theta: float
        radian, 0 along the major axis
    i: float
        in deg. Fix to 30. deg if we don't have better constraints

    """

    V_sys = 0.       # at the rest-frame of high-z gal
    # theoretically, a thin disk:
    # V_obs = V_sys + np.sin(i) * (V_tang * np.cos(theta) + V_radial * np.sin(theta))
    # in practice:
    V_obs = V_sys + np.sin(i) * (V_tang * np.cos(theta))
    return V_obs


def offset_to_physicalR(R_arcsec, z):
    """ Convert offset to physical radius"""
    arcsec_to_kpc = WMAP9.kpc_proper_per_arcmin(z).value / 60.
    return arcsec_to_kpc * R_arcsec


def arctang(R_kpc, vcsini, r_t, v_0):
    # hacky step to force r_t to be physical;
    if r_t > 0.1 and r_t < 15:
        return vcsini*(2./np.pi)*np.arctan(R_kpc/r_t) + v_0
    else:  # bascially reject this fit
        return 1e-5


def arctang2(p, R_kpc):
    vcsini, r_t, v_0 = p
    # hacky step to force r_t to be physical
    if r_t > 0.1 and r_t < 15:
        return vcsini*(2./np.pi)*np.arctan(R_kpc/r_t) + v_0
    else:  # bascially reject this fit
        return 1e-5


ydata = abs(np.array(z))
yerr = z_err
xdata = offset_to_physicalR(np.array(off), redshift)
xerr = offset_to_physicalR(np.array(offset_err), redshift)

# plot velocity v.s. physical radius
if plotRot:
    plt.errorbar(xdata, ydata, yerr=yerr, fmt='ko', xerr=xerr)
    plt.ylabel('V sin i')
    plt.xlabel('Physical distance from line center [kpc]')
    plt.xlim(-1, 8.5)
    plt.ylim(-50, 450)
    plt.tight_layout()
    plt.minorticks_on()
    plt.show()


# fit tangent model to PV along major axis
xdataPV = list(xdata[:4])
_blue = -xdata[4:]
for i in range(len(_blue)): xdataPV.append(_blue[i])
xdataPV = np.array(xdataPV)

# ODR
x4 = [max(z), max(xdataPV), 20]
_model = Model(arctang2)
data = RealData(xdataPV, np.array(z),
                sx=xerr, sy=yerr)
odr = ODR(data, _model, beta0=x4)
out = odr.run()
param = out.beta
covar = out.cov_beta
res_variance = out.res_var
if res_variance < 1.0 :
    uncertainty = out.sd_beta/np.sqrt(res_variance)
else:
    uncertainty = out.sd_beta
pearsonR = covar/np.outer(uncertainty, uncertainty)
DoF = len(xdataPV)-len(out.beta)
# Print ODR Status and Results
print 'Return Reason:\n', out.stopreason, '\n'
print 'Estimated Parameters:\n', param, '\n'
print 'Parameter Standard Errors:\n', out.sd_beta, '\n'
# print 'Covariance Matrix:\n', covar, '\n'
print 'Degree of Freedom:\n', DoF, '\n'
print 'Pearson R coefficients:\n', pearsonR, '\n'

# Chi^2
chi2 = 0.
for i in range(len(np.array(z))):
    residual = np.array(z)[i] - arctang2(param, xdataPV[i])
    sigma = (xerr[i]**2. + yerr**2.)**0.5
    chi2 += (residual / sigma)**2.
# Reduced Chi^2
ndof = len(z) - len(param)
redchi2 = chi2 / ndof
print 'Degrees of Freedom:\t', ndof
print 'Chi-Square:\t\t', chi2
print 'Reduced Chi-Square:\t', redchi2

if plotODR_arctan:
    # Calculate initial residuals and the 'adjusted error' for each data point
    delta = out.delta  # estimated x-component of the residuals
    eps   = out.eps    # estimated y-component of the residuals
    # (xstar,ystar) is the point where the 'residual line' (in black)
    #   intersects the 'ellipse' created by xerr & yerr.
    xstar = (xerr*np.sqrt(((yerr*delta)**2) / ((yerr*delta)**2 + (xerr*eps)**2)))
    ystar = (yerr*np.sqrt(((xerr*eps)**2) / ((yerr*delta)**2 + (xerr*eps)**2)))
    adjusted_err = np.sqrt(xstar**2 + ystar**2)
    residual = np.sign(np.array(z)-arctang2(param, xdataPV))*np.sqrt(delta**2 + eps**2)
    ## Plot
    fig = plt.figure(facecolor="0.98")   # with light gray background
    fig.subplots_adjust(hspace=0)
    fit = fig.add_subplot(211)       #, adjustable='box', aspect=1.2)
    fit.set_xticklabels(())
    plt.ylabel("V obs")

    space = 0.2
    xspan = xdataPV.max() - xdataPV.min()
    yspan = np.array(z).max() - np.array(z).min()
    xspace = space * xspan
    yspace = space * yspan
    xarray = np.linspace(xdataPV.min()-xspace, xdataPV.max()+xspace, 500)
    fit.plot(xdataPV, np.array(z), 'ko', xarray, arctang2(param, xarray))
    fit.errorbar(xdataPV, np.array(z), xerr=xerr, yerr=yerr, fmt='k+')
    plt.errorbar(xdataPV, np.array(z), fmt='ko', xerr=xerr, yerr=yerr)
    fit.set_yscale('linear')
    a = np.array([out.xplus, xdataPV])
    b = np.array([out.y, np.array(z)])
    # to show where the fit is relative to the points
    fit.plot(np.array([a[0][0], a[1][0]]), np.array([b[0][0], b[1][0]]),
             'k-', label='Residuals')
    for i in range(1, len(np.array(z))):
        fit.plot(np.array([a[0][i], a[1][i]]), np.array([b[0][i], b[1][i]]),
             'k-')
    fit.legend(loc='upper left')
    fit.grid()
    residuals = fig.add_subplot(212)
    residuals.errorbar(x=xdataPV, y=residual, yerr=adjusted_err,
                            fmt="k+", label="Residuals")
    residuals.set_xlim(fit.get_xlim())
    plt.axhline(y=0, color='b')
    plt.xlabel("R")
    plt.ticklabel_format(style='plain', useOffset=False, axis='x')
    plt.ylabel("Residuals")
    residuals.grid()
    plt.minorticks_on()
    plt.tight_layout()
    plt.show()

# curve fit
pfit, pcov = optimize.curve_fit(arctang, xdataPV, np.array(z), p0=x4, sigma=yerr, absolute_sigma=True)
perr = np.sqrt(np.diag(pcov))
print("\nFit parameters and parameter errors from curve_fit method:")
print("pfit = ", pfit)
print("perr = ", perr)
chi2 = 0.
for i in range(len(np.array(z))):
    residual = np.array(z)[i] - arctang(xdataPV[i], *pfit)
    sigma = (xerr[i]**2. + yerr**2.)**0.5
    chi2 += (residual / sigma)**2.
# Reduced Chi^2
ndof = len(z) - len(param)
redchi2 = chi2 / ndof
print 'Degrees of Freedom:\t', ndof
print 'Chi-Square:\t\t', chi2
print 'Reduced Chi-Square:\t', redchi2
print("\n Best-fit asymptotic V: {} km/s \n").format(pfit[0])
v_rot = pfit[0]/np.sin(i_rad)
print("V_rot = V_obs/sin i: {} km/s \n").format(v_rot)

## Plot
fig = plt.figure(facecolor="0.98")   # with light gray background
fig.subplots_adjust(hspace=0)
fit = fig.add_subplot(111)
plt.ylabel("Velocity [km/s]")
plt.xlabel("Physical Radius from dynamical center [kpc]")
space = 0.2
xspan = xdataPV.max() - xdataPV.min()
yspan = np.array(z).max() - np.array(z).min()
xspace = space * xspan
yspace = space * yspan
xarray = np.linspace(xdataPV.min()-xspace, xdataPV.max()+xspace, 500)
fit.plot(xdataPV, np.array(z), 'ko', xarray, arctang(xarray, *pfit), 'b--')
fit.errorbar(xdataPV, np.array(z), xerr=xerr, yerr=yerr, fmt='k+')
plt.minorticks_on()
plt.tight_layout()
plt.show(block=False)
User_input = raw_input('Save figure? (Y/N): ')
if User_input == 'Y':
    filename = "bestfit_PV.eps"
    fig.savefig(Plotpath + filename, dpi=100,
                bbox_inches="tight", pad_inches=0.1)
    print "-- Saved figure as : %s --" % (Plotpath + filename)

# monte-carlo: generates random data points starting from the given data plus a random variation based on the systematic error, take into acct systematic uncertainties
nsam = 500
ps = []
debug = False

for i in range(nsam):
    randomDelta = np.random.normal(0., yerr, size=len(z))
    randomdataY = np.array(z) + randomDelta

    if debug:
        # verify it's drawing from distribution we expect
        randomDelta = np.random.normal(0., yerr, size=1000)
        count, bins, ignored = plt.hist(randomDelta, 30, normed=True)
        plt.plot(bins, 1/(yerr * np.sqrt(2 * np.pi)) *
                 np.exp(-(bins)**2 / (2 * yerr**2)), linewidth=2, color='r')
        plt.show()

    # based on the systematic errors on each physical separation
    randomDeltaX = np.array([np.random.normal(0., derr, 1)[0]
                             for derr in xerr])
    randomdataX = xdataPV + randomDeltaX
    randomfit, randomcov = optimize.curve_fit(arctang,
                                              randomdataX, randomdataY, p0=x4)
    ps.append(randomfit)

ps = np.array(ps)
mean_pfit = np.median(ps, 0)
Nsigma = 1.       # 68.3% CI
err_pfit_MC = Nsigma * np.std(ps, 0)
print "\n MC method :"
print "pfit = ", mean_pfit
print "perr = ", err_pfit_MC


