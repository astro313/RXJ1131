#!/usr/bin/env python
'''

plot source locations from lens model of various channel as markers on observed 1st moment map

Last Modified: 10 May 16

TODO:
get source size in old papers
get i from optical
add xerr bar on v_0, incl in fit

History:
10 May 16:
  - remove ODR fit to rot. curve, since xerr drives poor fit
  - fix R function, major axis is not affected by inclination angle in circular
  - update error bars on velocity to using velocity range of models
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

'''
import matplotlib
import matplotlib.pyplot as plt
import pyfits
import pywcs
import wcsaxes
import numpy as np

font = {'family': 'Arial Narrow',
        'weight': 'bold',
        'size': 18}
matplotlib.rc('font', **font)


Plotpath = '/Users/admin/Research/RXJ1131/Figures/'
firstMom_obs = '/Users/admin/Research/RXJ1131/PdBI/data/30Apr16/centralizedCube4GILDAS-momCLIP5sigma.velo.fits'


interval = 50.        # km/s between velocity contour
DecCentroid = -12.5328629
RACentroid = 172.96434563
CellSize = 0.5        # arcsec per pixel in 1D
redshift = 0.6357

# source positions from models of different channel, offset in arcsec
p_list = [(-0.20, 0.12), (-0.02, 0.11), (0.20, -0.13), (0.29, -0.01),
          (0.57, -0.76), (0.87, -0.57), (1.04, -0.56), (1.32, -0.77)]
# corresponding error in arcsec
p_list_err = [(0.39, 0.3), (0.25, 0.2), (0.06, 0.21), (0.14, 0.22),
              (0.12, 0.29), (0.09, 0.15), (0.06, 0.12), (0.34, 0.45)]
# corresponding velocity in km/s for each points in p_list
z = [344.46, 344.46, 236.82, 129.16, 21.54, -86.08, -193.76, -301.38]

plotMajorFit = True
plotMajor_FirstMom = True
plotPV = True
plotRot = True
plotMdyn = False
# -----------------------------------------------
#  Read fits file and Set up WCS
# -----------------------------------------------
a = pyfits.open(firstMom_obs)
hdr = a[0].header
wcs = pywcs.WCS(hdr).sub([1, 2])
plotwcs = wcsaxes.WCS(hdr, naxis=2)
LensRA_PX, LensDEC_PX = wcs.wcs_sky2pix(RACentroid, DecCentroid, 1)
dat = a[0].data.reshape(a[0].data.shape[2], a[0].data.shape[3])

# -----------------------------------------------
#  Set up Figure
# -----------------------------------------------
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

from matplotlib.pyplot import cm
color = iter(cm.rainbow_r(np.linspace(0, 1, len(p_list))))

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
    ax.errorbar(x=[ra_px], y=[dec_px], xerr=[ra_err], yerr=[dec_err],
                marker='+', ms=8, mec=c, mew=1.25, zorder=3.1, ecolor=c,
                capsize=5, capthick=2, lw=2, label=str(i))
#     print ra_px, dec_px

# matplotlib.rcParams['legend.handlelength'] = 0
# matplotlib.rcParams['legend.markerscale'] = 0
# ax.legend(loc=2, bbox_to_anchor=(0.9, 0.9), borderaxespad=0., fontsize=15, numpoints=1)
plt.tight_layout()
plt.show()
a.close()
# User_input = raw_input('Save figure? (Y/N): ')
# if User_input == 'Y':
#     filename = "veloGradient_markers.eps"
#     fig.savefig(Plotpath + filename, dpi=100,
#                 bbox_inches="tight", pad_inches=0.1)
#     print "-- Saved figure as : %s --" % (Plotpath + filename)

theta_rad = abs(np.arctan(dec_coord_ls[0] - dec_coord_ls[-1]) / (
    (ra_coord_ls[0] - ra_coord_ls[-1]) * np.cos(np.mean(dec_coord_ls))))
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
print best

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
    s_sq = (errfunc(pfit, xdata, ydata, yerr)**2).sum() / (len(ydata) - len(x0))
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
quasi_chisq = out.res_var
# if quasi_chisq < 1.0 :
#     uncertainty = uncertainty/np.sqrt(quasi_chisq)
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
    fit = fig.add_subplot(211)       #, adjustable='box', aspect=1.2)
    fit.set_xticklabels(())
    plt.ylabel("Dec")

    fit.plot(xdata, ydata, 'ro', xdata, fun2(param, xdata))
    fit.errorbar(xdata, ydata, xerr=xerr, yerr=yerr, fmt='r+')
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
    # separate plot to show residuals
    residuals = fig.add_subplot(212, adjustable='box', aspect=0.4)
    residuals.errorbar(x=xdata, y=residual, yerr=adjusted_err,
                                fmt="r+", label="Residuals")
    # make sure residual plot has same x axis as fit plot
    residuals.set_xlim(fit.get_xlim())
    plt.axhline(y=0, color='b')
    plt.xlabel("RA ")
    # These data look better if 'plain', not scientific, notation is used, and if
    #   the tick labels are not offset by a constant (as is done by default).
    plt.ticklabel_format(style='plain', useOffset=False, axis='x')
    plt.ylabel("Residuals")
    residuals.grid()
    plt.tight_layout()
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
RA_major_err = np.array([p_list_err[i][0] / CellSize for i in range(len(p_list_err))])
Dec_major_err = np.array([p_list_err[i][1] / CellSize for i in range(len(p_list_err))])

def fitted_err_on_y(deltaP, x):
    """ calc y_error of point on fit

    """
    deltam, deltab = deltaP
    dy = np.sqrt((x * deltam)**2 + (deltab)**2)
    return dy

def calc_dist(xloc, yloc, delta_x, delta_y):
    """ calc distance in arcsec offset along major axis

    Parameters
    ----------
    xloc: array
        RA in deg of source position
    yloc: array
        dec in deg of source position along the fitted major axis

    Return
    ------
    offset: array
        in arcsec along the major axis, w.r.t the central emission position

    """
    diff_x = np.array([xloc[i] - xloc[4] for i in range(len(xloc))])
    diff_y = np.array([yloc[i] - yloc[4] for i in range(len(yloc))])
    offset_sq = (diff_x * np.cos(np.mean(yloc/3600. * np.pi/180.)))**2 + diff_y**2
    R = np.sqrt(offset_sq)
    offset_err = np.sqrt((diff_x*delta_x)**2 + (diff_y*delta_y)**2) / R

    # put error bar on center point as well, in arcsec
    xcenter_err, ycenter_err = p_list_err[4]
    offset_err[np.where(np.isnan(offset_err))[0][0]] = np.sqrt(xcenter_err**2 +ycenter_err**2)

    return R, offset_err


# youfitErr = fitted_err_on_y(uncertainty, xdata)     # very big..., because of x
# use uncertainty on ydata instead
# reasonble because people slice along their map, they don't fit for the major axis

offset, offset_err = calc_dist(RA_major*deg_to_arcsec, Dec_major*deg_to_arcsec, RA_major_err, Dec_major_err)

# to show blue, red on left and right of central
offset[5:] = -offset[5:]
off = list(offset)
offset_err = list(offset_err)
off.pop(1)
z.pop(1)
offset_err.pop(1)
z_err = 21.5*5/2.     # range of channels combined to make models/2. --> error bar
if plotPV:
    plt.errorbar(off, z, yerr=z_err, fmt='r+', xerr=offset_err)
    plt.ylabel('v_r = v sin i')
    plt.xlabel(' offset from line center position ["] ')
    plt.title('PV along major axis at PA {:.2f} deg'.format(PA_deg))
    plt.tight_layout()
    plt.show()

offset, offset_err = calc_dist(RA_major*deg_to_arcsec, Dec_major*deg_to_arcsec, RA_major_err, Dec_major_err)
off = list(offset)
offset_err = list(offset_err)
off.pop(1)
offset_err.pop(1)
if plotRot:
    plt.errorbar(off[3], z[3], yerr=z_err, fmt='go', label='line center')
    plt.errorbar(off[:3], z[:3], yerr=z_err, fmt='ro', label='red', xerr=offset_err[:3])

    plt.ylabel('|v_r = v sin i|')
    plt.errorbar(off[4:], abs(np.array(z[4:])), yerr=z_err, fmt='bo', label='blue', xerr=offset_err[4:])
    plt.xlabel(' radial offset from line center position ["] ')
    plt.title('"Rotation Curve" taken along major axis at PA {:.2f} deg'.format(PA_deg))
    plt.ylim(-50, 450)
    plt.xlim(-0.2, 1.2)
    plt.legend()
    plt.tight_layout()
    plt.show()

# --------------------------------------------------------------
#  disk model
# ---------------------------------------------------------------
kpc_to_m = 3.086e+19


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
    from astropy.cosmology import WMAP9
    arcsec_to_kpc = WMAP9.kpc_proper_per_arcmin(z).value / 60.
    return arcsec_to_kpc * R_arcsec


def arctang2(p, R_kpc):
    vcsini, r_s, c = p
    # hacky step to force r_s to be physical
    if r_s > 0.1:
        return vcsini*(2./np.pi)*np.arctan(r_s*R_kpc) + c
    else:  # bascially reject this fit
        return 1e-5

ydata = abs(np.array(z))
yerr = z_err
xdata = offset_to_physicalR(np.array(off), redshift)     # np.array(off)
xerr = offset_to_physicalR(np.array(offset_err), redshift)      # np.array(offset_err)

# plot velocity v.s. physical radius
plt.errorbar(xdata, ydata, yerr=yerr, fmt='ko', xerr=xerr)
plt.ylabel('V sin i')
plt.xlabel('Physical distance from line center [kpc]')
plt.xlim(-1, 8.5)
plt.ylim(-50, 450)
plt.tight_layout()
plt.show()

# Least sq fit tangent model to PV along major axis
x4 = [30, 1, 1]
xdataPV = list(xdata[:4])
_blue = -xdata[4:]
for i in range(len(_blue)): xdataPV.append(_blue[i])
xdataPV = np.array(xdataPV)
errfunc = lambda p, x, y, err: (y - arctang2(p, x)) / err
_pfit4 = optimize.leastsq(errfunc, x4, args=(xdataPV, np.array(z), yerr), full_output=1, maxfev=10000)
plt.errorbar(xdataPV, np.array(z), xerr=xerr, yerr=yerr, fmt='ko', linestyle='None', label='data')
# extrapolate for prettier looking
plt.plot(np.linspace(xdataPV.min(), xdataPV.max(), 100), arctang2(_pfit4[0], np.linspace(xdataPV.min(), xdataPV.max(), 100)), 'r--', label='no fit to sin i, tangent 2')
plt.xlabel("Radial offset from line center position kpc")
plt.ylabel("V [km/s]")
plt.legend(loc='best')
plt.xlim(-8, 8)
plt.ylim(-425, 425)
plt.tight_layout()
plt.show()
print("Bestfit Circ velocity * sin i: {} km/s").format(_pfit4[0][0])

# ----
# fit i altogether with least-sq
def inner_rot_incl_iter(p, x, y):
    """ solve for i on-the-fly """

    a, b, i = p

    def _calc_dist_only(xloc, yloc, i):
        """ calc distance in arcsec offset along major axis, taking into acct inclination
        """
        diff_x = np.array([xloc[j] - xloc[3] for j in range(len(xloc))])
        diff_y = np.array([yloc[j] - yloc[3] for j in range(len(yloc))])
        offset_sq = (diff_x * np.cos(np.mean(yloc/3600. * np.pi/180.)))**2 + diff_y**2
        R = np.sqrt(offset_sq)
        return R

    R = offset_to_physicalR(_calc_dist_only(x, y, i=i), redshift)
    V_obs = a * R * np.sin(i * np.pi/180.) + b
    return V_obs


def calc_dist_only(xloc, yloc, i):
    """ calc distance in arcsec offset along major axis, taking into acct inclination
    """
    diff_x = np.array([xloc[j] - xloc[3] for j in range(len(xloc))])
    diff_y = np.array([yloc[j] - yloc[3] for j in range(len(yloc))])
    offset_sq = (diff_x * np.cos(np.mean(yloc/3600. * np.pi/180.)))**2 + diff_y**2
    R = np.sqrt(offset_sq)
    return R


# copy list
_ra = list(RA_major*deg_to_arcsec)
_dec = list(Dec_major*deg_to_arcsec)
_ra.pop(1)                # remove source two in red
_dec.pop(1)
xidata = np.array(_ra)

p0 = [32, 0.5, 30.]
errfunc = lambda p, x, dec, y, err: (y - inner_rot_incl_iter(p, x, dec)) / err
pfit, pcov, infodict, errmsg, success = optimize.leastsq(errfunc, p0, args=(xidata, np.array(_dec), ydata, yerr), full_output=1)

R = offset_to_physicalR(calc_dist_only(xidata, np.array(_dec), pfit[-1]), redshift)
print("\n  best-fit inclination: {} deg \n").format(pfit[-1])
plt.errorbar(R, ydata, yerr, fmt='ko', linestyle='None', label='data')
plt.plot(R, inner_rot_incl_iter(pfit, xidata, np.array(_dec)), 'r--', label='fit i on the fly')
plt.xlabel("Radial offset from line center position kpc")
plt.ylabel("V [km/s]")
plt.legend(loc='best')
plt.xlim(-1, 8)
plt.ylim(-50, 425)
plt.tight_layout()
plt.show()

# ---
def M_encl(R_kpc, v_sini):
    """ calc mass enclosed """     # v in km/s
    from astropy.constants import G

    R = kpc_to_m * R_kpc
    Msun = 1.989e30      # kg
    v2 = (v_sini * 1e3)**2
    v2oG = v2 / G.value
    v2oGoMsun = v2oG / Msun
    M = R * v2oGoMsun           # per Msun
    M10 = M / 1e10              # in units of x 1e10 Msun
    return M10

# calc. mass within some physical dist from line center
R_in_kpc2 = R.max()              # using best-fit incl
M_model_iter_inc = M_encl(R_in_kpc2, inner_rot_incl_iter(pfit, xidata, np.array(_dec)).max())               # could be dangerous to just use .max() on velocity, but don't want to mess with functions now.
print("Mass enclosed within {:.2f} kpc: {:.2f} x 1e+10 Msun".format(R_in_kpc2, M_model_iter_inc))

# plot M rise with R
if plotMdyn:
    plt.plot(R, M_encl(R, inner_rot_incl_iter(pfit, xidata, np.array(_dec))), label='R dep. on i')
    plt.legend()
    plt.tight_layout()
    plt.show()
