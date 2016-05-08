'''

plot source locations from lens model of various channel as markers on observed 1st moment map

Last Modified: 08 May 16

History:
08 May 16:
  - calc. enclosed mass using velocity from fit
  - fit to rotation curve, just solid body
  - plot PV, and rotation curve along major axis
07 May 16:
  - fit through data to get `major axis` using fmn, leastsq, ODR
05 May 16:
  - added error bars to positions
04 May 16:
  - created code

Note:
- two things:
    1. rotation curve looks like solid body, even up to 10 kpc, seems very unexpected, recent stellar light is 8 kpc across, if the redmost and bluemost channel are truely ~10 kpc from center, should see more than just solid body part of the rotation curve? am I calculating the separation incorrectly (larger than true?)
    2. on the other hand, physical dist. from center also seem reasonable, kind of comforting that our model is giving some information
- any derived mass in literature to compare with, Mbh ~ 2e8 Msun
- can we fit SED to get stellar mass?

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
print("PA along the `major axis`: {} deg.").format(PA_deg)

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
def fun2(p, x):
    m, b = p
    return m * x + b
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
    fit = fig.add_subplot(211) #, adjustable='box', aspect=1.2)
    fit.set_xticklabels(())
    plt.ylabel("Dec")

    fit.plot(xdata, ydata, 'ro', xdata, fun2(param, xdata))
    fit.errorbar(xdata, ydata, xerr=xerr, yerr=yerr, fmt='r+')
    fit.set_yscale('linear')
    #   draw starting guess as dashed green line ('r-')
    # fit.plot(xdata, fun2(x0, xdata), 'g-', label="Start", linestyle="--")
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
    plt.show()

youtfit = fun2(out.beta, xdata)
# replaced by the above
# plt.figure()
# plt.errorbar(ra_px_ls, dec_px_ls, xerr=ra_err_ls, yerr=dec_err_ls,
#              marker='+', ms=8, mew=1.25, capsize=5, capthick=2, linestyle='None')
# plt.plot(xdata, youtfit, 'k--', lw=1.2)
# plt.show()


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
    plt.show()


# -------------------------------------------------------
# PV and rot. curve exclude 2nd component in redmost channel
# -------------------------------------------------------
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
    offset_sq = (diff_x * np.cos(np.mean(yloc * np.pi/180.)))**2 + diff_y**2
    offset_err = 2.*np.sqrt((diff_x*delta_x)**2 + (diff_y*delta_y)**2)
    return np.sqrt(offset_sq), offset_err

# youfitErr = fitted_err_on_y(uncertainty, xdata)     # very big..., probably because of the error on offset, which probably means we shouldn't have just fit a linear func to the coordinates

# this is hacky, but reasonable: use uncertainty on ydata instead
# reasonble because people slice along their map, they don't fit for the major axis
offset, offset_err = calc_dist(xdata, youtfit, ra_err_ls, dec_err_ls)
print offset_err
# to show blue, red on left and right of central
offset[5:] = -offset[5:]
off = list(offset)
offset_err = list(offset_err)
off.pop(1)
z.pop(1)
offset_err.pop(1)
# error bar?
z_err = 21.5     # bandwidth
if plotPV:
    plt.errorbar(off, z, xerr=offset_err, yerr=z_err, fmt='r+')
    plt.ylabel('v_r = v sin i')
    plt.xlabel(' offset from line center position ["] ')
    plt.title('PV along major axis at PA {:.2f} deg'.format(PA_deg))
    plt.show()

offset, offset_err = calc_dist(xdata, youtfit, ra_err_ls, dec_err_ls)
off = list(offset)
offset_err = list(offset_err)
off.pop(1)
offset_err.pop(1)
if plotRot:
    plt.errorbar(off[3], z[3], yerr=z_err, fmt='go', label='line center')
    plt.errorbar(off[:3], z[:3], xerr=offset_err[:3], yerr=z_err, fmt='ro', label='red')
    # errorbar
    plt.ylabel('|v_r = v sin i|')
    plt.errorbar(off[4:], abs(np.array(z[4:])), xerr=offset_err[:3], yerr=z_err, fmt='bo', label='blue')
    plt.xlabel(' radial offset from line center position ["] ')
    plt.title('"Rotation Curve" taken along major axis at PA {:.2f} deg'.format(PA_deg))
    plt.ylim(-50, 450)
    plt.xlim(-0.5, 4.0)
    plt.legend()
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


def inner_rot(p, R_kpc):
    """ solid body """

    a, b = p
    R = R_kpc   # * kpc_to_m
    V_sini = a * R + b
    return V_sini


ydata = abs(np.array(z))
yerr = z_err
xdata = offset_to_physicalR(np.array(off), 0.6357)     # np.array(off)
xerr = offset_to_physicalR(np.array(offset_err), 0.6357)      # np.array(offset_err)

# plot velocity v.s. physical radius
plt.errorbar(xdata, ydata, xerr=xerr, yerr=yerr, fmt='ko')
plt.ylabel('V sin i')
plt.xlabel('Physical distance from line center [kpc]')
plt.xlim(-1, 13)
plt.ylim(-50, 450)
plt.show()
print(" ok, the physical distance is somewhat consistent with HST")

# fit to rotation curve
data = RealData(xdata, ydata, sy=yerr)
x0 = [32, 0.5]

_model_1 = Model(inner_rot)

def run_odr_for_model(model, func, guess, x, y, xerr, yerr):
    # getting v sin i
    odr = ODR(data, model, beta0=guess)
    out = odr.run()
    param = out.beta
    cov = out.cov_beta
    uncertainty = out.sd_beta

    # Calculate initial residuals and the 'adjusted error' for each data point
    delta = out.delta  # estimated x-component of the residuals
    eps   = out.eps    # estimated y-component of the residuals
    xstar = (xerr*np.sqrt(((yerr*delta)**2) / ((yerr*delta)**2 + (xerr*eps)**2)))
    ystar = (yerr*np.sqrt(((xerr*eps)**2) / ((yerr*delta)**2 + (xerr*eps)**2)))
    adjusted_err = np.sqrt(xstar**2 + ystar**2)
    residual = np.sign(ydata-func(param, xdata))*np.sqrt(delta**2 + eps**2)

    ## Plot
    fig = plt.figure(facecolor="0.98")   # with light gray background
    fig.subplots_adjust(hspace=0)
    fit = fig.add_subplot(211)      # , adjustable='box', aspect=1.2)
    fit.set_xticklabels(())
    plt.ylabel(" |V_r sin i| ")

    fit.plot(xdata, ydata, 'ro', xdata, inner_rot(param, xdata))
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
    fit.set_xlim(-1, 12)
    fit.set_ylim(-50, 450)
    fit.grid()
    # separate plot to show residuals
    residuals = fig.add_subplot(212)   # , adjustable='box', aspect=0.4)
    residuals.errorbar(x=xdata, y=residual, yerr=adjusted_err,
                                fmt="r+", label="Residuals")
    # make sure residual plot has same x axis as fit plot
    residuals.set_xlim(fit.get_xlim())
    plt.axhline(y=0, color='b')
    plt.xlabel("Radial offset from line center position kpc")
    # These data look better if 'plain', not scientific, notation is used, and if
    #   the tick labels are not offset by a constant (as is done by default).
    plt.ticklabel_format(style='plain', useOffset=False, axis='x')
    plt.ylabel("Residuals")
    residuals.grid()
    plt.show()
    return param, uncertainty

# run
param1, p1_err = run_odr_for_model(_model_1, inner_rot, x0, xdata, ydata, xerr, yerr)


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
R_in_kpc = 10.
M = M_encl(R_in_kpc, inner_rot(param1, R_in_kpc))   # use velo from fit
print("Mass enclosed within {:.2f} kpc: {:.2f} x 1e+10 Msun".format(R_in_kpc, M))