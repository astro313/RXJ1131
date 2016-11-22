import aplpy       # http://aplpy.readthedocs.org/en/stable/howto_subplot.html
import numpy as np


def get_RA_dec(f):
    """
    Get RA, dec from fits

    Input:
    f: str
        filename with .fits extension
    """
    import atpy
    t = atpy.Table(f)
    return t.ra, t.dec

def sigma_contour_VLA(sigma):
    """
    Return a list of contour levels, (4, 8, 16, 32..)*sigma

    +/- 2^(n)*sigma; n from 2

    Parameters
    ----------
    sigma: float
        map rms

    Returns
    -------
    arr: list
    contour levels

    """
    arr = range(11)
    arr.remove(0)
    arr.remove(1)
#     arr = [2 ** i * np.sqrt(2) * sigma for i in arr]
    arr = [2 ** i * sigma for i in arr]
    arr.append(list(np.array(arr) * -1)[0])
    return arr



def sigma_contour_tenpercent(Smax):
    """
    return list of sigma multiples, +/- 10% max * n from 1 ....4

    increment in 15% of max
    """
    arr = [0.10 * Smax]
    for i in np.arange(1, 7):
        arr.append(0.15 * Smax * i + arr[0])
    Narr = -1 * arr
    arrE = np.hstack((Narr, arr))
    return arrE



def sigma_contour_array(sigma):
    """
    return list of sigma multiples, +/- 3*n*sigma from 1 ....14
    """
    arr = np.arange(15) * 3
    arr = arr[1:] * sigma
    Narr = -1 * arr
    arrE = np.hstack((Narr, arr))
    return arrE


def sigma_contour_CARMA(sigma):
    """
    return list of sigma multiples, -4,-3,-2,2,3...20  * sigma
    """
    arr = range(-4, 20)
    arr.remove(0)
    arr.remove(-1)
    arr.remove(1)
#    arr = [i * np.sqrt(2) * sigma for i in arr]
    arr = [i * sigma for i in arr]
    return arr


def standard_plot_setup(sp, ra_center, dec_center, size, tickc='black'):
    sp.set_frame_color('black')
    sp.frame.set_linewidth(1)
    sp.set_system_latex(True)
    sp.recenter(x=ra_center, y=dec_center, radius=size)   # width = blah, height = blah     # Scaling and panning

    sp.tick_labels.set_font(size='x-large', weight='bold')      # size='10'
#    sp.tick_labels.set_style('colons')
#     sp.tick_labels.set_xformat('hh:mm:ss')
#     sp.tick_labels.set_yformat('dd:mm:ss')
    sp.ticks.set_color(tickc)
    sp.ticks.set_linewidth(2)
    sp.ticks.set_length(10)
#     sp.ticks.set_minor_frequency(4)
    # sp.ticks.set_xspacing(45*15/3600.)        # deg

    sp.axis_labels.set_font(size='xx-large', weight='bold')      # (size='12')
    sp.axis_labels.set_xtext('Right Ascension (J2000)')
#    sp.axis_labels.set_xtext(u'$\\mathrm{Right}$'+' '+u'$\\mathrm{Ascension}$' + ' ' + u'$\\mathrm{(J2000)}$')
    sp.axis_labels.set_ytext('Declination (J2000)')
#    sp.axis_labels.set_ytext(u'$\\mathrm{Declination}$'+' '+u'$\\mathrm{(J2000)}$')
    sp.axis_labels.set_xpad(5)
    sp.axis_labels.set_ypad(-25)
    sp.axis_labels.set_xposition('bottom')


def inset_frmt(sp):
    """
    Hide the ticks, tick labels and axis labels (e.g. RA & Dec.)
    """
    sp.ticks.hide()
    sp.tick_labels.hide()
    sp.axis_labels.hide()


def setup_beam(sp, loc='bottom left', c='black', hatch=None, a=0.8, lw=3):
    """
    Auto fetch beam info from header
    Inputs:
    --------
    loc: str
    hatch: str
        one of /,|,-,+,x,o,O,.,*
    """
    sp.show_beam(corner=loc, color=c, fill=True, edgecolor='grey', facecolor='black', linestyle='solid', linewidth=lw, frame=True, alpha=a)
    if hatch is not None:
        sp.beam.set_hatch(hatch)


def setup_beam2(sp, major, minor, angle, idx, loc='top right', c='black', hatch=None, a=0.8, lw=3):
    """
    Manual input beam info

    Input:
    -------
    major: float
        deg
    minor: float
        deg
    PA: float
        deg
    idx: int
        since more than 1 beam, sp.beam becomes a list
    """
    idx = int(idx)
    sp.add_beam(0, 0, 0)
    sp.beam[idx].show(major=major, minor=minor, angle=angle, corner=loc, color=c, fill=True, edgecolor='grey', facecolor='black', linestyle='solid', linewidth=lw, frame=True, alpha=a)
    sp.beam[idx].set_frame(True)


def markers_cross(sp, ra, dec, layer=None, ec='yellow', fc='none', mk='+', s=500, a=1.0, lw=1.2, zorder=50):
    """
    Inputs:
    -------
    ra: float
    dec: float
    layer: str
    """
    sp.show_markers(ra, dec, layer=layer, edgecolor=ec, facecolor=fc, marker=mk, s=s, alpha=a, linewidth=lw, zorder=50)


def put_label(sp, x, y, text, layer, c='yellow', s='x-large', w='bold'):
    """
    Inputs:
    -------
    x: float
    y: float
    text: str
    layer: str
    """
    sp.add_label(x, y, text, relative=True, color=c, size=s, layer=layer, weight=w)


def frmt_colorbar(sp):
    """
    format colorbar
    """
    sp.colorbar.set_axis_label_font(size=15)
    sp.colorbar.set_width(0.1)
    sp.colorbar.set_axis_label_pad(15)
    sp.colorbar.set_axis_label_text('$Flux Density [mJy]$')


def setup_scalebar(sp, lg, label, unit=None, loc='bottom right', a=0.9, c='white', lw=3):
    """
    unit: str
        unit for label
    lg: float
        length [degrees]
    color: str
        white, lime, Red, DogderBlue, yellow, black, etc
    label: str
    """
    sp.show_scalebar(lg, color=c, corner=loc, alpha=a, linestyle='solid', linewidth=lw)
    sp.scalebar.set_font(size='large', weight='bold')
    # if unit == None:
    #     unit = r'$kpc / {\prime\prime}$'
    sp.scalebar.set_label(label)