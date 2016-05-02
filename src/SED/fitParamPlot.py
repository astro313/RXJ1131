'''
Plot mbb_emcee output

Last Modified: 01 May 2016

Original: /Users/admin/Research/3C220.3/3C220.3_CARMA/writeup/Code

History:
01 May 2016
    get rid of np.flipud() in ax.imshow in order for triplot to work correctly.. Don't know why it works for 3C220.3, but has issues here -- maybe MPL version different?
30 Apr 2016
    comment out code for ALMA prediction, set ylim for best-fit SED from min=0 to avoid display -ve due to upper limit data points
16 Nov 2015
    added function PlotSED, changed savePlot location

Note:

'''


import mbb_emcee
import matplotlib.pyplot as plt
import matplotlib
import sys
from pylab import savefig
import numpy as np
import scipy.ndimage.filters as filter

# matplotlib.rc('text', usetex=True)
matplotlib.rc('xtick', labelsize=18)
matplotlib.rc('ytick', labelsize=18)
font = {'family': 'Arial Narrow',
        'weight': 'bold',
        'size': 12}        # 15
matplotlib.rc('font', **font)
# monospace sans-serif
axes = {
    'titlesize': 'xx-large',
    'labelsize': 'x-large',
    'labelweight': 'normal',
    'linewidth': 1.5
}
matplotlib.rc('axes', **axes)


def convergence(plotTitle='ln prob', saveFig=False):
    """
    Plot the convergence profile.  I.e., Max(lnprob) - lnprob as a function of
    iteration number.
    """

    lnprob = res.lnprobability.flatten()
    Plotlnprob = max(lnprob) - lnprob
    plt.clf()
    plt.plot(Plotlnprob)  # , ',', alpha=0.5)
    plt.xlabel('iteration')
    plt.ylabel('max(lnprob) - lnprob')
    plt.title(plotTitle)
    plt.semilogy()
    if saveFig == True:
        outfile = 'convergence__'
        savefig('../Figures/' + outfile + filename.replace('.h5', '.png'))
    else:
        plt.show()


def chain_iter(nburn=100, saveFig=False):
    """
    plot the chains to visually assess convergence
    nburn default = 100, mbb_emcee default=50; modify as needed
    """
    par = ('T', 'beta', 'lambda0', 'alpha', 'fnorm')
    nwalkers = res._nwalkers
    plt.figure(figsize=[20, 10])
    for i, p in enumerate(par):
        ncols = (len(par)/2 + 1) if len(par)%2 != 0 else len(par)/2
        plt.subplot(len(par) / 2, ncols, i + 1)
        for w in range(nwalkers):
            plt.plot(
                np.arange(res.chain.shape[1]), res.chain[w, :, i], 'b-', alpha=0.1)
        plt.ylabel(p)
        plt.xlabel('Iter')
        aymin, aymax = plt.ylim()
        plt.vlines(nburn, aymin, aymax, linestyle=':')
        plt.ylim(aymin, aymax)
    plt.suptitle('performance of each paramater as a function of iter', fontsize=14, y=0.90)
    if saveFig == True:
        outfile = 'steps_convergence__'
        savefig('../Figures/' + outfile + filename.replace('.h5', '.eps'), dvi=600)
    else:
        plt.show()


def Ramdon():
    """
    Grab a set of values for the parameters randomly from the sampled result

    """
    return res.choice()


def extrapolateFlux(wave_um):
    """
    Purpose
    -------
        return best-fit SED flux density at some given wavelength

    Input
    -----
    wave_um: float
        wavelength of interest in micron
    """
    return res.best_fit_sed(wave_um)


def get_Paramvalues(percentile=68.3):
    """
    Input:
    ------
    percentile: +/- return of all parameters


    Print:
    --------
    5 element list of parameters (T, beta, alpha, fnorm, lambda0)
    4 elements (T, beta, alpha, fnorm) if op thin
    list of 3, [mean = central confidence, +errorUp, -errorLow]
    default percentile =  68.3
    The returned is not the max. ln prob value, only the mean
    For best value: see BestFit()

    Returns:
    --------
    None at the moment
    """

    paramsDict = {'T': 'K',
                  'beta': ' ',
                  'alpha': ' ',
                  'fnorm': 'mJy',
                  'lambda0': 'um'}
    if res._opthin:
        paramsDict.pop('lambda0')
    if res._noalpha:
        paramsDict.pop('alpha')

    for k, v in paramsDict.iteritems():
        param = k
        param_val = res.par_cen(param, percentile=percentile)
        unit = v
        print("{0:s}: Mean {2:0.2f}+{3:0.2f}-{4:0.2f} [{1:s}]".format(param, unit, *param_val))


def get_computation_value(FIR=False, percentile=68.3):
    """
    Purpose:
    --------
    Print peak lambda (p_val [mean, +errorUp, -errorLow])
    Returns, LIR, (LFIR), dustMass

    Input:
    ------
    Keywords:
    FIR for computing FIR luminosity; default: wavemin=42.5, wavemax=122.5
    percentile for computing the value at a different confidence interval


    Returns:
    ---------
    LIR
    LFIR    # if FIR == True
    Dust Mass


    Note:
    ------
    Most quantities are computed before saving results to .h5

    However, if one wish to get an integrate Luminosity of different range of wavelengths, that will require computing in the function, not available from result.
    """

    p_val = res.peaklambda_cen(percentile=percentile)
    print("Peak Obs wavelength: {:0.1f}+{:0.1f}-{:0.1f} [um]".format(*p_val))
    # fetch from computed, integrated
    lir = res.lir_cen(percentile=percentile)
    args = (res._lir_min, res._lir_max) + tuple(lir)
    lirstr = "L_IR({:0.1f} to {:0.1f}um): {:0.2f} "\
        "+{:0.2f} -{:0.2f} [10^12 L_sun]\n"
    print(lirstr).format(*args)
    IR_ = {'lir': lir, 'lirChain': res.lir_chain}    # Avoid over-written by FIR

    if FIR == True:
        res.compute_lir(wavemin=42.5, wavemax=122.5)    # compute from chain
        RangeMSG = "Current range of wavelength for computing LIR Luminosity: {0:.2f} - {1:.2f} um"
        print RangeMSG.format(*res.lir_wavelength)
        FIR_val = res.lir_cen(percentile=percentile)
        # _lir_min changed to 42.5
        args = (res._lir_min, res._lir_max) + tuple(FIR_val)
        lirstr = "L_FIR({:0.1f} to {:0.1f}um): {:0.2f} "\
            "+{:0.2f} -{:0.2f} [10^12 L_sun]\n"
        print(lirstr).format(*args)
        IR_['FIRchain'], IR_['LFIR'] = res.lir_chain, FIR_val
    return IR_


def Plot_Standard(saveFig=False):
    """
    Purpose:
    --------
    Make 4 panels plot of
    - Data + errorbar + best-fit SED (unnormalized)
    - Likelihood of T [Temperature rest Frame]
    - Likelihood of beta
    - Likelihood of M_dust
    """
    wave, flux, flux_unc = res.data
    redshiftZ = res.redshift
    f, axs = plt.subplots(nrows=2, ncols=2, figsize=(12, 8))

    ax1 = axs[0, 0]
    p_data = ax1.errorbar(wave, flux, yerr=flux_unc, fmt='ro')
    p_wave = np.linspace(wave.min() * 0.5, wave.max() * 1.5, 200)   # range of x
    # plot best fit using the params in res
    p_fit = ax1.plot(p_wave, res.best_fit_sed(p_wave), color='blue')
    ax1.set_xlabel('Wavelength [um]')
    ax1.set_ylabel('Flux Density [mJy]')
    ax1.set_title('Best-fitted SED')
#     ax1.set_ylim([-100, max(flux_unc)*1.25])

    ax2 = axs[0, 1]
    h1 = ax2.hist(res.parameter_chain('T/(1+z)') * (1 + redshiftZ))
    ax2.set_xlabel(r'$T_{d, rest}$')
    ax2.set_ylabel(r'$\cal L$')
    ax2.set_yticks([])

    ax3 = axs[1, 0]
    h2 = ax3.hist(res.parameter_chain('beta'))
    ax3.set_xlabel(r'$\beta$')
    ax3.set_ylabel(r'$\cal L$')
    ax3.set_yticks([])

    ax4 = axs[1, 1]
    h4 = ax4.hist(res.dustmass_chain)
    ax4.set_xlabel('Dust Mass  ' + r'$[10^8 M_{\odot}]$')
    ax4.set_ylabel(r'$\cal L$')
    ax4.set_yticks([])

    if saveFig == True:
        outfile = 'SED_ParamProb__'
        savefig('../Figures/' + outfile + filename.replace('.h5', '.eps'), dvi=600)
    else:
        f.show()


def PlotSED(saveFig=False):

    """
    Purpose
    --------
    Make a plot of
    - Data + errorbar + best-fit SED (unnormalized)
    """

    wave, flux, flux_unc = res.data
    redshiftZ = res.redshift
    fig = plt.figure(figsize=(5, 5))
    fig.subplots_adjust(left=0.18, bottom=0.12, top=0.90, right=0.95)

    p_data = plt.errorbar(wave, flux, yerr=flux_unc, fmt='ro')
    p_wave = np.linspace(wave.min() * 0.5, wave.max() * 1.25, 200)   # range of x
    # plot best fit using the params in res
    p_fit = plt.plot(p_wave, res.best_fit_sed(p_wave), color='blue')
    plt.ylim(-100, 500)
    plt.xlabel('Observed Frame Wavelength [um]')
    plt.ylabel('Flux Density [mJy]')
    plt.title('Best-fitted SED')

    if saveFig == True:
        outfile = 'SED_ParamProb__'
        savefig('../Figures/' + outfile + filename.replace('.h5', '.eps'), dvi=600)
    else:
        plt.show()


def Plot_Chain(IR_chain, saveFig=False, FIR=False):

    f, axs = plt.subplots(nrows=2, ncols=2, figsize=(12, 8))

    ax1 = axs[0, 0]
    h1 = ax1.hist(res.peaklambda_chain)
    ax1.set_xlabel('Peak Obs Wavelength [um]')
    ax1.set_ylabel(r'$\cal L$')
    ax1.set_yticks([])

    ax2 = axs[0, 1]
    h2 = ax2.hist(IR_chain['lirChain'])
    ax2.set_xlabel(r'$L_{IR} [10^{12}L_{\odot}]$')
    ax2.set_ylabel(r'$\cal L$')
    ax2.set_yticks([])

    # ax3 = axs[1, 0]
    # h3 = ax3.hist(IR_chain['lirChain'])
    # ax3.set_xlabel(r'$L_{IR}\ \ [10^{12}L_{\odot}]$')
    # ax3.set_ylabel(r'$\cal L$')
    # ax3.set_yticks([])

    if FIR == True:
        ax4 = axs[1, 1]
        h3 = ax4.hist(IR_chain['FIRchain'])
        ax4.set_xlabel(r'$L_{FIR} [10^{12}L_{\odot}]$')
        ax4.set_ylabel(r'$\cal L$')
        ax4.set_yticks([])

    if saveFig == True:
        outfile = 'ParamCHAINProb__'
        savefig('../Figures/' + outfile + filename.replace('.h5', '.png'))
    else:
        f.show()


def MCMC_scatterDots(Xparam='T', Yparam='beta', saveFig=False):
    """
    Purpose:
    --------
    Plot scattered chain samples of 2 parameters from MCMC results
    Similar to hist2D

    Input:
    ------
    Xparam: Observed frame dust temeprature
    Yparam: beta
    options: 'fnorm', 'beta', 'alpha', 'lambda_0', 'T'

    """
    X = res.parameter_chain('T')
    Y = res.parameter_chain('beta')

    # Define plot ranges
    Xrange = [min(X) - 0.2, max(X) + 0.2]
    Yrange = [min(Y) - 0.2, max(Y) + 0.2]

    # Define figure size and formatting
    fig = plt.figure(figsize=(7, 7))
    fig.subplots_adjust(left=0.10, bottom=0.09, top=0.98, right=0.98)
    plt.plot(X, Y, 'o', ms=4, alpha=.1, color='b')

    # set plot range, axes ticks, and axes labels
    plt.xlim(Xrange)
    plt.ylim(Yrange)
    plt.xlabel(Xparam)
    plt.ylabel(Yparam)
    plt.xticks(fontsize=16)
    plt.yticks(fontsize=16)
    # plt.gcf()
    if saveFig == True:
        outfile = 'ScatterPlot' + Xparam + Yparam + '__'
        savefig('../Figures/' + outfile + filename.replace('.h5', '.png'))
    else:
        plt.show()


def getsigmalevels(hist2d):
    """
    get sigma for samples contour plot

    Input:
    ------
    hist2d formed by the 2 parameters
    """
    sigma1 = 0.68268949
    level1 = 0
    sigma2 = 0.95449974
    level2 = 0
    sigma3 = 0.99730024
    level3 = 0

    lik = hist2d.reshape(hist2d.size)
    sortlik = np.sort(lik)

    # 1sigma level
    dTotal = np.sum(sortlik)
    nIndex = sortlik.size
    dSum = 0
    while (dSum < dTotal * sigma1):
        nIndex -= 1
        dSum += sortlik[nIndex]
    level1 = sortlik[nIndex]

    # 2 sigma level
    nIndex = sortlik.size
    dSum = 0
    while (dSum < dTotal * sigma2):
        nIndex -= 1
        dSum += sortlik[nIndex]
    level2 = sortlik[nIndex]

    # 3 sigma level
    nIndex = sortlik.size
    dSum = 0
    while (dSum < dTotal * sigma3):
        nIndex -= 1
        dSum += sortlik[nIndex]
    level3 = sortlik[nIndex]
    return level1, level2, level3


def makesubplot1d(ax, samples, weights=None, interpolate=False, smooth=True,
                  label=None, bins=30, range=None, color='k'):
    """
    Make histogram of samples
    """
    import scipy.interpolate as interp

    if range is None:
        hist, xedges = np.histogram(samples, bins, normed=True, weights=weights)
    else:
        hist, xedges = np.histogram(
            samples, bins, normed=True, range=range, weights=weights)

    xedges = np.delete(xedges, -1) + 0.5 * (xedges[1] - xedges[0])

    # gaussian smoothing
    if smooth:
        hist = filter.gaussian_filter(hist, sigma=0.75)
        if interpolate:
            f = interp.interp1d(xedges, hist, kind='cubic')
            xedges = np.linspace(xedges.min(), xedges.max(), 10000)
            hist = f(xedges)

    # make plot
    if label is not None:
        ax.plot(xedges, hist, color=color, lw=1.5, label=label)
    else:
        ax.plot(xedges, hist, color=color, lw=1.5)


def makesubplot2d(ax, samples1, samples2, color=True, weights=None, smooth=True,
                  bins=[40, 40], contours=True, x_range=None, y_range=None,
                  logx=False, logy=False, logz=False):

    if x_range is None:
        xmin = np.min(samples1)
        xmax = np.max(samples1)
    else:
        xmin = x_range[0]
        xmax = x_range[1]

    if y_range is None:
        ymin = np.min(samples2)
        ymax = np.max(samples2)
    else:
        ymin = y_range[0]
        ymax = y_range[1]

    if logx:
        bins[0] = np.logspace(np.log10(xmin), np.log10(xmax), bins[0])

    if logy:
        bins[1] = np.logspace(np.log10(ymin), np.log10(ymax), bins[1])

    hist2d, xedges, yedges = np.histogram2d(samples1, samples2,
                                            weights=weights,
                                            bins=bins,
                                            range=[[xmin, xmax], [ymin, ymax]])
    extent = [xedges[0], xedges[-1], yedges[0], yedges[-1]]

    if logz:
        for ii in range(hist2d.shape[0]):
            for jj in range(hist2d.shape[1]):
                if hist2d[ii, jj] <= 0:
                    hist2d[ii, jj] = 1

    xedges = np.delete(xedges, -1) + 0.5 * (xedges[1] - xedges[0])
    yedges = np.delete(yedges, -1) + 0.5 * (yedges[1] - yedges[0])

    # gaussian smoothing
    if smooth:
        hist2d = filter.gaussian_filter(hist2d, sigma=0.75)

    if contours:

        level1, level2, level3 = getsigmalevels(hist2d)
        contourlevels = (level1, level2, level3)

        #contourcolors = ('darkblue', 'darkblue', 'darkblue')
        contourcolors = ('black', 'black', 'black')
        contourlinestyles = ('-', '--', ':')
        contourlinewidths = (1.5, 1.5, 1.5)
        contourlabels = [r'1 $\sigma$', r'2 $\sigma$', r'3 $\sigma$']

        contlabels = (contourlabels[0], contourlabels[1], contourlabels[2])

        c1 = ax.contour(xedges, yedges, hist2d.T, contourlevels,
                        colors=contourcolors, linestyles=contourlinestyles,
                        linewidths=contourlinewidths, zorder=2)
    if color:
        if logz:
            c2 = ax.imshow(np.flipud(hist2d.T), extent=extent, aspect=ax.get_aspect(), interpolation='gaussian', norm=matplotlib.colors.LogNorm())
        else:
            c2 = ax.imshow(hist2d.T, extent=extent, aspect=ax.get_aspect(), interpolation='gaussian', alpha=0.8)

    if logx:
        ax.set_xscale('log')
    if logy:
        ax.set_yscale('log')


def makeChainTriPlot(covar=False):
    """
    take chain of each parameter and convert into N length array
    Combing all parameters chains into an array (N,dim); dim = # of params
    The ouput is then use for the input to plot triplot

    Parameters
    ----------
    covar: Boolean
        If true, print out the covariance matrix between T and beta, Don't do True if either of the parameter is fixed

    """
    _T = res.parameter_chain('T').flatten()
    _B = res.parameter_chain('beta').flatten()
    _alpha = res.parameter_chain('alpha').flatten()
    _f = res.parameter_chain('fnorm').flatten()
    if res._opthin != True:
        _l = res.parameter_chain('lambda0').flatten()
        if res._noalpha != True:
            stack = np.vstack((_T, _B, _alpha, _f, _l))
            mark_BestTriPlot = [res._best_fit[0][0], res._best_fit[0][1], res._best_fit[0][3], res._best_fit[0][-1], res._best_fit[0][2]]
        else:
            stack = np.vstack((_T, _B, _f, _l))
            mark_BestTriPlot = [res._best_fit[0][0], res._best_fit[0][1], res._best_fit[0][-1], res._best_fit[0][2]]
    else:
        if res._noalpha != True:
            stack = np.vstack((_T, _B, _alpha, _f))
            mark_BestTriPlot = [res._best_fit[0][0], res._best_fit[0][1], res._best_fit[0][3], res._best_fit[0][-1]]
        else:
            stack = np.vstack((_T, _B, _f))
            mark_BestTriPlot = [res._best_fit[0][0], res._best_fit[0][1], res._best_fit[0][-1]]

    if covar:
        print"Correlation Matrix between T, Beta: \n{0} \n".format(np.corrcoef(_T, _B))

    return stack.T, mark_BestTriPlot


def triplot(chain, color=True, weights=None, interpolate=False, smooth=True,
            labels=None, figsize=(11, 8.5), title=None, inj=None, saveFig=False):
    """
    Make Triangle plot of marginalized posterior distribution
    will dependent on
    - makesubplot2d()
    - makesubplot1d()
    - getsigmalevels()

    """
    from matplotlib.ticker import FormatStrFormatter, NullFormatter, NullLocator, MaxNLocator

    # rcParams settings
#    plt.rcParams['ytick.labelsize'] = 10.0
#    plt.rcParams['xtick.labelsize'] = 10.0
#    plt.rcParams['text.usetex'] = True
    plt.rcParams['figure.figsize'] = figsize

    # get number of parameters
    ndim = chain.shape[1]
    parameters = np.linspace(0, ndim - 1, ndim)

    f, axarr = plt.subplots(
        nrows=len(parameters), ncols=len(parameters), figsize=figsize)

    for i in range(len(parameters)):
        # for j in len(parameters[np.where(i <= parameters)]:

        for j in range(len(parameters)):

            ii = i
            jj = len(parameters) - j - 1
#            print "jj %d: " %(jj)

            if (jj == len(parameters) - 1) and (ii == jj):
#                 print ii
#                 print jj
                # the following line added to make the xtick labels of lambda_0 (1+z) to show clearly
                xmajorLocator = MaxNLocator(nbins=3, prune='both')
            else:
                xmajorLocator = MaxNLocator(nbins=4, prune='both')
            ymajorLocator = MaxNLocator(nbins=4, prune='both')

            if j <= len(parameters) - i - 1:
                axarr[jj][ii].xaxis.set_minor_locator(NullLocator())
                axarr[jj][ii].yaxis.set_minor_locator(NullLocator())
                axarr[jj][ii].xaxis.set_major_locator(NullLocator())
                axarr[jj][ii].yaxis.set_major_locator(NullLocator())

                axarr[jj][ii].xaxis.set_minor_formatter(NullFormatter())
                axarr[jj][ii].yaxis.set_minor_formatter(NullFormatter())
                axarr[jj][ii].xaxis.set_major_formatter(NullFormatter())
                axarr[jj][ii].yaxis.set_major_formatter(NullFormatter())
                xmajorFormatter = FormatStrFormatter('%g')
                ymajorFormatter = FormatStrFormatter('%g')

                if ii == jj:
                    # Make a 1D plot
                    makesubplot1d(axarr[ii][ii], chain[:, parameters[ii]],
                                  weights=weights, interpolate=interpolate,
                                  smooth=smooth)
                    axarr[ii][jj].set_ylim(ymin=0)

                    if inj is not None:
                        axarr[ii][ii].axvline(inj[ii], lw=2, color='k')
                else:
                    # Make a 2D plot
                    makesubplot2d(axarr[jj][ii], chain[:, parameters[ii]],
                                  chain[
                                      :, parameters[jj]], color=color, weights=weights,
                                  smooth=smooth)

                    if inj is not None:
                        axarr[jj][ii].plot(inj[ii], inj[jj], 'x', color='k', markersize=12,
                                           mew=2, mec='k')

                axarr[jj][ii].xaxis.set_major_locator(xmajorLocator)
                axarr[jj][ii].yaxis.set_major_locator(ymajorLocator)
            else:
                axarr[jj][ii].set_visible(False)
                # axarr[jj][ii].axis('off')

            if jj == len(parameters) - 1:
                axarr[jj][ii].xaxis.set_major_formatter(xmajorFormatter)
                if labels:
                    axarr[jj][ii].set_xlabel(labels[ii], labelpad=10)

            if ii == 0:
                if jj == 0:
                    axarr[jj][ii].yaxis.set_major_locator(NullLocator())
                    # axarr[jj][ii].set_ylabel('Post.')
                else:
                    axarr[jj][ii].yaxis.set_major_formatter(ymajorFormatter)
                    if labels:
                        axarr[jj][ii].set_ylabel(labels[jj], labelpad=10)

    # overall plot title
    if title:
        f.suptitle(title, fontsize=22, y=0.90)

    # make plots closer together
    f.subplots_adjust(hspace=0.1)
    f.subplots_adjust(wspace=0.1)
    if saveFig == True:
        outfile = 'CorrelationPlot__'
#        savefig('../Figures/' + outfile + filename.replace('.h5', '.png'), dvi=600)
        savefig('../Figures/' + outfile + filename.replace('.h5', '.pdf'), dvi=600)
    else:
        plt.show()


def PanelPlot(Xparam='beta', Yparam='T/(1+z)', maxLev=False, confid=False, saveFig=False):
    """
    Take MCMC samples parameters
    and make plots with sides panels = projected distribution of params

    Main panel as Contours of confidence levels and Make 2D histogram using colormesh and hist2d

    For single points , see MCMC_scatterDots()

    Input:
    ------
    Xparam: beta    (default, can be others)
    Yparam: T/(1+z)     the obs frame Temp (default, can be others)
        keyword:
        - maxLev = using percentage of the maximum value instead of confidence interval
        - confid = using percentile 68.3, 95.4, 99.7
        - else: use sigma levels

    """
    from matplotlib.gridspec import GridSpec
    X = res.parameter_chain(Xparam)
    Y = res.parameter_chain(Yparam)
    if Yparam.find('z') != -1:
        z = res.redshift
        Y *= (1 + z)
    Xrange = [min(X), max(X)]
    Yrange = [min(Y), max(Y)]

    fig = plt.figure(figsize=(7, 7))
    fig.subplots_adjust(
        hspace=0.001, wspace=0.001, left=0.10, bottom=0.095, right=0.98)
    gs = GridSpec(2, 2, width_ratios=[1, 4], height_ratios=[4, 1])

    plt.subplot(gs[1])
    Bins = 20
    hist2D, xedges, yedges = np.histogram2d(
        X, Y, bins=[Bins, Bins], range=[Xrange, Yrange], normed=False)
    # x in histogram2d on the abscissa, y on the ordinate axis
    hist2D = np.transpose(hist2D)
    # Use bin edges to restore extent
    extent = [xedges.min(), xedges.max(), yedges.min(), yedges.max()]

    # color hist2d
    plt.pcolormesh(xedges, yedges, hist2D, cmap=plt.cm.rainbow)
    # samples
    # plt.plot(X, Y, 'bo', markersize=3, alpha=.05)

    # Oplot with error contours
    if maxLev == True:
        # Use a percentage of the maximum
        maximum = np.max(hist2D)
        [L1, L2, L3] = [0.5 * maximum, 0.25 * maximum, 0.125 * maximum]
        fmtdict = {L1: '0.5Max', L2: '0.25Max', L3: '0.125Max'}
    elif confid == True:
        # using confidence interval as levels
        conf = np.percentile(hist2D, 68.3)
        conf2 = np.percentile(hist2D, 95.4)
        conf3 = np.percentile(hist2D, 99.7)
        [L1, L2, L3] = [conf, conf2, conf3]
        fmtdict = {L1: r'$68.3$', L2: r'$95.4 $', L3: r'$99.7$'}
    else:
        # use sigma
        sig1, sig2, sig3 = getsigmalevels(hist2D)
        [L1, L2, L3] = [sig1, sig2, sig3]
        fmtdict = {L1: r'$\sigma$', L2: r'$2\ \sigma $', L3: r'$3\ \sigma$'}

    cs = plt.contour(hist2D, extent=extent, levels=[L1, L2, L3], linestyles=[
        '--', '--', '--'], colors=['orange', 'orange', 'orange'], linewidths=3)

    # colour label for contour levels
    plt.clabel(cs, fmt=fmtdict, inline=True, fontsize=20)
    plt.xlim(Xrange)
    plt.ylim(Yrange)

    # Add side panels showing the projected distributions of X and Y:
    # Bin X, Y separately. In 1D bin, can use more bins now...
    S = 120
    LX = np.histogram(X, bins=S, range=Xrange, normed=True)[0]
    LY = np.histogram(Y, bins=S, range=Yrange, normed=True)[0]

    # restore positions lost by binning, similiar to making arrays in IDL
    X = Xrange[0] + (Xrange[1] - Xrange[0]) * \
        np.array(range(0, len(LX))) / float(len(LX) - 1)
    Y = Yrange[0] + (Yrange[1] - Yrange[0]) * \
        np.array(range(0, len(LY))) / float(len(LY) - 1)

    # bottom right panel: projected density of x.
    plt.subplot(gs[3])
    plt.plot(X, LX, '-', lw=3, color='black')
    plt.xticks(fontsize=16)
    plt.yticks([])
    plt.xlabel(Xparam, fontsize=24)
    plt.ylabel(r'$\cal L$', fontsize=24)
    plt.xlim(Xrange)
    plt.ylim(0.0, 1.1 * np.max(LX))

    # Top left panel: Projected density of Y
    plt.subplot(gs[0])
    # Likelihood on horizontal axis
    plt.plot(LY, Y, '-', lw=3, color='k')

    plt.yticks(fontsize=16)
    plt.xticks([])
    plt.xlabel(r'$\cal L$', fontsize=24)
    plt.ylabel(Yparam, fontsize=24)
    plt.xlim(0.0, 1.1 * np.max(LY))
    plt.ylim(Yrange)
    if saveFig == True:
        outfile = 'mcmcSamplesPanels' + Xparam + 'T' + '__'
        savefig('../Figures/' + outfile + filename.replace('.h5', '.png'))
    else:
        plt.show()


def BestFit(verbose=True):
    """
    Grab best fit parameters, e.g. chi_square, best_set of param values
    """
    BestFitIndex = res._best_fit[2]
    best_fit_chisq = -2.0 * res._best_fit[1]
    if res._opthin:
        res._best_fit[0][2] = None
    if verbose == True:
        print "BestFitParams {:s}".format(res._best_fit[0])
        print "BestFitLogLike %.2f" % (res._best_fit[1])
        print "Best Chi Square: %.2f" % best_fit_chisq
    return best_fit_chisq


def S_100Rest(alpha=0.80):
    """
    06- 25 2015: SHOULD REMOVE THIS, no reason to use 100 GHz for deriving qIR if we have different beta than 1.5 (condon 91 or 92)
    Compute flux density at rest frame 100 GHz

    Inputs:
    -------
    z: float
        redshift
    alpha: float
        spectral index in radio

    """
    z = res.redshift
    freq_obs100_Hz = 100.e9 / (1 + z)
    from scipy.constants import c
    wavelg_obs100_um = (c / freq_obs100_Hz) * 1.e6
    S100_obs_mJy = res.best_fit_sed(wavelg_obs100_um)
    S_1dot4_obs_mJy = S100_obs_mJy * (1.4e9 / 100.e9) ** (- np.abs(alpha))
    print "flux density of 1.4Ghz at observed frame: ", S_1dot4_obs_mJy
    return S_1dot4_obs_mJy


if __name__ == '__main__':
    """
    run script.py filename True False
    sys.argv[2] determines whether or not to compute values and plots for FIR
    sys.argv[3] determines whether to save all the plots or not
    """
    if len(sys.argv) < 4:
        errmsg = "Invalid number of arguments: {0:d}\n  run script.py filename FIR_True Save_True"
        raise IndexError(errmsg.format(len(sys.argv)))

    global filename
    fir_op = True if sys.argv[2].lower() == 'true' else False

    try:
        filename = sys.argv[1]
        if not filename.endswith('.h5'):
            filename += '.h5'
    except IOError:
        print"'{0:s}' is not a valid file".format(filename)

    print "... Retriving data from {0:s} ...".format(filename)
    global res
    res = mbb_emcee.mbb_results(h5file=filename)
    print "... Done reading file ..."

    if sys.argv[3].lower() == 'true':
        save_op = True
        import os
        dirStr = '../Figures'
        if not os.path.isdir(dirStr):
            print "... Making directory {:s}...".format(dirStr)
            os.mkdir(dirStr)
        else:
            print "... Directory {:s} exists...".format(dirStr)
    else:
        save_op = False

# ------- for ALMA proposal continuum -----------
#     nu_band6 = 214.2027e9
#     nu_band4 = 160.6565e9
#     nu_band7 = 348.20e9

#    from astropy import constants as const
#    wave_6 = const.c.value / nu_band6 * 1.e6
#    wave_4 = const.c.value / nu_band4 * 1.e6
#    wave_7 = const.c.value / nu_band7 * 1.e6
#
#     print("! -- ALMA Band 4 continuum: {} mJy").format(extrapolateFlux(wave_4)[0])
#    print("! -- ALMA Band 6 continuum: {} mJy").format(extrapolateFlux(wave_6)[0])
#    print("! -- ALMA Band 7 continuum: {} mJy").format(extrapolateFlux(wave_7)[0])
# --------------------------------------------------

#    convergence(saveFig=save_op)
#    chain_iter(nburn=100, saveFig=save_op)
#    Ramdon()
#    get_Paramvalues(percentile=68.3)
#    S_100Rest(alpha=0.8)
    IR = get_computation_value(FIR=fir_op)
#    Plot_Standard(saveFig=save_op)
    PlotSED(saveFig=save_op)
#    Plot_Chain(IR, FIR=fir_op, saveFig=save_op)
#    MCMC_scatterDots(saveFig=save_op)
#    PanelPlot(maxLev=False, confid=False, saveFig=save_op)
#    chi2 = BestFit()

    chainPlot, inj_best = makeChainTriPlot()
    triplotLabel = [r'$T/(1+z)$', r'$\beta$', r'$\alpha$', r'$f_{\rm norm}$', r'$\lambda_0 (1+z)$']
    if res._opthin == True:
        triplotLabel.remove(r'$\lambda_0 (1+z)$')
    if res._noalpha == True:
        triplotLabel.remove(r'$\alpha$')
    triplot(chainPlot, title=(filename.replace('.h5', ' ') + 'Parameters'),            labels=triplotLabel, inj=inj_best, saveFig=save_op)
