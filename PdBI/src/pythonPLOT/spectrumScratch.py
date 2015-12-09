'''
Build a spectrum from scratch (and fit a gaussian)

Can't figure out how to interactively show. For now, save it and check.

Output: spectrum with fit

'''
import numpy as np
import pyspeckit

xaxis = np.linspace(-50,150,100.)
sigma = 10.
center = 50.
synth_data = np.exp(-(xaxis-center)**2/(sigma**2 * 2.))

# Add noise
stddev = 0.1
noise = np.random.randn(xaxis.size)*stddev
error = stddev*np.ones_like(synth_data)
data = noise+synth_data

# this will give a "blank header" warning, which is fine
sp = pyspeckit.Spectrum(data=data, error=error, xarr=xaxis,
                        xarrkwargs={'unit':'km/s'},
                        unit='erg/s/cm^2/AA')

sp.plotter(figure=1)
sp.plotter(errstyle='fill')


# Fit with automatic guesses
sp.specfit(fittype='gaussian')
sp.plotter.savefig('spectrumScratch.png')

# Fit with input guesses
# The guesses initialize the fitter
# This approach uses the 0th, 1st, and 2nd moments
amplitude_guess = data.max()
center_guess = (data*xaxis).sum()/data.sum()
width_guess = data.sum() / amplitude_guess / np.sqrt(2*np.pi)
guesses = [amplitude_guess, center_guess, width_guess]
sp.specfit(fittype='gaussian', guesses=guesses)
sp.specfit.plot_fit()

sp.plotter.savefig('spectrumScratchGuess.png')