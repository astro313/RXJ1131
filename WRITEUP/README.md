
# HST photometry:
- on NED
- or Claeskens+06

# IR photometry
1. Herschel:
- only SPIRE, no PACS
2. WISE 
- w1 3.4 µm:
  - F_nu = 7.027 mJy +/- 0.1424 mJy
  - ang. res = 6.1"
  - 5sigma point source sensitivity = 0.08 mJy
  - 5sigma photometry sensitivity = 0.068 mJy 
- w2 4.6 µm: 
  - F_nu = 8.872 mJy +/- 0.1634 mJy
  - ang. res = 6.4"
  - 5sigma point source sensitivity = 0.11 mJy 
  - 5sigma photometry sensitivity = 0.098 mJy 
- w3 12 µm: 
  - F_nu = 21.96 +/- 0.4247 mJy
  - ang. res = 6.5"
  - 5sigma point source sensitivity = 1 mJy
  - 5sigma photometry sensitivity = 0.86 mJy
- w4 22µm: 
  - F_nu = 55.11 +/- 1.878 mJy
  - ang.res = 12.0"
  - 5sigma point source sensitivity = 6 mJy
  - 5 sigma photometry sensitivity = 5.4 mJy
- cite: WISE: Wright+10
3. 2MASS
- there are also values on NED, but with much larger aperture
- The ones listed below ~ enclosing the Einstein Ring
- J band (1.25 µm; instrumental FWHM = 2.9"):
  - F_nu = 1.009 +/- 0.090 mJy
- H band (1.65 µm; instrumental FWMH = 2.8"): 
  - F_nu = 1.448 +/- 0.1214 mJy
- Ks band (2.17 µm; instrumental FWHM = 2.9"): 
  - F_nu = 2.064 +/- 0.1597 mJy
4. Spitzer/IRAC 
- (3.8" diameter --> 5.8" diameter):
- 3.6 (PSF FWHM = 1.7"; use 3.8" aperpho: )
  * need aperture correction, --> get flux uncertainties 
- listed below are from db, applied aperture correction 
  - 4.5: 6.241 +/- 0.00207 mJy --> 7.803 +/- 2.082e-3 mJy
  - 5.8: 9.354 +/- 0.005694 mJy --> 10.72 +/- 5.124e-3 mJy
  - 8.9: 13.56 +/- 0.004518 mJy --> 14.47 +/- 4.122e-3 mJy
5. Spitzer/MIPS
- 24: 47.18 +/- 0.02621 mJy (PSF fit with native FWHM = 6"; use this)
  + for aperture fitting (12" diameter): 46.66 +/- 0.0181 mJy
    * aperture correction not applied, need to multiply by 1.488
    * error: 1.48*MAD, MAD = median of the absolute deviation from the median data
  + since source size < FWHM --> ~ point source --> aperphot ~ PSFfit 
- 70 (nope)
- 160 (nope)
6. IRAS
- upper limits, detection limit:
- 12: 0.4 Jy
- 25: 0.5 Jy
- 60: 0.6 Jy
- 100: 1.0 Jy
- ref: http://heasarc.gsfc.nasa.gov/W3Browse/iras/iraspsc.html
- I believe is the 1-sigma level


# intensity, flux density
### VLA C Band data @ 4.88 GHz
1. Core: 8.662040e-04 +/- 1.3e-5*sqrt(86/19.2738)
2. Arc: 1.273277e-03 +/- 4.2085778291312834e-05


### PdBI data
- natural weighting beam size: 4.44 x 1.95, PA 13deg

### Not PB corrected
### see 10Nov15/README.md:
#### continuum
1. logfile: CASA2Dfit_region1.log
   * 2D fit
   * region1
   * integrated flux: 1.23 +/- 0.22 mJy
   * Peak: 754 +/- 88 uJy/beam
   * Image component size (convolved with beam) ---
       - major axis FWHM:     4.17 +/- 0.53 arcsec
       - minor axis FWHM:     3.34 +/- 0.36 arcsec
       - position angle: 166 +/- 19 deg
       * deconvolved source size (N/A)
2. imstat()
   - region1
   - **continuum: 1.130 mJy**; # of pixels = 277
   - max value = 799.269 µJy/beam

#### CO(2-1)
##### Continuum-subtracted
CO (2-1) velocity ranges used for different purpose:
- line flux Intensity: 124-156
- moment 0 map: 127-155
- RGB: 126-160
- RedBlue: 134-150

# See 14Oct15/README.md (Deprecated)
- cell size =0.5"
- I = 20.3 Jy Km/s 
- sigma_ch ~ 1.451 mJy/ Beam
- velocity resolution per bin: 21.528154 km/s
- channel range 
    + [134, 142] (red) and [142, 150] (blue)
- sigma_I = 1.451 / sqrt(32) * 32 * 21.528154 = 176.6 mJy km/s / Beam
- sigma_I (spatially integrated) = 176.6 * 2.5 = 441.5 mJy km/s << wrong, because *region* is >2.5 beams, altho source is resovled over ~2.2 beams

# updated, see 15May16/README.md
- based on moment 0 map based with channel range [124, 156]
- updated CO21 I = 24.14 +/- 0.62 Jy kms
- sigma_ch: same as above
- velocity resolution per bin: same as above
- region ~ 13 beams
- sigma_I = 1.451 * sqrt(156-124+1) * 21.5 * 13 beam = 2.33 Jy km/s/beam
- peak flux from moment 0 map = 7.23 +/- 0.46 Jy km/s /Beam
- source size, deconvolved from beam: 5.1"+/-0.72" x 3.72"+/-0.66", PA: 158+/-23
- this corresponds to 85.5 px
--> resolved over ~2.2 beams
- note that our source is not quite like a 2D gaussian, but since it recovers the intensity from region, it's probably a good enough estimate for source size

# 0th moment map (all channels):
- 14Oct15/mom0.fits
- Range By eye in go view:
  - chan: 155, 127
  - vel: -530.33, 72.73 (z=0 obs. frame)
- no clipping
- 0th moments
- 14Oct15/README.md

-0th moment RGB:
- see 12Nov15/README.md
- blue: channel [160, 146]; velocity [-651.90, -342.664]
Noise:
    from line-free:
    calculate:
- Central: channel [141, 145]; velocity [-220.335, -342.664]
Map Noise:
    from line-free:
    calculate:
- red: channel [126, 140]; velocity [75.81, -220.335]

-1st, 2nd moments (clipped @ 5sigma):
- see 26Dec15/README.md & 30Apr16/README.md


### CARMA data, combined D2, D3
#### ~1.4mm Continuum
- @ spatial position of CO emission is burried in noise - 2.2 +/- 3.0 mJy
- **upper limit: line-free channels mfs image rms = 8.30807E-04 Jy/B**

#### CO (3-2):
- natural weighting beam size: **3.19 x 1.86 PA: 8 deg**
- **channel width = 2x17.935 = 35.87 km/s**
- Theoretical rms noise: 8.834E-03 Jy/B per channel
- map off region rms: 1.31957E-02 Jy/B per channel (before clean);  1.32997E-02 Jy/B per channel (after clean)
- image cubes: /Users/admin/Research/RXJ1131/CARMA/imagingD23/
- Line intensity FWZI from Gauss Fit = sqrt(2*pi)*a*c = **37.509 Jy km/s**, see /Users/admin/Research/RXJ1131/CARMA/imagingD23/GILDAS/Gaussfit_bin2.greg
- ~ chan [49,70]

##### Intensity from summing up fluxes within velo range guided by CO21
see CARMA/specOPlot/README.md
- I = 37.1 Jy km/s 
- in 21 channels [49, 69]
- velocity [-430.94, 305.26] km/s
- --> # of beams ~ 10
- Isigma ~ 13.3 mJy/B/channel * 21 channels * 10 beams ~ 2.793 Jy km/s

