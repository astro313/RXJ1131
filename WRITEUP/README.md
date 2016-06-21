
# HST photometry:
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
### Not PB corrected
### see 10Nov15/README.md:
#### continuum
- averaged over chan 1 120 165 360
- number of channels = 316
- sigma = 1.451 mJy/B / sqrt(316) = 0.081625 mJy/B
- go noise on map (whole map) after cleaning: 89.274 uJy/B
- F = 1.2 mJy (integrated spatially, from ~ 1 sigma contour; 150 pixels)
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
   - continuum: 1.130 mJy, over 277 px
   - max value = 799.269 µJy/beam

#### CO(2-1)
##### Continuum-subtracted
CO (2-1) velocity ranges used for different purpose:
- line flux Intensity: 124-156
- moment 0 map: 127-155 (deprecated as of 15May16)
  + --> use chan 126-160 (highest SNR)
- RGB: 126-160
- RedBlue: 126-160 (17May16)

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
- made moment 0 map with channel range [124, 156] (deprecated because this gives highest flux, but not highest SNR, 126-160 has highest SNR)

- updated CO21 I = 24.14 +/- 0.62 Jy kms, but report 2.44 Jy km/s
  - chan: 124-156, most flux
- sigma_ch: same as above
- velocity resolution per bin: same as above
- intensity region ~ 13 beams
- sigma_I = 1.451 * sqrt(156-124+1) * 21.5 * 13 beam = 2.33 Jy km/s
- peak flux from moment 0 map = 7.23 +/- 0.46 Jy km/s /Beam (Deprecated)
- source size, deconvolved from beam: 5.1"+/-0.72" x 3.72"+/-0.66", PA: 158+/-23
- this corresponds to 85.5 px
--> resolved over ~2.2 beams
- note that our source is not quite like a 2D gaussian, but since it recovers the intensity from region, it's probably a good enough estimate for source size

20 June 2016:
- drive intrinsic line flux from the updated CO21 I = 24.14 Jy km/s
  + sum up flux in chan 124 - 156, divide by mag. from respective models
  + 124-156 <-> -258.3, 430.6 km/s (z=0 frame, after shifting velo in ASCII file; see 22May16/)
  + assume no mag. factor for channel 124 (suggested by channel map) & no mag. for 125
  + 126-130: derive Mgas separately for RXJ and companion
  + Companion: M_gas = (0.192 +/- 0.090)e10 Msun
  + delensed flux_CO(2-1) for just RXJ1131 from 124, 125, 131-156: -0.825 + 5.157 + (72.845 + 59.771 + 48.868 + 41.295 + 31.943) / 7.6 + (38.128 + 25.169 + 29.686 + 14.481 + 17.699) / 8.7 + (13.464 + 13.718 + 14.077 + 6.991 + 9.303) / 4.1 + (14.191 + 14.675 + 21.192 + 32.956 + 22.292) / 4.2   + (18.555 + 24.522 + 23.988 + 16.965 + 17.155)/4.3 + 6.809 / 3.1 = 117.0726 mJy
  + --> intrinsic I(2-1) for just RXJ1131 from chan 124, 125, 131-156 = 117.0726 * 21.5/1e3 = 2.517 Jy km/s
  + ---> total intrinsic I(2-1) for just RXJ1131 = 2.517 + 0.417 = 2.934 Jy km/s; 
  + ---> unc. on intrinsic I(2-1) on RXJ: assuming overall uncertainty with all the magnification factors ~ sum((unc. on each/mag of each)**2)/5 = 1.00, where I am using error from 126-160 and assuming they contribute to the intrinsic I calculation equally (which is not true beucase chan 124,125  dont have mag. factor) = sqrt((1.0/5.5)^2 + (1.451e-3 * 2.5 * 21.5 * sqrt(156-124+1)/2.934)^2) * 2.934 = 0.697 Jy km/s
  + ---> RXJ M_gas = (1.375 +/- 0.327) e10 Msun 
- for diff. lensing section: derive intrinsic line flux for red and blue wing respectively. 
  + spectrum for corrected for lensing: de-lensed.spec in 22May16/
  + --> spectrum look symmetric
  + --> fit double Gauss to refine channel ranges for the red and blue-shifted kinematic components
    * --> peak flux densities are consistent with being same 
    * ----> symmetric disk
  + ---> red: channel [124-142]; blue: [143-156]
  + Ired = 21.5 * (-0.825 + 5.157 + (72.845 + 59.771 + 48.868 + 41.295 + 31.943) / 7.6 + (38.128 + 25.169 + 29.686 + 14.481 + 17.699) / 8.7 + (13.464 + 13.718) / 4.1) + 0.417 = 1266 +/- (sqrt((1.0/5.5)^2 + (1.451e-3 * 2.5 * 21.5 * sqrt(142-124+1)/1266)^2) * 1266 = 230.2) mJy km/s
  + Iblue = 21.5 * ((14.077 + 6.991 + 9.303) /4.1 + (14.191 + 14.675 + 21.192 + 32.956 + 22.292) / 4.2 + (18.555 + 24.522 + 23.988 + 16.965 + 17.155)/4.3 + 6.809 / 3.1) = 1251.5 +/- (sqrt((1.0/5.5)^2 + (1.451e-3 * 2.5 * 21.5 * sqrt(156-143+1)/1251.5)^2) * 1251.5 = 227.5) mJy km/s
  + --> approximately a symmetric disk

# 0th moment map (all channels):
- 15May16/centralizedCube4GILDAS-CASA_ch126-160_mom0.fits
- highest SNR: 
  + CASA chan=125~159 <=> GILDAS 126-160 <=> python 125-160
- no clipping
- sigma = 0.305 Jy km/s/Beam
  - note that theoretical sigma is lower = 1.5 * sqrt(160-126+1) * 21.5 ~ 0.2 Jy km/s /beam
  - due to higher noise in some channels with emission

- 14Oct15/mom0.fits (Deprecated because noise properties is messed up)
- Range By eye in go view:
  - chan: 155, 127
  - vel: -530.33, 72.73 (z=0 obs. frame)
- no clipping
- 14Oct15/README.md

-0th moment RGB: (Deprecated)
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

# 1st, 2nd moments (clipped):
- see 26Dec15/README.md & 30Apr16/README.md


### CARMA data, combined D2, D3
#### ~1.4mm Continuum
- spatial position of CO emission is burried in noise - 2.2 +/- 3.0 mJy
- consistent with the Gaussian + poly fit to the spectrum where continuum: 2.2 $\pm$ 3.0 mJy.
- line free channels:1, 118, 160, 600
- nchan =(117+600-160+1) = 558
- sigma =  0.83 mJy/B
- **upper limit: line-free channels mfs image rms = 8.30807E-04 Jy/B**

#### CO (3-2):
- natural weighting beam size: **3.19 x 1.86 PA: 8 deg**
- **channel width = 2x17.935 = 35.87 km/s**
- Theoretical rms noise: 8.834E-03 Jy/B per channel
- map off region rms: 13.1957 mJy/B per channel (before clean);  1.32997E-02 Jy/B per channel (after clean)
- image cubes: /Users/admin/Research/RXJ1131/CARMA/imagingD23/
- Line intensity FWZI from Gauss Fit = sqrt(2*pi)*a*c = **37.509 Jy km/s**, see /Users/admin/Research/RXJ1131/CARMA/imagingD23/GILDAS/Gaussfit_bin2.greg
- ~ chan [49,70]

##### Intensity from summing up fluxes within velo range guided by CO21
- I = 35.7 Jy km/s from `go view`
- Update channel range to be consistent with the velocity range used to get the new CO21 Intensity (124-156)
- --> integrated over velocity width ~ 710 km/s
- Channel:     48-68
- velocity:   -441.57, 278.79
- I_sigma: see below

see CARMA/specOPlot/README.md (Deprecated)
- I = 37.1 Jy km/s from `go view`
- in 21 channels [49, 69] 
- velocity [-430.94, 305.26] km/s
- vertices of polygon (Xi, Yi) = (133,150), (118,138), (122,104), (153, 114), (150,141)
- CASA region file of the above polygon: Flux_CO32_37.2Jy_region.txt
- pts area within region = 1159 px
- # of pix in a beam = 112.867 PX
- --> # of beams ~ 10
- Isigma ~ 13.3 mJy/B/channel * np.sqrt(21) channels * 35.87 km/s * 10 beams ~ 21.8621 Jy km/s

