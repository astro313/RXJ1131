copied centralizedCube4GILDAS.fits and centralizedCube4GILDAS.lmv-clean from 30Apr16.

Want to generate a moment 0 map without clipping.
Previously use 127-155, where in 30Apr16, I found that using just 127-155 have excluded significant flux.
Should use 126-160 for highest SNR, intensity use 124-156.

*Need to get source size, # of beams, and uncertainties on the integrated flux.*

-----------------------
go moments
! Moment.init
TASK\FILE "Input file name" IN$ "/Users/admin/Research/RXJ1131/PdBI/data/15May16/centralizedCube4GILDAS.lmv-clean"
TASK\FILE "Output files name (no extension)" OUT$ "/Users/admin/Research/RXJ1131/PdBI/data/15May16/centralizedCube4GILDAS-mom"
TASK\REAL "Velocity range" VELOCITY$[2]  -258.3 430.6
TASK\REAL "Detection threshold" THRESHOLD$ 0
TASK\LOGICAL "Smooth before detection ?" SMOOTH$ YES
TASK\GO

fits centralizedCube4GILDAS-mom0.fits from centralizedCube4GILDAS-mom.mean 

let name centralizedCube4GILDAS-mom
let type mean
go nice

poly
mean
--> rms = 0.1818 Jy/B, too high?
go noise
--> histogram skewed..

# why histogram skewed? (see below)
let name centralizedCube4GILDAS
let type lmv-clean
let first 10
let last 110
go noise
--> histogram looks fine, rms = 1.41 mJy/B

let first 120
let last 160 
go noise
--> histogram also looks fine, rms = 1.55 mJy/B, which is not surprising since there are emssions in these channels

let name centralizedCube4GILDAS
let type lmv-clean
for i 120 to 158
    let first i
    let last i
    go nice
    poly
    mean
    go noise
    sic wait 1
next
--> histogram of each plane still center around zero, and rms ~ 1.5 mJy/B, which is still ok

let first 130
let last 130
let size 100
go nice
poly
mean
--> rms = 1.57 mJy/B, still ok
go noise
--> looks ok
Don't know why. let's assume it's ok for now.

----------------------------

vertice to get I = 20.3 Jy km/s
133,109
141,115
142,130
138,140
128,140
122,131
123,118
128,110
--> saved as CASA region file: Flux_CO32_20.3Jy_region.txt

inp importfits
importfits(fitsimage='centralizedCube4GILDAS.fits',imagename='centralizedCube4GILDAS.image')
inp imstat
imstat(imagename='centralizedCube4GILDAS.image', region='Flux_CO32_20.3Jy_region.txt', logfile='Flux_CO32_20.3Jy_region.log')
--> no.. because it cal. N point in many channels, but I just want one plane

importfits(fitsimage='centralizedCube4GILDAS-mom0.fits',imagename='centralizedCube4GILDAS-mom0.image')
imstat(imagename='centralizedCube4GILDAS-mom0.image', region='Flux_CO32_20.3Jy_region.txt', logfile='Flux_CO32_20.3Jy_region.log')
- I = 24.14 Jy km/s
- pts area within region = 504 px
- number of pix in a beam = numpy.pi * 4.43 * 1.93 / 4 / numpy.log(2) / 0.5**2 = 38.8 px
- --> # of beams ~ 13 within the region
- theoretical Isigma ~ 1.451 mJy/B/channel * sqrt(156-124+1) channels * 21.5 * 13 beams ~ 2.330 Jy km/s
- --> should be 1.451 mJy/B/channel * sqrt(156-124+1) channels * 21.5 * 13^.5 beams ~ 0.647 Jy km/s

Altho region is 13 beams, if we do imfit():
imfit(imagename='centralizedCube4GILDAS-mom0.image', region='Flux_CO32_20.3Jy_region.txt', logfile='centralizedCube4GILDAS-mom0.imfit.log')
- I = 24.6 +/- 2.0 Jy km/s
- peak S = 7.23 +/- 0.46 Jy km/s /Beam
- source size, deconvolved from beam = 5.1"+/-0.72" x 3.72"+/-0.66", PA: 158+/-23
- this corresponds to 85.5 px
--> ~2.2 beams
- note that our source is not quite like a 2D gaussian, but since it recovers the intensity from region, it's probably a good enough estimate for source size
- hence, it's resolved over ~ 2.2 beams

imstat(imagename='centralizedCube4GILDAS-mom0.image', region='offsource_region.crtf', logfile='offsource_region.log')
--> rms offsource = 0.484 Jy km/s/B; but STD = 0.173 Jy km/s/B
--> the histogram also shows the mean is off from zero.
--> note rms higher than calculated from 1.451 mJy/B/channel


mapping
(try no smoothing)
go moments
! Moment.init
TASK\FILE "Input file name" IN$ "/Users/admin/Research/RXJ1131/PdBI/data/15May16/centralizedCube4GILDAS.lmv-clean"
TASK\FILE "Output files name (no extension)" OUT$ "/Users/admin/Research/RXJ1131/PdBI/data/15May16/centralizedCube4GILDAS-mom_nosmooth"
TASK\REAL "Velocity range" VELOCITY$[2]  -258.3 430.6
TASK\REAL "Detection threshold" THRESHOLD$ 0
TASK\LOGICAL "Smooth before detection ?" SMOOTH$ NO
TASK\GO

let name centralizedCube4GILDAS-mom_nosmooth
let type mean
go nice
--> same problem


# To solve this unexpected noise properties in the moment0 map made using GILDAS
1. see script RXJ1131/PdBI/src/python_mom0/makeMom0.py, where I made the moment0 manually, with its noise properties consistent with what we expect.
--> use this for moment0 in APLpy, this also uses channel 126-160 --> highest SNR



