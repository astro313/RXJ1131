'''
Last Modified: 16 Aug 2016


History:
16Aug16:
  - copied from J0939 VLA CO1-0
  - modified up to making .ms for all SBs, USBs, and LSBs; then making clean & dirty images of the full SB cont

Note:
- only ran up to making a clean image using all SBs & LSB
- in the future, if we want to generate separate images for the USB, will have to modified the remaning code that are currently commented out

Make dirty and clean continuum images

1. LSB+USB
2. USB     spw: 0,1,4,5
3. LSB     spw: 2,3,6,7    <--> PdBI
'''

import os

scriptmode = False

# The prefix to use for all output files
# default nterms
prefix = '/Users/admin/Research/RXJ1131/ALMA/C2/member.uid___A001_X13b_Xf3/Daisy_PL_inspectCont/RXJ1131'

splitms = prefix + '.src.split.ms'
outputvis = splitms


# split out science
default('split')

vis = "/Users/admin/Research/RXJ1131/ALMA/C2/member.uid___A001_X13b_Xf3/Daisy_PL_inspectCont/calibrated.ms"
splitms = prefix + '.src.split.ms'
outputvis = splitms
field = '5'
spw = ''
datacolumn = 'data'

saveinputs('split', prefix + '.split.RXJ1131.saved')

split()

#=====================================================================
# Create an Averaged Continuum MS by averaging line-free channels
# in both basebands
#
# use this for continuum imaging, not the one from uvcontsub
#

casalog.post('-- Create an Averaged Continuum MS of Both Basebands--', 'INFO')
print '-- Create an Averaged Continuum MS of Both Basebands--'


# Use plotms to identify line and continuum spw
plotms(vis=outputvis, xaxis='channel', yaxis='amp', ydatacolumn='data', avgtime='1e8s', avgscan=True, avgchannel='2', iteraxis='spw')


plotms(vis=outputvis, xaxis='channel', yaxis='amp', ydatacolumn='data', avgtime='1e8s', avgscan=True, spw='2:100~300,6:100~300')
# chan ~170-210 is line

plotms(vis=outputvis,xaxis='velocity', yaxis='amp', avgtime='1e7', avgscan=T, restfreq='139.40738949023GHz', xselfscale=T, freqframe='LSRK', transform=T, correlation='XX,YY', spw='2,6', avgbaseline=True)

plotms(vis=outputvis, xaxis="frequency", yaxis="amp",
       avgtime='1e7', avgscan=True, avgbaseline=True)

# "line channels"
linechans = '2:170~210,6:170~210'

# Set spws to be used to form continuum
contspws = '0,1,2,3,4,5,6,7'

# Need to flag the line channels prior to averaging.
flagmanager(vis=outputvis, mode='save', versionname='before_line_flagging')

flagdata(vis=outputvis, mode='manual', spw=linechans, flagbackup=False)


# check that flags are as expected, NOTE must check reload on plotms
# gui if its still open.
plotms(vis=outputvis, yaxis='amp', xaxis='channel',
       avgchannel='5', avgtime='1e8', avgscan=True, spw='2,6')

# Average the channels within spws
contvis = '/Users/admin/Research/RXJ1131/ALMA/C2/member.uid___A001_X13b_Xf3/Daisy_PL_inspectCont/RXJ1131_cont.ms'
rmtables(contvis)

listobs(outputvis)
# spw native width?

default('split')
outputvis=splitms
split(vis=outputvis,
      spw=contspws,
      outputvis=contvis,
      width=[128,128,480,128,128,128,480,128],
      datacolumn='data')           # s.t. not affected by uvcontsub

# Inspect continuum for any problems
plotms(vis=contvis, xaxis='uvdist', yaxis='amp', coloraxis='spw')


#=====================================================================
# Create an Averaged Continuum MS by averaging line-free channels
# in low frequency baseband
#
print '-- Create an Averaged Continuum MS of Lower Freq Baseband--'

# Set spws to be used to form continuum
contspwsL = '2~3,6~7'

# Average the channels within spws
contvisL = '/Users/admin/Research/RXJ1131/ALMA/C2/member.uid___A001_X13b_Xf3/Daisy_PL_inspectCont/RXJ1131_contL.ms'

rmtables(contvisL)

default('split')
outputvis = splitms
split(vis=outputvis,
      spw=contspwsL,
      outputvis=contvisL,
      width=[480, 128, 480, 128],
      datacolumn='data')

# Inspect continuum for any problems
plotms(vis=contvisL, xaxis='uvdist', yaxis='amp', coloraxis='spw')


#=====================================================================
# Create an Averaged Continuum MS by averaging line-free channels
# in high frequency baseband
#
print '-- Create an Averaged Continuum MS of High Freq Baseband--'

# Set spws to be used to form continuum
contspwsU = '0~1,4~5'

# Average the channels within spws
contvisU = '/Users/admin/Research/RXJ1131/ALMA/C2/member.uid___A001_X13b_Xf3/Daisy_PL_inspectCont/RXJ1131_contU.ms'
rmtables(contvisU)

default('split')
outputvis = splitms
split(vis=outputvis,
      spw=contspwsU,
      outputvis=contvisU,
      width=[128,128,128,128],
      datacolumn='data')

# Inspect continuum for any problems
plotms(vis=contvisU, xaxis='uvdist', yaxis='amp', coloraxis='spw')
#
#
#
# restore the previous flagged line channels
flagmanager(vis=outputvis, mode='restore',
            versionname='before_line_flagging')
#
#
#=====================================================================
# Make dirty continuum map of RXJ1131
#

print '-- Clean (make dirty continuum image) --'

contimagename = contvis[:contvis.find('.ms')]
imname = contimagename + '.dirty'
for ext in ['.flux', '.image', '.mask', '.model', '.pbcor', '.psf', '.residual', '.flux.pbcoverage']:
    rmtables(imname + ext)

default('clean')

vis = contvis
imagename = imname

mode = 'mfs'
psfmode = 'clark'
imsize = [640]
cell = ['0.07arcsec']
niter = 0
threshold = 0
stokes = 'I'
interactive = False

# weighting = 'briggs'
# robust = 0.5

saveinputs('clean', prefix + '.cont.invert.saved')

# Pause script if you are running in scriptmode
if scriptmode:
    inp()
    user_check = raw_input('Return to continue script\n')

clean()

dirtyimage = imname + '.image'
viewer(dirtyimage)
# get offline rms
rms = 1.55e-5

#=====================================================================
#
# Make and clean continuum map of RXJ1131
#
print '-- Clean (clean) Continuum image --'

imname = contimagename + '.clean'
for ext in ['.flux', '.image', '.mask', '.model', '.pbcor', '.psf', '.residual', '.flux.pbcoverage']:
    rmtables(imname + ext)

default('clean')

vis = contvis
imagename = imname

mode = 'mfs'
psfmode = 'clark'
imsize = [640]
cell = ['0.07arcsec']
niter = 10000
threshold = rms
stokes = 'I'
interactive = True
threshold = threshold

# Set up the weighting
# Use Briggs weighting (a moderate value, on the uniform side)
# weighting = 'briggs'
# robust = 0.5

saveinputs('clean', prefix + '.cont.clean.saved')

# Pause script if you are running in scriptmode
if scriptmode:
    inp()
    user_check = raw_input('Return to continue script\n')


clean()

clnimage = imname + '.image'
#=====================================================================
#
# Done with imaging
# Now view the image cube of RXJ1131
#
if scriptmode:
    print '--View image--'
    viewer(clnimage)
    user_check = raw_input('Return to continue script\n')

viewer(clnimage)
# rms = ?       # Jy /beam

# --> viewer olay Our PdBI continuum

# #
# # Alternatively, you can use the scripting "imview" approach.
# #
# imview(raster={'file': clnimage,
#                'range': [-1.5e-4, 1.5e-2],
#                'colormap': 'Rainbow 2', 'scaling': 0.0, 'colorwedge': True},
#        contour={'file': clnimage,
#                 'levels': [-6, -3, 3, 4, 5, 6, 7, 8, 9, 12, 15, 16],
#                 'unit': rms},       # Jy/beam
#        zoom=3)


# #=====================================================================
# #
# # export the Final CLEAN Image as FITS
# # Run asynchronously so as not to interfere with other tasks
# #
# print '--Final Export CLEAN FITS--'
# default('exportfits')
# #
# clnfits = prefix + '.cont.clean.fits'
# #
# imagename = clnimage
# fitsimage = clnfits
# async = True
# #
# saveinputs('exportfits', prefix+'cont.clean.exportfits.saved')
# #
# exportfits()

# #=====================================================================
# #
# # Print the image header
# #
# print '--Imhead--'
# default('imhead')

# imagename = clnimage

# mode = 'summary'

# imhead()

# # A summary of the map will be seen in the logger

# #=====================================================================
# #
# # Get map statistics
# #
# print '--Imstat --'
# default('imstat')

# imagename = clnimage

# box_on = '110,110,140,140'
# box_off = '30,30,218,100'

# box = box_off
# imgstat = imstat()               # dictionary; imgstat['key'][axis]
# rms = (imgstat['rms'][0])
# print '>> rms: '+str(rms)

# box = box_on
# imgstat_on = imstat()
# peak = (imgstat_on['max'][0])
# print '>> Peak: '+str(peak)

# print '>> Dynamic range: '+str(peak/rms)


# ------------ Repeat for basebnad continuum image --------------
#=====================================================================
# Make dirty continuum map of RXJ1131
#

print '-- Clean (make dirty continuum image) --'

contimagename = contvisL[:contvisL.find('.ms')]
imname = contimagename + '.dirty'
for ext in ['.flux', '.image', '.mask', '.model', '.pbcor', '.psf', '.residual', '.flux.pbcoverage']:
    rmtables(imname + ext)

default('clean')

vis = contvisL
imagename = imname

mode = 'mfs'
psfmode = 'clark'
imsize = [640]
cell = ['0.07arcsec']
niter = 0
threshold = 0
stokes = 'I'
interactive = False

weighting = 'briggs'
robust = 0.5

saveinputs('clean', prefix + '.contL.invert.saved')

# Pause script if you are running in scriptmode
if scriptmode:
    inp()
    user_check = raw_input('Return to continue script\n')

clean()

dirtyimage = imname + '.image'
viewer(dirtyimage)
# get offline rms
rms = 1.9e-5

#=====================================================================
#
# Make and clean continuum map of RXJ1131
#
print '-- Clean (clean) Continuum image --'

imname = contimagename + '.clean'
for ext in ['.flux', '.image', '.mask', '.model', '.pbcor', '.psf', '.residual', '.flux.pbcoverage']:
    rmtables(imname + ext)

default('clean')

vis = contvisL
imagename = imname

mode = 'mfs'
psfmode = 'clark'
imsize = [640]
cell = ['0.07arcsec']
niter = 10000
threshold = rms
stokes = 'I'
interactive = True
threshold = threshold

# Set up the weighting
# Use Briggs weighting (a moderate value, on the uniform side)
# weighting = 'briggs'
# robust = 0.5

saveinputs('clean', prefix + '.contL.clean.saved')

# Pause script if you are running in scriptmode
if scriptmode:
    inp()
    user_check = raw_input('Return to continue script\n')


clean()

clnimage = imname + '.image'


# convolve w/ the beam size of our PdBI obs.
BMJ = 4.44
BMN = 1.95
BPAa = 13.0

# smooth images
inputimage = clnimage
outputimage = inputimage + '.smooth'
imsmooth(imagename=inputimage,
         outfile=outputimage,
         kernel='gauss',
         major=str(BMJ)+'arcsec',
         minor=str(BMN)+'arcsec',
         pa=str(BPAa)+'deg',
         targetres=True
        )

viewer(outputimage)


# #=====================================================================
# #
# # Done with imaging
# # Now view the image cube of RXJ1131
# #
# if scriptmode:
#     print '--View image--'
#     viewer(clnimage)
#     user_check = raw_input('Return to continue script\n')

# viewer(clnimage)
# rms = 4.5e-5
# #
# # Alternatively, you can use the scripting "imview" approach.
# #
# imview(raster={'file': clnimage,
#                'range': [-0.001, 1e-2],
#                'colormap': 'Rainbow 2', 'scaling': 0.0, 'colorwedge': True},
#        contour={'file': clnimage,
#                 'levels': [-6, -3, 3, 4, 5, 6, 7, 8, 9, 12, 15, 16],
#                 'unit': rms},       # Jy/beam
#        zoom=3)


# #=====================================================================
# #
# # export the Final CLEAN Image as FITS
# # Run asynchronously so as not to interfere with other tasks
# #
# print '--Final Export CLEAN FITS--'
# default('exportfits')
# #
# clnfits = prefix + '.contL.clean.fits'
# #
# imagename = clnimage
# fitsimage = clnfits
# async = True
# #
# saveinputs('exportfits', prefix+'.contL.clean.exportfits.saved')
# #
# exportfits()

# #=====================================================================
# #
# # Print the image header
# #
# print '--Imhead--'
# default('imhead')

# imagename = clnimage

# mode = 'summary'

# imhead()

# # A summary of the cube will be seen in the logger

# #=====================================================================
# #
# # Get map statistics
# #
# print '--Imstat --'
# default('imstat')

# imagename = clnimage

# box_on = '110,110,140,140'
# box_off = '30,30,218,100'

# box = box_off
# imgstat = imstat()               # dictionary; imgstat['key'][axis]
# rms = (imgstat['rms'][0])
# print '>> rms: '+str(rms)

# box = box_on
# imgstat_on = imstat()
# peak = (imgstat_on['max'][0])
# print '>> Peak: '+str(peak)

# print '>> Dynamic range: '+str(peak/rms)

# #=====================================================================
# # Make dirty continuum map of RXJ1131, high frequency baseband
# #

# print '-- Clean (make dirty continuum image) --'

# contimagename = prefix.replace('RXJ1131', 'RXJ1131_contU')
# imname = contimagename + '.dirty'
# for ext in ['.flux', '.image', '.mask', '.model', '.pbcor', '.psf', '.residual', '.flux.pbcoverage']:
#     rmtables(imname + ext)

# default('clean')

# vis = contvisU
# imagename = imname

# mode = 'mfs'
# psfmode = 'clark'
# imsize = [256]
# cell = ['0.75arcsec']
# niter = 0
# threshold = 0
# stokes = 'I'
# interactive = False

# # weighting = 'briggs'
# # robust = 0.5

# saveinputs('clean', prefix + '.contU.invert.saved')

# # Pause script if you are running in scriptmode
# if scriptmode:
#     inp()
#     user_check = raw_input('Return to continue script\n')

# clean()

# dirtyimage = imname + '.image'
# viewer(dirtyimage)
# # get offline rms
# rms = 0.24        # mJy/beam

# #=====================================================================
# #
# # Make and clean continuum map of RXJ1131
# #
# print '-- Clean (clean) Continuum image --'

# contimagename = prefix.replace('RXJ1131', 'RXJ1131_contU')
# imname = contimagename + '.clean'
# for ext in ['.flux', '.image', '.mask', '.model', '.pbcor', '.psf', '.residual', '.flux.pbcoverage']:
#     rmtables(imname + ext)

# default('clean')

# vis = contvisU
# imagename = imname

# mode = 'mfs'
# psfmode = 'clark'
# imsize = [256]
# cell = ['0.75arcsec']
# niter = 10000
# threshold = rms
# stokes = 'I'
# interactive = True
# threshold = threshold

# # Set up the weighting
# # Use Briggs weighting (a moderate value, on the uniform side)
# # weighting = 'briggs'
# # robust = 0.5

# saveinputs('clean', prefix + '.contU.clean.saved')

# # Pause script if you are running in scriptmode
# if scriptmode:
#     inp()
#     user_check = raw_input('Return to continue script\n')


# clean()

# clnimage = imname + '.image'
# #=====================================================================
# #
# # Done with imaging
# # Now view the image cube of RXJ1131
# #
# if scriptmode:
#     print '--View image--'
#     viewer(clnimage)
#     user_check = raw_input('Return to continue script\n')

# viewer(clnimage)
# #
# # Alternatively, you can use the scripting "imview" approach.
# #
# imview(raster={'file': clnimage,
#                'range': [-0.001, 3e-3],
#                'colormap': 'Rainbow 2', 'scaling': 0.0, 'colorwedge': True},
#        contour={'file': clnimage,
#                 'levels': [-6, -3, 3, 4, 5, 6, 7, 8, 9, 12, 15, 16],
#                 'unit': threshold / 1e3},       # Jy/beam
#        zoom=3)


# #=====================================================================
# #
# # export the Final CLEAN Image as FITS
# # Run asynchronously so as not to interfere with other tasks
# #
# print '--Final Export CLEAN FITS--'
# default('exportfits')
# #
# clnfits = prefix + '.contU.clean.fits'
# #
# imagename = clnimage
# fitsimage = clnfits
# async = True
# #
# saveinputs('exportfits', prefix+'contU.clean.exportfits.saved')
# #
# myhandle2 = exportfits()

# #=====================================================================
# #
# # Print the image header
# #
# print '--Imhead--'
# default('imhead')

# imagename = clnimage

# mode = 'summary'

# imhead()

# # A summary of the cube will be seen in the logger

# #=====================================================================
# #
# # Get map statistics
# #
# print '--Imstat --'
# default('imstat')

# imagename = clnimage

# box_on = '110,110,140,140'
# box_off = '30,30,218,100'

# box = box_off
# imgstat = imstat()               # dictionary; imgstat['key'][axis]
# rms = (imgstat['rms'][0])
# print '>> rms: '+str(rms)

# box = box_on
# imgstat_on = imstat()
# peak = (imgstat_on['max'][0])
# print '>> Peak: '+str(peak)

# print '>> Dynamic range: '+str(peak/rms)


# Maybe we didn't exclude all the line channels when making LSB cont., maybe we just excluded the red wing CO..?
# if we just make the LSB cont. using the line-free channels, is the background still brighter than the fg cont?

dirtyimage = '/Users/admin/Research/RXJ1131/ALMA/C2/member.uid___A001_X13b_Xf3/Daisy_PL_inspectCont/RXJ1131_cont_spw3_7.dirty'

clean(vis=splitms,
      imagename=dirtyimage,
      spw='3,7',
      mode='mfs',
      interactive=False,
      imsize=[640, 640],
      cell='0.07arcsec',
      niter=0)

dirtyim = dirtyimage + '.image'
viewer(dirtyim)

# --> still.. bg is brighter