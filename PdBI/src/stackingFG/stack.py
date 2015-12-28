'''

Stack line with CASA


To combine them is to use clean to put each spectral line cube on the same velocity scale. Then the cubes are smoothed to the same resolution as the cube with the coarsest resolution using imsmooth.
Finally, the spectral line cubes are averaged together using immath.

This script assumes that you have already continuum subtracted the data using uvcontsub.

If you need to continuum subtract in the image domain using imcontsub, it should be done after the clean, not after the spectral line cubes are averaged together. Be careful with your continuum subtraction -- small errors in continuum subtraction can lead to spurious lines.

Instead of smoothing after cleaning, another option would be to have all the images use the same restoring beam size. You can set this parameter in clean.


>> execfile('stack.py')


'''

import os

fixtable = False
cleanfile = True


cPath = os.path.abspath(os.curdir)
path = '/Users/admin/Research/RXJ1131/PdBI/data/28Dec15/'
print('-- change to directory: {:} ---').format(path)
os.chdir(path)

# run from casapy
uvAIPS = ['co21_trial1-noCont-ext69_83AIPS.uvfits',
          'co21_trial1-noCont-ext175_189AIPS.uvfits',
          'co21_trial1-noCont-ext220_234AIPS.uvfits']

uvFITS = ['co21_trial1-noCont-ext69_83.uvfits',
          'co21_trial1-noCont-ext175_189.uvfits',
          'co21_trial1-noCont-ext220_234.uvfits']

if fixtable:
    for uv1, uv2 in zip(uvAIPS, uvFITS):
        print uv1, uv2
        importuvfits(fitsfile=uv1, vis=uv1.replace('.uvfits', '.ms.contsub'))
        importuvfits(fitsfile=uv2, vis=uv2.replace('.uvfits', '.ms.contsub'))

        # handle the uvfits filler
        # copy antenna table

        print ("--- Copying tables --- ")
        os.system('cp -rf ' + uv1.replace('.uvfits', '.ms.contsub') + '/ANTENNA/* ' + uv2.replace('.uvfits', '.ms.contsub') + '/ANTENNA/')
        os.system('cp -rf ' + uv1.replace('.uvfits', '.ms.contsub') + '/FEED/* ' + uv2.replace('.uvfits', '.ms.contsub') + '/FEED/')


# put the spectral lines in different ms
linems = ['co21_trial1-noCont-ext69_83.ms.contsub',
          'co21_trial1-noCont-ext175_189.ms.contsub',
          'co21_trial1-noCont-ext220_234.ms.contsub']

# line rest freq (i.e. redshifted @ z=0)
restfreq = ['138.72265GHz',
            '139.78174GHz',
            '140.22446GHz']

velstart = '-150km/s'     # velocity to start the cube at
chanwidth = '30km/s'      # width of channels in the cube
nchan = 10                # number of channels in the cube
imsize = [1568, 1568]
cell = ['0.05arcsec']


# shallow clean if line is present, but for us, we want to get upper limit / don't have strong line. (better not clean any)
# threshold = '5mJy'
# niter=50
threshold = '80mJy'
niter = 0

# output image names
output_line_image = [linems[i].replace('.ms.contsub', '.contsub') for i in range(len(linems))]

if cleanfile:
    # create images for each spectral line data set on the same velocity scale
    for i in range(len(linems)):
        os.system('rm -rf ' + output_line_image[i] + '.*')
        clean(vis=linems[i],
              imagename=output_line_image[i],
              mode='velocity',
              start=velstart,
              nchan=nchan,
              width=chanwidth,
              interactive=False,
              restfreq=restfreq[i],
              niter=niter,
              imsize=imsize,
              cell=cell
              )

# getting the sizes of the beam
BMAJOR = []
BMINOR = []
BPA = []

for image in output_line_image:
    myimagename = image+'.image'
    BMAJOR.append(imhead(imagename=myimagename, mode='get', hdkey='beammajor'))
    BMINOR.append(imhead(imagename=myimagename, mode='get', hdkey='beamminor'))
    BPA.append(imhead(imagename=myimagename, mode='get', hdkey='beampa'))


# determine smooth
BMJ = max(BMAJOR)
idx = BMAJOR.index(BMJ)
BMN = BMINOR[idx]
BPAa = BPA[idx]

# smooth images
for image in output_line_image:
    inputimage = image+'.image'
    outputimage = image+'.image'+'.smooth'
    if image != output_line_image[idx]:
        imsmooth(imagename=inputimage,
                 outfile=outputimage,
                 kernel='gauss',
                 major=str(BMJ['value'])+BMJ['unit'],
                 minor=str(BMN['value'])+BMN['unit'],
                 pa=str(BPAa['value'])+BPAa['unit'],
                 targetres=True
                 )

    else:
        os.system('cp -ir ' + inputimage + ' ' + outputimage)

# immath
output_line_image = [output_line_image[i] + '.image' for i in range(len(output_line_image))]
myim = ['IM'+str(i) for i in range(len(output_line_image))]
myexp = '(' + ('+'.join(myim)) + ')/' + str(float(len(output_line_image)))

# combining the images
output_stacked_image = '_stack.image'
immath(imagename=output_line_image,
       outfile=output_stacked_image,
       mode='evalexpr',
       expr=myexp
       )

print('-- change back to code directory: {:} ---').format(cPath)
os.chdir(cPath)
