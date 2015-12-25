# In order to importuvfits from a GILDAS-exported fits file, must supply /style CASA in GILDAS
#
# Purpose: Look into what are the difference in the header of the uvfits that uvmcmcfit is not understanding
#

from astropy.io import fits
from os.path import join

datapath = '/Users/admin/Research/RXJ1131/PdBI/data/10Sep15'
uvname = 'noCont_line_blue_oneBin.uvfits'
uvname_CASA = 'noCont_line_blue_oneBin_CASA.uvfits'
uvname_AIPS = 'noCont_line_blue_oneBin_styleAIPS.uvfits'
uvname_UVFITS = 'noCont_line_blue_oneBin_styleUV.uvfits'


def check(uvfitsname):
    im = fits.open(join(datapath, uvfitsname))
    im.info()
    im[0].data['DATA'].shape
    print im[0].header
    print("*"*40)

#     visfreq = im[1].data
    # check Frequency
    try:
        IF = im['AIPS FQ'].data['IF FREQ']
    except:
        print "No IF"
    try:
        nspw = im[0].data['DATA'][0, 0, 0, :, 0, 0, 0].size
        nfreq = im[0].data['DATA'][0, 0, 0, 0, :, 0, 0].size
        npol = im[0].data['DATA'][0, 0, 0, 0, 0, :, 0].size
    except:
        nfreq = im[0].data['DATA'][0, 0, 0, :, 0, 0].size
        npol = im[0].data['DATA'][0, 0, 0, 0, :, 0].size
    freq0 = im[0].header['CRVAL4']
    dfreq = im[0].header['CDELT4']
    cfreq = im[0].header['CRPIX4']

    im.close()
    print('Done')
check(uvname_UVFITS)
check(uvname_AIPS)
check(uvname_CASA)



def checkHDRwts(uvfitsname):
    """
    Check the weight column in header of uvfits file exported from GILDAS
    """
    uv = fits.open(join(datapath, uvfitsname))
    uv.info()
    print uv[0].header
    uv.close()
    print('Done')

# seems like the two below are the same
checkHDRwts(uvname_AIPS)
checkHDRwts(uvname_CASA)


checkHDRwts(uvname_UVFITS)


def pcd_uvfits(uvfitsName):
    """
    Get phase center
    """
    uvAfter = fits.open(join(datapath, uvfitsName))
    uvAfter.info()
    uvAfter[0].data['DATA'].shape
    hdr = uvAfter[0].header
    try:
        print uvAfter['AIPS SU '].data['RAEPO'][0]
        print uvAfter['AIPS SU '].data['DECEPO'][0]
        print 'AIPS'
    except:
        print hdr['CRVAL5']
        print hdr['CRVAL6']

    freq0_J = uvAfter[0].header['CRVAL4']
    dfreq_J = uvAfter[0].header['CDELT4']
    cfreq_J = uvAfter[0].header['CRPIX4']
    nfreq_J = uvAfter[0].data['DATA'][0, 0, 0, :, 0, 0].size
    npol_J = uvAfter[0].data['DATA'][0, 0, 0, 0, :, 0].size
    uvAfter.close()
    print('Done')

pcd_uvfits(uvname_AIPS)
pcd_uvfits(uvname_CASA)

# image fits

