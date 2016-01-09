#! /bin/csh -f
##############################################################################
#
# Combining D2 and D3 dataset, using the mir/ files that have `uvputhd` applied, shifting restfreq to desired
#
# Author : Daisy Leung
# Last Modified: Jan 08 2016
#
#
# Note
# ----
# Line in LSB for CO(3-2) at 0.655
#
#
#
# History
# -------
# Jan 08 2016:
#   - change cell size 0.75 -> 0.25
#   - make MFS image to get limit on 1.5-mm continuum
# Dec 17 2015:
#   - Created
#
#
###############################################################################
#
#
#

###################################
######## USER DEFINE ##############
####################################
set starting_dir   = `pwd`
set dataPATH       = $HOME/Research/RXJ1131/CARMA
set dir_im         = "imagingD23"    # combined
if (!(-e $dataPATH/"$dir_im"))     mkdir $dataPATH/"$dir_im"

cd $dataPATH/"$dir_im"

set D3             = "/Users/admin/Research/RXJ1131/CARMA/reduced_pretty/bin2/RXJ1131only_Calibrated_shift_vel.mir"
set D2             = "../reduced_pretty_D2/bin2/RXJ1131only_Calibrated_shift_vel.mir"
# output name for combined vis
set D23            = "RXJ1131only_shiftVel_D23.mir"

# Imaging
set lineops        = systemp,double #,mosaic
set contops        = systemp,mfs,double #,mosaic
set imsize         = 256
set cell           = 0.25            # FOV:
set imsize_mfs     = 256
set cellsize_mfs   = 0.25
set sig_two        = 2.0
set sig_conserve   = 3

set redshift        = 0.655
set Linefreq        = 345.7959899

set LSBbeam       = "RXJ1131_LSB.beam"
set LSBmap        = "RXJ1131_LSB.map"
set sourcefits    = "dirtyLSB_RXJ1131.fits"
set LSB_cln       = "LSB_cln.fits"

if ($lineops =~ *'mosaic'* ) then
    set LSBbeam       = $LSBbeam.mosaic
    set LSBmap        = $LSBmap.mosaic
    set sourcefits    = "dirtyLSB_RXJ1131.mosaic.fits"
    set LSB_cln       = "LSB_cln.mosaic.fits"
endif

if ($contops =~ *'mosaic'* ) then
    set src      = RXJ1131
    set contmap  = $src.contmap.mosaic
    set contbeam = $src.contbeam.mosaic
    set linmap   = $src.linmap.mosaic
    set linbeam  = $src.linbeam.mosaic
    set offRegion_slab_bin2  = "boxes(30,30,110,110)(80,140)"
    set offRegion            = "boxes(30,30,110,110)"         # for continuum
    set offRegion_lineNoCont = "boxes(30,30,110,110)(80,140)"
else
    set src      = RXJ1131
    set contmap  = $src.contmap
    set contbeam = $src.contbeam
    set linmap   = $src.linmap
    set linbeam  = $src.linbeam
    set offRegion_slab_bin2  = "boxes(30,30,110,110)(80,140)"
    set offRegion            = "boxes(30,30,110,110)"         # for continuum
    set offRegion_lineNoCont = "boxes(30,30,110,110)(80,140)"
endif

# on source regions
set region_LSB     = "boxes(128,111,150,137)(55,65)"
set region_flat    = "boxes(128,111,150,137)"

cd $dataPATH/$dir_im
echo "*** Current directory ***"
echo `pwd`
echo -n "*** HIT RETURN TO CONTINUE ***"
set ans = "$<"

echo "*** Change cleaning mask as needed ****"
echo -n "Cick Enter"
set ans ="$<"

goto cont

# --------------------------------
# combine the two sets
# --------------------------------
uvcat vis=$D2,$D3 out=$D23

# --------------------------------
#  LSB
# ---------------------------------
LSB:

if ($lineops =~ *'mosaic'* ) then
    set src     = RXJ1131.LSB_linewithCont.mosaic
else
    set src     = RXJ1131.LSB_linewithCont
endif

# goto LSBcont

set bin      = 2
set onebin_v = 17.935
set bluest   = -2192.49
set redest   = 3134.27
set edge     = 3

set width     = `calc -f f0.3 "$bin*$onebin_v"`
set nchan     = `calc -i "($redest+abs($bluest)-2*($edge*$onebin_v))/$width"`
set start     = `calc -f f0.3 "$bluest+$edge*$onebin_v"`
echo "*** Making image with imsize=$imsize, cell=$cell, bin=$bin, velocity width=$width ***"
rm -rf $LSBmap $LSBbeam

invert vis=$D23 map=$LSBmap beam=$LSBbeam imsize=$imsize cell=$cell slop=1 robust=+2 options=$lineops line=vel,$nchan,$start,$width,$width
# Visibilities accepted: 21959
### Warning [invert]:  Visibilities rejected: 16576
### Warning [invert]:  Number of channels with no good data: 7
# Theoretical rms noise: 8.834E-03

rm -rf $sourcefits
fits in=$LSBmap out=$sourcefits op=xyout
# CASAviewer to adjust region

histo in=$LSBmap region=$offRegion_slab_bin2  # off source
#1.31957E-02
echo -n "Cick Enter"
set ans ="$<"


echo ""
echo "*** Cleaning Image"
set cutoff      = `histo in=$LSBmap region=$offRegion_slab_bin2 | grep Rms | awk '{printf "%.3e", $4}'`
echo $cutoff
set cutoff      = `calc "$sig_two*$cutoff"`
set threshold   = `calc "$sig_conserve*$cutoff"`

if ($lineops =~ *'mosaic'* ) then
    echo " "
    echo "*** Making psf for imfit"
    ### create one synth beam out of ovro-ovro, bima-bima, ovro-bima
    # Get FWHM beam size of the mosaic
    # create psf using beam file from invert
    set srcB = $LSBbeam:r
    rm -rf $srcB.psf
    mospsf beam=$LSBbeam out=$srcB.psf

    echo " "
    echo "*** imfit a beam"
    set log_lsb = $srcB.psf.log
    rm -rf $log_lsb
    imfit in=$srcB.psf object=beam region="relcenter,boxes(-5,-5,5,5)" > $log_lsb

    set bmaj        = `grep "Major axis" $log_lsb | awk '{print $4}'`
    set bmin        = `grep "Minor axis" $log_lsb | awk '{print $4}'`
    set bpa         = `grep "  Position angle" $log_lsb | awk '{print $4}'`
    echo "Beam size = $bmaj x $bmin arcsec at PA = $bpa deg"

    echo ""
    echo "*** Cleaning Image"
    rm -rf $src.cc
    mossdi map=$LSBmap beam=$LSBbeam out=$src.cc niters=10000 cutoff=$cutoff region=$region_LSB
    echo ""
    echo "*** Restoring Cleaned Image"
    rm -rf $src.cm
    restor map=$LSBmap beam=$LSBbeam out=$src.cm model=$src.cc fwhm=$bmaj,$bmin pa=$bpa
else
    echo ""
    echo "*** Cleaning Image"
    rm -rf $src.cc
    clean map=$LSBmap beam=$LSBbeam out=$src.cc niters=10000 cutoff=$cutoff region=$region_LSB
    rm -rf $src.cm
    # restore image
    restor model=$src.cc map=$LSBmap beam=$LSBbeam out=$src.cm
endif

### restore residual with clean beam
rm -rf $src.res
restor model=$src.cc map=$LSBmap beam=$LSBbeam out=$src.res mode=residual

echo ""
echo " *** # check distribution plotting histogram on line residual "
echo ""
imhist in=$src.res region=$offRegion_slab_bin2 device=/xs options=nbin,100
histo in=$src.cm region=$offRegion_slab_bin2
echo -n "*** HIT RETURN TO CONTINUE ***"
set ans = "$<"

rm -rf $LSB_cln
fits in=$src.cm out=$LSB_cln op=xyout

###############################################################################
#
# No continuum at position of CO emission, see if we can get a limit from cont. map
#
#############################################################################
cont:
if ($contops =~ *'mosaic'* ) then
    set src=RXJ1131.mosaic
else
    set src=RXJ1131
endif

echo " "
echo "*** Continuum subtraction (ISSUE) check number of uv points ***"
rm -rf $src.uvcont
uvlin vis=$D23 chans=1,118,160,600 out=$src.uvcont order=1 mode=cont options=nowindow
echo -n "Cick Enter"
set ans ="$<"

echo " "
echo "*** Make a dirty image cont data "
echo " "
rm -rf $contmap $contbeam
invert vis=$src.uvcont map=$contmap beam=$contbeam imsize=$imsize_mfs cell=$cellsize_mfs robust=2 options=$contops sup=0
# Making MFS images
# Visibilities accepted: 10616694
### Warning [invert]:  Visibilities rejected: 13429146
# Mean Frequency(GHz):     216.
#Theoretical rms noise: 5.582E-04

echo " "
echo "** export dirty cont. map"
rm -rf $contmap.dirty.fits
fits in=$contmap out=$contmap.dirty.fits op=xyout
echo -n "Cick Enter"
set ans ="$<"

echo " "
echo "*** Determine noise level"
echo " "
histo in=$contmap region=$offRegion       # rms: 8.30807E-04
echo -n "Cick Enter"
set ans ="$<"


echo ""
echo "*** Cleaning Image"
set cutoff_mfs    = `histo in=$contmap region=$offRegion | grep Rms | awk '{printf "%.3e", $4}'`
set threshold_mfs = `calc "$cutoff_mfs"`

if ($contops =~ *'mosaic'* ) then   # should make it check with $lineops as well
    echo " "
    echo "*** Making psf for imfit for continuum"
    ### create one synth beam out of ovro-ovro, bima-bima, ovro-bima
    # Get FWHM beam size of the mosaic
    # create psf using beam file from invert
    rm -rf $src.psf
    mospsf beam=$contbeam out=$src.psf

    echo " "
    echo "*** imfit a beam"
    set log_cont = $src.psf.log
    rm -rf $log_cont
    imfit in=$src.psf object=beam region="relcenter,boxes(-5,-5,5,5)" > $log_cont

    set bmaj        = `grep "Major axis" $log_cont | awk '{print $4}'`
    set bmin        = `grep "Minor axis" $log_cont | awk '{print $4}'`
    set bpa         = `grep "  Position angle" $log_cont | awk '{print $4}'`
    echo "Beam size = $bmaj x $bmin arcsec at PA = $bpa deg"

    echo " "
    echo "*** Clean the dirty images"
    echo " "
    rm -rf $contmap.cc
    mossdi map=$contmap beam=$contbeam out=$contmap.cc cutoff=$threshold_mfs region=$region_flat


    echo " "
    echo "*** Restore cleaned continuum image"
    echo " "
    rm -rf $contmap.res $contmap.cm
    restor model=$contmap.cc map=$contmap beam=$contbeam out=$contmap.cm fwhm=$bmaj,$bmin pa=$bpa
    restor model=$contmap.cc map=$contmap beam=$contbeam out=$contmap.res mode=residual
else
    echo " "
    echo "*** Clean the dirty image"
    echo " "
    rm -rf $contmap.cc
    clean map=$contmap beam=$contbeam out=$contmap.cc cutoff=$threshold_mfs region=$region_flat

    echo " "
    echo "*** Restore cleaned continuum image"
    echo " "
    rm -rf $contmap.res $contmap.cm
    restor model=$contmap.cc map=$contmap beam=$contbeam out=$contmap.cm
    restor model=$contmap.cc map=$contmap beam=$contbeam out=$contmap.res mode=residual
endif

echo ""
echo " *** # check distribution plotting histogram on continuum residual "
echo ""
imhist in=$contmap.res region=$offRegion device=/xs options=nbin,100
echo " "
echo "*** Check continuum cm rms"
echo " "
imlist options=stat in=$contmap.cm region=$offRegion
echo -n "*** HIT RETURN TO CONTINUE ***"
set ans = "$<"

rm -rf $contmap.cm.fits
fits in=$contmap.cm out=$contmap.cm.fits op=xyout

exit
end:
exit 0
end