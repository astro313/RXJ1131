#! /bin/csh -f

# Author : Daisy Leung
# Last Modified: Dec 13 2015
#
#
# Todo:
# ----
# - modify script to make moment maps, continuum map, and cont-sub. line maps
#
#
# History:
# --------
# Dec-13-2015:
#   - change USB to LSB (where line is)
#   - used script to make line cube
#   - customized for making line cube binned by 2 channels
# 08-30-2015
#   - Created, copied from RXJ1131
#
#
# Purpose:
# --------
# - select out source data and put on same velocity scale
# - Shift Velocity axis for the Source
# - run after calibrated data with script calibrateD3.csh
# - make SNR-conserved image WITHOUT using options=mosaic
# - make PB-corrected, flux-conserving image using options=mosaic
#
#
# Note:
# -----
# Calibrated set : cf0098.1D_215RXJ113.3_wide.mir in reduced_testflag/
# usb.co32.mom0.fits is line moment0 map without continuum subtraction
############################################################



####################################
######## USER DEFINE ##############
####################################
set starting_dir   = `pwd`
set dataPATH       = $HOME/Research/RXJ1131/CARMA
set dir_reduced   = "reduced_testflag" # Directory with reduced data
set cal_file       = "cf0098.1D_215RXJ113.3_wide.mir"
set vis            = $dataPATH/$dir_reduced/$cal_file
set dir_cal        = "reduced_pretty/bin2"           # output
if (!(-e $dataPATH/"$dir_cal"))     mkdir $dataPATH/$dir_cal

set outscience     = "RXJ1131_calibrated.mir"
set shiftedscience = "RXJ1131only_Calibrated_shift_vel.mir"


set lineops        = systemp,double #,mosaic
set contops        = systemp,mfs,double #,mosaic
set imsize         = 256
set cell           = 0.75
set imsize_mfs     = 256
set cellsize_mfs   = 0.75
set sig_two        = 2.0
set sig_conserve   = 3

set redshift        = 0.655
set Linefreq        = 345.7959899

# Imaging
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
    set offRegion_slab_bin2       = "boxes(30,30,110,110)(80,140)"
    set offRegion            = "boxes(30,30,110,110)"         # for continuum
    set offRegion_lineNoCont = "boxes(30,30,110,110)(80,140)"
else
    set src      = RXJ1131
    set contmap  = $src.contmap
    set contbeam = $src.contbeam
    set linmap   = $src.linmap
    set linbeam  = $src.linbeam
    set offRegion_slab_bin2       = "boxes(30,30,110,110)(80,140)"
    set offRegion            = "boxes(30,30,110,110)"         # for continuum
    set offRegion_lineNoCont = "boxes(30,30,110,110)(80,140)"
endif

set region_LSB           = "boxes(125,125,130,130)(55,65)"
# #set region_lineNoCont    = "boxes(124,129,132,135)(65,85)"
# set region_lineNoCont    = "boxes(124,129,132,135)"
# set region_mfs           = "boxes(149,153,151,155)"


cd $dataPATH/$dir_cal
echo "*** Current directory ***"
echo `pwd`
echo -n "*** HIT RETURN TO CONTINUE ***"
set ans = "$<"

echo "*** Change cleaning mask as needed ****"
echo -n "Cick Enter"
set ans ="$<"

# goto LSB
# ===============================
# Take out Science object ONLY
rm -rf $outscience
uvcat vis=$vis out=$outscience select="source(RXJ1131)"
uvlist options=spec vis=$outscience

# *******  Shift velocity of SOURCE UV  ********
#  Calculated rest_velocity, see aftercal_uv2fits.csh
set restfreq=`calc -f f0.10 "$Linefreq/($redshift+1)"`
rm -rf $shiftedscience
uvputhd vis=$outscience hdvar=restfreq varval=$restfreq out=$shiftedscience
uvlist options=spec vis=$shiftedscience  # to check

# goto cont
#####################################
echo " "
echo "*** Check for spectra of RXJ1131"
smauvspec vis=$shiftedscience interval=1000 hann=5 axis=chan,amp device=/xs nxy=3,5
echo -n "Cick Enter"
set ans ="$<"


#  ========= LSB ============
LSB:

if ($lineops =~ *'mosaic'* ) then
    set src     = RXJ1131.LSB_linewithCont.mosaic
else
    set src     = RXJ1131.LSB_linewithCont
endif

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

invert vis=$shiftedscience map=$LSBmap beam=$LSBbeam imsize=$imsize cell=$cell slop=1 robust=+2 options=$lineops line=vel,$nchan,$start,$width,$width
# theoretical rms noise: 1.190E-02

rm -rf $sourcefits
fits in=$LSBmap out=$sourcefits op=xyout
# CASAviewer to adjust region

histo in=$LSBmap region=$offRegion_slab_bin2  # off source
# rms: 1.60419E-02
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
exit
# goto cont

# ############################################
# ###########  MOMENT Maps for CO(3-2) ########
# BLAH contains line in LSB without continuum subtraction
#################################################
mom0:
set src=lsb
echo ""
echo "*** Inspect USB spectrum near source center and channels before subtracting cont."
echo ""
imspect in=$src.cm device=/xs hann=1 region="relcenter,boxes(-5,-5,5,5)(130,220)"

echo ""
echo "*** Moment 0 map of CO(3-2) in the inner quarter map, no continuum subtraction"
echo ""
rm -rf $src.co32_quarter.cm
imsub in=$src.cm region="quarter(150,200)" out=$src.co32_quarter.cm
cgdisp in=$src.co32_quarter.cm type=p xybin=1,1 device=/xs nxy=1,1 options=full,beambr,wedge,trlab,3val labtyp=arcsec,arcsec csize=1,1,1  range=0,0.045,lin,2
rm -rf $src.co32.mom0
moment mom=0 in=$src.co32_quarter.cm out=$src.co32.mom0
cgdisp in=$src.co32.mom0,$src.co32.mom0 type=p,c xybin=1,1 device=/xs nxy=1,1 options=full,beambr,wedge,trlab,3val labtyp=arcsec,arcsec csize=1,1,1  range=0,20,lin,2 slev=a,1.350 levs1=-6,-3,3,5,7,10
rm -rf $src.co32.mom0.fits
fits in=$src.co32.mom0 out=$src.co32.mom0.fits op=xyout

# #############################################################################
# Image plane continuum subtraction and make mom0 map of CO(3-2)
##############################################################################
# =======
# Alternative way to subtract continuum in image plane DO in casa
# 1: get average map with channels excluding the line
# 2: get mfs image with line channels
# 3: subtract the entire map with line and continuum from the image from step 1


echo ""
echo "*** Image plane continuum subtraction"
echo ""
rm -rf $src.imcontsub.line.cm $src.cont.cm
contsub in=$src.cm contchan="(1,150)(200,355)" out=$src.imcontsub.line.cm cont=$src.cont.cm region="images(150,200)" # mode="poly,1"  verbose=true

echo ""
echo "*** Inspect USB spectrum near source center and channels after subtracting cont."
echo ""
imspect in=$src.imcontsub.line.cm device=/xs hann=1 region="relcenter,boxes(-5,-5,5,5)"

echo ""
echo "*** Inspect USB Continuum RMS"
histo in=$src.cont.cm region="relcenter,boxes(-25,-25,-10,-10)"
# 3.33418E-3
cgdisp in=$src.cont.cm,$src.cont.cm type=p,c xybin=1,1 device=/xs nxy=1,1 options=full,beambr,wedge,3val labtyp=hms,dms csize=1,1,1 range=0,0.02,log,2 slev=a,0.003341 levs1=-6,-3,3,5,7,10 region="relcenter,boxes(-25,-25,25,25)"

echo ""
echo "*** USB Continuum by poly-fit of first order ***"
echo ""
rm -rf $src.cont.cm.fits
fits in=$src.cont.cm out=$src.cont.cm.fits op=xyout

echo ""
echo "*** Make Moment 0 map of CO(3-2) in the inner quarter map"
echo ""
rm -rf $src.150_200_co32_contsub.mom0 $src.imcontsub.line.quarter.cm
imsub in=$src.imcontsub.line.cm out=$src.imcontsub.line.quarter.cm # region="quarter"
histo in=$src.imcontsub.line.cm region="relcenter,boxes(-25,-25,-10,-10)"
echo -n "Cick Enter"
set ans ="$<"
moment mom=0 in=$src.imcontsub.line.quarter.cm out=$src.150_200_co32_contsub.mom0  		# << making weird image

echo ""
echo "*** RMS of Moment 0 map of CO(3-2) -- Off emission for contour"
echo ""
histo in=$src.150_200_co32_contsub.mom0 region="relcenter,boxes(-25,-25,-5,-5)"
echo -n "Cick Enter"
set ans ="$<"
echo ""
echo "*** MAX of Moment 0 map of CO(3-2) -- on emission for cgdisp range"
echo ""
# imlist in=$src.150_200_co32_contsub.mom0 options=statistics
echo -n "Cick Enter"
set ans ="$<"
histo in=$src.150_200_co32_contsub.mom0 # region="relcenter,boxes(-5,-5,5,5)"
echo -n "Cick Enter"
set ans ="$<"

rm -rf $src.co32_imcontsub.quarter.mom0
imsub in=$src.150_200_co32_contsub.mom0 out=$src.co32_imcontsub.quarter.mom0 region="relcenter,boxes(-15,-15,15,15)"

rm -rf $src.co32_imcontsub.quarter.mom0.fits
fits in=$src.co32_imcontsub.quarter.mom0  out=$src.co32_imcontsub.quarter.mom0.fits op=xyout
echo -n "Cick Enter"
set ans ="$<"


exit
# ===================== CONTINUUM =========================
# ======= chan BLAH of lsb.cm is line emission
#====== Separate out the Line-free and line channels in LSB =========
# ========== and clean them again separately ============
# =========================================================
cont:
if ($contops =~ *'mosaic'* ) then
    set src=RXJ1131.mosaic
else
    set src=RXJ1131
endif

# 1. In uv-plane: uvlin
echo " "
echo "*** Determine line-free spectrum"
echo " "
# smauvspec vis=$shiftedscience options=nobase,avall interval=1000 hann=5 axis=chan,amp device=1/xs nxy=1,1 line=channel,156,1000,2,2

# smauvspec vis=$shiftedscience options=nobase,avall interval=1000 hann=5 axis=freq,amp device=2/xs nxy=1,1 line=channel,156,1000,2,2
#smauvspec vis=$shiftedscience options=nobase,avall interval=1000 hann=5 axis=vel,amp device=3/xs nxy=1,1 line=channel,313,1000
# echo -n "Cick Enter"
# set ans ="$<"
# 1-761, 761-920, 1371-1521 are line free in $shiftedscience (LSB and USB combined uvfile)

echo " "
echo "*** Continuum subtraction"
rm -rf $src.uvcont $src.uvlin
uvlin vis=$shiftedscience chans=1,920,1371,1521 out=$src.uvcont mode=cont options=nowindow
uvlin vis=$shiftedscience chans=1,920,1371,1521 out=$src.uvlin line=channel,156,1000,2,2 mode=line options=nowindow 	# the line data after continuum subtraction			line=channel,313,1000

# to compare noise with USB cube before removing continuum
# which is consistent,
# will continue using line=channel,156,100,2,2 for mom0 (post-cont.sub) because the way I coded it uses the channels that corresponds to that cube (i.e., cube generated with line=channel,156,100,2,2)
rm -rf $src.uvlin.vel
uvlin vis=$shiftedscience chans=1,920,1371,1521 out=$src.uvlin.vel line=vel,355,-5042.6200,29.088,29.088 mode=line options=nowindow

echo " "
echo "*** Make a dirty image of line data and cont data, Binned in uvlin line=channel..."
echo " "
rm -rf $contmap $contbeam $linmap $linbeam $linmap.vel $linbeam.vel
invert vis=$src.uvlin map=$linmap beam=$linbeam imsize=$imsize cell=$cell sup=0 robust=2 options=$lineops
invert vis=$src.uvlin.vel map=$linmap.vel beam=$linbeam.vel imsize=$imsize cell=$cell sup=0 robust=2 options=$lineops
invert vis=$src.uvcont map=$contmap beam=$contbeam imsize=$imsize_mfs cell=$cellsize_mfs robust=2 options=$contops

echo " "
echo "*** Determine noise level"
echo " "
histo in=$linmap region=$offRegion
histo in=$linmap.vel region=$offRegion
histo in=$contmap region=$offRegion
echo -n "Cick Enter"
set ans ="$<"

echo ""
echo "*** Cleaning Image"
set cutoff        = `histo in=$linmap region=$offRegion_lineNoCont | grep Rms | awk '{printf "%.3e", $4}'`
set threshold     = `calc "$sig*$cutoff"`
set cutoff_vel        = `histo in=$linmap.vel region=$offRegion_lineNoCont | grep Rms | awk '{printf "%.3e", $4}'`
set threshold_vel     = `calc "$sig*$cutoff_vel"`
set cutoff_mfs    = `histo in=$contmap region=$offRegion | grep Rms | awk '{printf "%.3e", $4}'`
set threshold_mfs = `calc "$sig_conserve*$cutoff_mfs"`

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
    echo "*** Making psf for imfit for line without cont."
    ### create one synth beam out of ovro-ovro, bima-bima, ovro-bima
    # Get FWHM beam size of the mosaic
    # create psf using beam file from invert
    rm -rf $src.line.psf
    mospsf beam=$linbeam out=$src.line.psf
    echo " "
    echo "*** imfit a beam"
    set log_line = $src.line.psf.log
    rm -rf $log_line
    imfit in=$src.line.psf object=beam region="relcenter,boxes(-5,-5,5,5)" > $log_line
    set bmaj_line   = `grep "Major axis" $log_line | awk '{print $4}'`
    set bmin_line   = `grep "Minor axis" $log_line | awk '{print $4}'`
    set bpa_line    = `grep "  Position angle" $log_line | awk '{print $4}'`
    echo "Beam size = $bmaj_line x $bmin_line arcsec at PA = $bpa_line deg"

    rm -rf $src.line.vel.psf
    mospsf beam=$linbeam.vel out=$src.line.vel.psf
    echo " "
    echo "*** imfit a beam"
    set log_linevel = $src.line.vel.psf.log
    rm -rf $log_linevel
    imfit in=$src.line.vel.psf object=beam region="relcenter,boxes(-5,-5,5,5)" > $log_linevel
    set bmaj_linevel   = `grep "Major axis" $log_linevel | awk '{print $4}'`
    set bmin_linevel   = `grep "Minor axis" $log_linevel | awk '{print $4}'`
    set bpa_linevel    = `grep "  Position angle" $log_linevel | awk '{print $4}'`
    echo "Beam size = $bmaj_linevel x $bmin_linevel arcsec at PA = $bpa_linevel deg"


    echo " "
    echo "*** Clean the dirty images"
    echo " "
    rm -rf $linmap.cc  $contmap.cc $linmap.vel.cc
    mossdi map=$linmap beam=$linbeam out=$linmap.cc cutoff=$threshold region=$region_lineNoCont
    mossdi map=$linmap.vel beam=$linbeam.vel out=$linmap.vel.cc cutoff=$threshold_vel region=$region_lineNoCont
    mossdi map=$contmap beam=$contbeam out=$contmap.cc cutoff=$threshold_mfs region=$region_mfs


    echo " "
    echo "*** Restore cleaned line emission image"
    echo " "
    rm -rf $linmap.cm $linmap.res $linmap.vel.cm $linmap.vel.res
    restor model=$linmap.cc map=$linmap beam=$linbeam out=$linmap.cm fwhm=$bmaj_line,$bmin_line pa=$bpa_line
    restor model=$linmap.cc map=$linmap beam=$linbeam out=$linmap.res mode=residual
    restor model=$linmap.vel.cc map=$linmap.vel beam=$linbeam.vel out=$linmap.vel.cm fwhm=$bmaj_linevel,$bmin_linevel pa=$bpa_linevel
    restor model=$linmap.vel.cc map=$linmap.vel beam=$linbeam.vel out=$linmap.vel.res mode=residual

    echo " "
    echo "*** Restore cleaned continuum image"
    echo " "
    rm -rf $contmap.res $contmap.cm
    restor model=$contmap.cc map=$contmap beam=$contbeam out=$contmap.cm fwhm=$bmaj,$bmin pa=$bpa
    restor model=$contmap.cc map=$contmap beam=$contbeam out=$contmap.res mode=residual
else
    echo " "
    echo "*** Clean the dirty images"
    echo " "
    rm -rf $linmap.cc  $linmap.vel.cc $contmap.cc
    clean map=$linmap beam=$linbeam out=$linmap.cc cutoff=$threshold region=$region_lineNoCont
    clean map=$linmap.vel beam=$linbeam.vel out=$linmap.vel.cc cutoff=$threshold_vel region=$region_lineNoCont
    clean map=$contmap beam=$contbeam out=$contmap.cc cutoff=$threshold_mfs region=$region_mfs
    echo " "
    echo "*** Restore cleaned line emission image"
    echo " "
    rm -rf $linmap.cm $linmap.res $linmap.vel.cm $linmap.vel.res
    restor model=$linmap.cc map=$linmap beam=$linbeam out=$linmap.cm
    restor model=$linmap.cc map=$linmap beam=$linbeam out=$linmap.res mode=residual
    restor model=$linmap.vel.cc map=$linmap.vel beam=$linbeam.vel out=$linmap.vel.cm
    restor model=$linmap.vel.cc map=$linmap.vel beam=$linbeam.vel out=$linmap.vel.res mode=residual


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

echo ""
echo " *** # check distribution plotting histogram on line residual, inverted with chan "
echo ""
imhist in=$linmap.res region=$offRegion_lineNoCont device=/xs options=nbin,100
echo " "
echo "*** Check line cm rms"
echo " "
histo in=$linmap.cm region=$offRegion_lineNoCont
# imlist options=stat in=$linmap.cm region=$offRegion_lineNoCont
echo -n "*** HIT RETURN TO CONTINUE ***"
set ans = "$<"

echo ""
echo " *** # check distribution plotting histogram on line residual, inverted with vel "
echo ""
imhist in=$linmap.vel.res region=$offRegion_lineNoCont device=/xs options=nbin,100
echo " "
echo "*** Check line cm rms"
echo " "
histo in=$linmap.vel.cm region=$offRegion_lineNoCont
# imlist options=stat in=$linmap.vel.cm region=$offRegion_lineNoCont
echo -n "*** HIT RETURN TO CONTINUE ***"
set ans = "$<"


# set plotMax = `histo in=$lin.cm | grep Max | awk '{printf "%d",$3+0.5}'`
# set Plotrms = `histo in=$lin.cm region=$offRegion | grep Rms | awk '{printf "%.3e", $4}'`

# cgdisp in=$lin.cm,$lin.cm type=p,c xybin=1,1 device=/xs nxy=1,1 options=full,beambr,wedge,3val labtyp=hms,dms csize=1,1,1 range=0,plotMax,log,2 slev=a,Plotrms levs1=-6,-3,3,5,7,10 region="relcenter,boxes(-25,-25,25,25)"

echo " "
echo "*** Fit continuum source size in uv-plane"
echo " "
uvfit vis=$src.uvcont object=gaussian spar=7,-.6,4.67,11,7,-55 line=chan,15,1,100,100
uvfit vis=$src.uvcont object=point spar=7,-.6,4.67 line=chan,15,1,100,100
echo " "
echo "*** Fit continuum source size in image plane"
echo " "
rm -rf $contmap.cm.imfit.res
imfit in=$contmap.cm object=gaussian region="relcenter,boxes(-10,-10,10,10)" out=$contmap.cm.imfit.res options=residual
imstat in=$contmap.cm region="relcenter,boxes(-5,-5,5,5)" options=stat # not gaussin summed flux
echo -n "Cick Enter"
set ans ="$<"

rm -rf $linmap.cm.fits $linmap.vel.cm.fits $contmap.cm.fits
fits in=$linmap.cm out=$linmap.cm.fits op=xyout
fits in=$linmap.cm out=$linmap.vel.cm.fits op=xyout
fits in=$contmap.cm out=$contmap.cm.fits op=xyout

echo -n "Cick Enter"
set ans ="$<"

exit

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Momemt 0 of line after subtracting continuum in uvplane (chan invert)
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
momUV:
if ($lineops =~ *'mosaic'* ) then
    set src = RXJ1131.mosaic
else
    set src = RXJ1131
endif

echo "*** Inspect channels that contains CO(3-2) after subtracting continuum in UV"
echo ""
rm -rf $src.lin_55_100.cm $src.lin_65_85.cm $src.lin_72_90.cm
imspect in=$linmap.cm device=/xs hann=5 region="relcenter,boxes(-5,-5,5,5)" xaxis='channel'
echo -n "Cick Enter"
set ans ="$<"

echo "*** Extract channels for moment0 of CO(3-2) after subtracting continuum in UV"
echo ""
imsub in=$linmap.cm out=$src.lin_55_100.cm region="image(55,100)"
# make it tighter! 55-100 contains noise! got chan 65-85 in GILDAS of $linmap.cm, use image instead of quarter to avoid cropping
imsub in=$linmap.cm out=$src.lin_65_85.cm region="image(65,85)"
imsub in=$linmap.cm out=$src.lin_72_90.cm region="image(72,90)"

echo "*** Making moment0 of CO(3-2) after subtracting continuum in UV"
echo ""
rm -rf $linmap.cm.mom0 $src.lin65_85.cm.mom0 $src.lin72_90.cm.mom0
moment mom=0 in=$src.lin_55_100.cm out=$linmap.cm.mom0
moment mom=0 in=$src.lin_65_85.cm out=$src.lin65_85.cm.mom0
moment mom=0 in=$src.lin_72_90.cm out=$src.lin72_90.cm.mom0

rm -rf $linmap.cm.mom0.fits $src.lin65_85.cm.mom0.fits $src.lin72_90.cm.mom0.fits
fits in=$linmap.cm.mom0 out=$src.lin.cm.mom0.fits op=xyout
fits in=$src.lin65_85.cm.mom0 out=$src.lin65_85.cm.mom0.fits op=xyout
fits in=$src.lin72_90.cm.mom0 out=$src.lin72_90.cm.mom0.fits op=xyout
echo -n "Cick Enter"
set ans ="$<"

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Momemt 0 of line after subtracting continuum in uvplane (vel invert)
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
momVel:
if ($lineops =~ *'mosaic'* ) then
    set src = RXJ1131.mosaic
else
    set src = RXJ1131
endif

echo "*** Inspect channels that contains CO(3-2) after subtracting continuum in UV"
echo ""
rm -rf $src.lin_55_100.cm $src.lin_65_85.cm $src.lin_72_90.cm
imspect in=$linmap.cm device=/xs hann=5 region="relcenter,boxes(-5,-5,5,5)" xaxis='channel'
echo -n "Cick Enter"
set ans ="$<"

echo "*** Extract channels for moment0 of CO(3-2) after subtracting continuum in UV"
echo ""
imsub in=$linmap.cm out=$src.lin_55_100.cm region="image(55,100)"
imsub in=$linmap.cm out=$src.lin_65_85.cm region="image(65,85)"
imsub in=$linmap.cm out=$src.lin_72_90.cm region="image(72,90)"

echo "*** Making moment0 of CO(3-2) after subtracting continuum in UV"
echo ""
rm -rf $linmap.cm.mom0 $src.lin65_85.cm.mom0 $src.lin72_90.cm.mom0
moment mom=0 in=$src.lin_55_100.cm out=$linmap.cm.mom0
moment mom=0 in=$src.lin_65_85.cm out=$src.lin65_85.cm.mom0
moment mom=0 in=$src.lin_72_90.cm out=$src.lin72_90.cm.mom0

rm -rf $linmap.cm.mom0.fits $src.lin65_85.cm.mom0.fits $src.lin72_90.cm.mom0.fits
fits in=$linmap.cm.mom0 out=$src.lin.cm.mom0.fits op=xyout
fits in=$src.lin65_85.cm.mom0 out=$src.lin65_85.cm.mom0.fits op=xyout
fits in=$src.lin72_90.cm.mom0 out=$src.lin72_90.cm.mom0.fits op=xyout
echo -n "Cick Enter"
set ans ="$<"

#=========================
cd $starting_dir

goto end

end:
exit 0
end