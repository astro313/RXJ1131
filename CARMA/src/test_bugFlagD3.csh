#! /bin/csh -f


# Author: Daisy
# Last Update :
# Purpose: Calibrate RX J1131 Data using Miriad
# Test how well it is calibrated without fixing the original flag, also run with old script that I have applied to other data
#
#

#################################
### USER DEFINED VARIABLES #####
#################################
set vis          = "cf0098.1D_215RXJ113.3.mir"

set ant          = 9
#set antpos      = $MIRCAT/baselines/carma/antpos.130404
set antpos       = $HOME/Research/RXJ1131/CARMA/raw/antpos.140124

set science      = RXJ1131
set pass         = 3C279
set gain         = 3C273

# The flux calibration procedure is as follows:
#    1) If flux_gaincal is set, that flux is used during gain calibration.
#    2) If flux_gaincal is not set, then flux_planet
#       is used to calibrate the fluxes.
#    3) If neither flux_planet and flux_gaincal are set, then the flux
#       for flux_gaincal is determined from the miriad flux table in
#       $MIR/cat/FluxSource.cat.
#
set flux_planet  = MARS      # fluxcal
set flux_gaincal = ""         # flux in Jy for gain calibrator
set flux_fluxcal = ""         # Flux in Jy for flux_planet. Set to "" if a planet.
set flux_object  = 3C273     # secondary

set edgechan     = "3"
set flagant      = "7"
set WIDE         = "channel,1,1,15,1"
# set win          = 4,5

# Time intervals (in minutes)
set interval_flux = 1.0       # [minutes]  Interval for flux calibration
set interval_gain = 14.6      # [minutes]  Interval for gain calibration
set interval_pb   = 1.0       # [minutes]  Interval for passband calibration


# Image parameters
set cell          = 0.75
set imsize        = 256
#set offRegion    = "boxes(-70,-70,-20,-20)"
set offRegion     = "abspixel,boxes(50,50,100,100)"
set robust        = 2
set dotsize       = 15
set badres        = 30
set nb_polyfit    = 0

#File outputs
set dataPATH      = $HOME/Research/RXJ1131/CARMA
set dir_vis       = "raw"    # Directory with raw data
set dir_reduced   = "reduced_testflag" # Directory with reduced data
set dir_tmp       = "tmp_D3config_testflag"     # Scratch directory
# Make directories
if (!(-e $dataPATH/"$dir_tmp"))     mkdir $dataPATH/$dir_tmp
if (!(-e $dataPATH/"$dir_reduced")) mkdir $dataPATH/$dir_reduced

# Move into temporary directory
set starting_dir = `pwd`
cd $dataPATH/$dir_tmp
echo "*** Current directory ***"
echo `pwd`
echo -n "*** HIT RETURN TO CONTINUE ***"
set ans = "$<"

# Get root name of raw visibility data
if ($vis:e == 'mir' || $vis:e == 'miriad') then
   set root = $vis:r
else
   set root = $vis
endif

set mir_raw        = "raw.mir"
set mir_cal_wide   = $dataPATH/$dir_reduced/${root}_wide.mir
set mir_cal_narrow = $dataPATH/$dir_reduced/${root}_narrow.mir
set out_base       = "base.vis"
set fname          = "base_flagMars.mir"

# Gain-calibrator flux used in mfcal
set mfcal_flux = ""
if ($flux_gaincal != "") set mfcal_flux = "flux=$flux_gaincal"

goto base
# ***********************
# **** COPY RAW DATA ****
# ***********************
# Only need to run once
set out = $mir_raw
# rm -rf $dataPATH/$dir_tmp/$out

if (!(-e $dataPATH/$dir_tmp/$out)) then
  echo ""
  echo "*** Making copy of raw data  (vis=$dataPATH/$dir_vis/$vis out=$dataPATH/$dir_tmp/$out)"
  uvcat vis=$dataPATH/$dir_vis/$vis out=$dataPATH/$dir_tmp/$out options=nowide
endif
set vis=$out

##############################
#### INSPECT DATA LOG ########
###############################
set obs_log 		   = "obsD.log"

rm -rf $obs_log
listobs vis=$vis log=$obs_log
uvindex vis=$vis
echo "*** Check auto flagged data, how many active flags you have now ***"
uvflag vis=$vis options=noapply flagval=flag
echo -n "*** HIT RETURN TO CONTINUE ***"
set ans = "$<"
# Good:        8243274
# Bad:        13247286


echo "*** Check what would happen if unflag data ***"
uvflag vis=$vis options=noapply flagval=unflag
echo -n "*** HIT RETURN TO CONTINUE, and unflag all data ***"
set ans = "$<"
echo " "
echo " *** Unflagging all data "
uvflag vis=$vis flagval=unflag
echo -n "*** HIT RETURN TO CONTINUE ***"
set ans = "$<"

uvlist vis=$vis options=spectra


########################
# look at data quality #
########################
echo -n "*** Do you want to inspect data before flagging (Y/N)? "
set ans = "$<"
if ($ans == "Y" || $ans == "y" || $ans == "1") then

  # look at the phase vs time: it should have a smooth trend
  uvplt vis=$vis select="-auto,source($gain)" device=/xs line=$WIDE axis=time,phase
  echo -n "*** HIT RETURN TO CONTINUE ***"
  set ans = "$<"

  # look at the amplitude vs. time: it should be constant and the same
  # for all antennas
  uvplt vis=$vis select="-auto,source($gain)" device=/xs line=$WIDE axis=time,amp
  echo -n "*** HIT RETURN TO CONTINUE ***"
  set ans = "$<"

  # look at the source phase vs time: if it has no continuum it
  # should be just random. Any trends would be due to false fringes
  uvplt vis=$vis select="-auto,source($science)" device=/xs line=$WIDE axis=time,phase
  echo -n "*** HIT RETURN TO CONTINUE ***"
  set ans = "$<"
endif

##################################
#### BASELINE CALIBRATION ########
##################################
rm -rf $out_base

base:
uvedit vis=$vis out=$out_base apfile=$antpos
set vis=$out_base

##########################
#### FLAGGING ############
##########################
# Flag shadowed baselines
echo ""
echo "*** Flagging shadowed baselines  (vis=$vis)"
csflag vis=$vis

echo ""
echo "*** Flag high elevation data, user-defined ant (vis=$vis)"
uvflag vis=$vis flagval=f select="el(85,90)"

# Flag used-specified bad antennas
if ($flagant != "") then
   echo ""
   echo "*** Flagging antennas: $flagant"
   uvflag vis=$vis flagval=f select="ant($flagant)"
endif

echo "*** Although `invert` weigh data by T_sys (options=systemp, it may also be advantageous to throw away data with large T_sys. "
echo " *** Flagging all records where T_sys > 1000K "
uvflag vis=$vis tsys=1000 flagval=flag
echo -n "*** HIT RETURN TO CONTINUE ***"
set ans = "$<"


echo ""
echo "*** Flagging $edgechan edge channels in all windows (vis=$vis)"
uvflag vis=$vis flagval=flag edge=$edgechan

#FLAG
uvflag vis=$vis flagval=flag select="ant(1)(4),time(08:10:00.0,08:20:00.0)"
uvflag vis=$vis flagval=flag select="ant(3)(4),time(09:30:00.0,09:56:47.0)"
uvflag vis=$vis flagval=flag select="ant(4)(5),time(07:30:00.0,09:00:00.0)"
uvflag vis=$vis flagval=flag select="ant(5)(11),time(07:30:00.0,07:40:00.0)"
uvflag vis=$vis flagval=flag select="ant(5)(12),time(08:10:00.0,08:20:00.0)"
echo -n "*** HIT RETURN TO CONTINUE ***"
set ans = "$<"

##################################
#### PLOT TRACK ########
##################################
#T_sys
echo -n "*** Do you want to plot the system temperatures (Y/N)? "
set ans = "$<"
if ($ans == "Y" || $ans == "y" || $ans == "1") then
    echo ""
    echo "*** Plotting tsys"
    echo "*** Look for antennas with high or variable tsys."
    echo "*** Note that all CARMA antennas are shown, even if they "
    echo "*** were not online during the track."
    varplt device=/xs vis=$vis nxy=5,3 yaxis=systemp
    echo -n "*** HIT RETURN TO CONTINUE ***"
    set ans = "$<"
endif

# Amplitudes
echo -n "*** Do you want to plot the Amplitudes (Y/N)? "
set ans = "$<"
if ($ans == "Y" || $ans == "y" || $ans == "1") then
  set win_track = 1   # Default window to plot
  echo ""
  echo "*** Raw amplitudes in window $win_track ***"
  echo "*** Look for integrations with abnormally high or low amplitudes."
  smauvplt device=/xs vis=$vis select="-source(noise),win($win_track),-auto"
  echo -n "*** HIT RETURN TO CONTINUE ***"
  set ans = "$<"
endif

# Phases
echo -n "*** Do you want to plot the Phases (Y/N)? "
set ans = "$<"
if ($ans == "Y" || $ans == "y" || $ans == "1") then
  set win_track = 1
  echo ""
  echo "*** Raw phases in window $win_track ***"
  echo "*** Look for incoherent phases on calibrators."
  smauvplt device=/xs vis=$vis select="-source(noise),-source($science),win($win_track),-auto" axis=time,phase
  echo -n "*** HIT RETURN TO CONTINUE ***"
  set ans = "$<"
endif

# Spectra , takes forever to finished plotting -- Don't set default on
echo -n "*** Do you want to plot the spectrum (Y/N)? "
set ans = "$<"
if ($ans == "Y" || $ans == "y" || $ans == "1") then
  echo ""
  echo "*** Plotting raw spectra for $pass ***"
  echo "*** Look for windows with poor signal to noise."
  smauvspec device=/xs vis=$vis select="source($pass),-auto" interval=10000 axis=chan,both
  echo -n "*** HIT RETURN TO CONTINUE ***"
  set ans = "$<"
endif

# first attempt to map the data
# rm -rf map0 beam0
# invert vis=$vis map=map0 beam=beam0 robust=$robust cell=$cell options=mosaic,mfs,systemp select="source($science)" imsize=$imsize
# cgdisp device=/xs in=map0 labtyp=arcsec options=full,beambr,wedge,3val region=quarter labtyp=hms,dms csize=1,1,1
# echo -n "*** HIT RETURN TO CONTINUE ***"
# set ans = "$<"


##############################################################################
# Flag data based on poor linecal phase v.s. time plot (see below)
##############################################################################
uvflag vis=$vis flagval=flag select="ant(10),time(09:00:00.0,09:30:00.0)"

##############################################################################
# Flag some data of Ant 12,
# phase noisy in some channels for some baselines after passband calibration (below) for the gain calibrator
# flag here, and re-calibrate the passband below
##############################################################################
uvflag vis=$vis flagval=flag select="ant(12)(5,14)"


####################################
# Flag some data of MARS
# Amp too high for some baselines
####################################
rm -rf $fname
cp -r $vis $fname

uvflag vis=$fname flagval=flag select="source($flux_planet),ant(2)(4)"
uvflag vis=$fname flagval=flag select="source($flux_planet),ant(3)(10)"
uvflag vis=$fname flagval=flag select="source($flux_planet),ant(5)(6)"
uvflag vis=$fname flagval=flag select="source($flux_planet),ant(6)(14)"
uvflag vis=$fname flagval=flag select="source($flux_planet),ant(8)(9,12,13,14,15)"
uvflag vis=$fname flagval=flag select="source($flux_planet),ant(9)(13,14)"
uvflag vis=$fname flagval=flag select="source($flux_planet),ant(10)(11,13)"
uvflag vis=$fname flagval=flag select="source($flux_planet),ant(11)(13)"
uvflag vis=$fname flagval=flag select="source($flux_planet),ant(12)(13,15)"

####################################
# Choose which vis file to use
# with MARS data flagged or not
####################################
linecal:
set vis = $fname
# set vis = out_base

##############################
#### SEPARATE NOISE ########
###############################

# This is needed since linecal should not be applied to noise source data.
# select out only the data (cut out the NOISE source and auto correlations
# leaving only the cross correlations)
set vis_astro = "no_noise_no_auto.vis"
rm -rf $vis_astro
echo ""
echo "*** Creating file with astronomical data only  (vis=$vis out=$vis_astro)"
uvcat vis=$vis select='-source(NOISE),-auto' out=$vis_astro



#################################################################
#### CALIBRATION ##############################################
# No auto correlation passband calibration (amplitude only)
# No noise source passband calibration (amplitude and phase)
##################################################################

# line length calibration
set out_line = "line_cal.vis"
rm -rf $out_line

echo ""
echo "*** Deriving linecal  (vis=$vis)"
linecal vis=$vis_astro

echo ""
echo "*** Plotting linecal phases ***"
smagpplt device=/xs vis=$vis_astro yaxis=phase options=wrap,nofit \
         yrange=-200,200 nxy=5,3
echo -n "*** HIT RETURN TO CONTINUE ***"
set ans = "$<"

echo ""
echo "*** Applying linecal  (vis=$vis_astro out=$out_line)"
uvcat vis=$vis_astro out=$out_line options=nopass
set vis_astro = $out_line


######## ######################## #######
# NOISE SOURCE CALIBRATION GOES HERE #####
# AUTO CORR CALIBRATION GOES HERE #####
# N/A in this data
##########################################

# ***********************************************************
# **** SPLIT ASTRONOMICAL DATA INTO WIDE AND NARROW BAND ****
# ***********************************************************
# Create file with wideband windows
echo ""
echo " *** NO NARROWBAND for RXJ1131 *** "
set vis_wide = $vis_astro

# ******************************************
# **** ASTRONOMICAL PASSBAND - WIDEBAND ****
# ******************************************
# This section will remove phase and amplitude offsets between windows.
# For the wide bands, we assume channel-by-channel passband is ok.

# Create copy of wide-band passband calibrator
set wb = "tmp_wb.mir"
rm -rf $wb
uvcat vis=$vis_wide out=$wb select="source($pass)" \
    options=nocal,nopass

# Derive wideband passband
echo ""
echo "*** Passband calibration for wideband channels  (vis=$wb source=$pass)"
set flux_fake = 1.0   # This can be an arbitrary number
mfcal vis=$wb select="source($pass)" \
      interval=$interval_pb refant=$ant flux=$flux_fake
gpcopy vis=$wb out=$vis_wide options=nocal      # meaning no copy antenna gain

# Show sources after passband
# Plot wide-band solutions - amplitude
echo ""
echo "*** Plotting mfcal passband amplitudes for wideband data  (vis=$vis_wide)"
smagpplt device=/xs vis=$vis_wide nxy=2,2 options=bandpass,nofit,dot dotsize=$dotsize
echo -n "*** HIT RETURN TO CONTINUE ***"
set ans = "$<"

# Plot wide-band solutions - phase
echo ""
echo "*** Plotting mfcal passband phases for wideband data  (vis=$vis_wide) "
smagpplt device=/xs vis=$vis_wide nxy=2,2 yrange=-200,200 \
         yaxis=phase options=wrap,bandpass,nofit,dot dotsize=$dotsize
echo -n "*** HIT RETURN TO CONTINUE ***"
set ans = "$<"

# Spectra on passband calibrator - wideband
echo ""
echo -n "*** Do you want to plot wide-band spectra for passband calibrator (Y/N)? "
set ans = "$<"
if ($ans == "Y" || $ans == "y" || $ans == "1") then
   echo ""
   echo "*** Plotting wideband data for $pass after passband calibration  (vis=$vis_wide)"
   echo "*** Amplitudes and phases should be flat, and there "
   echo "*** should be no amplitude/phase jumps across windows."
   smauvspec device=/xs vis=$vis_wide nxy=3,3 select="source($pass),-auto" axis=chan,both interval=10000
   echo -n "*** HIT RETURN TO CONTINUE ***"
   set ans = "$<"
endif

# Spectra on gain calibrator
echo ""
echo -n "*** Do you want to plot wide-band spec for gain calibrator (Y/N)? "
set ans = "$<"
if ($ans == "Y" || $ans == "y" || $ans == "1") then
   echo ""
   echo "*** Plotting wideband data for $gain after passband calibration  (vis=$vis_wide)"
   echo "*** Amplitudes and phases should be flat, and there "
   echo "*** should be no amplitude/phase jumps across windows."
   smauvspec device=/xs vis=$vis_wide nxy=3,3 select="source($gain),-auto" axis=chan,both interval=10000
   echo -n "*** HIT RETURN TO CONTINUE ***"
   set ans = "$<"
endif

# Apply passband - wideband
set mir_pb_wide = "astro_pb_wide.mir"
set out = $mir_pb_wide
rm -rf $out
echo ""
echo "*** Applying passband to wideband data  (vis=$vis_wide out=$out)"
uvcat vis=$vis_wide out=$out options=nocal
set vis_wide = $out

#######################################################
### Check passband calibration with point source ######
#######################################################
if ("$flux_object" != "") then
    set vis_wide = "astro_pb_wide.mir"
    set second = $flux_object
    # Continuum image - MFS
    echo ""
    echo "*** Making continum image of $second  (vis=$vis_wide)"
    echo "*** Units are in Janskys"
    echo ""
    echo "*** Should look like a Point source "
    set map="$second.cont.temp.map"
    set beam="$second.cont.temp.beam"
    set model="$second.cont.temp.cc"
    set cm="$second.cont.temp.cm"
    set res="$second.cont.temp.res"
    rm -rf $map $beam $model $res $cm
    invert vis=$vis_wide map=$map beam=$beam robust=$robust cell=$cell options=mosaic,mfs,systemp select="source($second)" imsize=$imsize
    cgdisp device=/xs in=$map labtyp=arcsec options=full,beambr,wedge,3val region=quarter labtyp=hms,dms csize=0.5,0.5,0.5
    echo -n "*** HIT RETURN TO CONTINUE ***"
    set ans = "$<"
endif

####################################
### FLUX CALIBRATION #########
####################################
fluxcal:
set vis_wide = "astro_pb_wide.mir"
if ("$flux_gaincal" == "" && "$flux_planet" != "") then
    echo ""
    echo "*** Starting flux calibration"
    set out_primary = "flux_primary.mir"
    rm -rf $out_primary
    echo ""
    echo "*** Creating file for $flux_planet  (vis=$vis_wide out=$out_primary)"
    uvcat vis=$vis_wide out=$out_primary options=nocal,nopass \
      select="source($flux_planet),-auto"

    # Create file with gain calibrator
    set out_gaincal = "flux_gaincal.mir"
    rm -rf $out_gaincal
    echo ""
    echo "*** Creating file for $gain  (vis=$vis_wide out=$out_gaincal)"
    uvcat vis=$vis_wide out=$out_gaincal options=nocal,nopass \
      select="source($gain),-auto"


    #Phase-only selfcal on primary flux calibrator
    set out = $out_primary:r_ph.mir
    rm -rf $out
    echo ""
    echo "*** Phase-only selfcal on flux calibrator  (vis=$out_primary out=$out)"
    selfcal vis=$out_primary refant=$ant interval=$interval_flux options=phase,apriori,noscale
    uvcat vis=$out_primary out=$out
    set out_primary = $out


    #plot amp vs uvdist using smauvamp
    echo ""
    echo "*** Plotting amplitudes vs. uvdistance for $flux_planet  (vis=$out_primary)"
    smauvamp device=/xs vis=$out_primary axis=uvd,amp options=zero nbin=20
    echo -n "*** HIT RETURN TO CONTINUE ***"
    set ans = "$<"


    #Phase only selfcal on gain calibrator
    set out = $out_gaincal:r_ph.mir
    rm -rf $out
    echo ""
    echo "*** Phase-only selfcal on gain calibrator  (vis=$out_gaincal out=$out)"
    selfcal vis=$out_gaincal interval=$interval_flux refant=$ant options=apriori,noscale,phase
    uvcat vis=$out_gaincal out=$out
    set out_gaincal = $out

    # Determine number of channels in spectrometer setup.
    # This is needed since bootflux requires only ony data point.
    set tmplog = tmp.uvlist
    rm -rf $tmplog
    uvlist vis=$out_primary options=var,full log=$tmplog
    set nchan=`grep nchan $tmplog | sed 's/.*nchan[ ]*://' | awk '{print $1}'`

    #Measure flux
    echo ""
    echo "*** Measuring flux using bootflux"
    set log = 'bootflux.log'
    set out = $gain.flux
    rm -rf $log $out
    set flux = ""
    set badres = 30
    if ($flux_fluxcal != "") set flux = "flux=$flux_fluxcal"
    bootflux vis=$out_primary,$out_gaincal select="source($flux_planet,$gain)" taver=1.01 primary=$flux_planet badres=$badres log=$log line=chan,1,1,$nchan,1 $flux
    grep Average $log | grep Flux | cut -c54-59 > $out
    echo ""
    echo -n "*** Measured flux in Jansky's on $gain is: "
    cat $out
    set mfcal_flux = "flux=@$out"
endif


####################################
### GAIN CALIBRATION #########
####################################
# apply gain calibration
echo "*** Gain calibration wideband data  (vis=$vis_wide source=$gain)"
mfcal vis=$vis_wide select="source($gain)" refant=$ant options=nopass $mfcal_flux interval=$interval_gain

echo "*** Plotting wideband amplitudes (vis=$vis_wide)"
echo "*** Amplitudes should vary smoothly in time and should have a value of ~ 1."
#plot result amp
smagpplt device=4/xs vis=$vis_wide nxy=5,3 options=dot,nofit dotsize=$dotsize yrange=0.5,1.2
echo -n "*** HIT RETURN TO CONTINUE ***"
set ans = "$<"

#plot result phase
echo "*** Plotting wideband phases  (vis=$vis_wide)"
echo "*** Phases should vary smoothly in time."
smagpplt device=/xs vis=$vis_wide nxy=5,3 yrange=-200,200 yaxis=phase options=wrap,dot,nofit dotsize=$dotsize
echo -n "*** HIT RETURN TO CONTINUE ***"
set ans = "$<"

# *******************
# **** COPY DATA ****
# *******************
# Copy calibrated data to reduced directory
echo ""
echo "*** Copying calibrated wide-band data ($vis_wide) to $dir_reduced"
puthd in=$vis_wide/senmodel value='GSV' type=ascii
set cal_wide = $mir_cal_wide
rm -rf $cal_wide
echo "*** Calibrated wideband data in $cal_wide"
cp -r $vis_wide $cal_wide


# ******************************
# **** PLOT CALIBRATED DATA ****
# ******************************
sanity:
set vis_wide = "astro_pb_wide.mir"
# Check data
# Phase vs. uv distance on phase gain calibrator
echo ""
echo "*** Plotting phase vs uvdistance on gain calibrator"
echo "*** Phases should be centered on zero. The smaller scatter the better!"
smauvplt device=/xs vis=$vis_wide select="source($gain)" options=nobase axis=uvdist,phase
echo -n "*** HIT RETURN TO CONTINUE ***"
set ans = "$<"

# Amplitude vs. uv distance on phase gain calibrator
echo ""
echo "*** Plotting amplitude vs uvdistance on gain calibrator"
echo "*** Amplitudes should be centered on flux calibrator flux."
smauvplt device=/xs vis=$vis_wide select="source($gain)" options=nobase axis=uvdist,amp
echo -n "*** HIT RETURN TO CONTINUE ***"
set ans = "$<"

# Self-calibration of the science target
echo "*** This will go into part of another imaging script *** "
echo "*** We'll continue making rough images for now"

temp_Im:
#Create Images in tmp
#set sources = `echo $science | sed 's/,/ /g'`
set vis_wide = "astro_pb_wide.mir"
set sources = $science
# Continuum image - MFS
echo ""
echo "*** Making continum image of $sources  (vis=$vis_wide)"
echo "*** Units are in Janskys"
set map="$sources.cont.temp.map"
set beam="$sources.cont.temp.beam"
set model="$sources.cont.temp.cc"
set cm="$sources.cont.temp.cm"
set res="$sources.cont.temp.res"
rm -rf $map $beam $model $res $cm
invert vis=$vis_wide map=$map beam=$beam robust=$robust cell=$cell options=mosaic,mfs,systemp select="source($sources)" imsize=$imsize
cgdisp device=/xs in=$map labtyp=arcsec options=full,beambr,wedge,3val region=quarter labtyp=hms,dms csize=0.5,0.5,0.5
#cgdisp device=2/xs in=$beam labtyp=arcsec options=full,beambr,wedge,3val labtyp=hms,dms csize=0.5,0.5,0.5   # dirty beam
echo -n "*** HIT RETURN TO CONTINUE ***"
set ans = "$<"
set fname="dirty_Dmfs.fits"
rm -rf $fname
fits in=$map out=$fname op=xyout


echo -n "*** HIT RETURN TO CONTINUE ***"
set ans = "$<"
set cutoff=`histo in=$map region=$offRegion | grep Rms | awk '{printf "%.3e", $4}'`
set minval=`histo in=$map | grep Min | awk '{printf "%3.1e",$3}'`
echo $cutoff
echo -n "*** HIT RETURN TO CONTINUE ***"
set ans = "$<"

# Get beam size
echo " "
echo "*** Making psf for imfit"
rm -rf $sources.psf
mospsf beam=$beam out=$sources.psf
echo " "
echo "*** imfit a beam"
set psf_log_mfs = $sources.psf.log
rm -rf $psf_log_mfs
imfit in=$sources.psf object=beam region=quarter > $psf_log_mfs
set bmaj=`grep "Major axis" $psf_log_mfs | awk '{print $4}'`
set bmin=`grep "Minor axis" $psf_log_mfs | awk '{print $4}'`
set bpa=`grep "  Position angle" $psf_log_mfs | awk '{print $4}'`
echo "Beam size = $bmaj x $bmin arcsec at PA = $bpa deg"
# Beam size = 6.333 x 4.166 arcsec at PA = 75.76 deg

echo ""
echo "*** Cleaning mfs Image with mossdi  **"
mossdi map=$map beam=$beam out=$model niters=10000 cutoff=$cutoff region=quarter
echo ""
echo "*** Restoring cleaned mfs Image **"
restor model=$model map=$map beam=$beam out=$cm fwhm=$bmaj,$bmin pa=$bpa
restor out=$res map=$map beam=$beam model=$model mode=residual

echo ""
echo "*** Plotting residuals with mossdi cleaning"
cgdisp device=/xs in=$res labtyp=arcsec options=full,beambr,wedge,3val region=quarter labtyp=hms,dms csize=1,1,1
echo -n "*** HIT RETURN TO CONTINUE ***"
set ans = "$<"

set f_mossdi = "$sources.cont.temp.mossdi"
rm -rf $f_mossdi.fits
fits in=$cm out=$f_mossdi.fits op=xyout


# echo ""
# echo "*** Cleaning mfs Image with CLEAN  **"
# set model="$sources.cont.tempCLEAN.cc"
# set cm="$sources.cont.tempCLEAN.cm"
# set res="$sources.cont.tempCLEAN.res"

# rm -rf $model $cm $res
# clean map=$map beam=$beam out=$model niters=10000 cutoff=$cutoff region=quarter
# echo ""
# echo "*** Restoring cleaned mfs Image **"
# restor map=$map beam=$beam model=$model out=$cm fwhm=$bmaj,$bmin pa=$bpa
# restor out=$res map=$map beam=$beam model=$model mode=residual
# echo -n "*** HIT RETURN TO CONTINUE ***"
# set ans = "$<"

# echo ""
# echo "*** Plotting residuals with CLEAN cleaning"
# cgdisp device=/xs in=$res labtyp=arcsec options=full,beambr,wedge,3val region=quarter labtyp=hms,dms csize=1,1,1
# echo -n "*** HIT RETURN TO CONTINUE ***"
# set ans = "$<"


# set f_clean  = "$sources.cont.temp.CLEAN"
# rm -rf $f_clean.fits
# fits in=$cm out=$f_clean.fits op=xyout

cd $starting_dir


