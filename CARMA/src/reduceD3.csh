#! /bin/csh -f


# Author: Daisy
# Last Update :
# Purpose: Calibrate RXJ1131s CARMA Data using Miriad, for the 3rd obs
# Last Edit: Sept 2 2015 --> need more flagging?
#
# Purpose: Fully calibrate the dataset
#
# Note:
# - may want to change interval_gain?
# - The phase calibrator is far away for this set
#
#
#
# make dirty and clean mfs image of the whole image cube
#
# Outputs: Everything is under dir_tmp except
# for the fully calibrated .mir which is under dir_reduced
#
# Further analysis such as shifting header rest freq and imaging
# is in another script
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
set win          = 4,5

# Time intervals (in minutes)
set interval_flux    = 1.0       # [minutes]  Interval for flux calibration
set int_fluxgain_phs = 0.5       # phase-selfcal with short interval, to ensure phs wrapping in time not wipe out the passband
set interval_gain    = 14.6      # [minutes]  Interval for gain calibration, amount spent on each cal-sci taget cycle
set interval_pb      = 1.0       # [minutes]  Interval for passband calibration


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
set starting_dir  = `pwd`
cd $dir_tmp
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

# Gain-calibrator flux used in mfcal
set mfcal_flux = ""
if ($flux_gaincal != "") set mfcal_flux = "flux=$flux_gaincal"

goto fluxcal
# ***********************
# **** COPY RAW DATA ****
# ***********************
# Only need to run once
set out = $mir_raw
rm -rf $dataPATH/$dir_tmp/$out
# echo ""
echo "*** Making copy of raw data  (vis=$dataPATH/$dir_vis/$vis out=$dataPATH/$dir_tmp/$out)"
uvcat vis=$dataPATH/$dir_vis/$vis out=$dataPATH/$dir_tmp/$out options=nowide
set vis=$out

##############################
#### INSPECT DATA LOG ########
###############################
# set obs_log 		   = "obsD.log"

# rm -rf $obs_log
# listobs vis=$vis log=$obs_log
# uvindex vis=$vis
# uvlist vis=$vis options=spectra
echo "*** Check auto flagged data, how many active flags you have now ***"
uvflag vis=$vis options=noapply flagval=flag
# Bad:        13247286
echo -n "*** HIT RETURN TO CONTINUE ***"
set ans = "$<"

echo " "
echo " *** Unflagging all data "
uvflag vis=$vis flagval=unflag

uvlist vis=$vis options=spectra

########################
# look at data quality #
########################
# # look at the phase vs time: it should have a smooth trend
uvplt vis=$vis select="-auto,source($gain)" device=/xs line=$WIDE axis=time,phase options=unwrap
echo -n "*** HIT RETURN TO CONTINUE ***"
set ans = "$<"

# # look at the amplitude vs. time: it should be constant and the same
# # for all antennas
uvplt vis=$vis select="-auto,source($gain)" device=/xs line=$WIDE axis=time,amp
echo -n "*** HIT RETURN TO CONTINUE ***"
set ans = "$<"

# # look at the source phase vs time: if it has no continuum it
# # should be just random. Any trends would be due to false fringes
uvplt vis=$vis select="-auto,source($science)" device=/xs line=$WIDE axis=time,phase
echo -n "*** HIT RETURN TO CONTINUE ***"
set ans = "$<"

echo "*** Plot amplitude, phase versus time for lsb wideband"
uvplt vis=$vis select="-auto,source($gain)" device=/xs axis=time,ampl line=wide,1,1
uvplt vis=$vis select="-auto,source($gain)" device=/xs axis=time,phase line=wide,1,1
echo -n "*** HIT RETURN TO CONTINUE ***"
set ans = "$<"

# all sources T_sys
varplt vis=$vis device=/xs yaxis=systemp nxy=5,3 yrange=0,2000 options=compress


echo " *** Inspect why flag so many data for this set? "
# inspect all objects even for the science target: it probably be random, unless there are false fringes
# or it is a very strong source
foreach i ($gain $science $pass $flux_planet)
  # phase vs. time
  echo " *** Phase v.s. time WIDEBAND"
  smauvplt vis=$vis select="-auto,source("$i")" device=/xs line=$WIDE axis=time,ph interval=3
  echo -n "*** HIT RETURN TO CONTINUE ***"
  set ans = "$<"
  # amplitude vs. time
  echo " *** Amplitude v.s. time window=$win"
  smauvplt vis=$vis select="-auto,source("$i")" device=/xs axis=time,amp "select=win($win)"
  echo -n "*** HIT RETURN TO CONTINUE ***"
  set ans = "$<"
  echo " *** WIDEBAND"
  smauvplt vis=$vis select="-auto,source("$i")" device=/xs line=$WIDE axis=time,amp interval=3
  echo -n "*** HIT RETURN TO CONTINUE ***"
  set ans = "$<"


  echo " *** amp,ph vs. channel"
  smauvspec vis=$vis select="-auto,source("$i")" axis=chan,phase interval=999 nxy=5,3 device=/xs
  echo -n "*** HIT RETURN TO CONTINUE ***"
  set ans = "$<"
  smauvspec vis=$vis select="-auto,source("$i")" axis=chan,amp interval=999 nxy=5,3 device=/xs
  echo -n "*** HIT RETURN TO CONTINUE ***"
  set ans = "$<"
  smauvspec vis=$vis select="-auto,source("$i")" axis=chan,phase line="channel,282,1,1,1" device=/xs interval=999 yrange=-180,180 nxy=5,3
  echo -n "*** HIT RETURN TO CONTINUE ***"
  set ans = "$<"

  echo " *** in velocity space... notice with uvlist-spectra, do we need to patch? "
  smauvspec vis=$vis select="-auto,source("$i")" device=1/xs interval=999 nxy=5,3 axis=velocity,amp "select=win(4,5,6)"
end
##################################
#### BASELINE CALIBRATION ########
##################################
base:
set out_base='base.vis'
rm -rf $out_base
uvedit vis=$vis out=$out_base apfile=$antpos
set vis=$out_base

# goto linecal
##############################################################################
############################## FLAGGING ######################################
# uvplt, varplt, uvflag, blflag, uvimage, uvlist option=bfmask (identify why data was flagged)
##############################################################################
# Flag shadowed baselines
echo ""
echo "*** Flagging shadowed baselines  (vis=$vis)"
csflag vis=$vis   # designed for CARMA, works for heterogenous array
echo -n "*** HIT RETURN TO CONTINUE ***"
set ans = "$<"
echo "*** check flags "
uvflag vis=$vis flagval=flag options=noapply
echo -n "*** HIT RETURN TO CONTINUE ***"
set ans = "$<"

echo ""
echo "*** Flag high elevation data (vis=$vis)"
uvflag vis=$vis flagval=flag select="el(85,90)"
echo -n "*** HIT RETURN TO CONTINUE ***"
set ans = "$<"
echo "*** check flags "
uvflag vis=$vis flagval=flag options=noapply
echo -n "*** HIT RETURN TO CONTINUE ***"
set ans = "$<"

echo ""
echo "*** Creating image cube from (vis=$vis), to guide flagging"
echo "*** Look for a change in noise, regions of pure 0s, vertical spikes (a.k.a. birdies), horizontal spikes (bad baselines or antennae). These will potentially all have to be flagged."
echo "*** Overall noise increase that is the result of a higher T_sys will be accounted for though (see invert options=systemp) "
# useful modes: amplitude (default), channel (x-axis), baselines(Y), time (Z)
rm -rf BLAH180
uvimage vis=$vis out=BLAH180
echo "  You do need to start ds9 manually yourself."
set xpasets=(`which xpaset`)
if ($#xpasets > 1) then
  echo 'It does not look like your $path knows about xpaset.'
  echo 'You need the xpatools for this from SAO.'
  exit 0
endif
ds9 &
sleep 5  # wait for ds9 to load
mirds9 BLAH180

# example of birdie (often antenna based single chan with high amp) can be flagged with `uvflag line=chan,n,start,width,step select="ant(blah)" flagval=f`

set flagchan=""
if ($flagchan != "") then
  echo ""
  echo " *** Flagging these channels in all antennas "
  # change the input before running...
  uvflag line=chan,n,start,width,step flagval=f
endif
echo -n "*** HIT RETURN TO CONTINUE ***"
set ans = "$<"

echo "*** Although `invert` weigh data by T_sys (options=systemp, it may also be advantageous to throw away data with large T_sys. "
echo " *** Flagging all records where T_sys > 1000K "
uvflag vis=$vis tsys=1000 flagval=flag
echo -n "*** HIT RETURN TO CONTINUE ***"
set ans = "$<"

# echo ""
# echo "*** TESTING: Flag based on tracking errors (vis=$vis)"
# # axisrms contains tracking error (in arcsec, in Az and El) for each antenna
# varplt vis=$vis device=/xs yaxis=axisrms options=overlay yrange=0,30
# echo -n "*** HIT RETURN TO CONTINUE ***"
# set ans = "$<"

# Flag used-specified bad antennas
if ($flagant != "") then
   echo ""
   echo "*** Flagging antennas: $flagant"
   uvflag vis=$vis flagval=flag select="ant($flagant)"
endif
echo -n "*** HIT RETURN TO CONTINUE ***"
set ans = "$<"


echo ""
echo "*** Flagging $edgechan edge channels in all windows (vis=$vis)"
uvflag vis=$vis flagval=flag edge=$edgechan
echo -n "*** HIT RETURN TO CONTINUE ***"
set ans = "$<"

#FLAG
uvflag vis=$vis flagval=flag select="ant(12),time(07:50:00.0,08:05:00)"
uvflag vis=$vis flagval=flag select="ant(1)(4),time(08:10:00,08:20:00)"
uvflag vis=$vis flagval=flag select="ant(1)(12),time(07:27:39.5,07:45:00)"
uvflag vis=$vis flagval=flag select="ant(2)(11),time(07:27:39.5,07:30:00.0)"
uvflag vis=$vis flagval=flag select="ant(3)(4),time(09:30:00,09:56:47.0)"
uvflag vis=$vis flagval=flag select="ant(4)(5),time(07:40:00.0,08:00:00)"


echo ""
echo "*** Copy unflagged vis after basic flagging to $vis.basicflags "
rm -rf $vis.basicflags
uvcat vis=$vis out=$vis.basicflags options=unflagged
echo -n "*** HIT RETURN TO CONTINUE ***"
set ans = "$<"

echo ""
echo "*** Copy unflagged vis after basic flagging to $vis.basicflags.cont to continue flagging those with wide spread phase"
rm -rf $vis.basicflags.cont
uvcat vis=$vis out=$vis.basicflags.cont options=unflagged
echo -n "*** HIT RETURN TO CONTINUE ***"
set ans = "$<"


# wide sperad in phase
uvflag vis=$vis.basicflags.cont flagval=flag select="ant(1)(14),time(07:27:39.5,07:40:00)"
uvflag vis=$vis.basicflags.cont flagval=flag select="ant(4)(15)"
uvflag vis=$vis.basicflags.cont flagval=flag select="ant(5)(11)"
 uvflag vis=$vis.basicflags.cont flagval=flag select="ant(5)(12),time(07:30:00,07:40:00)"
uvflag vis=$vis.basicflags.cont flagval=flag select="ant(5)(15)"
uvflag vis=$vis.basicflags.cont flagval=flag select="ant(6)(15)"
uvflag vis=$vis.basicflags.cont flagval=flag select="ant(9)(11)"
uvflag vis=$vis.basicflags.cont flagval=flag select="ant(11)(14)"


echo ""
echo "*** Creating image cube from (vis=$vis), to guide flagging"
echo "*** Look for a change in noise, regions of pure 0s, vertical spikes (a.k.a. birdies), horizontal spikes (bad baselines or antennae). These will potentially all have to be flagged."
echo "*** Overall noise increase that is the result of a higher T_sys will be accounted for though (see invert options=systemp) "
# useful modes: amplitude (default), channel (x-axis), baselines(Y), time (Z)
rm -rf BLAH180_afterFlags
uvimage vis=$vis.basicflags.cont out=BLAH180_afterFlags
echo "  You do need to start ds9 manually yourself."
ds9 &
sleep 5
mirds9 BLAH180_afterFlags

echo -n " *** Remove flagged data for plotting ****"
echo ""
echo "*** Copy unflagged vis after flagging wide spread to $vis.flagWidephs4Plot to plot only unflagged data"
rm -rf $vis.flagWidephs4Plot
uvcat vis=$vis.basicflags.cont out=$vis.flagWidephs4Plot options=unflagged
echo -n "*** HIT RETURN TO CONTINUE ***"
set ans = "$<"

# exit

#################################################
################ Validate Flags ##############
#################################################
verify:
foreach i ($vis.basicflags $vis.flagWidephs4Plot)
  echo -n "*** HIT RETURN TO CONTINUE ***"
  set ans = "$<"
  # vs. channel and time.  average over baselines
  smauvspec vis=$i device=/xs select='ant(7)' axis=chan,bo options=avall
  echo -n "*** HIT RETURN TO CONTINUE ***"
  set ans = "$<"
  # vs. channel and baseline.  average over time
  smauvspec vis=$i device=/xs select='ant(7)' axis=ch,bo interval=1000
  echo -n "*** HIT RETURN TO CONTINUE ***"
  set ans = "$<"
  # vs. baseline.  average over channel and time
  smauvplt vis=$i device=/xs axis=uvdistance,amp options=nobase select='-auto'
  echo -n "*** HIT RETURN TO CONTINUE ***"
  set ans = "$<"
end
##################################
#### PLOT TRACK ########
##################################
#T_sys
# set vis=$vis.basicflags
# set vis=$vis.flagWidephs4Plot
track:
foreach name ($vis.basicflags $vis.flagWidephs4Plot)
  echo -n "*** Do you want to plot the system temperatures (Y/N)? "
  set ans = "$<"
  if ($ans == "Y" || $ans == "y" || $ans == "1") then
      echo ""
      echo "*** Plotting tsys"
      echo "*** Look for antennas with high or variable tsys."
      echo "*** Note that all CARMA antennas are shown, even if they "
      echo "*** were not online during the track."
      smavarplt device=/xs vis=$name nxy=5,3 yaxis=systemp
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
    smauvplt device=/xs vis=$name select="-source(noise),win($win_track),-auto"
    echo -n "*** HIT RETURN TO CONTINUE ***"
    set ans = "$<"
  endif

  # Phases
  echo -n "*** Do you want to plot the Phases (Y/N)? "
  set ans = "$<"
  if ($ans == "Y" || $ans == "y" || $ans == "1") then
    echo ""
    echo "*** Raw phases in window $win_track ***"
    echo "*** Look for incoherent phases on calibrators."
    smauvplt device=/xs vis=$name select="-source(noise),win($win_track),-auto"  axis=time,phase
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
    smauvspec device=/xs vis=$name select="source($pass),-auto" interval=10000 axis=chan,both
    echo -n "*** HIT RETURN TO CONTINUE ***"
    set ans = "$<"
  endif
end

# first attempt to map the data
# rm -rf map0 beam0
# invert vis=$vis map=map0 beam=beam0 robust=$robust cell=$cell options=mosaic,mfs,systemp select="source($science)" imsize=$imsize
# cgdisp device=/xs in=map0 labtyp=arcsec options=full,beambr,wedge,3val region=quarter labtyp=hms,dms csize=1,1,1
# echo -n "*** HIT RETURN TO CONTINUE ***"
# set ans = "$<"
# exit

##############################
#### SEPARATE NOISE ########
###############################
linecal:
echo $vis
# This is needed since linecal should not be applied to noise source data.
# select out only the data (cut out the NOISE source and auto correlations
# leaving only the cross correlations)
set vis_astro = "no_noise_no_auto.vis"
set vis_noise = "noise_no_auto.vis"
rm -rf $vis_astro $vis_noise
echo ""
echo "*** Creating file with astronomical data only  (vis=$vis out=$vis_astro)"
uvcat vis=$vis select='-source(NOISE),-auto' out=$vis_astro
echo ""
echo "*** Creating file with noise source data  (vis=$vis out=$vis_noise)"
uvcat vis=$vis select="source(noise),-auto" out=$vis_noise


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
         yrange=-180,180 nxy=5,3
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
set out = "split_wide.mir"
rm -rf $out
echo ""
echo "*** Creating wideband data file  (vis=$vis_astro out=$out)"
uvcat vis=$vis_astro out=$out select="-auto"
set vis_wide = $out

# NO NARROWBAND for RXJ1131

# ************************************************
# **** ASTRONOMICAL PASSBAND - WIDEBAND **********
# **** use short int, esp for 1mm (default is 5 mm)
# ************************************************
# This section will remove phase and amplitude offsets between windows.
# For the wide bands, we assume channel-by-channel passband is ok.

# Create copy of wide-band passband calibrator
set wb = "tmp_wb.mir"
rm -rf $wb
uvcat vis=$vis_wide out=$wb select="source($pass)" \
    options=nocal,nopass

# Inspect passband calibrator
smauvspec vis=$wb device=/xs interval=999 nxy=5,3 axis=channel,phase
smauvspec vis=$wb device=/xs interval=999 nxy=5,3 axis=channel,amp
# and in vel space... notice with uvlist-spectra, do we need to patch?
smauvspec vis=$wb device=/xs interval=999 nxy=5,3 axis=velocity,amp select="win(4,5,6)"

# inspect amps of source
echo ""
echo "*** Check for bad points on source "
smauvplt vis=$vis_wide select="source($science)" device=/xs axis=time,amp select="win($win)"


# Derive wideband passband
echo ""
echo "*** Passband calibration for wideband channels  (vis=$wb source=$pass)"
# Antenna-based
set flux_fake = 1.0   # This can be an arbitrary number
mfcal vis=$wb select="source($pass)" \
      interval=$interval_pb refant=$ant flux=$flux_fake
gpcopy vis=$wb out=$vis_wide options=nocal      # meaning no copy antenna gain

# Show sources after passband calibration
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
   smauvspec device=/xs vis=$vis_wide nxy=3,3 select="source($gain),-auto" axis=ch,bo interval=10000
   echo -n "*** HIT RETURN TO CONTINUE ***"
   set ans = "$<"
  echo ""
   echo "*** Plotting wideband data for $gain after passband calibration in time (vis=$vis_wide)"
   echo "*** Amplitudes and phases should be flat, and there "
   echo "*** should be no amplitude/phase jumps across windows."
   smauvplt device=/xs vis=$vis_wide nxy=3,3 select="source($gain),-auto" axis=time,phase
   echo -n "*** HIT RETURN TO CONTINUE ***"
   set ans = "$<"
endif

# Apply passband - wideband
set mir_pb_wide = "astro_pb_wide.mir"
set out = $mir_pb_wide
rm -rf $out
echo ""
echo "*** Applying passband to wideband data  (vis=$vis_wide out=$out)"
rm -rf $vis_wide.bp.tmp
uvcat vis=$vis_wide out=$vis_wide.bp.tmp
gpcopy vis=$vis_wide out=$vis_wide.bp.tmp options=nopol,nocal
uvcat vis=$vis_wide.bp.tmp out=$out options=nocal
set vis_wide = $out


####################################
### FLUX CALIBRATION #########
# To determine antenna gains
####################################
fluxcal:
echo ""
echo "*** Starting flux calibration"

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
    echo "*** Phase-only selfcal on flux calibrator  (vis=$out_primary out=$out interval=$int_fluxgain_phs)"
    selfcal vis=$out_primary refant=$ant interval=$int_fluxgain_phs options=phase,apriori,noscale
    uvcat vis=$out_primary out=$out
    set out_primary = $out


    #plot amp vs uvdist using smauvamp
    echo ""
    echo "*** Plotting amplitudes vs. uvdistance for $flux_planet  (vis=$out_primary)"
    smauvamp device=/xs vis=$out_primary axis=uvd,amp options=zero nbin=20
    echo -n "*** HIT RETURN TO CONTINUE ***"
    set ans = "$<"


    # inspect phase (gain) calibrator
    echo ""
    echo "*** Plotting amplitudes, phas vs. time for $gain GAIN calibrator  (vis=$out_gaincal)"
    echo " *** Inspect Gain calibrator"
    uvplt vis=$out_gaincal device=/xs axis=time,phase "select=win($win)"
    uvplt vis=$out_gaincal device=/xs axis=time,amp "select=win($win)"

    #Phase only selfcal on gain calibrator
    set out = $out_gaincal:r_ph.mir
    rm -rf $out
    echo ""
    echo "*** Phase-only selfcal on gain calibrator  (vis=$out_gaincal out=$out interval=$int_fluxgain_phs)"
    echo " Note we may want interval to be greater than $int_fluxgain_phs"
    echo "Maybe $interval_gain <- use for gain cal amplitude? "
    selfcal vis=$out_gaincal interval=$int_fluxgain_phs refant=$ant options=apriori,noscale,phase
    echo " *** Inspect gain calibrator, assuming we trust the flux"
    echo " Phase gain should be smooth"
    gpplt vis=$out_gaincal device=/xs yaxis=phase nxy=5,3
    echo -n "*** HIT RETURN TO CONTINUE ***"
    set ans = "$<"
    echo " Inspect, every chan should have 0 phase"
    uvspec vis=$out_gaincal axis=chan,phase line="channel,200,1,10,1" device=/xs interval=999 yrange=-180,180
    echo -n "*** HIT RETURN TO CONTINUE ***"
    set ans = "$<"
    echo " Time series of self-calibrated wideband channels"
    uvplt vis=$out_gaincal axis=time,phase device=/xs line="channel,1,16,15,1"
    echo -n "*** HIT RETURN TO CONTINUE ***"
    set ans = "$<"
    echo " *** applying gain tables"
    rm -rf $out_gaincal.tmp
    uvcat vis=$out_gaincal out=$out_gaincal.tmp
    gpcopy vis=$out_gaincal out=$out_gaincal.tmp options=nopass,nopol
    uvcat vis=$out_gaincal.tmp out=$out
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
    # badres: Used to drop baselines on which planet flux is resolved out.
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
echo ""
echo "*** Checking time variance of the phase calibrator"
gplist vis=$vis_wide select="source($gain)" options=zeropha,amp > $vis_wide.gains
cat $vis_wide.gains
echo -n "*** HIT RETURN TO CONTINUE ***"
set ans = "$<"

echo "Sanity Check"
echo " if gain amp > 1 == phase calibrator was weaker than expected, maybe pointing issues"
gpplt vis=$vis_wide select="source($gain)" device=/xs yaxis=amp nxy=5,3
echo -n "*** HIT RETURN TO CONTINUE ***"
set ans = "$<"



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

#######################################################
### Check calibration with point source (second cal) ######
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
    set map   = "$second.cont.temp.map"
    set beam  = "$second.cont.temp.beam"
    set model = "$second.cont.temp.cc"
    set cm    = "$second.cont.temp.cm"
    set res   = "$second.cont.temp.res"
    rm -rf $map $beam $model $res $cm
    invert vis=$vis_wide map=$map beam=$beam robust=$robust cell=$cell options=mosaic,mfs,systemp select="source($second)" imsize=$imsize
    cgdisp device=/xs in=$map labtyp=arcsec options=full,beambr,wedge,3val region=quarter labtyp=hms,dms csize=0.5,0.5,0.5
    echo -n "*** HIT RETURN TO CONTINUE ***"
    set ans = "$<"
endif

#######################################################
### Image gain calibrator
#######################################################
if ($gain != $flux_object) then
    set vis_wide = "astro_pb_wide.mir"
    # Continuum image - MFS
    echo ""
    echo "*** Making continum image of $gain  (vis=$vis_wide)"
    echo "*** Units are in Janskys"
    echo ""
    echo "*** Should look like a Point source "
    set map   = "$gain.cont.temp.map"
    set beam  = "$gain.cont.temp.beam"
    set model = "$gain.cont.temp.cc"
    set cm    = "$gain.cont.temp.cm"
    set res   = "$gain.cont.temp.res"
    rm -rf $map $beam $model $res $cm
    invert vis=$vis_wide map=$map beam=$beam robust=$robust cell=$cell options=mosaic,mfs,systemp select="source($gain)" imsize=$imsize

    cgdisp device=/xs in=$map labtyp=arcsec options=full,beambr,wedge,3val region=quarter labtyp=hms,dms csize=0.5,0.5,0.5
    echo -n "*** HIT RETURN TO CONTINUE ***"
    set ans = "$<"
endif


# *******************
# **** COPY DATA ****
# *******************
# Copy calibrated data to reduced directory
# Add GSV keyword so that invert will take into account relative antenna gains.
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


