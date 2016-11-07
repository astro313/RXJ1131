# to get SNR separately for red wind and blue wing in repsonse to referee's point 4

# in paper, we used chan 48-68 to form 0th moment <-> FWZI and not FWHM
# spectrum extracted from (133,150), (118,138), (122,104), (153, 114), (150,141)
# blue wing ~ channel 48-56
# red ~ channel 56-69

echo ""
echo "*** Moment 0 map of CO(3-2) in the inner quarter map, no continuum subtraction"
echo ""

# stat before doing anything
histo in=RXJ1131.LSB_linewithCont.cm region="relcenter,boxes(-100,-100,100,-30)"

# blue wing
rm -rf LSB_withCont.co32_quarter_chan48_56.cm
imsub in=RXJ1131.LSB_linewithCont.cm region="quarter(48,56)" out=LSB_withCont.co32_quarter_chan48_56.cm

# noise
imlist in=LSB_withCont.co32_quarter_chan48_56.cm
histo in=LSB_withCont.co32_quarter_chan48_56.cm region="relcenter,boxes(-50,-50,50,-30)"
# rms ~ 1.26E-02 Jy/B per channel

rm -rf LSB_withCont.co32_quarter_chan48_56.mom0
moment mom=0 in=LSB_withCont.co32_quarter_chan48_56.cm out=LSB_withCont.co32_quarter_chan48_56.mom0

# noise in mom0
histo in=LSB_withCont.co32_quarter_chan48_56.mom0 region="relcenter,boxes(-50,-50,50,-30)"
# rms = 1.43 Jy/B
histo in=LSB_withCont.co32_quarter_chan48_56.mom0
# max = 5.13 Jy

cgdisp in=LSB_withCont.co32_quarter_chan48_56.mom0,LSB_withCont.co32_quarter_chan48_56.mom0 type=p,c xybin=1,1 device=/xs nxy=1,1 options=full,beambr,wedge,trlab,3val labtyp=arcsec,arcsec csize=1,1,1 range=0,5.5,log,2 slev=a,1.4 levs1=-3,-1,1,2,3,4,5
# highest: 3 sigma contours

fits in=LSB_withCont.co32_quarter_chan48_56.mom0 out=LSB_withCont.co32_quarter_chan48_56.mom0.fits op=xyout


# Red wing
rm -rf LSB_withCont.co32_quarter_chan56_69.cm
imsub in=RXJ1131.LSB_linewithCont.cm region="quarter(56,69)" out=LSB_withCont.co32_quarter_chan56_69.cm

# noise
imlist in=LSB_withCont.co32_quarter_chan56_69.cm
histo in=LSB_withCont.co32_quarter_chan56_69.cm region="relcenter,boxes(-50,-50,50,-30)"
# rms ~ 1.19e-2 Jy/B per channel

rm -rf LSB_withCont.co32_quarter_chan56_69.mom0
moment mom=0 in=LSB_withCont.co32_quarter_chan56_69.cm out=LSB_withCont.co32_quarter_chan56_69.mom0

# noise in mom0
histo in=LSB_withCont.co32_quarter_chan56_69.mom0 region="relcenter,boxes(-50,-50,50,-30)"
# rms = 1.99 Jy/B
histo in=LSB_withCont.co32_quarter_chan56_69.mom0
# max = 11.45

cgdisp in=LSB_withCont.co32_quarter_chan56_69.mom0,LSB_withCont.co32_quarter_chan56_69.mom0 type=p,c xybin=1,1 device=/xs nxy=1,1 options=full,beambr,wedge,trlab,3val labtyp=arcsec,arcsec csize=1,1,1 range=0,12,log,2 slev=a,2 levs1=-3,-1,1,3,4,5,6,7,8,9
# --> max: 5 sigma

fits in=LSB_withCont.co32_quarter_chan56_69.mom0 out=LSB_withCont.co32_quarter_chan56_69.mom0.fits op=xyout

#
# Next, look at it with GILDAS --> how many sigmas? # let noise, go nice
#
######################################################################

# GILDAS:
# Blue wing:
# let name LSB_withCont.co32_quarter_chan48_56.mom0
# let type fits
# let size 35
# go nice

# poly
# mean
# --> rms: 1.45 Jy/B
# let spacing 1.45

# extrema /compute
# --> max: 5.13 Jy/B
# - --> go nice show max 3*sigma contour, where sigma = 2 Jy km/s /B
# go nice
#


# Red wing:
# let name LSB_withCont.co32_quarter_chan56_69.mom0
# let type fits
# let size 35
# go nice

# poly
# mean
# --> rms: 1.94 Jy/B

# extrema /compute
# --> max: 11.455 Jy/B
# - --> go nice show max 5*sigma contour, where sigma = 2 Jy km/s /B
# consistency check: Zmax / noise ~ 5*sigma
