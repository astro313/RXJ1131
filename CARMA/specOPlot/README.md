fitted Gauss to CO32 fixing FWHM --> since it's more of a double Gaussian in principle but we don't have enough SNR to fit double gauss.

Just going to integrate line flux using FWHM as in CO21. 

Use CO21 & coOverlay.eps to guide the range of velocity to include.
-----------

CO (3-2):
May 17 2016:
- Update channel range to be consistent with the velocity range used to get the new CO21 Intensity (124-156)
- I = 35.7 Jy km/s from `go view`
- --> integrated over velocity width ~ 710 km/s
- Channel:     48-68
- velocity:   -441.57, 278.79
- I_sigma: see below

May 04 2016 (Deprecated):
- I = 37.1 Jy km/s from `go view` 
- from channel [49, 69]
- velocity [-430.94, 305.26] km/s
- vertices of polygon (Xi, Yi) = (133,150), (118,138), (122,104), (153, 114), (150,141)
- CASA region file of the above polygon: Flux_CO32_37.2Jy_region.txt
- pts area within region = 1159 px
- # of pix in a beam = 112.867 PX
- --> # of beams ~ 10
- sigma = 13.3mJy /B/channel of 35.87 km/s per channel
- Isigma ~ 13.3 mJy/B/channel * 21 channels * 10 beams ~ 2.793 Jy km/s
