### Copied continuum clean map (lmv-clean & fits) from 04Sep15

#### Determine source size, continuum flux using GILDAS and CASA.
-  frequency:        139.8139 GHz (2.144225 mm)

## UV plane with GILDAS
#### UV_fit result for source size:
- script: /Users/admin/Research/RXJ1131/PdBI/src/postIm_analysis/Continuum.map
- co21_trial1_cont.uvt
- one function, E_GAUSS
	- R.A.        =    -2.11931 (   .18227)  11:31:51.4553
	- Dec.        =    -1.29913 (   .22635) -12:31:58.2991
	- **Flux        =     1.36570 (   .20869)  milliJy**
	- **Major       =     3.03174 (594.09415)**
	- **Minor       =     0.00362 (   .53387)**
	- Pos.Ang.    =   -72.76842 ( 15.39338)
	- --> does not appear that this fits the 3*sigma blob seen in image plane, since the source size along minor axis appears unresolved, but the image seems to suggest it's extended as beam broadens it along that axis.
- tried `uv_fit` (E_GAUSS) after `uv_shift` --> similar results
	- R.A.        =    -0.00383 (   .18160)  11:31:51.4550
	- Dec.        =     0.00111 (   .22602) -12:31:58.2980
	- **Flux        =     1.36260 (   .20820)  milliJy**
	- **Major       =     3.01647 (365.39810)**
	- **Minor       =     0.00587 (   .53198)**
	- Pos.Ang.    =   -72.77988 ( 15.48605)

- try fitting other functions 
	- C_GAUSS (after shifting phase center):
		- R.A.        =    -0.10163 (   .14901)  11:31:51.4484
		- Dec.        =     0.06940 (   .23643) -12:31:58.2297
		- **Flux        =     1.28496 (   .18507)  milliJy**
		- **F.W.H.P.    =     2.14503 (   .40583)**
	- Double E_Gauss (**bad fit**):
		- had to use the UVT before shifting phase center, else `uv_fit` throws an error
		- R.A.        =    -6.68012 (      Inf)  11:31:51.1438
		- Dec.        =    42.77865 (      Inf) -12:31:14.2213
		- Flux        =   324.12138 (      Inf)  Jy      
		- Major       =   616.19074 (      Inf)
		- Minor       =   244.74343 (      Inf)
		- Pos.Ang.    =     6.61935 (      Inf)
		- R.A.        =     8.67984 (      Inf)  11:31:52.1928
		- Dec.        =   -40.77913 (      Inf) -12:32:37.7791
		- Flux        =  -324.12062 (      Inf)  Jy      
		- Major       =   624.19198 (      Inf)
		- Minor       =   238.74247 (      Inf)
		- Pos.Ang.    =    -0.55870 (      Inf)


## UV plane with CASA:
#### Attemped uvmodelfit()
- in GILDAS: 
    - `fits co21_trial1_cont_AIPSstyle.uvfits from co21_trial1_cont.uvt /style AIPS`
    - `fits co21_trial1_cont_UVstyle.uvfits from co21_trial1_cont.uvt /style UV`

- in CASA:
    - `importuvfits()` both uvfits file
    - copied ANTENNA and FEED from co21_trial1_cont_AIPSstyle.ms/ to co21_trial1_cont_UVstyle.ms/

    - `uvmodelfit(vis='co21_trial1_cont_UVstyle.ms', selectdata=False)`
      + throws an **error**... "Caught exception: Found no data to fit!"


## Image plane with GILDAS
- `go nice` --> `go flux`
- 1.2 mJy over ~150 pixels (region slightly larger than 3sigma contours)
- 1.1 mJy over ~155 pixels (differet poly, also slightly larger than 3sigma)
- 1.2 mJy over ~237 pixels (region larger than 3sigma, and larger than above)
- 0.458 mJy over ~29 pixels (along 6 sigma contours)
- 0.136 mJy over ~7 pixels (along 9 sigma contours)
- --> source appears slightly extended
- script: /Users/admin/Research/RXJ1131/PdBI/src/postIm_analysis/Continuum.map


## Image Plane with CASA
- CASA_region1 = CASA region file
- CASA_region2 = CASA region file

#### 2D fit in viewer:
1. logfile: CASA2Dfit_region1.log
	* integrated flux: 1.23 +/- 0.22 mJy
	* Image component size (convolved with beam) ---
       - major axis FWHM:     4.17 +/- 0.53 arcsec
       - minor axis FWHM:     3.34 +/- 0.36 arcsec
       - position angle: 166 +/- 19 deg
	* deconvolved source size (N/A)
2. logfile: CASA2Dfit_region2.log
	* integrated flux: 1.22 +/- 0.24 mJy
	* Image component size (convolved with beam) ---
		- major axis FWHM:     4.16 +/- 0.60 arcsec
		- minor axis FWHM:     3.33 +/- 0.40 arcsec
		- position angle: 165 +/- 22 deg
 	* deconvolved source size (N/A)

#### `imstat()` on co21_trial1_cont.fits
- `importfits(fitsimage='co21_trial1_cont.fits', imagename='co21_trial1_cont')`
- `imstat(region='CASA_region1', logfile='CASAimstat_region1.log')`
- **Integrated line flux: 1.130 mJy**; # of pixels = 277


if the continuum is resolved along one axis, then maybe the continuum is from the foreground, given the position of the emission
-> subtracted a point source model in UV plane (to remove the unresolved emission from the foreground): co21_trial1_cont.uv_shift-res.uvt
- confirmed a low SNR emission from the foreground that conincide with the arc in the pre-model subtracted cont. map
- the peak in that (co21_trial1_cont.uv_shift-res.lmv-clean or equivalently co21_trial1_cont.uv_shift-res.cln.fits) ~ difference between integrated and peak in the pre-model subtracted cont. map
