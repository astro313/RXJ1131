copied sup127_155_2ndcln_noCont.spec from 30Apr16/; But the central velocity is not shifted to line center.. (spectrum extracted from cube: sup127_155_2ndcln_noCont.lmv-clean in 30Apr16/) 
--> thus, the velocity range corresponding to chan 126-130 in this spectrum file != velocity range after shifted 0km/s to line center (i.e. centralizedCube4GILDAS.lmv-clean in 30Apr16/ & channel map).
*If the spectrum were written out in channel numbers instead of velocity range: it would work fine (since channel numbers didn't change after shifting the line center to the ref. velocity, only the velocity change).*
--> to obtain the intrinsic flux, I need to sum up fluxes from a known velocity range consistent with the channel map, I need to first write out a spectrum where its velocity col. is consistent with the channel map range. 

from 30Apr16/README.md: 
    note that we shifted 0 km/s from chan 130 to chan 145
    i.e. velocity is shifted by ~ 15 * 21.5 km/s

--> shift the velocity spectrum as well
transformVelo.map: write out ASCII spectrum with velocity shifted to ref. velocity (i.e. as if extracted from centralizedCube4GILDAS.lmv-clean)
output ASCII spectrum: veloLineCenterShifted.spec


For plotting spectrum, we would use the velocity scale w/o shifting --> GILDAS script will output axis w.r.t. given z. (sup127_155_2ndcln_noCont.spec, de-lensed.spec)
Note: the velocity label in channel map is also the shifted ones in z=0 frame. 

de-lensed.spec:
- line profile for RXJ1131 only
- for chan 131-160: directly divide by mag. factor
- for chan 126-130: decompose flux contribution from companion and correct for lensing using mag. factor for RXJ1131
- for other channels: copied from lensed spectrum so that Gauss fit can be done; but is of course much noiser compared to the intrinsic fluxes in chan (126-160).

singleFit_delensed.map:
- plot the de-lensed Spectrum for RXJ1131 only and fit single gauss to get FWHM

singleFit_lensed.map:
- fit single Gauss to get FWHM
    - does it differ from that of fit w/o de-lensed? (singleFit_delensed.map)
        + --> problem is: the fit is completely dominated by the redshifted emission.

