We thank the referee for providing constructive and thoughtful
comments to improve the manuscript. We have addressed all the concerns
raised by the referee and details can be found in the revised version of the manuscript -- please see attached. A brief description on how we incorporated the referee's suggestions can be found below.


1a) We have changed the title to "Molecular gas kinematics of the strongly lensed quasar host galaxy RXS\,J1131-1231".
We agree that we study more than just the CO gas of the science target in this paper; however, a large fraction of the paper is based on our new CO data and the analysis (incl. lens modeling). Therefore, we think that the title is suitable for the paper. 

The literature/archival data of the continuum would not be as meaningful without our CO data, e.g. their emission does not convey inforamtion about the kinematics, spatial distribution, and does not give a broader picture regarding the f_gas evolution, which is one of the main conclusions found in this study. 

1b) We took this suggestion. To minimize potential confusion, we have changed "dynamically resolved" to "kinematically resolved".

1c) We took this suggestion and have updated our redshift (and re-plotted the figure) based on the results from fitting a double-Gaussian to the de-lensed intrinsic line profile of our target (Fig. 5). 

1d) We took this suggestion and have changed "dynamical information" to "kinematic information", "reconstruct the intrinsic gas dynamics" to "reconstruct the intrinsic line profile and source-plane velocity structure".

We have rephrased the 2nd paragraph of Sec 4.1.2 to clarify what we meant.



2a) We took this suggestion and have corrected nu_LO --> nu_obs.

2b) We took this suggestion and re-fitted a doubled Gaussian to the de-lensed profile (Fig. 5) and updated the redshift and replotted all the relavent figures (and updated their captions correspondingly), updated column 1 in Table 3, the labeling in Fig 7, refitted the model shown in Fig. 9 and updated parameters in Sec 4.3 correspondingly, SED fitting, etc.


3) 
We found a bug in generating the map used for previous line-integrated flux extraction. In the resubmitted version, we have updated the intensity reported for CO(2-1) by summing the fluxes in the relavent channels, which is in agreement with the double-Gaussian fit to the spatially extracted spectrum (Fig. 1). Since our lensing-corrected M_gas for both RXJ1131 and the companion are derived from the flux in each channel, this doesn't change our conclusion.

We tried adding the companion's spectrum to plot (Fig. 5) but found the figure to be over-crowded/too busy, as there are already showing various things with different various linestyles, colors. This main message for this plot is to show that we find that the intrinsic profile of the main galaxy is symmetric, which is informative about the kinematics of this galaxy. On the other hand, the useful info for the companion from the de-lensed spectrum is the gas mass, which we report in the text (and in the table, as suggested by the referee).

Source size: We fitted a single Guassian in uv-plane to support the reported CO size based on image plane fitting, of which the convolved size is 6.6+/-0.5 x 4.5+/-0.3. We additionally report the source sizes by fitting two Gaussians to provide further evidence that the CO(2-1) emission is resolved. We have updated Section 3.1 to clarify the discussion.

CO physical size:
The rough estimate of de-magnified CO diameter provided in the referee report is flawed since it assumes the magnification factor to be constant along the arc. The fact that 
1) the lensed emission is not uniformly stretched along the arc (and contains complex structures within the arc)
2) CO is immune to dust extinction, therefore likely to be more extended than the optical emission, as supported by our deconvolved CO source size and in Fig. 7 (where the blue wing centroid is further out than the lensed arc), and from the results of our lens models 
indicates the CO source size is unlikely to be R ≤ 3.5 kpc.

The reported R_CO = 6 kpc is obtained based on the most extended kinematic component is consistent, see Sec 4.3.


4) While the most extreme peak to peak noise is ~100 mJy/beam, the noise varies across channels, and therefore, on average, sigma is lower across all channels. Thus, the uncertainty is more representative using those obtained from the integrated map. 

To demonstrate that our S/N are real regardless of the fainter blue wing, we have independently obtained SNR map to confirm our results, where the peak flux in the red peak is 11.46+/-1.95 Jy/B and the blue peak is 5.13 +/- 1.45 Jy/B. 
I think what was confusing was the reported I_CO(3-2)=35.7+/-21.3Jy km/s, for which I missed a sqrt, and it should have been 35.7+/-6.91 Jy km/s, which has been corrected. 


Since both CO lines are observed towards the same source, using prior knowledge about the blue wing from the CO(2-1) data to inform about the CO(3-2) emission is reasonable. 


While the additioanl figure in I_CO is not significant, we keep this 
to enable reproducible results and to minimize cumulative errors (e.g. when compiled for a larger number of samples). We only report significant figures when reporting the physical parameters.


5) 
We updated the 2mm flux in Table 1. We updated the flux uncertainty to include those from the uv-fitting process.

The peak emission in the residual map is 0.39 +/- 0.08 mJy. While the uncertainty associated with the flux measurement of the background source needs to propagate the errors from subtracting the fitted model flux, the significance of the emission in the residual doesn't not. It simply is Data-Model, which is > 4 sigma than the noise. This together with the co-spatial emission supports that flux contributiion from the background source to the total integrated 2mm continuum.

The term "deblending" is a broadly used in photometry, and it is clear that it refers to the flux subtraction within the same paragraph. We therefore make no changes to the wording.



6) While the large uncertainty on the flux density of the background source at 2mm results in a spectral index with a large uncertainty, this nevertheless is a constraint on the radio-mm spectral slope: our value of -0.345 is flatter than the typical spectral index for pure synchrotron emission, suggesting that the 2mm has a non-negligible contribution from the underlying free-free emission. We have also added an additional statement stating that if there is a non-negligible contribution from the background source contributing to the peak flux, then the derived spectral index for the background source would be shallower, and therefore increasing the fraction of free-free emission over synchrotron at 2mm. This would also explain the atypical flat spectrum found for the foreground galaxy.

The flatter than pure synchrotron spectral slope found for the background source indicates that the 2mm continuum emission is not dominated by synchrotron. Assuming a synchrotron index of alpha = 0.7, we find 0.12 mJy at 2mm for the background source. To account for it, we have updated our SED fitting and the resulting parameters throughout the paper using S = 0.27 mJy at 2mm.

We note that 0.27 mJy is consistent with our previous flux measurement used for SED fitting (0.39+/-0.23 mJy) within the unceratinties, where this uncertainty included in the fit included those from absolute flux calibration. Therefore changes to the physical parameters are only at most a few percent,  we have updated the numbers but this does not alter the conclusions of this paper.


7) We adopted the suggestions.


8) We adopted the suggestions and have added additional tables summarizing the observed CO line properties and the derived physical parameters.


9) It is of standard practice to report flux measurmenets from interferometric observations excluding those from the absolute flux calibration. To address the referee's concern, we have added a remark in Section 3.3, table 1, and Sec 4.5 para 1 to avoid possible confusion for readers.

We have fixed the number of decimals in spectral indices in Sec 3.3.


10) We have replotted the figure with different contour intervals with thinner lines and added a sentence in the caption to improve clarify of Fig 4.

11) We have corrected "blueshifted" to "redshifted", udpated the velocity offset from 780 km/s to 715 km/s using z = 0.65406.

12) 
a) We have removed the statement regarding "AGN dynamically offset".

Our CO measurement confirms that an offset between BLR and the systemic redshift, which is not known previously. We have updated this section, adding additional possible interpretation for the offset (e.g. gas inflows, geometric effects) in the discussion. 

We agree that it is unlikely that all objects with offset AGN lines are black hole recoils. Here, we are not asserting that the offset arises due to a recoiling black hole, we only discuss the possibility of such interpretation given the data (showing evidence of large offset between BLR and the systemic redshift as traced by our CO data), the high spin parameter, and theoretical expectation from hierarchical model of galaxy evolution, among other possible scenarios, as listed by the referee. 


We would greatly appreciate the referee if they can point us to a reference where it shows evidence of why the recoiling black hole scneario is highly unlikely and will gladly make the appropriate changes.


3) Assuming that effects (e.g., inflow, geometric effect, micro-lensing) are accountable for part of the velocity offset, then the velocity offset of ~715 km/s would then be at most an upper limit on the line-of-sight velocity. Then  the true recoil velocity is thus most likely to be << 1000 km/s. Ignoring dynamical friction against stellar background and assuming t_recoil ~ 100 Myr, then the crude BH displacement estimate given in the referee report (~few 10kpc) would be possible (as pointed out to be theoretically possible 
by Bonning+07, Loeb07). However, Loeb07 also pointed out that at low ejection speeds (few hundred km/s), the deceleration of the BH (damped) by the external gravitational potential of the galaxy and by dynamical friction (Chandrasekhar 1943; Apj 97, 255) on the background stars, dark matter, and gas would limit the displacement of the slow-moving BH and tend to bring it back to the center of its host galaxy within a few dynamical times (see also Madau & Quataert 04 for analytic results and Boylan-Kolchin+04 for numerical results). Therefore the displacement given in the referee report is unlikely to hold (see also Guedes+11).

Finally, again, we acknowledge that offsets between BLR and NLR lines are quite common in Seyferts and quasars, but typically not as big of an offset. In this paper, we only discuss possible scenarios and do not assert that any one of mentioned scenario is necessary in RXJ1131.



b & c) The details of the model (e.g. how we reconstruct the likelihood functions) are already referenced in other papers. We did not discuss any of that in this paper. 

However, description of our MCMC sampling has been retained given that conveying the priors and stating how we ensure convergence is necessary for reproducible results and are crucial information for readers to decide on how trust-worthy our inference are from a statistical viewpoint. Explicitly stating the priors is important to demonstrate that our choice are reasonable and physically motivated since they could affect the marginal posteriors strongly and thus the inferred physical parameters. 


13) We updated f_gas and f_mol using M_gas and M_dyn of both RXJ1131+companion , i.e., as a merger system. 

For other comments, please refer to our response to point 11 and 12 (since 
they are repeated here).

14) 
Indeed, the plot on the left is to demonstrate that despite the observed 1st moment map is distorted due to a combination of differential lensing (in which the red velocity are more heavily weighed) and beam smearing effect, our dynamical lens modeling approach, which finds the best-fit models by taking into account the beam convolution indirectly, has enabled us to reconstruct the source-plane kienamtic components, as shown atop of the observed map. The reconstructed information displays a clear velocity gradient. Taken together with previous results from literature and the intrinsically symmetric line profile, we proceed the next section with dynamical modeling. 

To address the referee's concern, we have added in a few lines in the caption to clarify what is being shown in the left panel.

PV: 
As discussed in Sec 4.1.1, we have detailed the evidence and reasons for which we interpret RXJ1131 as a disk galaxy. Under this assumption, the markers on the right plot are based on extracting along a fitted major axis of a disk (see first paragraphy of Sec 4.3). 

To address the referee's concern and to improve clarify, we have re-iterate in the beginning of the first paragraph of Sec 4.3 that we interpret RXJ1131 as disk galaxy.


Size:
As explained in para. 3 of sec 4.3, we use the maximum observed velocity (i.e.  corresponding to the redmost channel along the disk major axis) at 6+/-3kpc as a proxy to the rotation velocity. Thus the quoted CO disk radius corresponds to the extent of this velocity component, as measured from the line center. 


Clumps:
With the evidence and reasons listed in the begining of Sec 4.3, the reconstructed source-plane distribution is unlikely to be just two kinematically-ordered clumps rather than a single rotating disk.


