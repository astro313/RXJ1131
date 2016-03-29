#!/usr/bin/env csh
#
# Attempted to subtract a point source model in the UV plane with MIRIAD to get flux from arc
#
# Note:
# - didn't work, error in uvmodel: "Warning [uvmodel]:  Error accessing pixel data of co21_trial1_cont_UVstyle.fixvis.model, in XYOPEN"
#
#
##########################################################################
#

fits op=uvin in=co21_trial1_cont_UVstyle.fixvis.uvfits out=co21_trial1_cont_UVstyle.fixvis.mir

uvfit vis=co21_trial1_cont_UVstyle.fixvis.mir object=point spar=8.0E-6,0.,0. out=co21_trial1_cont_UVstyle.fixvis.model

uvmodel vis=co21_trial1_cont_UVstyle.fixvis.mir model=co21_trial1_cont_UVstyle.fixvis.model options=subtract out=co21_trial1_cont_UVstyle.fixvis.residual

# Get remaining flux here, and image to see if see arc

