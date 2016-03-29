'''

Continuum is mostly unresolved, can just take peak as the foreground emission.
subtract a point source model as the foreground, then extract the flux for the arc

Note
----
- have attempted this in CASA before in 10Nov15/, which has same problem as here -- uvmodelfit throws an error

'''
import os

uvfiles = ['co21_trial1_cont_UVstyle.uvfits', 'co21_trial1_cont_AIPSstyle.uvfits']
for i in uvfiles:
    importuvfits(fitsfile=i, vis=i.replace('.uvfits', '.ms'))

# copy over antenna tables
os.system('cp -r ' + uvfiles[1].replace('.uvfits', '.ms') + '/ANTENNA/* ' + uvfiles[0].replace('.uvfits', '.ms')+'/ANTENNA/.')
os.system('cp -r ' + uvfiles[1].replace('.uvfits', '.ms') + '/POINTING/* ' + uvfiles[0].replace('.uvfits', '.ms')+'/POINTING/.')
os.system('cp -r ' + uvfiles[1].replace('.uvfits', '.ms') + '/FEED/* ' + uvfiles[0].replace('.uvfits', '.ms')+'/FEED/.')
os.system('cp -r ' + uvfiles[1].replace('.uvfits', '.ms') + '/HISTORY/* ' + uvfiles[0].replace('.uvfits', '.ms')+'/HISTORY/.')

# coord of G (astrometric corrected)
ra, dec = 172.9643456388, -12.532862944
fixvis(vis=uvfiles[0].replace('.uvfits', '.ms'), outputvis=uvfiles[0].replace('.uvfits', '.fixvis.ms'), phasecenter='J2000 11h31m51.44 -12d31m58.3')

vis = 'co21_trial1_cont_UVstyle.fixvis.ms'
plotms(vis=vis, xaxis='uvwave', yaxis='real', ydatacolumn='data')
exportuvfits(vis=vis, fitsfile=vis.replace('.ms', '.uvfits'), datacolumn='data')

newvis = vis.replace('.ms', '.pointModel.ms')
os.system('cp -r ' + vis + ' ' + newvis)
default('uvmodelfit')
comptype = 'P'
sourcepar = [800e-6, 0, 0]
vis = newvis
field = '0'
spw = ''
outfile = newvis.replace('.ms', '.uvmodelfit.cl')
inp()
uvmodelfit()

# 2016-03-28 22:05:07 SEVERE  uvmodelfit::Calibrater::selectvis   Caught exception: Found no data to fit!

#
default('ft')
vis = newvis
complist = newvis.replace('.ms', '.uvmodelfit.cl')
usescratch = True    # ensured model visibilities are stored in the MS
inp()
ft()
plotms(vis=vis, xaxis='uvwave', yaxis='real', ydatacolumn='model')


# subtract the FT of the model (added with ft()), and write output in CORRECTED_DATA
vis = newvis
uvsub()

# Get flux for arc ...
