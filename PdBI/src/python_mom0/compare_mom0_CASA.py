'''

Issue: can't plot my mom0 with APLpy...
Compare mom0 I got from python versus using immoments from CASA to make sure the data array is fine and the issue comes from the header.
- compare stats on source region & off source region

Last Modified: 16 May 2016

History:
--------
16 May 2016:
    - created

Note:
-----
    - got similar RMS off source, and flux on source between mine and from CASA
    - CASA immoments chan=125~159 <=> GILDAS 126-160 <=> python 125-160 compare with that from makeMom0.py

'''
from os.path import join
from os import chdir, getcwd

ori_dir = getcwd()
chdir('/Users/admin/Research/RXJ1131/PdBI/data/15May16/')

# make a mom0 straight from CASA
image = 'centralizedCube4GILDAS.image'
rmtables('testmom0.image')
immoments(imagename=image, moments=[0], chans='125~159', outfile='testmom0.image')

# ----------------------------------------------------------
# compare with histograms in makeMom0.py
ia.open('testmom0.image')
data = ia.getregion()
ia.close()
data.shape
rms = np.std(data[15:220, 15:100])
import matplotlib.pylab as plt
import matplotlib.mlab as mlab
from scipy.stats import norm
n, bins, patches = plt.hist(data[15:220, 15:100].flatten(), bins=100, normed=True)
(mu, sigma) = norm.fit(data[15:220, 15:100].flatten())
y = mlab.normpdf(bins, mu, sigma)
l = plt.plot(bins, y, 'r--', linewidth=2)
print mu, sigma
plt.show()
# ----------------------------------------------------------

new_mom0 = 'centralizedCube4GILDAS-python_ch126-160_mom0.fits'
importfits(fitsimage=new_mom0, imagename=new_mom0.replace('.fits', '.image'))

rmtables(new_mom0.replace('.fits', '.image'))
imstat(imagename=new_mom0.replace('.fits', '.image'), region='offsource_region.crtf')
imstat(imagename='testmom0.image', region='offsource_region.crtf')
# --> rms and std same as above, but flux values are different.

# check on source region
imstat(imagename=new_mom0.replace('.fits', '.image'), region='Flux_CO32_20.3Jy_region.txt')
imstat(imagename='testmom0.image', region='Flux_CO32_20.3Jy_region.txt')
# --> overall consistent.., so for the most part I believe they are consistent, except for somehow off source region have some offset.


# -------------------------
#   Export, for APLpy
# --------------------------
exportfits(imagename='testmom0.image', fitsimage=new_mom0.replace('python', 'CASA'), overwrite=True)

# chdir(ori_dir)
