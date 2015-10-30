from os import path
import requests
import Image, StringIO

import numpy as np
import healpy as hp

from astropy.wcs import WCS

url = 'http://www.fvalk.com/images/MaptoGeo/world-view-total.jpg'

r = requests.get(url)

earth =  np.array(Image.open(StringIO.StringIO(r.content)))

NAXIS2, NAXIS1, colors = earth.shape

# Build a fake header for the supposed Platte Carre projection
wcs = WCS(naxis=2)
wcs.wcs.crpix = [NAXIS1/2, NAXIS2/2]
# -alph to have the same convetion as astro
# -delt to inverse the map (north is up)
wcs.wcs.cdelt = np.array([-360./NAXIS1, -180./NAXIS2])
wcs.wcs.crval = [0, 0]
wcs.wcs.ctype = ["RA---CAR", "DEC--CAR"]
#wcs.wcs.ctype = ["RA---MER", "DEC--MER"] diverge at the pole, so limits must be presents...

# import matplotlib.pyplot as plt
# from wcsaxes import WCS
# plt.ion()
# fig = plt.figure()
# ax = fig.add_axes([0.15, 0.1, 0.8, 0.8], projection=WCS(wcs.to_header()))
# ax.imshow(earth)
# ax.grid(color='white')

# Size of the corresponding healpix map :

# - find out the surface of a typical pixel on the equator in sr
pixel_size = np.radians(np.abs(wcs.wcs.cdelt)).prod()

# healpix goes to 4 pi over the sphere and
# npix = 12 * nside **2 with nside a power of two
# so, if take an healpix resolution just below the original resolution
nside = 2**np.int(np.ceil(np.log2((np.sqrt((4*np.pi / pixel_size)/12)))))


nest = False
npix = hp.nside2npix(nside)
ipix = np.arange(npix)

# Map all healpix pixel to lon and lat
lat, lon = hp.pix2ang(nside, ipix, nest=nest)

# Healpix convention
ra = lon
dec = np.pi/2 - lat

# back to earth maps, nearest neighbor
ind_x, ind_y = wcs.wcs_world2pix(np.degrees(ra), np.degrees(dec), 0)
ind_x = np.ceil(ind_x).astype(np.int)
ind_y = np.ceil(ind_y).astype(np.int)

print earth.shape
print ind_x.min(), ind_x.max()
print ind_y.min(), ind_y.max()

# extract from earth maps to build heapix maps
healpix_maps = [earth[ind_y, ind_x, i] for i in range(3)]

outfilename = path.splitext(path.basename(url))[0] + '_healpix_%s.fits'
for healpix_map, color in zip(healpix_maps, ['r','g','b']):
    hp.write_map(outfilename%color, healpix_map, nest=nest, coord='E')
