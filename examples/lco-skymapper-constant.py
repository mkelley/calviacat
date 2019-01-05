# Licensed with the MIT License, see LICENSE for details

import os
import requests
import numpy as np
import matplotlib.pyplot as plt
from astropy.table import Table
from astropy.io import fits
from astropy.coordinates import SkyCoord
import calviacat as cvc

if os.path.exists('lco.fits'):
    hdu = fits.open('lco.fits')
else:
    r = requests.get('https://archive-api.lco.global/frames/6913550/').json()
    hdu = fits.open(r['url'])
    hdu.writeto('lco.fits')

im = hdu['sci'].data
h = hdu['sci'].header
phot = Table(hdu['cat'].data)

phot = phot[phot['FLAG'] == 0]  # clean LCO catalog
lco = SkyCoord(phot['RA'], phot['DEC'], unit='deg')

# initialize catalog
skym = cvc.SkyMapper('cat.db')

# download catalog?
if len(skym.search(lco)[0]) < 500:
    skym.fetch_field(lco)

# crossmatch LCO photometry table with catalog
objids, distances = skym.xmatch(lco)

# Calibrate this g-band image
g_inst = -2.5 * np.log10(phot['FLUX'])
g_err = phot['FLUXERR'] / phot['FLUX'] * 1.0857

# no color info, so make g-i unlimited
zp_mean, zp_median, unc, g, gmi = skym.cal_constant(
    objids, g_inst, 'g', gmi_lim=None)

# plot results
fig = plt.figure(1)
fig.clear()
ax = fig.gca()
ax.scatter(g, g - g_inst, marker='.', color='k')
ax.axhline(zp_mean)
ax.axhline(zp_median)
plt.setp(ax, xlabel='$g$ (mag)', ylabel=r'$g - g_{\rm inst}$ (mag)')
plt.tight_layout()
plt.savefig('lco-skymapper-constant.png', dpi=150)
