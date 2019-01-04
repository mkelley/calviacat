# calviacat
Calibrate star photometry by comparison to a catalog.  PanSTARRS 1 catalog currently supported.

## Requirements

* sqlite3
* requests
* astropy

## Example

Calibrate a Las Cumbres Observatory image of 41P/Tuttle-Giacobini-Kresak:

``` python
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
    r = requests.get('https://archive-api.lco.global/frames/6143031/').json()
    hdu = fits.open(r['url'])
    hdu.writeto('lco.fits')

im = hdu['sci'].data
h = hdu['sci'].header
phot = Table(hdu['cat'].data)

phot = phot[phot['FLAG'] == 0]  # clean LCO catalog
lco = SkyCoord(phot['RA'], phot['DEC'], unit='deg')

# initialize catalog
ps1 = cvc.PanSTARRS1('cat.db')

# download PS1 catalog?
if len(ps1.search(lco)[0]) < 500:
    ps1.fetch_field(lco)

# crossmatch LCO photometry table with catalog
objids, distances = ps1.xmatch(lco)

# Calibrate this g-band image, include a color correction
g_inst = -2.5 * np.log10(phot['FLUX'])
g_err = phot['FLUXERR'] / phot['FLUX'] * 1.0857

zp, C, unc, g, gmr = ps1.cal_color(objids, g_inst, 'g', 'g-r')

# plot results
fig = plt.figure(1)
fig.clear()
ax = fig.gca()
ax.scatter(gmr, g - g_inst, marker='.', color='k')
x = np.linspace(0, 1.5)
ax.plot(x, C * x + zp, 'r-')
plt.setp(ax, xlabel='$g-r$ (mag)', ylabel='$g$ (mag)')
plt.tight_layout()

```

![alt text](lco-example.png "Best-fit calibration")
