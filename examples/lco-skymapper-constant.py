# Licensed with the MIT License, see LICENSE for details

import os

import matplotlib.pyplot as plt
import numpy as np
import requests
from astropy.coordinates import SkyCoord
from astropy.io import fits
from astropy.table import Table

import calviacat as cvc

if os.path.exists("lco.fits"):
    hdu = fits.open("lco.fits")
else:
    r = requests.get("https://archive-api.lco.global/frames/6913550/").json()
    hdu = fits.open(r["url"])
    hdu.writeto("lco.fits")

im = hdu["sci"].data
h = hdu["sci"].header
phot = Table(hdu["cat"].data)

phot = phot[phot["FLAG"] == 0]  # clean LCO catalog
lco = SkyCoord(phot["RA"], phot["DEC"], unit="deg")

# initialize catalog
skym = cvc.SkyMapper("cat.db")

# check the local catalog database for stars in the vicinity of the LCO catalog
catalog_ids, coordinates = skym.search(lco)

# If there are no objects in the local database, then we need to retrieve them
# from PS1's online database.  For these LCO data, finding anything fewer than
# 500 catalog objects suggests there may only be partial overlap between the
# local database and the LCO data.  This number will depend on image depth and
# location on the sky.
if len(catalog_ids) < 500:
    # search the PS1 online database for objects in the vicinity of the LCO
    # catalog
    skym.fetch_field(lco)

# crossmatch LCO photometry table with catalog
objids, distances = skym.xmatch(lco)

# Calibrate this g-band image
g_inst = -2.5 * np.log10(phot["FLUX"])
g_err = phot["FLUXERR"] / phot["FLUX"] * 1.0857

# no color info, so make g-i unlimited
zp_mean, zp_median, unc, g, gmi = skym.cal_constant(objids, g_inst, "g", gmi_lim=None)

# plot results
fig = plt.figure(1)
fig.clear()
ax = fig.gca()
ax.scatter(g, g - g_inst, marker=".", color="k")
ax.axhline(zp_mean)
ax.axhline(zp_median)
plt.setp(ax, xlabel="$g$ (mag)", ylabel=r"$g - g_{\rm inst}$ (mag)")
plt.tight_layout()
plt.savefig("lco-skymapper-constant.png", dpi=150)
