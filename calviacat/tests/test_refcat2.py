# Licensed with the MIT License, see LICENSE for details

import os

import astropy.units as u
import numpy as np
from astropy.coordinates import SkyCoord

from ..data.refcat2 import add
from ..refcat2 import RefCat2


def data_path(filename):
    data_dir = os.path.join(os.path.dirname(__file__), "data")
    return os.path.join(data_dir, filename)


def test_refcat2(tmpdir):
    fn = data_path("mast-atlas-refcat2.csv")
    dbfn = tmpdir.join("test.db")
    add(fn, dbfn)

    rc2 = RefCat2(dbfn, offline=True)

    # first, try to get all the test catalog data
    sources = SkyCoord([0, 0, 90, 180, 270], [-90, -89, -89, -89, -89], unit="deg")
    objids, _ = rc2.search(sources)
    assert len(objids) == 100

    # match two specific sources
    rc2.min_matches = 1
    sources = SkyCoord(
        [179.45364084, 193.41319168], [-89.95472658, -89.95036938], unit="deg"
    )
    objids, distances = rc2.xmatch(sources)
    assert u.allclose(distances, 0)

    # test calibrate to Gaia G and r
    m_inst = [-10.96, -10.67]
    zp, C, _, G, bmr, gmi = rc2.cal_color(objids, m_inst, "G", "BP-RP", mlim=[0, 20])
    assert np.isclose(zp, 30.218, atol=0.001)
    assert np.isclose(C, -0.144, atol=0.001)
    assert np.allclose(G, [19.011, 19.388])
    assert np.allclose(bmr, [1.723, 1.117])
    assert np.allclose(gmi, [1.587, 0.833])

    zp, C, _, r, gmr, gmi = rc2.cal_color(objids, m_inst, "r", "g-r", mlim=[0, 20])
    assert np.isclose(zp, 29.995, atol=0.001)
    assert np.isclose(C, 0.009, atol=0.001)
    assert np.allclose(r, [19.044, 19.330])
    assert np.allclose(gmr, [1.039, 0.590])
