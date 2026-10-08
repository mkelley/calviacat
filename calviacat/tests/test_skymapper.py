# Licensed with the MIT License, see LICENSE for details

import astropy.units as u
from astropy.coordinates import SkyCoord

from calviacat.skymapper import SkyMapper


class TestSkyMapper:
    def test_fetch_field(self):
        corners = SkyCoord([1, 1.2, 1.2, 1] * u.deg, [1, 1, 1.2, 1.2] * u.deg)
        skym = SkyMapper(":memory:", dr=2)
        skym.fetch_field(corners)
        catalog_ids, coords = skym.search(corners)
        assert len(catalog_ids) == 225

    def test_dr4(self):
        """DR4 has a slightly different set of columns"""
        corners = SkyCoord([1, 1.2, 1.2, 1] * u.deg, [1, 1, 1.2, 1.2] * u.deg)
        skym = SkyMapper(":memory:", dr=4)
        skym.fetch_field(corners)
        catalog_ids, coords = skym.search(corners)
        assert len(catalog_ids) == 337
