# Licensed with the MIT License, see LICENSE for details

"""Offline query data may be retrieved from:

https://catalogs.mast.stsci.edu/api/v0.1/panstarrs/dr1/mean?ra=210.802429&dec=54.348750&radius=0.00833333&pagesize=10&format=json&columns=[objName,objID,raMean,decMean,raMeanErr,decMeanErr,nDetections,randomID,projectionID,skyCellID,objInfoFlag,qualityFlag,raStack,decStack,raStackErr,decStackErr,epochMean,nStackDetections,ng,gQfPerfect,gMeanPSFMag,gMeanPSFMagErr,gMeanKronMag,gMeanKronMagErr,gMeanApMag,gMeanApMagErr,gFlags,nr,rQfPerfect,rMeanPSFMag,rMeanPSFMagErr,rMeanKronMag,rMeanKronMagErr,rMeanApMag,rMeanApMagErr,rFlags,ni,iQfPerfect,iMeanPSFMag,iMeanPSFMagErr,iMeanKronMag,iMeanKronMagErr,iMeanApMag,iMeanApMagErr,iFlags,nz,zQfPerfect,zMeanPSFMag,zMeanPSFMagErr,zMeanKronMag,zMeanKronMagErr,zMeanApMag,zMeanApMagErr,zFlags,ny,yQfPerfect,yMeanPSFMag,yMeanPSFMagErr,yMeanKronMag,yMeanKronMagErr,yMeanApMag,yMeanApMagErr,yFlags]

save as data/ps1.json

"""

import os

import astropy.units as u
import numpy as np
import pytest
import requests
from astropy.coordinates import SkyCoord

from calviacat.panstarrs1 import PanSTARRS1


def mocked_requests_get(*args, **kwargs) -> requests.Response:
    """Return saved PS1 request."""
    response = requests.Response()
    with open(f"{os.path.dirname(__file__)}/data/ps1.json") as inf:
        response._content = inf.read(-1).encode()

    return response


class TestPanSTARRS1:
    @pytest.mark.remote_data
    def test_fetch_field_remote(self):
        corners = SkyCoord([1, 1.05, 1.05, 1] * u.deg, [1, 1, 1.05, 1.05] * u.deg)
        ps1 = PanSTARRS1(":memory:")
        ps1.fetch_field(corners)
        catalog_ids, _ = ps1.search(corners)
        assert len(catalog_ids) == 144

    def test_fetch_field(self, monkeypatch):
        monkeypatch.setattr(requests, "get", mocked_requests_get)

        corners = SkyCoord(
            [210.79, 210.81, 210.81, 210.79] * u.deg,
            [54.34, 54.34, 54.36, 54.36] * u.deg,
        )

        ps1 = PanSTARRS1(":memory:")
        ps1.fetch_field(corners)
        catalog_ids, _ = ps1.search(corners)
        assert len(catalog_ids) == 10

        row = ps1.lookup(catalog_ids[:1], ",".join(ps1.table.columns))[0]

        expected = [
            "PSO J140310.964+542031.310",
            173212107956881105,
            210.79565067,
            54.34202381,
            0.03886000066995621,
            0.03886000066995621,
            1,
            0.567481060109549,
            2381,
            53,
            310456320,
            37,
            -999.0,
            -999.0,
            -999.0,
            -999.0,
            56842.2708912,
            0,
            0,
            -999.0,
            -999.0,
            -999.0,
            -999.0,
            -999.0,
            -999.0,
            -999.0,
            0,
            0,
            -999.0,
            -999.0,
            -999.0,
            -999.0,
            -999.0,
            -999.0,
            -999.0,
            0,
            1,
            0.9707229733467102,
            18.699100494384766,
            0.03585200011730194,
            12.125699996948242,
            0.0010189999593421817,
            17.134700775146484,
            0.03585200011730194,
            16777496,
            0,
            -999.0,
            -999.0,
            -999.0,
            -999.0,
            -999.0,
            -999.0,
            -999.0,
            0,
            0,
            -999.0,
            -999.0,
            -999.0,
            -999.0,
            -999.0,
            -999.0,
            -999.0,
            0,
        ]
        assert row[0] == expected[0]
        for i in range(1, len(row)):
            assert np.isclose(row[i], expected[i])
