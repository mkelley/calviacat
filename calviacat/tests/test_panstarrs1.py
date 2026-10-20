# Licensed with the MIT License, see LICENSE for details

"""Offline query data may be retrieved from:

https://catalogs.mast.stsci.edu/api/v0.1/panstarrs/dr1/mean?ra=200&dec=50&radius=0.01&nDetections.gte=2&pagesize=10&format=json&columns=[objName,objID,raMean,decMean,raMeanErr,decMeanErr,nDetections,randomID,projectionID,skyCellID,objInfoFlag,qualityFlag,raStack,decStack,raStackErr,decStackErr,epochMean,nStackDetections,ng,gQfPerfect,gMeanPSFMag,gMeanPSFMagErr,gMeanKronMag,gMeanKronMagErr,gMeanApMag,gMeanApMagErr,gFlags,nr,rQfPerfect,rMeanPSFMag,rMeanPSFMagErr,rMeanKronMag,rMeanKronMagErr,rMeanApMag,rMeanApMagErr,rFlags,ni,iQfPerfect,iMeanPSFMag,iMeanPSFMagErr,iMeanKronMag,iMeanKronMagErr,iMeanApMag,iMeanApMagErr,iFlags,nz,zQfPerfect,zMeanPSFMag,zMeanPSFMagErr,zMeanKronMag,zMeanKronMagErr,zMeanApMag,zMeanApMagErr,zFlags,ny,yQfPerfect,yMeanPSFMag,yMeanPSFMagErr,yMeanKronMag,yMeanKronMagErr,yMeanApMag,yMeanApMagErr,yFlags]

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

        sources = SkyCoord(
            [199.9, 200.1, 200.1, 199.9] * u.deg,
            [49.9, 50.1, 50.1, 49.9] * u.deg,
        )

        ps1 = PanSTARRS1(":memory:")
        ps1.fetch_field(sources)
        catalog_ids, _ = ps1.search(sources)
        assert len(catalog_ids) == 10

        row = ps1.lookup([168002000088058883], ",".join(ps1.table.columns))[0]

        expected = [
            "PSO J132002.116+500024.642",
            168002000088058883,
            200.00889097,
            50.00682584,
            0.013030000030994415,
            0.014449999667704105,
            63,
            0.0899922261852699,
            2322,
            51,
            444915712,
            53,
            200.0088184,
            50.00684503,
            0.0010000000474974513,
            0.0010000000474974513,
            56023.74878472,
            5,
            14,
            0.9997140169143677,
            21.28660011291504,
            0.027884000912308693,
            21.096200942993164,
            0.058371998369693756,
            21.015100479125977,
            0.09455599635839462,
            16892216,
            18,
            0.9996370077133179,
            20.04360008239746,
            0.023608999326825142,
            19.75149917602539,
            0.01755400002002716,
            19.702600479125977,
            0.014662000350654125,
            16892216,
            10,
            0.9997000098228455,
            19.59149932861328,
            0.028703000396490097,
            19.299100875854492,
            0.021524999290704727,
            19.26289939880371,
            0.019453000277280807,
            16892216,
            11,
            0.9996650218963623,
            19.3262996673584,
            0.030858999118208885,
            19.109699249267578,
            0.04419099912047386,
            19.038000106811523,
            0.02351599931716919,
            16892216,
            10,
            0.9991440176963806,
            19.010400772094727,
            0.072113998234272,
            18.770299911499023,
            0.045180998742580414,
            18.678800582885742,
            0.024421000853180885,
            16892216,
        ]

        assert row[0] == expected[0]
        for i in range(1, len(row)):
            assert np.isclose(row[i], expected[i])
