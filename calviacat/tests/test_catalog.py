# Licensed with the MIT License, see LICENSE for details

import astropy.units as u
import pytest
from astropy.coordinates import SkyCoord

from ..exceptions import XMatchError
from .mock_catalog import MockCatalog


def test_fetch_field_center_and_radius():
    cat = MockCatalog(None)

    # center is a "scalar" coordinate
    center = SkyCoord(1 * u.deg, 2 * u.deg)
    cat.fetch_field(center, 15 * u.arcmin)
    rows = cat.lookup(range(10), "r")
    assert len(rows) == 10

    # center is a list of 1
    center = SkyCoord([1] * u.deg, [2] * u.deg)
    cat.fetch_field(center, 15 * u.arcmin)
    rows = cat.lookup(range(10), "r")
    assert len(rows) == 10


def test_fetch_field_coordinate_list():
    cat = MockCatalog(None)

    # center is a list of coordinates
    center = SkyCoord([1, 2] * u.deg, [2, 3] * u.deg)
    cat.fetch_field(center)
    rows = cat.lookup(range(10), "r")
    assert len(rows) == 10


def test_offline():
    cat = MockCatalog(None, offline=True)
    sources = SkyCoord(
        [179.45364084, 193.41319168], [-89.95472658, -89.95036938], unit="deg"
    )
    with pytest.warns(UserWarning):
        cat.fetch_field(sources)


def test_no_catalog_sources():
    cat = MockCatalog(None)
    sources = SkyCoord(
        [179.45364084, 193.41319168], [-89.95472658, -89.95036938], unit="deg"
    )
    with pytest.raises(XMatchError, match="No catalog sources to match."):
        cat.xmatch(sources)


def test_not_enough_matches():
    cat = MockCatalog(None)
    sources = SkyCoord(
        [179.45364084, 193.41319168], [-89.95472658, -89.95036938], unit="deg"
    )

    cat.db.execute("INSERT INTO test VALUES(1, 179.45364084, -89.95472658, 15.0, 0.01)")
    cat.db.commit()

    with pytest.raises(
        XMatchError, match=f"Fewer than {cat.min_matches} sources matched: 1."
    ):
        cat.xmatch(sources)
