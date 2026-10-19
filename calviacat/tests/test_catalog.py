# Licensed with the MIT License, see LICENSE for details

import astropy.units as u
import numpy as np
import pytest
from astropy.coordinates import SkyCoord

from ..catalog import Catalog
from ..exceptions import XMatchError
from ..table import TableDefinition

COLUMN_DEFS = (
    ("source_id", "INTEGER PRIMARY KEY"),
    ("ra", "FLOAT"),
    ("dec", "FLOAT"),
    ("r", "FLOAT"),
    ("rerr", "FLOAT"),
)
table_def = TableDefinition(
    "test", COLUMN_DEFS, "source_id", "ra", "dec", {"r": {"mag": "r", "err": "rerr"}}
)


class MockCatalog(Catalog):
    def __init__(self, dbfile, **kwargs):
        super().__init__(dbfile, table_def, **kwargs)

    def _fetch_field(self, ra, dec, sr):  # pragma: no cover
        np.random.seed(23)
        n = 10
        coord = SkyCoord(ra * u.deg, dec * u.deg)
        pa = np.random.rand(n) * 360 * u.deg
        offset = np.random.rand(n) * sr * u.deg
        coords = coord.directional_offset_by(pa, offset)
        r = 10 + np.arange(n)
        rerr = np.ones(n) * 0.03

        self.db.executemany(
            "INSERT OR IGNORE INTO {} VALUES({})".format(
                self.table.name, ",".join("?" * len(self.table.columns))
            ),
            zip(range(n), coords.ra.value, coords.dec.value, r, rerr),
        )
        self.db.commit()


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
