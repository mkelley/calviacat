# Licensed with the MIT License, see LICENSE for details

import pytest
from astropy.coordinates import SkyCoord

from ..catalog import Catalog
from ..exceptions import XMatchError
from ..table import TableDefinition


class TestCatalog(Catalog):
    def _fetch_field(self, sources):  # pragma: no cover
        pass


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


def test_offline():
    cat = TestCatalog(None, table_def, offline=True)
    sources = SkyCoord(
        [179.45364084, 193.41319168], [-89.95472658, -89.95036938], unit="deg"
    )
    with pytest.warns(UserWarning):
        cat.fetch_field(sources)


def test_no_catalog_sources():
    cat = TestCatalog(None, table_def, offline=True)
    sources = SkyCoord(
        [179.45364084, 193.41319168], [-89.95472658, -89.95036938], unit="deg"
    )
    with pytest.raises(XMatchError, match="No catalog sources to match."):
        cat.xmatch(sources)


def test_not_enough_matches():
    cat = TestCatalog(None, table_def, offline=True)
    sources = SkyCoord(
        [179.45364084, 193.41319168], [-89.95472658, -89.95036938], unit="deg"
    )

    cat.db.execute("INSERT INTO test VALUES(1, 179.45364084, -89.95472658, 15.0, 0.01)")
    cat.db.commit()

    with pytest.raises(
        XMatchError, match=f"Fewer than {cat.min_matches} sources matched: 1."
    ):
        cat.xmatch(sources)
