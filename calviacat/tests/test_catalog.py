# Licensed with the MIT License, see LICENSE for details

import pytest
from astropy.coordinates import SkyCoord

from ..catalog import Catalog
from ..table import TableDefinition


class TestCatalog(Catalog):
    def _fetch_field(self, sources):
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
