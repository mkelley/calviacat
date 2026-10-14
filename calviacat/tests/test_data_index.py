# Licensed with the MIT License, see LICENSE for details

from sqlite3 import OperationalError

import pytest

from ..catalog import Catalog
from ..data.index import create_index, drop_index
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


def test_valid_table():
    cat = TestCatalog(None, table_def, offline=True)
    drop_index(cat.db, "test")
    create_index(cat.db, "test")

    # index already exists
    with pytest.raises(OperationalError):
        create_index(cat.db, "test")


def test_invalid_table():
    cat = TestCatalog(None, table_def, offline=True)
    with pytest.raises(LookupError):
        create_index(cat.db, "doesnotexist")
