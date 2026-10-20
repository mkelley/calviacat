# Licensed with the MIT License, see LICENSE for details


import pytest

from ..data.index import create_index, drop_index
from ..exceptions import DatabaseError
from .mock_catalog import MockCatalog


def test_valid_table():
    cat = MockCatalog(None, offline=True)
    drop_index(cat.db, "test")
    create_index(cat.db, "test")

    # OK if index already exists
    create_index(cat.db, "test")


def test_invalid_table():
    cat = MockCatalog(None, offline=True)
    with pytest.raises(LookupError):
        create_index(cat.db, "doesnotexist")


def test_missing_trigger():
    cat = MockCatalog(None)
    cat.db.execute("DROP TRIGGER test_insert")
    with pytest.raises(DatabaseError):
        create_index(cat.db, "test")
