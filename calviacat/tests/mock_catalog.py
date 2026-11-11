# Licensed with the MIT License, see LICENSE for details

import astropy.units as u
import numpy as np
from astropy.coordinates import SkyCoord

from ..catalog import Catalog
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
