# Licensed with the MIT License, see LICENSE for details

__all__ = ["SkyMapper"]

from copy import copy
import numpy as np
import pyvo as vo
from .catalog import Catalog, TableDefinition

# column names and SQLite type
COLUMN_DEFS_DR12 = (
    ("object_id", "INTEGER PRIMARY KEY"),
    ("raj2000", "FLOAT"),
    ("dej2000", "FLOAT"),
    ("e_raj2000", "FLOAT"),
    ("e_dej2000", "FLOAT"),
    ("smss_j", "TEXT"),
    ("flags", "INTEGER"),
    ("ngood", "INTEGER"),
    ("ngood_min", "INTEGER"),
    ("u_flags", "INTEGER"),
    ("u_ngood", "INTEGER"),
    ("v_flags", "INTEGER"),
    ("v_ngood", "INTEGER"),
    ("g_flags", "INTEGER"),
    ("g_ngood", "INTEGER"),
    ("r_flags", "INTEGER"),
    ("r_ngood", "INTEGER"),
    ("i_flags", "INTEGER"),
    ("i_ngood", "INTEGER"),
    ("z_flags", "INTEGER"),
    ("z_ngood", "INTEGER"),
    ("class_star", "FLOAT"),
    ("u_psf", "FLOAT"),
    ("e_u_psf", "FLOAT"),
    ("v_psf", "FLOAT"),
    ("e_v_psf", "FLOAT"),
    ("g_psf", "FLOAT"),
    ("e_g_psf", "FLOAT"),
    ("r_psf", "FLOAT"),
    ("e_r_psf", "FLOAT"),
    ("i_psf", "FLOAT"),
    ("e_i_psf", "FLOAT"),
    ("z_psf", "FLOAT"),
    ("e_z_psf", "FLOAT"),
    ("prox", "FLOAT"),
    ("prox_id", "INTEGER"),
)

COLUMN_DEFS_DR4 = (
    ("object_id", "INTEGER PRIMARY KEY"),
    ("raj2000", "FLOAT"),
    ("dej2000", "FLOAT"),
    ("e_raj2000", "FLOAT"),
    ("e_dej2000", "FLOAT"),
    ("smss_j", "TEXT"),
    ("flags", "INTEGER"),
    ("ngood", "INTEGER"),
    ("u_flags", "INTEGER"),
    ("u_ngood", "INTEGER"),
    ("v_flags", "INTEGER"),
    ("v_ngood", "INTEGER"),
    ("g_flags", "INTEGER"),
    ("g_ngood", "INTEGER"),
    ("r_flags", "INTEGER"),
    ("r_ngood", "INTEGER"),
    ("i_flags", "INTEGER"),
    ("i_ngood", "INTEGER"),
    ("z_flags", "INTEGER"),
    ("z_ngood", "INTEGER"),
    ("class_star", "FLOAT"),
    ("u_psf", "FLOAT"),
    ("e_u_psf", "FLOAT"),
    ("v_psf", "FLOAT"),
    ("e_v_psf", "FLOAT"),
    ("g_psf", "FLOAT"),
    ("e_g_psf", "FLOAT"),
    ("r_psf", "FLOAT"),
    ("e_r_psf", "FLOAT"),
    ("i_psf", "FLOAT"),
    ("e_i_psf", "FLOAT"),
    ("z_psf", "FLOAT"),
    ("e_z_psf", "FLOAT"),
)


class SkyMapper(Catalog):
    """Calibrate to SkyMapper catalog data.


    Parameters
    ----------
    dbfile : string
        Sqlite3 database file name.  If a database at the provided file name
        does not exist, a new one will be created.

    max_records : int, optional
        Maximum number of records to return from online queries to SkyMapper.

    dr : int
        Use this SkyMapper data release number: 1, 2, or 4.

    logger : Logger
        Use this python logger for logging.

    match_limit : astropy Quantity
        Plane of sky tolerance for catalog matches.

    min_matches : int
        Throw an error if fewer than this many matches are found.

    """

    def __init__(self, dbfile, max_records=2000, dr=4, **kwargs):
        filter2col = {}
        for f in "uvgriz":
            filter2col[f] = {"mag": f + "_psf", "err": "e_" + f + "_psf"}

        if dr in [1, 2]:
            column_defs = COLUMN_DEFS_DR12
        elif dr == 4:
            column_defs = COLUMN_DEFS_DR4
        else:
            raise ValueError(
                f"SkyMapper data releases 1, 2, and 4 are supported (DR{dr} requested)"
            )

        skym = TableDefinition(
            f"skymapper_dr{dr}",
            column_defs,
            "object_id",
            "raj2000",
            "dej2000",
            filter2col,
        )
        self.dr = dr
        super().__init__(dbfile, skym, max_records=max_records, **kwargs)

    def fetch_field(self, ra: float, dec: float, sr: float) -> None:
        self.logger.debug(
            (
                "Fetching SkyMapper catalog from ASVO over {:.2g}" " field-of-view."
            ).format(sr)
        )

        query = """
        SELECT TOP {max}
        {columns}
        FROM dr{dr}.master
        WHERE 1=CONTAINS(POINT('ICRS', raj2000, dej2000),
                         CIRCLE('ICRS', {ra}, {dec}, {sr}))
        ORDER BY ngood DESC
        """.format(
            dr=self.dr,
            max=self.max_records,
            columns=",".join(self.table.columns),
            ra=ra,
            dec=dec,
            sr=sr,
        )

        skym = vo.dal.TAPService("https://api.skymapper.nci.org.au/public/tap")
        tab = skym.search(query).to_table()

        self.logger.debug(
            "Updating {} with {} sources.".format(self.table.name, len(tab))
        )

        self.db.executemany(
            """
        INSERT OR IGNORE INTO {}
          VALUES({})
        """.format(
                self.table.name, ",".join("?" * len(self.table.columns))
            ),
            self._masked_to_null(tab),
        )
        self.db.commit()

    @staticmethod
    def _masked_to_null(tab):
        """Replace masked values with ``None``."""
        for row in tab:
            yield [None if val is np.ma.masked else val for val in row]
