# Licensed with the MIT License, see LICENSE for details

__all__ = ["PanSTARRS1"]

import io
import requests
import numpy as np
from astropy.io import votable
from astropy.coordinates import SkyCoord
import astropy.units as u
from .catalog import Catalog, TableDefinition

# column names and SQLite type
COLUMN_DEFS = (
    ("objname", "TEXT"),
    ("objid", "INTEGER PRIMARY KEY"),
    ("ramean", "FLOAT"),
    ("decmean", "FLOAT"),
    ("rameanerr", "FLOAT"),
    ("decmeanerr", "FLOAT"),
    ("ndetections", "INTEGER"),
    ("randomid", "FLOAT"),
    ("projectionid", "INTEGER"),
    ("skycellid", "INTEGER"),
    ("objinfoflag", "INTEGER"),
    ("qualityflag", "INTGER"),
    ("rastack", "TEXT"),
    ("decstack", "TEXT"),
    ("rastackerr", "FLOAT"),
    ("decstackerr", "FLOAT"),
    ("epochmean", "FLOAT"),
    ("nstackdetections", "INTEGER"),
    ("ng", "INTEGER"),
    ("gqfperfect", "FLOAT"),
    ("gmeanpsfmag", "FLOAT"),
    ("gmeanpsfmagerr", "FLOAT"),
    ("gmeankronmag", "FLOAT"),
    ("gmeankronmagerr", "FLOAT"),
    ("gmeanapmag", "FLOAT"),
    ("gmeanapmagerr", "FLOAT"),
    ("gflags", "INTEGER"),
    ("nr", "INTEGER"),
    ("rqfperfect", "FLOAT"),
    ("rmeanpsfmag", "FLOAT"),
    ("rmeanpsfmagerr", "FLOAT"),
    ("rmeankronmag", "FLOAT"),
    ("rmeankronmagerr", "FLOAT"),
    ("rmeanapmag", "FLOAT"),
    ("rmeanapmagerr", "FLOAT"),
    ("rflags", "INTEGER"),
    ("ni", "INTEGER"),
    ("iqfperfect", "FLOAT"),
    ("imeanpsfmag", "FLOAT"),
    ("imeanpsfmagerr", "FLOAT"),
    ("imeankronmag", "FLOAT"),
    ("imeankronmagerr", "FLOAT"),
    ("imeanapmag", "FLOAT"),
    ("imeanapmagerr", "FLOAT"),
    ("iflags", "INTEGER"),
    ("nz", "INTEGER"),
    ("zqfperfect", "FLOAT"),
    ("zmeanpsfmag", "FLOAT"),
    ("zmeanpsfmagerr", "FLOAT"),
    ("zmeankronmag", "FLOAT"),
    ("zmeankronmagerr", "FLOAT"),
    ("zmeanapmag", "FLOAT"),
    ("zmeanapmagerr", "FLOAT"),
    ("zflags", "INTEGER"),
    ("ny", "INTEGER"),
    ("yqfperfect", "FLOAT"),
    ("ymeanpsfmag", "FLOAT"),
    ("ymeanpsfmagerr", "FLOAT"),
    ("ymeankronmag", "FLOAT"),
    ("ymeankronmagerr", "FLOAT"),
    ("ymeanapmag", "FLOAT"),
    ("ymeanapmagerr", "FLOAT"),
    ("yflags", "INTEGER"),
)


class PanSTARRS1(Catalog):
    def __init__(self, dbfile, phot_type="meanpsfmag", **kwargs):
        filter2col = {}
        for f in "grizy":
            filter2col[f] = {"mag": f + phot_type, "err": f + phot_type + "err"}
        ps1 = TableDefinition(
            "panstarrs1", COLUMN_DEFS, "objid", "ramean", "decmean", filter2col
        )
        super().__init__(dbfile, ps1, **kwargs)

    def _fetch_field(self, ra: float, dec: float, sr: float) -> None:
        self.logger.debug(
            f"Fetching PS1 catalog from STScI over {sr:.2f} deg field-of-view."
        )

        params = {
            "RA": ra,
            "DEC": dec,
            "SR": sr,
            "max_records": self.max_records,
            "ordercolumn1": "ndetections",
            "descending1": "on",
            "selectedColumnsCsv": ",".join(self.table.columns),
        }

        q = requests.get(
            "https://archive.stsci.edu/panstarrs/search.php", params=params
        )

        with io.BytesIO(q.text.encode()) as xml:
            try:
                tab = votable.parse_single_table(xml).to_table()
            except Exception as e:
                self.logger.error(q.text)
                return

        self.logger.debug(f"Updating {self.table.name} with {len(tab)} sources.")

        self.db.executemany(
            """
        INSERT OR IGNORE INTO {}
          VALUES({})
        """.format(self.table.name, ",".join("?" * len(self.table.columns))),
            tab,
        )
        self.db.commit()
