# Licensed with the MIT License, see LICENSE for details

__all__ = ["PanSTARRS1"]

import io

import requests
from astropy.io import votable

from .catalog import Catalog
from .table import TableDefinition

# column names and SQLite type
COLUMN_DEFS = (
    ("objName", "TEXT"),
    ("objID", "INTEGER PRIMARY KEY"),
    ("raMean", "FLOAT"),
    ("decMean", "FLOAT"),
    ("raMeanErr", "FLOAT"),
    ("decMeanErr", "FLOAT"),
    ("nDetections", "INTEGER"),
    ("randomID", "FLOAT"),
    ("projectionID", "INTEGER"),
    ("skyCellID", "INTEGER"),
    ("objInfoFlag", "INTEGER"),
    ("qualityFlag", "INTGER"),
    ("raStack", "TEXT"),
    ("decStack", "TEXT"),
    ("raStackErr", "FLOAT"),
    ("decStackErr", "FLOAT"),
    ("epochMean", "FLOAT"),
    ("nStackDetections", "INTEGER"),
    ("ng", "INTEGER"),
    ("gQfPerfect", "FLOAT"),
    ("gMeanPSFMag", "FLOAT"),
    ("gMeanPSFMagErr", "FLOAT"),
    ("gMeanKronMag", "FLOAT"),
    ("gMeanKronMagErr", "FLOAT"),
    ("gMeanApMag", "FLOAT"),
    ("gMeanApMagErr", "FLOAT"),
    ("gFlags", "INTEGER"),
    ("nr", "INTEGER"),
    ("rQfPerfect", "FLOAT"),
    ("rMeanPSFMag", "FLOAT"),
    ("rMeanPSFMagErr", "FLOAT"),
    ("rMeanKronMag", "FLOAT"),
    ("rMeanKronMagErr", "FLOAT"),
    ("rMeanApMag", "FLOAT"),
    ("rMeanApMagErr", "FLOAT"),
    ("rFlags", "INTEGER"),
    ("ni", "INTEGER"),
    ("iQfPerfect", "FLOAT"),
    ("iMeanPSFMag", "FLOAT"),
    ("iMeanPSFMagErr", "FLOAT"),
    ("iMeanKronMag", "FLOAT"),
    ("iMeanKronMagErr", "FLOAT"),
    ("iMeanApMag", "FLOAT"),
    ("iMeanApMagErr", "FLOAT"),
    ("iFlags", "INTEGER"),
    ("nz", "INTEGER"),
    ("zQfPerfect", "FLOAT"),
    ("zMeanPSFMag", "FLOAT"),
    ("zMeanPSFMagErr", "FLOAT"),
    ("zMeanKronMag", "FLOAT"),
    ("zMeanKronMagErr", "FLOAT"),
    ("zMeanApMag", "FLOAT"),
    ("zMeanApMagErr", "FLOAT"),
    ("zFlags", "INTEGER"),
    ("ny", "INTEGER"),
    ("yQfPerfect", "FLOAT"),
    ("yMeanPSFMag", "FLOAT"),
    ("yMeanPSFMagErr", "FLOAT"),
    ("yMeanKronMag", "FLOAT"),
    ("yMeanKronMagErr", "FLOAT"),
    ("yMeanApMag", "FLOAT"),
    ("yMeanApMagErr", "FLOAT"),
    ("yFlags", "INTEGER"),
)

filter2col = {}
for f in "grizy":
    filter2col[f] = {"mag": f + "meanpsfmag", "err": f + "meanpsfmag" + "err"}

table_def = TableDefinition(
    "panstarrs1", COLUMN_DEFS, "objid", "ramean", "decmean", filter2col
)

del filter2col


class PanSTARRS1(Catalog):
    def __init__(self, dbfile, **kwargs):
        super().__init__(dbfile, table_def, **kwargs)

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
            except Exception:
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
