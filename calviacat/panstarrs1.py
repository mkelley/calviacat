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
    """Calibrate photometry with the PanSTARRS 1 DR2 catalog.

    Calibration uses the meanPSFMag data columns.


    Parameters
    ----------
    {parameters}

    """.format(parameters=Catalog._inherited_parameters)

    URL_BASE: str = "https://catalogs.mast.stsci.edu/api/v0.1/panstarrs/dr2/mean"

    def __init__(self, dbfile, **kwargs):
        super().__init__(dbfile, table_def, **kwargs)

    def _fetch_field(self, ra: float, dec: float, radius: float) -> None:
        self.logger.debug(
            f"Fetching PS1 DR2 catalog from STScI over {radius:.2f} deg field-of-view."
        )

        params = {
            "ra": ra,
            "dec": dec,
            "radius": radius,
            "max_records": self.max_records,
            "nDetections.gte": 2,
            "sort_by": "nDetections.DESC",
            "columns": "[" + ",".join(self.table.columns) + "]",
            "format": "json",
        }

        q = requests.get(self.URL_BASE, params=params)

        with io.BytesIO(q.text.encode()) as xml:
            breakpoint()
            try:
                tab = votable.parse_single_table(xml).to_table()
            except Exception:
                self.logger.error(q.text)
                return

        self.logger.debug(f"Updating {self.table.name} with {len(tab)} sources.")

        self.db.executemany(
            "INSERT OR IGNORE INTO {} VALUES({})".format(
                self.table.name, ",".join("?" * len(self.table.columns))
            ),
            tab,
        )
        self.db.commit()
