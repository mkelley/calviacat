# Licensed with the MIT License, see LICENSE for details

__all__ = ["RefCat2"]

import os
import random
import string

try:
    from mastcasjobs import MastCasJobs
except ImportError:
    MastCasJobs = None

from .catalog import Catalog
from .table import TableDefinition


class CasJobsFailed(Exception):
    pass


# column names and SQLite type
COLUMN_DEFS = (
    ("objid", "INTEGER PRIMARY KEY"),
    ("RA", "FLOAT"),
    ("Dec", "FLOAT"),
    ("Gaia", "FLOAT"),
    ("dGaia", "FLOAT"),
    ("BP", "FLOAT"),
    ("dBP", "FLOAT"),
    ("RP", "FLOAT"),
    ("dRP", "FLOAT"),
    ("g", "FLOAT"),
    ("dg", "FLOAT"),
    ("r", "FLOAT"),
    ("dr", "FLOAT"),
    ("i", "FLOAT"),
    ("di", "FLOAT"),
    ("z", "FLOAT"),
    ("dz", "FLOAT"),
    ("J", "FLOAT"),
    ("dJ", "FLOAT"),
    ("H", "FLOAT"),
    ("dH", "FLOAT"),
    ("K", "FLOAT"),
    ("dK", "FLOAT"),
)

filter2col = {}
for f in ["BP", "RP", "g", "r", "i", "z", "J", "H", "K"]:
    filter2col[f] = {"mag": f, "err": "d{}".format(f)}
filter2col["G"] = {"mag": "Gaia", "err": "dGaia"}

# change table name and objid here and at calviacat.data.index
table_def = TableDefinition("refcat2", COLUMN_DEFS, "objid", "RA", "Dec", filter2col)

del filter2col


class RefCat2(Catalog):
    """ATLAS-RefCat2 catalog queries.

    Catalog includes Gaia G, BP, RP; g, r, i, z; and 2MASS J, H, K.

    Note, the Gaia G magnitudes are labeled "Gaia".


    Parameters
    ----------

    dbfile : string
        Database file name.

    wsid : int, optional
        CasJobs WSID.  If not given, then the environment variable
        ``CASJOBS_WSID`` will be used.

    password : string, optional
        CasJobs password.  If not given, then the environment variable
        ``CASJOBS_PW`` will be used.

    **kwargs
        Any `Catalog` keyword arguments.

    """

    def __init__(self, dbfile, wsid=None, password=None, max_records=5000, **kwargs):
        super().__init__(dbfile, table_def, max_records=max_records, **kwargs)

        self.credentials = {}
        self.credentials["userid"] = (
            os.environ.get("CASJOBS_WSID") if wsid is None else wsid
        )
        self.credentials["password"] = (
            os.environ.get("CASJOBS_PW") if password is None else password
        )

    def _fetch_field(self, ra: float, dec: float, sr: float) -> None:
        if MastCasJobs is None:
            raise ImportError(
                "Attemping online query, but mastcasjobs is not installed"
            )

        name = "calviacat_{}".format(
            "".join(random.choices(string.ascii_uppercase, k=5))
        )

        self.logger.debug(
            (
                "Fetching ATLAS-RefCat2 catalog from MAST over {:.2g} field-of-view."
            ).format(sr)
        )

        q = """
        SELECT TOP {max} {columns}
        INTO MyDB.{name}
        FROM fGetNearbyObjEq({ra},{dec},{sr}) as n
        INNER JOIN refcat2 AS r ON (n.objid=r.objid)
        WHERE r.dr < 0.1
        ORDER BY n.distance
        """.format(
            max=self.max_records,
            columns="r.{}".format(",r.".join(self.table.columns)),
            name=name,
            ra=ra,
            dec=dec,
            sr=sr,
        )

        job = MastCasJobs(context="HLSP_ATLAS_REFCAT2", **self.credentials)
        jobid = job.submit(
            q, task_name=("calviacat refcat2 cone search {:.5f} {:.5f}".format(ra, dec))
        )
        status = job.monitor(jobid)
        if status[0] in (3, 4):
            raise CasJobsFailed("status={}, {}".format(status[0], status[1]))
        tab = job.get_table(name, format="CSV")
        job.drop_table_if_exists(name)

        self.logger.debug(
            "Updating {} with {} sources.".format(self.table.name, len(tab))
        )

        self.db.executemany(
            """
        INSERT OR IGNORE INTO {}
          VALUES({})
        """.format(self.table.name, ",".join("?" * len(self.table.columns))),
            tab,
        )
        self.db.commit()
