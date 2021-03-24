# Licensed with the MIT License, see LICENSE for details

__all__ = [
    'RefCat2'
]

import io
import random
import string
import requests
import numpy as np
from astropy.io import votable
from mastcasjobs import MastCasJobs
from .catalog import Catalog, TableDefinition


class CasJobsFailed(Exception):
    pass


# column names and SQLite type
COLUMN_DEFS = (
    ('objid', 'INTEGER PRIMARY KEY'),
    ('RA', 'FLOAT'),
    ('Dec', 'FLOAT'),
    ('g', 'FLOAT'),
    ('dg', 'FLOAT'),
    ('gchi', 'FLOAT'),
    ('gcontrib', 'INTEGER'),
    ('r', 'FLOAT'),
    ('dr', 'FLOAT'),
    ('rchi', 'FLOAT'),
    ('rcontrib', 'INTEGER'),
    ('i', 'FLOAT'),
    ('di', 'FLOAT'),
    ('ichi', 'FLOAT'),
    ('icontrib', 'INTEGER'),
    ('z', 'FLOAT'),
    ('dz', 'FLOAT'),
    ('zchi', 'FLOAT'),
    ('zcontrib', 'INTEGER'),
    ('nstat', 'INTEGER'),
    ('J', 'FLOAT'),
    ('dJ', 'FLOAT'),
    ('H', 'FLOAT'),
    ('dH', 'FLOAT'),
    ('K', 'FLOAT'),
    ('dK', 'FLOAT')
)


class RefCat2(Catalog):
    """ATLAS-RefCat2 catalog queries.


    Parameters
    ----------

    dbfile : string
        Database file name.

    wsid : int, optional
        CasJobs WSID.  If not given, then the enviornment variable
        ``CASJOBS_WSID`` will be used.

    password : string, optional
        CasJobs password.  If not given, then the enviornment variable
        ``CASJOBS_PW`` will be used.

    **kwargs
        Any `Catalog` keyword arguments.

    """

    def __init__(self, dbfile, wsid=None, password=None, max_records=5000,
                 **kwargs):
        filter2col = {}
        for f in 'griz':
            filter2col[f] = {
                'mag': f,
                'err': 'd{}'.format(f)
            }
        skym = TableDefinition('refcat2', COLUMN_DEFS, 'objid',
                               'RA', 'Dec', filter2col)
        super().__init__(dbfile, skym, max_records=max_records, **kwargs)
        self.credentials = {}
        self.credentials['userid'] = wsid
        self.credentials['password'] = password

    def fetch_field(self, sources, scale=1.25):
        """Fetch catalog sources for this field and save to database.

        Search radius and center are derived from the source list.

        Parameters
        ----------
        sources : SkyCoord
            Sources to be matched.

        scale : float, optional
            Search radius scale factor.

        """
        sr = max((sources.separation(c).max() for c in sources)) * scale / 2

        name = 'calviacat_{}'.format(
            ''.join(random.choices(string.ascii_uppercase, k=5)))

        self.logger.debug(
            ('Fetching ATLAS-RefCat2 catalog from MAST over {:.2g}'
             ' field-of-view.').format(sr))

        q = '''
        SELECT TOP {max} {columns}
        INTO MyDB.{name}
        FROM fGetNearbyObjEq({ra},{dec},{sr}) as n
        INNER JOIN refcat2 AS r ON (n.objid=r.objid)
        WHERE r.dr < 0.1
        ORDER BY n.distance
        '''.format(
            max=self.max_records,
            columns='r.{}'.format(',r.'.join(self.table.columns)),
            name=name,
            ra=np.mean(sources.ra.deg),
            dec=np.mean(sources.dec.deg),
            sr=sr.deg
        )

        job = MastCasJobs(context="HLSP_ATLAS_REFCAT2", **self.credentials)
        jobid = job.submit(q, task_name=('calviacat refcat2 cone search {:.5f} {:.5f}'
                                         .format(np.mean(sources.ra.deg),
                                                 np.mean(sources.dec.deg))))
        status = job.monitor(jobid)
        if status[0] in (3, 4):
            raise CasJobsFailed('status={}, {}'.format(status[0], status[1]))
        tab = job.get_table(name, format='CSV')
        job.drop_table_if_exists(name)

        self.logger.debug('Updating {} with {} sources.'.format(
            self.table.name, len(tab)))

        self.db.executemany('''
        INSERT OR IGNORE INTO {}
          VALUES({})
        '''.format(self.table.name, ','.join('?' * len(self.table.columns))),
            tab)
        self.db.commit()
