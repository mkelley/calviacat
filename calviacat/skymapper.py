# Licensed with the MIT License, see LICENSE for details

__all__ = [
    'SkyMapper'
]

import io
import requests
import numpy as np
from astropy.io import votable
from astroquery.utils.tap.core import Tap
from .catalog import Catalog, TableDefinition

# column names and SQLite type
COLUMN_DEFS = (
    ('object_id', 'INTEGER PRIMARY KEY'),
    ('raj2000', 'FLOAT'),
    ('dej2000', 'FLOAT'),
    ('e_raj2000', 'FLOAT'),
    ('e_dej2000', 'FLOAT'),
    ('smss_j', 'TEXT'),
    ('flags', 'INTEGER'),
    ('ngood', 'INTEGER'),
    ('ngood_min', 'INTEGER'),
    ('u_flags', 'INTEGER'),
    ('u_ngood', 'INTEGER'),
    ('v_flags', 'INTEGER'),
    ('v_ngood', 'INTEGER'),
    ('g_flags', 'INTEGER'),
    ('g_ngood', 'INTEGER'),
    ('r_flags', 'INTEGER'),
    ('r_ngood', 'INTEGER'),
    ('i_flags', 'INTEGER'),
    ('i_ngood', 'INTEGER'),
    ('z_flags', 'INTEGER'),
    ('z_ngood', 'INTEGER'),
    ('class_star', 'FLOAT'),
    ('a', 'FLOAT'),
    ('e_a', 'FLOAT'),
    ('b', 'FLOAT'),
    ('e_b', 'FLOAT'),
    ('u_psf', 'FLOAT'),
    ('e_u_psf', 'FLOAT'),
    ('v_psf', 'FLOAT'),
    ('e_v_psf', 'FLOAT'),
    ('g_psf', 'FLOAT'),
    ('e_g_psf', 'FLOAT'),
    ('r_psf', 'FLOAT'),
    ('e_r_psf', 'FLOAT'),
    ('i_psf', 'FLOAT'),
    ('e_i_psf', 'FLOAT'),
    ('z_psf', 'FLOAT'),
    ('e_z_psf', 'FLOAT'),
    ('prox', 'FLOAT'),
    ('prox_id', 'INTEGER')
)


class SkyMapper(Catalog):
    def __init__(self, dbfile, max_records=2000, **kwargs):
        filter2col = {}
        for f in 'uvgriz':
            filter2col[f] = {
                'mag': f + '_psf',
                'err': 'e_' + f + '_psf'
            }
        skym = TableDefinition('skymapper', COLUMN_DEFS, 'object_id',
                               'raj2000', 'dej2000', filter2col)
        super().__init__(dbfile, skym, max_records=max_records, **kwargs)

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

        self.logger.debug(
            ('Fetching SkyMapper catalog from ASVO over {:.2g}'
             ' field-of-view.').format(sr))

        q = '''
        SELECT TOP {max}
        {columns}
        FROM dr1.master
        WHERE 1=CONTAINS(POINT('ICRS', raj2000, dej2000),
                         CIRCLE('ICRS', {ra}, {dec}, {sr}))
        ORDER BY ngood DESC
        '''.format(
            max=self.max_records,
            columns=','.join(self.table.columns),
            ra=np.mean(sources.ra.deg),
            dec=np.mean(sources.dec.deg),
            sr=sr.deg
        )
        # self.logger.debug(q)

        skym = Tap(url='http://skymappertap.asvo.nci.org.au/ncitap/tap/')
        job = skym.launch_job(q)
        tab = job.get_results()

        self.logger.debug('Updating {} with {} sources.'.format(
            self.table.name, len(tab)))

        self.db.executemany('''
        INSERT OR IGNORE INTO {}
          VALUES({})
        '''.format(self.table.name, ','.join('?' * len(self.table.columns))),
            tab)
        self.db.commit()
