# Licensed with the MIT License, see LICENSE for details

__all__ = [
    'Gaia'
]

import io
import requests
import numpy as np
from astropy.io import votable
from astroquery.utils.tap.core import TapPlus
from .catalog import Catalog, TableDefinition

# column names and SQLite type
COLUMN_DEFS = (
    ('source_id', 'INTEGER PRIMARY KEY'),
    ('ra', 'FLOAT'),
    ('dec', 'FLOAT'),
    ('ra_error', 'FLOAT'),
    ('dec_error', 'FLOAT'),
    ('designation', 'TEXT'),
    ('phot_variable_flag', 'TEXT'),
    ('astrometric_n_good_obs_al', 'INTEGER'),
    ('phot_g_n_obs', 'INTEGER'),
    ('phot_bp_n_obs', 'INTEGER'),
    ('phot_rp_n_obs', 'INTEGER'),
    ('phot_g_mean_mag', 'FLOAT'),
    ('phot_g_mean_flux_over_error', 'FLOAT'),
    ('phot_bp_mean_mag', 'FLOAT'),
    ('phot_bp_mean_flux_over_error', 'FLOAT'),
    ('phot_rp_mean_mag', 'FLOAT'),
    ('phot_rp_mean_flux_over_error', 'FLOAT'),
    ('astrometric_gof_al', 'FLOAT'),
    ('astrometric_excess_noise', 'FLOAT'),
    ('astrometric_excess_noise_sig', 'FLOAT')
)


class Gaia(Catalog):
    """Calibrate to Gaia catalog data.


    Parameters
    ----------
    dbfile : string
        Sqlite3 database file name.  If a database at the provided file name
        does not exist, a new one will be created.

    max_records : int, optional
        Maximum number of records to return from online queries to SkyMapper.

    dr : string
        Use this Gaia data release name.

    logger : Logger
        Use this python logger for logging.

    match_limit : astropy Quantity
        Plane of sky tolerance for catalog matches.

    min_matches : int
        Throw an error if fewer than this many matches are found.    

    """

    def __init__(self, dbfile, max_records=2000, dr='dr2', **kwargs):
        table_mapping = { 'dr2' : 'gaiadr2.gaia_source',
                          'edr3' : 'gaiaedr3.gaia_source',
                        }
        filter2col = {}
        for f in ['g', 'bp', 'rp']:
            filter2col[f] = {
                'mag': f'phot_{f:}_mean_mag',
                'err': f'phot_{f:}_flux_over_error'
            }
        gaia = TableDefinition('gaia', COLUMN_DEFS, 'source_id',
                               'ra', 'dec', filter2col)
        self.dr = dr
        self.gaia_table = table_mapping.get(dr.lower, table_mapping['dr2'])
        super().__init__(dbfile, gaia, max_records=max_records, **kwargs)

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
            ('Fetching Gaia {} catalog from ESA over {:.2g}'
             ' field-of-view.').format(self.dr.upper(), sr))

        q = '''
        SELECT TOP {max}
        {columns}
        FROM {table}
        WHERE 1=CONTAINS(POINT('ICRS', ra, dec),
                         CIRCLE('ICRS', {ra}, {dec}, {sr}))
        ORDER BY astrometric_n_good_obs_al DESC
        '''.format(
            table=self.gaia_table,
            max=self.max_records,
            columns=','.join(self.table.columns),
            ra=np.mean(sources.ra.deg),
            dec=np.mean(sources.dec.deg),
            sr=sr.deg
        )
        # self.logger.debug(q)

        gaia = TapPlus(url='https://gea.esac.esa.int/', server_context='tap-server', tap_context='tap')
        job = gaia.launch_job(q)
        tab = job.get_results()

        self.logger.debug('Updating {} with {} sources.'.format(
            self.table.name, len(tab)))

        self.db.executemany('''
        INSERT OR IGNORE INTO {}
          VALUES({})
        '''.format(self.table.name, ','.join('?' * len(self.table.columns))),
            self._masked_to_null(tab))
        self.db.commit()

    @staticmethod
    def _masked_to_null(tab):
        """Replace masked values with ``None``."""
        for row in tab:
            yield [None if val is np.ma.masked else val
                   for val in row]
