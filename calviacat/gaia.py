# Licensed with the MIT License, see LICENSE for details

__all__ = [
    'Gaia'
]

import io
import requests
import numpy as np
import astropy.units as u
from astropy.io import votable
from astropy.table import Column
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
    ('astrometric_excess_noise_sig', 'FLOAT'),
    ('_gmag', 'FLOAT'),
    ('_err_gmag', 'FLOAT'),
    ('_rmag', 'FLOAT'),
    ('_err_rmag', 'FLOAT'),
    ('_imag', 'FLOAT'),
    ('_err_imag', 'FLOAT'),
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
        Use this Gaia data release name. (One of 'dr2', 'edr3', dr3'; 'dr2' is
        assumed if not specified)

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
                          'dr3' : 'gaiadr3.gaia_source',
                        }
        filter2col = {}
        for f in ['G', 'bp', 'rp']:
            filter2col[f] = {
                'mag': f'phot_{f.lower():}_mean_mag',
                'err': f'phot_{f.lower():}_mean_flux_over_error'
            }
        # Add synthesized SDSS mags
        for f in ['g', 'r', 'i']:
            filter2col[f] = {
                'mag': f'_{f:}mag',
                'err': f'_err_{f:}mag'
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
            columns=','.join([col for col in self.table.columns if col[0] != '_']),
            ra=np.mean(sources.ra.deg),
            dec=np.mean(sources.dec.deg),
            sr=sr.deg
        )
        # self.logger.debug(q)

        gaia = TapPlus(url='https://gea.esac.esa.int/', server_context='tap-server', tap_context='tap')
        job = gaia.launch_job(q)
        tab = job.get_results()

        # Compute mag. errors from fluxes
        tab = self._add_mag_errors(tab)

        # Transform Gaia magnitudes into SDSS g, r, i ones
        tab = self._transform_filters(tab)
        self.logger.debug('Updating {} with {} sources.'.format(
            self.table.name, len(tab)))

        self.db.executemany('''
        INSERT OR IGNORE INTO {}
          VALUES({})
        '''.format(self.table.name, ','.join('?' * len(self.table.columns))),
            self._masked_to_null(tab))
        self.db.commit()

    def _add_mag_errors(self, tab):
        """Adds magnitude errors to Gaia table <tab> and returns it
        This assumes a symmetric error distribution which won't be correct at
        very low flux levels but you shouldn't be using those for calibration.
        """

        FLUX2MAG = 2.5/np.log(10)
        for filt in ['G', 'bp', 'rp']:
            col_name = self.table.filter2col[filt]['err']
            tab[col_name] = FLUX2MAG/tab[col_name]
        return tab

    def _transform_filters(self, tab):
        """Transform the Gaia G, G_BP, G_RP magnitudes to SDSS g, r, i
        using the equations in Table 5.7 of (DR2)
        https://gea.esac.esa.int/archive/documentation/GDR2/Data_processing/chap_cu5pho/sec_cu5pho_calibr/ssec_cu5pho_PhotTransf.html
        and Table 5.6 of (EDR3)
        https://gea.esac.esa.int/archive/documentation/GEDR3/Data_processing/chap_cu5pho/cu5pho_sec_photSystem/cu5pho_ssec_photRelations.html
        and Table 5.8 of (DR3;appears identical, may not have been updated?)
        https://gea.esac.esa.int/archive/documentation/GDR3/Data_processing/chap_cu5pho/cu5pho_sec_photSystem/cu5pho_ssec_photRelations.html
        """

        g = tab[self.table.filter2col['G']['mag']].data
        e_g = tab[self.table.filter2col['G']['err']].data
        bp = tab[self.table.filter2col['bp']['mag']].data
        rp = tab[self.table.filter2col['rp']['mag']].data

        if self.dr == 'dr2':
            g_sdss = g - (0.13518 - 0.46245*(bp-rp) -
                          0.25171*(bp-rp)**2 + 0.021349*(bp-rp)**3)
            e_g_sdss = np.sqrt(e_g**2 + 0.16497**2)
            r_sdss = g - (-0.12879 + 0.24662*(bp-rp) -
                          0.027464*(bp-rp)**2 - 0.049465*(bp-rp)**3)
            e_r_sdss = np.sqrt(e_g**2 + 0.066739**2)
            i_sdss = g - (-0.29676 + 0.64728*(bp-rp) - 0.10141*(bp-rp)**2)
            e_i_sdss = np.sqrt(e_g**2 + 0.098957**2)
        elif self.dr == 'edr3' or self.dr == 'dr3':
            g_sdss = g - (0.2199 - 0.6365*(bp-rp) -
                          0.1548*(bp-rp)**2 + 0.0064*(bp-rp)**3)
            e_g_sdss = np.sqrt(e_g**2 + 0.0745**2)
            r_sdss = g - (-0.09837 + 0.08592*(bp-rp) +
                          0.1907*(bp-rp)**2 - 0.1701*(bp-rp)**3 +
                          0.02263*(bp-rp)**4)
            e_r_sdss = np.sqrt(e_g**2 + 0.03776**2)
            i_sdss = g - (-0.293 + 0.6404*(bp-rp) -
                          0.09609*(bp-rp)**2 - 0.002104*(bp-rp)**3)
            e_i_sdss = np.sqrt(e_g**2 + 0.04092**2)
        else:
            raise ValueError("Transformations not defined yet for other than DR2 and EDR3/DR3")

        tab.add_column(Column(data=g_sdss, name='_gmag',
                                    unit=u.mag))
        tab.add_column(Column(data=e_g_sdss, name='_err_gmag',
                                    unit=u.mag))
        tab.add_column(Column(data=r_sdss, name='_rmag',
                                    unit=u.mag))
        tab.add_column(Column(data=e_r_sdss, name='_err_rmag',
                                    unit=u.mag))
        tab.add_column(Column(data=i_sdss, name='_imag',
                                    unit=u.mag))
        tab.add_column(Column(data=e_i_sdss, name='_err_imag',
                                    unit=u.mag))
        return tab

    @staticmethod
    def _masked_to_null(tab):
        """Replace masked values with ``None``."""
        for row in tab:
            yield [None if val is np.ma.masked else val
                   for val in row]
