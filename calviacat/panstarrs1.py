# Licensed with the MIT License, see LICENSE for details

__all__ = [
    'PanSTARRS1'
]

import io
import requests
import numpy as np
from astropy.io import votable
from astropy.coordinates import SkyCoord
import astropy.units as u
from .catalog import Catalog, TableDefinition

# column names and SQLite type
COLUMN_DEFS = (
    ('objname', 'TEXT'),
    ('objid', 'INTEGER PRIMARY KEY'),
    ('ramean', 'FLOAT'),
    ('decmean', 'FLOAT'),
    ('rameanerr', 'FLOAT'),
    ('decmeanerr', 'FLOAT'),
    ('ndetections', 'INTEGER'),
    ('randomid', 'FLOAT'),
    ('projectionid', 'INTEGER'),
    ('skycellid', 'INTEGER'),
    ('objinfoflag', 'INTEGER'),
    ('qualityflag', 'INTGER'),
    ('rastack', 'TEXT'),
    ('decstack', 'TEXT'),
    ('rastackerr', 'FLOAT'),
    ('decstackerr', 'FLOAT'),
    ('epochmean', 'FLOAT'),
    ('nstackdetections', 'INTEGER'),
    ('ng', 'INTEGER'),
    ('gqfperfect', 'FLOAT'),
    ('gmeanpsfmag', 'FLOAT'),
    ('gmeanpsfmagerr', 'FLOAT'),
    ('gmeankronmag', 'FLOAT'),
    ('gmeankronmagerr', 'FLOAT'),
    ('gmeanapmag', 'FLOAT'),
    ('gmeanapmagerr', 'FLOAT'),
    ('gflags', 'INTEGER'),
    ('nr', 'INTEGER'),
    ('rqfperfect', 'FLOAT'),
    ('rmeanpsfmag', 'FLOAT'),
    ('rmeanpsfmagerr', 'FLOAT'),
    ('rmeankronmag', 'FLOAT'),
    ('rmeankronmagerr', 'FLOAT'),
    ('rmeanapmag', 'FLOAT'),
    ('rmeanapmagerr', 'FLOAT'),
    ('rflags', 'INTEGER'),
    ('ni', 'INTEGER'),
    ('iqfperfect', 'FLOAT'),
    ('imeanpsfmag', 'FLOAT'),
    ('imeanpsfmagerr', 'FLOAT'),
    ('imeankronmag', 'FLOAT'),
    ('imeankronmagerr', 'FLOAT'),
    ('imeanapmag', 'FLOAT'),
    ('imeanapmagerr', 'FLOAT'),
    ('iflags', 'INTEGER'),
    ('nz', 'INTEGER'),
    ('zqfperfect', 'FLOAT'),
    ('zmeanpsfmag', 'FLOAT'),
    ('zmeanpsfmagerr', 'FLOAT'),
    ('zmeankronmag', 'FLOAT'),
    ('zmeankronmagerr', 'FLOAT'),
    ('zmeanapmag', 'FLOAT'),
    ('zmeanapmagerr', 'FLOAT'),
    ('zflags', 'INTEGER'),
    ('ny', 'INTEGER'),
    ('yqfperfect', 'FLOAT'),
    ('ymeanpsfmag', 'FLOAT'),
    ('ymeanpsfmagerr', 'FLOAT'),
    ('ymeankronmag', 'FLOAT'),
    ('ymeankronmagerr', 'FLOAT'),
    ('ymeanapmag', 'FLOAT'),
    ('ymeanapmagerr', 'FLOAT'),
    ('yflags', 'INTEGER')
)


class PanSTARRS1(Catalog):
    def __init__(self, dbfile, phot_type='meanpsfmag', **kwargs):
        filter2col = {}
        for f in 'grizy':
            filter2col[f] = {
                'mag': f + phot_type,
                'err': f + phot_type + 'err'
            }
        ps1 = TableDefinition('panstarrs1', COLUMN_DEFS, 'objid',
                              'ramean', 'decmean', filter2col)
        super().__init__(dbfile, ps1, **kwargs)

    def fetch_field(self, sources, scale=1.25, radius=0.5):
        """Fetch catalog sources for this field and save to database.

        Search radius and center are derived from the source list if
        len(sources) > 1 or from the first entry and radius if there is only
        a single source.

        Parameters
        ----------
        sources : SkyCoord or tuple of (RA center, Dec center) as angle Quantities
                  or in degrees.
            Sources to be matched or center to search around.

        scale : float, optional
            Search radius scale factor.

        radius: float or Quantity, optional, default=0.5
            Search radius in degrees or a Quantity that can be transformed
            to degrees (`scale` is not applied)
        """

        if type(sources) == SkyCoord:
            try:
                sr = max((sources.separation(c).max() for c in sources)) * scale / 2
                sr = sr.to(u.deg).value
                ra_center = sources.ra.mean().deg
                dec_center = sources.dec.mean().deg
            except TypeError:
                # single center position given
                ra_center = sources.ra.deg
                dec_center = sources.dec.deg
                try:
                    sr = radius.to(u.deg).value
                except AttributeError:
                    sr = radius
        elif type(sources) == tuple:
            ra_center = sources[0]
            dec_center = sources[1]
            try:
                ra_center = ra_center.to(u.deg).value
                dec_center = dec_center.to(u.deg).value
            except AttributeError:
                self.logger.debug("Assuming coordinates are in degrees already")
                pass
            except u.UnitConversionError:
                self.logger.error(("Could not convert {} to center coordinates".format(sources)))
            try:
                sr = radius.to(u.deg).value
            except AttributeError:
                sr = radius
        self.logger.debug(
            ('Fetching PS1 catalog from STScI over {:.2g} deg field-of-view.')
            .format(sr))

        params = dict(RA=ra_center, DEC=dec_center, SR=sr,
                      max_records=self.max_records,
                      ordercolumn1='ndetections',
                      descending1='on',
                      selectedColumnsCsv=','.join(self.table.columns))

        q = requests.get('https://archive.stsci.edu/panstarrs/search.php',
                         params=params)

        with io.BytesIO(q.text.encode()) as xml:
            try:
                tab = votable.parse_single_table(xml).to_table()
            except Exception as e:
                self.logger.error(q.text)
                return

        tab['objname'] = np.array(
            [str(x.decode()) for x in tab['objname']])
        tab['objid'] = np.array(
            [str(x.decode()) for x in tab['objid']], int)
        tab['rastack'] = np.array(
            [str(x.decode()) for x in tab['rastack']])
        tab['decstack'] = np.array(
            [str(x.decode()) for x in tab['decstack']])

        self.logger.debug('Updating {} with {} sources.'.format(
            self.table.name, len(tab)))

        self.db.executemany('''
        INSERT OR IGNORE INTO {}
          VALUES({})
        '''.format(self.table.name, ','.join('?' * len(self.table.columns))),
            tab)
        self.db.commit()
