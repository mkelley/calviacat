import os
import numpy as np
import astropy.units as u

__all__ = [
    'Catalog'
]

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

COLUMNS = tuple((c[0] for c in COLUMN_DEFS))


class Catalog:
    def __init__(self, dbfile=None, logger=None, max_records=3000,
                 match_limit=1.5 * u.arcsec, min_matches=10):
        import logging

        if logger is None:
            self.logger = logging.getLogger('calviacat')
            self.logger.setLevel(logging.DEBUG)
        else:
            self.logger = logger

        self._open_db(dbfile)
        self.max_records = max_records
        self.match_limit = match_limit
        self.min_matches = min_matches

    def _open_db(self, dbfile):
        """Open and, if needed, initialize database."""
        import sqlite3

        if dbfile is None:
            self.db = sqlite3.connect(':memory:')
        else:
            self.db = sqlite3.connect(dbfile)

        sqlite3.register_adapter(np.int64, int)
        sqlite3.register_adapter(np.int32, int)
        sqlite3.register_adapter(np.float64, float)
        sqlite3.register_adapter(np.float32, float)

        self.db.create_function('sin', 1, np.sin)

        defs = ',\n  '.join(['{} {}'.format(*c) for c in COLUMN_DEFS])
        self.db.execute('''
        CREATE TABLE IF NOT EXISTS
        panstarrs(
          ''' + defs + '''
        )''')

        self.db.execute('''
        CREATE VIRTUAL TABLE IF NOT EXISTS
        skytree USING RTREE(
          objid INTEGER PRIMARY KEY,
          ra0, ra1, sra0, sra1, dec0, dec1
        )''')

        self.db.execute('''
        CREATE TRIGGER IF NOT EXISTS
        panstarrs_insert AFTER INSERT ON panstarrs
        BEGIN
          INSERT INTO skytree VALUES (
            new.objid,
            new.ramean,new.ramean,
            sin(new.ramean),sin(new.ramean),
            new.decmean,new.decmean
          );
        END''')

    def fetch_field(self, sources):
        """Fetch catalog sources for this field.

        Parameters
        ----------
        sources : SkyCoord
            Sources to be matched.

        """
        import io
        import requests
        from astropy.io import votable

        sr = max((sources.separation(c).max() for c in sources)) * 5 / 8

        self.logger.debug(
            'Fetching PS1 catalog from STScI over {} field-of-view.'.format(
                sr))

        params = dict(RA=sources.ra.mean().deg, DEC=sources.dec.mean().deg,
                      SR=sr.deg, max_records=self.max_records,
                      ordercolumn1='ndetections',
                      descending1='on',
                      selectedColumnsCsv=','.join(COLUMNS))

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

        self.logger.debug('Updating panstarrs with {} sources.'.format(
            len(tab)))

        self.db.executemany('''
        INSERT OR IGNORE INTO panstarrs
          VALUES(''' + ','.join('?' * len(COLUMNS)) + ''')
        ''', tab)
        self.db.commit()

    def search(self, sources, col='rmeanpsfmag'):
        """Search catalog for objects near sources.

        Parameters
        ----------
        sources : SkyCoord
            Sources to be matched.
        col : string
          Column to return.

        Returns
        -------
        objid : array
          Catalog ID.
        cat : SkyCoord
          Catalog coordinates.

        """

        from astropy.coordinates import SkyCoord

        mean_ra = np.mean(sources.ra)
        select = '''
        SELECT panstarrs.objid,ramean,decmean,{} FROM panstarrs
        INNER JOIN skytree ON panstarrs.rowid = skytree.rowid
        '''.format(col)

        if mean_ra < 10 * u.deg:
            # Near RA=0, use sin(ra) for catalog search
            rows = self.db.execute(select + '''
            WHERE skytree.sra0 > ?
              AND skytree.sra1 < ?
              AND skytree.dec0 > ?
              AND skytree.dec1 < ?
            ''', (np.sin(sources.ra.min().rad), np.sin(sources.ra.max().rad),
                  sources.dec.min().deg, sources.dec.max().deg)).fetchall()
        else:
            rows = self.db.execute(select + '''
            WHERE skytree.ra0 > ?
              AND skytree.ra1 < ?
              AND skytree.dec0 > ?
              AND skytree.dec1 < ?
            ''', (sources.ra.min().deg, sources.ra.max().deg,
                  sources.dec.min().deg, sources.dec.max().deg)).fetchall()

        if len(rows) == 0:
            return np.array([]), SkyCoord([], {}, unit='deg'), np.array([])

        objid, ra, dec, c = [np.array(x) for x in zip(*rows)]
        cat = SkyCoord(ra, dec, unit='deg')
        return objid, cat, c

    def xmatch(self, sources, col):
        """Cross-match sources to catalog.

        Parameters
        ----------
        sources : SkyCoord
            Sources to cross-match.
        col : string
            Column to return.

        Returns
        -------
        matched : masked array
            Matched object ID of each source.
        d : masked array
            Distance between match and source in arcsec.
        c : masked array
            The requested column.

        """

        objid, cat, c = self.search(sources, col)

        if len(cat) == 0:
            self.logger.error('No catalog sources to match.')
            return

        idx, d2d = sources.match_to_catalog_sky(cat)[:2]
        i = d2d < self.match_limit
        n = i.sum()
        if n < self.min_matches:
            self.logger.error('Fewer than 10 sources matched to PS1.')
            return

        self.logger.info(
            'Matched {} sources to photometric catalog.'.format(n))

        matched = np.ma.MaskedArray(objid[idx], int)
        matched.mask = ~i

        d = np.ma.MaskedArray(d2d.arcsec, float)
        d.mask = ~i

        c = np.ma.MaskedArray(c[idx], float)
        c.mask = ~i

        return matched, d, c
