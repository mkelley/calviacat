# Licensed with the MIT License, see LICENSE for details

__all__ = [
    'Catalog'
]

import logging
import sqlite3
from abc import ABC, abstractmethod

import numpy as np
import astropy.units as u
from astropy.coordinates import SkyCoord
from astropy.stats import sigma_clipped_stats, sigma_clip
from astropy.modeling import models, fitting

try:
    from astropy.version import version_info as astropy_version
except ImportError:
    import astropy.version
    astropy_version = [int(x) for x in astropy.version.version.split('.')]


class CalibrationError(Exception):
    pass


class TableDefinition:
    def __init__(self, name, column_definitions, objid, ra, dec,
                 filter2col):
        self.name = name
        self.column_definitions = column_definitions
        self.objid = objid
        self.ra = ra
        self.dec = dec
        self.filter2col = filter2col

    @property
    def columns(self):
        return tuple((c[0] for c in self.column_definitions))


class Catalog(ABC):
    def __init__(self, dbfile, table, logger=None,
                 max_records=3000, match_limit=1.5 * u.arcsec,
                 min_matches=10):

        if logger is None:
            self.logger = logging.getLogger(__name__)
            self.logger.setLevel(logging.DEBUG)
            if len(self.logger.handlers) == 0:
                self.logger.addHandler(logging.StreamHandler())
        else:
            self.logger = logger

        self.table = table
        self.max_records = max_records
        self.match_limit = match_limit
        self.min_matches = min_matches

        self._open_db(dbfile)

    def _open_db(self, file_name):
        """Open and, if needed, initialize database.

        Parameters
        ----------
        file_name : string
            Name of the data base file, or `None` to use a temporary
            memory file.

        """

        if file_name is None:
            self.db = sqlite3.connect(':memory:')
        else:
            self.db = sqlite3.connect(file_name, timeout=30)

        sqlite3.register_adapter(np.int64, int)
        sqlite3.register_adapter(np.int32, int)
        sqlite3.register_adapter(np.float64, float)
        sqlite3.register_adapter(np.float32, float)

        self.db.create_function('SIN', 1, np.sin)
        self.db.create_function('COS', 1, np.cos)

        defs = ',\n  '.join(['{} {}'.format(*c)
                             for c in self.table.column_definitions])
        self.db.execute('''
        CREATE TABLE IF NOT EXISTS
        {}(
          {}
        )'''.format(self.table.name, defs))

        self.db.execute('''
        CREATE VIRTUAL TABLE IF NOT EXISTS
        {}_skytree USING RTREE(
          {} INTEGER PRIMARY KEY,
          x0 FLOAT, x1 FLOAT,
          y0 FLOAT, y1 FLOAT,
          z0 FLOAT, z1 FLOAT
        )'''.format(self.table.name, self.table.objid))

        self.db.execute('''
        CREATE TRIGGER IF NOT EXISTS
        {table}_insert AFTER INSERT ON {table}
        BEGIN
          INSERT INTO {table}_skytree VALUES (
            new.{objid},
            COS(new.{dec} * 0.017453292519943295)
            * COS(new.{ra} * 0.017453292519943295),
            COS(new.{dec} * 0.017453292519943295)
            * COS(new.{ra} * 0.017453292519943295),
            COS(new.{dec} * 0.017453292519943295)
            * SIN(new.{ra} * 0.017453292519943295),
            COS(new.{dec} * 0.017453292519943295)
            * SIN(new.{ra} * 0.017453292519943295),
            SIN(new.{dec} * 0.017453292519943295),
            SIN(new.{dec} * 0.017453292519943295)
          );
        END'''.format(table=self.table.name, objid=self.table.objid,
                      ra=self.table.ra, dec=self.table.dec))

    @abstractmethod
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
        pass

    def search(self, sources):
        """Search catalog for objects near sources.

        Parameters
        ----------
        sources : SkyCoord
            Sources to be matched.


        Returns
        -------
        objids : array
            Catalog IDs.

        cat : SkyCoord
            Coordinates for `objids`.

        """

        ra = sources.ra.rad
        dec = sources.dec.rad
        mean_ra = np.mean(ra)
        mean_dec = np.mean(dec)
        mean_xyz = (np.cos(mean_dec) * np.cos(mean_ra),
                    np.cos(mean_dec) * np.sin(mean_ra),
                    np.sin(mean_dec))
        xyz = np.array((np.cos(dec) * np.cos(ra),
                        np.cos(dec) * np.sin(ra),
                        np.sin(dec)))
        box = (max(xyz[0].max(), mean_xyz[0]),
               min(xyz[0].min(), mean_xyz[0]),
               max(xyz[1].max(), mean_xyz[1]),
               min(xyz[1].min(), mean_xyz[1]),
               max(xyz[2].max(), mean_xyz[2]),
               min(xyz[2].min(), mean_xyz[2]))

        rows = self.db.execute('''
        SELECT {objid},{ra},{dec} FROM {table}
        INNER JOIN {table}_skytree USING ({objid})
        WHERE x0 < ?
          AND x1 > ?
          AND y0 < ?
          AND y1 > ?
          AND z0 < ?
          AND z1 > ?
        '''.format(table=self.table.name, objid=self.table.objid,
                   ra=self.table.ra, dec=self.table.dec), box
        ).fetchall()

        if len(rows) == 0:
            return [], SkyCoord([], [], unit='deg')

        objids, ra, dec = [np.array(x) for x in zip(*rows)]
        cat = SkyCoord(ra, dec, unit='deg')
        return objids, cat

    def xmatch(self, sources):
        """Cross-match sources to catalog.

        Parameters
        ----------
        sources : SkyCoord
            Sources to cross-match.


        Returns
        -------
        matched : masked array
            Matched object ID of each source.

        d : masked array
            Distance between match and source in arcsec.

        """

        objids, cat = self.search(sources)

        if len(cat) == 0:
            self.logger.error('No catalog sources to match.')
            return

        idx, d2d = sources.match_to_catalog_sky(cat)[:2]
        i = d2d < self.match_limit
        n = i.sum()
        if n < self.min_matches:
            self.logger.error(
                'Fewer than {} sources matched: {}.'.format(self.min_matches, n))
            return

        self.logger.info(
            'Matched {} sources to photometric catalog.'.format(n))

        matched = np.ma.MaskedArray(objids[idx], int)
        matched.mask = ~i

        d = np.ma.MaskedArray(d2d.arcsec, float)
        d.mask = ~i

        return matched, d

    def lookup(self, objids, columns, allow_null=False):
        """Lookup these columns for these objects.

        Parameters
        ----------
        objids : array-like
            Catalog object IDs.  May be a `~numpy.ma.MaskedArray`.

        columns : string
            Database columns to return.

        allow_null : bool, optional
            Set to ``True`` if null values are allowed.

        Returns
        -------
        rows : list
            List of rows from the catalog.  If an object was not
            matched, the row will be an empty list.

        """

        statement = '''
        SELECT {columns} from {table}
        WHERE {objid}=?
        '''.format(columns=columns, table=self.table.name,
                   objid=self.table.objid)

        rows = []
        for i, objid in enumerate(objids):
            if objid is np.ma.masked or objid is None:
                row = []
            else:
                row = self.db.execute('''
                SELECT {columns} FROM {table}
                WHERE {objid}=?
                '''.format(
                    objid=self.table.objid,
                    table=self.table.name,
                    columns=columns
                ), [objid]).fetchone()

            if None in row and not allow_null:
                rows.append([])
            else:
                rows.append(row)

        return rows

    def cal_constant(self, matched, m_inst, filt, mlim=[14, 18],
                     gmi_lim=None):
        """Estimate calibration constant without color correction.

        Parameters
        ----------
        matched : array-like
            Object IDs matched to sources.  May be a masked array.

        m_inst : array-like
            Instrumental magnitudes for each source in matched.

        filt : string
            Filter to calibrate to, e.g., 'r'.

        mlim : list, optional
            Only fit stars with this magnitude range in filter `filt`.

        gmi_lim : list, optional
            Only fit stars with this g-i color range, or `None` to
            skip test.

        Returns
        -------
        zp_mean, zp_median, unc : float
            Zero-point magnitude mean, median, and uncertainty.

        m : float
            Catalog magnitude.

        gmi : ndarray
            g-i color for each source.

        """

        if filt not in self.table.filter2col:
            raise ValueError('Filter must be one of {}.'.format(
                self.table.filter2col.keys()))

        if gmi_lim is None:
            limits = [-np.inf, np.inf, min(mlim), max(mlim)]
        else:
            limits = [min(gmi_lim), max(gmi_lim), min(mlim), max(mlim)]

        columns = ("{filt[mag]},{filt[err]},{g[mag]}-{i[mag]}").format(
            filt=self.table.filter2col[filt],
            g=self.table.filter2col['g'],
            i=self.table.filter2col['i'])
        cat = self.lookup(matched, columns)

        m = np.ma.MaskedArray(np.zeros(len(matched)),
                              mask=np.ones(len(matched), bool))
        gmi = np.zeros_like(m.data)
        for i in range(len(cat)):
            if len(cat[i]) > 0:
                m[i], merr, gmi[i] = cat[i]
                if all((gmi[i] >= limits[0], gmi[i] <= limits[1],
                        m[i] >= limits[2], m[i] <= limits[3],
                        m[i] / merr > 2)):
                    m.mask[i] = False
                else:
                    m.mask[i] = True

        if np.all(m.mask):
            raise CalibrationError(
                'No data returned from database.  Check `matched` and catalog '
                'coverage of requested field.'
            )

        dm = m - m_inst
        i = np.isfinite(dm) * ~m.mask
        mms = sigma_clipped_stats(dm[i])
        return mms[0], mms[1], mms[2], m, gmi

    def cal_color(self, matched, m_inst, filt, color, C=None,
                  mlim=[14, 18], gmi_lim=[0.2, 3.0]):
        """Estimate calibration constant with color correction.

        Parameters
        ----------
        matched : array-like
            Object IDs matched to sources.  May be a
            `~numpy.ma.MaskedArray`.

        m_inst : array-like
            Instrumental magnitudes for each source in matched.

        filt : string
            Filter to calibrate to, e.g., 'r'.

        color : string
            Color to consider, e.g., 'g-r'.

        C : float, optional
            Set to a value to hold the color correction fixed.

        mlim : list, optional
            Only fit stars with this magnitude range in filter ``filt``.

        gmi_lim : list, optional
            Only fit stars with this g-i color range, or `None` to disable.

        Returns
        -------
        zp, C, unc : float
            Zero-point magnitude, color slope, and uncertainty.
                m - m_inst = C * color + zp

        m, cindex : MaskedArray
            Catalog magnitude and color index.

        gmi : ndarray
            g-i color for each source.

        """

        if filt not in self.table.filter2col:
            raise ValueError('Filter must be one of {}.'.format(
                self.table.filter2col.keys()))

        blue, red = color.split('-')
        if gmi_lim is None:
            limits = [-np.inf, np.inf, min(mlim), max(mlim)]
        else:
            limits = [min(gmi_lim), max(gmi_lim), min(mlim), max(mlim)]

        columns = ("{filt[mag]},{filt[err]},{b[mag]}-{r[mag]},"
                   "{g[mag]}-{i[mag]}").format(
                       filt=self.table.filter2col[filt],
                       b=self.table.filter2col[blue],
                       r=self.table.filter2col[red],
                       g=self.table.filter2col['g'],
                       i=self.table.filter2col['i'])
        cat = self.lookup(matched, columns)

        m = np.ma.MaskedArray(np.zeros(len(matched)),
                              mask=np.ones(len(matched), bool))
        gmi = np.zeros_like(m.data)
        cindex = m.copy()
        for i in range(len(cat)):
            if len(cat[i]) > 0:
                m[i], merr, cindex[i], gmi[i] = cat[i]
                if all((gmi[i] >= limits[0], gmi[i] <= limits[1],
                        m[i] >= limits[2], m[i] <= limits[3],
                        m[i] / merr > 2)):
                    m.mask[i] = False
                    cindex.mask[i] = False
                else:
                    m.mask[i] = True
                    cindex.mask[i] = True

        if np.all(m.mask):
            raise CalibrationError(
                'No data returned from database.  Check `matched` and catalog '
                'coverage of requested field.'
            )
        dm = m - m_inst

        model = models.Linear1D(slope=0, intercept=28)
        if C is not None:
            model.slope.value = C
            model.slope.fixed = True

        fitter = fitting.FittingWithOutlierRemoval(
            fitting.LinearLSQFitter(), sigma_clip)

        i = np.isfinite(dm) * ~dm.mask
        if sum(i) == 0:
            raise CalibrationError('All sources masked.')

        if (astropy_version[0] > 3 or
                (astropy_version[0] == 3 and astropy_version[1] >= 1)):
            # Return order changed in astropy 3.1
            # (http://docs.astropy.org/en/stable/changelog.html#id10)
            # Also now returns a boolean mask array rather than a
            # MaskedArray of the data which could be applied back to
            # reconstuct 'line' if needed (not currently used)
            fit, mask = fitter(model, cindex[i], dm[i])
        else:
            line, fit = fitter(model, cindex[i], dm[i])
        C = fit.slope.value
        zp = fit.intercept.value
        cal_unc = sigma_clipped_stats((dm - fit(cindex))[i])[2]

        return zp, C, cal_unc, m, cindex, gmi
