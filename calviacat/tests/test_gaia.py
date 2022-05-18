import os
from unittest.mock import patch

import numpy as np
import pytest
import warnings
warnings.simplefilter("ignore", pytest.PytestUnknownMarkWarning)

from requests import HTTPError
import astropy.units as u
from astropy.table import Table
from astropy.tests.helper import assert_quantity_allclose

from calviacat.gaia import Gaia

FLUX2MAG = 2.5/np.log(10)

def data_path(filename):
    data_dir = os.path.join(os.path.dirname(__file__), 'data')
    return os.path.join(data_dir, filename)

class TestGaia:

    def test_mag_errors(self, tmpdir):

        db_path = tmpdir.join('test.db')
        gaia = Gaia(db_path)
        job_data_file = data_path('dr2_1.vot')
        table = Table.read(job_data_file)

        assert len(table.colnames) == 20
        flux_errors = u.Quantity([976.93274, 493.7522, 97.520424])
        assert_quantity_allclose(table['phot_g_mean_flux_over_error'], flux_errors)

        new_table = gaia._add_mag_errors(table)
        assert_quantity_allclose(table['phot_g_mean_flux_over_error'], FLUX2MAG/flux_errors)
