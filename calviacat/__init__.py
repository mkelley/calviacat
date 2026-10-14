# Licensed with the MIT License, see LICENSE for details
from importlib.metadata import PackageNotFoundError
from importlib.metadata import version as _version

from .catalog import CalibrationError, Catalog
from .gaia import Gaia
from .panstarrs1 import PanSTARRS1
from .refcat2 import RefCat2
from .skymapper import SkyMapper

try:
    __version__ = _version(__name__)
except PackageNotFoundError:  # pragma: no cover
    pass

__all__ = ["CalibrationError", "Catalog", "Gaia", "PanSTARRS1", "RefCat2", "SkyMapper"]
