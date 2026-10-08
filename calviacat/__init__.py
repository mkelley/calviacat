# Licensed with the MIT License, see LICENSE for details
from importlib.metadata import PackageNotFoundError
from importlib.metadata import version as _version

from .catalog import CalibrationError, Catalog
from .panstarrs1 import PanSTARRS1
from .skymapper import SkyMapper

try:
    from .refcat2 import RefCat2
except ImportError:
    pass

from .gaia import Gaia

try:
    __version__ = _version(__name__)
except PackageNotFoundError:
    pass

__all__ = ["CalibrationError", "Catalog", "Gaia", "PanSTARRS1", "RefCat2", "SkyMapper"]
