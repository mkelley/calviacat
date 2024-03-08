# Licensed with the MIT License, see LICENSE for details
from importlib.metadata import version as _version, PackageNotFoundError

from .catalog import *
from .panstarrs1 import PanSTARRS1
from .skymapper import SkyMapper
from .refcat2 import RefCat2
from .gaia import Gaia

try:
    __version__ = _version(__name__)
except PackageNotFoundError:
    pass
