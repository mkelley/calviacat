class CVCException(Exception):
    """Calviacat exception."""


class CalibrationError(CVCException):
    """Photometric calibration error."""


class XMatchError(CVCException):
    """Catalog cross-matching error."""
