class CVCException(Exception):
    """Calviacat exception."""


class CalibrationError(CVCException):
    """Photometric calibration error."""


class DatabaseError(CVCException):
    """Error with the database."""


class XMatchError(CVCException):
    """Catalog cross-matching error."""
