try:
    from enum import StrEnum  # Python 3.11+
except ImportError:
    from backports.strenum import StrEnum  # Python 3.8 - 3.10
from collections import namedtuple
from types import MappingProxyType


R_STR = "([A-Z][a-z]?)"
R_INT = "(\d+)"
R_FLOAT = "(-?\d+\.\d+)"


FormatParams = namedtuple("FormatParams", ["format_trigger", "read_trigger", "regex_", "group_", "types_"])


FormatGroup = namedtuple("FormatGroup", ["name", "x", "y", "z"])


FORMATS = MappingProxyType({
    "orca": FormatParams(
        "* O   R   C   A *",
        "CARTESIAN COORDINATES (ANGSTROEM)",
        rf"^\s*{R_STR}\s+{R_FLOAT}\s+{R_FLOAT}\s+{R_FLOAT}\s*$",
        FormatGroup(name=1, x=2, y=3, z=4),
        FormatGroup(name=str, x=float, y=float, z=float)
    ),
    "gaussian": FormatParams(
        "Gaussian",
        "Standard orientation",
        rf"^\s*{R_INT}\s+{R_INT}\s+{R_INT}\s+{R_FLOAT}\s+{R_FLOAT}\s+{R_FLOAT}\s*$",
        FormatGroup(name=2, x=4, y=5, z=6),
        FormatGroup(name=int, x=float, y=float, z=float)
    ),
    "xyz": FormatParams(
        None,
        None,
        rf"^\s*{R_STR}\s+{R_FLOAT}\s+{R_FLOAT}\s+{R_FLOAT}\s*$",
        FormatGroup(name=1, x=2, y=3, z=4),
        FormatGroup(name=str, x=float, y=float, z=float),
    )
})


class FormatEnum(StrEnum):
    @property
    def format_trigger(self):
        return FORMATS[self].format_trigger
    @property
    def read_trigger(self):
        return FORMATS[self].read_trigger
    @property
    def regex_(self):
        return FORMATS[self].regex_
    @property
    def group_(self):
        return FORMATS[self].group_
    @property
    def types_(self):
        return FORMATS[self].types_


FileFormat = FormatEnum("FileFormat", ["XYZ", "GAUSSIAN", "ORCA"])


